#include "time_stepper.hpp"
#include "trajectory_recorder.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// Clamp the adaptive proposal to the configured admissible interval (only for adaptive timestepping)
Real clamp_dt(Real value, Real dtMin, Real dtMax)
{
    return std::max(dtMin, std::min(value, dtMax));
}

// Contact solves can demand a stricter time step than free flight
Real contact_target_dt(const TimeStepper& stepper, Real baseDt)
{
    return std::min(baseDt, stepper.contactSolveDt);
}

// adaptive dt controller used after each accepted step
// The proposal is based on:
// - translational travel relative to DEM cell size
// - rotational increment per step
// - proximity to first contact
// - active contact state
Real propose_next_dt(const TimeStepper& stepper,
    const Terrain& dem,
    Real minPositiveGap,
    std::uint32_t contactCount,
    const Vec6& velocity,
    Real baseDt)
{
    if (!stepper.adaptive) {
        return baseDt;
    }

    const Real tiny = Real(1e-12);
    const Real dtMin = std::max(baseDt * stepper.minDtRatio, tiny);
    const Real dtContact = std::max(contact_target_dt(stepper, baseDt), dtMin);
    Real nextDt = baseDt;

    const Real vNorm = velocity.head<3>().norm();
    if (vNorm > tiny) {
        nextDt = std::min(nextDt,
            (stepper.maxTranslationFrac * dem.cell) / vNorm);
    }

    const Real omegaNorm = velocity.tail<3>().norm();
    if (omegaNorm > tiny) {
        nextDt = std::min(nextDt,
            stepper.maxRotationPerStep / omegaNorm);
    }

    if (std::isfinite(minPositiveGap) && vNorm > tiny) {
        const Real nearBand = stepper.contactBandFrac * dem.cell;
        if (minPositiveGap < nearBand) {
            nextDt = std::min(nextDt,
                stepper.nearContactDtScale * minPositiveGap / vNorm);
            if (stepper.fixedNearContactDt) {
                nextDt = std::min(nextDt, dtContact);
            }
        }
    }

    if (contactCount > 0 && stepper.fixedContactDt) {
        nextDt = std::min(nextDt, dtContact);
    }

    return clamp_dt(nextDt, dtMin, baseDt);
}

// transforms 6-dimensional generalised velocities (vx,vy,vz,wx,wy,wz) 
// into 7-dimensional derivative of q d/dt(x,y,z,e0,e1,e2,e3)
static Eigen::Matrix<Real, 7, 1> qdot(const Eigen::Matrix<Real, 7, 1>& q,
    const Vec6& u)
{
    Eigen::Matrix<Real, 7, 1> qprime;

    const Real wx = u(3), wy = u(4), wz = u(5);

    // quaternion
    const Real e0 = q(3);
    const Real e1 = q(4);
    const Real e2 = q(5);
    const Real e3 = q(6);

    qprime(0) = u(0);
    qprime(1) = u(1);
    qprime(2) = u(2);
    qprime(3) = Real(-0.5) * (e1 * wx + e2 * wy + e3 * wz);
    qprime(4) = Real(0.5) * (e0 * wx - e3 * wy + e2 * wz);
    qprime(5) = Real(0.5) * (e3 * wx + e0 * wy - e1 * wz);
    qprime(6) = Real(0.5) * (-e2 * wx + e1 * wy + e0 * wz);

    return qprime;
}

// Shared stepping kernel used by both the in-memory and recorder-based APIs.
// The scheme follows the same high-level structure throughout:
// 1. build midpoint state qM
// 2. query gap at the midpoint
// 3. solve contact impulses with JOR-prox if needed
// 4. advance to the end state
// 5. accept if converged, record, and continue to the next dt
TimeStepperRunSummary run_core(const TimeStepper& stepper,
    Real mu,
    const Terrain& dem,
    const Rock& rock,
    Real dt,
    TimeStepperWorkspace& workspace,
    const TimeStepperStopCriteria& stop,
    TimeStepperResult* result,
    TrajectoryRecorder* recorder)
{
    const Real tStart = Real(0);
    const Real tEnd = stop.maxTime;

    TimeStepperRunSummary summary;
    summary.tFinal = tStart;
    summary.qFinal = rock.q0();
    summary.uFinal = rock.u0();
    summary.minGapFinal.reset();
    summary.stopReason = TimeStepperStopReason::Unknown;

    // reset output containers before run starts
    if (result) {
        result->t.clear();
        result->q.clear();
        result->u.clear();
        result->stopReason = TimeStepperStopReason::Unknown;
    }

    // check for invalid input
    if (!(dt > Real(0)) || !(tEnd > tStart)) {
        summary.stopReason = TimeStepperStopReason::InvalidInput;
        if (result) {
            result->t.push_back(tStart);
            result->q.push_back(rock.q0());
            result->u.push_back(rock.u0());
            result->stopReason = summary.stopReason;
        }
        return summary;
    }

    const Real baseDt = dt;

    // only important if adaptive timestepping is chosen
    const Real dtMin = stepper.adaptive
        ? std::max(baseDt * stepper.minDtRatio, Real(1e-12))
        : baseDt;
    const Real dtContact = std::max(contact_target_dt(stepper, baseDt), dtMin);
    const Real allowedPenetration = stepper.maxPenetrationFrac * dem.cell;

    // preallocate approximate space to avoid repeated reallocations
    if (result) {
        const std::uint32_t reserveGuess = std::max<std::uint32_t>(
            static_cast<std::uint32_t>(std::ceil(stop.maxTime / baseDt)) + 1,
            16u);
        result->t.reserve(reserveGuess);
        result->q.reserve(reserveGuess);
        result->u.reserve(reserveGuess);
        // store initial state
        result->t.push_back(tStart);
        result->q.push_back(rock.q0());
        result->u.push_back(rock.u0());
    }

    Real t = tStart;
    std::uint32_t stepIdx = 0;
    // estimate how many steps we expect if we roughly advance with baseDt until maxTime
    const std::uint32_t reserveGuess = std::max<std::uint32_t>(
        static_cast<std::uint32_t>(std::ceil(stop.maxTime / baseDt)) + 1,
        16u);
    // safety limit for adaptive stepping
    const std::uint32_t fallbackMaxSteps = std::max<std::uint32_t>(
        reserveGuess,
        reserveGuess * stepper.maxAdaptiveStepsFactor);
    const std::uint32_t maxSteps = stop.maxSteps > 0
        ? stop.maxSteps
        : fallbackMaxSteps;

    Vec6 uA = rock.u0();
    Eigen::Matrix<Real, 7, 1> qA = rock.q0();

    ArrX PN = ArrX::Zero(static_cast<Eigen::Index>(rock.nVert()));
    MatX PT = MatX::Zero(static_cast<Eigen::Index>(rock.nVert()), 2);

    const Mat6 MM_inv = rock.massMatrix().inverse();
    Real stepDt = baseDt;

    workspace.terrain.ensureVertexCapacity(rock.nVert());
    workspace.terrain.ensureContactCapacity(rock.nVert());
    workspace.solver.ensureContactCapacity(rock.nVert());

    while (t < tEnd && stepIdx < maxSteps)
    {
        // Never step beyond the requested final simulation time.
        stepDt = std::min(stepDt, tEnd - t);
        bool accepted = false;
        bool terminate = false;

        for (std::uint32_t refine = 0; refine <= stepper.maxStepRefinements; ++refine)
        {
            // Evaluate contact at the midpoint of the current tentative step.
            const auto qM = qA + qdot(qA, uA) * (stepDt * Real(0.5));
            Vec6 hM = Vec6::Zero();
            hM(2) = -rock.mass() * Real(9.81);
            const Vec6 rhs = hM * stepDt;

            Quat quatK(qM(3), qM(4), qM(5), qM(6));
            quatK.normalize();
            const Mat3 AIK = quatK.toRotationMatrix();
            const Vec3 xM = qM.head<3>();

            // queryPose returns:
            // - bounds / nodata classification
            // - all vertex gaps
            // - midpoint contact set
            // - DEM / roughness gradients only for contact vertices
            const TerrainQueryResult query = dem.queryPose(rock, AIK, xM, workspace.terrain);

            // out of bounds check
            if (query.state != TerrainContainment::Inside) {
                if (query.state == TerrainContainment::OutOfBounds
                    && stepper.adaptive && stepDt > dtMin * Real(1.01) && (t + stepDt) < tEnd) {
                    stepDt = std::max(stepDt * Real(0.5), dtMin);
                    continue;
                }
                summary.stopReason = (query.state == TerrainContainment::NoData)
                    ? TimeStepperStopReason::NoData
                    : TimeStepperStopReason::OutOfBounds;
                terminate = true;
                break;
            }

            // adaptive time stepping logic
            const Real maxPenetration = std::max(Real(0), -query.minGap);
            if (stepper.adaptive && maxPenetration > allowedPenetration && stepDt > dtMin * Real(1.01)) {
                stepDt = std::max(stepDt * Real(0.5), dtMin);
                continue;
            }
            const bool nearContact = std::isfinite(query.minPositiveGap)
                && query.minPositiveGap < stepper.contactBandFrac * dem.cell;
            if (stepper.adaptive && stepper.fixedNearContactDt && nearContact && stepDt > dtContact * Real(1.01)) {
                stepDt = dtContact;
                continue;
            }
            if (stepper.adaptive && stepper.fixedContactDt && query.contactCount > 0 && stepDt > dtContact * Real(1.01)) {
                stepDt = dtContact;
                continue;
            }

            Vec6 uE;
            // solve contact problem
            if (query.contactCount > 0)
            {
                // Build normal / tangential contact operators for the active
                // contact vertices only, then solve the proximal contact step.
                make_Wmat(rock,
                    workspace.terrain.contactIdx.data(),
                    query.contactCount,
                    AIK,
                    workspace.terrain.dhdx.data(),
                    workspace.terrain.dhdy.data(),
                    workspace.solver.WN,
                    workspace.solver.WT1,
                    workspace.solver.WT2);

                const Eigen::Index nC = static_cast<Eigen::Index>(query.contactCount);
                for (Eigen::Index j = 0; j < nC; ++j) {
                    const Eigen::Index v = static_cast<Eigen::Index>(workspace.terrain.contactIdx(j));
                    workspace.solver.PNsub(j) = PN(v);
                    workspace.solver.PTsub(j, 0) = PT(v, 0);
                    workspace.solver.PTsub(j, 1) = PT(v, 1);
                }

                uE = jor_project(workspace.solver, query.contactCount, MM_inv, rhs, uA,
                    mu, stepper.alpha, stepper.atol, stepper.rtol, stepper.maxiter);

                for (Eigen::Index j = 0; j < nC; ++j) {
                    const Eigen::Index v = static_cast<Eigen::Index>(workspace.terrain.contactIdx(j));
                    PN(v) = workspace.solver.PNsub(j);
                    PT(v, 0) = workspace.solver.PTsub(j, 0);
                    PT(v, 1) = workspace.solver.PTsub(j, 1);
                }
            }
            else {
                // all contacts are open
                uE = uA + MM_inv * rhs;
                PN.setZero();
                PT.setZero();
            }

            // Trapezoidal-style state update using midpoint kinematics and the
            // average of old/new generalized velocities.
            Eigen::Matrix<Real, 7, 1> qE =
                (qA + qdot(qM, (uA + uE) * Real(0.5)) * stepDt).eval();
            Quat qc(qE(3), qE(4), qE(5), qE(6));
            qc.normalize();
            qE.segment<4>(3) << qc.w(), qc.x(), qc.y(), qc.z();

            t += stepDt;
            ++stepIdx;

            summary.tFinal = t;
            summary.acceptedSteps = stepIdx;
            summary.qFinal = qE;
            summary.uFinal = uE;
            summary.minGapFinal = query.minGap;

            if (result) {
                result->t.push_back(t);
                result->q.push_back(qE);
                result->u.push_back(uE);
            }
            // dump time step result to binary file
            if (recorder) {
                recorder->onAcceptedStep(RecordedStepData{ t, qE, uE, query.minGap });
            }

            qA.swap(qE);
            uA.swap(uE);
            stepDt = propose_next_dt(stepper, dem, query.minPositiveGap, query.contactCount, uA, baseDt);
            accepted = true;
            break;
        }

        if (terminate || !accepted) {
            break;
        }

        // Low-energy settling stop. This is intentionally based on the norms
        // of the accepted generalized velocity, not on contact state.
        if (stepIdx >= stop.settleMinSteps
            && uA.head<3>().norm() < stop.settleLinearSpeed
            && uA.tail<3>().norm() < stop.settleAngularSpeed) {
            summary.stopReason = TimeStepperStopReason::Energy;
            break;
        }
    }

    if (summary.stopReason == TimeStepperStopReason::Unknown) {
        if (stepIdx >= maxSteps) {
            summary.stopReason = TimeStepperStopReason::MaxSteps;
        }
        else if (t >= tEnd) {
            summary.stopReason = TimeStepperStopReason::MaxTime;
        }
    }

    if (result) {
        result->stopReason = summary.stopReason;
    }

    return summary;
}

} // namespace

TimeStepperResult
TimeStepper::run(Real mu, const Terrain& dem,
    const Rock& rock, Real dt,
    const TimeStepperStopCriteria& stop) const
{
    TimeStepperWorkspace workspace;
    return run(mu, dem, rock, dt, workspace, stop);
}

TimeStepperResult
TimeStepper::run(Real mu, const Terrain& dem,
    const Rock& rock, Real dt,
    TimeStepperWorkspace& workspace,
    const TimeStepperStopCriteria& stop) const
{
    TimeStepperResult result;
    run_core(*this, mu, dem, rock, dt, workspace, stop, &result, nullptr);
    return result;
}

TimeStepperRunSummary
TimeStepper::runRecorded(Real mu, const Terrain& dem,
    const Rock& rock, Real dt,
    TrajectoryRecorder& recorder,
    const TimeStepperStopCriteria& stop) const
{
    TimeStepperWorkspace workspace;
    return runRecorded(mu, dem, rock, dt, workspace, recorder, stop);
}

TimeStepperRunSummary
TimeStepper::runRecorded(Real mu, const Terrain& dem,
    const Rock& rock, Real dt,
    TimeStepperWorkspace& workspace,
    TrajectoryRecorder& recorder,
    const TimeStepperStopCriteria& stop) const
{
    return run_core(*this, mu, dem, rock, dt, workspace, stop, nullptr, &recorder);
}
