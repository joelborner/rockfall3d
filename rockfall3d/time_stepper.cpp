#include "time_stepper.hpp"
#include "trajectory_recorder.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// Transforms 6-dimensional generalised velocities (vx,vy,vz,wx,wy,wz)
// into the 7-dimensional derivative d/dt(x,y,z,e0,e1,e2,e3).
static Eigen::Matrix<Real, 7, 1> qdot(const Eigen::Matrix<Real, 7, 1>& q,
    const Vec6& u)
{
    Eigen::Matrix<Real, 7, 1> qprime;

    const Real wx = u(3), wy = u(4), wz = u(5);

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

// Shared fixed-step kernel used by both the in-memory and recorder-based APIs.
// The scheme follows the same high-level structure throughout:
// 1. build midpoint state qM
// 2. query gap at the midpoint
// 3. solve contact impulses with the configured prox solver if needed
// 4. advance to the end state
// 5. accept, record, and continue with the next fixed step
TimeStepperRunSummary run_core(const TimeStepper& stepper,
    Real mu,
    const Terrain& dem,
    const Rock& rock,
    Real dt,
    TimeStepperWorkspace& workspace,
    const TimeStepperStopCriteria& stop,
    TimeStepperResult* result,
    TrajectoryRecorder* recorder,
    FragmentationShapeWriter* shapeWriter)
{
    const Real tStart = Real(0);
    const Real tEnd = stop.maxTime;

    TimeStepperRunSummary summary;
    summary.tFinal = tStart;
    summary.qFinal = rock.q0();
    summary.uFinal = rock.u0();
    summary.minGapFinal.reset();
    summary.massFinal = rock.mass();
    summary.inertiaFinal = rock.inertiaK();
    summary.stopReason = TimeStepperStopReason::Unknown;

    if (result) {
        result->t.clear();
        result->q.clear();
        result->u.clear();
        result->stopReason = TimeStepperStopReason::Unknown;
    }

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

    const std::uint32_t reserveGuess = std::max<std::uint32_t>(
        static_cast<std::uint32_t>(std::ceil(stop.maxTime / dt)) + 1,
        16u);
    const std::uint32_t maxSteps = stop.maxSteps > 0
        ? stop.maxSteps
        : reserveGuess;

    if (result) {
        result->t.reserve(reserveGuess);
        result->q.reserve(reserveGuess);
        result->u.reserve(reserveGuess);
        result->t.push_back(tStart);
        result->q.push_back(rock.q0());
        result->u.push_back(rock.u0());
    }

    Real t = tStart;
    std::uint32_t stepIdx = 0;

    Vec6 uA = rock.u0();
    Eigen::Matrix<Real, 7, 1> qA = rock.q0();

	// initiator for the active rock, which may be replaced by a fragmented rock during the simulation
    const Rock* activeRock = &rock;
    std::optional<Rock> fragmentedRock;
    std::uint32_t shapeIndex = 0;

    ArrX PN = ArrX::Zero(static_cast<Eigen::Index>(activeRock->nVert()));
    MatX PT = MatX::Zero(static_cast<Eigen::Index>(activeRock->nVert()), 2);

    Mat6 MM_inv = activeRock->massMatrix().inverse();

    workspace.terrain.ensureVertexCapacity(rock.nVert());
    workspace.terrain.ensureContactCapacity(rock.nVert());
    workspace.solver.ensureContactCapacity(rock.nVert());

    while (t < tEnd && stepIdx < maxSteps)
    {
        const Real stepDt = std::min(dt, tEnd - t);

        // Evaluate contact at the midpoint of the fixed step.
        const auto qM = qA + qdot(qA, uA) * (stepDt * Real(0.5));
        Vec6 hM = Vec6::Zero();
        hM(2) = -activeRock->mass() * Real(9.81);
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
        const TerrainQueryResult query = dem.queryPose(*activeRock, AIK, xM, workspace.terrain);

        if (query.state != TerrainContainment::Inside) {
            summary.stopReason = (query.state == TerrainContainment::NoData)
                ? TimeStepperStopReason::NoData
                : TimeStepperStopReason::OutOfBounds;
            break;
        }

        Vec6 uE;
        if (query.contactCount > 0)
        {
            // Build normal / tangential contact operators for the active
            // contact vertices only, then solve the proximal contact step.
            make_Wmat(*activeRock,
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

            bool converged = false;
            if (stepper.solverType == ContactSolverType::SorProx) {
                const SorProjectResult solve = sor_project(workspace.solver, query.contactCount, MM_inv, rhs, uA,
                    mu, stepper.alpha, stepper.epsilonN, stepper.epsilonT,
                    stepper.atol, stepper.rtol, stepper.maxiter);
                converged = solve.converged;
                uE = solve.velocity;
            }
            else {
                const JorProjectResult solve = jor_project(workspace.solver, query.contactCount, MM_inv, rhs, uA,
                    mu, stepper.alpha, stepper.epsilonN, stepper.epsilonT,
                    stepper.atol, stepper.rtol, stepper.maxiter);
                converged = solve.converged;
                uE = solve.velocity;
            }

            if (!converged) {
                summary.stopReason = TimeStepperStopReason::SolverNonConvergence;
                break;
            }

            for (Eigen::Index j = 0; j < nC; ++j) {
                const Eigen::Index v = static_cast<Eigen::Index>(workspace.terrain.contactIdx(j));
                PN(v) = workspace.solver.PNsub(j);
                PT(v, 0) = workspace.solver.PTsub(j, 0);
                PT(v, 1) = workspace.solver.PTsub(j, 1);
            }
        }
        else {
            uE = uA + MM_inv * rhs;
            PN.setZero();
            PT.setZero();
        }

        if (summary.stopReason == TimeStepperStopReason::SolverNonConvergence) {
            break;
        }

        // Trapezoidal-style state update using midpoint kinematics and the
        // average of old/new generalized velocities.
        Eigen::Matrix<Real, 7, 1> qE =
            (qA + qdot(qM, (uA + uE) * Real(0.5)) * stepDt).eval();
        Quat qc(qE(3), qE(4), qE(5), qE(6));
        qc.normalize();
        qE.segment<4>(3) << qc.w(), qc.x(), qc.y(), qc.z();

		// Fragmentation check. If fragmentation is enabled and there are active contacts, compute the dissipated work at each contact and determine if fragmentation should occur.
        if (stepper.fragmentation && query.contactCount > 0) {
            std::vector<FragmentationContact> fragmentationContacts;
            fragmentationContacts.reserve(query.contactCount);
            const Eigen::Index nC = static_cast<Eigen::Index>(query.contactCount);
            for (Eigen::Index j = 0; j < nC; ++j) {
                const Real gammaMinus = workspace.solver.WN.col(j).dot(uA);
                const Real gammaPlus = workspace.solver.WN.col(j).dot(uE);
                const Real impulse = workspace.solver.PNsub(j);
                const Real work = std::max(
                    Real(0), -Real(0.5) * impulse * (gammaMinus + gammaPlus));
                if (work > Real(0)) {
                    const Eigen::Index vertex = static_cast<Eigen::Index>(
                        workspace.terrain.contactIdx(j));
                    fragmentationContacts.push_back(FragmentationContact{
                        activeRock->verticesK().col(vertex),
                        work
                    });
                }
            }

			// If fragmentation occurs, update the active rock to the new fragmented rock and adjust the state accordingly.
            if (auto fragmented = fragment_rock(
                *activeRock, fragmentationContacts, *stepper.fragmentation)) {
                Quat endQuaternion(qE(3), qE(4), qE(5), qE(6));
                endQuaternion.normalize();
                const Mat3 endRotation = endQuaternion.toRotationMatrix();
                const Vec3 centerShiftK = fragmented->centerShiftK;
                qE.head<3>() += endRotation * centerShiftK;
                uE.head<3>() += endRotation * uE.tail<3>().cross(centerShiftK);

                fragmentedRock = fragmented->rock.withInitialState(uE, qE);
                activeRock = &*fragmentedRock;
                ++shapeIndex;
                if (shapeWriter) {
                    shapeWriter->write(shapeIndex, t + stepDt, *activeRock);
                }

                PN = ArrX::Zero(static_cast<Eigen::Index>(activeRock->nVert()));
                PT = MatX::Zero(static_cast<Eigen::Index>(activeRock->nVert()), 2);
                MM_inv = activeRock->massMatrix().inverse();
                workspace.terrain.ensureVertexCapacity(activeRock->nVert());
                workspace.terrain.ensureContactCapacity(activeRock->nVert());
                workspace.solver.ensureContactCapacity(activeRock->nVert());
            }
        }

        t += stepDt;
        ++stepIdx;

        summary.tFinal = t;
        summary.acceptedSteps = stepIdx;
        summary.qFinal = qE;
        summary.uFinal = uE;
        summary.minGapFinal = query.minGap;
        summary.massFinal = activeRock->mass();
        summary.inertiaFinal = activeRock->inertiaK();
        summary.shapeIndexFinal = shapeIndex;

        if (result) {
            result->t.push_back(t);
            result->q.push_back(qE);
            result->u.push_back(uE);
        }
        if (recorder) {
            recorder->onAcceptedStep(RecordedStepData{
                t, qE, uE, query.minGap,
                activeRock->mass(), activeRock->inertiaK(), shapeIndex
            });
        }

        qA.swap(qE);
        uA.swap(uE);

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
    run_core(*this, mu, dem, rock, dt, workspace, stop, &result, nullptr, nullptr);
    return result;
}

TimeStepperRunSummary
TimeStepper::runRecorded(Real mu, const Terrain& dem,
    const Rock& rock, Real dt,
    TrajectoryRecorder& recorder,
    const TimeStepperStopCriteria& stop,
    FragmentationShapeWriter* shapeWriter) const
{
    TimeStepperWorkspace workspace;
    return runRecorded(mu, dem, rock, dt, workspace, recorder, stop, shapeWriter);
}

TimeStepperRunSummary
TimeStepper::runRecorded(Real mu, const Terrain& dem,
    const Rock& rock, Real dt,
    TimeStepperWorkspace& workspace,
    TrajectoryRecorder& recorder,
    const TimeStepperStopCriteria& stop,
    FragmentationShapeWriter* shapeWriter) const
{
    return run_core(*this, mu, dem, rock, dt, workspace, stop, nullptr, &recorder, shapeWriter);
}
