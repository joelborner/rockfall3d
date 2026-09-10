#include "time_stepper.hpp"
#include "trajectory_recorder.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// Shared fixed-step kernel used by both the in-memory and recorder-based APIs.
// The scheme follows the same high-level structure throughout:
// 1. build midpoint state qM
// 2. query gap at the midpoint
// 3. solve contact impulses with the configured prox solver if needed
// 4. advance to the end state
// 5. accept, record, and continue with the next fixed step
template <bool RollingEnabled, bool SlippageDependent, typename RotationalIntegrator>
TimeStepperRunSummary run_core(const TimeStepper& stepper,
    const RotationalIntegrator& integrator,
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
    TerrainInteractionState interactionState;

    const CoulombFrictionSettings& friction =
        stepper.terrainInteractions.coulombFriction;
    const Real constantMu = friction.mu;

    // The active rock is replaced after each accepted fragmentation event.
    const Rock* activeRock = &rock;
    std::optional<Rock> fragmentedRock;
    std::uint32_t shapeIndex = 0;

    ArrX PN = ArrX::Zero(static_cast<Eigen::Index>(activeRock->nVert()));
    MatX PT = MatX::Zero(static_cast<Eigen::Index>(activeRock->nVert()), 2);
    Real muR = Real(0);
    MatX PR;
    if constexpr (RollingEnabled) {
        muR = stepper.terrainInteractions.rollingFrictionCoefficient
            * activeRock->equivalentRadius();
        PR = MatX::Zero(static_cast<Eigen::Index>(activeRock->nVert()), 2);
        workspace.solver.ensureRollingCapacity(activeRock->nVert());
    }

    Mat3 inertiaInv = activeRock->inertiaK().inverse();

    workspace.terrain.ensureVertexCapacity(activeRock->nVert());
    workspace.terrain.ensureContactCapacity(activeRock->nVert());
    workspace.solver.ensureContactCapacity(activeRock->nVert());

    while (t < tEnd && stepIdx < maxSteps)
    {
        const Real stepDt = std::min(dt, tEnd - t);

        // Evaluate contact at the midpoint of the fixed step.
        const auto qM = make_midpoint_configuration(qA, uA, stepDt);

        Quat quatK(qM(3), qM(4), qM(5), qM(6));
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

        const IntegratorFreeStep freeStep = integrator.prepareFreeStep(
            activeRock->mass(), activeRock->inertiaK(), inertiaInv, uA, stepDt);
        if (!freeStep.converged) {
            summary.stopReason = TimeStepperStopReason::IntegratorNonConvergence;
            break;
        }

        Real mu = constantMu;
        if constexpr (SlippageDependent) {
            mu = friction.coefficient(interactionState.slippage);
        }

        Vec6 uE;
        bool hasActiveTerrainContact = false;
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

            if constexpr (RollingEnabled) {
                make_rolling_Wmat(query.contactCount,
                    AIK,
                    workspace.terrain.dhdx.data(),
                    workspace.terrain.dhdy.data(),
                    workspace.solver.WR1,
                    workspace.solver.WR2);
            }

            const Eigen::Index nC = static_cast<Eigen::Index>(query.contactCount);
            for (Eigen::Index j = 0; j < nC; ++j) {
                const Eigen::Index v = static_cast<Eigen::Index>(workspace.terrain.contactIdx(j));
                workspace.solver.PNsub(j) = PN(v);
                workspace.solver.PTsub(j, 0) = PT(v, 0);
                workspace.solver.PTsub(j, 1) = PT(v, 1);
                if constexpr (RollingEnabled) {
                    workspace.solver.PRsub(j, 0) = PR(v, 0);
                    workspace.solver.PRsub(j, 1) = PR(v, 1);
                }
            }

            bool converged = false;
            if (stepper.solverType == ContactSolverType::SorProx) {
                const auto solve = sor_project<RollingEnabled>(workspace.solver, query.contactCount,
                    freeStep.mobility, freeStep.velocity, uA,
                    mu, muR, stepper.alpha, stepper.epsilonN, stepper.epsilonT,
                    stepper.atol, stepper.rtol, stepper.maxiter);
                converged = solve.converged;
                uE = solve.velocity;
            }
            else {
                const auto solve = jor_project<RollingEnabled>(workspace.solver, query.contactCount,
                    freeStep.mobility, freeStep.velocity, uA,
                    mu, muR, stepper.alpha, stepper.epsilonN, stepper.epsilonT,
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
                hasActiveTerrainContact = hasActiveTerrainContact
                    || workspace.solver.PNsub(j) > Real(0);
                PT(v, 0) = workspace.solver.PTsub(j, 0);
                PT(v, 1) = workspace.solver.PTsub(j, 1);
                if constexpr (RollingEnabled) {
                    PR(v, 0) = workspace.solver.PRsub(j, 0);
                    PR(v, 1) = workspace.solver.PRsub(j, 1);
                }
            }
        }
        else {
            uE = freeStep.velocity;
            PN.setZero();
            PT.setZero();
            if constexpr (RollingEnabled) PR.setZero();
        }

        if constexpr (SlippageDependent) {
            // The old RAMMS implementation evaluates the midpoint COM speed
            // from bsm.p = bs.p, which is the accepted start velocity here.
            update_slippage_state(interactionState, friction,
                hasActiveTerrainContact, uA.head<3>().norm(), stepDt);
        }

        Eigen::Matrix<Real, 7, 1> qE = integrator.finishConfiguration(
            qA, qM, uA, uE, activeRock->inertiaK(), freeStep, stepDt);

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
                inertiaInv = activeRock->inertiaK().inverse();
                workspace.terrain.ensureVertexCapacity(activeRock->nVert());
                workspace.terrain.ensureContactCapacity(activeRock->nVert());
                workspace.solver.ensureContactCapacity(activeRock->nVert());
                if constexpr (RollingEnabled) {
                    muR = stepper.terrainInteractions.rollingFrictionCoefficient
                        * activeRock->equivalentRadius();
                    PR = MatX::Zero(static_cast<Eigen::Index>(activeRock->nVert()), 2);
                    workspace.solver.ensureRollingCapacity(activeRock->nVert());
                }
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

template <bool RollingEnabled, bool SlippageDependent>
TimeStepperRunSummary run_integrator_dispatched(const TimeStepper& stepper,
    const Terrain& dem,
    const Rock& rock,
    Real dt,
    TimeStepperWorkspace& workspace,
    const TimeStepperStopCriteria& stop,
    TimeStepperResult* result,
    TrajectoryRecorder* recorder,
    FragmentationShapeWriter* shapeWriter)
{
    switch (stepper.integrator.type) {
    case RotationalIntegratorType::Schweizer:
        return run_core<RollingEnabled, SlippageDependent>(
            stepper, SchweizerIntegrator{}, dem, rock, dt,
            workspace, stop, result, recorder, shapeWriter);
    case RotationalIntegratorType::Capobianco:
        return run_core<RollingEnabled, SlippageDependent>(stepper,
            CapobiancoIntegrator{
                stepper.integrator.atol,
                stepper.integrator.rtol,
                stepper.integrator.maxiter
            },
            dem, rock, dt, workspace, stop, result, recorder, shapeWriter);
    }

    TimeStepperRunSummary summary;
    summary.stopReason = TimeStepperStopReason::InvalidInput;
    return summary;
}

TimeStepperRunSummary run_dispatched(const TimeStepper& stepper,
    const Terrain& dem,
    const Rock& rock,
    Real dt,
    TimeStepperWorkspace& workspace,
    const TimeStepperStopCriteria& stop,
    TimeStepperResult* result,
    TrajectoryRecorder* recorder,
    FragmentationShapeWriter* shapeWriter)
{
    const bool rolling = stepper.terrainInteractions.rollingFrictionEnabled();
    const bool slippage =
        stepper.terrainInteractions.coulombFriction.slippageDependent();
    if (rolling) {
        if (slippage) {
            return run_integrator_dispatched<true, true>(stepper, dem, rock, dt,
                workspace, stop, result, recorder, shapeWriter);
        }
        return run_integrator_dispatched<true, false>(stepper, dem, rock, dt,
            workspace, stop, result, recorder, shapeWriter);
    }
    if (slippage) {
        return run_integrator_dispatched<false, true>(stepper, dem, rock, dt,
            workspace, stop, result, recorder, shapeWriter);
    }
    return run_integrator_dispatched<false, false>(stepper, dem, rock, dt,
        workspace, stop, result, recorder, shapeWriter);
}

} // namespace

TimeStepperResult
TimeStepper::run(const Terrain& dem,
    const Rock& rock, Real dt,
    const TimeStepperStopCriteria& stop) const
{
    TimeStepperWorkspace workspace;
    return run(dem, rock, dt, workspace, stop);
}

TimeStepperResult
TimeStepper::run(const Terrain& dem,
    const Rock& rock, Real dt,
    TimeStepperWorkspace& workspace,
    const TimeStepperStopCriteria& stop) const
{
    TimeStepperResult result;
    run_dispatched(*this, dem, rock, dt, workspace, stop, &result, nullptr, nullptr);
    return result;
}

TimeStepperRunSummary
TimeStepper::runRecorded(const Terrain& dem,
    const Rock& rock, Real dt,
    TrajectoryRecorder& recorder,
    const TimeStepperStopCriteria& stop,
    FragmentationShapeWriter* shapeWriter) const
{
    TimeStepperWorkspace workspace;
    return runRecorded(dem, rock, dt, workspace, recorder, stop, shapeWriter);
}

TimeStepperRunSummary
TimeStepper::runRecorded(const Terrain& dem,
    const Rock& rock, Real dt,
    TimeStepperWorkspace& workspace,
    TrajectoryRecorder& recorder,
    const TimeStepperStopCriteria& stop,
    FragmentationShapeWriter* shapeWriter) const
{
    return run_dispatched(*this, dem, rock, dt, workspace, stop,
        nullptr, &recorder, shapeWriter);
}
