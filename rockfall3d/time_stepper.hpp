#pragma once
#include "math_alias.hpp"
#include "rock.hpp"
#include "terrain.hpp"
#include "contact_wmat.hpp"
#include "jor_solver.hpp"
#include "rotational_integrator.hpp"
#include "sor_solver.hpp"

#include <optional>

class TrajectoryRecorder;

// Contact solver variants available in the fixed-step midpoint loop.
enum class ContactSolverType : std::uint8_t
{
    JorProx,
    SorProx
};

// Reasons why a trajectory stopped. These values are propagated to
// manifest.csv and therefore form part of the public run contract.
enum class TimeStepperStopReason : std::uint8_t
{
    Unknown = 0,
    Energy = 1,
    OutOfBounds = 2,
    NoData = 3,
    MaxTime = 4,
    MaxSteps = 5,
    InvalidInput = 6,
    SolverNonConvergence = 7,
    IntegratorNonConvergence = 8
};

// returns strings for the reached stop criterion
inline const char* stop_reason_name(TimeStepperStopReason reason) noexcept
{
    switch (reason) {
    case TimeStepperStopReason::Energy: return "energy";
    case TimeStepperStopReason::OutOfBounds: return "out_of_bounds";
    case TimeStepperStopReason::NoData: return "nodata_encountered";
    case TimeStepperStopReason::MaxTime: return "max_time";
    case TimeStepperStopReason::MaxSteps: return "max_steps";
    case TimeStepperStopReason::InvalidInput: return "invalid_input";
    case TimeStepperStopReason::SolverNonConvergence: return "solver_nonconvergence";
    case TimeStepperStopReason::IntegratorNonConvergence: return "integrator_nonconvergence";
    default: return "unknown";
    }
}

// Accepted fixed-step states.
struct TimeStepperResult
{
    std::vector<Real>                           t;// time
    std::vector<Eigen::Matrix<Real, 7, 1>>     q;// generalised coordinates
	std::vector<Vec6>                           u;// generalised velocities
    TimeStepperStopReason stopReason{ TimeStepperStopReason::Unknown };
};

struct TimeStepperStopCriteria
{
    Real maxTime = Real(20.0);                 // hard upper bound on simulated time
    std::uint32_t maxSteps = 200000;           // hard upper bound on accepted steps
    std::uint32_t settleMinSteps = 100;        // ignore settle checks before this many steps
    Real settleLinearSpeed = Real(0.1);        // stop once translation falls below this value
    Real settleAngularSpeed = Real(0.1);       // stop once rotation falls below this value
};

// Lightweight end-of-run summary for recorder-based output modes.
struct TimeStepperRunSummary
{
    Real tFinal{ Real(0) };// total simulated time
	std::uint32_t acceptedSteps{ 0 };// total accepted steps
    Eigen::Matrix<Real, 7, 1> qFinal = (Eigen::Matrix<Real, 7, 1>() << 0, 0, 0, 1, 0, 0, 0).finished();// generalized coordinates at stop
    Vec6 uFinal = Vec6::Zero();// generalized velocities at stop
    std::optional<Real> minGapFinal;// minimum gap at stop
    TimeStepperStopReason stopReason{ TimeStepperStopReason::Unknown };// reason for stopping
};

// Per-thread scratch storage reused across trajectories. This avoids
// repeated heap allocations in the hot terrain-query and contact-solve
// path.
struct TimeStepperWorkspace
{
    TerrainQueryWorkspace terrain;
    ContactSolverWorkspace solver;
};

class TimeStepper
{
public:
    RotationalIntegratorSettings integrator;

    // Contact solver variant and shared scalar iteration parameters.
    ContactSolverType solverType = ContactSolverType::JorProx;
    Real alpha = Real(0.5);
    Real atol = Real(1e-8);
    Real rtol = Real(1e-6);
    Real epsilonN = Real(0);
    Real epsilonT = Real(0);
    std::uint32_t maxiter = 10000;

    // Run and keep all dumped states in memory.
    TimeStepperResult run(Real mu,
        const Terrain& dem,
        const Rock& rock,
        Real dt,
        const TimeStepperStopCriteria& stop = {}) const;

    // Run and keep all dumped states in memory, caller provides reusable scratch storage.
    TimeStepperResult run(Real mu,
        const Terrain& dem,
        const Rock& rock,
        Real dt,
        TimeStepperWorkspace& workspace,
        const TimeStepperStopCriteria& stop = {}) const;

    // Run and stream accepted states directly to a recorder.
    TimeStepperRunSummary runRecorded(Real mu,
        const Terrain& dem,
        const Rock& rock,
        Real dt,
        TrajectoryRecorder& recorder,
        const TimeStepperStopCriteria& stop = {}) const;

    // Recorder-based run with caller-owned scratch storage.
    TimeStepperRunSummary runRecorded(Real mu,
        const Terrain& dem,
        const Rock& rock,
        Real dt,
        TimeStepperWorkspace& workspace,
        TrajectoryRecorder& recorder,
        const TimeStepperStopCriteria& stop = {}) const;
};

