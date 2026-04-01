#pragma once
#include "math_alias.hpp"
#include "rock.hpp"
#include "terrain.hpp"
#include "contact_wmat.hpp"
#include "jor_solver.hpp"

#include <optional>

class TrajectoryRecorder;

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
    InvalidInput = 6
};

inline const char* stop_reason_name(TimeStepperStopReason reason) noexcept
{
    switch (reason) {
    case TimeStepperStopReason::Energy: return "energy";
    case TimeStepperStopReason::OutOfBounds: return "out_of_bounds";
    case TimeStepperStopReason::NoData: return "nodata_encountered";
    case TimeStepperStopReason::MaxTime: return "max_time";
    case TimeStepperStopReason::MaxSteps: return "max_steps";
    case TimeStepperStopReason::InvalidInput: return "invalid_input";
    default: return "unknown";
    }
}

struct TimeStepperResult
{
    // Accepted states only. Rejected adaptive substeps never appear here.
    std::vector<Real>                           t;
    std::vector<Eigen::Matrix<Real, 7, 1>>     q;
    std::vector<Vec6>                           u;
    TimeStepperStopReason stopReason{ TimeStepperStopReason::Unknown };
};

struct TimeStepperStopCriteria
{
    // Default stop criteria. usually overwritten by parsing an input XML
    Real maxTime = Real(20.0);                 // hard upper bound on simulated time
    std::uint32_t maxSteps = 200000;           // hard upper bound on accepted steps
    std::uint32_t settleMinSteps = 100;        // ignore settle checks before this many steps
    Real settleLinearSpeed = Real(0.1);        // stop once translation falls below this value
    Real settleAngularSpeed = Real(0.1);       // stop once rotation falls below this value
};

struct TimeStepperRunSummary
{
    // Lightweight end-of-run summary for recorder-based output modes.
    Real tFinal{ Real(0) };
    std::uint32_t acceptedSteps{ 0 };
    Eigen::Matrix<Real, 7, 1> qFinal = (Eigen::Matrix<Real, 7, 1>() << 0, 0, 0, 1, 0, 0, 0).finished();
    Vec6 uFinal = Vec6::Zero();
    std::optional<Real> minGapFinal;
    TimeStepperStopReason stopReason{ TimeStepperStopReason::Unknown };
};

struct TimeStepperWorkspace
{
    // Per-thread scratch storage reused across trajectories. This avoids
    // repeated heap allocations in the hot terrain-query and contact-solve
    // path.
    TerrainQueryWorkspace terrain;
    ContactSolverWorkspace solver;
};

class TimeStepper
{
public:
    // Parameters of the JOR-Prox contact solver
    Real alpha = Real(0.5);
    Real atol = Real(1e-8);
    Real rtol = Real(1e-6);
    std::uint32_t maxiter = 10000;

    // Adaptive timestep controls. dt in run(...) is interpreted as the
    // nominal/free-flight step size, which may be reduced when needed.
    bool adaptive = true;
    Real minDtRatio = Real(0.05);              // lower bound relative to dt
    Real maxTranslationFrac = Real(0.35);      // max traveled cell fraction per step
    Real maxRotationPerStep = Real(0.15);      // rad, based on |omega| * dt
    Real contactBandFrac = Real(0.75);         // near-contact band in DEM cell units
    Real nearContactDtScale = Real(0.6);       // extra shrink factor close to contact
    Real contactSolveDt = Real(0.002);         // maximum dt while solving contact
    bool fixedContactDt = true;                // cap dt by contactSolveDt once contact is active
    bool fixedNearContactDt = true;            // cap dt by contactSolveDt in the near-contact band
    Real maxPenetrationFrac = Real(0.10);      // reject a step if midpoint penetration is too deep
    std::uint32_t maxStepRefinements = 8;      // max halvings for one attempted step
    std::uint32_t maxAdaptiveStepsFactor = 20; // fallback guard when stop.maxSteps is not tightened

    // Run and keep all dumped states in memory.
    TimeStepperResult run(Real mu,
        const Terrain& dem,
        const Rock& rock,
        Real dt,
        const TimeStepperStopCriteria& stop = {}) const;

    // Same as above, but caller provides reusable scratch storage.
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

