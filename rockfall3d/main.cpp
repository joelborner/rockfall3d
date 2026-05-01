#include "config_io.hpp"
#include "fast_io.hpp"
#include "raster_output.hpp"
#include "rock.hpp"
#include "terrain.hpp"
#include "terrain_export.hpp"
#include "terrain_io.hpp"
#include "time_stepper.hpp"
#include "trajectory_io.hpp"
#include "trajectory_recorder.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <memory>
#include <iomanip>
#include <vector>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <Windows.h>
#include <Psapi.h>
#endif

namespace fs = std::filesystem;

namespace {

// One job corresponds to one rock / start point / orientation combination.
struct SimulationJob
{
    std::size_t runIndex{};
    std::size_t rockIndex{};
    std::size_t startIndex{};
    std::uint32_t orientationIndex{};
};

std::string output_format_name(OutputFormat format)
{
    switch (format) {
    case OutputFormat::LegacyText: return "legacy-text";
    case OutputFormat::BinarySharded: return "binary-sharded";
    }
    return "unknown";
}

std::string output_format_name(const OutputSettings& settings)
{
    if (!settings.enabled) {
        return "none";
    }
    return output_format_name(settings.format);
}

const char* run_status_name(bool failed) noexcept
{
    return failed ? "FAILED" : "COMPLETED";
}

// maps a caught C++ exception to the string we want to write into the simulation report as stop criterion
std::string stop_reason_from_exception(const std::exception& e)
{
    const std::string message = e.what();
    if (message.find("DEM nodata") != std::string::npos) {
        return "nodata_encountered";
    }
    if (message.find("outside the DEM extent") != std::string::npos) {
        return "out_of_bounds";
    }
    return "initialization_failed";
}

// converts arbitrary text into valid quoted CSV field
std::string csv_quoted(const std::string& value)
{
    std::string escaped;
    escaped.reserve(value.size() + 2);
    escaped.push_back('"');
    for (char c : value) {
        if (c == '"') {
            escaped.push_back('"');
        }
        escaped.push_back(c);
    }
    escaped.push_back('"');
    return escaped;
}

std::string csv_optional_real(const std::optional<Real>& value)
{
    if (!value.has_value()) {
        return std::string();
    }
    std::ostringstream out;
    out << *value;
    return out.str();
}

double bytes_to_gib(std::uint64_t bytes) noexcept
{
    return static_cast<double>(bytes) / (1024.0 * 1024.0 * 1024.0);
}

// returns number of available CPU threads on the machine
std::uint32_t available_parallel_threads() noexcept
{
#ifdef _WIN32
    // sum of logical processors across all processor groups (for machines with more than 64 logical CPUs)
    const WORD groupCount = GetActiveProcessorGroupCount();
    std::uint32_t total = 0;
    for (WORD group = 0; group < groupCount; ++group) {
        total += static_cast<std::uint32_t>(GetActiveProcessorCount(group));
    }
    return std::max<std::uint32_t>(1u, total);
#else
    return std::max(1u, std::thread::hardware_concurrency());
#endif
}

struct WorkerGroupBinding
{
#ifdef _WIN32
    WORD group{ 0 };
    KAFFINITY affinityMask{ 0 };
#endif
};

// maps a worker index to a specific logical CPU
#ifdef _WIN32
WorkerGroupBinding choose_worker_binding(std::uint32_t workerIndex) noexcept
{
    WorkerGroupBinding binding;
    const WORD groupCount = GetActiveProcessorGroupCount();
    if (groupCount == 0) {
        binding.affinityMask = 1;
        return binding;
    }

    std::uint32_t remaining = workerIndex;
    for (WORD group = 0; group < groupCount; ++group) {
        const DWORD logicalCount = GetActiveProcessorCount(group);
        if (logicalCount == 0) {
            continue;
        }
        if (remaining < logicalCount) {
            binding.group = group;
            binding.affinityMask = (KAFFINITY(1) << remaining);
            return binding;
        }
        remaining -= logicalCount;
    }

    // Fallback to the last non-empty group if workerIndex exceeds the total.
    for (WORD group = groupCount; group-- > 0;) {
        const DWORD logicalCount = GetActiveProcessorCount(group);
        if (logicalCount == 0) {
            continue;
        }
        binding.group = group;
        binding.affinityMask = (KAFFINITY(1) << (logicalCount - 1));
        return binding;
    }

    binding.affinityMask = 1;
    return binding;
}

void bind_thread_to_worker_group(std::thread& thread, std::uint32_t workerIndex)
{
    const WorkerGroupBinding binding = choose_worker_binding(workerIndex);
    GROUP_AFFINITY affinity{};
    affinity.Group = binding.group;
    affinity.Mask = binding.affinityMask;
    SetThreadGroupAffinity(thread.native_handle(), &affinity, nullptr);
}
#else
void bind_thread_to_worker_group(std::thread&, std::uint32_t) {}
#endif

// for tracking memory during simulation
std::uint64_t current_process_peak_working_set_bytes() noexcept
{
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS counters{};
    if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters))) {
        return static_cast<std::uint64_t>(counters.PeakWorkingSetSize);
    }
#endif
    return 0;
}

// Timing/performance/memory report written only when requested in <debug/>
void write_timing_report(const fs::path& path,
    std::uint32_t threadCount,
    std::size_t totalRuns,
    std::size_t completedRuns,
    std::size_t failedRuns,
    std::size_t skippedRuns,
    double setupSeconds,
    double simulationSeconds,
    double rasterSeconds,
    double manifestSeconds,
    double totalSeconds,
    std::uint64_t peakWorkingSetBytes,
    std::size_t rasterLayerCount,
    std::size_t rasterOrderContributionCount)
{
    std::ofstream out(path.string());
    if (!out) {
        throw std::runtime_error("cannot open timing report file: " + path.string());
    }

    out << std::fixed << std::setprecision(6);
    out << "threads=" << threadCount << '\n';
    out << "total_runs=" << totalRuns << '\n';
    out << "completed_runs=" << completedRuns << '\n';
    out << "failed_runs=" << failedRuns << '\n';
    out << "skipped_runs=" << skippedRuns << '\n';
    out << "setup_seconds=" << setupSeconds << '\n';
    out << "simulation_seconds=" << simulationSeconds << '\n';
    out << "raster_seconds=" << rasterSeconds << '\n';
    out << "manifest_seconds=" << manifestSeconds << '\n';
    out << "total_seconds=" << totalSeconds << '\n';
    out << "peak_working_set_bytes=" << peakWorkingSetBytes << '\n';
    out << "peak_working_set_gib=" << bytes_to_gib(peakWorkingSetBytes) << '\n';
    out << "raster_layers=" << rasterLayerCount << '\n';
    out << "raster_order_contributions=" << rasterOrderContributionCount << '\n';
    if (totalSeconds > 0.0) {
        out << "throughput_trajectories_per_second="
            << (static_cast<double>(completedRuns) / totalSeconds) << '\n';
    }
}

// turns a file name into a safe output name by replacing unsupported characters with underscore
std::string sanitize_stem(const fs::path& path)
{
    std::string stem = path.stem().string();
    for (char& c : stem) {
        const bool ok = (c >= 'a' && c <= 'z')
            || (c >= 'A' && c <= 'Z')
            || (c >= '0' && c <= '9')
            || c == '_' || c == '-';
        if (!ok) {
            c = '_';
        }
    }
    return stem.empty() ? std::string("rock") : stem; //fall back to rock if result is empty
}

// generates unique file name for trajectory
std::string make_case_name(const fs::path& rockPath,
    std::size_t rockIndex,
    std::size_t startIndex,
    std::size_t orientationIndex)
{
    std::ostringstream out;
    out << sanitize_stem(rockPath)
        << "_r" << std::setw(3) << std::setfill('0') << rockIndex
        << "_s" << std::setw(4) << std::setfill('0') << startIndex
        << "_o" << std::setw(4) << std::setfill('0') << orientationIndex;
    return out.str();
}

// Small deterministic mixer used for reproducible per-trajectory seeds
std::uint64_t splitmix64(std::uint64_t x) noexcept
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

struct InitialPoseResult
{
    Eigen::Matrix<Real, 7, 1> pose;
    std::optional<Real> initialMinGap;
};

class CompositeTrajectoryRecorder final : public TrajectoryRecorder
{
public:
    CompositeTrajectoryRecorder(TrajectoryRecorder& primary,
        TrajectoryRecorder& secondary)
        : primary_(primary), secondary_(secondary)
    {
    }

    void beginTrajectory(const TrajectoryMetadata& meta,
        const Rock& rock,
        const RecordedStepData& initialStep) override
    {
        primary_.beginTrajectory(meta, rock, initialStep);
        secondary_.beginTrajectory(meta, rock, initialStep);
    }

    void onAcceptedStep(const RecordedStepData& step) override
    {
        primary_.onAcceptedStep(step);
        secondary_.onAcceptedStep(step);
    }

    RecordedTrajectoryLocation endTrajectory(const RecordedStepData& finalStep,
        double simSeconds) override
    {
        const RecordedTrajectoryLocation location = primary_.endTrajectory(finalStep, simSeconds);
        secondary_.endTrajectory(finalStep, simSeconds);
        return location;
    }

private:
    TrajectoryRecorder& primary_;
    TrajectoryRecorder& secondary_;
};

// Build the initial pose either from an absolute XYZ release point or from an
// X/Y location plus a dz clearance above the terrain after orientation.
InitialPoseResult make_initial_pose(const Terrain& dem,
    const Rock& rock,
    const StartPoint& startPoint,
    std::uint64_t orientationSeed,
    std::size_t globalRunIndex,
    TerrainQueryWorkspace& terrainWorkspace)
{
    // pseudo-random initial orientation. Reproducible with same seed.
    std::mt19937_64 rng(splitmix64(orientationSeed ^ (globalRunIndex + 1)));
    const auto q = make_random_quaternion_wxyz(rng);

    Quat quatK(q(0), q(1), q(2), q(3));
    quatK.normalize();
    const Mat3 AIK = quatK.toRotationMatrix();

    // Calculate the terrain gap for z = 0 so we can either convert a dz
    // clearance into an absolute release height or derive the actual initial
    // jump height for an absolute-Z release point.
    const Vec3 startXY0(startPoint.x, startPoint.y, Real(0));
    const TerrainQueryResult query = dem.queryPose(rock, AIK, startXY0, terrainWorkspace);
    if (query.state == TerrainContainment::OutOfBounds) {
        throw std::runtime_error("initial rock placement extends outside the DEM extent at start point");
    }
    if (query.state == TerrainContainment::NoData) {
        throw std::runtime_error("initial rock placement encounters DEM nodata at start point");
    }

    const Real z0 = startPoint.dz.has_value()
        ? (-query.minGap + *startPoint.dz)
        : startPoint.z;
    const Real initialMinGap = query.minGap + z0;

    InitialPoseResult result;
    result.pose << startPoint.x, startPoint.y, z0, q(0), q(1), q(2), q(3);
    result.initialMinGap = initialMinGap;
    return result;
}

// Translate XML configuration into the runtime timestepper object.
TimeStepper make_time_stepper(const BatchSimulationConfig& cfg)
{
    TimeStepper step;
    step.solverType = cfg.solver.type;
    step.alpha = cfg.solver.alpha;
    step.atol = cfg.solver.atol;
    step.rtol = cfg.solver.rtol;
    step.maxiter = cfg.solver.maxiter;
    return step;
}

// Cartesian product of rocks, release points, and sampled orientations.
std::vector<SimulationJob> build_jobs(const BatchSimulationConfig& cfg)
{
    std::vector<SimulationJob> jobs;
    jobs.reserve(cfg.rockPaths.size()
        * cfg.startPoints.size()
        * static_cast<std::size_t>(cfg.randomOrientations));

    std::size_t runIndex = 0;
    for (std::size_t rockIndex = 0; rockIndex < cfg.rockPaths.size(); ++rockIndex) {
        for (std::size_t startIndex = 0; startIndex < cfg.startPoints.size(); ++startIndex) {
            for (std::uint32_t orientationIndex = 0; orientationIndex < cfg.randomOrientations; ++orientationIndex) {
                jobs.push_back(SimulationJob{
                    runIndex,
                    rockIndex,
                    startIndex,
                    orientationIndex
                    });
                ++runIndex;
            }
        }
    }
    return jobs;
}

int run_batch_threaded(const BatchSimulationConfig& cfg,
    std::uint32_t threadCount,
    std::uint64_t orientationSeed,
    std::uint64_t roughnessSeed)
{
    // Setup phase: load terrain, configure roughness, build reusable rock templates, and assemble the full batch of jobs once.
    const auto totalPhaseStart = std::chrono::high_resolution_clock::now();
    const auto setupPhaseStart = totalPhaseStart;
    fs::create_directories(cfg.outputDir);

    if (cfg.output.enabled && cfg.output.format == OutputFormat::LegacyText && !cfg.output.fields.empty()) {
        throw std::runtime_error("variable-selectable output is currently implemented only for format=\"binary-sharded\"");
    }
    if (cfg.rasterOutput.has_value() && cfg.output.enabled && cfg.output.format == OutputFormat::LegacyText) {
        throw std::runtime_error("raster output currently supports raster-only runs or raster combined with format=\"binary-sharded\"");
    }

    Terrain dem = load_terrain(cfg.demPath.string());
    if (cfg.roughnessRasterPath) {
        dem.addRoughness(
            cfg.roughnessRasterPath->string(),
            cfg.roughnessStepFunc,
            roughnessSeed);
    }

    // export terrain with roughness in a specific resolution if defined in <debug/>
    if (cfg.terrainExportPath && cfg.terrainExportResolution) {
        fs::create_directories(cfg.terrainExportPath->parent_path());
        export_sampled_terrain_asc(dem, cfg.terrainExportPath->string(), *cfg.terrainExportResolution);
        std::cout << "Exported terrain raster to " << cfg.terrainExportPath->string() << "\n";
    }

    // read rock PTS files
    std::vector<Rock> rockTemplates;
    std::vector<fs::path> rockTemplatePaths;
    rockTemplates.reserve(cfg.rockPaths.size());
    rockTemplatePaths.reserve(cfg.rockPaths.size());
    for (const auto& rockPath : cfg.rockPaths) {
        MatX verts = load_xyz(rockPath.string());
        Eigen::Matrix<Real, 7, 1> identityPose;
        identityPose << 0, 0, 0, 1, 0, 0, 0;
        rockTemplates.emplace_back(verts, cfg.rockDensity, Vec6::Zero(), identityPose);
        rockTemplatePaths.push_back(rockPath);
    }

    const std::vector<SimulationJob> jobs = build_jobs(cfg);
    const std::size_t totalRuns = jobs.size();
    const std::optional<RasterGridGeometry> rasterGrid = cfg.rasterOutput.has_value()
        ? std::optional<RasterGridGeometry>(make_raster_grid_geometry(dem, cfg.rasterOutput->resolution))
        : std::nullopt;
    std::vector<std::string> manifestLines(totalRuns);
    std::vector<MemoryRasterWorkerResult> rasterBuffers(threadCount);
    std::atomic<std::size_t> nextJob{ 0 };
    std::atomic<std::size_t> finishedJobs{ 0 };
    std::atomic<std::size_t> completedJobs{ 0 };
    std::atomic<std::size_t> failedJobs{ 0 };
    std::atomic<std::size_t> skippedJobs{ 0 };
    std::atomic<bool> progressDone{ false };
    const Real dumpStep = cfg.dumpStep.value_or(cfg.dt);
    std::mutex coutMutex;
    const auto setupPhaseEnd = std::chrono::high_resolution_clock::now();

    // lightweight progress reporter keeps console output human-readable without adding per-trajectory synchronization noise
    auto progressReporter = [&]() {
        using namespace std::chrono_literals;
        while (!progressDone.load()) {
            const std::size_t finished = finishedJobs.load();
            const double percent = totalRuns > 0
                ? (100.0 * static_cast<double>(finished) / static_cast<double>(totalRuns))
                : 100.0;
            {
                std::lock_guard<std::mutex> lock(coutMutex);
                std::cout << '\r'
                    << "Progress: " << finished << " / " << totalRuns
                    << " trajectories (" << std::fixed << std::setprecision(1) << percent << "%)"
                    << std::flush;
            }
            if (finished >= totalRuns) {
                break;
            }
            // only report once every second
            std::this_thread::sleep_for(1s);
        }

        const std::size_t finished = finishedJobs.load();
        const double percent = totalRuns > 0
            ? (100.0 * static_cast<double>(finished) / static_cast<double>(totalRuns))
            : 100.0;
        std::lock_guard<std::mutex> lock(coutMutex);
        std::cout << '\r'
            << "Progress: " << finished << " / " << totalRuns
            << " trajectories (" << std::fixed << std::setprecision(1) << percent << "%)"
            << '\n';
    };

    auto worker = [&](std::uint32_t threadIndex) {
        // One workspace and one optional recorder shard per worker thread.
        TimeStepper step = make_time_stepper(cfg);
        TimeStepperWorkspace workspace;
        std::unique_ptr<BinaryShardRecorder> binaryRecorder;
        if (cfg.output.enabled && cfg.output.format == OutputFormat::BinarySharded) {
            binaryRecorder = std::make_unique<BinaryShardRecorder>(
                cfg.outputDir,
                threadIndex,
                dumpStep,
                cfg.output);
        }
        std::unique_ptr<MemoryRasterRecorder> rasterRecorder;
        if (cfg.rasterOutput.has_value()) {
            rasterRecorder = std::make_unique<MemoryRasterRecorder>(*rasterGrid, *cfg.rasterOutput, dumpStep);
        }

        while (true) {
            // automatically pull available job
            const std::size_t jobIdx = nextJob.fetch_add(1);
            if (jobIdx >= jobs.size()) {
                break;
            }

            const SimulationJob& job = jobs[jobIdx];
            const Rock& rockTemplate = rockTemplates[job.rockIndex];
            const fs::path& rockPath = rockTemplatePaths[job.rockIndex];
            const StartPoint& startPoint = cfg.startPoints[job.startIndex];
            std::optional<Real> manifestStartZ = startPoint.dz.has_value()
                ? std::nullopt
                : std::optional<Real>(startPoint.z);

            const std::string caseName = make_case_name(
                rockPath,
                job.rockIndex + 1,
                job.startIndex + 1,
                static_cast<std::size_t>(job.orientationIndex + 1));
            try {
                // generate initial conditions
                const InitialPoseResult initialPlacement = make_initial_pose(
                    dem,
                    rockTemplate,
                    startPoint,
                    orientationSeed,
                    job.runIndex,
                    workspace.terrain);

                Rock runRock = rockTemplate.withInitialState(
                    Vec6::Zero(),
                    initialPlacement.pose);
                manifestStartZ = initialPlacement.pose(2);

                const auto t0 = std::chrono::high_resolution_clock::now();
                std::size_t writtenSamples = 0;
                std::string trajectoryStorage;
                std::uint64_t byteOffset = 0;
                TimeStepperRunSummary recordedSummary;
                bool hasRecordedSummary = false;
                TimeStepperStopReason finalStopReason = TimeStepperStopReason::Unknown;
                std::size_t acceptedStepCount = 0;
                std::string stopCriterion = "unknown";

                if (cfg.output.enabled && cfg.output.format == OutputFormat::BinarySharded) {
                    // Stream dumped samples directly into a per-thread binary
                    // shard instead of materializing the full trajectory.
                    const TrajectoryMetadata meta{
                        job.runIndex,
                        job.rockIndex,
                        job.startIndex,
                        job.orientationIndex,
                        caseName
                    };
                    const RecordedStepData initialStep{
                        Real(0),
                        runRock.q0(),
                        runRock.u0(),
                        initialPlacement.initialMinGap
                    };
                    if (cfg.rasterOutput.has_value()) {
                        CompositeTrajectoryRecorder recorder(*binaryRecorder, *rasterRecorder);
                        recorder.beginTrajectory(meta, runRock, initialStep);
                        recordedSummary = step.runRecorded(
                            cfg.mu,
                            dem,
                            runRock,
                            cfg.dt,
                            workspace,
                            recorder,
                            cfg.stop);
                    }
                    else {
                        binaryRecorder->beginTrajectory(meta, runRock, initialStep);
                        recordedSummary = step.runRecorded(
                            cfg.mu,
                            dem,
                            runRock,
                            cfg.dt,
                            workspace,
                            *binaryRecorder,
                            cfg.stop);
                    }
                    hasRecordedSummary = true;
                    finalStopReason = recordedSummary.stopReason;
                    acceptedStepCount = recordedSummary.acceptedSteps;
                    stopCriterion = stop_reason_name(recordedSummary.stopReason);
                }
                else if (cfg.rasterOutput.has_value()) {
                    const TrajectoryMetadata meta{
                        job.runIndex,
                        job.rockIndex,
                        job.startIndex,
                        job.orientationIndex,
                        caseName
                    };
                    rasterRecorder->beginTrajectory(meta, runRock, RecordedStepData{
                        Real(0),
                        runRock.q0(),
                        runRock.u0(),
                        initialPlacement.initialMinGap
                        });

                    recordedSummary = step.runRecorded(
                        cfg.mu,
                        dem,
                        runRock,
                        cfg.dt,
                        workspace,
                        *rasterRecorder,
                        cfg.stop);
                    hasRecordedSummary = true;
                    finalStopReason = recordedSummary.stopReason;
                    acceptedStepCount = recordedSummary.acceptedSteps;
                    stopCriterion = stop_reason_name(recordedSummary.stopReason);
                }
                else {
                    // Legacy text output still stores the full accepted
                    // trajectory in memory first, then writes it afterward.
                    const fs::path trajectoryPath = cfg.outputDir / (caseName + ".txt");
                    TimeStepperResult result = step.run(cfg.mu, dem, runRock, cfg.dt, workspace, cfg.stop);
                    finalStopReason = result.stopReason;
                    acceptedStepCount = result.t.empty() ? 0 : (result.t.size() - 1);
                    writtenSamples = write_trajectory_txt(
                        trajectoryPath.string(),
                        result,
                        6,
                        dumpStep);
                    trajectoryStorage = trajectoryPath.filename().string();
                    stopCriterion = stop_reason_name(result.stopReason);
                }
                const auto t1 = std::chrono::high_resolution_clock::now();
                const double simSeconds = std::chrono::duration<double>(t1 - t0).count();

                if (hasRecordedSummary) {
                    const RecordedStepData finalStep{
                        recordedSummary.tFinal,
                        recordedSummary.qFinal,
                        recordedSummary.uFinal,
                        recordedSummary.minGapFinal
                    };
                    if (cfg.output.enabled && cfg.output.format == OutputFormat::BinarySharded) {
                        if (cfg.rasterOutput.has_value()) {
                            CompositeTrajectoryRecorder recorder(*binaryRecorder, *rasterRecorder);
                            const RecordedTrajectoryLocation location = recorder.endTrajectory(finalStep, simSeconds);
                            writtenSamples = static_cast<std::size_t>(location.sampleCount);
                            trajectoryStorage = location.storageFile;
                            byteOffset = location.byteOffset;
                        }
                        else {
                            const RecordedTrajectoryLocation location = binaryRecorder->endTrajectory(finalStep, simSeconds);
                            writtenSamples = static_cast<std::size_t>(location.sampleCount);
                            trajectoryStorage = location.storageFile;
                            byteOffset = location.byteOffset;
                        }
                    }
                    else if (cfg.rasterOutput.has_value()) {
                        rasterRecorder->endTrajectory(finalStep, simSeconds);
                    }
                }

                const bool failed = (finalStopReason == TimeStepperStopReason::SolverNonConvergence);

                // write trajectory report
                std::ostringstream line;
                line << (job.runIndex + 1) << ','
                    << (job.rockIndex + 1) << ','
                    << csv_quoted(rockPath.string()) << ','
                    << (job.startIndex + 1) << ','
                    << startPoint.x << ',' << startPoint.y << ',' << csv_optional_real(manifestStartZ) << ','
                    << (job.orientationIndex + 1) << ','
                    << csv_quoted(caseName) << ','
                    << run_status_name(failed) << ','
                    << csv_quoted(stopCriterion) << ','
                    << csv_quoted(output_format_name(cfg.output)) << ','
                    << csv_quoted(trajectoryStorage) << ','
                    << byteOffset << ','
                    << simSeconds << ','
                    << writtenSamples;
                manifestLines[job.runIndex] = line.str();
                finishedJobs.fetch_add(1);
                if (failed) {
                    failedJobs.fetch_add(1);
                }
                else {
                    completedJobs.fetch_add(1);
                }

                if (failed) {
                    std::lock_guard<std::mutex> lock(coutMutex);
                    std::cout << "Warning: solver did not converge for " << caseName
                        << " on thread " << (threadIndex + 1)
                        << "; trajectory marked FAILED after "
                        << acceptedStepCount << " accepted steps.\n";
                }
            }
            catch (const std::exception& e) {
                if (binaryRecorder) {
                    binaryRecorder->abortActiveTrajectory();
                }
                if (rasterRecorder) {
                    rasterRecorder->abortActiveTrajectory();
                }
                std::ostringstream line;
                line << (job.runIndex + 1) << ','
                    << (job.rockIndex + 1) << ','
                    << csv_quoted(rockPath.string()) << ','
                    << (job.startIndex + 1) << ','
                    << startPoint.x << ',' << startPoint.y << ',' << csv_optional_real(manifestStartZ) << ','
                    << (job.orientationIndex + 1) << ','
                    << csv_quoted(caseName) << ','
                    << "SKIPPED" << ','
                    << csv_quoted(stop_reason_from_exception(e)) << ','
                    << csv_quoted(output_format_name(cfg.output)) << ','
                    << csv_quoted("") << ','
                    << 0 << ','
                    << 0 << ','
                    << 0;
                manifestLines[job.runIndex] = line.str();
                finishedJobs.fetch_add(1);
                skippedJobs.fetch_add(1);

                std::lock_guard<std::mutex> lock(coutMutex);
                std::cout << "Warning: skipped " << caseName
                    << " on thread " << (threadIndex + 1)
                    << " because initialization failed: " << e.what() << '\n';
            }
            catch (...) {
                if (binaryRecorder) {
                    binaryRecorder->abortActiveTrajectory();
                }
                if (rasterRecorder) {
                    rasterRecorder->abortActiveTrajectory();
                }
                std::ostringstream line;
                line << (job.runIndex + 1) << ','
                    << (job.rockIndex + 1) << ','
                    << csv_quoted(rockPath.string()) << ','
                    << (job.startIndex + 1) << ','
                    << startPoint.x << ',' << startPoint.y << ',' << csv_optional_real(manifestStartZ) << ','
                    << (job.orientationIndex + 1) << ','
                    << csv_quoted(caseName) << ','
                    << "SKIPPED" << ','
                    << csv_quoted("unknown") << ','
                    << csv_quoted(output_format_name(cfg.output)) << ','
                    << csv_quoted("") << ','
                    << 0 << ','
                    << 0 << ','
                    << 0;
                manifestLines[job.runIndex] = line.str();
                finishedJobs.fetch_add(1);
                skippedJobs.fetch_add(1);

                std::lock_guard<std::mutex> lock(coutMutex);
                std::cout << "Warning: skipped " << caseName
                    << " on thread " << (threadIndex + 1)
                    << " because initialization failed with an unknown exception.\n";
            }
        }

        if (rasterRecorder) {
            rasterBuffers[threadIndex] = rasterRecorder->releaseResult();
        }
    };

    const auto simulationPhaseStart = std::chrono::high_resolution_clock::now();
    std::thread progressThread(progressReporter);

    // launch all worker threads
    std::vector<std::thread> threads;
    threads.reserve(threadCount);
    for (std::uint32_t i = 0; i < threadCount; ++i) {
        threads.emplace_back(worker, i);
        bind_thread_to_worker_group(threads.back(), i);
    }
    for (auto& thread : threads) {
        thread.join();
    }
    progressDone.store(true);
    if (progressThread.joinable()) {
        progressThread.join();
    }
    const auto simulationPhaseEnd = std::chrono::high_resolution_clock::now();
    const auto rasterPhaseStart = simulationPhaseEnd;

    RasterBuildSummary rasterSummary;
    if (cfg.rasterOutput.has_value()) {
        rasterSummary = build_rasters_from_workers(
            cfg.outputDir,
            *cfg.rasterOutput,
            *rasterGrid,
            std::move(rasterBuffers));
    }
    const auto rasterPhaseEnd = std::chrono::high_resolution_clock::now();

    // write simulation report
    const auto manifestPhaseStart = rasterPhaseEnd;
    std::ofstream manifest((cfg.outputDir / "manifest.csv").string());
    if (!manifest) {
        throw std::runtime_error("cannot open manifest.csv in output directory");
    }
    manifest << "run_index,rock_index,rock_path,start_index,start_x,start_y,start_z,orientation_index,case_name,status,stop_criterion,output_format,trajectory_file,byte_offset,sim_seconds,trajectory_samples\n";
    for (const auto& line : manifestLines) {
        manifest << line << '\n';
    }
    const auto manifestPhaseEnd = std::chrono::high_resolution_clock::now();

    if (cfg.writeTimingReport) {
        const auto totalPhaseEnd = std::chrono::high_resolution_clock::now();
        const double setupSeconds = std::chrono::duration<double>(setupPhaseEnd - setupPhaseStart).count();
        const double simulationSeconds = std::chrono::duration<double>(simulationPhaseEnd - simulationPhaseStart).count();
        const double rasterSeconds = std::chrono::duration<double>(rasterPhaseEnd - rasterPhaseStart).count();
        const double manifestSeconds = std::chrono::duration<double>(manifestPhaseEnd - manifestPhaseStart).count();
        const double totalSeconds = std::chrono::duration<double>(totalPhaseEnd - totalPhaseStart).count();
        write_timing_report(
            cfg.outputDir / "timing_report.txt",
            threadCount,
            totalRuns,
            completedJobs.load(),
            failedJobs.load(),
            skippedJobs.load(),
            setupSeconds,
            simulationSeconds,
            rasterSeconds,
            manifestSeconds,
            totalSeconds,
            current_process_peak_working_set_bytes(),
            rasterSummary.layerCount,
            rasterSummary.orderContributionCount);
    }

    std::cout << "\nCompleted " << totalRuns << " simulations using "
        << threadCount << " thread(s). Output written to "
        << cfg.outputDir.string() << "\n";
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "usage: rockfall3d config.xml\n";
        return 1;
    }

    try {
        const BatchSimulationConfig cfg = load_batch_simulation_config(argv[1]);
        const std::uint64_t orientationSeed = cfg.orientationSeed.value_or(
            splitmix64(static_cast<std::uint64_t>(std::random_device{}())));
        const std::uint64_t roughnessSeed = cfg.roughnessSeed.value_or(
            splitmix64(static_cast<std::uint64_t>(std::random_device{}()) ^ 0xa5a5a5a5ULL));

        const std::uint32_t availableCores = available_parallel_threads();
        const std::uint32_t requestedThreads = std::max(1u, cfg.parallelThreads);
        const std::uint32_t threadCount = std::min(requestedThreads, availableCores);

        if (requestedThreads > availableCores) {
            std::cout << "Requested " << requestedThreads << " parallel threads, limiting to "
                << availableCores << " available logical processors.\n";
        }

        return run_batch_threaded(cfg, threadCount, orientationSeed, roughnessSeed);
    }
    catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
