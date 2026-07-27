#pragma once

#include "terrain.hpp"
#include "trajectory_recorder.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef WITH_GDAL
#include "gdal.h"
#endif

namespace fs = std::filesystem;

// Possible metrics for raster output
enum class RasterMetric : std::uint8_t
{
    None = 0,
    X,
    Y,
    Z,
    Vx,
    Vy,
    Vz,
    Vtot,
    Wx,
    Wy,
    Wz,
    Wtot,
    EkinTrans,
    EkinRot,
    EkinTot,
    JumpH
};

// Possible statistics to compute for each raster cell
enum class RasterStatistic : std::uint8_t
{
    Count = 0,
    Min,
    Max,
    Sum,
    Mean,
    Median,
    Percentile
};

// Order-based statistics require storing all values for a cell, while streaming statistics can be computed in one pass without storing all values.
inline bool raster_statistic_is_order_based(RasterStatistic statistic) noexcept
{
    return statistic == RasterStatistic::Median
        || statistic == RasterStatistic::Percentile;
}

// Order-based statistics require storing all values for a cell, while streaming statistics can be computed in one pass without storing all values.
inline bool raster_statistic_is_streaming(RasterStatistic statistic) noexcept
{
    return !raster_statistic_is_order_based(statistic);
}

inline std::string raster_metric_name(RasterMetric metric)
{
    switch (metric) {
    case RasterMetric::X: return "x";
    case RasterMetric::Y: return "y";
    case RasterMetric::Z: return "z";
    case RasterMetric::Vx: return "vx";
    case RasterMetric::Vy: return "vy";
    case RasterMetric::Vz: return "vz";
    case RasterMetric::Vtot: return "vtot";
    case RasterMetric::Wx: return "wx";
    case RasterMetric::Wy: return "wy";
    case RasterMetric::Wz: return "wz";
    case RasterMetric::Wtot: return "wtot";
    case RasterMetric::EkinTrans: return "ekin_trans";
    case RasterMetric::EkinRot: return "ekin_rot";
    case RasterMetric::EkinTot: return "ekin_tot";
    case RasterMetric::JumpH: return "jump_h";
    default: return "rocks";
    }
}

inline std::string raster_statistic_name(RasterStatistic statistic)
{
    switch (statistic) {
    case RasterStatistic::Count: return "count";
    case RasterStatistic::Min: return "min";
    case RasterStatistic::Max: return "max";
    case RasterStatistic::Sum: return "sum";
    case RasterStatistic::Mean: return "mean";
    case RasterStatistic::Median: return "median";
    case RasterStatistic::Percentile: return "percentile";
    }
    return "unknown";
}

struct RasterLayerSpec
{
    RasterMetric metric{ RasterMetric::None };
    RasterStatistic statistic{ RasterStatistic::Max };
    std::optional<Real> percentile;
    std::string outputStem;
};

struct RasterOutputConfig
{
    Real resolution{ Real(2) };
    std::vector<RasterLayerSpec> layers;
};

struct RasterGridGeometry
{
    Real xll{ Real(0) };
    Real yll{ Real(0) };
    Real cell{ Real(2) };
    std::uint32_t ncols{ 0 };
    std::uint32_t nrows{ 0 };
    std::string projectionWkt;
};

// snaps raster origin down to a multiple of the resolution, so that the output grid is nicely aligned and not offset by a small fraction of the resolution
inline Real snap_origin_down(Real value, Real resolution) noexcept
{
    return std::floor(value / resolution) * resolution;
}

// Compute the geometry of the output raster grid based on the DEM geometry and the requested output resolution. 
// The grid is snapped down to a multiple of the resolution to ensure nice alignment. 
// The grid covers at least the DEM extent, but may be larger due to snapping and rounding.
inline RasterGridGeometry make_raster_grid_geometry(const Terrain& dem, Real resolution)
{
    if (!(resolution > Real(0))) {
        throw std::runtime_error("raster output resolution must be positive");
    }

    RasterGridGeometry grid;
    grid.xll = snap_origin_down(dem.xll, resolution);
    grid.yll = snap_origin_down(dem.yll, resolution);
    grid.cell = resolution;
    grid.projectionWkt = dem.projectionWkt;

    const Real width = dem.xmax - grid.xll;
    const Real height = dem.ymax - grid.yll;
    grid.ncols = std::max<std::uint32_t>(
        1u,
        static_cast<std::uint32_t>(std::ceil(width / resolution)));
    grid.nrows = std::max<std::uint32_t>(
        1u,
        static_cast<std::uint32_t>(std::ceil(height / resolution)));
    if (grid.ncols == 0 || grid.nrows == 0) {
        throw std::runtime_error("raster output grid is empty");
    }
    return grid;
}

// Resolve a world-space query point to a raster cell index. 
// Returns nullopt if the point is outside the grid.
inline std::optional<std::uint32_t> raster_cell_index(const RasterGridGeometry& grid,
    Real x,
    Real y) noexcept
{
    if (x < grid.xll || y < grid.yll) {
        return std::nullopt;
    }

    const Real colf = std::floor((x - grid.xll) / grid.cell);
    const Real rowf = std::floor((y - grid.yll) / grid.cell);
    if (colf < Real(0) || rowf < Real(0)) {
        return std::nullopt;
    }

    std::uint32_t col = static_cast<std::uint32_t>(colf);
    std::uint32_t row = static_cast<std::uint32_t>(rowf);
    if (col >= grid.ncols) {
        if (colf == Real(grid.ncols)) {
            col = grid.ncols - 1;
        }
        else {
            return std::nullopt;
        }
    }
    if (row >= grid.nrows) {
        if (rowf == Real(grid.nrows)) {
            row = grid.nrows - 1;
        }
        else {
            return std::nullopt;
        }
    }

    return row * grid.ncols + col;
}

// Generate a default output stem for a raster layer based on its metric and statistic.
inline std::string default_raster_output_stem(const RasterLayerSpec& layer)
{
    if (!layer.outputStem.empty()) {
        return layer.outputStem;
    }
    if (layer.statistic == RasterStatistic::Percentile) {
        const double p = layer.percentile.has_value() ? static_cast<double>(*layer.percentile) : 0.0;
        const long long rounded = static_cast<long long>(std::llround(p));
        return "q" + std::to_string(rounded) + "_" + raster_metric_name(layer.metric);
    }
    if (layer.statistic == RasterStatistic::Count && layer.metric == RasterMetric::None) {
        return "count_rocks";
    }
    return raster_statistic_name(layer.statistic) + "_" + raster_metric_name(layer.metric);
}

// Sample the value of the specified metric from the given step data. 
// Returns nullopt if the metric is not available.
inline std::optional<Real> sample_metric_value(RasterMetric metric,
    const RecordedStepData& step)
{
    switch (metric) {
    case RasterMetric::X: return step.q(0);
    case RasterMetric::Y: return step.q(1);
    case RasterMetric::Z: return step.q(2);
    case RasterMetric::Vx: return step.u(0);
    case RasterMetric::Vy: return step.u(1);
    case RasterMetric::Vz: return step.u(2);
    case RasterMetric::Vtot: return step.u.head<3>().norm();
    case RasterMetric::Wx: return step.u(3);
    case RasterMetric::Wy: return step.u(4);
    case RasterMetric::Wz: return step.u(5);
    case RasterMetric::Wtot: return step.u.tail<3>().norm();
    case RasterMetric::EkinTrans:
        return Real(0.5) * step.mass * step.u.head<3>().squaredNorm();
    case RasterMetric::EkinRot:
        return Real(0.5) * step.u.tail<3>().dot(step.inertiaK * step.u.tail<3>());
    case RasterMetric::EkinTot:
        return Real(0.5) * step.mass * step.u.head<3>().squaredNorm()
            + Real(0.5) * step.u.tail<3>().dot(step.inertiaK * step.u.tail<3>());
    case RasterMetric::JumpH:
        if (step.minGap.has_value()) {
            return *step.minGap;
        }
        return std::nullopt;
    default:
        return std::nullopt;
    }
}

// Contribution of a single trajectory step to an order-based statistic for a specific raster cell.
struct OrderContribution
{
    std::uint32_t cellIndex{ 0 };
    float value{ 0.0f };
};

// Accumulator for streaming statistics for a single raster layer.
class StreamingLayerAccumulator
{
public:
    StreamingLayerAccumulator() = default;
	// Initialize the accumulator for the specified statistic and number of cells.
    StreamingLayerAccumulator(RasterStatistic statistic,
        std::size_t cellCount)
        : statistic_(statistic),
        count_(cellCount, 0u)
    {
        switch (statistic_) {
        case RasterStatistic::Min:
            min_.assign(cellCount, std::numeric_limits<float>::infinity());
            break;
        case RasterStatistic::Max:
            max_.assign(cellCount, -std::numeric_limits<float>::infinity());
            break;
        case RasterStatistic::Sum:
        case RasterStatistic::Mean:
            sum_.assign(cellCount, 0.0);
            break;
        case RasterStatistic::Count:
            break;
        default:
            throw std::runtime_error("streaming layer accumulator received non-streaming statistic");
        }
    }
	// Accumulate a count contribution for the specified cell index.
    void accumulate_count(std::uint32_t cellIndex)
    {
        ++count_[cellIndex];
    }

	// Accumulate a value contribution for the specified cell index, updating the count and the relevant statistic.
    void accumulate_value(std::uint32_t cellIndex, float value)
    {
        ++count_[cellIndex];
        switch (statistic_) {
        case RasterStatistic::Min:
            min_[cellIndex] = std::min(min_[cellIndex], value);
            break;
        case RasterStatistic::Max:
            max_[cellIndex] = std::max(max_[cellIndex], value);
            break;
        case RasterStatistic::Sum:
        case RasterStatistic::Mean:
            sum_[cellIndex] += value;
            break;
        case RasterStatistic::Count:
            break;
        default:
            break;
        }
    }

	// Merge another accumulator into this one, combining their counts and statistics.
    void merge_from(const StreamingLayerAccumulator& other)
    {
        if (statistic_ != other.statistic_ || count_.size() != other.count_.size()) {
            throw std::runtime_error("cannot merge incompatible streaming layer accumulators");
        }

        for (std::size_t i = 0; i < count_.size(); ++i) {
            count_[i] += other.count_[i];
        }

        switch (statistic_) {
        case RasterStatistic::Min:
            for (std::size_t i = 0; i < min_.size(); ++i) {
                min_[i] = std::min(min_[i], other.min_[i]);
            }
            break;
        case RasterStatistic::Max:
            for (std::size_t i = 0; i < max_.size(); ++i) {
                max_[i] = std::max(max_[i], other.max_[i]);
            }
            break;
        case RasterStatistic::Sum:
        case RasterStatistic::Mean:
            for (std::size_t i = 0; i < sum_.size(); ++i) {
                sum_[i] += other.sum_[i];
            }
            break;
        case RasterStatistic::Count:
            break;
        default:
            break;
        }
    }

	// Finalize the accumulated statistics into a raster layer, applying the specified nodata value for cells with no contributions.
    std::vector<float> finalize(float nodataValue) const
    {
        std::vector<float> raster(count_.size(), nodataValue);
        for (std::size_t i = 0; i < count_.size(); ++i) {
            if (count_[i] == 0u) {
                continue;
            }
            switch (statistic_) {
            case RasterStatistic::Count:
                raster[i] = static_cast<float>(count_[i]);
                break;
            case RasterStatistic::Min:
                raster[i] = min_[i];
                break;
            case RasterStatistic::Max:
                raster[i] = max_[i];
                break;
            case RasterStatistic::Sum:
                raster[i] = static_cast<float>(sum_[i]);
                break;
            case RasterStatistic::Mean:
                raster[i] = static_cast<float>(sum_[i] / static_cast<double>(count_[i]));
                break;
            default:
                break;
            }
        }
        return raster;
    }

private:
    RasterStatistic statistic_{ RasterStatistic::Count };
    std::vector<std::uint32_t> count_;
    std::vector<double> sum_;
    std::vector<float> min_;
    std::vector<float> max_;
};

struct MemoryRasterWorkerResult
{
    std::vector<StreamingLayerAccumulator> streamingLayers;
    std::vector<std::vector<OrderContribution>> orderLayers;
    std::size_t orderContributionCount{ 0 };
};

struct RasterBuildSummary
{
    std::size_t layerCount{ 0 };
    std::size_t orderContributionCount{ 0 };
};

struct RasterCropBounds
{
    std::uint32_t colMin{ 0 };
    std::uint32_t colMax{ 0 };
    std::uint32_t rowMin{ 0 };
    std::uint32_t rowMax{ 0 };
    bool hasData{ false };
};

// Trajectory recorder that accumulates raster contributions in memory for a single trajectory.
// Produces a result that can be merged with other trajectories' results to build the final raster layers.
class MemoryRasterRecorder final : public TrajectoryRecorder
{
public:
    MemoryRasterRecorder(const RasterGridGeometry& grid,
        const RasterOutputConfig& config,
        Real dumpStep)
        : grid_(grid),
        config_(config),
        dumpStep_(dumpStep),
        layerCount_(config.layers.size())
    {
        if (!(dumpStep_ > Real(0))) {
            throw std::runtime_error("memory raster output requires a positive dump step");
        }

        const std::size_t cellCount = static_cast<std::size_t>(grid_.ncols) * grid_.nrows;
        streamingSlotByLayer_.assign(layerCount_, -1);
        orderSlotByLayer_.assign(layerCount_, -1);
        metricSlotByLayer_.assign(layerCount_, -1);

		// Initialize accumulators for each layer and track unique metrics.
        for (std::size_t layerIdx = 0; layerIdx < layerCount_; ++layerIdx) {
            const auto& layer = config_.layers[layerIdx];
            if (raster_statistic_is_streaming(layer.statistic)) {
                streamingSlotByLayer_[layerIdx] = static_cast<int>(result_.streamingLayers.size());
                result_.streamingLayers.emplace_back(layer.statistic, cellCount);
            }
            else {
                orderSlotByLayer_[layerIdx] = static_cast<int>(result_.orderLayers.size());
                result_.orderLayers.emplace_back();
            }

            if (layer.statistic != RasterStatistic::Count) {
                int metricSlot = find_metric_slot(layer.metric);
                if (metricSlot < 0) {
                    metricSlot = static_cast<int>(uniqueMetrics_.size());
                    uniqueMetrics_.push_back(layer.metric);
                }
                metricSlotByLayer_[layerIdx] = metricSlot;
            }
        }

        metricValues_.resize(uniqueMetrics_.size());
    }

	// Begin a new trajectory, resetting the internal state and accumulating the initial step.
    void beginTrajectory(const TrajectoryMetadata&,
        const Rock& rock,
        const RecordedStepData& initialStep) override
    {
        if (active_) {
            throw std::runtime_error("memory raster recorder beginTrajectory called while another trajectory is still active");
        }

        (void)rock;
        active_ = true;
        trajectoryCellIndex_.clear();
        trajectoryCells_.clear();
        trajectoryLayerMax_.clear();
        lastWrittenTime_ = std::numeric_limits<Real>::quiet_NaN();
        nextDumpTime_ = initialStep.t + dumpStep_;
        accumulate(initialStep);
    }

	// On each accepted step, check if it falls onto the next dump time. If so, accumulate its contribution and update the next dump time.
    void onAcceptedStep(const RecordedStepData& step) override
    {
        requireActive();

        const Real eps = Real(1e-12);
        if (step.t + eps < nextDumpTime_) {
            return;
        }

        accumulate(step);
        nextDumpTime_ = step.t + dumpStep_;
    }

	// Finalize the trajectory by flushing the accumulated contributions into the result and resetting the internal state.
    RecordedTrajectoryLocation endTrajectory(const RecordedStepData& finalStep,
        double) override
    {
        requireActive();

        const Real eps = Real(1e-12);
        if (!std::isfinite(lastWrittenTime_) || std::abs(finalStep.t - lastWrittenTime_) > eps) {
            accumulate(finalStep);
        }

        flush_trajectory();
        active_ = false;
        trajectoryCellIndex_.clear();
        trajectoryCells_.clear();
        trajectoryLayerMax_.clear();

        return RecordedTrajectoryLocation{};
    }

	// If a trajectory is aborted, reset the internal state without flushing contributions, effectively discarding any accumulated data for that trajectory.
    void abortActiveTrajectory() noexcept
        override
    {
        active_ = false;
        trajectoryCellIndex_.clear();
        trajectoryCells_.clear();
        trajectoryLayerMax_.clear();
        nextDumpTime_ = Real(0);
        lastWrittenTime_ = std::numeric_limits<Real>::quiet_NaN();
    }

    MemoryRasterWorkerResult releaseResult()
    {
        return std::move(result_);
    }

private:
	// Find the slot index for the specified metric in the uniqueMetrics_ vector.
    int find_metric_slot(RasterMetric metric) const noexcept
    {
        for (std::size_t i = 0; i < uniqueMetrics_.size(); ++i) {
            if (uniqueMetrics_[i] == metric) {
                return static_cast<int>(i);
            }
        }
        return -1;
    }

    void requireActive() const
    {
        if (!active_) {
            throw std::runtime_error("memory raster recorder received step data without an active trajectory");
        }
    }

	// Evaluate all unique metrics for the current step and store their values in metricValues_ for later use in accumulating layer contributions.
    void evaluate_metrics(const RecordedStepData& step)
    {
        for (std::size_t i = 0; i < uniqueMetrics_.size(); ++i) {
            const std::optional<Real> value = sample_metric_value(uniqueMetrics_[i], step);
            metricValues_[i] = (value.has_value() && std::isfinite(*value))
                ? static_cast<float>(*value)
                : std::numeric_limits<float>::quiet_NaN();
        }
    }
	// Accumulate the contribution of the current step to the relevant raster cell and layers, updating the maximum value for order-based statistics and the streaming accumulators for streaming statistics.
    void accumulate(const RecordedStepData& step)
    {
        lastWrittenTime_ = step.t;

        const auto cellIndex = raster_cell_index(grid_, step.q(0), step.q(1));
        if (!cellIndex.has_value()) {
            return;
        }

        auto it = trajectoryCellIndex_.find(*cellIndex);
        std::size_t cellSlot = 0;
        if (it == trajectoryCellIndex_.end()) {
            cellSlot = trajectoryCells_.size();
            trajectoryCellIndex_.emplace(*cellIndex, cellSlot);
            trajectoryCells_.push_back(*cellIndex);
            trajectoryLayerMax_.resize((cellSlot + 1u) * layerCount_, -std::numeric_limits<float>::infinity());
        }
        else {
            cellSlot = it->second;
        }

        evaluate_metrics(step);
        const std::size_t base = cellSlot * layerCount_;
        for (std::size_t layerIdx = 0; layerIdx < layerCount_; ++layerIdx) {
            if (config_.layers[layerIdx].statistic == RasterStatistic::Count) {
                continue;
            }

            const int metricSlot = metricSlotByLayer_[layerIdx];
            const float value = metricValues_[static_cast<std::size_t>(metricSlot)];
            if (!std::isfinite(value)) {
                continue;
            }
            trajectoryLayerMax_[base + layerIdx] = std::max(trajectoryLayerMax_[base + layerIdx], value);
        }
    }

	// Flush the accumulated contributions for the trajectory into the result, updating the streaming accumulators and order-based contribution lists for each layer based on the maximum metric values observed for each cell.
    void flush_trajectory()
    {
        for (std::size_t cellSlot = 0; cellSlot < trajectoryCells_.size(); ++cellSlot) {
            const std::uint32_t cellIndex = trajectoryCells_[cellSlot];
            const std::size_t base = cellSlot * layerCount_;
            for (std::size_t layerIdx = 0; layerIdx < layerCount_; ++layerIdx) {
                const auto& layer = config_.layers[layerIdx];
                if (layer.statistic == RasterStatistic::Count) {
                    result_.streamingLayers[static_cast<std::size_t>(streamingSlotByLayer_[layerIdx])].accumulate_count(cellIndex);
                    continue;
                }

                const float value = trajectoryLayerMax_[base + layerIdx];
                if (!std::isfinite(value)) {
                    continue;
                }

                if (streamingSlotByLayer_[layerIdx] >= 0) {
                    result_.streamingLayers[static_cast<std::size_t>(streamingSlotByLayer_[layerIdx])]
                        .accumulate_value(cellIndex, value);
                }
                else {
                    result_.orderLayers[static_cast<std::size_t>(orderSlotByLayer_[layerIdx])]
                        .push_back(OrderContribution{ cellIndex, value });
                    ++result_.orderContributionCount;
                }
            }
        }
    }

private:
    RasterGridGeometry grid_;
    RasterOutputConfig config_;
    Real dumpStep_{ Real(0) };
    std::size_t layerCount_{ 0 };
    bool active_{ false };
    Real nextDumpTime_{ Real(0) };
    Real lastWrittenTime_{ std::numeric_limits<Real>::quiet_NaN() };
    std::vector<RasterMetric> uniqueMetrics_;
    std::vector<float> metricValues_;
    std::vector<int> metricSlotByLayer_;
    std::vector<int> streamingSlotByLayer_;
    std::vector<int> orderSlotByLayer_;
    std::unordered_map<std::uint32_t, std::size_t> trajectoryCellIndex_;
    std::vector<std::uint32_t> trajectoryCells_;
    std::vector<float> trajectoryLayerMax_;
    MemoryRasterWorkerResult result_;
};

// Compute the index of the value corresponding to the specified percentile in a sorted list of the given count.
inline std::size_t percentile_index(std::size_t count, double percentile) noexcept
{
    if (count == 0) {
        return 0;
    }
    return static_cast<std::size_t>(std::ceil(percentile * static_cast<double>(count))) - 1u;
}

// Compute the value of an order-based statistic (median or percentile) for a raster cell by partially sorting the list of values and selecting the appropriate index.
inline float compute_order_statistic_in_place(const RasterLayerSpec& layer,
    std::vector<float>& values)
{
    if (values.empty()) {
        return std::numeric_limits<float>::quiet_NaN();
    }

    std::size_t idx = 0;
    if (layer.statistic == RasterStatistic::Median) {
        idx = percentile_index(values.size(), 0.5);
    }
    else {
        const double percentile = static_cast<double>(layer.percentile.value_or(Real(95))) / 100.0;
        idx = percentile_index(values.size(), percentile);
    }
    std::nth_element(values.begin(), values.begin() + idx, values.end());
    return values[idx];
}

// Write a single-band float32 GeoTIFF raster file with the specified grid geometry, cell values, and nodata value. 
// The input data is expected to be in south-up order and will be flipped to north-up order for GeoTIFF output.
inline void write_float32_geotiff(const fs::path& path,
    const RasterGridGeometry& grid,
    const std::vector<float>& southWestData,
    float nodataValue)
{
#ifndef WITH_GDAL
    (void)path;
    (void)grid;
    (void)southWestData;
    (void)nodataValue;
    throw std::runtime_error("GeoTIFF raster output requires a build with GDAL");
#else
    if (southWestData.size() != static_cast<std::size_t>(grid.ncols) * grid.nrows) {
        throw std::runtime_error("raster output size does not match grid geometry");
    }

    GDALAllRegister();
    GDALDriverH driver = GDALGetDriverByName("GTiff");
    if (!driver) {
        throw std::runtime_error("GTiff driver is not available");
    }

    GDALDatasetH ds = GDALCreate(driver, path.string().c_str(),
        static_cast<int>(grid.ncols),
        static_cast<int>(grid.nrows),
        1,
        GDT_Float32,
        nullptr);
    if (!ds) {
        throw std::runtime_error("failed to create GeoTIFF: " + path.string());
    }

    double gt[6]{
        static_cast<double>(grid.xll),
        static_cast<double>(grid.cell),
        0.0,
        static_cast<double>(grid.yll + Real(grid.nrows) * grid.cell),
        0.0,
        -static_cast<double>(grid.cell)
    };
    if (GDALSetGeoTransform(ds, gt) != CE_None) {
        GDALClose(ds);
        throw std::runtime_error("failed to set geotransform for GeoTIFF: " + path.string());
    }

    if (!grid.projectionWkt.empty()) {
        if (GDALSetProjection(ds, grid.projectionWkt.c_str()) != CE_None) {
            GDALClose(ds);
            throw std::runtime_error("failed to set projection for GeoTIFF: " + path.string());
        }
    }

    GDALRasterBandH band = GDALGetRasterBand(ds, 1);
    GDALSetRasterNoDataValue(band, static_cast<double>(nodataValue));

    std::vector<float> northUpRow(grid.ncols);
    for (std::uint32_t row = 0; row < grid.nrows; ++row) {
        const std::uint32_t srcRow = grid.nrows - 1u - row;
        const float* src = southWestData.data() + static_cast<std::size_t>(srcRow) * grid.ncols;
        std::copy(src, src + grid.ncols, northUpRow.begin());
        if (GDALRasterIO(band, GF_Write, 0, static_cast<int>(row),
            static_cast<int>(grid.ncols), 1,
            northUpRow.data(),
            static_cast<int>(grid.ncols), 1,
            GDT_Float32, 0, 0) != CE_None) {
            GDALClose(ds);
            throw std::runtime_error("failed to write GeoTIFF row: " + path.string());
        }
    }

    GDALClose(ds);
#endif
}

// Compute the bounding box of the union of the occupied cells across multiple raster layers, where a cell is considered occupied if it has a finite value different from the nodata value in any of the layers.
inline RasterCropBounds find_union_crop_bounds(const std::vector<std::vector<float>>& rasters,
    const RasterGridGeometry& grid,
    float nodataValue) noexcept
{
    RasterCropBounds bounds;
    if (rasters.empty() || grid.ncols == 0 || grid.nrows == 0) {
        return bounds;
    }

    for (std::uint32_t row = 0; row < grid.nrows; ++row) {
        for (std::uint32_t col = 0; col < grid.ncols; ++col) {
            const std::size_t idx = static_cast<std::size_t>(row) * grid.ncols + col;
            bool occupied = false;
            for (const auto& raster : rasters) {
                const float value = raster[idx];
                if (std::isfinite(value) && value != nodataValue) {
                    occupied = true;
                    break;
                }
            }

            if (!occupied) {
                continue;
            }

            if (!bounds.hasData) {
                bounds.colMin = bounds.colMax = col;
                bounds.rowMin = bounds.rowMax = row;
                bounds.hasData = true;
            }
            else {
                bounds.colMin = std::min(bounds.colMin, col);
                bounds.colMax = std::max(bounds.colMax, col);
                bounds.rowMin = std::min(bounds.rowMin, row);
                bounds.rowMax = std::max(bounds.rowMax, row);
            }
        }
    }

    return bounds;
}

// Crop the raster grid geometry to the specified bounds, returning a new geometry that corresponds to the subset of the original grid defined by the bounds.
inline RasterGridGeometry crop_raster_grid(const RasterGridGeometry& grid,
    const RasterCropBounds& bounds) noexcept
{
    if (!bounds.hasData) {
        return grid;
    }

    RasterGridGeometry cropped = grid;
    cropped.xll = grid.xll + Real(bounds.colMin) * grid.cell;
    cropped.yll = grid.yll + Real(bounds.rowMin) * grid.cell;
    cropped.ncols = bounds.colMax - bounds.colMin + 1u;
    cropped.nrows = bounds.rowMax - bounds.rowMin + 1u;
    return cropped;
}

// Crop the raster values to the specified bounds, returning a new vector that contains only the values corresponding to the subset of the original grid defined by the bounds.
inline std::vector<float> crop_raster_values(const std::vector<float>& raster,
    const RasterGridGeometry& grid,
    const RasterCropBounds& bounds)
{
    if (!bounds.hasData) {
        return raster;
    }

    const std::uint32_t croppedCols = bounds.colMax - bounds.colMin + 1u;
    const std::uint32_t croppedRows = bounds.rowMax - bounds.rowMin + 1u;
    std::vector<float> cropped(static_cast<std::size_t>(croppedCols) * croppedRows);

    for (std::uint32_t row = 0; row < croppedRows; ++row) {
        const std::uint32_t srcRow = bounds.rowMin + row;
        const std::size_t srcOffset = static_cast<std::size_t>(srcRow) * grid.ncols + bounds.colMin;
        const std::size_t dstOffset = static_cast<std::size_t>(row) * croppedCols;
        std::copy_n(raster.data() + srcOffset, croppedCols, cropped.data() + dstOffset);
    }

    return cropped;
}

// Build the final raster layers by merging the results from multiple workers,
// finalizing the streaming statistics, computing the order-based statistics, finding the union crop bounds, cropping the rasters, 
// and writing them to GeoTIFF files in the specified output directory. 
// Returns a summary of the build process, including the number of layers and order contributions.
inline RasterBuildSummary build_rasters_from_workers(const fs::path& outputDir,
    const RasterOutputConfig& config,
    const RasterGridGeometry& grid,
    std::vector<MemoryRasterWorkerResult> workerResults)
{
    const std::size_t rasterCellCount = static_cast<std::size_t>(grid.ncols) * grid.nrows;
    const float nodataValue = -9999.0f;

    RasterBuildSummary summary;
    summary.layerCount = config.layers.size();

    std::vector<std::size_t> streamingSlotByLayer(config.layers.size(), static_cast<std::size_t>(-1));
    std::vector<std::size_t> orderSlotByLayer(config.layers.size(), static_cast<std::size_t>(-1));
    std::size_t streamingCount = 0;
    std::size_t orderCount = 0;
    for (std::size_t layerIdx = 0; layerIdx < config.layers.size(); ++layerIdx) {
        if (raster_statistic_is_streaming(config.layers[layerIdx].statistic)) {
            streamingSlotByLayer[layerIdx] = streamingCount++;
        }
        else {
            orderSlotByLayer[layerIdx] = orderCount++;
        }
    }

    std::vector<std::vector<float>> finalizedRasters(config.layers.size());
    for (std::size_t layerIdx = 0; layerIdx < config.layers.size(); ++layerIdx) {
        const auto& layer = config.layers[layerIdx];

        if (raster_statistic_is_streaming(layer.statistic)) {
            StreamingLayerAccumulator merged(layer.statistic, rasterCellCount);
            const std::size_t slot = streamingSlotByLayer[layerIdx];
            for (const auto& worker : workerResults) {
                merged.merge_from(worker.streamingLayers[slot]);
            }
            finalizedRasters[layerIdx] = merged.finalize(nodataValue);
        }
        else {
            const std::size_t slot = orderSlotByLayer[layerIdx];
            std::size_t totalRecords = 0;
            for (const auto& worker : workerResults) {
                totalRecords += worker.orderLayers[slot].size();
            }
            summary.orderContributionCount += totalRecords;

            std::vector<OrderContribution> records;
            records.reserve(totalRecords);
            for (auto& worker : workerResults) {
                auto& source = worker.orderLayers[slot];
                records.insert(records.end(),
                    std::make_move_iterator(source.begin()),
                    std::make_move_iterator(source.end()));
            }

            std::vector<float> raster(rasterCellCount, nodataValue);
            if (!records.empty()) {
                std::sort(records.begin(), records.end(), [](const OrderContribution& a, const OrderContribution& b) {
                    return a.cellIndex < b.cellIndex;
                    });

                std::vector<float> cellValues;
                std::size_t begin = 0;
                while (begin < records.size()) {
                    const std::uint32_t cellIndex = records[begin].cellIndex;
                    std::size_t end = begin;
                    while (end < records.size() && records[end].cellIndex == cellIndex) {
                        ++end;
                    }

                    cellValues.clear();
                    cellValues.reserve(end - begin);
                    for (std::size_t i = begin; i < end; ++i) {
                        cellValues.push_back(records[i].value);
                    }

                    if (cellIndex < rasterCellCount) {
                        raster[cellIndex] = compute_order_statistic_in_place(layer, cellValues);
                    }
                    begin = end;
                }
            }

            finalizedRasters[layerIdx] = std::move(raster);
        }
    }

    const RasterCropBounds cropBounds = find_union_crop_bounds(finalizedRasters, grid, nodataValue);
    const RasterGridGeometry croppedGrid = crop_raster_grid(grid, cropBounds);
    for (std::size_t layerIdx = 0; layerIdx < config.layers.size(); ++layerIdx) {
        const auto& layer = config.layers[layerIdx];
        const fs::path outPath = outputDir / (default_raster_output_stem(layer) + ".tif");
        write_float32_geotiff(
            outPath,
            croppedGrid,
            crop_raster_values(finalizedRasters[layerIdx], grid, cropBounds),
            nodataValue);
    }

    return summary;
}
