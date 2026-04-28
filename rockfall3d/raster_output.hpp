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

struct RasterGridGeometry
{
    Real xll{ Real(0) };
    Real yll{ Real(0) };
    Real cell{ Real(2) };
    std::uint32_t ncols{ 0 };
    std::uint32_t nrows{ 0 };
    std::string projectionWkt;
};

inline RasterGridGeometry make_raster_grid_geometry(const Terrain& dem, Real resolution)
{
    if (!(resolution > Real(0))) {
        throw std::runtime_error("raster output resolution must be positive");
    }

    RasterGridGeometry grid;
    grid.xll = dem.xll;
    grid.yll = dem.yll;
    grid.cell = resolution;
    grid.projectionWkt = dem.projectionWkt;

    const Real width = Real(dem.nx - 1) * dem.cell;
    const Real height = Real(dem.ny - 1) * dem.cell;
    grid.ncols = static_cast<std::uint32_t>(std::floor(width / resolution)) + 1u;
    grid.nrows = static_cast<std::uint32_t>(std::floor(height / resolution)) + 1u;
    if (grid.ncols == 0 || grid.nrows == 0) {
        throw std::runtime_error("raster output grid is empty");
    }
    return grid;
}

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

struct RasterContribution
{
    std::uint32_t cellIndex{ 0 };
    float ekinTot{ 0.0f };
    float jumpH{ 0.0f };
};

struct TrajectoryCellAggregate
{
    Real maxEkinTot{ -std::numeric_limits<Real>::infinity() };
    Real maxJumpH{ -std::numeric_limits<Real>::infinity() };
};

class MemoryRasterRecorder final : public TrajectoryRecorder
{
public:
    MemoryRasterRecorder(const RasterGridGeometry& grid,
        Real dumpStep)
        : grid_(grid),
        dumpStep_(dumpStep)
    {
        if (!(dumpStep_ > Real(0))) {
            throw std::runtime_error("memory raster output requires a positive dump step");
        }
    }

    void beginTrajectory(const TrajectoryMetadata&,
        const Rock& rock,
        const RecordedStepData& initialStep) override
    {
        if (active_) {
            throw std::runtime_error("memory raster recorder beginTrajectory called while another trajectory is still active");
        }

        currentRock_ = &rock;
        active_ = true;
        visited_.clear();
        lastWrittenTime_ = std::numeric_limits<Real>::quiet_NaN();
        nextDumpTime_ = initialStep.t + dumpStep_;
        accumulate(initialStep);
    }

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
        currentRock_ = nullptr;
        visited_.clear();

        return RecordedTrajectoryLocation{};
    }

    void abortActiveTrajectory() noexcept
    {
        active_ = false;
        currentRock_ = nullptr;
        visited_.clear();
        nextDumpTime_ = Real(0);
        lastWrittenTime_ = std::numeric_limits<Real>::quiet_NaN();
    }

    const std::vector<RasterContribution>& contributions() const noexcept
    {
        return contributions_;
    }

    std::vector<RasterContribution> releaseContributions()
    {
        return std::move(contributions_);
    }

private:
    void requireActive() const
    {
        if (!active_) {
            throw std::runtime_error("memory raster recorder received step data without an active trajectory");
        }
    }

    Real total_kinetic_energy(const RecordedStepData& step) const
    {
        return Real(0.5) * currentRock_->mass() * step.u.head<3>().squaredNorm()
            + Real(0.5) * step.u.tail<3>().dot(currentRock_->inertiaK() * step.u.tail<3>());
    }

    void accumulate(const RecordedStepData& step)
    {
        lastWrittenTime_ = step.t;

        const auto cellIndex = raster_cell_index(grid_, step.q(0), step.q(1));
        if (!cellIndex.has_value()) {
            return;
        }

        auto& slot = visited_[*cellIndex];
        slot.maxEkinTot = std::max(slot.maxEkinTot, total_kinetic_energy(step));
        if (step.minGap.has_value() && std::isfinite(*step.minGap)) {
            slot.maxJumpH = std::max(slot.maxJumpH, *step.minGap);
        }
    }

    void flush_trajectory()
    {
        contributions_.reserve(contributions_.size() + visited_.size());
        for (const auto& entry : visited_) {
            contributions_.push_back(RasterContribution{
                entry.first,
                static_cast<float>(entry.second.maxEkinTot),
                static_cast<float>(entry.second.maxJumpH)
                });
        }
    }

private:
    RasterGridGeometry grid_;
    Real dumpStep_{ Real(0) };
    bool active_{ false };
    const Rock* currentRock_{ nullptr };
    Real nextDumpTime_{ Real(0) };
    Real lastWrittenTime_{ std::numeric_limits<Real>::quiet_NaN() };
    std::unordered_map<std::uint32_t, TrajectoryCellAggregate> visited_;
    std::vector<RasterContribution> contributions_;
};

inline std::size_t percentile_index_q95(std::size_t count) noexcept
{
    if (count == 0) {
        return 0;
    }
    return static_cast<std::size_t>(std::ceil(0.95 * static_cast<double>(count))) - 1u;
}

inline float compute_q95_in_place(std::vector<float>& values)
{
    if (values.empty()) {
        return std::numeric_limits<float>::quiet_NaN();
    }
    const std::size_t idx = percentile_index_q95(values.size());
    std::nth_element(values.begin(), values.begin() + idx, values.end());
    return values[idx];
}

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

inline void build_q95_rasters_from_memory(const fs::path& outputDir,
    const RasterGridGeometry& grid,
    std::vector<RasterContribution> records)
{
    const std::size_t rasterCellCount = static_cast<std::size_t>(grid.ncols) * grid.nrows;
    const float nodataValue = -9999.0f;
    std::vector<float> q95Ekin(rasterCellCount, nodataValue);
    std::vector<float> q95JumpH(rasterCellCount, nodataValue);

    if (!records.empty()) {
        std::sort(records.begin(), records.end(), [](const RasterContribution& a, const RasterContribution& b) {
            return a.cellIndex < b.cellIndex;
            });

        std::vector<float> ekinValues;
        std::vector<float> jumpValues;
        std::size_t begin = 0;
        while (begin < records.size()) {
            const std::uint32_t cellIndex = records[begin].cellIndex;
            std::size_t end = begin;
            while (end < records.size() && records[end].cellIndex == cellIndex) {
                ++end;
            }

            ekinValues.clear();
            jumpValues.clear();
            ekinValues.reserve(end - begin);
            jumpValues.reserve(end - begin);
            for (std::size_t i = begin; i < end; ++i) {
                ekinValues.push_back(records[i].ekinTot);
                jumpValues.push_back(records[i].jumpH);
            }

            if (cellIndex < rasterCellCount) {
                q95Ekin[cellIndex] = compute_q95_in_place(ekinValues);
                q95JumpH[cellIndex] = compute_q95_in_place(jumpValues);
            }
            begin = end;
        }
    }

    write_float32_geotiff(outputDir / "q95_ekin_tot.tif", grid, q95Ekin, nodataValue);
    write_float32_geotiff(outputDir / "q95_jump_h.tif", grid, q95JumpH, nodataValue);
}
