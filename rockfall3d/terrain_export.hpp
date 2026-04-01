#pragma once
#include "terrain.hpp"
#include "roughness.hpp"
#include "terrain_tapered.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct ExportStepEval
{
    Real value{};
};

// evaluates interpolation kernel value used when exporting roughness (linear, smoothstep, smootherstep)
inline ExportStepEval export_eval_step(StepFunc step, Real t) noexcept
{
    switch (step) {
    case StepFunc::Linear:
        return { t };
    case StepFunc::Smoothstep:
        return { t * t * (Real(3) - Real(2) * t) };
    case StepFunc::Smootherstep: {
        const Real t2 = t * t;
        const Real t3 = t2 * t;
        const Real t4 = t3 * t;
        const Real t5 = t4 * t;
        return { Real(6) * t5 - Real(15) * t4 + Real(10) * t3 };
    }
    }
    return { t };
}

// returns roughness elevation value for a given coordinate (legacy roughness model)
inline Real export_roughness_height(const Roughness& layer, Real x, Real y)
{
    const auto nx = layer.nx();
    const auto ny = layer.ny();
    if (nx < 2 || ny < 2 || layer.res <= Real(0)) return Real(0);

    const Real xmax = layer.xr(nx - 1) - layer.res;
    const Real ymax = layer.yr(ny - 1) - layer.res;
    if (x < layer.xll || y < layer.yll || x > xmax || y > ymax) return Real(0);

    const Real invRes = Real(1) / layer.res;
    const Real kf = std::floor((x - layer.xll) * invRes);
    const Real lf = std::floor((y - layer.yll) * invRes);
    if (kf < Real(0) || lf < Real(0)) return Real(0);

    const auto k = static_cast<std::uint32_t>(kf);
    const auto l = static_cast<std::uint32_t>(lf);
    if (k >= nx - 1 || l >= ny - 1) return Real(0);

    const Real dx1 = x - layer.xr(k);
    const Real dx2 = dx1 - layer.res;
    const Real dy1 = y - layer.yr(l);
    const Real dy2 = dy1 - layer.res;
    const Real eta = dx1 * invRes;
    const Real xi = dy1 * invRes;

    const Real dot1 = dx1 * layer.ax(l, k) + dy1 * layer.ay(l, k);
    const Real dot2 = dx2 * layer.ax(l, k + 1) + dy1 * layer.ay(l, k + 1);
    const Real dot3 = dx1 * layer.ax(l + 1, k) + dy2 * layer.ay(l + 1, k);
    const Real dot4 = dx2 * layer.ax(l + 1, k + 1) + dy2 * layer.ay(l + 1, k + 1);

    const Real stepEta = export_eval_step(layer.stepFunc, eta).value;
    const Real stepXi = export_eval_step(layer.stepFunc, xi).value;
    const Real line1 = dot1 + stepEta * (dot2 - dot1);
    const Real line2 = dot3 + stepEta * (dot4 - dot3);
    return line1 + stepXi * (line2 - line1);
}

} // namespace

inline void export_sampled_terrain_asc(const Terrain& dem,
    const std::string& path,
    Real resolution)
{
    if (resolution <= Real(0)) {
        throw std::runtime_error("terrain export resolution must be positive");
    }

    const Real width = Real(dem.nx - 1) * dem.cell;
    const Real height = Real(dem.ny - 1) * dem.cell;
    // Match MATLAB's xll:drawsize:x(end) and yll:drawsize:y(end) behavior:
    // keep stepping by resolution until the next point would overshoot.
    const std::uint32_t ncols = static_cast<std::uint32_t>(std::floor(width / resolution)) + 1u;
    const std::uint32_t nrows = static_cast<std::uint32_t>(std::floor(height / resolution)) + 1u;
    if (ncols < 2 || nrows < 2) {
        throw std::runtime_error("terrain export grid is too small");
    }

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot open terrain export file: " + path);
    }
    out << std::setprecision(std::numeric_limits<Real>::max_digits10);

    out << "ncols " << ncols << '\n';
    out << "nrows " << nrows << '\n';
    out << "xllcorner " << dem.xll << '\n';
    out << "yllcorner " << dem.yll << '\n';
    out << "cellsize " << resolution << '\n';
    out << "NODATA_value -9999\n";

    for (std::int64_t j = static_cast<std::int64_t>(nrows) - 1; j >= 0; --j) {
        const Real y = dem.yll + Real(j) * resolution;
        for (std::uint32_t i = 0; i < ncols; ++i) {
            const Real x = dem.xll + Real(i) * resolution;
            Real z = Real(0);
            if (!dem_height_at_xy(dem, x, y, z)) {
                throw std::runtime_error("terrain export query left DEM extent");
            }
            if (dem.roughnessBackend == RoughnessBackend::LegacyLayers) {
                for (const auto& layer : dem.roughness) {
                    z += export_roughness_height(layer, x, y);
                }
            }
            else if (dem.roughnessBackend == RoughnessBackend::TaperedCells) {
                Real h = Real(0);
                tapered_roughness::eval_tapered_cell_height(dem, x, y, h);
                z += h;
            }
            out << z;
            if (i + 1 < ncols) out << ' ';
        }
        out << '\n';
    }
}
