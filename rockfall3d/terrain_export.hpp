#pragma once
#include "terrain.hpp"
#include "roughness.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

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
            if (dem.hasRoughness()) {
                Real h = Real(0);
                roughness_model::eval_roughness_height(dem, x, y, h);
                z += h;
            }
            out << z;
            if (i + 1 < ncols) out << ' ';
        }
        out << '\n';
    }
}
