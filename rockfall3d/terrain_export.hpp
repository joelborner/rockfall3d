#pragma once
#include "terrain.hpp"
#include "roughness.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

// Exports the terrain to an ASCII raster file with the given resolution. 
// The output grid is aligned with the DEM grid and covers the same extent, but may have fewer rows and columns. 
// The output height values include roughness if present. 
// The output file is in a simple format with a header followed by rows of space-separated height values, starting from the top row (maximum y) to the bottom row (minimum y).
inline void export_sampled_terrain_asc(const Terrain& dem,
    const std::string& path,
    Real resolution)
{
    if (resolution <= Real(0)) {
        throw std::runtime_error("terrain export resolution must be positive");
    }
    const Real width = Real(dem.nx - 1) * dem.cell;
    const Real height = Real(dem.ny - 1) * dem.cell;
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
