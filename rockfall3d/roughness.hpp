// Per-cell Perlin-noise roughness layer that can cover the whole DEM or be
// clipped to polygons loaded from a shapefile.
//
// Only declarations live here.  Implementations go in roughness.cpp.
#pragma once
#include "math_alias.hpp"
#include "polygon_utils.hpp"

#include <vector>
#include <string>
#include <random>
#include <cstdint>

class Terrain;                      // forward declaration

// Interpolation kernel to be used later when sampling Perlin noise
enum class StepFunc : std::uint8_t
{
    Linear = 0,
    Smoothstep = 1,
    Smootherstep = 2
};

class Roughness
{
public:
    std::vector<PolyXY> bound;      // polygons that delimit the layer
    Real xll{}, yll{};              // lower-left world coordinate of layer
    Real xur{}, yur{};              // upper-right world coordinate
    Real res{};                     // grid spacing of this roughness layer

    MatX X, Y;                      // node coordinates (ny x nx)
    ArrX xr, yr;                    // 1D coordinate vectors (cols / rows)

    MatX ax, ay;                    // unit gradient vectors at nodes
    StepFunc stepFunc{ StepFunc::Linear };

    std::uint32_t index{ 0 };         // layer index

    Roughness() = default;

    // Overall layer that covers the whole DEM
    Roughness(const Terrain& dem,
        Real res,
        StepFunc step,
        std::uint32_t layerIndex,
        std::mt19937* rng = nullptr);

    // Polygon-bounded layer loaded from a shapefile
    Roughness(const Terrain& dem,
        Real res,
        const std::string& shapefilePath,
        StepFunc step,
        std::uint32_t layerIndex,
        std::mt19937* rng = nullptr);

    // Zero-out tuples that are overwritten by a later layer
    void removeLocalRoughness(const Roughness& over);

    // Rasterise the layer's polygons on an arbitrary node grid (utility)
    ArrXXb checkPolygonInside(const MatX& Xg, const MatX& Yg) const;

    // Node counts
    std::uint32_t nx() const noexcept { return static_cast<std::uint32_t>(xr.size()); }
    std::uint32_t ny() const noexcept { return static_cast<std::uint32_t>(yr.size()); }

private:
    // Shared initialisation logic
    void init_common(const Terrain& dem,
        Real res,
        StepFunc step,
        std::uint32_t layerIndex);

    // Random unit-vector seeding
    void seed_random_unit_vectors(std::mt19937& gen,
        const ArrXXb& activeMask);
};
