#pragma once
#include "math_alias.hpp"
#include <cstdint>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

class Rock;

// Interpolation kernel used by the roughness Perlin field and taper.
enum class StepFunc : std::uint8_t
{
    Linear = 0,
    Smoothstep = 1,
    Smootherstep = 2
};

enum class TerrainContainment : std::uint8_t
{
    Inside = 0,
    OutOfBounds = 1,
    NoData = 2
};

// Summary returned by the hot terrain-query kernel used during stepping.
struct TerrainQueryResult
{
    TerrainContainment state{ TerrainContainment::Inside };
    Real minGap{ std::numeric_limits<Real>::infinity() };
    Real minPositiveGap{ std::numeric_limits<Real>::infinity() };
    std::uint32_t contactCount{ 0 };
};

struct RoughnessClassParams
{
    Real gridSize{ 0 };
    Real invGridSize{ 0 };
    Real taperWidth{ 0 };
    Real invTaperWidth{ 0 };
    StepFunc stepFunc{ StepFunc::Smoothstep };
};

struct TerrainQueryWorkspace
{
    ArrX gN;
    ArrX dhdx;
    ArrX dhdy;
    ArrX contactX;
    ArrX contactY;
    ArrX contactEta;
    ArrX contactXi;
    IdxV contactIdx;
    IdxV contactRow;
    IdxV contactCol;

    // ensures workspace has enough storage for n entries. If too small, we resize once and reuse later
    // avoids repeated allocations in simulation loop
    void ensureVertexCapacity(std::size_t n)
    {
        if (static_cast<std::size_t>(gN.size()) < n) {
            gN.resize(n);
        }
    }

    void ensureContactCapacity(std::size_t n)
    {
        if (static_cast<std::size_t>(dhdx.size()) < n) {
            dhdx.resize(n);
            dhdy.resize(n);
            contactX.resize(n);
            contactY.resize(n);
            contactEta.resize(n);
            contactXi.resize(n);
            contactIdx.resize(n);
            contactRow.resize(n);
            contactCol.resize(n);
        }
    }
};

class Terrain
{
public:
    // Hot simulation query:
    // - transform/query all vertices at the midpoint qM
    // - collect the contact set
    // - compute terrain gradients for contact vertices only
    TerrainQueryResult queryPose(const Rock& rock,
        const Mat3& AIK,
        const Vec3& xI,
        TerrainQueryWorkspace& ws,
        Real hardLayerDepth = Real(0)) const;

    void addRoughness(const std::string& rasterPath,
        StepFunc stepFunc,
        std::uint64_t seed);

    bool hasRoughness() const noexcept
    {
        return roughnessCellSize > Real(0)
            && roughnessNx > 0
            && roughnessNy > 0
            && !roughnessClassParams.empty();
    }

    Real          xll{}, yll{}, cell{};
    Real          xmax{}, ymax{}, invCell{};
    std::uint32_t nx{}, ny{};
    std::string   projectionWkt;
    MatX          z;
    MatX          gradx;
    MatX          grady;

    // Internal compact representation of the roughness raster. Input raster
    // values are physical grid sizes; repeated grid sizes are stored once in
    // roughnessClassParams and referenced here by class ID. 0 means no roughness.
    Eigen::Array<std::uint16_t, Eigen::Dynamic, Eigen::Dynamic> roughnessCellClass;
    std::vector<RoughnessClassParams> roughnessClassParams;
    Real                   roughnessCellSize{ 0 };
    Real                   invRoughnessCellSize{ 0 };
    Real                   roughnessXll{ 0 }, roughnessYll{ 0 };
    std::uint32_t          roughnessNx{ 0 }, roughnessNy{ 0 };
    std::uint64_t          roughnessSeed{ 0 };
};

// Bilinear DEM interpolation once the cell and local coordinates are already
// known. This is the leanest form and is therefore used in hot paths that
// already resolved the DEM cell separately.
inline Real dem_height_at_cell(const Terrain& dem,
    std::uint32_t row,
    std::uint32_t col,
    Real eta,
    Real xi) noexcept
{
    const Real z00 = dem.z(row, col);
    const Real z10 = dem.z(row, col + 1);
    const Real z01 = dem.z(row + 1, col);
    const Real z11 = dem.z(row + 1, col + 1);
    const Real line1 = z00 + eta * (z10 - z00);
    const Real line2 = z01 + eta * (z11 - z01);
    return line1 + xi * (line2 - line1);
}

// Resolve a world-space query point and evaluate the DEM elevation with one
// shared helper. Optional outputs allow callers in the hot loop to reuse the
// resolved cell information without redoing the lookup.
inline bool dem_height_at_xy(const Terrain& dem,
    Real x,
    Real y,
    Real& h,
    std::uint32_t* rowOut = nullptr,
    std::uint32_t* colOut = nullptr,
    Real* etaOut = nullptr,
    Real* xiOut = nullptr) noexcept
{
    if (x < dem.xll || y < dem.yll || x > dem.xmax || y > dem.ymax) {
        return false;
    }

    const Real colf = std::floor((x - dem.xll) * dem.invCell);
    const Real rowf = std::floor((y - dem.yll) * dem.invCell);
    if (colf < Real(0) || rowf < Real(0)) {
        return false;
    }

    std::uint32_t col = static_cast<std::uint32_t>(colf);
    std::uint32_t row = static_cast<std::uint32_t>(rowf);
    if (col >= dem.nx - 1) col = dem.nx - 2;
    if (row >= dem.ny - 1) row = dem.ny - 2;

    const Real eta = (x - (dem.xll + Real(col) * dem.cell)) * dem.invCell;
    const Real xi = (y - (dem.yll + Real(row) * dem.cell)) * dem.invCell;
    h = dem_height_at_cell(dem, row, col, eta, xi);

    if (rowOut) *rowOut = row;
    if (colOut) *colOut = col;
    if (etaOut) *etaOut = eta;
    if (xiOut) *xiOut = xi;
    return true;
}

