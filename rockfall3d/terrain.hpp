#pragma once
#include "math_alias.hpp"
#include <cmath>
#include <limits>
#include <random>
#include <vector>

class Roughness;                 // forward declaration
class Rock;
enum class StepFunc : std::uint8_t;   // from roughness.hpp

// Legacy layered roughness and tapered-cell roughness use different fast
// lookup paths. The enum is kept explicit so the simulation code can dispatch
// once and then stay in a specialized hot path.
enum class RoughnessLookupMode : std::uint8_t
{
    Optimized = 0,
    Legacy = 1,
    BruteForce = 2
};

enum class RoughnessBackend : std::uint8_t
{
    None = 0,
    LegacyLayers = 1,
    TaperedCells = 2
};

enum class TerrainContainment : std::uint8_t
{
    Inside = 0,
    OutOfBounds = 1,
    NoData = 2
};

struct TerrainQueryResult
{
    // Summary returned by the hot terrain-query kernel used during stepping.
    TerrainContainment state{ TerrainContainment::Inside };
    Real minGap{ std::numeric_limits<Real>::infinity() };
    Real minPositiveGap{ std::numeric_limits<Real>::infinity() };
    std::uint32_t contactCount{ 0 };
};

struct TaperedClassParams
{
    Real gridSize{ 0 };
    Real invGridSize{ 0 };
    Real taperWidth{ 0 };
    Real invTaperWidth{ 0 };
};

struct TerrainQueryWorkspace
{
    // IrOP is kept for the compatibility path that still routes through the
    // older matrix-based gap() implementation.
    MatX IrOP;
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
        if (IrOP.cols() < static_cast<Eigen::Index>(n)) {
            IrOP.resize(3, n);
        }
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
    // Classify only checks geometric admissibility of a transformed rock:
    // rectangular DEM bounds and nodata-free interpolation cells. It does not
    // compute gap values or gradients. Only for legacy layered roughness model
    template <typename Derived>
    TerrainContainment classify(const Eigen::MatrixBase<Derived>& IrOP) const noexcept
    {
        const auto x = IrOP.row(0).array();
        const auto y = IrOP.row(1).array();

        const Real xmin = x.minCoeff();
        const Real xmax = x.maxCoeff();
        const Real ymin = y.minCoeff();
        const Real ymax = y.maxCoeff();

        if (!(xmin >= xll &&
            xmax <= this->xmax &&
            ymin >= yll &&
            ymax <= this->ymax)) {
            return TerrainContainment::OutOfBounds;
        }

        for (Eigen::Index p = 0; p < IrOP.cols(); ++p) {
            const Real xp = IrOP(0, p);
            const Real yp = IrOP(1, p);

            const Real colf = std::floor((xp - xll) * invCell);
            const Real rowf = std::floor((yp - yll) * invCell);
            if (colf < Real(0) || rowf < Real(0)) {
                return TerrainContainment::OutOfBounds;
            }

            std::uint32_t col = static_cast<std::uint32_t>(colf);
            std::uint32_t row = static_cast<std::uint32_t>(rowf);
            if (col >= nx - 1) col = nx - 2;
            if (row >= ny - 1) row = ny - 2;

            const Real z00 = z(row, col);
            const Real z10 = z(row, col + 1);
            const Real z01 = z(row + 1, col);
            const Real z11 = z(row + 1, col + 1);
            if (!std::isfinite(z00) || !std::isfinite(z10)
                || !std::isfinite(z01) || !std::isfinite(z11)) {
                return TerrainContainment::NoData;
            }
        }

        return TerrainContainment::Inside;
    }

    template <typename Derived>
    bool inside(const Eigen::MatrixBase<Derived>& IrOP) const noexcept
    {
        return classify(IrOP) == TerrainContainment::Inside;
    }

    void gap(const Eigen::Ref<const MatX>& IrOP,
        ArrX& gN,
        ArrX& dhdx,
        ArrX& dhdy,
        IdxV& I_N) const;

    // Hot simulation query:
    // - transform/query all vertices at the midpoint qM
    // - collect the contact set
    // - compute terrain gradients for contact vertices only
    TerrainQueryResult queryPose(const Rock& rock,
        const Mat3& AIK,
        const Vec3& xI,
        TerrainQueryWorkspace& ws) const;

    // legacy layered roughness
    void addRoughness(const std::vector<std::string>& paths,
        const std::vector<Real>& res,
        StepFunc stepFunc,
        std::mt19937* rng = nullptr);

    // new optimised roughness
    void addTaperedCellRoughness(const std::string& classRasterPath,
        const std::vector<Real>& gridSizes,
        std::uint64_t seed);

    // bitMask used for legacy roughness
    template <class UInt = std::uint64_t>
    Eigen::Array<UInt, Eigen::Dynamic, Eigen::Dynamic>
        buildBitMask(const std::vector<ArrXXb>& layers) const;

    Real          xll{}, yll{}, cell{};
    Real          xmax{}, ymax{}, invCell{};
    std::uint32_t nx{}, ny{};
    MatX          z;
    MatX          gradx;
    MatX          grady;

    std::vector<Roughness> roughness;           // legacy roughness
    Real                   minrough{ 0 };       // legacy roughness
    std::uint32_t roughNx{ 0 }, roughNy{ 0 };   // legacy roughness
    Eigen::Array<std::uint64_t, Eigen::Dynamic, Eigen::Dynamic> roughIndices;   // legacy roughness

    // Coarse roughness classification grid: each cell stores a class ID, not
    // the physical grid size itself.
    Eigen::Array<std::uint16_t, Eigen::Dynamic, Eigen::Dynamic> taperedCellClass;
    std::vector<TaperedClassParams> taperedClassParams;
    std::vector<Real>      taperedGridSizes;
    Real                   taperedCellSize{ 0 };
    Real                   invTaperedCellSize{ 0 };
    Real                   taperedXll{ 0 }, taperedYll{ 0 };
    std::uint32_t          taperedNx{ 0 }, taperedNy{ 0 };
    std::uint64_t          taperedSeed{ 0 };

    bool shapes{ false };
    bool hasRoughness{ false };
    RoughnessBackend roughnessBackend{ RoughnessBackend::None };
    RoughnessLookupMode roughnessLookupMode{ RoughnessLookupMode::Optimized };
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

