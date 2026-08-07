#include "terrain.hpp"
#include "rock.hpp"
#include "roughness.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

// A DEM query is invalid if any of the four bilinear interpolation corners is
// nodata / non-finite.
inline bool dem_cell_is_valid(const Terrain& dem,
    std::uint32_t row,
    std::uint32_t col) noexcept
{
    const Real z00 = dem.z(row, col);
    const Real z10 = dem.z(row, col + 1);
    const Real z01 = dem.z(row + 1, col);
    const Real z11 = dem.z(row + 1, col + 1);
    return std::isfinite(z00)
        && std::isfinite(z10)
        && std::isfinite(z01)
        && std::isfinite(z11);
}

// DEM gradient interpolation reuses the precomputed node gradients and only
// needs to be done for contact vertices.
inline Real dem_grad_x_at(const Terrain& dem,
    std::uint32_t row,
    std::uint32_t col,
    Real xi) noexcept
{
    const Real gx00 = dem.gradx(row, col);
    const Real gx01 = dem.gradx(row + 1, col);
    return (Real(1) - xi) * gx00 + xi * gx01;
}

// DEM gradient interpolation reuses the precomputed node gradients and only
// needs to be done for contact vertices.
inline Real dem_grad_y_at(const Terrain& dem,
    std::uint32_t row,
    std::uint32_t col,
    Real eta) noexcept
{
    const Real gy00 = dem.grady(row, col);
    const Real gy10 = dem.grady(row, col + 1);
    return (Real(1) - eta) * gy00 + eta * gy10;
}

// Specialized simulation query for the plain DEM and the roughness model. The
// loop is intentionally scalar and explicit: it follows
// the physical query order closely and avoids temporary matrix batches in the
// hot stepping path.
template <bool WithRoughness>
TerrainQueryResult query_pose_impl(const Terrain& dem,
    const Rock& rock,
    const Mat3& AIK,
    const Vec3& xI,
    TerrainQueryWorkspace& ws,
    Real hardLayerDepth)
{
    const std::size_t nVert = rock.nVert();
    ws.ensureVertexCapacity(nVert);
    ws.ensureContactCapacity(nVert);

    TerrainQueryResult result;
    result.state = TerrainContainment::Inside;
    result.minGap = std::numeric_limits<Real>::infinity();
    result.minPositiveGap = std::numeric_limits<Real>::infinity();
    result.contactCount = 0;

    const auto& vertsK = rock.verticesK();
    const Real a00 = AIK(0, 0), a01 = AIK(0, 1), a02 = AIK(0, 2);
    const Real a10 = AIK(1, 0), a11 = AIK(1, 1), a12 = AIK(1, 2);
    const Real a20 = AIK(2, 0), a21 = AIK(2, 1), a22 = AIK(2, 2);

    // First pass:
    // - transform each vertex to world coordinates
    // - classify DEM cell validity
    // - compute gap
    // - collect the contact set
    for (std::size_t v = 0; v < nVert; ++v) {
        const Real rx = vertsK(0, v);
        const Real ry = vertsK(1, v);
        const Real rz = vertsK(2, v);

        const Real x = xI.x() + a00 * rx + a01 * ry + a02 * rz;
        const Real y = xI.y() + a10 * rx + a11 * ry + a12 * rz;
        const Real z = xI.z() + a20 * rx + a21 * ry + a22 * rz;

        std::uint32_t row = 0;
        std::uint32_t col = 0;
        Real eta = Real(0);
        Real xi = Real(0);
        Real h = Real(0);
        if (!dem_height_at_xy(dem, x, y, h, &row, &col, &eta, &xi)) {
            result.state = TerrainContainment::OutOfBounds;
            return result;
        }
        if (!dem_cell_is_valid(dem, row, col)) {
            result.state = TerrainContainment::NoData;
            return result;
        }

        if constexpr (WithRoughness) {
            Real roughH = Real(0);
            roughness_model::eval_roughness_height(dem, x, y, roughH);
            h += roughH;
        }

        const Real gap = z - h;
        ws.gN(static_cast<Eigen::Index>(v)) = gap;
        result.minGap = std::min(result.minGap, gap);
        if (gap > Real(0) && gap < result.minPositiveGap) {
            result.minPositiveGap = gap;
        }

        if (gap + hardLayerDepth <= Real(0)) {
            const Eigen::Index cIdx = static_cast<Eigen::Index>(result.contactCount);
            ws.contactIdx(cIdx) = static_cast<std::uint32_t>(v);
            ws.contactRow(cIdx) = row;
            ws.contactCol(cIdx) = col;
            ws.contactEta(cIdx) = eta;
            ws.contactXi(cIdx) = xi;
            ws.contactX(cIdx) = x;
            ws.contactY(cIdx) = y;
            ++result.contactCount;
        }
    }

    // Second pass over contact vertices only: now the DEM / roughness
    // derivatives are needed to build the contact frame.
    for (std::uint32_t j = 0; j < result.contactCount; ++j) {
        const Eigen::Index idx = static_cast<Eigen::Index>(j);
        const std::uint32_t row = ws.contactRow(idx);
        const std::uint32_t col = ws.contactCol(idx);
        const Real eta = ws.contactEta(idx);
        const Real xi = ws.contactXi(idx);

        Real gx = dem_grad_x_at(dem, row, col, xi);
        Real gy = dem_grad_y_at(dem, row, col, eta);
        if constexpr (WithRoughness) {
            Real roughH = Real(0);
            Real roughDx = Real(0);
            Real roughDy = Real(0);
            roughness_model::eval_roughness_height_grad(
                dem,
                ws.contactX(idx),
                ws.contactY(idx),
                roughH,
                roughDx,
                roughDy);
            gx += roughDx;
            gy += roughDy;
        }
        ws.dhdx(idx) = gx;
        ws.dhdy(idx) = gy;
    }

    return result;
}

} // namespace

TerrainQueryResult Terrain::queryPose(const Rock& rock,
    const Mat3& AIK,
    const Vec3& xI,
    TerrainQueryWorkspace& ws,
    Real hardLayerDepth) const
{
    if (hasRoughness()) {
        return query_pose_impl<true>(*this, rock, AIK, xI, ws, hardLayerDepth);
    }
    return query_pose_impl<false>(*this, rock, AIK, xI, ws, hardLayerDepth);
}
