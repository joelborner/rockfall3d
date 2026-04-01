#include "terrain.hpp"
#include "rock.hpp"
#include "roughness.hpp"
#include "terrain_tapered.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

struct StepEval
{
    Real value{};
    Real deriv{};
};

inline StepEval eval_step(StepFunc step, Real t) noexcept
{
    switch (step) {
    case StepFunc::Linear:
        return { t, Real(1) };
    case StepFunc::Smoothstep:
        return {
            t * t * (Real(3) - Real(2) * t),
            Real(6) * t * (Real(1) - t)
        };
    case StepFunc::Smootherstep: {
        const Real t2 = t * t;
        const Real t3 = t2 * t;
        const Real t4 = t3 * t;
        const Real t5 = t4 * t;
        return {
            Real(6) * t5 - Real(15) * t4 + Real(10) * t3,
            Real(30) * t4 - Real(60) * t3 + Real(30) * t2
        };
    }
    }
    return { t, Real(1) };
}

// Evaluate one legacy roughness layer on a batch of query points.
// This path is kept for compatibility with the older multi-layer roughness
// model where gradients are pre-generated and stored explicitly.
void eval_roughness_layer(const Roughness& layer,
    const Eigen::Ref<const MatX>& IrOP,
    ArrX& dh,
    ArrX* ddhdx = nullptr,
    ArrX* ddhdy = nullptr)
{
    const std::size_t P = IrOP.cols();
    dh.resize(P);
    dh.setZero();
    if (ddhdx) {
        ddhdx->resize(P);
        ddhdx->setZero();
    }
    if (ddhdy) {
        ddhdy->resize(P);
        ddhdy->setZero();
    }

    const auto nx = layer.nx();
    const auto ny = layer.ny();
    if (nx < 2 || ny < 2 || layer.res <= Real(0)) return;

    const Real x0 = layer.xll;
    const Real y0 = layer.yll;
    const Real xmax = x0 + Real(nx - 1) * layer.res;
    const Real ymax = y0 + Real(ny - 1) * layer.res;
    const Real invRes = Real(1) / layer.res;

    for (std::size_t p = 0; p < P; ++p) {
        const Real x = IrOP(0, p);
        const Real y = IrOP(1, p);
        if (x < x0 || y < y0 || x >= xmax || y >= ymax) continue;

        const Real kf = std::floor((x - x0) * invRes);
        const Real lf = std::floor((y - y0) * invRes);
        if (kf < Real(0) || lf < Real(0)) continue;

        const auto k = static_cast<std::uint32_t>(kf);
        const auto l = static_cast<std::uint32_t>(lf);
        if (k >= nx - 1 || l >= ny - 1) continue;

        const Real xk = x0 + Real(k) * layer.res;
        const Real yl = y0 + Real(l) * layer.res;
        const Real dx1 = x - xk;
        const Real dx2 = dx1 - layer.res;
        const Real dy1 = y - yl;
        const Real dy2 = dy1 - layer.res;
        const Real eta = dx1 * invRes;
        const Real xi = dy1 * invRes;

        const StepEval sEta = eval_step(layer.stepFunc, eta);
        const StepEval sXi = eval_step(layer.stepFunc, xi);

        const Real ax00 = layer.ax(l, k);
        const Real ay00 = layer.ay(l, k);
        const Real ax10 = layer.ax(l, k + 1);
        const Real ay10 = layer.ay(l, k + 1);
        const Real ax01 = layer.ax(l + 1, k);
        const Real ay01 = layer.ay(l + 1, k);
        const Real ax11 = layer.ax(l + 1, k + 1);
        const Real ay11 = layer.ay(l + 1, k + 1);

        const Real dot1 = ax00 * dx1 + ay00 * dy1;
        const Real dot2 = ax10 * dx2 + ay10 * dy1;
        const Real dot3 = ax01 * dx1 + ay01 * dy2;
        const Real dot4 = ax11 * dx2 + ay11 * dy2;

        const Real line1 = dot1 + sEta.value * (dot2 - dot1);
        const Real line2 = dot3 + sEta.value * (dot4 - dot3);
        dh(p) = line1 + sXi.value * (line2 - line1);

        if (!ddhdx || !ddhdy) continue;

        const Real di1dx = ax00
            + sEta.deriv * invRes * (dot2 - dot1)
            + sEta.value * (ax10 - ax00);
        const Real di2dx = ax01
            + sEta.deriv * invRes * (dot4 - dot3)
            + sEta.value * (ax11 - ax01);
        const Real di1dy = ay00 + sEta.value * (ay10 - ay00);
        const Real di2dy = ay01 + sEta.value * (ay11 - ay01);

        (*ddhdx)(p) = di1dx + sXi.value * (di2dx - di1dx);
        (*ddhdy)(p) = di1dy + sXi.value * (di2dy - di1dy)
            + sXi.deriv * invRes * (line2 - line1);
    }
}

std::uint64_t full_layer_mask(std::size_t layerCount) noexcept
{
    if (layerCount == 0) return 0;
    if (layerCount >= 64) return std::numeric_limits<std::uint64_t>::max();
    return (std::uint64_t(1) << layerCount) - 1;
}

std::vector<std::uint64_t> roughness_point_masks(
    const Terrain& dem,
    const Eigen::Ref<const MatX>& IrOP)
{
    std::vector<std::uint64_t> masks(static_cast<std::size_t>(IrOP.cols()), 0);
    if (!dem.hasRoughness || dem.roughnessBackend != RoughnessBackend::LegacyLayers || dem.roughness.empty()) return masks;

    if (!dem.shapes || dem.roughIndices.size() == 0 || dem.minrough <= Real(0)) {
        const std::uint64_t mask = full_layer_mask(dem.roughness.size());
        std::fill(masks.begin(), masks.end(), mask);
        return masks;
    }

    for (Eigen::Index p = 0; p < IrOP.cols(); ++p) {
        const Real x = IrOP(0, p);
        const Real y = IrOP(1, p);
        const Real iflt = std::floor((x - dem.xll) / dem.minrough);
        const Real jflt = std::floor((y - dem.yll) / dem.minrough);
        if (iflt < Real(0) || jflt < Real(0)) continue;

        const auto i = static_cast<Eigen::Index>(iflt);
        const auto j = static_cast<Eigen::Index>(jflt);
        if (i < 0 || j < 0
            || i >= dem.roughIndices.cols()
            || j >= dem.roughIndices.rows()) {
            continue;
        }

        masks[static_cast<std::size_t>(p)] = dem.roughIndices(j, i);
    }

    return masks;
}

std::vector<std::size_t> active_roughness_layers(
    const Terrain& dem,
    const std::vector<std::uint64_t>& pointMasks)
{
    std::vector<std::size_t> layers;
    if (!dem.hasRoughness || dem.roughness.empty()) return layers;

    std::uint64_t unionMask = 0;
    for (const std::uint64_t mask : pointMasks) unionMask |= mask;

    for (std::size_t bit = 0; bit < dem.roughness.size() && bit < 64; ++bit) {
        if (unionMask & (std::uint64_t(1) << bit)) layers.push_back(bit);
    }

    return layers;
}

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

inline Real dem_grad_y_at(const Terrain& dem,
    std::uint32_t row,
    std::uint32_t col,
    Real eta) noexcept
{
    const Real gy00 = dem.grady(row, col);
    const Real gy10 = dem.grady(row, col + 1);
    return (Real(1) - eta) * gy00 + eta * gy10;
}

// Specialized simulation query for the plain DEM and the tapered-cell
// roughness backend. The loop is intentionally scalar and explicit: it follows
// the physical query order closely and avoids temporary matrix batches in the
// hot stepping path.
template <bool WithTaperedRoughness>
TerrainQueryResult query_pose_impl(const Terrain& dem,
    const Rock& rock,
    const Mat3& AIK,
    const Vec3& xI,
    TerrainQueryWorkspace& ws)
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

        if constexpr (WithTaperedRoughness) {
            Real roughH = Real(0);
            tapered_roughness::eval_tapered_cell_height(dem, x, y, roughH);
            h += roughH;
        }

        const Real gap = z - h;
        ws.gN(static_cast<Eigen::Index>(v)) = gap;
        result.minGap = std::min(result.minGap, gap);
        if (gap > Real(0) && gap < result.minPositiveGap) {
            result.minPositiveGap = gap;
        }

        if (gap <= Real(0)) {
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
        if constexpr (WithTaperedRoughness) {
            Real roughH = Real(0);
            Real roughDx = Real(0);
            Real roughDy = Real(0);
            tapered_roughness::eval_tapered_cell_height_grad(
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
    TerrainQueryWorkspace& ws) const
{
    // Legacy layered roughness still reuses the older matrix-based gap()
    // implementation because that path already contains its own optimized
    // masking logic.
    if (roughnessBackend == RoughnessBackend::LegacyLayers) {
        const std::size_t nVert = rock.nVert();
        ws.ensureVertexCapacity(nVert);
        ws.IrOP.leftCols(static_cast<Eigen::Index>(nVert)).noalias() = AIK * rock.verticesK();
        ws.IrOP.block(0, 0, 1, static_cast<Eigen::Index>(nVert)).array() += xI.x();
        ws.IrOP.block(1, 0, 1, static_cast<Eigen::Index>(nVert)).array() += xI.y();
        ws.IrOP.block(2, 0, 1, static_cast<Eigen::Index>(nVert)).array() += xI.z();

        TerrainQueryResult result;
        result.state = classify(ws.IrOP.leftCols(static_cast<Eigen::Index>(nVert)));
        if (result.state != TerrainContainment::Inside) {
            return result;
        }

        ArrX dhdxTmp, dhdyTmp;
        IdxV I_N;
        gap(ws.IrOP.leftCols(static_cast<Eigen::Index>(nVert)), ws.gN, dhdxTmp, dhdyTmp, I_N);
        result.minGap = ws.gN.head(static_cast<Eigen::Index>(nVert)).minCoeff();
        result.minPositiveGap = std::numeric_limits<Real>::infinity();
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(nVert); ++i) {
            const Real gapValue = ws.gN(i);
            if (gapValue > Real(0) && gapValue < result.minPositiveGap) {
                result.minPositiveGap = gapValue;
            }
        }
        result.contactCount = static_cast<std::uint32_t>(I_N.size());
        ws.ensureContactCapacity(result.contactCount);
        for (std::uint32_t j = 0; j < result.contactCount; ++j) {
            const Eigen::Index idx = static_cast<Eigen::Index>(j);
            ws.contactIdx(idx) = I_N(idx);
            ws.dhdx(idx) = dhdxTmp(idx);
            ws.dhdy(idx) = dhdyTmp(idx);
        }
        return result;
    }

    if (roughnessBackend == RoughnessBackend::TaperedCells) {
        return query_pose_impl<true>(*this, rock, AIK, xI, ws);
    }
    return query_pose_impl<false>(*this, rock, AIK, xI, ws);
}

// legacy gap function for old roughness model
void Terrain::gap(const Eigen::Ref<const MatX>& IrOP,
    ArrX& gN, ArrX& dhdx, ArrX& dhdy, IdxV& I_N) const
{
    // Matrix-based compatibility path used by:
    // - legacy layered roughness
    // - non-hot-path helpers such as export / initialization
    const std::size_t P = IrOP.cols();

    // grid indices
    ArrX colf = ((IrOP.row(0).array() - xll) * invCell).floor();   // k-1
    ArrX rowf = ((IrOP.row(1).array() - yll) * invCell).floor();   // l-1
    IdxV col = colf.cast<std::uint32_t>();
    IdxV row = rowf.cast<std::uint32_t>();

    // keep cell inside DEM
    col = (col < nx - 1).select(col, nx - 2);
    row = (row < ny - 1).select(row, ny - 2);

    // barycentric coords inside the cell
    ArrX eta = (IrOP.row(0).transpose().array() - (xll + col.cast<Real>() * cell)) * invCell;
    ArrX xi = (IrOP.row(1).transpose().array() - (yll + row.cast<Real>() * cell)) * invCell;

    // bilinear height
    ArrX h(P);
    for (std::size_t p = 0; p < P; ++p)
    {
        const unsigned r = row(p), c = col(p);

        const Real z00 = z(r, c);
        const Real z10 = z(r, c + 1);
        const Real z01 = z(r + 1, c);
        const Real z11 = z(r + 1, c + 1);

        h(p) = dem_height_at_cell(*this, r, c, eta(p), xi(p));
    }

    ArrX roughH = ArrX::Zero(P);
    std::vector<std::uint64_t> pointMasks;
    std::vector<std::size_t> activeLayers;
    if (roughnessBackend == RoughnessBackend::LegacyLayers) {
        pointMasks = roughness_point_masks(*this, IrOP);
        activeLayers = active_roughness_layers(*this, pointMasks);
        if (roughnessLookupMode == RoughnessLookupMode::Optimized
            || roughnessLookupMode == RoughnessLookupMode::Legacy) {
            ArrX layerH;
            for (const std::size_t layerIdx : activeLayers) {
                eval_roughness_layer(roughness[layerIdx], IrOP, layerH);
                roughH += layerH;
            }
        }
        else {
            for (std::size_t layerIdx = 0; layerIdx < roughness.size(); ++layerIdx) {
                ArrX layerH;
                eval_roughness_layer(roughness[layerIdx], IrOP, layerH);
                roughH += layerH;
            }
        }
    }
    else if (roughnessBackend == RoughnessBackend::TaperedCells) {
        tapered_roughness::eval_tapered_cell_batch(*this, IrOP, roughH);
    }

    h += roughH;
    gN = IrOP.row(2).transpose().array() - h;

    // contact indices
    std::vector<std::uint32_t> contactIdx;
    contactIdx.reserve(P);
    for (std::size_t p = 0; p < P; ++p)
        if (gN(p) <= Real(0))
            contactIdx.push_back(static_cast<std::uint32_t>(p));

    I_N.resize(contactIdx.size());
    for (std::size_t i = 0; i < contactIdx.size(); ++i)
        I_N(i) = contactIdx[i];

    const std::size_t nC = I_N.size();
    dhdx.resize(nC);
    dhdy.resize(nC);
    if (nC == 0) return;

    ArrX roughDx = ArrX::Zero(nC);
    ArrX roughDy = ArrX::Zero(nC);
    MatX contactPts(3, nC);
    for (std::size_t j = 0; j < nC; ++j) {
        contactPts.col(j) = IrOP.col(I_N(j));
    }

    if (roughnessBackend == RoughnessBackend::LegacyLayers) {
        if (roughnessLookupMode == RoughnessLookupMode::Optimized
            || roughnessLookupMode == RoughnessLookupMode::Legacy) {
            ArrX layerH, layerDx, layerDy;
            for (const std::size_t layerIdx : activeLayers) {
                eval_roughness_layer(roughness[layerIdx], contactPts, layerH, &layerDx, &layerDy);
                roughDx += layerDx;
                roughDy += layerDy;
            }
        }
        else {
            for (std::size_t layerIdx = 0; layerIdx < roughness.size(); ++layerIdx) {
                ArrX layerH, layerDx, layerDy;
                eval_roughness_layer(roughness[layerIdx], contactPts, layerH, &layerDx, &layerDy);
                roughDx += layerDx;
                roughDy += layerDy;
            }
        }
    }
    else if (roughnessBackend == RoughnessBackend::TaperedCells) {
        ArrX dummyH;
        tapered_roughness::eval_tapered_cell_batch(*this, contactPts, dummyH, &roughDx, &roughDy);
    }

    /* ---- gradients (bilinear) ---------------------------------------- */
    for (std::size_t j = 0; j < nC; ++j)
    {
        const unsigned r = row(I_N(j));
        const unsigned c = col(I_N(j));

        const Real gx00 = gradx(r, c);
        const Real gx01 = gradx(r + 1, c);
        dhdx(j) = (Real(1) - xi(I_N(j))) * gx00 + xi(I_N(j)) * gx01 + roughDx(j);

        const Real gy00 = grady(r, c);
        const Real gy10 = grady(r, c + 1);
        dhdy(j) = (Real(1) - eta(I_N(j))) * gy00 + eta(I_N(j)) * gy10 + roughDy(j);
    }
}



