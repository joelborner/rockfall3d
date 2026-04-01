#pragma once

#include "terrain.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

namespace tapered_roughness {

struct StepEval
{
    Real value{};
    Real deriv{};
};

inline StepEval smoothstep01(Real t) noexcept
{
    t = std::clamp(t, Real(0), Real(1));
    return {
        t * t * (Real(3) - Real(2) * t),
        Real(6) * t * (Real(1) - t)
    };
}

// Small integer mixer used to derive pseudo-random but deterministic lattice
// gradients from (seed, ix, iy).
inline std::uint64_t mix64(std::uint64_t x) noexcept
{
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ULL;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebULL;
    x ^= x >> 31;
    return x;
}

// Gradient selection for one Perlin lattice node. Continuity comes from the
// fact that the same lattice coordinates always map to the same gradient.
inline void hashed_unit_gradient(std::uint64_t seed,
    std::int64_t ix,
    std::int64_t iy,
    Real& gx,
    Real& gy) noexcept
{
    static constexpr Real invSqrt2 = Real(0.70710678118654752440);
    static constexpr std::array<std::array<Real, 2>, 8> gradients{ {
        // instead of generating a hashed random unit gradient in every time step for every vertex, we randomly select a predefined unit gradient from a list
        // this is way faster, not truly random but variable enough for the roughness to look authentic
        { Real(1), Real(0) },
        { Real(-1), Real(0) },
        { Real(0), Real(1) },
        { Real(0), Real(-1) },
        { invSqrt2, invSqrt2 },
        { invSqrt2, -invSqrt2 },
        { -invSqrt2, invSqrt2 },
        { -invSqrt2, -invSqrt2 }
    } };

    std::uint64_t h = seed;
    h ^= mix64(static_cast<std::uint64_t>(ix) + 0x9e3779b97f4a7c15ULL);
    h ^= mix64(static_cast<std::uint64_t>(iy) + 0xbf58476d1ce4e5b9ULL);
    const auto& grad = gradients[mix64(h) & 7ULL];
    gx = grad[0];
    gy = grad[1];
}

// Evaluate one globally defined Perlin field and, optionally, its spatial
// derivatives. "Global" here means that the field is a function of world
// coordinates, class parameters, and seed; it is not restarted per coarse
// roughness cell.
inline void eval_global_perlin_impl(Real x,
    Real y,
    const TaperedClassParams& params,
    std::uint64_t seed,
    Real& value,
    Real* ddx,
    Real* ddy) noexcept
{
    const Real res = params.gridSize;
    if (res <= Real(0)) {
        value = Real(0);
        if (ddx) *ddx = Real(0);
        if (ddy) *ddy = Real(0);
        return;
    }

    const Real invRes = params.invGridSize;
    const Real xf = x * invRes;
    const Real yf = y * invRes;
    const std::int64_t ix = static_cast<std::int64_t>(std::floor(xf));
    const std::int64_t iy = static_cast<std::int64_t>(std::floor(yf));
    const Real eta = xf - Real(ix);
    const Real xi = yf - Real(iy);

    Real g00x, g00y, g10x, g10y, g01x, g01y, g11x, g11y;
    hashed_unit_gradient(seed, ix, iy, g00x, g00y);
    hashed_unit_gradient(seed, ix + 1, iy, g10x, g10y);
    hashed_unit_gradient(seed, ix, iy + 1, g01x, g01y);
    hashed_unit_gradient(seed, ix + 1, iy + 1, g11x, g11y);

    const Real dx1 = eta * res;
    const Real dx2 = dx1 - res;
    const Real dy1 = xi * res;
    const Real dy2 = dy1 - res;

    const Real dot1 = g00x * dx1 + g00y * dy1;
    const Real dot2 = g10x * dx2 + g10y * dy1;
    const Real dot3 = g01x * dx1 + g01y * dy2;
    const Real dot4 = g11x * dx2 + g11y * dy2;

    const StepEval sEta = smoothstep01(eta);
    const StepEval sXi = smoothstep01(xi);
    const Real line1 = dot1 + sEta.value * (dot2 - dot1);
    const Real line2 = dot3 + sEta.value * (dot4 - dot3);
    value = line1 + sXi.value * (line2 - line1);

    if (!ddx || !ddy) return;

    const Real di1dx = g00x
        + sEta.deriv * invRes * (dot2 - dot1)
        + sEta.value * (g10x - g00x);
    const Real di2dx = g01x
        + sEta.deriv * invRes * (dot4 - dot3)
        + sEta.value * (g11x - g01x);
    const Real di1dy = g00y + sEta.value * (g10y - g00y);
    const Real di2dy = g01y + sEta.value * (g11y - g01y);

    *ddx = di1dx + sXi.value * (di2dx - di1dx);
    *ddy = di1dy + sXi.value * (di2dy - di1dy)
        + sXi.deriv * invRes * (line2 - line1);
}

inline void eval_global_perlin(Real x,
    Real y,
    const TaperedClassParams& params,
    std::uint64_t seed,
    Real& value,
    Real* ddx = nullptr,
    Real* ddy = nullptr) noexcept
{
    eval_global_perlin_impl(x, y, params, seed, value, ddx, ddy);
}

struct CellLookupResult
{
    std::uint32_t col{ 0 };
    std::uint32_t row{ 0 };
    const TaperedClassParams* params{ nullptr };
};

// Resolve a world-space query point to one coarse roughness class cell and the
// corresponding precomputed class parameters. Returns the Perlin grid size of the active roughness layer
inline bool cell_lookup(const Terrain& dem,
    Real x,
    Real y,
    CellLookupResult& cell) noexcept
{
    if (dem.invTaperedCellSize <= Real(0) || dem.taperedNx == 0 || dem.taperedNy == 0) {
        return false;
    }

    const Real colf = std::floor((x - dem.taperedXll) * dem.invTaperedCellSize);
    const Real rowf = std::floor((y - dem.taperedYll) * dem.invTaperedCellSize);
    if (colf < Real(0) || rowf < Real(0)) {
        return false;
    }

    cell.col = static_cast<std::uint32_t>(colf);
    cell.row = static_cast<std::uint32_t>(rowf);
    if (cell.col >= dem.taperedNx || cell.row >= dem.taperedNy) {
        return false;
    }

    const std::uint16_t cls = dem.taperedCellClass(cell.row, cell.col);
    if (cls == 0) {
        return false;
    }

    const std::size_t idx = static_cast<std::size_t>(cls - 1);
    if (idx >= dem.taperedClassParams.size()) {
        return false;
    }

    cell.params = &dem.taperedClassParams[idx];
    return cell.params->gridSize > Real(0);
}

// Height-only tapered roughness evaluation. This is the dominant path because
// all vertices need heights, while only contact vertices need derivatives.
inline void eval_tapered_cell_height(const Terrain& dem,
    Real x,
    Real y,
    Real& value) noexcept
{
    CellLookupResult cell;
    if (!cell_lookup(dem, x, y, cell)) {
        value = Real(0);
        return;
    }

    const TaperedClassParams& params = *cell.params;
    const Real x0 = dem.taperedXll + Real(cell.col) * dem.taperedCellSize;
    const Real x1 = x0 + dem.taperedCellSize;
    const Real y0 = dem.taperedYll + Real(cell.row) * dem.taperedCellSize;
    const Real y1 = y0 + dem.taperedCellSize;

    Real perlin = Real(0);
    eval_global_perlin(x, y, params, dem.taperedSeed, perlin);

    if (params.taperWidth <= Real(0)) {
        value = Real(0);
        return;
    }

    const Real dx1 = x - x0;
    const Real dx2 = x1 - x;
    const Real dy1 = y - y0;
    const Real dy2 = y1 - y;
    const Real dxMin = std::min(dx1, dx2);
    const Real dyMin = std::min(dy1, dy2);

    if (dxMin >= params.taperWidth && dyMin >= params.taperWidth) {
        value = perlin;
        return;
    }

    const StepEval sx = smoothstep01(dxMin * params.invTaperWidth);
    const StepEval sy = smoothstep01(dyMin * params.invTaperWidth);
    const Real taper = sx.value * sy.value;
    value = taper * perlin;
}

// Height + derivatives for contact vertices. Only called for contact set. The taper derivative is only
// active inside the edge band. in the cell interior the result reduces to the
// raw Perlin gradient.
inline void eval_tapered_cell_height_grad(const Terrain& dem,
    Real x,
    Real y,
    Real& value,
    Real& ddx,
    Real& ddy) noexcept
{
    CellLookupResult cell;
    if (!cell_lookup(dem, x, y, cell)) {
        value = Real(0);
        ddx = Real(0);
        ddy = Real(0);
        return;
    }

    const TaperedClassParams& params = *cell.params;
    const Real x0 = dem.taperedXll + Real(cell.col) * dem.taperedCellSize;
    const Real x1 = x0 + dem.taperedCellSize;
    const Real y0 = dem.taperedYll + Real(cell.row) * dem.taperedCellSize;
    const Real y1 = y0 + dem.taperedCellSize;

    Real perlin = Real(0);
    Real perlinDx = Real(0);
    Real perlinDy = Real(0);
    eval_global_perlin(x, y, params, dem.taperedSeed, perlin, &perlinDx, &perlinDy);

    if (params.taperWidth <= Real(0)) {
        value = Real(0);
        ddx = Real(0);
        ddy = Real(0);
        return;
    }

    const Real dx1 = x - x0;
    const Real dx2 = x1 - x;
    const Real dy1 = y - y0;
    const Real dy2 = y1 - y;
    const Real dxMin = std::min(dx1, dx2);
    const Real dyMin = std::min(dy1, dy2);

    if (dxMin >= params.taperWidth && dyMin >= params.taperWidth) {
        value = perlin;
        ddx = perlinDx;
        ddy = perlinDy;
        return;
    }

    const StepEval sx = smoothstep01(dxMin * params.invTaperWidth);
    const StepEval sy = smoothstep01(dyMin * params.invTaperWidth);
    const Real taper = sx.value * sy.value;
    value = taper * perlin;

    Real taperDx = Real(0);
    Real taperDy = Real(0);
    if (dx1 <= dx2 && dx1 < params.taperWidth) taperDx = sx.deriv * params.invTaperWidth * sy.value;
    else if (dx2 < dx1 && dx2 < params.taperWidth) taperDx = -sx.deriv * params.invTaperWidth * sy.value;

    if (dy1 <= dy2 && dy1 < params.taperWidth) taperDy = sy.deriv * params.invTaperWidth * sx.value;
    else if (dy2 < dy1 && dy2 < params.taperWidth) taperDy = -sy.deriv * params.invTaperWidth * sx.value;

    ddx = taperDx * perlin + taper * perlinDx;
    ddy = taperDy * perlin + taper * perlinDy;
}

// Convenience wrapper retained for call sites that want the older optional
// derivative signature.
inline void eval_tapered_cell_point(const Terrain& dem,
    Real x,
    Real y,
    Real& value,
    Real* ddx = nullptr,
    Real* ddy = nullptr) noexcept
{
    if (!ddx && !ddy) {
        eval_tapered_cell_height(dem, x, y, value);
        return;
    }

    Real gx = Real(0);
    Real gy = Real(0);
    eval_tapered_cell_height_grad(dem, x, y, value, gx, gy);
    if (ddx) *ddx = gx;
    if (ddy) *ddy = gy;
}

// Batch helper kept for compatibility with the older matrix-based gap() path.
inline void eval_tapered_cell_batch(const Terrain& dem,
    const Eigen::Ref<const MatX>& pts,
    ArrX& h,
    ArrX* dhdx = nullptr,
    ArrX* dhdy = nullptr)
{
    const std::size_t P = static_cast<std::size_t>(pts.cols());
    h.resize(P);
    h.setZero();
    if (dhdx) {
        dhdx->resize(P);
        dhdx->setZero();
    }
    if (dhdy) {
        dhdy->resize(P);
        dhdy->setZero();
    }
    for (std::size_t p = 0; p < P; ++p) {
        if (dhdx && dhdy) {
            Real gx = Real(0);
            Real gy = Real(0);
            eval_tapered_cell_height_grad(dem, pts(0, p), pts(1, p), h(p), gx, gy);
            (*dhdx)(p) = gx;
            (*dhdy)(p) = gy;
        }
        else {
            eval_tapered_cell_height(dem, pts(0, p), pts(1, p), h(p));
        }
    }
}

} // namespace tapered_roughness
