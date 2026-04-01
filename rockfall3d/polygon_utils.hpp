// Minimal polygon containers and rasterisation helpers (header-only).

#pragma once
#include "math_alias.hpp"
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <stdexcept>

// XY coordinates of polygon
struct PolyXY {
    std::vector<Real> x; 
    std::vector<Real> y;
};

// bounding-box of polygon
struct BBox {
    Real xmin{ std::numeric_limits<Real>::infinity() };
    Real ymin{ std::numeric_limits<Real>::infinity() };
    Real xmax{ -std::numeric_limits<Real>::infinity() };
    Real ymax{ -std::numeric_limits<Real>::infinity() };
    bool valid() const noexcept { return xmin <= xmax && ymin <= ymax; }
};


// returns bounding-box of polygon
inline BBox bbox_of(const PolyXY& p)
{
    BBox b; if (p.x.empty()) return b;
    for (std::size_t i = 0; i < p.x.size(); ++i) {
        b.xmin = std::min(b.xmin, p.x[i]); b.xmax = std::max(b.xmax, p.x[i]);
        b.ymin = std::min(b.ymin, p.y[i]); b.ymax = std::max(b.ymax, p.y[i]);
    }
    return b;
}


// Returns whether point is in polygon or not
// Works for simple or self-intersecting polygons
inline bool pip_evenodd(Real px, Real py, const PolyXY& poly) noexcept
{
    const std::size_t n = poly.x.size();
    if (n < 3) return false;
    bool inside = false;
    std::size_t j = n - 1;
    for (std::size_t i = 0; i < n; ++i) {
        const Real xi = poly.x[i];
        const Real yi = poly.y[i];
        const Real xj = poly.x[j];
        const Real yj = poly.y[j];
        const bool intersect = ((yi > py) != (yj > py)) &&
            (px < (xj - xi) * (py - yi) / (yj - yi + Real(1e-300)) + xi);
        if (intersect) inside = !inside;
        j = i;
    }
    return inside;
}


// Rasterise polygon set onto node grid 
// Returns ny x nx logical array where each element corresponds to a node at
// (xll + i*res, yll + j*res) that lies inside any polygon
inline ArrXXb rasterize_polygons_nodes(
    Real xll, Real yll, Real res,
    std::uint32_t nx, std::uint32_t ny,
    const std::vector<PolyXY>& polys)
{
    ArrXXb mask(ny, nx); mask.setConstant(false);
    if (polys.empty()) return mask;

    // Precompute bounding boxes for acceleration
    std::vector<BBox> boxes; boxes.reserve(polys.size());
    for (auto const& p : polys) boxes.push_back(bbox_of(p));

    for (std::uint32_t j = 0; j < ny; ++j) {
        const Real y = yll + Real(j) * res;
        for (std::uint32_t i = 0; i < nx; ++i) {
            const Real x = xll + Real(i) * res;
            bool in = false;
            for (std::size_t k = 0; k < polys.size(); ++k) {
                auto const& b = boxes[k];
                if (x < b.xmin || x > b.xmax || y < b.ymin || y > b.ymax) continue;
                if (pip_evenodd(x, y, polys[k])) { in = true; break; }
            }
            mask(j, i) = in;
        }
    }
    return mask;
}
