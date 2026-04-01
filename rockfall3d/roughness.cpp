// Implementation of the Roughness class declared in roughness.hpp.
#include "roughness.hpp"
#include "terrain.hpp"
#include <memory>
#include <iostream>
#include <limits>
#include <fstream>
#include <cmath>

#ifdef WITH_GDAL
#  include "ogrsf_frmts.h"
#  include "gdal_alg.h"
#endif

/*
inline void save_mask_pgm(const ArrXXb& mask, const std::string& path, bool flipY = false)
{
    const int H = static_cast<int>(mask.rows());
    const int W = static_cast<int>(mask.cols());

    std::ofstream out(path, std::ios::binary);
    if (!out) throw std::runtime_error("save_mask_pgm: cannot open " + path);

    out << "P5\n" << W << " " << H << "\n255\n";

    for (int j = 0; j < H; ++j) {
        const int row = flipY ? (H - 1 - j) : j;
        for (int i = 0; i < W; ++i) {
            unsigned char v = mask(row, i) ? 255 : 0;
            out.write(reinterpret_cast<char*>(&v), 1);
        }
    }
}*/

namespace {

    // -------------- Sutherland–Hodgman polygon clip in image space -------------
    struct Pt { double x, y; };
    using PolyD = std::vector<Pt>;

    static PolyD clip_poly_rect(const PolyD& in,
        double xmin, double xmax,
        double ymin, double ymax)
    {
        auto clip_edge = [](const PolyD& poly,
            auto inside, auto intersect) -> PolyD
            {
                PolyD out;
                if (poly.empty()) return out;

                Pt S = poly.back();
                bool S_in = inside(S);
                for (const Pt& E : poly) {
                    bool E_in = inside(E);
                    if (S_in && E_in) {
                        // S -> E fully inside
                        out.push_back(E);
                    }
                    else if (S_in && !E_in) {
                        // leaving – add intersection
                        out.push_back(intersect(S, E));
                    }
                    else if (!S_in && E_in) {
                        // entering – add intersection + E
                        out.push_back(intersect(S, E));
                        out.push_back(E);
                    }
                    S = E; S_in = E_in;
                }
                return out;
            };

        PolyD p = in;

        // left
        p = clip_edge(p,
            [xmin](const Pt& P) { return P.x >= xmin; },
            [xmin](const Pt& S, const Pt& E) {
                double t = (xmin - S.x) / (E.x - S.x);
                return Pt{ xmin, S.y + t * (E.y - S.y) };
            });

        // right
        p = clip_edge(p,
            [xmax](const Pt& P) { return P.x <= xmax; },
            [xmax](const Pt& S, const Pt& E) {
                double t = (xmax - S.x) / (E.x - S.x);
                return Pt{ xmax, S.y + t * (E.y - S.y) };
            });

        // top (remember: MATLAB img space is top-down, row 0 at the top)
        p = clip_edge(p,
            [ymin](const Pt& P) { return P.y >= ymin; },
            [ymin](const Pt& S, const Pt& E) {
                double t = (ymin - S.y) / (E.y - S.y);
                return Pt{ S.x + t * (E.x - S.x), ymin };
            });

        // bottom
        p = clip_edge(p,
            [ymax](const Pt& P) { return P.y <= ymax; },
            [ymax](const Pt& S, const Pt& E) {
                double t = (ymax - S.y) / (E.y - S.y);
                return Pt{ S.x + t * (E.x - S.x), ymax };
            });

        return p;
    }

} // unnamed namespace
// Helper: load polygons from a shapefile  (GDAL / OGR)
static std::vector<PolyXY> load_polygons_from_shapefile(const std::string& path)
{
    std::vector<PolyXY> polys;

#ifndef WITH_GDAL
    throw std::runtime_error(
        "Shapefile support requested but the code was built without WITH_GDAL.");
#else
    GDALAllRegister();
    std::unique_ptr<GDALDataset, decltype(&GDALClose)>
        ds(static_cast<GDALDataset*>(GDALOpenEx(path.c_str(),
            GDAL_OF_VECTOR, nullptr,
            nullptr, nullptr)),
            &GDALClose);
    if (!ds)
        throw std::runtime_error("GDALOpenEx failed for shapefile: " + path);

    OGRLayer* layer = ds->GetLayer(0);
    if (!layer)
        throw std::runtime_error("Shapefile has no layers: " + path);

    layer->ResetReading();
    while (OGRFeature* feat = layer->GetNextFeature()) {
        OGRGeometry* geom = feat->GetGeometryRef();
        if (!geom) { OGRFeature::DestroyFeature(feat); continue; }

        auto add_ring = [&polys](OGRLinearRing* ring)
            {
                if (!ring) return;
                PolyXY p;
                p.x.reserve(ring->getNumPoints());
                p.y.reserve(ring->getNumPoints());
                for (int i = 0; i < ring->getNumPoints(); ++i) {
                    p.x.push_back(static_cast<Real>(ring->getX(i)));
                    p.y.push_back(static_cast<Real>(ring->getY(i)));
                }
                polys.emplace_back(std::move(p));
            };

        switch (wkbFlatten(geom->getGeometryType())) {
        case wkbPolygon: {
            auto* pg = geom->toPolygon();
            add_ring(pg->getExteriorRing());
            break;
        }
        case wkbMultiPolygon: {
            auto* mpg = geom->toMultiPolygon();
            for (auto&& g : *mpg) {
                add_ring(g->toPolygon()->getExteriorRing());
            }
            break;
        }
        default:
            std::cerr << "Skipping non-polygon geometry in " << path << '\n';
            break;
        }

        OGRFeature::DestroyFeature(feat);
    }
    return polys;
#endif
}

void Roughness::init_common(const Terrain& dem,
    Real           res_,
    StepFunc       step,
    std::uint32_t  layerIndex)
{
    res = res_;
    stepFunc = step;
    index = layerIndex;

    // Build 1D coordinate vectors
    const std::size_t nxLoc =
        static_cast<std::size_t>(std::floor((xur - xll) / res)) + 2;
    const std::size_t nyLoc =
        static_cast<std::size_t>(std::floor((yur - yll) / res)) + 2;

    xr.resize(nxLoc);
    yr.resize(nyLoc);
    for (std::size_t i = 0; i < nxLoc; ++i) xr(i) = xll + Real(i) * res;
    for (std::size_t j = 0; j < nyLoc; ++j) yr(j) = yll + Real(j) * res;

    // Mesh grid
    X.resize(nyLoc, nxLoc);
    Y.resize(nyLoc, nxLoc);
    for (std::size_t j = 0; j < nyLoc; ++j)
        for (std::size_t i = 0; i < nxLoc; ++i) {
            X(j, i) = xr(i);
            Y(j, i) = yr(j);
        }

    // Allocate gradient arrays
    ax.setZero(nyLoc, nxLoc);
    ay.setZero(nyLoc, nxLoc);
}

void Roughness::seed_random_unit_vectors(std::mt19937& gen,
    const ArrXXb& activeMask)
{
    std::uniform_real_distribution<Real> dist(-0.5, 0.5);

    const Eigen::Index nyLoc = ax.rows();
    const Eigen::Index nxLoc = ax.cols();

    for (Eigen::Index j = 0; j < nyLoc; ++j)
        for (Eigen::Index i = 0; i < nxLoc; ++i) {
            if (!activeMask(j, i)) continue;          // leave at zero
            Real vx = dist(gen);
            Real vy = dist(gen);
            Real n = std::sqrt(vx * vx + vy * vy);
            if (n <= Real(0)) n = Real(1);
            ax(j, i) = vx / n;
            ay(j, i) = vy / n;
        }
}

// Constructors
Roughness::Roughness(const Terrain& dem,
    Real           res_,
    StepFunc       step,
    std::uint32_t  layerIndex,
    std::mt19937* rng)
{
    // Extent full DEM
    xll = dem.xll;
    yll = dem.yll;
    xur = dem.xll + dem.cell * (dem.nx - 1);
    yur = dem.yll + dem.cell * (dem.ny - 1);

    init_common(dem, res_, step, layerIndex);

    // Active everywhere
    ArrXXb mask(ax.rows(), ax.cols());
    mask.setConstant(true);

    std::mt19937 gen = rng ? *rng : std::mt19937{ std::random_device{}() };
    seed_random_unit_vectors(gen, mask);
}

Roughness::Roughness(const Terrain& dem,
    Real              res_,
    const std::string& shpPath,
    StepFunc          step,
    std::uint32_t     layerIndex,
    std::mt19937* rng)
{
    bound = load_polygons_from_shapefile(shpPath);
    if (bound.empty())
        std::cerr << "Warning: shapefile has no polygons: " << shpPath << '\n';

    // Bounding box of all polygons (fallback to DEM if none)
    BBox bb;
    for (const auto& p : bound) {
        auto b = bbox_of(p);
        if (!b.valid()) continue;
        bb.xmin = std::min(bb.xmin, b.xmin);
        bb.xmax = std::max(bb.xmax, b.xmax);
        bb.ymin = std::min(bb.ymin, b.ymin);
        bb.ymax = std::max(bb.ymax, b.ymax);
    }
    if (!bb.valid()) {
        bb.xmin = dem.xll;
        bb.ymin = dem.yll;
        bb.xmax = dem.xll + dem.cell * (dem.nx - 1);
        bb.ymax = dem.yll + dem.cell * (dem.ny - 1);
    }

    // Snap to DEM and add one cell margin
    const Real snap = dem.minrough > Real(0) ? dem.minrough : dem.cell;
    xll = std::floor((dem.xll-std::max(dem.xll, bb.xmin - res_)) / snap) * snap+dem.xll;
    yll = std::floor((dem.yll-std::max(dem.yll, bb.ymin - res_)) / snap) * snap+dem.yll;
    xur = std::min(dem.xll + dem.cell * dem.nx, bb.xmax) + res_;
    yur = std::min(dem.yll + dem.cell * dem.ny, bb.ymax) + res_;

    init_common(dem, res_, step, layerIndex);

    // Rasterise polygons onto this grid
    ArrXXb mask = checkPolygonInside(X, Y);

    std::mt19937 gen = rng ? *rng : std::mt19937{ std::random_device{}() };
    seed_random_unit_vectors(gen, mask);
}

ArrXXb Roughness::checkPolygonInside(const MatX& Xg, const MatX& Yg) const
{
    const int rows = static_cast<int>(Xg.rows());
    const int cols = static_cast<int>(Xg.cols());

    ArrXXb mask(rows, cols);
    mask.setConstant(false);

    if (bound.empty()) {
        mask.setConstant(true);
        return mask;
    }

    // ---- map to MATLAB image space (top-down) ------------------------------
    const Real xMin = Xg.minCoeff();
    const Real xMax = Xg.maxCoeff();
    const Real yMin = Yg.minCoeff();
    const Real yMax = Yg.maxCoeff();

    const double sx = (cols > 1) ? double(cols - 1) / double(xMax - xMin) : 1.0;
    const double sy = (rows > 1) ? double(rows - 1) / double(yMax - yMin) : 1.0;

    auto to_img_x = [&](Real px) {
        return double(std::llround((px - xMin) * sx));
        };
    auto to_img_y_topdown = [&](Real py) {
        return double(std::llround((py - yMin) * sy)); // 0 at top, like MATLAB
        };

    std::vector<double> xints;
    xints.reserve(256);

    const double xmin = -0.5;
    const double xmax = cols - 0.5;
    const double ymin = -0.5;
    const double ymax = rows - 0.5;

    for (const auto& P : bound)
    {
        const std::size_t n = P.x.size();
        if (n < 3) continue;

        // Original polygon in image space (continuous, unclamped)
        PolyD poly;
        poly.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            poly.push_back({ to_img_x(P.x[i]), to_img_y_topdown(P.y[i]) });

        // Clip to the raster rectangle FIRST
        PolyD C = clip_poly_rect(poly, xmin, xmax, ymin, ymax);
        if (C.size() < 3) continue;

        // Integer scanline span (top-down)
        int minRow = rows - 1, maxRow = 0;
        for (auto& v : C) {
            int yi = int(std::floor(v.y)); // conservative
            minRow = std::min(minRow, yi);
            maxRow = std::max(maxRow, yi);
        }
        minRow = std::max(minRow, 0);
        maxRow = std::min(maxRow, rows - 1);
        if (minRow > maxRow) continue;

        // Scanline fill on the clipped polygon
        for (int y = minRow; y <= maxRow; ++y)
        {
            xints.clear();
            const double yc = double(y) + 0.5;   // pixel center

            for (std::size_t i = 0, j = C.size() - 1; i < C.size(); j = i++)
            {
                const double y0 = C[j].y, y1 = C[i].y;
                const double x0 = C[j].x, x1 = C[i].x;

                // canonical even-odd crossing test
                if (((y0 > yc) != (y1 > yc))) {
                    const double xi = x0 + (yc - y0) * (x1 - x0) / (y1 - y0);
                    xints.push_back(xi);
                }
            }

            if (xints.empty()) continue;
            std::sort(xints.begin(), xints.end());

            for (std::size_t k = 0; k + 1 < xints.size(); k += 2) {
                int xStart = int(std::ceil(std::min(xints[k], xints[k + 1])));
                int xEnd = int(std::floor(std::max(xints[k], xints[k + 1])));
                if (xEnd < xStart) continue;

                xStart = std::max(0, std::min(cols - 1, xStart));
                xEnd = std::max(0, std::min(cols - 1, xEnd));

                for (int x = xStart; x <= xEnd; ++x)
                    mask(y, x) = true;
            }
        }
    }
    return mask;
}


void Roughness::removeLocalRoughness(const Roughness& over)
{
    ArrXXb mask = over.checkPolygonInside(X, Y);
    for (Eigen::Index j = 0; j < ax.rows(); ++j)
        for (Eigen::Index i = 0; i < ax.cols(); ++i)
            if (mask(j, i)) { ax(j, i) = Real(0); ay(j, i) = Real(0); }
}
