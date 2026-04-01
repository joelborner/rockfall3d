// Adds Perlin noise roughness layers to an existing Terrain instance and
// packs a fast lookup bit mask.
#include "terrain.hpp"          // extended Terrain declaration
#include "roughness.hpp"
#include "polygon_utils.hpp"
#include "math_alias.hpp"

#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>
#include <limits>

// Helper: choose an unsigned-integer type wide enough for k layers
static std::size_t choose_mask_bits(std::size_t k) noexcept
{
    if (k <= 8)   return 8;
    if (k <= 16)  return 16;
    if (k <= 32)  return 32;
    return 64;                      // clamp at 64 layers
}

// Turn a node mask (ny×nx) into a cell mask ((ny-1)×(nx-1))
// cell(i,j) = OR of its 4 node corners
static ArrXXb nodes_to_cells(const ArrXXb& node)
{
    const Eigen::Index H = node.rows() - 1;
    const Eigen::Index W = node.cols() - 1;
    ArrXXb cell(H, W);
    cell = node.block(0, 0, H, W)
        || node.block(1, 0, H, W)
        || node.block(0, 1, H, W)
        || node.block(1, 1, H, W);
    return cell;
}

// Nearest-neighbour upsampling by an *integer* factor `f >= 1`
static ArrXXb upsample_nearest(const ArrXXb& src, std::uint32_t f)
{
    if (f == 1) return src;
    const Eigen::Index Hs = src.rows();
    const Eigen::Index Ws = src.cols();
    ArrXXb dst(Hs * f, Ws * f);
    for (Eigen::Index j = 0; j < dst.rows(); ++j) {
        const Eigen::Index js = j / f;
        for (Eigen::Index i = 0; i < dst.cols(); ++i) {
            const Eigen::Index is = i / f;
            dst(j, i) = src(js, is);
        }
    }
    return dst;
}

// Paste `src` into `dst` (both bool arrays), with integer offset (ox,oy),
// cropping when src overflows dst. Returns nothing (in-place update).
static void blit_with_offset(const ArrXXb& src,
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>& dst,
    int ox, int oy)
{
    const int Wd = static_cast<int>(dst.cols());
    const int Hd = static_cast<int>(dst.rows());
    const int Ws = static_cast<int>(src.cols());
    const int Hs = static_cast<int>(src.rows());

    int minx = std::max(0, ox);
    int miny = std::max(0, oy);
    int maxx = std::min(Wd - 1, ox + Ws - 1);
    int maxy = std::min(Hd - 1, oy + Hs - 1);
    if (minx > maxx || miny > maxy) return;

    int u0 = std::max(0, -ox);
    int v0 = std::max(0, -oy);

    for (int y = 0; y <= (maxy - miny); ++y) {
        const int dy = miny + y;
        const int sy = v0 + y;
        for (int x = 0; x <= (maxx - minx); ++x) {
            const int dx = minx + x;
            const int sx = u0 + x;
            dst(dy, dx) = dst(dy, dx) || src(sy, sx);
        }
    }
}

void Terrain::addRoughness(const std::vector<std::string>& paths,
    const std::vector<Real>& res,
    StepFunc                        stepFunc,
    std::mt19937*                   rng)
{
    if (paths.size() != res.size())
        throw std::runtime_error("addRoughness: paths/res size mismatch");

    taperedCellClass.resize(0, 0);
    taperedGridSizes.clear();
    taperedCellSize = Real(0);
    taperedXll = taperedYll = Real(0);
    taperedNx = taperedNy = 0;
    taperedSeed = 0;

    const std::size_t L = res.size();          // no. of layers
    if (L == 0) {
        roughness.clear();
        roughIndices.resize(0, 0);
        hasRoughness = shapes = false;
        roughnessBackend = RoughnessBackend::None;
        minrough = cell;
        return;
    }

    //determine minrough
    minrough = *std::min_element(res.begin(), res.end());

    //build layers
    roughness.clear();
    roughness.reserve(L);

    std::mt19937 localRng = rng ? *rng : std::mt19937{ std::random_device{}() };
    std::mt19937& roughRng = rng ? *rng : localRng;

    for (std::size_t r = 0; r < L; ++r) {
        const std::string& p = paths[r];
        const Real         rr = res[r];

        std::cout << "Adding roughness layer " << r << " (" << rr << " m)\n";

        if (p == "overall" || p == "OVERALL" || p == "Overall")
            roughness.emplace_back(*this, rr, stepFunc,
                static_cast<std::uint32_t>(r), &roughRng);
        else
            roughness.emplace_back(*this, rr, p, stepFunc,
                static_cast<std::uint32_t>(r), &roughRng);

        // remove overlaps from earlier layers
        if (r > 0)
            for (std::size_t q = 0; q < r; ++q)
                roughness[q].removeLocalRoughness(roughness[r]);
    }

    // global minrough grid
    roughNx = static_cast<std::uint32_t>(
        std::round(((nx - 1) * cell) / minrough));
    roughNy = static_cast<std::uint32_t>(
        std::round(((ny - 1) * cell) / minrough));
    roughNx = std::max<std::uint32_t>(roughNx, 1u);
    roughNy = std::max<std::uint32_t>(roughNy, 1u);

    // Cube of logical layers at minrough resolution H x W x L
    std::vector<ArrXXb> cube;
    cube.reserve(L);

    // Size of the global minroughness cell grid
    roughNx = static_cast<std::uint32_t>(
        std::llround(((nx - 1) * cell) / minrough));
    roughNy = static_cast<std::uint32_t>(
        std::llround(((ny - 1) * cell) / minrough));
    roughNx = std::max<std::uint32_t>(roughNx, 1u);
    roughNy = std::max<std::uint32_t>(roughNy, 1u);

    // Pre-allocate cube layers to false
    for (std::size_t r = 0; r < L; ++r) {
        cube.emplace_back(roughNy, roughNx);
        cube.back().setConstant(false);
    }

    for (std::size_t r = 0; r < L; ++r)
    {
        const Roughness& R = roughness[r];

        // node mask: non-zero tuples active
        ArrXXb node = (R.ax.array() != Real(0)) || (R.ay.array() != Real(0));

        // cells mask (Hc x Wc)
        ArrXXb cells = nodes_to_cells(node);

        // upsample to global minrough grid using nearest neighbour
        //    (assume R.res is an integer multiple of minrough)
        const double ratio = static_cast<double>(R.res / minrough);
        const std::uint32_t f = static_cast<std::uint32_t>(std::llround(ratio));
        if (std::abs(ratio - double(f)) > 1e-6) {
            throw std::runtime_error("R.res must be an integer multiple of minrough");
        }
        ArrXXb im = upsample_nearest(cells, f);

        // compute integer offset (in cell indices at minrough resolution)
        const int offx = static_cast<int>(std::llround((R.xll - xll) / minrough));
        const int offy = static_cast<int>(std::llround((R.yll - yll) / minrough));

        // paste into the r-th layer of the cube
        blit_with_offset(im, cube[r], offx, offy);
    }


    // build packed mask
    shapes = std::any_of(paths.begin(), paths.end(),
        [](const std::string& s)
        { return s != "overall" && s != "OVERALL" && s != "Overall"; });

    hasRoughness =
        std::any_of(res.begin(), res.end(), [](Real v) { return v > Real(0); });
    roughnessBackend = hasRoughness
        ? RoughnessBackend::LegacyLayers
        : RoughnessBackend::None;

    // Pick smallest integer class that fits
    const std::size_t bits = choose_mask_bits(L);
    if (bits == 8)
        roughIndices = buildBitMask<std::uint8_t>(cube).template cast<std::uint64_t>();
    else if (bits == 16)
        roughIndices = buildBitMask<std::uint16_t>(cube).template cast<std::uint64_t>();
    else if (bits == 32)
        roughIndices = buildBitMask<std::uint32_t>(cube).template cast<std::uint64_t>();
    else
        roughIndices = buildBitMask<std::uint64_t>(cube);
}

// Terrain::buildBitMask – template definition + explicit instantiations
template <class UInt>
Eigen::Array<UInt, Eigen::Dynamic, Eigen::Dynamic>
Terrain::buildBitMask(const std::vector<ArrXXb>& layers) const
{
    if (layers.empty())
        return Eigen::Array<UInt, Eigen::Dynamic, Eigen::Dynamic>();

    const Eigen::Index H = layers[0].rows();
    const Eigen::Index W = layers[0].cols();

    Eigen::Array<UInt, Eigen::Dynamic, Eigen::Dynamic> mask(H, W);
    mask.setZero();

    for (std::size_t n = 0; n < layers.size(); ++n) {
        const UInt bit = UInt(1) << n;         // assume n < bit-width(UInt)
        const auto& L = layers[n];
        for (Eigen::Index j = 0; j < H; ++j)
            for (Eigen::Index i = 0; i < W; ++i)
                if (L(j, i)) mask(j, i) |= bit;
    }
    return mask;
}

// explicit instantiations to keep the template out of headers
template Eigen::Array<std::uint8_t, Eigen::Dynamic, Eigen::Dynamic>
Terrain::buildBitMask<std::uint8_t >(const std::vector<ArrXXb>&) const;

template Eigen::Array<std::uint16_t, Eigen::Dynamic, Eigen::Dynamic>
Terrain::buildBitMask<std::uint16_t>(const std::vector<ArrXXb>&) const;

template Eigen::Array<std::uint32_t, Eigen::Dynamic, Eigen::Dynamic>
Terrain::buildBitMask<std::uint32_t>(const std::vector<ArrXXb>&) const;

template Eigen::Array<std::uint64_t, Eigen::Dynamic, Eigen::Dynamic>
Terrain::buildBitMask<std::uint64_t>(const std::vector<ArrXXb>&) const;













