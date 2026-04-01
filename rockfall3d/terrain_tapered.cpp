#include "terrain.hpp"
#include "roughness.hpp"
#include "terrain_io.hpp"

#include <stdexcept>

void Terrain::addTaperedCellRoughness(const std::string& classRasterPath,
    const std::vector<Real>& gridSizes,
    std::uint64_t seed)
{
    if (gridSizes.empty()) {
        throw std::runtime_error("tapered-cell roughness requires at least one grid size");
    }

    Terrain cls = load_terrain(classRasterPath); // load roughness raster

    taperedCellClass = Eigen::Array<std::uint16_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(
        cls.z.rows(), cls.z.cols()); // allocate internal class-ID grid
    taperedGridSizes = gridSizes;
    taperedClassParams.resize(gridSizes.size());
    for (std::size_t i = 0; i < gridSizes.size(); ++i) {
        // loops through all roughness magnitudes and precomputes static values
        const Real gridSize = gridSizes[i];
        const Real invGridSize = gridSize > Real(0) ? Real(1) / gridSize : Real(0);
        const Real taperWidth = Real(0.4) * gridSize;
        const Real invTaperWidth = taperWidth > Real(0) ? Real(1) / taperWidth : Real(0);
        taperedClassParams[i] = TaperedClassParams{
            gridSize,
            invGridSize,
            taperWidth,
            invTaperWidth
        };
    }

    // convert raster values to integer class IDs for faster access
    // 0 -> no roughness
    // 1 -> first grid size
    // 2 -> second grid size
    for (Eigen::Index j = 0; j < cls.z.rows(); ++j) {
        for (Eigen::Index i = 0; i < cls.z.cols(); ++i) {
            const Real raw = cls.z(j, i);
            const auto idx = static_cast<std::int64_t>(std::llround(raw));
            if (idx >= 1 && idx <= static_cast<std::int64_t>(gridSizes.size())) {
                taperedCellClass(j, i) = static_cast<std::uint16_t>(idx);
            }
        }
    }

    taperedCellSize = cls.cell;
    invTaperedCellSize = taperedCellSize > Real(0) ? Real(1) / taperedCellSize : Real(0);
    taperedXll = cls.xll;
    taperedYll = cls.yll;
    taperedNx = static_cast<std::uint32_t>(cls.nx);
    taperedNy = static_cast<std::uint32_t>(cls.ny);
    taperedSeed = seed;

    // ensures legacy roughness model is disabled and does not cause trouble
    roughness.clear();
    roughIndices.resize(0, 0);
    minrough = Real(0);
    roughNx = roughNy = 0;
    shapes = false;
    hasRoughness = true;
    roughnessBackend = RoughnessBackend::TaperedCells;
}
