#include "terrain.hpp"
#include "terrain_io.hpp"

#include <stdexcept>

void Terrain::addRoughness(const std::string& classRasterPath,
    const std::vector<Real>& gridSizes,
    StepFunc stepFunc,
    std::uint64_t seed)
{
    if (gridSizes.empty()) {
        throw std::runtime_error("roughness requires at least one grid size");
    }

    Terrain cls = load_terrain(classRasterPath); // load roughness class raster

    roughnessCellClass = Eigen::Array<std::uint16_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(
        cls.z.rows(), cls.z.cols()); // allocate internal class-ID grid
    roughnessGridSizes = gridSizes;
    roughnessClassParams.resize(gridSizes.size());
    for (std::size_t i = 0; i < gridSizes.size(); ++i) {
        // Precompute values that are constant for each roughness class.
        const Real gridSize = gridSizes[i];
        const Real invGridSize = gridSize > Real(0) ? Real(1) / gridSize : Real(0);
        const Real taperWidth = Real(0.4) * gridSize;
        const Real invTaperWidth = taperWidth > Real(0) ? Real(1) / taperWidth : Real(0);
        roughnessClassParams[i] = RoughnessClassParams{
            gridSize,
            invGridSize,
            taperWidth,
            invTaperWidth,
            stepFunc
        };
    }

    // Convert raster values to integer class IDs for faster access:
    // 0 -> no roughness, 1 -> first grid size, 2 -> second grid size, ...
    for (Eigen::Index j = 0; j < cls.z.rows(); ++j) {
        for (Eigen::Index i = 0; i < cls.z.cols(); ++i) {
            const Real raw = cls.z(j, i);
            const auto idx = static_cast<std::int64_t>(std::llround(raw));
            if (idx >= 1 && idx <= static_cast<std::int64_t>(gridSizes.size())) {
                roughnessCellClass(j, i) = static_cast<std::uint16_t>(idx);
            }
        }
    }

    roughnessCellSize = cls.cell;
    invRoughnessCellSize = roughnessCellSize > Real(0) ? Real(1) / roughnessCellSize : Real(0);
    roughnessXll = cls.xll;
    roughnessYll = cls.yll;
    roughnessNx = static_cast<std::uint32_t>(cls.nx);
    roughnessNy = static_cast<std::uint32_t>(cls.ny);
    roughnessSeed = seed;
}
