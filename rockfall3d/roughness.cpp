#include "terrain.hpp"
#include "terrain_io.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

void Terrain::addRoughness(const std::string& rasterPath,
    StepFunc stepFunc,
    std::uint64_t seed)
{
    Terrain rough = load_terrain(rasterPath); // load raster with grid sizes in map units

    roughnessCellClass = Eigen::Array<std::uint16_t, Eigen::Dynamic, Eigen::Dynamic>::Zero(
        rough.z.rows(), rough.z.cols());
    roughnessClassParams.clear();

    // Convert physical grid sizes to compact class IDs for the hot loop:
    // non-finite values and values <= 0 mean no roughness.
    std::map<Real, std::uint16_t> classByGridSize;
    for (Eigen::Index j = 0; j < rough.z.rows(); ++j) {
        for (Eigen::Index i = 0; i < rough.z.cols(); ++i) {
            const Real gridSize = rough.z(j, i);
            if (!std::isfinite(gridSize) || gridSize <= Real(0)) {
                continue;
            }

            auto it = classByGridSize.find(gridSize);
            if (it == classByGridSize.end()) {
                if (roughnessClassParams.size() >= std::numeric_limits<std::uint16_t>::max()) {
                    throw std::runtime_error("roughness raster contains too many distinct grid sizes");
                }

                const auto classId = static_cast<std::uint16_t>(roughnessClassParams.size() + 1);
                classByGridSize.emplace(gridSize, classId);

                const Real invGridSize = Real(1) / gridSize;
                const Real taperWidth = Real(0.4) * gridSize;
                const Real invTaperWidth = taperWidth > Real(0) ? Real(1) / taperWidth : Real(0);
                roughnessClassParams.push_back(RoughnessClassParams{
                    gridSize,
                    invGridSize,
                    taperWidth,
                    invTaperWidth,
                    stepFunc
                });
                roughnessCellClass(j, i) = classId;
            }
            else {
                roughnessCellClass(j, i) = it->second;
            }
        }
    }

    roughnessCellSize = rough.cell;
    invRoughnessCellSize = roughnessCellSize > Real(0) ? Real(1) / roughnessCellSize : Real(0);
    roughnessXll = rough.xll;
    roughnessYll = rough.yll;
    roughnessNx = static_cast<std::uint32_t>(rough.nx);
    roughnessNy = static_cast<std::uint32_t>(rough.ny);
    roughnessSeed = seed;
}
