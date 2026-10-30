#pragma once

#include "math_alias.hpp"
#include "rock.hpp"

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

// Fragmentation configuration parameters. These are used to determine whether a rock should fragment and how much mass should be removed.
struct FragmentationConfig
{
	Real thresholdSpecificEnergy{ Real(0) }; // minimum specific energy required to trigger fragmentation [J/kg]
	Real halfSpecificEnergy{ Real(0) };      // saturation specific energy; legacy XML name retained [J/kg]
    Real exponent{ Real(2) };                // exponent for the fragmentation mass fraction curve
    Real maxRemovedFraction{ Real(0.5) };   // maximum fraction of mass that can be removed
    std::filesystem::path shapeOutputDirectory{ "fragment_shapes" }; // directory to save fragment shapes
};

// Fragmentation contact information. Each contact contributes to the total dissipated work that may trigger fragmentation.
struct FragmentationContact
{
	Vec3 pointK = Vec3::Zero(); // contact point in the body frame K
	Real dissipatedNormalWork{ Real(0) }; // work dissipated at the contact in the normal direction [J]
};

// Result of a fragmentation operation. Contains the new rock shape, the shift in center of mass, and the fraction of mass removed.
struct FragmentationResult
{
    Rock rock;
    Vec3 centerShiftK = Vec3::Zero(); // shift in the center of mass due to fragmentation
    Real removedMassFraction{ Real(0) }; // fraction of mass removed during fragmentation
};

// Fragment a rock based on the provided contacts and configuration. Returns an optional FragmentationResult if fragmentation occurs, or std::nullopt if no fragmentation is needed.
std::optional<FragmentationResult> fragment_rock(
    const Rock& rock,
    const std::vector<FragmentationContact>& contacts,
    const FragmentationConfig& config);

// Generate a filename for a fragmented rock shape based on its index. The filename is formatted as "shape_XXXXXX.pts" where XXXXXX is the zero-padded shape index.
class FragmentationShapeWriter
{
public:
    FragmentationShapeWriter(const std::filesystem::path& outputDirectory,
        const std::filesystem::path& shapeOutputDirectory,
        const std::string& caseName);

    void write(std::uint32_t shapeIndex, Real time, const Rock& rock);

private:
    std::filesystem::path trajectoryDirectory_;
    std::filesystem::path indexPath_;
};
