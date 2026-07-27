#include "fragmentation.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace {

	constexpr Real minimumTargetRatio = Real(0.9); // minimum fraction of the target volume loss that is considered acceptable
constexpr int maximumCutIterations = 32; // maximum number of iterations to find a cut that meets the target volume loss

// Result of a single cut operation, including the resulting geometry and the volume removed by the cut.
struct CutResult
{
	ConvexPolyhedron geometry;     // resulting geometry after the cut
	Real removedVolume{ Real(0) }; // volume removed by the cut
};

// Perform a cut on the given convex polyhedron along the specified normal direction to remove a target volume. 
// The function uses a binary search approach to find the appropriate cut offset that achieves the desired volume loss within a specified tolerance.
CutResult cut_to_volume_loss(
    const ConvexPolyhedron& geometry,
    const Vec3& normal,
    Real targetRemovedVolume)
{
	const Real initialVolume = convex_polyhedron_volume(geometry); // volume before the cut
	const Real targetFraction = targetRemovedVolume / initialVolume; // target fraction of volume to be removed
    const Eigen::Matrix<Real, 1, Eigen::Dynamic> projection =
        normal.transpose() * geometry.vertices;
	Real low = projection.minCoeff(); // the minimum projection value along the normal direction
	Real high = projection.maxCoeff(); // the maximum projection value along the normal direction

    // The outward endpoint retains the complete body and is therefore always
    // a safe candidate below the requested loss.
    CutResult best{ geometry, Real(0) };

	// Perform a binary search to find the cut offset that achieves the target volume loss.
    for (int iteration = 0; iteration < maximumCutIterations; ++iteration) {
		const Real middle = (low + high) * Real(0.5); // midpoint of the current search range
        ConvexPolyhedron candidate =
            clip_convex_polyhedron(geometry, normal, middle);
        const Real candidateVolume = convex_polyhedron_volume(candidate);
        const Real removedVolume = initialVolume - candidateVolume;
        const Real removedFraction = removedVolume / initialVolume;

        if (removedFraction > targetFraction) {
            low = middle;
            continue;
        }

        high = middle;
        best = CutResult{ std::move(candidate), removedVolume };
        if (removedFraction >= minimumTargetRatio * targetFraction) {
            break;
        }
    }
	// Check if the best cut found meets the minimum target ratio of the requested volume loss.
    if (targetFraction > Real(1e-10)
        && best.removedVolume < minimumTargetRatio * targetRemovedVolume) {
        throw std::runtime_error(
            "fragmentation cut did not reach its volume-loss tolerance");
    }
    return best;
}

// Generate a filename for the shape file based on the shape index. 
// The filename is formatted as "shape_XXXXXX.pts", where XXXXXX is a zero-padded six-digit number representing the shape index.
std::string shape_file_name(std::uint32_t shapeIndex)
{
    std::ostringstream name;
    name << "shape_" << std::setw(6) << std::setfill('0') << shapeIndex << ".pts";
    return name.str();
}

} // namespace

std::optional<FragmentationResult> fragment_rock(
    const Rock& rock,
    const std::vector<FragmentationContact>& contacts,
    const FragmentationConfig& config)
{
    Real totalWork = Real(0);
	// Sum the total dissipated normal work from all contacts. Only positive work contributes to fragmentation.
    for (const auto& contact : contacts) {
        totalWork += std::max(contact.dissipatedNormalWork, Real(0));
    }
    if (!(totalWork > Real(0)) || !(rock.mass() > Real(0))) {
        return std::nullopt;
    }

	// Compute the specific energy (work per unit mass) and check if it exceeds the threshold for fragmentation.
    const Real specificEnergy = totalWork / rock.mass();
    if (!(specificEnergy > config.thresholdSpecificEnergy)) {
        return std::nullopt;
    }

	// Compute the normalized specific energy and determine the target fraction of mass to be removed based on the configuration parameters.
    const Real normalized = std::clamp(
        (specificEnergy - config.thresholdSpecificEnergy)
        / (config.halfSpecificEnergy - config.thresholdSpecificEnergy),
        Real(0), Real(1));
    const Real targetRemovedFraction = config.maxRemovedFraction
        * std::pow(normalized, config.exponent);
    if (!(targetRemovedFraction > Real(0))) {
        return std::nullopt;
    }

    const Real initialVolume = convex_polyhedron_volume(rock.polyhedronK());
    ConvexPolyhedron retained = rock.polyhedronK();
    Real removedVolume = Real(0);

    // Each simultaneous contact receives its impact-work-weighted share of
    // the original-volume target. Applying the cuts to the evolving mesh
    // avoids double-counting overlapping chipped regions.
    for (const auto& contact : contacts) {
        if (!(contact.dissipatedNormalWork > Real(0))
            || contact.pointK.squaredNorm() == Real(0)) {
            continue;
        }

        const Real share = contact.dissipatedNormalWork / totalWork;
        const Real allocatedVolume = initialVolume * targetRemovedFraction * share;
        if (!(allocatedVolume > Real(0))) {
            continue;
        }

        CutResult cut = cut_to_volume_loss(
            retained, contact.pointK.normalized(), allocatedVolume);
        removedVolume += cut.removedVolume;
        retained = std::move(cut.geometry);
    }

    if (!(removedVolume > initialVolume * Real(1e-12))) {
        return std::nullopt;
    }

    Eigen::Matrix<Real, 7, 1> identityPose;
    identityPose << 0, 0, 0, 1, 0, 0, 0;
    Rock fragment(std::move(retained), rock.density(), Vec6::Zero(), identityPose);
    const Real actualRemovedFraction = Real(1) - fragment.mass() / rock.mass();
    const Real roundoffTolerance = Real(1e-10);
    if (actualRemovedFraction > targetRemovedFraction + roundoffTolerance) {
        throw std::runtime_error(
            "fragmentation exceeded its target mass-loss fraction");
    }
    const Vec3 centerShift = fragment.inputCenterOfMassK();

    return FragmentationResult{
        std::move(fragment),
        centerShift,
        actualRemovedFraction
    };
}

FragmentationShapeWriter::FragmentationShapeWriter(
    const std::filesystem::path& outputDirectory,
    const std::filesystem::path& shapeOutputDirectory,
    const std::string& caseName)
    : trajectoryDirectory_(outputDirectory / shapeOutputDirectory / caseName),
    indexPath_(trajectoryDirectory_ / "shapes.csv")
{
    std::filesystem::create_directories(trajectoryDirectory_);
    std::ofstream index(indexPath_, std::ios::trunc);
    if (!index) {
        throw std::runtime_error("cannot create fragmentation shape index: " + indexPath_.string());
    }
    index << "shape_index,time,mass,vertex_count,file\n";
}

void FragmentationShapeWriter::write(std::uint32_t shapeIndex, Real time, const Rock& rock)
{
    const std::string fileName = shape_file_name(shapeIndex);
    const std::filesystem::path path = trajectoryDirectory_ / fileName;
    std::ofstream output(path, std::ios::trunc);
    if (!output) {
        throw std::runtime_error("cannot write fragmentation point cloud: " + path.string());
    }
    output << std::setprecision(17);
    for (Eigen::Index i = 0; i < rock.verticesK().cols(); ++i) {
        output << rock.verticesK()(0, i) << ' '
            << rock.verticesK()(1, i) << ' '
            << rock.verticesK()(2, i) << '\n';
    }

    std::ofstream index(indexPath_, std::ios::app);
    if (!index) {
        throw std::runtime_error("cannot append fragmentation shape index: " + indexPath_.string());
    }
    index << shapeIndex << ',' << std::setprecision(17) << time << ',' << rock.mass() << ','
        << rock.nVert() << ',' << fileName << '\n';
}
