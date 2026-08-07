#include "scarring_drag.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace {

inline Vec3 edge_intersection(const Vec3& a, const Vec3& b, Real ga, Real gb) noexcept
{
    const Real denominator = ga - gb;
    const Real t = std::abs(denominator) > Real(1e-15) ? ga / denominator : Real(0.5);
    return a + std::clamp(t, Real(0), Real(1)) * (b - a);
}

void append_clipped_triangle(
    const std::array<Vec3, 3>& points,
    const std::array<Real, 3>& gaps,
    std::vector<Vec3>& areaVectors,
    std::vector<Vec3>& centroids,
    Vec3& capAreaVector,
    Vec3& capCentroidSum,
    std::uint32_t& capPointCount)
{
    std::array<Vec3, 5> clippedPoints;
    std::array<Real, 5> clippedGaps{};
    std::uint32_t clippedCount = 0;

    for (std::uint32_t i = 0; i < 3; ++i) {
        const std::uint32_t j = (i + 1) % 3;
        const bool insideA = gaps[i] <= Real(0);
        const bool insideB = gaps[j] <= Real(0);

        if (insideA) {
            clippedPoints[clippedCount] = points[i];
            clippedGaps[clippedCount] = gaps[i];
            ++clippedCount;
        }
        if (insideA != insideB) {
            const Vec3 intersection = edge_intersection(points[i], points[j], gaps[i], gaps[j]);
            clippedPoints[clippedCount] = intersection;
            clippedGaps[clippedCount] = Real(0);
            ++clippedCount;
            capCentroidSum += intersection;
            ++capPointCount;
        }
    }

    if (clippedCount < 3) {
        return;
    }

    for (std::uint32_t i = 1; i + 1 < clippedCount; ++i) {
        const Vec3& a = clippedPoints[0];
        const Vec3& b = clippedPoints[i];
        const Vec3& c = clippedPoints[i + 1];
        const Vec3 areaVector = Real(0.5) * (b - a).cross(c - a);
        if (areaVector.squaredNorm() <= Real(1e-24)) {
            continue;
        }
        areaVectors.push_back(areaVector);
        centroids.push_back((a + b + c) / Real(3));
        capAreaVector -= areaVector;
    }
}

inline Vec3 point_velocity(const Vec3& translationalVelocity,
    const Vec3& angularVelocityWorld,
    const Vec3& point,
    const Vec3& center) noexcept
{
    return translationalVelocity + angularVelocityWorld.cross(point - center);
}

} // namespace

Real calculate_scar_layer_depth(Real ME, Real mass, Real verticalImpactSpeed) noexcept
{
    if (!(ME > Real(0)) || !(mass > Real(0)) || !(verticalImpactSpeed > Real(0))) {
        return Real(0);
    }
    return Real(0.1)
        * std::cbrt(mass)
        * std::pow(ME, Real(-0.4))
        * std::pow(verticalImpactSpeed, Real(0.8));
}

ScarringDragResult evaluate_scarring_drag(
    const ScarringConfig& config,
    const Rock& rock,
    const Mat3& AIK,
    const Vec3& center,
    const Vec6& velocity,
    const TerrainQueryWorkspace& terrainQuery,
    MatX& worldVertices,
    std::vector<Vec3>& areaVectors,
    std::vector<Vec3>& centroids)
{
    ScarringDragResult result;
    if (!config.enabled() || !(config.dragCoefficient > Real(0))
        || !(config.density > Real(0))) {
        return result;
    }

    const std::size_t nVert = rock.nVert();
    worldVertices.resize(3, static_cast<Eigen::Index>(nVert));
    for (std::size_t i = 0; i < nVert; ++i) {
        worldVertices.col(static_cast<Eigen::Index>(i)) =
            center + AIK * rock.verticesK().col(static_cast<Eigen::Index>(i));
    }

    areaVectors.clear();
    centroids.clear();
    const std::size_t reserveCount = 2 * (rock.triangleFaces().size() / 3) + 1;
    if (areaVectors.capacity() < reserveCount) {
        areaVectors.reserve(reserveCount);
        centroids.reserve(reserveCount);
    }
    Vec3 capAreaVector = Vec3::Zero();
    Vec3 capCentroidSum = Vec3::Zero();
    std::uint32_t capPointCount = 0;

    const auto& faces = rock.triangleFaces();
    for (std::size_t i = 0; i + 2 < faces.size(); i += 3) {
        const std::array<std::uint32_t, 3> indices{ faces[i], faces[i + 1], faces[i + 2] };
        const std::array<Vec3, 3> points{
            worldVertices.col(indices[0]),
            worldVertices.col(indices[1]),
            worldVertices.col(indices[2])
        };
        const std::array<Real, 3> gaps{
            terrainQuery.gN(indices[0]),
            terrainQuery.gN(indices[1]),
            terrainQuery.gN(indices[2])
        };
        append_clipped_triangle(points, gaps, areaVectors, centroids,
            capAreaVector, capCentroidSum, capPointCount);
    }

    if (capPointCount > 0 && capAreaVector.squaredNorm() > Real(1e-24)) {
        areaVectors.push_back(capAreaVector);
        centroids.push_back(capCentroidSum / Real(capPointCount));
    }
    if (areaVectors.empty()) {
        return result;
    }

    const Vec3 vCom = velocity.head<3>();
    const Vec3 omegaWorld = AIK * velocity.tail<3>();
    const Real coefficient = Real(0.5) * config.dragCoefficient * config.density;
    constexpr Real speedEpsilon = Real(1e-12);

    if (config.mode == ScarringMode::Distributed) {
        for (std::size_t i = 0; i < areaVectors.size(); ++i) {
            const Vec3 localVelocity =
                point_velocity(vCom, omegaWorld, centroids[i], center);
            const Real speed = localVelocity.norm();
            if (speed <= speedEpsilon) {
                continue;
            }
            const Real projectedArea =
                std::max(Real(0), areaVectors[i].dot(localVelocity / speed));
            if (!(projectedArea > Real(0))) {
                continue;
            }
            const Vec3 force = -coefficient * projectedArea * speed * localVelocity;
            result.force += force;
            result.torque += (centroids[i] - center).cross(force);
            result.projectedArea += projectedArea;
        }
        return result;
    }

    const Real speed = vCom.norm();
    if (speed <= speedEpsilon) {
        return result;
    }
    const Vec3 direction = vCom / speed;
    Vec3 weightedCentroid = Vec3::Zero();
    for (std::size_t i = 0; i < areaVectors.size(); ++i) {
        const Real area = std::max(Real(0), areaVectors[i].dot(direction));
        result.projectedArea += area;
        weightedCentroid += area * centroids[i];
    }
    if (!(result.projectedArea > Real(0))) {
        return result;
    }

    if (config.mode == ScarringMode::CenterOfMass) {
        result.force = -coefficient * result.projectedArea * speed * vCom;
        return result;
    }

    const Vec3 applicationPoint = weightedCentroid / result.projectedArea;
    const Vec3 localVelocity =
        point_velocity(vCom, omegaWorld, applicationPoint, center);
    const Real localSpeed = localVelocity.norm();
    if (localSpeed <= speedEpsilon) {
        return result;
    }
    result.force = -coefficient * result.projectedArea * localSpeed * localVelocity;
    result.torque = (applicationPoint - center).cross(result.force);
    return result;
}
