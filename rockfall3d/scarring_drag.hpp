#pragma once

#include "math_alias.hpp"
#include "rock.hpp"
#include "terrain.hpp"
#include "time_stepper.hpp"

struct ScarringDragResult
{
    Vec3 force = Vec3::Zero();
    Vec3 torque = Vec3::Zero();
    Real projectedArea{ Real(0) };
};

// Werner-type maximum depth of the deformable layer for one impact episode.
Real calculate_scar_layer_depth(Real ME, Real mass, Real verticalImpactSpeed) noexcept;

// Evaluate one of the three drag models on the portion of the convex rock
// lying below the original terrain surface. Force and torque are in world axes.
ScarringDragResult evaluate_scarring_drag(
    const ScarringConfig& config,
    const Rock& rock,
    const Mat3& AIK,
    const Vec3& center,
    const Vec6& velocity,
    const TerrainQueryWorkspace& terrainQuery,
    MatX& worldVertices,
    std::vector<Vec3>& areaVectors,
    std::vector<Vec3>& centroids);
