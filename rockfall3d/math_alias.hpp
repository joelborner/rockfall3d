// Common numeric aliases and compact Eigen typedefs used across the rockfall
// simulation codebase
#pragma once
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <cstdint>

// Select scalar type (double by default, can be overridden at compile-time)
#ifndef ROCKFALL_REAL_TYPE
#define ROCKFALL_REAL_TYPE double            // set to float if desired
#endif

using Real = ROCKFALL_REAL_TYPE;


// Fixed-size matrices / vectors

template <int R, int C = 1>
using Mat = Eigen::Matrix<Real, R, C>;

using Vec3 = Mat<3>;                    // 3 - vector
using Vec6 = Mat<6>;                    // 6 - vector
using Mat3 = Mat<3, 3>;                 // 3 x 3 matrix
using Mat6 = Mat<6, 6>;                 // 6 x 6 matrix
using Quat = Eigen::Quaternion<Real>;   // Quaternion

// Dynamic-size shorthands

using MatX = Eigen::Matrix<Real,
    Eigen::Dynamic, Eigen::Dynamic,
    Eigen::RowMajor>;             // row-major matches DEM storage
using ArrX = Eigen::Array<Real, Eigen::Dynamic, 1>;
using IdxV = Eigen::Array<std::uint32_t, Eigen::Dynamic, 1>;
using ArrXXb = Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>;

// Utility functions
// Tilde operator - returns the skew-symmetric matrix of a 3 - vector.
inline Mat3 tilde(const Vec3& k) noexcept
{
    return (Mat3() << 0, -k.z(), k.y(),
        k.z(), 0, -k.x(),
        -k.y(), k.x(), 0).finished();
}
