#pragma once
#include "math_alias.hpp"
#include <vector>

// Rigid geometric/mechanical description of one rock shape together with
// an initial state used to start a trajectory. Geometry is stored in the body
// frame K with the center of mass at the origin.
class Rock
{
public:
    // verts is an arbitrary 3D point cloud describing the rock surface.
    // The constructor builds a convex hull, derives mass/inertia, and shifts
    // the vertices so the center of mass becomes the body-frame origin.
    explicit Rock(MatX verts, Real density,
        const Vec6& u0_in,
        const Eigen::Matrix<Real, 7, 1>& q0_in);

    // Reuse the same geometric data but replace the initial state.
    Rock withInitialState(const Vec6& u0_in,
        const Eigen::Matrix<Real, 7, 1>& q0_in) const;

    // Geometric data
    std::size_t   nVert() const noexcept { return m_nVert; }
    const MatX& verticesK() const noexcept { return m_KrSP; }
    Real          mass() const noexcept { return m_mass; }
    Real          equivalentRadius() const noexcept { return m_equivalentRadius; }
    const Mat3& inertiaK() const noexcept { return m_theta; }
    const Mat6& massMatrix() const noexcept { return m_MM; }
    const MatX& skewConcat() const noexcept { return m_KrSP_t; }

    // initial state
    const Vec6& u0() const noexcept { return m_u0; }
    const Eigen::Matrix<Real, 7, 1>& q0() const noexcept { return m_q0; }

private:
    MatX m_KrSP;                      // 3 x n body-frame vertices, COM at origin
    MatX m_KrSP_t;                    // concatenated skew matrices [tilde(r1) ... tilde(rn)]
    std::size_t m_nVert{};
    Real        m_density{};
    Real        m_mass{};
    Real        m_equivalentRadius{};    // radius of a sphere with the same volume
    Mat3        m_theta{};            // inertia tensor in body frame K
    Mat6        m_MM{};               // block-diagonal generalized mass matrix

    Vec6                      m_u0;
    Eigen::Matrix<Real, 7, 1> m_q0;

    // Build convex-hull-based volume, mass, center of mass, and inertia.
    void compute_mass_inertia_volume();

    // Precompute matrices reused throughout contact assembly and time stepping.
    void build_constant_matrices();
};
