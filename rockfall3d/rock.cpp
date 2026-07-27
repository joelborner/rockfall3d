#include "rock.hpp"

#include <utility>

Rock::Rock(MatX verts, Real rho, const Vec6& u0_in,
    const Eigen::Matrix<Real, 7, 1>& q0_in)
    : Rock(make_convex_polyhedron(std::move(verts)), rho, u0_in, q0_in)
{
}

Rock::Rock(ConvexPolyhedron geometry, Real rho, const Vec6& u0_in,
    const Eigen::Matrix<Real, 7, 1>& q0_in)
    : m_geometry(std::move(geometry)),
    m_nVert(static_cast<std::size_t>(m_geometry.vertices.cols())),
    m_density(rho),
    m_u0(u0_in),
    m_q0(q0_in)
{
    compute_mass_inertia_volume();
    build_constant_matrices();
}

Rock Rock::withInitialState(const Vec6& u0_in,
    const Eigen::Matrix<Real, 7, 1>& q0_in) const
{
    Rock copy(*this);
    copy.m_u0 = u0_in;
    copy.m_q0 = q0_in;
    return copy;
}

void Rock::compute_mass_inertia_volume()
{
    const PolyhedronMassProperties properties =
        convex_polyhedron_mass_properties(m_geometry);

    m_inputCenterOfMassK = properties.center;
    m_mass = properties.volume * m_density;
    m_theta = properties.inertiaAtCenter * m_density;
    m_geometry.vertices.colwise() -= properties.center;
}

void Rock::build_constant_matrices()
{
    m_KrSP_t.resize(3, 3 * static_cast<Eigen::Index>(m_nVert));
    for (std::size_t i = 0; i < m_nVert; ++i) {
        m_KrSP_t.block<3, 3>(0, 3 * static_cast<Eigen::Index>(i)) =
            tilde(m_geometry.vertices.col(static_cast<Eigen::Index>(i)));
    }

    m_MM.setZero();
    m_MM.block<3, 3>(0, 0) = Mat3::Identity() * m_mass;
    m_MM.block<3, 3>(3, 3) = m_theta;
}
