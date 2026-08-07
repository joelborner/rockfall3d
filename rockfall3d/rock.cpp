#include "rock.hpp"
#include "quickhull.hpp"
#include <Eigen/Geometry>


Rock::Rock(MatX verts, Real rho, const Vec6& u0_in,
    const Eigen::Matrix<Real, 7, 1>& q0_in)
    : m_KrSP(std::move(verts)),
    m_nVert(m_KrSP.cols()),
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

// Build convex-hull-based volume, mass, center of mass, and inertia.
void Rock::compute_mass_inertia_volume()
{

    // Build convex hull with quick_hull (header-only external library) Copyright(c) 2014 - 2015, Anatoliy V.Tomilov)
    using Point = std::array<Real, 3>;

    // quickhull works on a STL-like point container. 
    // The input vertex cloud is copied into that format first.
    std::vector<Point> pts(m_nVert);
    for (std::size_t i = 0; i < m_nVert; ++i) {
        const auto p = m_KrSP.col(i);
        pts[i] = { p.x(), p.y(), p.z() };
    }

    constexpr Real eps = std::numeric_limits<Real>::epsilon();
    quick_hull<std::vector<Point>::iterator, Real> qh{ 3, eps };

    qh.add_points(pts.begin(), pts.end());

    // The affine basis contains extreme points suitable for seeding the hull.
    // In 3D we need at least 4 non-coplanar points to build an initial simplex.
    auto basis = qh.get_affine_basis();
    if (basis.size() < 4)
        throw std::runtime_error("point set is not 3-D");

    // The initial simplex is the tetrahedron from which QuickHull starts
    // growing the final convex hull. After create_convex_hull(), facets_
    // stores the triangular boundary facets of that hull.
    qh.create_initial_simplex(basis.begin(), std::prev(basis.end()));
    qh.create_convex_hull();

    // 2.  Copy facet indices and enforce a consistent outward orientation
    m_triangleFaces.clear();
    m_triangleFaces.reserve(qh.facets_.size() * 3);

    // The mean of the input vertices is a convex combination of hull vertices
    // and therefore lies inside the convex hull. It can be used to decide
    // whether a facet normal points outward or inward.
    const Vec3 hullCtr = m_KrSP.rowwise().mean();

    for (const auto& f : qh.facets_) {
        std::array<std::size_t, 3> tri;
        for (int k = 0; k < 3; ++k)
            tri[k] = std::distance(pts.begin(), f.vertices_[k]);

        // orientation test
        const Vec3  A = m_KrSP.col(tri[0]),
            B = m_KrSP.col(tri[1]),
            C = m_KrSP.col(tri[2]);
        const Vec3  n = (B - A).cross(C - A);           // triangle normal for current winding
        if (n.dot(A - hullCtr) < Real(0))
            std::swap(tri[1], tri[2]);                  // flip inward-facing triangles

        for (const std::size_t index : tri)
            m_triangleFaces.push_back(static_cast<std::uint32_t>(index));
    }


    // 3.  Accumulate signed volume, first moment, and inertia
    Real  V = 0.0;
    Vec3  S1 = Vec3::Zero();
    Real  Ixx = 0, Iyy = 0, Izz = 0,
        Ixy = 0, Iyz = 0, Izx = 0;

    // Each oriented surface triangle together with the origin defines one
    // signed tetrahedron. Summing these tetrahedra gives the total polyhedron
    // volume and inertia. The closed-form monomials below are the standard
    // tetrahedral formulas written out explicitly.
    for (std::size_t i = 0; i < m_triangleFaces.size(); i += 3)
    {
        const Vec3& A = m_KrSP.col(m_triangleFaces[i]);
        const Vec3& B = m_KrSP.col(m_triangleFaces[i + 1]);
        const Vec3& C = m_KrSP.col(m_triangleFaces[i + 2]);

        const Real d = A.dot(B.cross(C));   // 6 x Vt
        const Real Vt = d / Real(6);

        const Vec3 Ct = (A + B + C) / Real(4);
        V += Vt;
        S1 += Vt * Ct;

        const Real Ax = A.x(), Ay = A.y(), Az = A.z();
        const Real Bx = B.x(), By = B.y(), Bz = B.z();
        const Real Cx = C.x(), Cy = C.y(), Cz = C.z();

        const Real
            xx = Ax * Ax + Bx * Bx + Cx * Cx + Ax * Bx + Bx * Cx + Cx * Ax,
            yy = Ay * Ay + By * By + Cy * Cy + Ay * By + By * Cy + Cy * Ay,
            zz = Az * Az + Bz * Bz + Cz * Cz + Az * Bz + Bz * Cz + Cz * Az,
            xy = 2 * (Ax * Ay + Bx * By + Cx * Cy) +
            (Ax * By + Bx * Ay + Bx * Cy + Cx * By + Cx * Ay + Ax * Cy),
            yz = 2 * (Ay * Az + By * Bz + Cy * Cz) +
            (Ay * Bz + By * Az + By * Cz + Cy * Bz + Cy * Az + Ay * Cz),
            zx = 2 * (Az * Ax + Bz * Bx + Cz * Cx) +
            (Az * Bx + Bz * Ax + Bz * Cx + Cz * Bx + Cz * Ax + Az * Cx);

        const Real kDiag = Vt / Real(10);
        const Real kOff = Vt / Real(20);

        Ixx += kDiag * (yy + zz);           // (y^2+z^2)
        Iyy += kDiag * (zz + xx);           // (z^2+x^2)
        Izz += kDiag * (xx + yy);           // (x^2+y^2)

        Ixy += kOff * xy;                 // xy
        Iyz += kOff * yz;                 // yz
        Izx += kOff * zx;                 // zx
    }

    // 4.  Make the aggregate quantities positive if the hull winding
    if (V < Real(0)) {
        V = -V;  S1 = -S1;
        Ixx = -Ixx; Iyy = -Iyy; Izz = -Izz;
        Ixy = -Ixy; Iyz = -Iyz; Izx = -Izx;
    }      

    // 5.  Shift inertia from the origin to the center of mass
    const Vec3 COM = S1 / V;       // vector from origin to COM
    m_mass = V * m_density;

    const Real cx = COM.x(), cy = COM.y(), cz = COM.z();

    Ixx -= V * (cy * cy + cz * cz);
    Iyy -= V * (cz * cz + cx * cx);
    Izz -= V * (cx * cx + cy * cy);

    // Parallel-axis correction for the product-of-inertia terms.
    Ixy -= V * cx * cy;
    Iyz -= V * cy * cz;
    Izx -= V * cz * cx;

    // 6.  Apply material density and assemble the inertia tensor
    const Real rho = m_density;
    Ixx *= rho;  Iyy *= rho;  Izz *= rho;
    Ixy *= rho;  Iyz *= rho;  Izx *= rho;

    m_theta <<
        Ixx, -Ixy, -Izx,
        -Ixy, Iyy, -Iyz,
        -Izx, -Iyz, Izz;

    // 7.  Translate the point cloud so the body-frame origin is the COM
    for (std::size_t i = 0; i < m_nVert; ++i)
        m_KrSP.col(i) -= COM;
}

void Rock::build_constant_matrices()
{
    // Concatenate the skew matrices of all body-frame vertices once. Contact
    // assembly reuses these blocks repeatedly when building W matrices.
    m_KrSP_t.resize(3, 3 * m_nVert);
    for (std::size_t i = 0; i < m_nVert; ++i)
        m_KrSP_t.block<3, 3>(0, 3 * i) = tilde(m_KrSP.col(i));

    // Generalized mass matrix [mI, 0; 0, Theta].
    m_MM.setZero();
    m_MM.block<3, 3>(0, 0) = Mat3::Identity() * m_mass;
    m_MM.block<3, 3>(3, 3) = m_theta;
}

