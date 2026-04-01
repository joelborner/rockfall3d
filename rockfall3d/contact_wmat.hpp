#pragma once
#include "math_alias.hpp"
#include "rock.hpp"

// Build the contact operators W_N, W_T1, W_T2 for the active contact set.
// Each column corresponds to one contact vertex and maps generalized velocity
// to normal / tangential point velocity at that contact.

inline void make_Wmat(const Rock& rock,
    const std::uint32_t* I_N,
    std::size_t nN,
    const Mat3& AIK,
    const Real* dhdx,
    const Real* dhdy,
    MatX& WN,
    MatX& WT1,
    MatX& WT2)
{
    WN.resize(6, nN);
    WT1.resize(6, nN);
    WT2.resize(6, nN);

    for (std::size_t j = 0; j < nN; ++j)
    {
        const std::size_t v = I_N[j];

        // Contact frame induced by the local terrain slope.
        Vec3 n = Vec3{ -dhdx[j], -dhdy[j], Real(1) }.normalized();
        Vec3 t1 = Vec3{ Real(1),  Real(0),  dhdx[j] }.normalized();
        Vec3 t2 = n.cross(t1);

        // Reuse the precomputed skew matrix block of vertex v.
        const Eigen::Block<const MatX, 3, 3> Sk = rock.skewConcat().block<3, 3>(0, 3 * v);

        Mat3 C = -(AIK * Sk).transpose();
        WN.topRows<3>().col(static_cast<Eigen::Index>(j)) = n;
        WN.bottomRows<3>().col(static_cast<Eigen::Index>(j)) = C * n;
        WT1.topRows<3>().col(static_cast<Eigen::Index>(j)) = t1;
        WT1.bottomRows<3>().col(static_cast<Eigen::Index>(j)) = C * t1;
        WT2.topRows<3>().col(static_cast<Eigen::Index>(j)) = t2;
        WT2.bottomRows<3>().col(static_cast<Eigen::Index>(j)) = C * t2;
    }
}
