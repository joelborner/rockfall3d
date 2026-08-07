#pragma once
#include "jor_solver.hpp"

#include <algorithm>
#include <cmath>

struct SorProjectResult
{
    Vec6 velocity = Vec6::Zero();
    bool converged{ true };
    std::uint32_t iterations{ 0 };
};

/* -------------------------------------------------------------------------
   Successive over-relaxation prox solver using the same caller-owned scratch
   storage as JOR-Prox. The only algorithmic difference is the Gauss-Seidel
   sweep order: each contact update immediately feeds back into the current
   generalized velocity state used by later contacts in the same iteration.
--------------------------------------------------------------------------- */
inline SorProjectResult sor_project(ContactSolverWorkspace& ws,
    std::size_t nN,
    const Mat6& MM_inv, const Vec6& rhs, const Vec6& uA,
    Real mu, Real alpha,
    Real epsilonN, Real epsilonT,
    Real atol, Real rtol, std::uint32_t maxiter)
{
    if (nN == 0) {
        return { uA + MM_inv * rhs, true, 0 };
    }

    const Eigen::Index nNI = static_cast<Eigen::Index>(nN);
    auto WN = ws.WN.leftCols(nNI);
    auto WT1 = ws.WT1.leftCols(nNI);
    auto WT2 = ws.WT2.leftCols(nNI);
    auto MWN = ws.MWN.leftCols(nNI);
    auto MWT1 = ws.MWT1.leftCols(nNI);
    auto MWT2 = ws.MWT2.leftCols(nNI);
    auto gammaNA = ws.gammaNA.head(nNI);
    auto gammaT1A = ws.gammaT1A.head(nNI);
    auto gammaT2A = ws.gammaT2A.head(nNI);
    auto PN = ws.PNsub.head(nNI);
    auto PT = ws.PTsub.topRows(nNI);

    MWN.noalias() = MM_inv * WN;
    MWT1.noalias() = MM_inv * WT1;
    MWT2.noalias() = MM_inv * WT2;

    Real maxDiagN = Real(0);
    Real maxDiagT = Real(0);
    for (Eigen::Index j = 0; j < nNI; ++j) {
        maxDiagN = std::max(maxDiagN, WN.col(j).dot(MWN.col(j)));
        maxDiagT = std::max(maxDiagT, WT1.col(j).dot(MWT1.col(j)));
        maxDiagT = std::max(maxDiagT, WT2.col(j).dot(MWT2.col(j)));
    }

    if (!(maxDiagN > Real(0))) maxDiagN = Real(1);
    if (!(maxDiagT > Real(0))) maxDiagT = Real(1);

    const Real rN = alpha / maxDiagN;
    const Real rT = alpha / maxDiagT;
    const Vec6 base = uA + MM_inv * rhs;
    const bool zeroRestitution = (epsilonN == Real(0) && epsilonT == Real(0));

    if (!zeroRestitution) {
        gammaNA.matrix().noalias() = WN.transpose() * uA;
        gammaT1A.matrix().noalias() = WT1.transpose() * uA;
        gammaT2A.matrix().noalias() = WT2.transpose() * uA;
    }

    Vec6 uE = base
        + MWN * PN.matrix()
        + MWT1 * PT.col(0).matrix()
        + MWT2 * PT.col(1).matrix();

    bool converged = false;
    std::uint32_t iterations = 0;
    for (std::uint32_t it = 0; it < maxiter; ++it)
    {
        iterations = it + 1;
        Real deltaPN2 = Real(0);
        Real deltaPT2 = Real(0);
        Real scalePN2 = Real(0);
        Real scalePT2 = Real(0);

        for (Eigen::Index j = 0; j < nNI; ++j)
        {
            const Real oldPN = PN(j);
            const Real oldPTx = PT(j, 0);
            const Real oldPTy = PT(j, 1);

            const Real gammaN = WN.col(j).dot(uE);
            const Real gammaT1 = WT1.col(j).dot(uE);
            const Real gammaT2 = WT2.col(j).dot(uE);
            const Real XiN = zeroRestitution ? gammaN : (gammaN + epsilonN * gammaNA(j));
            const Real XiT1 = zeroRestitution ? gammaT1 : (gammaT1 + epsilonT * gammaT1A(j));
            const Real XiT2 = zeroRestitution ? gammaT2 : (gammaT2 + epsilonT * gammaT2A(j));

            const Real newPN = std::max(oldPN - rN * XiN, Real(0));
            const Real zx = oldPTx - rT * XiT1;
            const Real zy = oldPTy - rT * XiT2;
            const Real nz = std::hypot(zx, zy);

            Real newPTx = zx;
            Real newPTy = zy;
            if (nz > mu * newPN && nz > Real(0)) {
                const Real scale = mu * newPN / nz;
                newPTx *= scale;
                newPTy *= scale;
            }

            const Real dPN = newPN - oldPN;
            const Real dPTx = newPTx - oldPTx;
            const Real dPTy = newPTy - oldPTy;

            PN(j) = newPN;
            PT(j, 0) = newPTx;
            PT(j, 1) = newPTy;

            if (dPN != Real(0) || dPTx != Real(0) || dPTy != Real(0)) {
                uE += MWN.col(j) * dPN
                    + MWT1.col(j) * dPTx
                    + MWT2.col(j) * dPTy;
            }

            deltaPN2 += dPN * dPN;
            deltaPT2 += dPTx * dPTx + dPTy * dPTy;
            scalePN2 += newPN * newPN;
            scalePT2 += newPTx * newPTx + newPTy * newPTy;
        }

        const Real delta = std::sqrt(deltaPN2) + std::sqrt(deltaPT2);
        const Real scale = std::sqrt(scalePN2) + std::sqrt(scalePT2);
        if (delta <= atol + rtol * scale) {
            converged = true;
            break;
        }
    }

    return { uE, converged, iterations };
}
