#pragma once
#include "math_alias.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

struct ContactSolverWorkspace
{
    // Contact operator blocks and solver iterates for the currently active
    // contact set. The arrays are resized only when the active contact count
    // exceeds previous capacity.
    MatX WN;
    MatX WT1;
    MatX WT2;
    MatX MWN;
    MatX MWT1;
    MatX MWT2;
    ArrX gammaNA;
    ArrX gammaT1A;
    ArrX gammaT2A;
    ArrX PNsub;
    MatX PTsub;

    void ensureContactCapacity(std::size_t n)
    {
        const Eigen::Index nI = static_cast<Eigen::Index>(n);
        if (WN.cols() < nI) {
            WN.resize(6, nI);
            WT1.resize(6, nI);
            WT2.resize(6, nI);
            MWN.resize(6, nI);
            MWT1.resize(6, nI);
            MWT2.resize(6, nI);
        }
        if (PNsub.size() < nI) {
            PNsub.resize(nI);
        }
        if (gammaNA.size() < nI) {
            gammaNA.resize(nI);
            gammaT1A.resize(nI);
            gammaT2A.resize(nI);
        }
        if (PTsub.rows() < nI) {
            PTsub.resize(nI, 2);
        }
    }
};

struct JorProjectResult
{
    Vec6 velocity = Vec6::Zero();
    bool converged{ true };
};

/* -------------------------------------------------------------------------
   Jacobi over-relaxation prox solver using caller-owned scratch storage.
   Inputs : WN, WT1, WT2 already filled in workspace, nN active contacts.
--------------------------------------------------------------------------- */
inline JorProjectResult jor_project(ContactSolverWorkspace& ws,
    std::size_t nN,
    const Mat6& mobility, const Vec6& freeVelocity, const Vec6& uA,
    Real mu, Real alpha,
    Real epsilonN, Real epsilonT,
    Real atol, Real rtol, std::uint32_t maxiter)
{
    if (nN == 0) {
        return { freeVelocity, true };
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

    // Precompute the generalized mobility times W once per solve so the iteration only sees
    // matrix-vector products in the current contact subspace.
    MWN.noalias() = mobility * WN;
    MWT1.noalias() = mobility * WT1;
    MWT2.noalias() = mobility * WT2;

    // The diagonal maxima define the Jacobi relaxation scales for normal and
    // tangential updates. These are equivalent to the diagonal entries of
    // W^T M^{-1} W but avoid materializing the full dense matrices.
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
    const Vec6& base = freeVelocity;
    const bool zeroRestitution = (epsilonN == Real(0) && epsilonT == Real(0));

    if (!zeroRestitution) {
        gammaNA.matrix().noalias() = WN.transpose() * uA;
        gammaT1A.matrix().noalias() = WT1.transpose() * uA;
        gammaT2A.matrix().noalias() = WT2.transpose() * uA;
    }

    Vec6 uE = base;
    bool converged = false;
    for (std::uint32_t it = 0; it < maxiter; ++it)
    {
        // Evaluate the current generalized velocity using the latest impulse
        // iterate.
        uE = base
            + MWN * PN.matrix()
            + MWT1 * PT.col(0).matrix()
            + MWT2 * PT.col(1).matrix();

        // Convergence metric intentionally mirrors the original formulation:
        // normal update norm + tangential update norm.
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

            PN(j) = newPN;
            PT(j, 0) = newPTx;
            PT(j, 1) = newPTy;

            const Real dPN = newPN - oldPN;
            const Real dPTx = newPTx - oldPTx;
            const Real dPTy = newPTy - oldPTy;
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

    return { uE, converged };
}
