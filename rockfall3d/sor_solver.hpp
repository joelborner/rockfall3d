#pragma once
#include "jor_solver.hpp"

#include <algorithm>
#include <cmath>

struct SorProjectResult
{
    Vec6 velocity = Vec6::Zero();
    bool converged{ true };
};

/// Successive over-relaxation proximal solver with a compile-time contact law.
///
/// Each contact immediately updates uE (Gauss-Seidel feedback). JOR instead
/// evaluates every contact from a common velocity snapshot.
/// RollingEnabled=false uses [PN, PT1, PT2], with separate normal and
/// tangential-disc projections. RollingEnabled=true adds two moment impulses
/// [PR1, PR2] and projects jointly onto the rolling-friction cone, with the
/// De Saxce correction required for non-associated friction.
///
/// mu is dimensionless; muR is the rolling-resistance length in metres.
/// The caller supplies initialized active operators/impulses and workspace
/// capacity, including rolling storage only for the enabled specialization.
/// The disabled specialization neither accesses nor allocates rolling storage.
///
/// Existing relaxation scales and stopping norms are retained. Convergence
/// measures impulse-iterate change, not an independent contact-law residual.
template <bool RollingEnabled>
inline SorProjectResult sor_project(ContactSolverWorkspace& ws,
    std::size_t nN,
    const Mat6& mobility, const Vec6& freeVelocity, const Vec6& uA,
    Real mu, Real muR, Real alpha,
    Real epsilonN, Real epsilonT,
    Real atol, Real rtol, std::uint32_t maxiter)
{
    if (nN == 0) return { freeVelocity, true };

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

    MWN.noalias() = mobility * WN;
    MWT1.noalias() = mobility * WT1;
    MWT2.noalias() = mobility * WT2;

    Real rN = Real(0), rT = Real(0);
    if constexpr (RollingEnabled) {
        auto WR1 = ws.WR1.leftCols(nNI);
        auto WR2 = ws.WR2.leftCols(nNI);
        auto MWR1 = ws.MWR1.leftCols(nNI);
        auto MWR2 = ws.MWR2.leftCols(nNI);
        MWR1.noalias() = mobility * WR1;
        MWR2.noalias() = mobility * WR2;

        // Preserve the rolling SOR path's per-contact diagonal scale.
        for (Eigen::Index j = 0; j < nNI; ++j) {
            Real maxDiag = std::max({
                WN.col(j).dot(MWN.col(j)),
                WT1.col(j).dot(MWT1.col(j)),
                WT2.col(j).dot(MWT2.col(j)),
                WR1.col(j).dot(MWR1.col(j)),
                WR2.col(j).dot(MWR2.col(j))
                });
            if (!(maxDiag > Real(0))) maxDiag = Real(1);
            ws.rollingStep(j) = alpha / maxDiag;
        }
    }
    else {
        // Preserve the Coulomb path's separate global normal/tangent scales.
        Real maxDiagN = Real(0), maxDiagT = Real(0);
        for (Eigen::Index j = 0; j < nNI; ++j) {
            maxDiagN = std::max(maxDiagN, WN.col(j).dot(MWN.col(j)));
            maxDiagT = std::max(maxDiagT, WT1.col(j).dot(MWT1.col(j)));
            maxDiagT = std::max(maxDiagT, WT2.col(j).dot(MWT2.col(j)));
        }
        if (!(maxDiagN > Real(0))) maxDiagN = Real(1);
        if (!(maxDiagT > Real(0))) maxDiagT = Real(1);
        rN = alpha / maxDiagN;
        rT = alpha / maxDiagT;
    }

    const bool zeroRestitution = epsilonN == Real(0) && epsilonT == Real(0);
    if (!zeroRestitution) {
        gammaNA.matrix().noalias() = WN.transpose() * uA;
        gammaT1A.matrix().noalias() = WT1.transpose() * uA;
        gammaT2A.matrix().noalias() = WT2.transpose() * uA;
    }

    const auto velocityFromImpulses = [&]() -> Vec6 {
        if constexpr (RollingEnabled) {
            return freeVelocity
                + MWN * PN.matrix()
                + MWT1 * PT.col(0).matrix()
                + MWT2 * PT.col(1).matrix()
                + ws.MWR1.leftCols(nNI) * ws.PRsub.topRows(nNI).col(0).matrix()
                + ws.MWR2.leftCols(nNI) * ws.PRsub.topRows(nNI).col(1).matrix();
        }
        else {
            return freeVelocity
                + MWN * PN.matrix()
                + MWT1 * PT.col(0).matrix()
                + MWT2 * PT.col(1).matrix();
        }
    };

    Vec6 uE = velocityFromImpulses();
    bool converged = false;
    for (std::uint32_t it = 0; it < maxiter; ++it) {
        Real deltaPN2 = Real(0), deltaPT2 = Real(0);
        Real scalePN2 = Real(0), scalePT2 = Real(0);
        Real delta2 = Real(0), scale2 = Real(0);
        for (Eigen::Index j = 0; j < nNI; ++j) {
            const Real oldPN = PN(j);
            const Real oldPT1 = PT(j, 0), oldPT2 = PT(j, 1);
            Real oldPR1 = Real(0), oldPR2 = Real(0);
            const Real gammaN = WN.col(j).dot(uE);
            const Real gammaT1 = WT1.col(j).dot(uE);
            const Real gammaT2 = WT2.col(j).dot(uE);
            const Real xiN = zeroRestitution ? gammaN : gammaN + epsilonN * gammaNA(j);
            const Real xiT1 = zeroRestitution ? gammaT1 : gammaT1 + epsilonT * gammaT1A(j);
            const Real xiT2 = zeroRestitution ? gammaT2 : gammaT2 + epsilonT * gammaT2A(j);

            Real newPN, newPT1, newPT2;
            Real newPR1 = Real(0), newPR2 = Real(0);
            if constexpr (RollingEnabled) {
                const Real r = ws.rollingStep(j);
                oldPR1 = ws.PRsub(j, 0);
                oldPR2 = ws.PRsub(j, 1);
                const Real gammaR1 = ws.WR1.col(j).dot(uE);
                const Real gammaR2 = ws.WR2.col(j).dot(uE);
                // Joint projection needs this shift to avoid normal dilatancy.
                newPN = oldPN - r * (xiN
                    + mu * std::hypot(xiT1, xiT2)
                    + muR * std::hypot(gammaR1, gammaR2));
                newPT1 = oldPT1 - r * xiT1;
                newPT2 = oldPT2 - r * xiT2;
                newPR1 = oldPR1 - r * gammaR1;
                newPR2 = oldPR2 - r * gammaR2;
                project_rolling_friction_cone(
                    newPN, newPT1, newPT2, newPR1, newPR2, mu, muR);
            }
            else {
                // Separate projections already encode non-associated Coulomb
                // friction: do NOT add the rolling path's De Saxce shift.
                newPN = std::max(oldPN - rN * xiN, Real(0));
                newPT1 = oldPT1 - rT * xiT1;
                newPT2 = oldPT2 - rT * xiT2;
                const Real normT = std::hypot(newPT1, newPT2);
                if (normT > mu * newPN && normT > Real(0)) {
                    const Real scale = mu * newPN / normT;
                    newPT1 *= scale;
                    newPT2 *= scale;
                }
            }

            const Real dPN = newPN - oldPN;
            const Real dPT1 = newPT1 - oldPT1, dPT2 = newPT2 - oldPT2;
            PN(j) = newPN;
            PT(j, 0) = newPT1;
            PT(j, 1) = newPT2;
            if constexpr (RollingEnabled) {
                const Real dPR1 = newPR1 - oldPR1, dPR2 = newPR2 - oldPR2;
                ws.PRsub(j, 0) = newPR1;
                ws.PRsub(j, 1) = newPR2;
                // Immediate feedback is the defining SOR/Gauss-Seidel step.
                if (dPN != Real(0) || dPT1 != Real(0) || dPT2 != Real(0)
                    || dPR1 != Real(0) || dPR2 != Real(0)) {
                    uE += MWN.col(j) * dPN
                        + MWT1.col(j) * dPT1
                        + MWT2.col(j) * dPT2
                        + ws.MWR1.col(j) * dPR1
                        + ws.MWR2.col(j) * dPR2;
                }
                delta2 += dPN * dPN + dPT1 * dPT1 + dPT2 * dPT2
                    + dPR1 * dPR1 + dPR2 * dPR2;
                scale2 += newPN * newPN + newPT1 * newPT1 + newPT2 * newPT2
                    + newPR1 * newPR1 + newPR2 * newPR2;
            }
            else {
                if (dPN != Real(0) || dPT1 != Real(0) || dPT2 != Real(0)) {
                    uE += MWN.col(j) * dPN
                        + MWT1.col(j) * dPT1
                        + MWT2.col(j) * dPT2;
                }
                deltaPN2 += dPN * dPN;
                deltaPT2 += dPT1 * dPT1 + dPT2 * dPT2;
                scalePN2 += newPN * newPN;
                scalePT2 += newPT1 * newPT1 + newPT2 * newPT2;
            }
        }

        if constexpr (RollingEnabled) {
            converged = std::sqrt(delta2) <= atol + rtol * std::sqrt(scale2);
        }
        else {
            const Real delta = std::sqrt(deltaPN2) + std::sqrt(deltaPT2);
            const Real scale = std::sqrt(scalePN2) + std::sqrt(scalePT2);
            converged = delta <= atol + rtol * scale;
        }
        if (converged) break;
    }
    return { uE, converged };
}

/// Compatibility entry point for the three-component Coulomb specialization.
inline SorProjectResult sor_project(ContactSolverWorkspace& ws,
    std::size_t nN,
    const Mat6& mobility, const Vec6& freeVelocity, const Vec6& uA,
    Real mu, Real alpha, Real epsilonN, Real epsilonT,
    Real atol, Real rtol, std::uint32_t maxiter)
{
    return sor_project<false>(ws, nN, mobility, freeVelocity, uA,
        mu, Real(0), alpha, epsilonN, epsilonT, atol, rtol, maxiter);
}

/// Compatibility entry point for the five-component rolling specialization.
inline SorProjectResult sor_project_rolling(ContactSolverWorkspace& ws,
    std::size_t nN,
    const Mat6& mobility, const Vec6& freeVelocity, const Vec6& uA,
    Real mu, Real muR, Real alpha, Real epsilonN, Real epsilonT,
    Real atol, Real rtol, std::uint32_t maxiter)
{
    return sor_project<true>(ws, nN, mobility, freeVelocity, uA,
        mu, muR, alpha, epsilonN, epsilonT, atol, rtol, maxiter);
}
