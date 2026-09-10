#pragma once

#include "math_alias.hpp"

#include <cmath>

// Euclidean projection of p = [pN, pT1, pT2, pR1, pR2] onto the rolling
// friction cone
//   ||pT|| <= mu*pN,  ||pR|| <= muR*pN.
// This is the closed-form projection derived by Acary and Bourrier (2021).
inline void project_rolling_friction_cone(
    Real& pN, Real& pT1, Real& pT2, Real& pR1, Real& pR2,
    Real mu, Real muR) noexcept
{
    const Real normT = std::hypot(pT1, pT2);
    const Real normR = std::hypot(pR1, pR2);

    // The trial point lies in the polar cone.
    if (mu * normT + muR * normR <= -pN) {
        pN = pT1 = pT2 = pR1 = pR2 = Real(0);
        return;
    }

    // The trial point is already admissible.
    if (normT <= mu * pN && normR <= muR * pN) {
        return;
    }

    const Real mu2 = mu * mu;
    const Real muR2 = muR * muR;

    // Both sliding and rolling constraints are active.
    Real projectedN = (pN + mu * normT + muR * normR)
        / (Real(1) + mu2 + muR2);
    if (normT > mu * projectedN && normR > muR * projectedN) {
        pN = projectedN;
        const Real scaleT = normT > Real(0) ? mu * pN / normT : Real(0);
        const Real scaleR = normR > Real(0) ? muR * pN / normR : Real(0);
        pT1 *= scaleT;
        pT2 *= scaleT;
        pR1 *= scaleR;
        pR2 *= scaleR;
        return;
    }

    // Sliding constraint only.
    projectedN = (pN + mu * normT) / (Real(1) + mu2);
    if (normT > mu * projectedN && normR <= muR * projectedN) {
        pN = projectedN;
        const Real scaleT = normT > Real(0) ? mu * pN / normT : Real(0);
        pT1 *= scaleT;
        pT2 *= scaleT;
        return;
    }

    // Rolling constraint only. The three cases are exhaustive apart from
    // roundoff exactly where two projection regions meet.
    projectedN = (pN + muR * normR) / (Real(1) + muR2);
    pN = projectedN;
    const Real scaleR = normR > Real(0) ? muR * pN / normR : Real(0);
    pR1 *= scaleR;
    pR2 *= scaleR;

    // Protect the shared boundary against a one-ulp tangential overshoot.
    const Real projectedNormT = std::hypot(pT1, pT2);
    if (projectedNormT > mu * pN && projectedNormT > Real(0)) {
        const Real scaleT = mu * pN / projectedNormT;
        pT1 *= scaleT;
        pT2 *= scaleT;
    }
}
