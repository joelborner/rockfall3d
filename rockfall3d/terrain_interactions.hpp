#pragma once

#include "math_alias.hpp"

#include <cmath>
#include <cstdint>

enum class CoulombFrictionModel : std::uint8_t
{
    Constant,
    SlippageDependent
};

struct CoulombFrictionSettings
{
    CoulombFrictionModel model{ CoulombFrictionModel::Constant };
    Real mu{ Real(0.5) };
    Real muMin{ Real(0) };
    Real muMax{ Real(0) };
    Real kappa{ Real(0) };
    Real beta{ Real(0) };

    bool slippageDependent() const noexcept
    {
        return model == CoulombFrictionModel::SlippageDependent;
    }

    Real coefficient(Real slippage) const noexcept
    {
        constexpr Real twoOverPi =
            Real(0.636619772367581343075535053490057448L);
        return muMin + (muMax - muMin) * twoOverPi
            * std::atan(kappa * slippage);
    }
};

// Dynamic terrain-interaction state belonging to one simulated rock trajectory.
// It must not be shared between parallel runs of the same Rock geometry.
struct TerrainInteractionState
{
    Real slippage{ Real(0) };
};

// Optional terrain-contact mechanisms are configured independently so they can
// be composed without changing the base contact formulation.
struct TerrainInteractionSettings
{
    CoulombFrictionSettings coulombFriction;

    // Dimensionless terrain coefficient. As in PlatRock, the solver multiplies
    // it by the equivalent rock radius to obtain the rolling-resistance length.
    Real rollingFrictionCoefficient{ Real(0) };

    bool rollingFrictionEnabled() const noexcept
    {
        return rollingFrictionCoefficient > Real(0);
    }
};

inline void update_slippage_state(
    TerrainInteractionState& state,
    const CoulombFrictionSettings& settings,
    bool hasActiveTerrainContact,
    Real midpointComSpeed,
    Real dt) noexcept
{
    if (hasActiveTerrainContact) {
        state.slippage += dt * midpointComSpeed;
    }
    else {
        state.slippage *= std::exp(-settings.beta * dt);
    }
}
