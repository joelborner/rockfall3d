#include "terrain_interactions.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close_relative(Real actual, Real expected, Real tolerance = Real(1e-13))
{
    return std::abs(actual - expected)
        <= tolerance * std::max({ Real(1), std::abs(actual), std::abs(expected) });
}

bool test_slippage_dependent_coefficient()
{
    CoulombFrictionSettings friction;
    friction.model = CoulombFrictionModel::SlippageDependent;
    friction.muMin = Real(0.2);
    friction.muMax = Real(1.3);
    friction.kappa = Real(0.7);

    return close_relative(friction.coefficient(Real(0)), friction.muMin)
        && close_relative(
            friction.coefficient(Real(1) / friction.kappa),
            Real(0.5) * (friction.muMin + friction.muMax))
        && friction.coefficient(Real(1e12)) < friction.muMax
        && close_relative(friction.coefficient(Real(1e12)), friction.muMax, Real(1e-11));
}

bool test_slippage_state_evolution()
{
    CoulombFrictionSettings friction;
    friction.model = CoulombFrictionModel::SlippageDependent;
    friction.beta = Real(20);

    TerrainInteractionState state;
    update_slippage_state(state, friction, true, Real(3), Real(0.1));
    if (!close_relative(state.slippage, Real(0.3))) return false;

    update_slippage_state(state, friction, false, Real(100), Real(0.05));
    return close_relative(state.slippage, Real(0.3) * std::exp(Real(-1)));
}

} // namespace

int main()
{
    if (!test_slippage_dependent_coefficient()) {
        std::cerr << "Slippage-dependent Coulomb coefficient test failed\n";
        return 1;
    }
    if (!test_slippage_state_evolution()) {
        std::cerr << "Slippage state evolution test failed\n";
        return 1;
    }
    return 0;
}
