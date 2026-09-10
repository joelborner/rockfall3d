#include "fragmentation.hpp"
#include "rotational_integrator.hpp"

#include <cmath>
#include <iostream>

namespace {

bool check(bool condition, const char* message)
{
    if (!condition) std::cerr << "FAILED: " << message << '\n';
    return condition;
}

} // namespace

int main()
{
    MatX vertices(3, 8);
    Eigen::Index column = 0;
    for (Real x : { Real(-1), Real(1) }) {
        for (Real y : { Real(-1), Real(1) }) {
            for (Real z : { Real(-1), Real(1) }) {
                vertices.col(column++) = Vec3(x, y, z);
            }
        }
    }

    Eigen::Matrix<Real, 7, 1> pose;
    pose << 0, 0, 0, 2, 0, 0, 0;
    const Rock original(std::move(vertices), Real(2500), Vec6::Zero(), pose);

    bool ok = true;
    ok &= check(std::abs(original.q0()(3) - Real(1)) < Real(1e-14),
        "rock construction normalizes its quaternion");

    FragmentationConfig config;
    config.thresholdSpecificEnergy = Real(0);
    config.halfSpecificEnergy = Real(100);
    config.exponent = Real(1);
    config.maxRemovedFraction = Real(0.1);

    const std::vector<FragmentationContact> contacts{
        FragmentationContact{ Vec3(1, 1, 1), original.mass() * Real(50) }
    };
    auto result = fragment_rock(original, contacts, config);
    ok &= check(result.has_value(), "impact fragments the test cube");
    if (!result) return 1;

    const Rock& fragment = result->rock;
    ok &= check(fragment.mass() < original.mass(), "fragment loses mass");
    ok &= check(fragment.equivalentRadius() < original.equivalentRadius(),
        "fragment updates its equivalent radius");
    ok &= check(fragment.inertiaK().allFinite(), "fragment inertia is finite");
    const auto fragmentProperties =
        convex_polyhedron_mass_properties(fragment.polyhedronK());
    ok &= check(fragmentProperties.center.norm() < Real(1e-12),
        "fragment geometry remains centered at its center of mass");

    const Mat3 inertiaInv = fragment.inertiaK().inverse();
    const Vec6 velocity = (Vec6() << 1, 2, 3, 4, 5, 6).finished();
    const auto schweizer = SchweizerIntegrator{}.prepareFreeStep(
        fragment.mass(), fragment.inertiaK(), inertiaInv, velocity, Real(0.002));
    const auto capobianco = CapobiancoIntegrator{}.prepareFreeStep(
        fragment.mass(), fragment.inertiaK(), inertiaInv, velocity, Real(0.002));
    ok &= check(schweizer.converged, "Schweizer accepts fragment properties");
    ok &= check(capobianco.converged, "Capobianco accepts fragment properties");

    return ok ? 0 : 1;
}
