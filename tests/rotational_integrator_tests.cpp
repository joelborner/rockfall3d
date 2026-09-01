#include "rotational_integrator.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close_relative(Real actual, Real expected, Real tolerance)
{
    return std::abs(actual - expected)
        <= tolerance * std::max({ Real(1), std::abs(actual), std::abs(expected) });
}

Real rotational_energy(const Mat3& inertia, const Vec3& omega)
{
    return Real(0.5) * omega.dot(inertia * omega);
}

bool test_schweizer_energy()
{
    Mat3 inertia;
    inertia << Real(2.0), Real(0.2), Real(-0.1),
        Real(0.2), Real(3.0), Real(0.15),
        Real(-0.1), Real(0.15), Real(5.0);
    Vec6 uA = Vec6::Zero();
    uA.tail<3>() << Real(1.3), Real(-0.7), Real(2.1);

    const IntegratorFreeStep step = SchweizerIntegrator{}.prepareFreeStep(
        Real(1), inertia, inertia.inverse(), uA, Real(0.01));
    return step.converged
        && close_relative(
            rotational_energy(inertia, step.velocity.tail<3>()),
            rotational_energy(inertia, uA.tail<3>()),
            Real(5e-13));
}

bool test_capobianco_invariants()
{
    Mat3 inertia;
    inertia << Real(2.0), Real(0.2), Real(-0.1),
        Real(0.2), Real(3.0), Real(0.15),
        Real(-0.1), Real(0.15), Real(5.0);
    const Mat3 inertiaInv = inertia.inverse();
    const Real dt = Real(0.002);
    const CapobiancoIntegrator integrator{ Real(1e-12), Real(1e-12), 50 };

    Vec6 u = Vec6::Zero();
    u.tail<3>() << Real(1.3), Real(-0.7), Real(2.1);
    Eigen::Matrix<Real, 7, 1> q;
    q << Real(0), Real(0), Real(0), Real(1), Real(0), Real(0), Real(0);

    const Real energy0 = rotational_energy(inertia, u.tail<3>());
    const Vec3 spatialSpin0 = quaternion_from_state(q) * (inertia * u.tail<3>());

    for (int stepIndex = 0; stepIndex < 5000; ++stepIndex) {
        const auto qM = make_midpoint_configuration(q, u, dt);
        const IntegratorFreeStep freeStep = integrator.prepareFreeStep(
            Real(1), inertia, inertiaInv, u, dt);
        if (!freeStep.converged) return false;
        q = integrator.finishConfiguration(q, qM, u, freeStep.velocity,
            inertia, freeStep, dt);
        u = freeStep.velocity;
        u.head<3>().setZero();
    }

    const Real energy = rotational_energy(inertia, u.tail<3>());
    const Vec3 spatialSpin = quaternion_from_state(q) * (inertia * u.tail<3>());
    return close_relative(energy, energy0, Real(5e-10))
        && (spatialSpin - spatialSpin0).norm()
            <= Real(5e-10) * std::max(Real(1), spatialSpin0.norm());
}

} // namespace

int main()
{
    if (!test_schweizer_energy()) {
        std::cerr << "Schweizer energy-conservation test failed\n";
        return 1;
    }
    if (!test_capobianco_invariants()) {
        std::cerr << "Capobianco energy/spin-conservation test failed\n";
        return 1;
    }
    return 0;
}

