#pragma once

#include "math_alias.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>

/// Available discretizations of the torque-free rigid-body rotation.
enum class RotationalIntegratorType : std::uint8_t
{
    /// Energy-consistent pseudo-implicit scheme of Schweizer (2015).
    Schweizer,
    /// Stability-preserving scheme of Leine et al. (2021), called CAPOBIANCO2
    /// in the predecessor code.
    Capobianco
};

/// User-configurable rotational integration parameters.
struct RotationalIntegratorSettings
{
    /// Selected rotational integration scheme.
    RotationalIntegratorType type{ RotationalIntegratorType::Schweizer };
    /// Absolute tolerance of the Capobianco free-rotation Newton solve.
    Real atol{ Real(1e-10) };
    /// Relative tolerance of the Capobianco free-rotation Newton solve.
    Real rtol{ Real(1e-10) };
    /// Maximum Newton iterations per Capobianco free-rotation step.
    std::uint32_t maxiter{ 50 };
};

/// Integrator-specific free motion supplied to the proximal contact solver.
///
/// Both contact solvers use the affine velocity relation
/// `uE = velocity + mobility * W * P`, where `P` contains contact impulses.
/// Generalized velocities are ordered as `[v_I; omega_K]`: translational
/// velocity in the inertial frame followed by angular velocity in the body
/// frame.
struct IntegratorFreeStep
{
    /// End-of-step generalized velocity in the absence of contact impulses.
    Vec6 velocity = Vec6::Zero();
    /// Algorithmic operator mapping a generalized impulse to a velocity change.
    Mat6 mobility = Mat6::Zero();
    /// Torque-free end angular velocity used by the Capobianco orientation map.
    Vec3 omegaFree = Vec3::Zero();
    /// Whether all integrator-specific nonlinear and linear operations succeeded.
    bool converged{ true };
};

/// Extract `(e0,e1,e2,e3)` from a configuration whose quaternion is normalized.
inline Quat quaternion_from_state(const Eigen::Matrix<Real, 7, 1>& q)
{
    Quat result(q(3), q(4), q(5), q(6));
    assert(result.coeffs().allFinite());
    assert(std::abs(result.squaredNorm() - Real(1))
        <= std::sqrt(std::numeric_limits<Real>::epsilon()));
    return result;
}

/// Store a normalized quaternion in the orientation part of a state vector.
inline void write_quaternion(Eigen::Matrix<Real, 7, 1>& q, const Quat& value)
{
    Quat normalized = value;
    normalized.normalize();
    q.segment<4>(3) << normalized.w(), normalized.x(), normalized.y(), normalized.z();
}

/// Advance a quaternion by an explicit Euler step.
///
/// Uses `qdot = 1/2 * q * (0,omegaK)` with body-frame angular velocity
/// `omegaK` and step length `step` in seconds. The caller normalizes the
/// quaternion when storing the resulting configuration.
inline Quat quaternion_euler_step(const Quat& q, const Vec3& omegaK, Real step)
{
    const Quat omega(Real(0), omegaK.x(), omegaK.y(), omegaK.z());
    Quat result;
    result.coeffs() = q.coeffs() + (q * omega).coeffs() * (step * Real(0.5));
    return result;
}

/// Construct the Moreau midpoint configuration from the beginning state.
///
/// Translation and orientation are advanced by `dt/2` using the beginning
/// generalized velocity `uA`.
inline Eigen::Matrix<Real, 7, 1> make_midpoint_configuration(
    const Eigen::Matrix<Real, 7, 1>& qA,
    const Vec6& uA,
    Real dt)
{
    Eigen::Matrix<Real, 7, 1> qM = qA;
    qM.head<3>().noalias() += uA.head<3>() * (dt * Real(0.5));
    write_quaternion(qM,
        quaternion_euler_step(quaternion_from_state(qA), uA.tail<3>(), dt * Real(0.5)));
    return qM;
}

/// Assemble `diag(massInv*I, rotationalMobility)` in generalized coordinates.
inline Mat6 make_generalized_mobility(Real massInv, const Mat3& rotationalMobility)
{
    Mat6 mobility = Mat6::Zero();
    mobility.block<3, 3>(0, 0) = Mat3::Identity() * massInv;
    mobility.block<3, 3>(3, 3) = rotationalMobility;
    return mobility;
}

/// Schweizer's energy-consistent pseudo-implicit rotational integrator.
///
/// The gyroscopic term is incorporated through the effective inertia
/// `inertiaHat = inertia + dt/2 * (inertia*omegaA~ + omegaA~*inertia)`.
/// Consequently, contact impulses change angular velocity through
/// `inertiaHat^-1`, rather than through the physical inverse inertia alone.
struct SchweizerIntegrator
{
    /// Compute torque-free motion and the impulse-to-velocity operator.
    ///
    /// @param mass Rock mass in kg.
    /// @param inertia Constant body-frame inertia tensor in kg m^2.
    /// @param inertiaInv Unused; retained for the common integrator contract.
    /// @param uA Beginning generalized velocity `[v_I; omega_K]`.
    /// @param dt Timestep in seconds.
    IntegratorFreeStep prepareFreeStep(Real mass,
        const Mat3& inertia,
        [[maybe_unused]] const Mat3& inertiaInv,
        const Vec6& uA,
        Real dt) const
    {
        const Real massInv = Real(1) / mass;
        const Vec3 omegaA = uA.tail<3>();
        const Mat3 gyro = inertia * tilde(omegaA) + tilde(omegaA) * inertia;
        const Mat3 inertiaHat = inertia + (dt * Real(0.5)) * gyro;
        const Mat3 rotationalMobility = inertiaHat.inverse();

        IntegratorFreeStep result;
        result.velocity = uA;
        result.velocity.z() -= Real(9.81) * dt;
        result.velocity.tail<3>().noalias() = rotationalMobility
            * (inertia - (dt * Real(0.5)) * gyro) * omegaA;
        result.mobility = make_generalized_mobility(massInv, rotationalMobility);
        result.omegaFree = result.velocity.tail<3>();
        result.converged = result.velocity.allFinite() && result.mobility.allFinite();
        return result;
    }

    /// Reconstruct the end configuration using the accepted contact velocity.
    ///
    /// @param qA Unused; retained for the common integrator contract.
    /// @param qM Moreau midpoint configuration.
    /// @param uA Unused; retained for the common integrator contract.
    /// @param uE End generalized velocity after contact impulses.
    /// @param inertia Unused; retained for the common integrator contract.
    /// @param freeStep Unused; retained for the common integrator contract.
    /// @param dt Timestep in seconds.
    Eigen::Matrix<Real, 7, 1> finishConfiguration(
        [[maybe_unused]] const Eigen::Matrix<Real, 7, 1>& qA,
        const Eigen::Matrix<Real, 7, 1>& qM,
        [[maybe_unused]] const Vec6& uA,
        const Vec6& uE,
        [[maybe_unused]] const Mat3& inertia,
        [[maybe_unused]] const IntegratorFreeStep& freeStep,
        Real dt) const
    {
        Eigen::Matrix<Real, 7, 1> qE = qM;
        qE.head<3>().noalias() += uE.head<3>() * (dt * Real(0.5));
        write_quaternion(qE,
            quaternion_euler_step(quaternion_from_state(qM), uE.tail<3>(), dt * Real(0.5)));
        return qE;
    }
};

/// Published stability-preserving Capobianco rotational integrator.
///
/// The torque-free angular velocity satisfies the implicit midpoint Euler
/// equation `inertia*(omegaE-omegaA) + dt*omegaM x inertia*omegaM = 0`.
/// A subsequent orientation correction preserves inertial angular momentum.
struct CapobiancoIntegrator
{
    /// Absolute Newton convergence tolerance.
    Real atol{ Real(1e-10) };
    /// Relative Newton convergence tolerance.
    Real rtol{ Real(1e-10) };
    /// Maximum Newton iterations per timestep.
    std::uint32_t maxiter{ 50 };

    /// Test the implicit midpoint Euler-equation residual against the configured tolerances.
    bool freeRotationConverged(const Mat3& inertia,
        const Vec3& omegaA,
        const Vec3& omegaE,
        Real dt) const
    {
        const Vec3 omegaM = (omegaA + omegaE) * Real(0.5);
        const Vec3 residual = inertia * (omegaE - omegaA)
            + dt * omegaM.cross(inertia * omegaM);
        const Real scale = std::max({
            (inertia * omegaA).norm(),
            (inertia * omegaE).norm(),
            dt * omegaM.cross(inertia * omegaM).norm(),
            Real(1)
            });
        return residual.norm() <= atol + rtol * scale;
    }

    /// Solve torque-free rotation and provide the physical inverse-mass mobility.
    ///
    /// @param mass Rock mass in kg.
    /// @param inertia Constant body-frame inertia tensor in kg m^2.
    /// @param inertiaInv Inverse body-frame inertia tensor.
    /// @param uA Beginning generalized velocity `[v_I; omega_K]`.
    /// @param dt Timestep in seconds.
    IntegratorFreeStep prepareFreeStep(Real mass,
        const Mat3& inertia,
        const Mat3& inertiaInv,
        const Vec6& uA,
        Real dt) const
    {
        const Vec3 omegaA = uA.tail<3>();
        Vec3 omegaE = omegaA;
        bool converged = false;

        for (std::uint32_t iter = 0; iter < maxiter; ++iter) {
            const Vec3 omegaM = (omegaA + omegaE) * Real(0.5);
            const Vec3 inertiaOmegaM = inertia * omegaM;
            const Vec3 residual = inertia * (omegaE - omegaA)
                + dt * omegaM.cross(inertiaOmegaM);

            if (freeRotationConverged(inertia, omegaA, omegaE, dt)) {
                converged = true;
                break;
            }

            const Mat3 gyroJacobian = -tilde(inertiaOmegaM) + tilde(omegaM) * inertia;
            const Mat3 jacobian = inertia + (dt * Real(0.5)) * gyroJacobian;
            const Vec3 correction = jacobian.partialPivLu().solve(residual);
            if (!correction.allFinite()) {
                break;
            }
            omegaE -= correction;
        }

        if (!converged) {
            converged = freeRotationConverged(inertia, omegaA, omegaE, dt);
        }

        IntegratorFreeStep result;
        result.velocity = uA;
        result.velocity.z() -= Real(9.81) * dt;
        result.velocity.tail<3>() = omegaE;
        result.mobility = make_generalized_mobility(Real(1) / mass, inertiaInv);
        result.omegaFree = omegaE;
        result.converged = converged && result.velocity.allFinite();
        return result;
    }

    /// Reconstruct the end configuration with the spin-preserving correction.
    ///
    /// Contact impulses affect `uE`, while the orientation map uses the
    /// torque-free angular velocity stored in `freeStep`, as prescribed by the
    /// published splitting.
    Eigen::Matrix<Real, 7, 1> finishConfiguration(
        const Eigen::Matrix<Real, 7, 1>& qA,
        const Eigen::Matrix<Real, 7, 1>& qM,
        const Vec6& uA,
        const Vec6& uE,
        const Mat3& inertia,
        const IntegratorFreeStep& freeStep,
        Real dt) const
    {
        Eigen::Matrix<Real, 7, 1> qE = qA;
        qE.head<3>() = qM.head<3>() + uE.head<3>() * (dt * Real(0.5));

        const Vec3 omegaM = (uA.tail<3>() + freeStep.omegaFree) * Real(0.5);
        const Real angle = omegaM.norm() * dt;
        Quat midpointRotation = Quat::Identity();
        if (angle > Real(0)) {
            midpointRotation = Quat(Eigen::AngleAxis<Real>(angle, omegaM.normalized()));
        }

        // Express both angular momenta in the body frame reached by the
        // midpoint rotation before constructing the correction.
        const Vec3 spinA = inertia * uA.tail<3>();
        const Vec3 spinFree = inertia * freeStep.omegaFree;
        const Vec3 transportedSpinA = midpointRotation.conjugate() * spinA;

        // The shortest rotation aligning these vectors supplies the published
        // B correction without evaluating a general matrix exponential.
        Quat correction = Quat::Identity();
        const Real spinScale = std::max(spinA.norm(), spinFree.norm());
        if (spinScale > std::numeric_limits<Real>::epsilon()) {
            correction = Quat::FromTwoVectors(spinFree, transportedSpinA);
            correction.normalize();
        }

        const Quat qEnd = quaternion_from_state(qA) * midpointRotation * correction;
        write_quaternion(qE, qEnd);
        return qE;
    }
};
