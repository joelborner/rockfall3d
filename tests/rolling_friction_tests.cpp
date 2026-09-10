#include "jor_solver.hpp"
#include "rolling_friction.hpp"
#include "sor_solver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace {

bool close(Real actual, Real expected, Real tolerance = Real(1e-10))
{
    return std::abs(actual - expected)
        <= tolerance * std::max({ Real(1), std::abs(actual), std::abs(expected) });
}

bool cone_projection_cases()
{
    {
        Real n = 2, t1 = 0.5, t2 = 0, r1 = 0.25, r2 = 0;
        project_rolling_friction_cone(n, t1, t2, r1, r2, Real(0.5), Real(0.25));
        if (!close(n, 2) || !close(t1, 0.5) || !close(r1, 0.25)) return false;
    }
    {
        Real n = -2, t1 = 1, t2 = 0, r1 = 1, r2 = 0;
        project_rolling_friction_cone(n, t1, t2, r1, r2, Real(0.5), Real(0.25));
        if (n != Real(0) || t1 != Real(0) || r1 != Real(0)) return false;
    }
    {
        Real n = 0, t1 = 2, t2 = 0, r1 = 1, r2 = 0;
        const Real mu = Real(0.5);
        const Real muR = Real(0.25);
        const Real expectedN = (mu * Real(2) + muR) / (Real(1) + mu * mu + muR * muR);
        project_rolling_friction_cone(n, t1, t2, r1, r2, mu, muR);
        if (!close(n, expectedN) || !close(t1, mu * expectedN)
            || !close(r1, muR * expectedN)) return false;
    }
    return true;
}

void setup_single_contact(ContactSolverWorkspace& ws)
{
    ws.ensureContactCapacity(1);
    ws.ensureRollingCapacity(1);
    ws.WN.leftCols(1).setZero();
    ws.WT1.leftCols(1).setZero();
    ws.WT2.leftCols(1).setZero();
    ws.WR1.leftCols(1).setZero();
    ws.WR2.leftCols(1).setZero();
    ws.WN(2, 0) = Real(1);
    ws.WT1(0, 0) = Real(1);
    ws.WT2(1, 0) = Real(1);
    ws.WR1(3, 0) = Real(1);
    ws.WR2(4, 0) = Real(1);
    ws.PNsub(0) = Real(0);
    ws.PTsub.topRows(1).setZero();
    ws.PRsub.topRows(1).setZero();
}

template <typename Solve>
bool rolling_solver_sticks(Solve&& solve)
{
    ContactSolverWorkspace ws;
    setup_single_contact(ws);

    Vec6 freeVelocity = Vec6::Zero();
    freeVelocity.z() = Real(-1);
    freeVelocity(3) = Real(0.2);
    const auto result = solve(ws, freeVelocity);
    return result.converged
        && close(result.velocity.z(), Real(0), Real(1e-8))
        && close(result.velocity(3), Real(0), Real(1e-8))
        && ws.PNsub(0) >= Real(0)
        && std::hypot(ws.PRsub(0, 0), ws.PRsub(0, 1))
            <= Real(0.5) * ws.PNsub(0) + Real(1e-10);
}

bool rolling_solvers()
{
    const Mat6 mobility = Mat6::Identity();
    const Vec6 uA = Vec6::Zero();
    const auto jor = [&](ContactSolverWorkspace& ws, const Vec6& freeVelocity) {
        return jor_project_rolling(ws, 1, mobility, freeVelocity, uA,
            Real(0.5), Real(0.5), Real(0.5), Real(0), Real(0),
            Real(1e-12), Real(1e-10), 10000);
    };
    const auto sor = [&](ContactSolverWorkspace& ws, const Vec6& freeVelocity) {
        return sor_project_rolling(ws, 1, mobility, freeVelocity, uA,
            Real(0.5), Real(0.5), Real(0.5), Real(0), Real(0),
            Real(1e-12), Real(1e-10), 10000);
    };
    return rolling_solver_sticks(jor) && rolling_solver_sticks(sor);
}

bool drilling_spin_is_unresisted()
{
    ContactSolverWorkspace ws;
    setup_single_contact(ws);
    Vec6 freeVelocity = Vec6::Zero();
    freeVelocity.z() = Real(-1);
    freeVelocity(5) = Real(2);

    const auto result = sor_project_rolling(ws, 1, Mat6::Identity(),
        freeVelocity, Vec6::Zero(), Real(0.5), Real(0.5), Real(1),
        Real(0), Real(0), Real(1e-12), Real(1e-10), 10000);
    return result.converged
        && close(result.velocity.z(), Real(0), Real(1e-8))
        && close(result.velocity(5), Real(2), Real(1e-12));
}

} // namespace

int main()
{
    if (!cone_projection_cases()) {
        std::cerr << "Rolling-friction cone projection test failed\n";
        return 1;
    }
    if (!rolling_solvers()) {
        std::cerr << "Rolling-friction solver test failed\n";
        return 1;
    }
    if (!drilling_spin_is_unresisted()) {
        std::cerr << "Rolling-friction drilling-spin test failed\n";
        return 1;
    }
    return 0;
}
