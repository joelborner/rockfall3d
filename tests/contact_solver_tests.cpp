#include "sor_solver.hpp"

#include <iostream>

namespace {

bool close(const Vec6& a, const Vec6& b, Real tol = Real(1e-8))
{
    return (a - b).norm() <= tol * std::max(Real(1), b.norm());
}

template <bool Rolling>
void initialize(ContactSolverWorkspace& ws)
{
    // Leave extra capacity to catch accidental use of inactive columns.
    ws.ensureContactCapacity(3);
    ws.WN.setZero(); ws.WT1.setZero(); ws.WT2.setZero();
    ws.WN(2, 0) = 1; ws.WT1(0, 0) = 1; ws.WT2(1, 0) = 1;
    ws.PNsub.setZero(); ws.PTsub.setZero();
    if constexpr (Rolling) {
        ws.ensureRollingCapacity(3);
        ws.WR1.setZero(); ws.WR2.setZero(); ws.PRsub.setZero();
        ws.WR1(3, 0) = 1; ws.WR2(4, 0) = 1;
    }
}

template <bool Sor, bool Rolling>
auto solve(ContactSolverWorkspace& ws, std::size_t n, const Vec6& free,
    const Vec6& before, Real muR, Real en, Real et,
    Real atol, Real rtol, std::uint32_t maxiter, Real alpha = Real(0.25))
{
    if constexpr (Sor) {
        return sor_project<Rolling>(ws, n, Mat6::Identity(), free, before,
            Real(0.5), muR, alpha, en, et, atol, rtol, maxiter);
    }
    else {
        return jor_project<Rolling>(ws, n, Mat6::Identity(), free, before,
            Real(0.5), muR, alpha, en, et, atol, rtol, maxiter);
    }
}

template <bool Rolling>
Vec6 reconstructed(const ContactSolverWorkspace& ws, const Vec6& free)
{
    Vec6 u = free + ws.MWN.col(0) * ws.PNsub(0)
        + ws.MWT1.col(0) * ws.PTsub(0, 0)
        + ws.MWT2.col(0) * ws.PTsub(0, 1);
    if constexpr (Rolling) {
        u += ws.MWR1.col(0) * ws.PRsub(0, 0)
            + ws.MWR2.col(0) * ws.PRsub(0, 1);
    }
    return u;
}

template <bool Sor, bool Rolling>
bool test_empty()
{
    ContactSolverWorkspace ws;
    Vec6 free; free << 1, 2, 3, 4, 5, 6;
    const auto r = solve<Sor, Rolling>(ws, 0, free, free, 0, 0, 0, 0, 0, 0);
    return r.converged && close(r.velocity, free, 0) && ws.WN.size() == 0;
}

template <bool Sor, bool Rolling>
bool test_return_consistency()
{
    Vec6 free = Vec6::Zero(); free.z() = -1;
    // Loose tolerance forces acceptance after a nonzero first update. This
    // exposes JOR's previous one-iterate lag even on a 'converged' return.
    for (int exitCase = 0; exitCase < 3; ++exitCase) {
        ContactSolverWorkspace ws;
        initialize<Rolling>(ws);
        ws.PNsub(0) = Real(0.1);
        const bool accepted = exitCase == 0;
        const auto r = solve<Sor, Rolling>(ws, 1, free, free, Real(0.3),
            Real(0.2), Real(0.1), accepted ? Real(10) : Real(0), 0,
            exitCase == 2 ? 0u : 1u);
        if (r.converged != accepted
            || !close(r.velocity, reconstructed<Rolling>(ws, free), Real(1e-13))) return false;
        if constexpr (!Rolling) {
            if (ws.WR1.size() || ws.WR2.size() || ws.PRsub.size()
                || ws.MWR1.size() || ws.MWR2.size() || ws.rollingStep.size()) return false;
        }
    }
    return true;
}

template <bool Sor, bool Rolling>
bool test_contact_law()
{
    for (Real en : { Real(0), Real(0.2) }) {
        for (Real muR : { Real(0), Real(0.3) }) {
            ContactSolverWorkspace ws;
            initialize<Rolling>(ws);
            Vec6 free = Vec6::Zero();
            free.x() = 2; free.z() = -1; free(3) = 2; free(5) = 3;
            const auto r = solve<Sor, Rolling>(ws, 1, free, free, muR,
                en, Real(0.1), Real(1e-12), Real(1e-10), 4000);
            Vec6 expected = free;
            expected.z() = en;
            expected.x() -= Real(0.5) * (Real(1) + en);
            if constexpr (Rolling) expected(3) -= muR * (Real(1) + en);
            if (!r.converged || !close(r.velocity, expected)
                || !close(r.velocity, reconstructed<Rolling>(ws, free), Real(1e-12))) return false;
        }
    }
    return true;
}

template <bool Sor, bool Rolling>
bool test_coupled_contacts()
{
    ContactSolverWorkspace ws;
    initialize<Rolling>(ws);
    ws.WN(4, 0) = 1;
    ws.WN(2, 1) = 1; ws.WN(4, 1) = -1;
    ws.WT1(0, 1) = 1; ws.WT2(1, 1) = 1;
    if constexpr (Rolling) { ws.WR1(3, 1) = 1; ws.WR2(4, 1) = 1; }
    Vec6 free = Vec6::Zero(); free.z() = -1;
    const auto r = solve<Sor, Rolling>(ws, 2, free, free, Real(0.3),
        0, 0, Real(1e-12), Real(1e-10), 4000);
    return r.converged && close(r.velocity, Vec6::Zero());
}

template <bool Sor, bool Rolling>
bool run_tests()
{
    return test_empty<Sor, Rolling>() && test_return_consistency<Sor, Rolling>()
        && test_contact_law<Sor, Rolling>() && test_coupled_contacts<Sor, Rolling>();
}

} // namespace

int main()
{
    if (!run_tests<false, false>()) { std::cerr << "JOR Coulomb regression failed\n"; return 1; }
    if (!run_tests<false, true>()) { std::cerr << "JOR rolling regression failed\n"; return 1; }
    if (!run_tests<true, false>()) { std::cerr << "SOR Coulomb regression failed\n"; return 1; }
    if (!run_tests<true, true>()) { std::cerr << "SOR rolling regression failed\n"; return 1; }
    return 0;
}
