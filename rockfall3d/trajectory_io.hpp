#pragma once
#include "time_stepper.hpp"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

// only for legacy text trajectory writing
inline TimeStepperResult load_trajectory_txt(const std::string& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open trajectory file: " + path);
    }

    TimeStepperResult result;
    Real t{};
    Eigen::Matrix<Real, 7, 1> q;
    Vec6 u;
    while (in >> t
        >> q(0) >> q(1) >> q(2) >> q(3) >> q(4) >> q(5) >> q(6)
        >> u(0) >> u(1) >> u(2) >> u(3) >> u(4) >> u(5)) {
        result.t.push_back(t);
        result.q.push_back(q);
        result.u.push_back(u);
    }

    if (result.t.empty()) {
        throw std::runtime_error("trajectory file is empty or malformed: " + path);
    }
    return result;
}

inline std::size_t write_trajectory_txt(const std::string& path,
    const TimeStepperResult& result,
    int precision = 6,
    Real dumpStep = Real(0))
{
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot open trajectory output file: " + path);
    }

    out << std::fixed << std::setprecision(precision);
    const bool useDecimation = dumpStep > Real(0);
    Real nextDumpTime = result.t.empty() ? Real(0) : result.t.front();
    const Real eps = Real(1e-12);
    std::size_t written = 0;

    for (std::size_t i = 0; i < result.t.size(); ++i) {
        const bool isLast = (i + 1 == result.t.size());
        if (useDecimation && !isLast && result.t[i] + eps < nextDumpTime) {
            continue;
        }

        out << result.t[i] << ' '
            << result.q[i].transpose() << ' '
            << result.u[i].transpose() << '\n';
        ++written;

        if (useDecimation) {
            nextDumpTime = result.t[i] + dumpStep;
        }
    }
    return written;
}
