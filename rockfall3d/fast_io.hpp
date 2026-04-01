#pragma once
#include "math_alias.hpp"
#include <fstream>

// Loads point-cloud from *.pts filepath
inline MatX load_xyz(const std::string& path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open rock file: " + path);
    }

    std::vector<Real> buf;
    Real x, y, z;
    while (in >> x >> y >> z) { buf.insert(buf.end(), { x,y,z }); }
    if (buf.empty()) {
        throw std::runtime_error("rock file is empty or malformed: " + path);
    }
    std::size_t n = buf.size() / 3;
    MatX P(3, n);
    for (std::size_t i = 0; i < n; ++i)
        P.col(i) << buf[3 * i], buf[3 * i + 1], buf[3 * i + 2];
    return P;
}
