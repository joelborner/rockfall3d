#pragma once

#include "math_alias.hpp"

#include <cstdint>
#include <vector>

// Convex polyhedron representation and mass properties computation
struct ConvexPolyhedron
{
    MatX vertices; // 3 x n body-frame vertices, COM at origin
    std::vector<std::vector<std::uint32_t>> faces; // list of faces, each face is a list of vertex indices
};

// Mass properties of a convex polyhedron
struct PolyhedronMassProperties
{
	Real volume{ Real(0) }; // volume of the polyhedron
	Vec3 center = Vec3::Zero(); // center of mass of the polyhedron
    Mat3 inertiaAtCenter = Mat3::Zero(); // unit-density inertia
};

ConvexPolyhedron make_convex_polyhedron(MatX points); // Build a convex polyhedron from a point cloud, removing duplicate and collinear vertices.

// Clip a convex polyhedron with a plane defined by an outward normal and an offset from the origin. The resulting polyhedron is the portion of the original polyhedron that lies on the side of the plane pointed to by the outward normal.
ConvexPolyhedron clip_convex_polyhedron(
	const ConvexPolyhedron& polyhedron,
    const Vec3& outwardNormal, 
    Real offset);

// Compute the volume of a convex polyhedron.
Real convex_polyhedron_volume(const ConvexPolyhedron& polyhedron);

// Compute the mass properties of a convex polyhedron, including volume, center of mass, and unit-density inertia tensor.
PolyhedronMassProperties convex_polyhedron_mass_properties(
    const ConvexPolyhedron& polyhedron);

