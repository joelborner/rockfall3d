#include "convex_polyhedron.hpp"
#include "quickhull.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace {

using Face = std::vector<std::uint32_t>;
using Edge = std::pair<std::uint32_t, std::uint32_t>;

// Compute a scale factor for the geometry based on the maximum extent of the points. This is used to normalize tolerances for geometric operations.
Real geometry_scale(const MatX& points)
{
    if (points.cols() == 0) return Real(1);
    const Vec3 extent = points.rowwise().maxCoeff() - points.rowwise().minCoeff();
    return std::max<Real>(Real(1), extent.maxCoeff());
}

// Compute the normal vector of a face defined by its vertex indices. The normal is computed as the sum of the cross products of consecutive edges of the face.
Vec3 face_normal(const MatX& vertices, const Face& face)
{
    Vec3 normal = Vec3::Zero();
    for (std::size_t i = 0; i < face.size(); ++i) {
        const Vec3 a = vertices.col(static_cast<Eigen::Index>(face[i]));
        const Vec3 b = vertices.col(static_cast<Eigen::Index>(face[(i + 1) % face.size()]));
        normal += a.cross(b);
    }
    return normal;
}

// Remove consecutive duplicate vertex indices from a face. If the first and last indices are the same, remove the last index to avoid a degenerate face.
void remove_consecutive_duplicates(Face& face)
{
    Face cleaned;
    cleaned.reserve(face.size());
    for (const std::uint32_t index : face) {
        if (cleaned.empty() || cleaned.back() != index) {
            cleaned.push_back(index);
        }
    }
    if (cleaned.size() > 1 && cleaned.front() == cleaned.back()) {
        cleaned.pop_back();
    }
    face = std::move(cleaned);
}

// Remove collinear vertices from a face. A vertex is considered collinear if it lies on the line segment between its two neighboring vertices, within a specified tolerance. 
// The function iteratively removes collinear vertices until no more can be removed or the face has only three vertices left.
void remove_collinear_vertices(Face& face, const MatX& vertices, Real tolerance)
{
    bool changed = true;
    while (changed && face.size() > 3) {
        changed = false;
        for (std::size_t i = 0; i < face.size(); ++i) {
            const Vec3 previous = vertices.col(static_cast<Eigen::Index>(
                face[(i + face.size() - 1) % face.size()]));
            const Vec3 current = vertices.col(static_cast<Eigen::Index>(face[i]));
            const Vec3 next = vertices.col(static_cast<Eigen::Index>(
                face[(i + 1) % face.size()]));
            const Vec3 a = current - previous;
            const Vec3 b = next - current;
            const Real scale = std::max<Real>(Real(1), a.norm() * b.norm());
            if (a.cross(b).norm() <= tolerance * scale && a.dot(b) >= Real(0)) {
                face.erase(face.begin() + static_cast<std::ptrdiff_t>(i));
                changed = true;
                break;
            }
        }
    }
}

// Compact a convex polyhedron by removing unused vertices and remapping the vertex indices in the faces. 
// The function checks that the resulting polyhedron has at least four vertices and four faces, throwing an exception if not.
ConvexPolyhedron compact_polyhedron(MatX vertices, std::vector<Face> faces)
{
    std::set<std::uint32_t> used;
    for (const Face& face : faces) {
        used.insert(face.begin(), face.end());
    }
    if (used.size() < 4 || faces.size() < 4) {
        throw std::runtime_error("fragmentation produced a non-3-D polyhedron");
    }

    std::vector<std::uint32_t> remap(static_cast<std::size_t>(vertices.cols()));
    MatX compact(3, static_cast<Eigen::Index>(used.size()));
    Eigen::Index next = 0;
    for (const std::uint32_t oldIndex : used) {
        remap[oldIndex] = static_cast<std::uint32_t>(next);
        compact.col(next++) = vertices.col(static_cast<Eigen::Index>(oldIndex));
    }
    for (Face& face : faces) {
        for (std::uint32_t& index : face) {
            index = remap[index];
        }
    }
    return ConvexPolyhedron{ std::move(compact), std::move(faces) };
}

struct PlaneGroup
{
    Vec3 normal = Vec3::Zero();
    Real offset{ Real(0) };
    std::vector<std::array<std::uint32_t, 3>> triangles;
};

// Merge coplanar triangles into convex polygonal faces. 
// The function groups triangles by their plane equations (normal and offset) within a specified tolerance,
std::vector<Face> merge_coplanar_triangles(
    const MatX& points,
    const std::vector<std::array<std::uint32_t, 3>>& triangles,
    Real tolerance)
{
    std::vector<PlaneGroup> groups;
    for (const auto& triangle : triangles) {
        const Vec3 a = points.col(triangle[0]);
        const Vec3 b = points.col(triangle[1]);
        const Vec3 c = points.col(triangle[2]);
        Vec3 normal = (b - a).cross(c - a);
        const Real norm = normal.norm();
        if (!(norm > tolerance)) continue;
        normal /= norm;
        const Real offset = normal.dot(a);

        PlaneGroup* group = nullptr;
        for (PlaneGroup& candidate : groups) {
            if ((candidate.normal - normal).norm() <= tolerance
                && std::abs(candidate.offset - offset) <= tolerance) {
                group = &candidate;
                break;
            }
        }
        if (group == nullptr) {
            groups.push_back(PlaneGroup{ normal, offset, {} });
            group = &groups.back();
        }
        group->triangles.push_back(triangle);
    }

    std::vector<Face> faces;
    faces.reserve(groups.size());
    for (const PlaneGroup& group : groups) {
        struct DirectedEdge { std::uint32_t from; std::uint32_t to; std::uint32_t count; };
        std::map<Edge, DirectedEdge> edges;
        for (const auto& triangle : group.triangles) {
            for (int i = 0; i < 3; ++i) {
                const std::uint32_t from = triangle[static_cast<std::size_t>(i)];
                const std::uint32_t to = triangle[static_cast<std::size_t>((i + 1) % 3)];
                const Edge key = std::minmax(from, to);
                auto [it, inserted] = edges.emplace(key, DirectedEdge{ from, to, 1 });
                if (!inserted) ++it->second.count;
            }
        }

        std::map<std::uint32_t, std::uint32_t> next;
        for (const auto& [key, edge] : edges) {
            (void)key;
            if (edge.count == 1) {
                next.emplace(edge.from, edge.to);
            }
        }
        if (next.size() < 3) {
            throw std::runtime_error("failed to reconstruct convex hull face boundary");
        }

        Face face;
        face.reserve(next.size());
        const std::uint32_t first = next.begin()->first;
        std::uint32_t current = first;
        do {
            face.push_back(current);
            const auto it = next.find(current);
            if (it == next.end()) {
                throw std::runtime_error("convex hull face boundary is not closed");
            }
            current = it->second;
        } while (current != first && face.size() <= next.size());
        if (current != first || face.size() != next.size()) {
            throw std::runtime_error("convex hull face boundary contains multiple loops");
        }
        remove_collinear_vertices(face, points, tolerance);
        faces.push_back(std::move(face));
    }
    return faces;
}

// Append a vertex to the list of vertices if it is not already present within a specified tolerance.
std::uint32_t append_vertex(std::vector<Vec3>& vertices, const Vec3& point, Real tolerance)
{
    const Real tolerance2 = tolerance * tolerance;
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        if ((vertices[i] - point).squaredNorm() <= tolerance2) {
            return static_cast<std::uint32_t>(i);
        }
    }
    vertices.push_back(point);
    return static_cast<std::uint32_t>(vertices.size() - 1);
}

// Convert a list of vertices into a matrix where each column represents a vertex.
MatX matrix_from_vertices(const std::vector<Vec3>& vertices)
{
    MatX result(3, static_cast<Eigen::Index>(vertices.size()));
    for (std::size_t i = 0; i < vertices.size(); ++i) {
        result.col(static_cast<Eigen::Index>(i)) = vertices[i];
    }
    return result;
}

} // namespace

// Create a convex polyhedron from a set of 3D points. 
// The function computes the convex hull of the points using the QuickHull algorithm, merges coplanar triangles into polygonal faces, and returns a compact representation of the polyhedron.
ConvexPolyhedron make_convex_polyhedron(MatX points)
{
    using Point = std::array<Real, 3>;
    if (points.cols() < 4) {
        throw std::runtime_error("point set contains fewer than four points");
    }

    std::vector<Point> input(static_cast<std::size_t>(points.cols()));
    for (Eigen::Index i = 0; i < points.cols(); ++i) {
        input[static_cast<std::size_t>(i)] = {
            points(0, i), points(1, i), points(2, i)
        };
    }

	const Real tolerance = Real(1e-10) * geometry_scale(points); // normalize tolerance based on geometry scale
    quick_hull<std::vector<Point>::iterator, Real> hull{ 3, tolerance };
    hull.add_points(input.begin(), input.end());
    auto basis = hull.get_affine_basis();
    if (basis.size() < 4) {
        throw std::runtime_error("point set is not 3-D");
    }
    hull.create_initial_simplex(basis.begin(), std::prev(basis.end()));
    hull.create_convex_hull();

	const Vec3 interior = points.rowwise().mean(); // Use the mean of the input points as an interior point for orientation checks
    std::vector<std::array<std::uint32_t, 3>> triangles;
    triangles.reserve(hull.facets_.size());

	// Ensure that the triangles are oriented consistently with respect to the interior point
    for (const auto& facet : hull.facets_) {
        std::array<std::uint32_t, 3> triangle{};
        for (int i = 0; i < 3; ++i) {
            triangle[static_cast<std::size_t>(i)] = static_cast<std::uint32_t>(
                std::distance(input.begin(), facet.vertices_[i]));
        }
        const Vec3 a = points.col(triangle[0]);
        const Vec3 b = points.col(triangle[1]);
        const Vec3 c = points.col(triangle[2]);
        if ((b - a).cross(c - a).dot(a - interior) < Real(0)) {
            std::swap(triangle[1], triangle[2]);
        }
        triangles.push_back(triangle);
    }

    std::vector<Face> faces = merge_coplanar_triangles(points, triangles, tolerance);
    return compact_polyhedron(std::move(points), std::move(faces));
}

ConvexPolyhedron clip_convex_polyhedron(
    const ConvexPolyhedron& polyhedron,
    const Vec3& outwardNormal,
    Real offset)
{
	const Vec3 normal = outwardNormal.normalized(); // normal unit vector
	const Real tolerance = Real(1e-10) * geometry_scale(polyhedron.vertices); // normalize tolerance based on geometry scale
    std::vector<Vec3> vertices;
    vertices.reserve(static_cast<std::size_t>(polyhedron.vertices.cols()) * 2);
    std::vector<Face> faces;
    faces.reserve(polyhedron.faces.size() + 1);
    std::set<std::uint32_t> cutVertices;

	// Clip each face of the polyhedron against the plane defined by the normal and offset.
    for (const Face& sourceFace : polyhedron.faces) {
        Face clippedFace;
        clippedFace.reserve(sourceFace.size() + 1);
        for (std::size_t i = 0; i < sourceFace.size(); ++i) {
            const Vec3 a = polyhedron.vertices.col(sourceFace[i]);
            const Vec3 b = polyhedron.vertices.col(sourceFace[(i + 1) % sourceFace.size()]);
            const Real da = normal.dot(a) - offset;
            const Real db = normal.dot(b) - offset;
            const bool insideA = da <= tolerance;
            const bool insideB = db <= tolerance;

            if (insideA) {
                clippedFace.push_back(append_vertex(vertices, a, tolerance));
            }
            if (insideA != insideB) {
                const Real fraction = da / (da - db);
                const std::uint32_t intersection = append_vertex(
                    vertices, a + fraction * (b - a), tolerance);
                clippedFace.push_back(intersection);
                cutVertices.insert(intersection);
            }
        }
        remove_consecutive_duplicates(clippedFace);
        if (clippedFace.size() >= 3) {
            faces.push_back(std::move(clippedFace));
        }
    }

	// If the clipping plane intersects the polyhedron, create a new face from the intersection vertices. 
    // The vertices are sorted in a counter-clockwise order around the center of the intersection polygon.
    if (cutVertices.size() >= 3) {
        Vec3 center = Vec3::Zero();
        for (const std::uint32_t index : cutVertices) center += vertices[index];
        center /= static_cast<Real>(cutVertices.size());
        const Vec3 axisU = normal.unitOrthogonal();
        const Vec3 axisV = normal.cross(axisU);

        Face cutFace(cutVertices.begin(), cutVertices.end());
        std::sort(cutFace.begin(), cutFace.end(), [&](std::uint32_t lhs, std::uint32_t rhs) {
            const Vec3 a = vertices[lhs] - center;
            const Vec3 b = vertices[rhs] - center;
            return std::atan2(axisV.dot(a), axisU.dot(a))
                < std::atan2(axisV.dot(b), axisU.dot(b));
        });
        const MatX vertexMatrix = matrix_from_vertices(vertices);
        remove_collinear_vertices(cutFace, vertexMatrix, tolerance);
        if (face_normal(vertexMatrix, cutFace).dot(normal) < Real(0)) {
            std::reverse(cutFace.begin(), cutFace.end());
        }
        faces.push_back(std::move(cutFace));
    }

    MatX vertexMatrix = matrix_from_vertices(vertices);
    for (Face& face : faces) {
        remove_consecutive_duplicates(face);
        remove_collinear_vertices(face, vertexMatrix, tolerance);
    }
    faces.erase(std::remove_if(faces.begin(), faces.end(),
        [](const Face& face) { return face.size() < 3; }), faces.end());
    return compact_polyhedron(std::move(vertexMatrix), std::move(faces));
}

// Compute the volume of a convex polyhedron using the divergence theorem. The volume is calculated as the sum of the signed volumes of tetrahedra formed by each face and the origin.
Real convex_polyhedron_volume(const ConvexPolyhedron& polyhedron)
{
    Real volume = Real(0);
    for (const Face& face : polyhedron.faces) {
        const Vec3 a = polyhedron.vertices.col(face[0]);
        for (std::size_t i = 1; i + 1 < face.size(); ++i) {
            const Vec3 b = polyhedron.vertices.col(face[i]);
            const Vec3 c = polyhedron.vertices.col(face[i + 1]);
            volume += a.dot(b.cross(c)) / Real(6);
        }
    }
    return std::abs(volume);
}

PolyhedronMassProperties convex_polyhedron_mass_properties(
    const ConvexPolyhedron& polyhedron)
{
    Real volume = Real(0);
    Vec3 firstMoment = Vec3::Zero();
    Real Ixx = Real(0), Iyy = Real(0), Izz = Real(0);
    Real Ixy = Real(0), Iyz = Real(0), Izx = Real(0);

	// Compute the volume, first moment, and inertia tensor of the polyhedron by summing contributions from each face. Each face is triangulated into tetrahedra with the origin, and the contributions are calculated using standard formulas for tetrahedra.
    for (const Face& face : polyhedron.faces) {
        const Vec3 A = polyhedron.vertices.col(face[0]);
        for (std::size_t i = 1; i + 1 < face.size(); ++i) {
            const Vec3 B = polyhedron.vertices.col(face[i]);
            const Vec3 C = polyhedron.vertices.col(face[i + 1]);
            const Real Vt = A.dot(B.cross(C)) / Real(6);
            volume += Vt;
            firstMoment += Vt * (A + B + C) / Real(4);

            const Real Ax = A.x(), Ay = A.y(), Az = A.z();
            const Real Bx = B.x(), By = B.y(), Bz = B.z();
            const Real Cx = C.x(), Cy = C.y(), Cz = C.z();
            const Real xx = Ax * Ax + Bx * Bx + Cx * Cx + Ax * Bx + Bx * Cx + Cx * Ax;
            const Real yy = Ay * Ay + By * By + Cy * Cy + Ay * By + By * Cy + Cy * Ay;
            const Real zz = Az * Az + Bz * Bz + Cz * Cz + Az * Bz + Bz * Cz + Cz * Az;
            const Real xy = Real(2) * (Ax * Ay + Bx * By + Cx * Cy)
                + Ax * By + Bx * Ay + Bx * Cy + Cx * By + Cx * Ay + Ax * Cy;
            const Real yz = Real(2) * (Ay * Az + By * Bz + Cy * Cz)
                + Ay * Bz + By * Az + By * Cz + Cy * Bz + Cy * Az + Ay * Cz;
            const Real zx = Real(2) * (Az * Ax + Bz * Bx + Cz * Cx)
                + Az * Bx + Bz * Ax + Bz * Cx + Cz * Bx + Cz * Ax + Az * Cx;

            Ixx += Vt / Real(10) * (yy + zz);
            Iyy += Vt / Real(10) * (zz + xx);
            Izz += Vt / Real(10) * (xx + yy);
            Ixy += Vt / Real(20) * xy;
            Iyz += Vt / Real(20) * yz;
            Izx += Vt / Real(20) * zx;
        }
    }

	// Ensure that the volume is positive. If the computed volume is negative, it indicates that the face orientations are inconsistent, and we need to negate the volume and first moment to correct the orientation. The inertia tensor components are also negated accordingly.
    if (volume < Real(0)) {
        volume = -volume;
        firstMoment = -firstMoment;
        Ixx = -Ixx; Iyy = -Iyy; Izz = -Izz;
        Ixy = -Ixy; Iyz = -Iyz; Izx = -Izx;
    }
    if (!(volume > Real(0))) {
        throw std::runtime_error("convex polyhedron has non-positive volume");
    }

    const Vec3 center = firstMoment / volume;
    Ixx -= volume * (center.y() * center.y() + center.z() * center.z());
    Iyy -= volume * (center.z() * center.z() + center.x() * center.x());
    Izz -= volume * (center.x() * center.x() + center.y() * center.y());
    Ixy -= volume * center.x() * center.y();
    Iyz -= volume * center.y() * center.z();
    Izx -= volume * center.z() * center.x();

    Mat3 inertia;
    inertia <<
        Ixx, -Ixy, -Izx,
        -Ixy, Iyy, -Iyz,
        -Izx, -Iyz, Izz;
    return PolyhedronMassProperties{ volume, center, inertia };
}
