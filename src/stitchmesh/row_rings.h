#pragma once
// Method A's rings: the level
// sets t = level of a vertex field on a triangle mesh, as contours whose points keep their bindings on the mesh. A
// component supports Method A at a level only when the level is exactly one closed contour; anything else (several
// contours, an open one ending on a boundary or hole) is a contour event, reported, never repaired here.
#include <string>
#include <vector>

#include "model/mesh.h"

namespace knit {

// One crossing of the level with a mesh edge.
struct RingPoint {
    Vec3 p;
    int a = -1, b = -1;   // the edge (as stored in the triangle that the contour leaves through it)
    double f = 0;         // p = (1 - f) pos[a] + f pos[b]
    int tri = -1;         // the triangle the contour runs through after this point (-1: the contour ends here)
};

struct Contour {
    std::vector<RingPoint> points;   // oriented so that larger values lie on the left, seen from outside
    bool closed = false;             // the last point connects back to the first
    double length = 0;               // polyline length, closing segment included [mm]
};

struct RingLevel {
    double level = 0;
    std::vector<Contour> contours;
    int crossedTriangles = 0;
};

// Values equal to the level count as above it (one consistent rule, so a level through a ring of vertices is still one
// closed contour; its points may then lie on those vertices). Segments are oriented from the corners' classes, never
// from positions.
RingLevel extractLevel(const TriMesh& m, const std::vector<double>& t, double level);

// The longest closed contour beside a row that is knitted over instead of being a contour event, in stitch widths w:
// below 2.5 w a round of its own would have fewer than 3 stitches (round(L / w) <= 2), and 3 is the smallest round
// (the pole start's minimum, guided_rows.h).
constexpr double kKnitOverStitches = 2.5;

// A bump or a pit of the surface too small for a round of its own closes a level into a tiny extra contour beside the
// row, so the row is knitted over it. Removes from `level` every closed contour shorter than minLength (callers pass
// kKnitOverStitches * w) when another contour remains, puts the longest remaining contour first, and returns the
// removed ones. A level whose contours are all short, or none of them, is left as it is.
std::vector<Contour> dropContoursShorterThan(RingLevel& level, double minLength);

// Where a closed contour crosses the half-plane {x : dot(x - origin, normal) = 0, dot(x - origin, side) > 0}
// (the cut guide of a surface of revolution: the meridian half-plane at one azimuth). Returns the parameter
// i + f of the crossing along the contour's points (point i to point i + 1, f in [0, 1)), or -1 when it does not
// cross exactly once; crossings receives the count.
double cutCrossing(const Contour& c, const Vec3& origin, const Vec3& normal, const Vec3& side, int* crossings);

// A drawn row-start path: the polyline of a marking path, each segment inside one triangle of
// the mesh it was drawn on. When the rings run on a trimmed copy of that mesh (surface_trim.h), parentTri maps every
// ring-mesh triangle to the triangle it was split from (nullptr: the rings run on the drawn mesh itself).
struct SurfacePathGuide {
    std::vector<Vec3> points;
    std::vector<int> segmentTri;                 // segment i (points i to i + 1) lies in this triangle
    const std::vector<int>* parentTri = nullptr;
};

// Where a closed contour on m crosses the path, tested per triangle in that triangle's plane (a ring segment and a path
// segment in the same piece triangle). Same return as cutCrossing: the parameter of the one crossing, else -1.
// pathParameters (may be null) receives every crossing's place along the path: segment index + fraction.
double pathCrossing(const TriMesh& m, const Contour& c, const SurfacePathGuide& g, int* crossings, std::vector<double>* pathParameters = nullptr);
// The contour parameter closest to the path (a boundary ring: the path may start or end just beside it).
double pathNearest(const Contour& c, const SurfacePathGuide& g);

} // namespace knit
