#pragma once
// Chart of a sphere-like piece: (a, b) -> surface point.
//
// a in [0, 1] runs from pole A (a = 0) to pole B (a = 1); b is periodic with period
// 1 around the pole axis. The chart is the piecewise-linear map given by the
// per-corner (a, b) of the triangles (TriMesh::uv): inside a triangle, the point is
// the barycentric combination of its corners, so the chart is exactly as smooth as
// the (subdivided) mesh. Normals are the barycentric interpolation of the vertex
// normals (smooth shading normals, not facet normals).
//
// Domain: the chart is a valid parameterization only on a in [aMin, aMax], the band
// between the first and the last ring of the original sphere grid (aMin = 1/(R+1),
// aMax = R/(R+1)). Inside the two pole caps, subdivision with linear UVs leaves the
// pole fan with b discontinuous along its radial edges (each original pole triangle
// carried its own pole-corner b), so (a, b) is not single-valued there; the caps are
// covered by the pole constructions, not by the chart.
//
// Checks (all on the stored data, no sampling tolerance beyond the stated ones):
//   orientation  every chart triangle is nondegenerate in (a, b) and has the same
//                orientation sign there, and that sign agrees with the outward 3D
//                orientation (d p/da x d p/db . n of one sign);
//   seam         across the seam, the two sides of each shared edge agree: per
//                vertex in the domain the corner b values agree modulo 1, and
//                evalRaw(a, 0) == evalRaw(a, 1) in position and normal on a sweep
//                of a over the domain;
//   monotone a   a is a vertex function (all corners agree), a = 0 / 1 only at the
//                two poles, and every other vertex has a neighbour with smaller and
//                one with larger a (no interior extremum of a).
#include <string>
#include <vector>

#include "model/status.h"
#include "model/mesh.h"
#include "model/spatial_hash.h"

namespace knit {

struct ChartPoint {
    bool ok = false;
    int tri = -1;
    double bary[3] = {0, 0, 0};
    Vec3 p, n;          // point, unit shading normal
    Vec3 dpda, dpdb;    // constant per triangle
    Vec2 ab;            // chart coordinates (b as stored on that triangle)
    double distance = 0;  // closestPoint only
};

struct ChartCheck {
    Status orientation = Status::Unresolved;
    Status seam = Status::Unresolved;
    Status monotoneA = Status::Unresolved;
    int degenerateTriangles = 0, flippedTriangles = 0, flipped3D = 0;
    int seamVertexMismatches = 0;
    double seamMaxPositionGap = 0, seamMaxNormalGap = 0;
    int seamSamples = 0;
    int aCornerMismatches = 0, aExtrema = 0, aPoles = 0;
    int domainGaps = 0;  // sweep points in the domain not covered by any triangle
    std::string detail;
    Status overall() const { return combine(combine(orientation, seam), monotoneA); }
};

class Chart {
public:
    // rings = R of the original sphere grid (manifest); it fixes the domain.
    Chart(const TriMesh& mesh, int rings);

    const TriMesh& mesh() const { return mesh_; }
    int rings() const { return rings_; }  // as constructed (relax-call dumps rebuild the chart)
    // Point at chart coordinates; b is reduced modulo 1 first. ok = false if (a, b)
    // lies outside every chart triangle (should not happen for a in [0, 1]).
    ChartPoint eval(double a, double b) const;
    // As eval, without reducing b (to address either side of the seam).
    ChartPoint evalRaw(double a, double b) const;
    // Closest point on the mesh, with interpolated normal and chart coordinates.
    ChartPoint closestPoint(const Vec3& p) const;

    double aMin() const { return aMin_; }
    double aMax() const { return aMax_; }
    Vec3 poleA() const { return mesh_.pos[poleA_]; }
    Vec3 poleB() const { return mesh_.pos[poleB_]; }

    // Length of the level set {a = const} on the mesh.
    double ringLength(double a) const;
    // Longest level set over a in [0, 1] (scan, then golden-section refinement).
    struct RingAt {
        double a, length;
    };
    RingAt longestRing() const;
    // Arc length of the chart curve b = const from a0 to a1 (n samples). Portions in a
    // pole cap are measured as the chord from the pole to the cap boundary.
    double meridianLength(double b, double a0, double a1, int n = 256) const;

    ChartCheck check() const;

    // Read-only view of the (a, b) triangle grid for the device's copy of evalRaw.
    struct GridView {
        int na = 0, nb = 0;
        double b0 = 0, b1 = 1;
        const std::vector<std::vector<int>>* cells = nullptr;  // na x nb, row-major in a
    };
    GridView gridView() const { return {na_, nb_, b0_, b1_, &grid_}; }

private:
    TriMesh mesh_;
    int rings_ = 0;
    int poleA_ = -1, poleB_ = -1;
    double aMin_ = 0, aMax_ = 1;
    // 2D grid over the chart triangles in (a, b).
    int na_ = 1, nb_ = 1;
    double b0_ = 0, b1_ = 1;
    std::vector<std::vector<int>> grid_;
    std::vector<Vec3> dpda_, dpdb_;
    SpatialHash3 hash_;
    double hashCell_ = 1;

    void buildTriangleGrid();
    void buildTriangleHash();
    void checkOrientation(ChartCheck& c) const;
    void checkVertexValues(ChartCheck& c) const;
    void checkSeam(ChartCheck& c) const;
    bool inTriangle(int t, double a, double b, double w[3], double tol) const;
    ChartPoint at(int t, const double w[3]) const;
};

// How far a chart's columns (levels of b) are from perpendicular to its rows (levels of a), in degrees: per
// triangle the angle between grad a and grad b, off 90, over the triangles that touch no pole, area-weighted.
// A hand-made UV-sphere grid has its columns on meridians; a skewed chart shears the stitches.
struct ChartSkew {
    double mean = 0, p95 = 0;
};
ChartSkew chartSkew(const TriMesh& charted);

} // namespace knit
