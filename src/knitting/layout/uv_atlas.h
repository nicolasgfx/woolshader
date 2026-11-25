#pragma once
// The UV atlas of a welded triangle mesh: its charts, each with a metric that makes uv isometric on
// average, and which triangles of a chart are usable.
//
// Charts: triangles joined across edges whose two corner uvs agree on both sides (within
// uvTol). Per triangle the Jacobian J_t (3 x 2) maps uv to 3D. Per chart:
//   metric      M = Gbar^(1/2), Gbar = sum A_t J_t^T J_t / sum A_t over usable triangles;
//               uv' = M (uv - c) is isometric ON AVERAGE (mm). A chart whose majority
//               orientation is mirrored gets x' flipped (frames right-handed w.r.t. the
//               outward normal).
//   distortion  both singular values s1 >= s2 of K_t = J_t M^-1 per triangle (reported:
//               s1, s2, s1 / s2, s1 s2); the bounds [0.75, 1.3] apply to s1 and s2
//               individually, so they catch residual LOCAL stretch, isotropic too. A
//               uniform scale of the whole chart cancels in K_t (J = 2 J0 -> M = 2 M0).
//   validity    degenerate (|det uv| < 1e-12 x chart mean, or rank-deficient J),
//               folded (orientation against the chart's area majority), overlapping
//               (uv' interiors of non-adjacent triangles intersect): unusable; a chart
//               with folded + overlap area > 1 % is UNRESOLVED (not tiled).
//   orientation world up projected into each triangle's tangent plane and pulled back
//               into uv' (d = M J_t^+ u_t), normalized, area-weighted mean in the
//               common chart coordinates; mean length < 0.2: the uv' v axis (recorded).
#include <array>
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/status.h"
#include "model/vec.h"

namespace knit {

struct Mat2 {
    double a = 1, b = 0, c = 0, d = 1;  // [[a, b], [c, d]]
    Vec2 operator*(const Vec2& v) const { return {a * v.x + b * v.y, c * v.x + d * v.y}; }
    Mat2 operator*(const Mat2& o) const { return {a * o.a + b * o.c, a * o.b + b * o.d, c * o.a + d * o.c, c * o.b + d * o.d}; }
    double det() const { return a * d - b * c; }
    Mat2 inverse() const { const double k = 1.0 / det(); return {d * k, -b * k, -c * k, a * k}; }
};
// Square root of a symmetric positive definite 2 x 2 matrix (closed form).
Mat2 sqrtSpd(const Mat2& g);

struct Distortion {  // area-weighted quantiles {min, p5, p50, p95, max} over the usable triangles
    double s1[5] = {}, s2[5] = {}, ratio[5] = {}, areaScale[5] = {};
    double outsideFraction = 0;  // area fraction with s1 or s2 outside [lo, hi]
};

struct UvChart {
    std::vector<int> tris;
    Mat2 M;                    // uv -> uv' (mirror included)
    Vec2 centre;               // uv centroid c
    bool mirrored = false;
    Status status = Status::Pass;  // UNRESOLVED: folded + overlap area > 1 %
    double area = 0, degenerateArea = 0, foldedArea = 0, overlapArea = 0;  // 3D mm^2
    Distortion distortion;
    Vec2 upDir;                // unit wale base direction in uv' (before the user angle)
    bool upFromWorld = true;   // false: the uv' v-axis fallback
    double upMeanLength = 0;
    std::string detail;
};

class UvAtlas {
public:
    // m: welded, per-corner uv; up: world up in the mesh's frame (library +Z).
    UvAtlas(const TriMesh& m, const Vec3& up = {0, 0, 1}, double uvTol = 1e-7, double lo = 0.75, double hi = 1.3);
    const TriMesh& mesh() const { return m_; }
    const std::vector<UvChart>& charts() const { return charts_; }
    bool usable(int t) const { return usable_[t]; }
    // uv' of the corners of triangle t (its chart's map).
    std::array<Vec2, 3> uvp(int t) const;
    // K_t = J_t M^-1 applied to a uv' vector (3D tangent vector, mm per uv' mm).
    Vec3 pushForward(int t, const Vec2& d) const;
    // Usable triangles of chart c containing uv' point p (0: outside, 2+: ambiguous).
    std::vector<int> locate(int chart, const Vec2& p, double bary[3]) const;
    // Edges of chart c's usable region (mesh boundary, uv seams, edges to unusable
    // triangles), in uv'.
    const std::vector<std::array<Vec2, 2>>& boundary(int chart) const { return boundary_[chart]; }
    // The part of boundary(chart) on open edges of the mesh (edges of one triangle: no surface
    // beyond them, unlike a seam), in uv'.
    const std::vector<std::array<Vec2, 2>>& openBoundary(int chart) const { return openBoundary_[chart]; }
    // Smallest s2 over the usable triangles of chart c whose uv' boxes meet [lo, hi].
    double minS2In(int chart, const Vec2& lo, const Vec2& hi) const;

    // Read-only views for the device's copy of the atlas (the lattice's node tests).
    struct GridView {
        Vec2 lo;
        double cell = 1;
        int nx = 0, ny = 0;
        const std::vector<std::vector<int>>* cells = nullptr;  // ny x nx, indices into charts()[c].tris
    };
    GridView gridView(int chart) const {
        const Grid& g = grid_[chart];
        return {g.lo, g.cell, g.nx, g.ny, &g.cells};
    }
    const std::vector<double>& s2() const { return s2_; }

private:
    // The constructor's stages, in order (chart by chart from chartMetric on).
    struct Construction;
    void joinCharts(Construction& b);
    void measureTriangles(Construction& b);
    bool chartMetric(int c, Construction& b);
    void buildChartGrid(int c, Construction& b);
    void chartValidity(int c, Construction& b);
    void chartDistortion(int c, Construction& b);
    void chartOrientation(int c, Construction& b);
    void collectBoundaries(int c, Construction& b);
    // Singular values s1 >= s2 of K_t.
    struct SingularValues {
        double s1 = 0, s2 = 0;
    };
    SingularValues singularValues(int t) const;

    const TriMesh& m_;
    std::vector<int> chartOf_;
    std::vector<char> usable_;
    std::vector<std::array<double, 6>> J_;  // column-major 3 x 2 (dP/du, dP/dv)
    std::vector<double> s1_, s2_;           // singular values of K_t (usable triangles)
    std::vector<UvChart> charts_;
    std::vector<std::vector<std::array<Vec2, 2>>> boundary_;
    std::vector<std::vector<std::array<Vec2, 2>>> openBoundary_;
    struct Grid { Vec2 lo; double cell = 1; int nx = 0, ny = 0; std::vector<std::vector<int>> cells; };
    std::vector<Grid> grid_;
};

} // namespace knit
