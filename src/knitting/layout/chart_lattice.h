#pragma once
// The stitch lattice on the UV charts of a welded triangle mesh: which stitches exist on a chart, how
// they group into yarns, and their program, layout and placed curves.
//
// Per chart (uv_atlas.h) a regular lattice of stitch nodes is laid in the isometric-on-average uv'
// plane: wale direction = the chart's base direction rotated by the user angle, spacing w x h in mm.
// A node is KEPT iff
//   - it lies in exactly one usable triangle of the chart (outside / ambiguous dropped),
//   - its 3D point is not inside another eligible closed component (TileParams::solids),
//   - its stitch footprint clears the chart's usable-region boundary (pre-filter, below),
// and it belongs to a run of >= minRun consecutive kept nodes of its row. A run is one open yarn
// (turn = run); rows are not continued across charts (unsewn seams). Parents, children and the
// pattern face are functions of the lattice coordinates (chart, i, r) only: the parent of (i, r) is
// (i, r - 1) of the same chart if kept, whichever run it belongs to.
//
// Footprint pre-filter (placement, NOT a guarantee): the template's in-plane extent (X +-0.5 w; run
// ends +-Xpad w from the padding controls; Y from the template's min to max height), inflated by
// m = R_s + c/2 converted into uv' with the smallest singular value s2 over the triangles under the
// inflated rectangle (iterated), must not meet the boundary of the chart's usable region. The
// guarantee is the check of settled yarn, B1 (certified) and B2 (sampled), see BoundaryCheck.
#include <string>
#include <vector>

#include "knitting/layout/base_surface.h"
#include "knitting/layout/stitch_layout.h"
#include "knitting/numerics.h"
#include "knitting/layout/uv_atlas.h"
#include "model/mesh_query.h"

namespace knit {

// A chart's candidate lattice nodes (r-major: (r, i) at index r * ni + i), classified in one batch
// (KnitNumerics::classifyLattice).
struct LatticeChartQuery {
    int chart = -1;
    Vec2 c2, v2;             // course / wale directions in uv'
    double xmin = 0, ymin = 0, w = 0, h = 0;
    int ni = 0, nr = 0;
    double margin = 0, Ymin = 0, Ymax = 0, Xpad = 0;  // the footprint test
    bool footprintCheck = true;
};
struct LatticeNodeClass {
    int code = 1;            // 0 kept, 1 outside, 2 ambiguous, 3 inside another component, 4 footprint
    int tri = -1;
    double bary[3] = {0, 0, 0};
    Vec3 p;
    // Kept nodes: the run-end footprint tests (padding extents at a run's start / end); -1: not
    // evaluated (tileMesh evaluates them when needed).
    signed char clearStart = -1, clearEnd = -1;
};

// One chart's lattice frame: course / wale directions in uv', the first node (x0, y0), at the
// chart's lattice-frame minimum, and the node counts covering the usable triangles. ok = false: no
// usable triangle.
struct LatticeFrame {
    bool ok = false;
    Vec2 c2, v2;
    double x0 = 0, y0 = 0;
    int ni = 0, nr = 0;
};
struct LatticeSpacing {
    double angleDeg;  // the user orientation
    double w, h;      // gauge [mm]
};
LatticeFrame latticeFrame(const UvAtlas& A, int chart, const LatticeSpacing& spacing);

// Closed components a stitch may not lie inside. Their frame is not the part's: library point +
// offset = the solids' frame. Tested in one batch by the classification, and on the host for the
// nodes that surface walks reach.
struct LatticeSolids {
    std::vector<const MeshQuery*> meshes;
    Vec3 offset;
    bool contains(const Vec3& p) const;
};

struct TileParams {
    PatternParams pp;          // D, c, A, pattern, template; spacings are set from g
    GaugeParams g;
    double angleDeg = 0;       // user orientation, added to each chart's base direction
    Vec3 up{0, 0, 1};
    int minRun = 4;
    double supportRadius = 0;  // R_s [mm] (FinalSurfaceSpec); <= 0: D / 2
    ShapingParams sh;
    ContinuationParams cp;     // lambdaW <= 0: 0.5 w
    LatticeSolids solids;      // none: no culling
    // false skips the footprint pre-filter: overlap parts and gap patches (scene/knit/mesh_part.cpp), whose
    // stitches run past their chart on purpose, and the tests' seam negative control. A zero inset alone
    // does not make two charts collide: each chart's lattice is anchored at its own boundary, so its first
    // stitch starts at least (1 - Xpad) w from a seam.
    bool footprintCheck = true;
    // Overlap coverage (knit.coverage = overlap): tile only chart onlyChart (-1: every chart), and
    // continue it `overhang` stitches / rows beyond its kept nodes by walks on the surface (course
    // steps of w, wale steps of h from the last node, projected by continueFrom; the frame is
    // transported into each new tangent plane). The walked nodes may lie on other charts: the
    // overhang overlaps the neighbouring charts' yarn (hidden under it), so the seams close without
    // sewing. Columns are continued first (rows above and below each column's ends), then every
    // run at both ends. A walk stops at an open boundary, a failed projection, an existing lattice
    // node or inside another component.
    int onlyChart = -1;
    int overhang = 0;
    // Overlap coverage, gap patches: a local lattice walked from a seed point instead of the charts
    // (tri >= 0): wale along `wale` projected into the seed's tangent plane, course = wale x n; rows
    // walked from the seed along the wale, stitches along each row, both ways, keeping the nodes
    // within `radius` [mm] of the seed. Independent of the atlas (no chart is tiled).
    struct WalkPatch {
        int tri = -1;
        double bary[3] = {0, 0, 0};
        Vec3 wale;
        double radius = 0;
    } patch;
};

struct ChartTiling {
    int chart = -1;
    Status status = Status::Pass;  // the atlas chart's status (UNRESOLVED: not tiled)
    double baseAngleDeg = 0;       // of the uv' base direction
    bool upFromWorld = true;
    int candidates = 0, kept = 0;
    int outside = 0, ambiguous = 0, boundary = 0, insideOther = 0, shortRun = 0, runEnd = 0;  // drops by reason
    int overhang = 0;              // walked nodes beyond the chart (TileParams::overhang)
    bool skipped = false;          // not this part's chart (TileParams::onlyChart)
    int runs = 0, droppedRuns = 0;
    double coverage = 0;           // kept stitches x w x h / chart area
};

struct LatticeCheck {
    Status status = Status::Pass;
    int stitches = 0, parentless = 0, uncaught = 0;  // counts only (cast-on / bind-off omitted)
    std::string detail;            // first violation or summary
};
// FAIL if a parent is not at (i, r - 1) of the same chart, a loop has more than one child,
// a run is not contiguous in i, or ids are not in yarn order.
LatticeCheck checkLatticeProgram(const StitchProgram& prog);

struct LatticeKnit {
    StitchProgram prog;
    PieceLayout lay;
    std::vector<StitchFaces> faces;
    CurveSet curves;
    SurfaceBinding bind;
    std::vector<ChartTiling> charts;
    LatticeCheck lattice;
    LayoutCheck layout;
    PatternParams pp;              // with the gauge's spacings (the placement's: nominal D and c)
    std::string detail;
};

// The component a lattice is laid on (welded mesh with uv, library mm).
struct ChartedMesh {
    const UvAtlas& atlas;
    const MeshSurface& surface;  // must outlive the result (LatticeKnit::bind points to it)
};
// Lattice, program, layout and placement on one component (buildRuns, buildStitchProgram,
// layOutNodes, placeYarn); numerics classifies the nodes and stores the placed curves.
bool tileMesh(const ChartedMesh& on, const TileParams& p, KnitNumerics& numerics, LatticeKnit& out, std::string* err);

// The footprint test of one node: no border edge of the chart in the box
// [x - left, x + right] x [y + Ymin h - m, y + Ymax h + m] of the lattice frame; m is the margin in
// uv', grown by the chart's smallest stretch in the box. In the overlap mode only the chart's
// open edges count: the seams are crossed, an open edge is not.
struct NodeExtent {
    double x, y;         // the node in the lattice frame
    double left, right;  // the stitch's reach along the course
};
bool footprintClear(const UvAtlas& atlas, const LatticeChartQuery& q, const NodeExtent& node);

// The template's extents over every face combination (tileMesh's footprint): its height range, and
// how far the padding controls reach beyond a run's end nodes.
struct TemplateExtent {
    double Ymin, Ymax, Xpad;
};
TemplateExtent templateExtent(const PatternParams& pp);

} // namespace knit
