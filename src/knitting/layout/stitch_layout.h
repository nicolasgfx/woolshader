#pragma once
// Stitch layout and yarn placement on a base surface: where each stitch goes, the frame around it,
// and the yarn control points that follow from them.
//
// Coordinate contract (see stockinette.cpp): a stitch template control (X, Y, Z) has
// its origin at the stitch's leg/foot level (X = 0, Y = 0; crown Y = 1.40, sinkers
// Y = -0.36, X = +-0.5). +Y points from the parent towards the stitch's own head (up
// the wale), +X along the course in knitting order, Z along the outward normal in
// units of A * D. The LAYOUT NODE of a stitch is its template origin; the parent's
// node sits one row (about h) below along the wale.
//
// Frame map of a node k (physical units, not chart units):
//   F_k(X, Y) = node_k + X * (X < 0 ? wL_k : wR_k) * c_k + Y * h_k * v_k
// (tangent plane at the node; c the unit course direction towards the next yarn
// neighbour, v the unit wale direction away from the parent, skew kept; wL, wR the
// distances to the previous / next node, so X = +-0.5 is the boundary point shared
// with the neighbour; h the row height). A control is placed at
//   s, n = closest point of p0 on the base surface,   control = s + Z * A * D * n.
//
// Every interlock lives in the CHILD's frame: a stitch's lower part (feet, legs:
// Y <= 0.3) is mapped in its own frame, its head (Y >= 0.95) in the frame of the
// stitch(es) drawn through it, one row down (child Y = Y - 1); in between (neck) the
// two are blended smoothly, so a loop whose children sit sideways leans instead of
// kinking. Heads therefore follow the child row's spacing and neighbouring heads meet
// at the children's shared boundaries:
//   K child        p0 = F_c(X, Y - 1)
//   INC children   the head is mapped onto both children's interval (2x stretch):
//                  X < 0 -> F_c0(2X + 0.5, Y - 1), X >= 0 -> F_c1(2X - 0.5, Y - 1);
//                  the crown sits on the boundary between the two children
//   DEC child      both parents' heads in F_c, shifted by -/+ decShift along X and
//                  layered in depth: the later parent (p1) in front by +decDepth, the
//                  earlier (p0) behind by -decDepth (k2tog: p1 lies on top); the DEC
//                  child's lower part is pushed apart in depth by decLegDepth so it
//                  clears both layers
// A stitch without children (the last turn) keeps its head in its own frame. On a
// regular planar layout every frame coincides with the flat generator's, so the
// mapping reproduces it exactly (test planar_mapping_equals_flat).
#include <string>
#include <vector>

#include "knitting/layout/base_surface.h"
#include "knitting/stitch/stitch_program.h"
#include "knitting/stitch/stockinette.h"
#include "model/params.h"
#include "model/status.h"
#include "model/yarn_path.h"

namespace knit {

struct GaugeParams {
    double w = 5.5, h = 3.4;  // stitch width and row height [mm]
};

struct ShapingParams {
    double decShift = 0.1;     // DEC parent heads: X shift in child widths
    double decDepth = 1.0;     // DEC parent heads: depth offset in units of A * D
    double decLegDepth = 0.6;  // DEC child lower part: relative depth stretch
    double headBlend0 = 0.3, headBlend1 = 0.6;  // own -> children's frame between these template Y
    double decBlend0 = 0.3, decBlend1 = 0.6;    // DEC parents: own -> child frame
    double decDepth0 = 0.15, decDepth1 = 0.45;  // DEC parents: depth layering ramp
    double closeY = 0.1;       // CLOSE strand height in the last turn's virtual child frames
    double closeDepth = 1.8;   // CLOSE strand depth swing, in units of A * D
};

// One node per stitch (template origin) plus the local frame and size.
struct LayoutNode {
    Vec3 p, n;       // node on the surface, unit normal
    Vec3 c, v;       // unit course and wale directions in the tangent plane (skew kept)
    double w = 0, h = 0;
    double wL = 0, wR = 0;  // distances to the previous / next node (X < 0 / X >= 0 scale)
    double a = 0, b = 0;  // chart coordinates of the initial placement
    double f = 0;         // fraction of the turn (arc-length fraction of the row when given)
};

struct PieceLayout {
    std::vector<LayoutNode> nodes;  // by stitch id
    LayoutNode before, after;       // virtual neighbours for the padding controls
    bool helical = false;           // spiral: course neighbours are yarn neighbours across turns
    double aSign = 1;               // direction of increasing row in a (turn-0 wale direction)
    // Virtual children of the last turn (one row further, not knitted): its heads are
    // mapped into them as into K children. virtualChildOf[k] indexes virtualNodes or -1.
    std::vector<LayoutNode> virtualNodes;
    std::vector<int> virtualChildOf;
    // Lattice layouts (chart_lattice.h); all empty for spiral and flat programs:
    std::vector<int> rowPrev, rowNext;  // course neighbours by stitch id (-1: run end); else same-turn ids
    std::vector<Vec3> waleHint;         // parentless stitches: wale direction (else the chart's at())
    std::vector<ProjState> nodeState;   // each node's accepted surface location (placement starts there)
};

// Position of each stitch along its turn, as a fraction f of the turn (1D layout
// relaxation, the course-direction part of Yuksel et al. 2012, §6.1): course springs
// towards uniform spacing 1/n_r, wale springs (weight `wale`) towards the stitch's
// parent(s) (K: above the parent, DEC: above the midpoint of its parents, INC: half a
// stitch to either side). Turn 0 is uniform.
std::vector<double> turnFractions(const StitchProgram& prog, double wale = 0.5, int iterations = 400);

// Where one stitch of a spiral sits on the chart.
struct HelixSpot {
    double a = 0, b = 0;  // chart coordinates
    double f = 0;         // fraction of the turn after any warp (LayoutNode::f)
};

// Where a spiral's stitches sit on the chart, worked out by the caller before the surface is asked:
// stitch k of turn r at row coordinate rho = r + f_k (f_k from turnFractions); the virtual neighbours
// before the first and after the last stitch, where the padding controls go; and the virtual row after
// the last turn, the frames the last turn's heads are mapped into (a gathered end's ring round the end
// pole, or one row on for a threaded end; see PieceLayout).
struct HelixNodes {
    std::vector<HelixSpot> stitches;     // by stitch id
    HelixSpot before, after;
    std::vector<HelixSpot> virtualRow;   // the last turn's virtual children, in its order; empty: no virtual row
    double aSign = 1;                    // +1: the rows rise in a (PieceLayout::aSign)
};

// The points and frames of a helix on the surface: the stitches and their two neighbours in one batch
// and computeFrames, then the virtual row in a second batch with its frames.
PieceLayout helicalLayout(const BaseSurface& s, const StitchProgram& prog, const GaugeParams& g, const HelixNodes& nodes);

// The plain helix of a surface of revolution: row rho at a = (rho + 1) rowHeight / length, b = f, with f
// from turnFractions; no virtual row.
HelixNodes evenHelix(const StitchProgram& prog, double rowHeight, double length, double walePull = 0.5);

// Frames (c, v, w, h) from the node positions: c towards the yarn neighbours, v away
// from the parent(s) (turn 0: along the row direction, aSign * increasing a).
void computeFrames(const BaseSurface& s, const StitchProgram& prog, const GaugeParams& g, PieceLayout& lay);

// One open yarn through all stitches in program order (16 controls per stitch, one
// padding control before and two after, as the flat generator), with stitch spans,
// topology (records with frames) and kind = Surface.
// With prog.closeDrawstring and a virtual row, the yarn continues after the last
// stitch's sinker as the CLOSE drawstring: through every loop of the last turn, in
// order, alternately back to front (even loops) and front to back (odd loops), at
// height closeY of the virtual child frame, 5 controls per loop (depth -z, -z, 0, +z,
// +z or reversed). The strand changes depth inside the loops, never in the narrow
// gaps between neighbouring heads. A CLOSE terminal record with the per-loop u-ranges
// is added.
// The yarn comes without its arc-length tables: the caller finishes it (KnitNumerics::store on
// the device, finalizeCurves on the CPU). Nothing here reads them; uEnd needs only the controls.
CurveSet placeSpiral(const BaseSurface& s, const PieceLayout& lay, const StitchProgram& prog, const PatternParams& pp,
                     const ShapingParams& sh = {});

struct RowPlacement {
    ShapingParams sh;
    CurveKind kind = CurveKind::FlatPatch;           // Surface for lattices on meshes: local signature
                                                     // instead of the global projection
    const std::vector<StitchFaces>* faces = nullptr;  // patterns; padding controls use the end stitch's
                                                     // faces with no wale transition
    SurfaceBinding* bind = nullptr;                  // the projection binding (sized and filled here)
};

// One open yarn PER TURN (flat test patches: no helix), all knitted in +X, kind =
// FlatPatch so the global projected checks apply. Padding uses the turn's end frames.
// Without arc-length tables, as placeSpiral.
CurveSet placeRows(const BaseSurface& s, const PieceLayout& lay, const StitchProgram& prog, const PatternParams& pp,
                   const RowPlacement& rows = {});

// Every yarn's arc-length tables on the CPU (fixtures and tests).
void finalizeCurves(CurveSet& cs);

// Layout check (part of gate a), per-op cell dimensions:
//   K, INC child, DEC child   width (wL + wR) / 2 vs w; height (to the parent, the INC
//                             slot or the DEC parents' midpoint) vs h
//   INC parent                head span (its two children's outer boundaries) vs 2 w
// Bounds: every ratio in [lo, hi], the angle between course and wale in [60, 120]
// degrees, and no inverted cell (c x v . n > 0). Turn 0 has no height (no parent).
struct LayoutCheck {
    Status status = Status::Pass;
    int cells = 0, outOfBounds = 0, skewed = 0, inverted = 0;
    double minW = 1e300, maxW = 0, minH = 1e300, maxH = 0, minAngle = 180, maxAngle = 0;
    int worst = -1;  // stitch with the largest relative deviation
    std::string detail;
};
LayoutCheck checkLayout(const StitchProgram& prog, const PieceLayout& lay, const GaugeParams& g, double lo = 0.75,
                        double hi = 1.3);

// Regular planar layout of a program: turn r at y = r h, stitch j at
// x = (j - (n_r - 1) / 2) w (rows centred on x = 0).
PieceLayout planarLayout(const StitchProgram& prog, const GaugeParams& g);

} // namespace knit
