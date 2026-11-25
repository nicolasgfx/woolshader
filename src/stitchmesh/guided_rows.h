#pragma once
// Method A of the stockinette experiment: the typed stitch mesh from
// the time field's rings.
//
//   rings     ring 0 is the cast-on loop, ring R the bind-off loop, ring k between them the level t = k h' (h' = T / R,
//             R = round(T / h), T = the field's smallest value on the bind-off loop), each one closed contour.
//   vertices  ring k gets n_k = round(L_k / w) vertices at arc fractions j / n_k from its crossing with the cut guide,
//             so vertex 0 of every ring is on the cut; each keeps its binding (triangle, barycentrics).
//   shaping   strip k (between ring k and ring k + 1) needs e = |n_{k+1} - n_k| shaping slots, INC when the count grows.
//             Their slots follow section 7: evenly spaced round the strip at fractions (i + sigma) / e, sigma staggered
//             between successive shaping strips by the golden ratio; never beside the shift-path (slot 0 and the last
//             slot); at least one K between two shaping slots; a slot that cannot be legalized within m = n_k / (2 e)
//             slots of its target, or a strip that needs shaping in the first or last strip, keeps its geometry and is
//             typed Unsupported with the reason code. Every placement is logged with its target and move.
//   faces     the slots between consecutive wale edges of a monotone matching that starts at vertex 0 on the cut; the
//             key (bottom count, top count, course roles, beside the path) selects the construction (lookupConstruction).
//   path      the shift-path is the chain of wale edges (ring k vertex 0, ring k + 1 vertex 0).
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/status.h"
#include "stitchmesh/row_rings.h"
#include "stitchmesh/typed_mesh.h"

namespace knit {

// The lookup's key (section 6) and its answer.
struct SlotKey {
    int bottom = 1, top = 1;                         // course edges shared with the ring below / above
    EdgeBoundary bottomRole = EdgeBoundary::Interior, topRole = EdgeBoundary::Interior;
    bool besidePath = false;                         // slot 0 or the last slot of its strip
    // Variant "boundary shaping": kp / d12k with a cast-on bottom or a bind-off
    // top edge are admitted ("INC kp + cast-on", "DEC d12k + bind-off", ...). Off: such a slot is Unsupported,
    // shaping-at-boundary-strip (the CP2 rule).
    bool boundaryShaping = false;
};
struct LookupResult {
    StitchFaceKind kind = StitchFaceKind::Unsupported;
    std::string reason;                              // the reason code when Unsupported
    std::string construction;                        // "K", "K + cast-on", "K + bind-off", "INC kp", "DEC d12k"
};
LookupResult lookupConstruction(const SlotKey& key);

// A row-start path point with the field's value there (pathOnField in marked_build.cpp).
struct MarkedPathPoint {
    Vec3 p;
    double t = 0;            // the field at p
    int segmentTri = -1;     // the triangle the segment from this point to the next lies in (-1 for the last point)
};

// Where each ring's vertices sit around its contour.
enum class RowPhase {
    // Vertex 0 where the ring crosses the row-start path, the others at equal arcs from it (CP2). The shaping rule puts
    // a strip's increases or decreases after slot 0, so every wale edge of the strip leans the same way: by 2 to 3 mm
    // per ring on the teddy body, which the relaxation then undoes by turning whole rings against each other (up to
    // 23 mm one way and 22 mm the other over the piece).
    OnPath,
    // As OnPath, then each ring from the second on is turned on its contour until its wale edges to the ring below lean
    // as much one way as the other. The row start leaves the path by the turns added up (GuidedRowsProduct::rowStartOffPath).
    WalesUpright,
};

// How many stitches the bind-off ring has.
enum class BindOffCount {
    // round(length of the bind-off loop / w), as every ring (CP2). Where the piece narrows towards its opening, the last
    // strip then holds a decrease: a pentagon two stitches wide whose top edge lies on the border. On both teddy legs that
    // one face is the only one over the face-deviation limit.
    FromLoopLength,
    // The count of the round below: the bind-off strip is plain, and its stitches are as much narrower or wider as the
    // loop is shorter or longer than that round. The local port of Wu 2018 (Method B) puts no decrease in the last row
    // of the same leg either, and all its decreases fit.
    OfRoundBelow,
};

struct GuidedRowsRequest {
    double w = 5.5, h = 5.5;
    Vec3 cutOrigin{0, 0, 0}, cutNormal{0, 1, 0}, cutSide{1, 0, 0};   // the cut guide's half-plane
    // A drawn row-start path instead of the half-plane (nullptr: the half-plane). It plays the same role: vertex 0 of
    // every ring is where the ring crosses it (boundary rings: where they come closest), so it never cuts the mesh.
    const SurfacePathGuide* path = nullptr;
    // The same path with the field's value at its points (nullptr: not given). When given, an interior ring's crossings
    // are the sign changes of t - level along it (levelCrossings, the rule of buildMarkedRows and of the contours: a value
    // equal to the level counts as above it) and its origin is the contour point nearest the crossing. This holds when
    // rows run exactly through vertices (a regular grid), where intersecting the path and the contour per triangle
    // touches only at end points and loses crossings. Without it, `path` is intersected per triangle (pathCrossing).
    const std::vector<MarkedPathPoint>* pathField = nullptr;
    double stagger = 0.6180339887498949;              // sigma_m = frac(0.5 + m * stagger)
    bool adjacentShaping = false;                     // the adjacent-shaping variant (MarkedRowsRequest), off: the frozen rule
    bool boundaryShaping = false;                     // the boundary-shaping variant (SlotKey)
    RowPhase rowPhase = RowPhase::OnPath;
    BindOffCount bindOffCount = BindOffCount::FromLoopLength;
};

struct RingSample {
    int k = 0;
    double level = 0, length = 0, cutParameter = -1;
    int contours = 0, n = 0;
    int guideCrossings = 0;                           // crossings with the cut guide or the row-start path
    bool boundary = false;
};

struct ShapingPlacement {
    int strip = -1, slot = -1, face = -1;
    StitchFaceKind kind = StitchFaceKind::Unsupported;
    int events = 0, slotsInStrip = 0, index = 0;      // the i-th of e events in the strip
    double sigma = 0, targetFraction = 0;
    int targetSlot = -1, moved = 0;
    bool placed = false;                              // false: Unsupported with `reason`
    bool adjacent = false;                            // placed by the adjacent-shaping variant next to another shaping slot
    std::string reason;                               // why here (and the code when not placed)
};

// A closed contour beside ring k too small for a round of its own (a bump or pit of the surface): knitted over, not a
// contour event (row_rings.h, dropContoursShorterThan and kKnitOverStitches). `at` is a point of it, in the field mesh's coordinates.
struct KnittedOverBump {
    int ring = -1;
    double length = 0;
    Vec3 at;
};

struct GuidedRowsProduct {
    TypedStitchMesh mesh;
    std::vector<RingSample> rings;
    std::vector<Contour> contours;                    // per ring k = 0 .. R: the contour tested (empty: a contour event, or the pole)
    std::vector<ShapingPlacement> shaping;
    int rows = 0;
    double hPrime = 0, extent = 0;
    int contourEvents = 0, limitBinds = 0;
    double rowStartOffPath = 0;                       // RowPhase::WalesUpright: the farthest a ring's vertex 0 was turned from the path [mm]
    std::vector<KnittedOverBump> knittedOver;
    Status status = Status::Unresolved;               // FAIL: no mesh (a contour event); UNRESOLVED: an unsupported slot; PASS
    std::string detail;
    std::vector<std::string> log;
};

// castOn and bindOff are the reference's boundary loops (vertex ids in either walking order). The time field t is per
// reference vertex. False only for unusable input; a contour event returns true with a FAIL product and no mesh.
bool buildGuidedRows(const TriMesh& m, const std::vector<double>& t, const std::vector<int>& castOn, const std::vector<int>& bindOff,
                     const GuidedRowsRequest& r, GuidedRowsProduct& out, std::string* err);

// Method A on a marked component: the same rings, counts, shaping rule and lookup, with the row
// origin and the ends taken from saved markings instead of a meridian half-plane and two boundary loops.
//   start     a pole at the start point: ring 0 is the pole vertex; strip 0 a gathered pole start, n_1 triangles (W, T, W)
//             n_1 = round(L_1 / w) with at least 3, raised until strip 1's increases fit the
//             placement rule in use, or the marking's ring1 setting
//   finish    the bind-off ring is the level t_f of the finish point: R = round(t_f / h), h' = t_f / R, ring k the level k h'
//   origin    vertex 0 of ring k is where the row-start path crosses that level. The path is a polyline on the surface
//             (each segment inside one triangle) with the field's value at its points; along a segment the field is linear,
//             so a level is crossed once per sign change. A ring the path crosses other than once is a guide problem: no
//             row origin, no mesh, the crossings listed
//   shaping   the frozen rule of buildGuidedRows, unchanged; the pole strip has no slots to shape (its count is n_1)
//   path      the shift-path runs from ring 1 to ring R (strips 1 .. R - 1); the pole round has none
struct MarkedRowsRequest {
    double w = 5.5, h = 3.4;
    double stagger = 0.6180339887498949;
    Vec3 pole;
    VertexBinding poleBinding;
    double finishLevel = 0;
    int ring1 = 0;                                   // 0: automatic
    std::vector<MarkedPathPoint> path;
    // Variant "adjacent shaping": a shaping slot that cannot keep one K from the other
    // shaping slots within its search window takes the nearest free slot between the ends instead (still never slot 0
    // or the last slot, never in a boundary strip). Off: the frozen rule exactly.
    bool adjacentShaping = false;
    bool boundaryShaping = false;                    // the boundary-shaping variant (SlotKey)
    RowPhase rowPhase = RowPhase::OnPath;
    BindOffCount bindOffCount = BindOffCount::FromLoopLength;
    // A finish curve instead of a finish level: the bind-off ring is this border loop of m (vertex ids, either walking
    // order), T its smallest field value, R = round(T / h); its row origin is its point closest to the path. Empty: the
    // finish level.
    std::vector<int> bindOffLoop;
};
struct PathCrossing {
    int ring = -1;
    double level = 0;
    int crossings = 0;                               // sign changes of t - level along the path
    std::vector<double> pathParameters;              // segment index + fraction of each crossing
    double contourParameter = -1;                    // the row origin on the contour (point i + fraction), -1 without one
    double locateError = 0;                          // [mm] distance of the crossing from the contour segment it was put on
};
struct MarkedRowsProduct {
    GuidedRowsProduct rows;                          // mesh, rings, shaping, log, status (FAIL without a mesh)
    std::vector<PathCrossing> crossings;             // per ring 1 .. R
    int guideProblems = 0;                           // rings the path does not cross exactly once
    int ring1Auto = 0;                               // what the automatic rule gives for n_1
};
bool buildMarkedRows(const TriMesh& m, const std::vector<double>& t, const MarkedRowsRequest& r, MarkedRowsProduct& out, std::string* err);

// How often the path crosses a level: sign changes of t - level between consecutive path points (a value equal to the level
// counts as above it, as the contour does); params (may be null) receives segment index + fraction of each crossing. The
// rows builder and the quick review (quick_review.h) count with this one rule.
int levelCrossings(const std::vector<MarkedPathPoint>& path, double level, std::vector<double>* params);
// The point of `path` at parameter u (segment index + fraction).
Vec3 pathPointAt(const std::vector<MarkedPathPoint>& path, double u);
// The contour parameter (point j + fraction) of the contour point nearest x, preferring the segment that runs through
// triangle tri (a linear field crosses a triangle in one segment); distance (may be null) receives how far x lies from it.
double contourParameterNear(const Contour& c, const Vec3& x, int tri, double* distance);

// The faces of one strip from the two rings' sizes and the shaping slots: what buildGuidedRows and the hand-built
// fixtures share. ringBelow / ringAbove are vertex ids in yarn order (vertex 0 on the path). shapingSlots[i] is a slot
// index; the strip has INC slots when ringAbove is longer, DEC slots when shorter. Appends faces, kinds, roles,
// boundary roles, strips and slots to `s`; false when the counts and slots disagree.
bool assembleStrip(int strip, const std::vector<int>& ringBelow, const std::vector<int>& ringAbove, const std::vector<int>& shapingSlots,
                   EdgeBoundary bottomRole, EdgeBoundary topRole, TypedStitchMesh& s, std::vector<int>* shapingFaces, std::string* err,
                   bool boundaryShaping = false);

} // namespace knit
