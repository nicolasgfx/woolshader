#pragma once
// The typed stitch mesh both methods of the stockinette experiment produce,
// and its independent check.
//
// Faces are polygons on the reference surface, corners counter-clockwise seen from outside. Every corner edge carries a
// role (Y12 section 4: a course edge, Bottom or Top, crossed by two yarns; a Wale edge crossed by one) and a boundary role
// (Interior, or CastOn / BindOff on a declared open boundary, or Side). A face's kind follows from its roles alone:
//   K            (B, W, T, W)        Regular
//   INC kp       (B, W, T, T, W)     Increase, one bottom and two top course edges
//   DEC d12k     (B, B, W, T, W)     Decrease, two bottom and one top course edge
//   special quad (B, T, W, W)        SpecialQuad: the course edges meet, the wale edges meet (W18 section 7, Method B only)
//   pole start   (W, T, W)           PoleStart: a triangle of a gathered pole start, its bottom edge collapsed into the
//                                   pole vertex (proposed); its strip has a one-vertex ring below
//   anything else, or a slot the lookup could not fill: Unsupported, with a reason code.
//
// The check re-derives everything from the faces and roles and keeps three verdicts apart:
//   connectivity  the surface is manifold and oriented, kinds agree with the roles, every shared course edge is an
//                 interlock (Top of one face, Bottom of the other), yarn agrees across wale edges, boundary roles sit
//                 exactly on the declared boundaries, the count invariant holds per strip. FAIL for a contradiction,
//                 UNRESOLVED for listed open items (unsupported faces, direction mismatches), PASS otherwise.
//   geometry      inverted, degenerate or intersecting faces, and the vertices' distance from the reference: what a
//                 relaxation cannot be asked to repair if it is wrong before it, reported apart from connectivity.
//   cells         the frozen cell bounds, per face: measured and listed, never a CP2 verdict, since
//                 relaxation (CP3) moves vertices and is expected to change them.
#include <array>
#include <string>
#include <utility>
#include <vector>

#include "model/mesh.h"
#include "model/status.h"
#include "stitchmesh/polygon_mesh.h"
#include "stitchmesh/reference_surface.h"
#include "stitchmesh/stitch_mesh.h"

namespace knit {

enum class EdgeBoundary : char { Interior, CastOn, BindOff, Side };
const char* toString(EdgeBoundary b);

// One wale edge of the shift-path: the edge between ring `strip` and ring `strip` + 1 on the cut guide.
struct ShiftPathEdge {
    int strip = -1;
    int bottom = -1, top = -1;   // vertex ids on ring strip and ring strip + 1
};

struct TypedStitchMesh {
    PolygonMesh mesh;
    std::vector<StitchFaceKind> kind;                 // per face
    std::vector<std::vector<StitchEdgeRole>> role;    // per face, per corner edge k (corner k to corner k + 1)
    std::vector<std::vector<EdgeBoundary>> boundary;  // per face, per corner edge k
    std::vector<std::string> reason;                  // per face: why it is Unsupported ("" otherwise)
    // Method A's row structure (empty for Method B): the strip of every face, its slot in the strip (slot 0 has its W-
    // edge on the shift-path), and every ring's vertices in order, vertex 0 on the cut guide.
    std::vector<int> strip, slot;
    std::vector<std::vector<int>> rings;
    std::vector<ShiftPathEdge> shiftPath;
    std::vector<VertexBinding> binding;               // per vertex (empty: not bound)
};

// The kind a face's roles support (cyclic corner order), independent of what was stored.
StitchFaceKind kindFromRoles(const std::vector<StitchEdgeRole>& roles);
// The corner edge through which the yarn leaves a face with a bottom edge: the wale edge right after the (last) bottom
// edge, counter-clockwise (W+). -1 for a special quad or an unsupported face.
int exitWaleCorner(const std::vector<StitchEdgeRole>& roles);
int entryWaleCorner(const std::vector<StitchEdgeRole>& roles);   // W-: the wale edge right before the (first) bottom edge

// One reference boundary the component declares, as a closed polyline, with the role its course edges take.
struct DeclaredBoundary {
    std::string id;
    EdgeBoundary role = EdgeBoundary::CastOn;
    std::vector<Vec3> loop;
};

struct TypedCheckOptions {
    double w = 5.5, h = 5.5;          // gauge [mm]
    double boundaryTolerance = 0.05;  // a boundary-role edge's vertices within this of their declared loop [mm]
    double bindingTolerance = 0.0055; // vertex to reference: 0.001 w
    double faceDeviation = 0.55;      // face centroid and edge midpoints to reference: 0.1 w
    double degenerateArea = 0.1;      // x w h
    bool requireStrips = false;       // Method A: every face has a strip and slot, and the rings are given
};

struct StrandTrace {
    std::vector<std::vector<int>> strands;   // faces in yarn order
    std::vector<char> closed;
    std::vector<int> strandOfFace;
    int yarnEnds = 0;                         // wale edges with no face across
    int waleEdgesTracedOnce = 0, waleEdgesNotTracedOnce = 0;
};
// Faces joined across their wale edges (the yarn passes through one wale edge and leaves through the other; a
// special quad passes it from one wale edge to the other as well, W18 section 8). Re-derived from the polygons.
StrandTrace traceYarnStrands(const PolygonMesh& m, const std::vector<std::vector<StitchEdgeRole>>& role);

struct StripCount {
    int strip = -1;
    int ringBelow = 0, ringAbove = 0;         // vertices of the rings (record)
    int bottomEdges = 0, topEdges = 0;        // course edges of the strip's faces
    int increases = 0, decreases = 0, unsupported = 0, poleStarts = 0;
    bool invariant = false;                   // bottomEdges == ringBelow, topEdges == ringAbove, top - bottom == INC - DEC; a pole strip
                                              // (ring below of one vertex): no bottom edges, every face a pole start, topEdges == ringAbove
};

// One face's geometry check, as checkTypedMesh measured it (the counts in TypedMeshCheck are the sums of these).
struct FaceCheck {
    float deviation = 0;           // face centroid and edge midpoints to the reference [mm]
    bool overDeviation = false;    // deviation above TypedCheckOptions::faceDeviation
    bool inverted = false, degenerate = false, untriangulable = false;
};

struct TypedMeshCheck {
    // connectivity
    Status connectivity = Status::Unresolved;
    PolygonMeshCheck polygon;
    int faces = 0, regular = 0, increase = 0, decrease = 0, special = 0, unsupported = 0, poleStart = 0;
    int kindMismatches = 0, facesWithoutTwoWaleEdges = 0;
    int interlocks = 0, bothTop = 0, bothBottom = 0, courseAgainstWale = 0;
    int waleAgree = 0, waleDisagree = 0;      // shared wale edges: W+ of one face and W- of the other, or not
    int castOnEdges = 0, bindOffEdges = 0, sideEdges = 0;
    int boundaryRoleOnInteriorEdge = 0;       // a boundary role on an edge with a face across
    int borderCourseWithoutRole = 0;          // a course edge on the mesh border with no boundary role
    int boundaryRoleWrongSide = 0;            // CastOn on a Top edge, BindOff on a Bottom edge
    int boundaryRoleOffDeclared = 0;          // a boundary-role edge away from every declared loop of that role
    int borderWaleEdges = 0;                  // yarn ends on an undeclared side
    StrandTrace strands;
    int closedStrands = 0, openStrands = 0;
    std::vector<StripCount> strips;
    int stripsBroken = 0, mixedShapingStrips = 0;
    std::vector<std::pair<std::string, int>> unsupportedReasons;
    // geometry
    Status geometry = Status::Unresolved;
    int invertedFaces = 0, degenerateFaces = 0, intersectingPairs = 0, untriangulableFaces = 0;
    Status intersections = Status::Unresolved;
    double maxBindingDistance = 0;
    double maxFaceDeviation = 0;              // face centroid and edge midpoints to the reference [mm]
    int facesOverDeviation = 0;
    int unboundVertices = 0;
    std::vector<FaceCheck> perFace;           // per face, in face order
    std::string detail;
    std::string text() const;
};

TypedMeshCheck checkTypedMesh(const TypedStitchMesh& s, const TriMesh& reference, const std::vector<DeclaredBoundary>& declared,
                              const TypedCheckOptions& o);

// One face's dimensions against the frozen bounds.
struct CellMeasure {
    int face = -1;
    StitchFaceKind kind = StitchFaceKind::Unsupported;
    double widthRatio = 0, heightRatio = 0;   // K: mean course / w, mean wale / h; pentagons: worst course / w, worst wale / slant
    double minCourse = 0, maxCourse = 0, minWale = 0, maxWale = 0;   // pentagons and special quads: edge ratios
    double minAngle = 0, maxAngle = 0;        // wale-course corners [deg]
    double courseCourseAngle = 0;             // pentagons: the corner between the two course edges [deg]
    double area = 0;                          // / (w h)
    bool withinBounds = false;
    std::string outside;                      // which bound fails, "" when within
};
std::vector<CellMeasure> measureCells(const TypedStitchMesh& s, double w, double h);

// The shaping-cell assessment, reported beside
// the frozen pentagon bounds above, which stay as they are. The relaxation's model (Yuksel et al. 2012 section 6.1, Fig. 15)
// gives a kp or d12k pentagon no trapezoid rest shape: it pairs the single course edge of one side with each course edge of
// the other into two sub-quads, each pulled towards a w x h rectangle, and its wale edges rest at h. So a shaping cell is
// judged by what the model and a later yarn embedding use:
//   interfaces  every course edge / w and every wale edge / h in [0.75, 1.30] (the frozen bounds of a K cell's edges)
//   sub-quads   each of the two sub-quads (bottom course edge, top course edge, as the solver pairs them) non-degenerate
//               (area >= 0.1 w h, the frozen degenerate bound), not folded (both diagonal splits turn like the face), and
//               every corner in [45, 135] degrees (the frozen pentagon wale-course corner range)
// The thresholds are the frozen ones of checkpoint 0, not chosen from any model result.
struct ShapingCellMeasure {
    int face = -1;
    StitchFaceKind kind = StitchFaceKind::Unsupported;
    double minCourse = 0, maxCourse = 0;   // course edge / w
    double minWale = 0, maxWale = 0;       // wale edge / h
    int subQuads = 0;
    double minCorner = 0, maxCorner = 0;   // over the sub-quads' corners [deg]
    double minSubArea = 0;                 // smallest sub-quad area / (w h)
    int foldedSubQuads = 0;
    bool interfacesWithin = false, subQuadsWithin = false, within = false;
    std::string outside;
};
std::vector<ShapingCellMeasure> measureShapingCells(const TypedStitchMesh& s, double w, double h);
struct ShapingSummary {
    int cells = 0, interfacesWithin = 0, subQuadWithin = 0, within = 0, invertedSubQuads = 0;
};
ShapingSummary summarizeShaping(const std::vector<ShapingCellMeasure>& cells);

// A face's triangles by ear clipping in the plane of its Newell normal (a non-convex pentagon fanned from one corner
// would fold over itself). False, with a fan from corner 0 in `out`, when no ear is found (a self-overlapping face).
bool earClipFace(const PolygonMesh& m, int f, std::vector<std::array<int, 3>>& out);

// The typed mesh as one digest: the polygon mesh, then kinds, roles and boundary roles.
uint64_t typedDigest(const TypedStitchMesh& s);

// Method B's product as a typed mesh: kinds and roles as built; border course edges take the role of the declared loop
// they lie on (CastOn for a Bottom edge, BindOff for a Top edge), anything else stays Interior and the check lists it.
TypedStitchMesh typedFromStitchMesh(const StitchMesh& s, const std::vector<DeclaredBoundary>& declared, double tolerance);

} // namespace knit
