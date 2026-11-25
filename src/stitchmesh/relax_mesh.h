#pragma once
// Checkpoint 3: mesh-based relaxation of a typed stitch mesh on its reference surface (Yuksel et al. 2012, section 6.1),
// before any yarn.
//
// The model of the paper:
//   dynamics      a first-order model x' = f(x), stepped by backward Euler, (I - h df/dx) dx = h f0 (Eq. 1), solved by the
//                 filtered Conjugate Gradients of Baraff and Witkin 1998; the normal component of dx is discarded and each
//                 vertex then moves ON the surface in the direction of dx. Vertices on the stitch-mesh border stay on the
//                 border of the surface.
//   stretch       every course and wale edge, f_i = k_stretch (l / |xi - xj| - 1) (xi - xj) / |xi - xj| (Eq. 2), with
//                 l = r_course or r_wale; also the two diagonals of every sub-quad, l = sqrt(r_course^2 + r_wale^2)
//   shear         every corner j of every sub-quad with neighbours i, k: f_j = -k_shear ((xi - xj) . (xk - xj)) (xj - (xi + xk) / 2),
//                 half of it the other way on i and on k (Eq. 3)
//   wale strut    consecutive wale edges i-j, j-k: f_i = -k_wale (|xi - xk| / r_ijk - 1) (xi - xk) / |xi - xk|, the opposite on k,
//                 r_ijk = max(r_wale, |xi - xj| + |xk - xj|) (Eq. 4)
//   sub-quads     a face's bottom course edges paired with its top course edges (paper Fig. 15): a quad has one, a kp or
//                 d12k pentagon two, a face whose bottom and top course edges share a vertex (a triangle, the helix's start
//                 step, a special quad) none
//   stiffness     k_stretch = 2, k_shear = 0.2, k_wale = 2 (the paper's values)
// Departures:
//   S1  rest lengths r_course = w, r_wale = h from the frozen gauge, not alpha A / N
//   S2  the surface is the reference triangle mesh (a directional walk across its triangles), not a subdivision surface
//   S3  lengths in units of r_course, so the paper's dimensionless stiffnesses apply at any scale (the paper gives no units)
//   S4  the implicit matrix holds the stretch and strut Jacobians without their transverse terms when those would make it
//       indefinite (the standard positive-definite approximation); shear is applied explicitly. Equilibrium (f = 0) is the
//       paper's; only the path to it differs.
//   S5  a per-step displacement limit (maxStep) as a safeguard, and a stop rule: converged, stalled, step cap or time cap.
#include <array>
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/status.h"
#include "stitchmesh/reference_surface.h"
#include "stitchmesh/stitch_mesh.h"

namespace knit {

// One face as the solver sees it: corners counter-clockwise from outside and the role of each corner edge. For Method A
// these are the virtual polygons of the shift (a start step repeats one corner: its collapsed edge carries no spring).
struct RelaxFace {
    std::vector<int> corners;
    std::vector<StitchEdgeRole> roles;
};

// What keeps the border vertices of one reference loop in their order around it.
enum class BorderOrder {
    // Nothing (CP3): each border vertex moves along the loop by its own step. A step can carry it past its neighbour; the
    // course spring between them is then at rest again on the other side. On the teddy body two vertices of the bind-off
    // ring passed the held first vertex this way.
    Free,
    // A step along the loop is at most 0.45 of the gap to the neighbour on that side, measured at the start of the step.
    // Two neighbours moving towards each other then keep a tenth of their gap, so the order the layout gave never changes.
    Kept,
    // All border vertices of a loop move by the same arc length per step, the mean of their own steps (none when one of
    // them is held): the ring keeps its spacing and turns around the loop with the rows below it. Measured reason: with
    // Kept, the wale springs of a bind-off strip laid out taller than h pull each border vertex to the point of the loop
    // nearest the ring below; on the teddy body eight of 22 gathered in two bends of the loop, 0.05 mm apart, and the
    // forces grew a hundredfold. Where the vertices sit on the loop for the fit is the border placement's part
    // (border_placement.h). Not the default: the ring cannot settle its own spacing, so a run that would converge with
    // Kept ends as stalled (the open tube of the tests: RMS force 0.074 after 61 steps against 1e-4 after 26).
    RingSlidesAsOne,
};

// What the solver does with a step that leaves larger forces than it found.
enum class StepControl {
    // Nothing (CP3): the time step is the request's throughout. The shear force is applied explicitly (S4); in rows laid
    // out strongly sheared (the bind-off strips of the teddy body, wale edges 7 to 9 mm for h = 3.4 mm) the RMS force
    // then grew from 0.29 to 13.9 over forty steps, and the run ended as stalled after its first fifty.
    FixedTimeStep,
    // The step is taken back and the time step halved (down to a 64th of the request's); after a step that lowers the
    // forces the time step grows by a quarter again.
    ForcesMayNotGrow,
};

// What keeps a face from turning over in one step.
enum class FaceCorners {
    // Nothing (CP3): every vertex takes its solved step (at most maxStep). Two ends of a short edge can pass each other
    // in one step; the springs are then at rest again in the folded state. On the teddy body a valid quad beside a
    // decrease in the bind-off strip was turned over this way.
    Free,
    // Per corner (a, b, c) of every solver face, the size (b - a) x (c - b) . n, with n the reference normal under b at
    // the start of the step, stays positive: a step that would take a corner of at least w h / 4 below w h / 20 is
    // halved for the corner's three vertices, up to six times, then left out. A smaller corner (a straight pentagon
    // corner, the thin corners at a pole, a face the layout laid out folded) is not judged: its sign changes with a
    // small bend of the face's outline, not with the face turning over.
    Kept,
};

struct MeshRelaxRequest {
    double w = 5.5, h = 5.5;                  // r_course, r_wale [mm]
    double kStretch = 2, kShear = 0.2, kWale = 2;
    double timestep = 1.0;                    // h of Eq. 1, in the normalized units (S3)
    double maxStep = 0.25;                    // largest vertex move per step / w (S5)
    int maxSteps = 2000;
    double capSeconds = 120;
    double forceTolerance = 1e-4;             // converged: RMS of the filtered force per free coordinate below this
    int stallWindow = 50;                     // stalled: RMS force not below stallRatio x its value stallWindow steps earlier
    double stallRatio = 0.99;
    int cgMaxIterations = 300;
    double cgTolerance = 1e-8;                // relative residual
    // Filter normal (solver note R1): false, the bound triangle's own normal (CP3); true, the vertex normals
    // interpolated across it, a tangent plane that does not jump where a vertex crosses a crease of a faceted reference. The
    // walk stays on the facets either way, so the surface is the same; only the projection of forces and steps changes.
    bool smoothNormals = false;
    // Solver note R2, kink constraints (integration report): the faceted surface and the guide polylines are not smooth at
    // their edges and corners, and a vertex resting against such a kink is an equilibrium of the constrained problem
    // (KKT) that a one-sided projection can only circle. Per step, a free vertex within kinkTolerance x w of a facet edge
    // whose force pushes into both facets is held on that edge and moves only along it; a vertex on a guide line or border
    // loop within that distance of a corner whose force pushes into the corner from both sides is held there. False:
    // the CP3 behaviour exactly.
    bool kinkConstraints = false;
    double kinkTolerance = 0.01;
    BorderOrder borderOrder = BorderOrder::Kept;
    FaceCorners faceCorners = FaceCorners::Kept;
    StepControl stepControl = StepControl::ForcesMayNotGrow;
};

struct RelaxStep {
    int step = 0;
    double rmsForce = 0, maxForce = 0, maxMove = 0;
    int maxForceVertex = -1;                  // where the largest filtered force acts
    int maxForceTri = -1;                     // the reference triangle that vertex is bound to at this step
    double maxForceEdgeDistance = 0;          // its distance from the nearest edge of that triangle [mm]
    int cgIterations = 0;
};

struct RelaxProduct {
    std::vector<Vec3> pos;
    std::vector<VertexBinding> binding;
    std::vector<RelaxStep> history;
    std::string stopReason;                   // converged | stalled | step cap | time cap | diverged
    Status status = Status::Unresolved;       // PASS converged; UNRESOLVED stalled or capped; FAIL diverged
    int steps = 0;
    double ms = 0;
    double rmsForceStart = 0, rmsForceEnd = 0;
    // the elements built
    int courseSprings = 0, waleSprings = 0, diagonalSprings = 0, shearCorners = 0, struts = 0;
    int addedCourseSprings = 0;               // of courseSprings: RelaxGuides::courseEdges that no face had
    int subQuads = 0, facesWithoutSubQuads = 0, waleJunctionsWithoutStrut = 0, conflictingEdgeRoles = 0;
    int borderVertices = 0, borderVerticesSnapped = 0, walkClamps = 0;
    int borderMovesLimited = 0;               // BorderOrder::Kept: steps along a loop cut short at 0.45 of the gap
    int cornerMovesLimited = 0;               // FaceCorners::Kept: vertex steps shortened, summed over the steps
    int stepsTakenBack = 0;                   // StepControl::ForcesMayNotGrow
    std::string firstCornerLimit;             // the first of them: step, corner, its size before and after the full step
    int borderOrderViolations = 0;            // BorderOrder::Free: times a border vertex passed its neighbour on the loop
    std::string firstBorderOrderViolation;    // the first of them: step, the two vertices, their gap before and after
    int pinnedVertices = 0, guidedVertices = 0;
    int kinkHeldLast = 0, kinkOnEdgeLast = 0;  // R2: vertices held at a corner / on a crease edge in the last step
    double maxGuideSnap = 0;                  // [mm] guided vertices moved onto their line before the first step
    double maxBorderSnap = 0;                 // [mm] border vertices moved onto the reference border before the first step
    std::string detail;
};

// Constraints from a component's markings: what the markings fix, the solver keeps. They restrict where
// vertices may go; the forces (the model) are unchanged.
//   pinned   a vertex that stays where it starts (the pole of a pole start, at the start marking)
//   lines    polylines on the reference, each segment inside one reference triangle; a vertex with lineOf >= 0 slides along
//            its line as a border vertex slides along its loop (the bind-off ring on its level contour, a row's first
//            vertex on the row-start path). An open line keeps its vertices between its ends.
struct GuideLine {
    std::vector<Vec3> points;
    std::vector<int> segmentTri;   // per segment i (points i, i + 1); a closed line has one segment more, back to point 0
    bool closed = false;
};
struct RelaxGuides {
    std::vector<GuideLine> lines;
    std::vector<int> lineOf;       // per stitch vertex, -1: none (then `border` applies)
    std::vector<char> pinned;      // per stitch vertex
    // Course edges that no solver face has, each with the stretch spring of a course edge (Eq. 2). Method A's solver
    // faces are the virtual (helical) ones: slot 0 of the last strip is read one ring lower, so the edge between the
    // first and the second vertex of the bind-off ring belongs to no face, and without its spring nothing holds the two
    // apart.
    std::vector<std::array<int, 2>> courseEdges;
};
// reference: the surface; refLoops: its boundary loops (vertex ids in walking order); border: per stitch vertex, 1 when it
// lies on the stitch mesh's border (it then slides along the nearest reference loop). start: the stitch vertices.
// guides (optional): the markings' constraints above; a vertex on a guide line or pinned ignores `border`.
bool relaxStitchMesh(const TriMesh& reference, const std::vector<std::vector<int>>& refLoops, const std::vector<RelaxFace>& faces,
                     const std::vector<Vec3>& start, const std::vector<char>& border, const MeshRelaxRequest& r, RelaxProduct& out, std::string* err,
                     const RelaxGuides* guides = nullptr);

// The faces of a polygon mesh with its roles (Method B, the fixtures without the shift).
std::vector<RelaxFace> relaxFacesOf(const PolygonMesh& m, const std::vector<std::vector<StitchEdgeRole>>& role);

// Stationarity under the actual faceted constraint (integration report, solver note R1): for each vertex, the size of the
// force component along the motions the reference really allows it, in the solver's units.
//   inside a facet      the force projected on the facet
//   on an edge (within edgeTolerance [mm] of one)   the force along the edge, plus any part pulling into either facet
//                       (a vertex held on a crease by forces pushing into both facets is at rest there)
//   on a mesh vertex    the largest pull into any incident facet (directions sampled)
//   border / guided     the force along the loop or guide line at the vertex; at a corner of the polyline (within
//                       edgeTolerance), only the parts leading out along either segment count (1-D, as a crease); pinned: 0
// A state is a faceted equilibrium when every value is ~0, whatever filter produced it.
struct WalkResidual {
    std::vector<double> perVertex;
    double rms = 0, max = 0;
    int maxVertex = -1, interior = 0, onEdge = 0, onVertex = 0, constrained = 0;
};
WalkResidual walkResidualOf(const TriMesh& reference, const std::vector<Vec3>& pos, const std::vector<Vec3>& forces, const std::vector<char>& border,
                            const std::vector<std::vector<int>>& refLoops, const RelaxGuides* guides, double edgeTolerance);
// The same with the model's forces at pos (relaxForces).
WalkResidual walkStationarity(const TriMesh& reference, const std::vector<RelaxFace>& faces, const std::vector<Vec3>& pos, const std::vector<char>& border,
                              const std::vector<std::vector<int>>& refLoops, const MeshRelaxRequest& r, const RelaxGuides* guides, double edgeTolerance);
// The total force of the model at x (normalized units), for tests: per vertex, unfiltered.
std::vector<Vec3> relaxForces(const std::vector<RelaxFace>& faces, const std::vector<Vec3>& x, const MeshRelaxRequest& r);

} // namespace knit
