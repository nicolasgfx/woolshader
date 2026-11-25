#pragma once
// Stage 1 of the stitch-mesh layout: the reference surface remeshed into
// a field-aligned, quad-dominant coarse mesh at twice the stitch width (Wu et al. 2018, section 4: a 2-RoSy orientation
// field, a position field, the extraction; one coarse quad becomes 2 x 2 stitches in section 7). The backend is RIM
// (third_party/rim), behind this stage's own checks:
//   validity      StructuralCheck: manifold, oriented, no inverted, collapsed or unbound element. What any later stage
//                 needs to read the mesh at all.
//   qualification GeometricQualification: edge lengths, angles, surface distance, field alignment, irregular vertices and
//                 triangles within the bounds written next to each field. What a stitch mesh needs from its coarse mesh.
// Up to three attempts (intrinsic field, more sweeps, extrinsic field) run only while the previous one is not qualified;
// the kept one ranks qualified before valid before invalid, then by distortion, surface fidelity, irregular vertices,
// then the earlier attempt. The product's status is PASS only for a qualified mesh; every attempt stays inspectable.
#include <string>
#include <vector>

#include "stitchmesh/quantiles.h"
#include "model/status.h"
#include "model/vec.h"
#include "stitchmesh/polygon_mesh.h"
#include "stitchmesh/reference_surface.h"
#include "stitchmesh/stage_report.h"

namespace knit {

// Refused before anything is allocated for the solve.
struct SizeLimits {
    int inputTriangles = 1'000'000;
    int finestVertices = 400'000;   // the input after RIM's subdivision (knit_rim::predictedVertices)
    int coarseFaces = 200'000;      // area / coarse edge^2
};

struct RemeshRequest {
    double w = 5.5, h = 5.5;        // the stitch gauge [mm]: coarse edges of coarseFactor x w along the field (the course) and
                                    // coarseFactor x h across it (the wale); h != w runs RIM's rectangular lattice (intrinsic 2-RoSy)
    double coarseFactor = 2;        // coarse edge = coarseFactor x w
    bool swapAxes = false;          // the course across the field instead of along it (after the roles came out exchanged)
    uint64_t seed = 0;              // 0: RIM's own default start
    int sweepsPerLevel = 20;
    bool extrinsic = false;         // the first attempt's field; the third attempt uses the other kind
    int rosy = 2;                   // the field's symmetry: 2 (Wu 2018 section 4) or 4 (RIM's own, a lab experiment)
    std::vector<Vec3> axisGuide;    // per reference vertex: a fixed 2-RoSy axis (tests, M5 guides); empty: solved
    int maxAttempts = 3;            // the initial attempt plus at most two retries
    SizeLimits limits;
    double attemptCapS = 120;       // one attempt's deadline [s]
    double coarseEdge() const { return coarseFactor * 0.5 * (w + h); }   // the thresholds' edge
    double coarseEdgeAlong() const { return coarseFactor * (swapAxes ? h : w); }    // along the field's axis
    double coarseEdgeAcross() const { return coarseFactor * (swapAxes ? w : h); }   // across it
};

// A global axis guide by name: "ring-x", "ring-y", "ring-z" put the course across that world axis (the rows run round
// the object), "along-x" .. "along-z" along it, "boundary" the rows along the open borders (boundaryAxisGuide); any
// other name ("", "none") gives an empty guide (the field is solved). Per vertex of `mesh` (its vertex normals), unit,
// in the tangent plane.
std::vector<Vec3> globalAxisGuide(const std::string& name, const TriMesh& mesh);

// The rows along the open borders of a sheet (M6, open boundaries): every boundary loop's dominant direction (its edges'
// 2-RoSy mean: the long sides of a rectangle) or, for a loop with no dominant direction (a round cuff), the ring round
// the loop's axis; every vertex takes its nearest loop's direction in its tangent plane. Empty for a closed mesh.
std::vector<Vec3> boundaryAxisGuide(const TriMesh& mesh);

struct SizePreflight {
    int inputTriangles = 0, predictedFinestVertices = 0, predictedCoarseFaces = 0;
    bool ok = false;
    std::string why;
};
SizePreflight preflightSize(const ReferenceSurface& s, const RemeshRequest& r);

// The 2-RoSy line field at the coarse vertices: axis and -axis are one direction.
struct LineField {
    std::vector<Vec3> axis;         // per coarse vertex, unit, in the tangent plane
    std::vector<Vec3> orientationSingularPoints, positionSingularPoints;  // on the input surface
};

struct StructuralCheck {
    PolygonMeshCheck mesh;
    int invertedFaces = 0;          // face normal against the reference normal under the face's centroid
    int collapsedFaces = 0;         // area below 1e-3 coarse edge^2
    int unboundVertices = 0;        // farther than half a coarse edge from the reference surface
    Status status = Status::Unresolved;  // PASS iff mesh PASS and the three counts are zero
    std::string detail;
    std::string text() const;
};

// Repository bounds, written next to each field; a bound that rejects a configuration the paper shows as regular is
// revised with a recorded reason, never switched off.
struct GeometricQualification {
    // Edges between two interior vertices, / coarse edge: p5 >= 0.85 and p95 <= 1.15. Edges touching an open boundary are
    // reported apart (boundaryEdgeLength) and not bounded: the remesher's lattice has no phase at a boundary, so the cells
    // next to it are stretched or shrunk by up to half a cell (the proposal: "apart from declared boundaries").
    Quantiles edgeLength;
    Quantiles boundaryEdgeLength;
    Quantiles minAngleDeg;          // per face                                  p5 >= 45
    Quantiles surfaceDistance;      // binding distance / coarse edge            p95 <= 0.10
    Quantiles fieldAlignmentDeg;    // interior edge against the axis, mod 90    p95 <= 20
    std::vector<int> valenceHistogram;  // interior vertices, index = valence
    int interiorVertices = 0, irregularVertices = 0;  // interior valence != 4   <= 10 % of the interior vertices
    double triangleFraction = 0;    //                                           <= 0.20
    // The local defects the quantiles hide, kept visible for subdivision and relaxation: edges outside [0.5, 2.0] coarse
    // edges and faces with a corner under 20 degrees. Reported, not bounded by the qualification (which reads the quantiles).
    struct EdgeDefect {
        int a = -1, b = -1;
        double length = 0;  // / coarse edge
    };
    struct FaceDefect {
        int face = -1;
        double minAngleDeg = 0;
    };
    std::vector<EdgeDefect> extremeEdges;
    std::vector<FaceDefect> sharpFaces;
    double distortionScore = 0;     // |p95 edge - 1| + |p5 edge - 1| + (90 - p5 minAngle) / 90
    double fidelityScore = 0;       // p95 surface distance
    Status status = Status::Unresolved;  // PASS iff every bound holds; FAIL names the first that does not
    std::string detail;
    std::string text() const;
};

struct RemeshAttempt {
    std::string variant;            // "intrinsic-20", "intrinsic-60", "extrinsic-20"
    PolygonMesh mesh;
    std::vector<Vec3> normal;       // per coarse vertex
    LineField field;
    std::vector<VertexBinding> binding;
    NormalizeLog normalization;
    StructuralCheck structural;
    GeometricQualification quality;
    StageReport report;
    std::string backendLog;
};

struct RemeshProduct {
    RemeshAttempt kept;             // the ranked first attempt (an empty mesh when none ran)
    std::vector<RemeshAttempt> rejected;
    SizePreflight preflight;
    StageReport report;             // PASS: kept is qualified; UNRESOLVED: valid but unqualified, or nothing to solve; FAIL: no valid attempt
};

// False only for an unusable request (no reference mesh); a FAIL or UNRESOLVED product returns true with its report.
bool remeshStage(const ReferenceSurface& s, const RemeshRequest& r, RemeshProduct& out, std::string* err);

// The attempts' order (qualified, then valid, then invalid; distortion, fidelity, irregular vertices, the earlier one): the
// index of the one to keep. For the lab's supervised attempts, which come back from child processes.
int bestAttempt(const std::vector<RemeshAttempt>& attempts);

// The checks alone, for fixtures and for meshes from elsewhere. `coarseEdge` scales the bounds.
StructuralCheck checkStructure(const PolygonMesh& m, const ReferenceSurface& s, const std::vector<VertexBinding>& binding, double coarseEdge);
// coarseEdgeAlong / coarseEdgeAcross: an interior edge is measured against the target of its own direction (within 45
// degrees of the axis at its origin: along; else across).
GeometricQualification qualifyGeometry(const PolygonMesh& m, const std::vector<Vec3>& axis, const std::vector<VertexBinding>& binding,
                                       double coarseEdgeAlong, double coarseEdgeAcross);

} // namespace knit
