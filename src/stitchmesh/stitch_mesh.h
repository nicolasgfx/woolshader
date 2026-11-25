#pragma once
// Stage 4 of the stitch-mesh layout: the typed stitch mesh (Wu et al. 2018, section 7) from the labelled, directed coarse
// mesh of stages 2 and 3.
//
// Every coarse edge gets its midpoint and every coarse face its centre; a coarse quad becomes four quads and a coarse
// triangle three. The new edge from an edge's midpoint to the face centre takes the role perpendicular to that edge, so
// a coarse row becomes two stitch rows. A coarse triangle with two course edges yields a SPECIAL QUAD at its apex (two
// wale edges meeting at the centre, a bottom and a top course edge meeting at the apex): the end of a row, where the
// yarn turns from the lower stitch row into the upper one and the loop below passes through to the row above. A coarse
// triangle with two wale edges yields a quad with two course edges of the same side (both top or both bottom) meeting at
// the centre; it is split along the wale diagonal from the apex and each half merged with the neighbouring stitch face
// across the apex's wale edge, giving two pentagons: INCREASE (one bottom, two top: the paper's kp) or DECREASE (two
// bottom, one top: d12k), consecutive in their stitch row. Nothing else is rewritten; a face that fits no type is kept
// and named UNSUPPORTED, never dropped and never accepted.
//
// Ports (section 2.2): a face's yarn enters through one wale edge and leaves through the other; its bottom course edges
// receive the loops of the row below, its top course edges offer loops to the row above. Faces joined across their wale
// edges form the STRANDS, each a closed yarn or an open one with two ends on the mesh border. A course edge offered by
// one face and received by the other is an interlock; one offered or received by both is a direction mismatch, listed.
// An independent validator (validateStitchMesh) re-derives all of this from the faces and roles alone and compares. It
// keeps two questions apart: is the SURFACE valid (manifold, one fan of faces around every vertex, every face wound with
// its coarse face and of positive area, no two faces intersecting, the topology kept: FAIL when not) and is the YARN
// SUPPORTED (every face of a known kind, every course edge an interlock or a declared end, no two faces sharing two edges:
// UNRESOLVED when not, listed). A degree-two interior vertex, two faces sharing two edges, is a valid surface (an edge
// split makes one); whether yarn can run through it is a question for the yarn stage, so it is listed, not refused.
#include <string>
#include <utility>
#include <vector>

#include "model/status.h"
#include "stitchmesh/polygon_mesh.h"
#include "stitchmesh/stage_report.h"

namespace knit {

// PoleStart (stockinette experiment): a triangle (W+, T, W-) of a gathered pole start; only
// the typed mesh of the stockinette experiment makes it. Last, so the other values (and stored digests) are unchanged.
enum class StitchFaceKind : char { Regular, Increase, Decrease, SpecialQuad, Unsupported, PoleStart };
enum class StitchEdgeRole : char { Wale, Bottom, Top };

const char* toString(StitchFaceKind k);
const char* toString(StitchEdgeRole r);

// One continuous yarn: the faces it passes in order, each entered through one wale edge and left through the other.
struct Strand {
    std::vector<int> faces;
    bool closed = false;   // the last face's exit edge is the first face's entry edge
};

// One course edge of the stitch mesh seen from its two faces.
struct LoopPort {
    enum class Kind : char { Interlock, BothTop, BothBottom, BoundaryBottom, BoundaryTop };
    int a = -1, b = -1;          // the edge (min, max vertex)
    int offering = -1;           // the face whose top edge it is (-1: none)
    int receiving = -1;          // the face whose bottom edge it is (-1: none)
    int faceA = -1, faceB = -1;  // its two faces, for the mismatch kinds (BothTop: both offer; BothBottom: both receive)
    Kind kind = Kind::Interlock;
};

struct StitchMesh {
    PolygonMesh mesh;                              // vertices: the coarse ones, then the edge midpoints, then the centres
    std::vector<StitchFaceKind> kind;              // per face
    std::vector<std::vector<StitchEdgeRole>> role; // per face, per corner edge k (the edge from corner k to corner k + 1)
    std::vector<int> coarseFace;                   // per face: the coarse face it came from (a pentagon: the quad's)
    // yarn
    std::vector<Strand> strands;
    std::vector<int> strandOfFace;                 // -1: an unsupported face (its wale edges end the strands they touch)
    std::vector<char> reversed;                    // per face: the yarn runs against the face's own course direction
    int yarnEnds = 0;                              // wale edges on the mesh border (each ends an open strand)
    int unsupportedEnds = 0;                       // wale edges at an unsupported face
    // loops
    std::vector<LoopPort> loops;
};

// What the validator finds, from the faces and roles alone, compared with the builder's strands and loops.
struct StitchMeshCheck {
    Status status = Status::Unresolved;  // FAIL on a surface fault (below) or a record that contradicts the faces and roles;
                                         // UNRESOLVED when the surface is valid but the yarn support is open (unsupported
                                         // faces, mismatches, shared edge pairs); PASS otherwise
    Status surface = Status::Unresolved; // the surface alone
    PolygonMeshCheck mesh;
    int brokenVertexLinks = 0;           // a vertex whose faces form more than one fan (a pinched surface)
    int invertedFaces = 0;               // wound against its coarse face's normal, or of no area
    std::vector<int> invertedFaceIds;    // the first of them (up to 64), for inspection
    int nonConvexFaces = 0;              // a corner turning the other way (reported; a pentagon may be)
    Status intersections = Status::Unresolved;   // model/mesh_check.h on the faces' validated ear-clipped triangles; UNRESOLVED
                                                 // (not certified) when a face has no valid triangulation in its plane
    int unvalidatedTriangulations = 0;
    std::vector<int> unvalidatedFaces;             // the first of them (up to 64), for inspection
    int intersectingPairs = 0, unresolvedPairs = 0;
    std::vector<std::pair<int, int>> intersectingFaces;   // the first pairs of faces that intersect (up to 64), for inspection
    int faces = 0, regular = 0, increase = 0, decrease = 0, special = 0, unsupported = 0;
    int facesWithoutTwoWaleEdges = 0, kindMismatches = 0;   // kinds the roles do not support
    int sharedEdgePairs = 0;             // two faces sharing two edges (a degree-two interior vertex between them): the
                                         // shared wale connection of a row that folds back; yarn support to be established
    int strands = 0, closedStrands = 0, openStrands = 0, yarnEnds = 0;
    int waleEdgesTracedOnce = 0, waleEdgesNotTracedOnce = 0;
    int interlocks = 0, bothTop = 0, bothBottom = 0, boundaryBottom = 0, boundaryTop = 0;
    int loopsNotClassified = 0;          // course edges the builder's list misses or classifies differently
    long eulerBefore = 0, eulerAfter = 0;
    int boundaryLoopsBefore = 0, boundaryLoopsAfter = 0;
    std::string detail;
    std::string text() const;
};
StitchMeshCheck validateStitchMesh(const StitchMesh& s, const PolygonMesh& coarse);

struct StitchMeshProduct {
    StitchMesh stitch;
    StitchMeshCheck check;
    int coarseQuads = 0, coarseTriangles = 0, rowEndTriangles = 0, apexTriangles = 0;
    int merges = 0, unmergedApexFaces = 0;
    StageReport report;   // PASS iff check PASS with no unsupported face and no mismatch; UNRESOLVED when those are listed;
                          // FAIL when the check fails
};

// top: per coarse face, per corner edge k, 1 top, 0 bottom, -1 wale (DirectedMesh::top). False only when the input cannot
// be read as a labelled coarse mesh (sizes disagree, a face of fewer than 3 or more than 4 sides, a non-manifold edge).
bool buildStitchMesh(const PolygonMesh& coarse, const std::vector<std::vector<char>>& top, StitchMeshProduct& out, std::string* err);

} // namespace knit
