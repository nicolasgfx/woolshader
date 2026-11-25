#pragma once
// Stage 2 of the stitch-mesh layout: one role per edge, course or wale (Wu et al. 2018, section 5), on the coarse mesh
// of stage 1. A regular quad alternates course and wale around its boundary; a triangle has one or two course edges.
//
// The roles come from a binary program on the HALF-edges (section 5.1, solved with HiGHS): each face's half-edges take
// a valid configuration, and the number of interior edges whose two half-edges disagree is minimized. The disagreeing
// edges are then settled by mesh operations (section 5.2: a quad split along a diagonal, a triangle pair rotated), and
// the triangles that can merge into valid quads do (section 5.3). An independent validator reads the result: it never
// trusts the solver's own status. Before the program, the quads around interior vertices of odd valence are
// triangulated (section 5.1, Fig. 5), since no valid alternation exists around them.
//
// The vertices never change; faces split, merge or exchange a diagonal. Each operation is counted and the topology
// (vertex set, boundary edges, Euler characteristic) compared before and after, so a hole cannot slip in.
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "model/status.h"
#include "model/vec.h"
#include "stitchmesh/polygon_mesh.h"
#include "stitchmesh/stage_report.h"

namespace knit {

enum class EdgeRole : char { Wale = 0, Course = 1 };

// One role per undirected edge (min vertex, max vertex) of a polygon mesh.
struct EdgeLabels {
    std::map<std::pair<int, int>, EdgeRole> role;
    bool has(int a, int b) const;
    EdgeRole at(int a, int b) const;  // the edge must exist
    int courseEdges() const;
};

// The mean length of the course edges over the mean length of the wale edges: on a rectangular remesh the roles must
// follow the metric (the ratio near w / h), else the coarse axes came out exchanged. 0 without edges of both roles.
double courseToWaleLength(const PolygonMesh& m, const EdgeLabels& labels);

struct LabelOptions {
    double timeLimitS = 120;       // the program; past it the incumbent is used and reported as feasible, not optimal
    bool triangulateOdd = true;    // section 5.1
    bool postProcess = true;       // section 5.3
    bool heuristicStart = true;    // a feasible start for HiGHS from the quads' parity classes
    // Section 5.1 allows a triangle one course edge (a shaping triangle, stage 4's kp or d12k pentagon) or two (a row end,
    // stage 4's special quad) and the paper's program does not distinguish them. A row end at an odd vertex has both
    // course edges against the same neighbouring row, so one of them is a wale-direction mismatch whatever stage 3
    // decides (sphere, 2026-11-13: 42 row ends, 18 mismatches, 49 unlinked loops). With this, the program charges each
    // course half-edge of a triangle one unit below the cost of one disagreeing edge, forbids the two-wale pairs that
    // section 5.3 would flip, and forbids a row end whose neighbours are one row (labels.cpp). Measured 2026-11-14 on
    // the 5.5 mm pipelines (direction mismatches, paper's program -> this): flat 20 -> 6, torus 72 -> 26, sphere 36 -> 30,
    // capsule 24 -> 44; the flat's stitch mesh then fails its trace. An experiment (stitchmesh_lab --prefer-shaping), off.
    bool preferShaping = false;
    // The start prefers course for edges more perpendicular than 45 degrees to `up` (rows climb `up`): a tie-break among
    // equally good labellings, a repository choice, not a term of the paper's program.
    Vec3 up{0, 0, 1};
};

// What the independent validator finds: a face whose roles are no valid configuration, an edge whose faces disagree
// (impossible once roles are per edge; counted from the half-edge solution before 5.2), and two triangles with two
// wale edges each sharing a wale edge (section 7 cannot subdivide those: a gate for stage 4).
struct LabelCheck {
    Status status = Status::Unresolved;   // PASS iff invalidFaces == 0 and every edge has a role
    int faces = 0, invalidFaces = 0, unlabeledEdges = 0, twoWalePairs = 0;
    std::string detail;
};
LabelCheck validateLabels(const PolygonMesh& m, const EdgeLabels& labels);

struct LabelProduct {
    PolygonMesh mesh;              // after 5.1 to 5.3: the same vertices, other faces
    EdgeLabels labels;
    LabelCheck check;
    // the program
    std::string solverStatus;      // "optimal", "time limit", or HiGHS's own word
    double objective = 0, bound = 0, gap = 0, msSolve = 0;
    bool optimal = false;
    // the operations
    int dartsRefused = 0;          // triangle pairs not merged because the quad would have a reflex corner
    int oddVertices = 0, triangulatedQuads = 0, untriangulatedQuads = 0;
    int rowEndTriangles = 0;       // triangles with two course edges in the program's solution (before 5.2 and 5.3)
    int pairRules = 0;             // triangle pairs constrained by section 5.3's two-wale rule inside the program
    int apexRules = 0;             // triangle corners at interior vertices constrained against the forced mismatch
    int halfEdges = 0, interiorEdges = 0, conflicts = 0, conflictQuadQuad = 0, conflictQuadTriangle = 0, conflictTriangleTriangle = 0;
    int splitQuads = 0, keptQuads = 0, unsplittableQuads = 0, rotatedPairs = 0, unrotatablePairs = 0;
    int merged = 0, flipMerged = 0, walePairFlips = 0;
    // topology before and after
    bool sameVertices = false, sameBoundary = false, sameEuler = false;
    long eulerBefore = 0, eulerAfter = 0;
    StageReport report;            // PASS iff check PASS and the topology is kept; UNRESOLVED when feasible but not optimal
};

// False only when the mesh cannot be labelled at all (a face of 5 or more sides, a non-manifold edge, an infeasible
// program); a product with a FAIL verdict returns true with its report.
bool labelEdges(const PolygonMesh& m, const LabelOptions& o, LabelProduct& out, std::string* err);

} // namespace knit
