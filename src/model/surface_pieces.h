#pragma once
// The pieces of an imported file: the connected parts of each component (a teddy's arm is one piece of the component
// "white_bear_fabric"), welded so that a piece is one surface across its texture seams. Markings are stored on the
// file's triangles (surface_marking.h) and resolved to a piece's own triangles here; a piece keeps every triangle's
// corner order, so barycentric weights carry over unchanged.
#include <cstdint>
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/surface_marking.h"

namespace knit {

struct SurfacePiece {
    int node = -1;                  // the component (ImportedMesh::nodes index)
    int ordinal = 0;                // its place among the component's pieces (1, 2, ...: by lowest file triangle)
    std::vector<int> fileTriangle;  // per piece triangle: the file triangle it is
    TriMesh mesh;                   // file units and axes; welded vertices, normals computed
};

struct SurfacePieces {
    std::vector<SurfacePiece> pieces;   // ordered by lowest file triangle
    std::vector<int> pieceOfTriangle;   // per file triangle (-1: degenerate after welding)
    std::vector<int> localTriangle;     // per file triangle: its index in its piece (-1: degenerate)
};

// Pieces are connected through shared (welded) vertices within one component. weldTolerance in file units.
SurfacePieces findSurfacePieces(const std::vector<Vec3>& pos, const std::vector<std::array<int, 3>>& tri, const std::vector<int>& nodeOfTriangle,
                                double weldTolerance);

// The weld tolerance every user of the pieces takes (the Mark step and the lab must find the same pieces): 1e-6 of
// the bounding box diagonal, far below any modelled detail and above a texture seam's rounding.
double pieceWeldTolerance(const std::vector<Vec3>& pos);

// meshDigest of the file's mesh as read: the identity a markings file records (MarkingSet::digest).
uint64_t surfaceDigest(const std::vector<Vec3>& pos, const std::vector<std::array<int, 3>>& tri);

// A file point as a point of its piece (false: the triangle is outside the file or degenerate).
bool pieceLocalPoint(const SurfacePieces& p, const SurfacePoint& filePoint, int& piece, SurfacePoint& local);

// The piece a marking lies on: the piece of its first point (-1: no points, or not on any piece).
int pieceOfMarking(const SurfacePieces& p, const SurfaceMarking& m);
// True when the marking lies on `piece` and belongs to its knitting component `knitComponent` (SurfaceMarking::knitComponent).
bool markingInKnitComponent(const SurfacePieces& p, const SurfaceMarking& m, int piece, int knitComponent);

// The identity of what a stitch-mesh build of `piece`'s knitting component `knitComponent` uses (stitchmesh/marked_build.h takes the same
// markings): the set's source digest and preparation, then every marking of that component that is neither offered nor
// rejected (a user marking or an accepted proposal), in set order, with its id, meaning, shape, settings and points bit for
// bit. Offered and rejected proposals, other knitting components and other pieces change nothing (component 0 hashes as
// before knitting components existed).
uint64_t effectiveMarkingsDigest(const MarkingSet& set, const SurfacePieces& pieces, int piece, int knitComponent = 0);

// What a proposal for `piece` was computed for: the file, the piece, the reference preparation, the unit and the gauge.
// An offered proposal whose basis differs from today's is outdated: shown, never acceptable (suggest again instead).
uint64_t proposalBasis(uint64_t fileDigest, int piece, const std::string& preparation, double mmPerUnit, double w, double h);

// What accepting (or previewing) a setup on `piece` sets aside: the markings in use on that piece and in the setup's knitting
// component (user markings, and accepted proposals of other setups) whose construction role (markingRole) one of the
// setup's markings takes. Finish and Boundary are the same role; cuts and attachments are never set aside.
struct SetupConflicts {
    std::vector<int> acceptedProposals;   // become Rejected on acceptance
    std::vector<int> userMarkings;        // removed on acceptance only when the user asked to replace them
};
SetupConflicts setupConflicts(const MarkingSet& s, const SurfacePieces& pieces, int piece, int setup);

// Accepts one setup as a whole: its markings Accepted, conflicting accepted proposals Rejected, conflicting user markings
// removed when replaceUser. False, and nothing changed, when the setup is empty, not offered on `piece`, made for another
// basis (outdated), or conflicts with user markings and replaceUser is false.
bool acceptSetup(MarkingSet& s, const SurfacePieces& pieces, int piece, int setup, uint64_t currentBasis, bool replaceUser);

// The set a preview of `setup` is reviewed and built from (never saved): a copy of s without the setup's conflicts, so the
// setup's start, termination and path stand alone for their roles while everything else on the piece (cuts,
// attachments, other roles) and on other pieces stays. The setup's markings stay Offered; previewIds receives those to
// evaluate as if accepted (MarkedBuildRequest::previewOffered: the builder labels the result a research preview).
MarkingSet previewMarkings(const MarkingSet& s, const SurfacePieces& pieces, int piece, int setup, std::vector<int>& previewIds);

// New offered proposals added to the set: every marking gets a fresh id (never reused, proposalId = id) and its setup a fresh
// number above any in the set (setups keep their grouping). newSetups (may be null) receives the set's number of each
// distinct setup of `offered`, in order of first appearance.
void appendProposals(MarkingSet& s, const std::vector<SurfaceMarking>& offered, std::vector<int>* newSetups = nullptr);
// New offered proposals for `piece` replace the piece's previous offered and rejected proposals of the same knitting
// component (the component the offered ones carry; neither is input nor the user's intent); user markings, accepted
// proposals, other components and other pieces stay; then appendProposals.
void mergeProposals(MarkingSet& s, const SurfacePieces& pieces, int piece, const std::vector<SurfaceMarking>& offered, std::vector<int>* newSetups = nullptr);

} // namespace knit
