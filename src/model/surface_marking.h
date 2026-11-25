#pragma once
// Surface markings: what the user (or a proposal) draws on an imported reference surface to tell the stitch
// mesher where to start, where to finish, where each circular row begins, where to trim, cut or attach.
// A marking stores only its control points, each as a
// triangle of the file's mesh as read and barycentric weights on it, so it survives any later remeshing; the
// path between two control points is derived (surface_path.h), never stored.
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace knit {

// What a marking tells the mesher. A row-start path and a cut stay distinct: only a cut changes topology.
enum class MarkingMeaning {
    Start,         // where the row field begins: a closed curve (cast-on) or a point (pole start, experimental)
    Finish,        // where knitting terminates: a closed curve (the bind-off loop), or a point (v2: the bind-off is the
                   // row field's level through that point; a different request from a curve)
    RowStartPath,  // where each circular row begins and the helix advances: an open path, never a cut
    Boundary,      // where covering stops (e.g. inside the shoulder): a closed curve
    Cut,           // a deliberate change of the surface's construction: an open path or a closed curve
    Attachment     // where two components connect later: a closed curve
};
enum class MarkingShape { Point, Path, ClosedCurve };
enum class MarkingOrigin { User, Proposal };
enum class ProposalState { None, Offered, Accepted, Rejected };  // a proposal is accepted, edited or rejected

// A point on the file's mesh: triangle index in the reader's order and weights on its three corners (sum 1).
struct SurfacePoint {
    int tri = -1;
    double bary[3] = {1, 0, 0};
};

struct SurfaceMarking {
    int id = 0;                                   // unique in its set, never reused
    MarkingMeaning meaning = MarkingMeaning::Start;
    MarkingShape shape = MarkingShape::Point;
    std::vector<SurfacePoint> points;             // control points (a closed curve joins the last to the first)
    MarkingOrigin origin = MarkingOrigin::User;
    ProposalState proposal = ProposalState::None;
    int proposalId = 0;        // a user edit of a proposal keeps the id of the proposal it came from
    int alternativeGroup = 0;  // proposals offered as alternatives to each other share a group (> 0)
    std::string note;          // the user's note, or a proposal's purpose
    std::string reason;        // why it was proposed
    std::string unsuitable;    // non-empty: the proposal is flagged unsuitable, and why
    // Construction settings of this marking (file v2, "setting <key> = <value>"), e.g. a start point's "ring1" (the
    // pole start's first ring: "auto" or a count). Unknown keys are kept and reported by the mesher, never applied.
    std::vector<std::pair<std::string, std::string>> settings;
    // File v4: proposals made together as one compatible setup (a start, a termination and a row-start path from the
    // same field) share a setup number > 0, and record the basis they were computed for (proposalBasis in
    // surface_pieces.h: file, piece, preparation, unit, gauge). A user's edit keeps both, so it can be compared with the
    // proposal; the mesher never reads them, and neither enters effectiveMarkingsDigest.
    int setup = 0;
    uint64_t basis = 0;
    // File v5: the knitting component of its piece the marking belongs to (0: the piece's first). Each component of a piece
    // is knitted on its own (its own start, termination and row-start path) and may overlap the others, the way a collar
    // overlaps a body; component n lies n yarn layers outside the piece's surface (scene/knit/guided_part.cpp).
    int knitComponent = 0;
};
// A setting's value ("" when the marking has none).
std::string markingSetting(const SurfaceMarking& m, const std::string& key);

// All markings of one imported file. The digest is meshDigest of the file's mesh as read: markings whose digest
// differs from the open file are stale (shown, never applied).
struct MarkingSet {
    int version = 5;            // the file version it was read from (1 to 5); the writer always writes the current one
    std::string source;
    // How the mesher prepares the marked pieces' reference (v3; stitchmesh/prepared_reference.h): "facets" (the imported
    // surface, the default and what v1 and v2 files mean) or "loop_limit" (the smooth-reference variant). The markings stay
    // bound to the file's triangles either way.
    std::string preparation = "facets";
    uint64_t digest = 0;
    std::vector<SurfaceMarking> markings;
    int nextId = 1;
};

const char* toString(MarkingMeaning m);  // "start", "finish", "row_start", "boundary", "cut", "attachment"
const char* toString(MarkingShape s);    // "point", "path", "closed"
const char* toString(MarkingOrigin o);   // "user", "proposal"
const char* toString(ProposalState p);   // "none", "offered", "accepted", "rejected"
bool parseMeaning(const std::string& text, MarkingMeaning& out);
bool parseShape(const std::string& text, MarkingShape& out);
bool parseOrigin(const std::string& text, MarkingOrigin& out);
bool parseProposalState(const std::string& text, ProposalState& out);

const char* meaningLabel(MarkingMeaning m);           // for people: "Start", "Row-start path", ...
MarkingShape defaultShape(MarkingMeaning m);          // the shape a new marking of this meaning starts with
bool shapeAllowed(MarkingMeaning m, MarkingShape s);  // Start, Finish: point or closed; RowStartPath: path; Cut: path
                                                      // or closed; Boundary, Attachment: closed
// Whether the marking has enough points for its shape: a point 1, a path 2, a closed curve 3.
bool markingComplete(const SurfaceMarking& m);

SurfaceMarking* findMarking(MarkingSet& s, int id);
const SurfaceMarking* findMarking(const MarkingSet& s, int id);

// The construction role a meaning plays for the mesher: one start, one termination and one row-start path per piece.
// Finish and Boundary are both a termination (the mesher takes one or the other); Cut and Attachment have no
// construction role, so a setup never sets them aside.
enum class MarkingRole { Start, Termination, RowStart, None };
MarkingRole markingRole(MarkingMeaning m);

// The ids of one setup's markings, in set order (setup > 0).
std::vector<int> setupMarkings(const MarkingSet& s, int setup);
// A setup rejected as a whole: its offered markings become Rejected (accepted or edited ones are left as they are).
void rejectSetup(MarkingSet& s, int setup);

} // namespace knit
