#pragma once
// Reference-bound marking records of the imported-arm experiment,
// namespace knit::refmark. Historical: the shared user-facing contract is the Mark step's file (model/surface_marking.h,
// io/marking_file.h); these records stay as fixtures and are migrated to it by stitchmesh_lab propose-marks (checked).
//
// A record set names exactly one component of an imported model and exactly one prepared reference mesh of it. Surface
// positions are (triangle id, barycentric coordinates) on that reference; the reference's identity is a digest of its
// vertices and triangles, so a changed triangulation no longer matches and its markings are refused until they are
// rebound by a checked closest-point rebinding (rebindMarkings). Coordinates of the reference are library mm (the import's
// space: library = scale * zUp(file) + offset), and the record states that frame explicitly.
//
// Semantics, kept apart on purpose:
//   start_point       where the component's rows start (a pole start, construction setting `construction`)
//   finish_boundary   the bind-off: the level contour of the row-time field through this point
//   growth_direction  from the start marking to the finish marking (settings from=, to=)
//   row_start_path    the computational path where every row starts (Method A's ring vertex 0 and the helix's
//                     shift-path); it does NOT cut the surface
//   cut               a physical cut of the fabric (an opening or a component split); not consumed by the mesher yet
//   attachment_region triangles where another component overlaps (setting other=); a candidate for a join, never a join
// Provenance is one of manual | automatic_proposal | user_edited. There is no "approved" value: a proposal stays a
// proposal until a person edits or re-enters it.
//
// The file is plain text, one statement per line (`stockinette-markings 1` first); see writeMarkings.
#include <string>
#include <utility>
#include <vector>

#include "model/mesh.h"
#include "model/status.h"

namespace knit::refmark {

struct SurfacePoint {
    int tri = -1;
    double bary[3] = {1, 0, 0};
};

enum class MarkingType : char { StartPoint, FinishBoundary, GrowthDirection, RowStartPath, Cut, AttachmentRegion };
enum class Provenance : char { Manual, AutomaticProposal, UserEdited };
const char* toString(MarkingType t);
const char* toString(Provenance p);

struct Marking {
    std::string id;
    MarkingType type = MarkingType::StartPoint;
    Provenance provenance = Provenance::AutomaticProposal;
    std::vector<SurfacePoint> points;    // a point marking: one; a path: in order, consecutive points sharing a triangle
    std::vector<int> triangles;          // a region marking
    bool closed = false;                 // a path or boundary that closes on itself
    std::vector<std::pair<std::string, std::string>> settings;   // construction and semantic settings, in order
    std::string note;                    // free text (why it was proposed, what was edited)
    std::string setting(const std::string& key, const std::string& fallback = "") const;
};

// The original component: the source file, the import that prepared it, the node (an OBJ object of the file), the
// material and the connected piece within that node and material (pieces numbered by their lowest file triangle).
struct ComponentIdentity {
    std::string name;                    // a label for people ("white bear, left arm")
    std::string sourceDigest;            // the converted OBJ the import read
    std::string importIdentity;          // importIdentity() of that import, 16 hex digits
    std::string node, material;
    int piece = -1;
    bool operator==(const ComponentIdentity&) const = default;
};

struct ReferenceIdentity {
    std::string digest;                  // referenceDigest() of the prepared reference mesh
    int vertices = 0, triangles = 0;
    std::string frame;                   // e.g. "library mm, Z up; library = 68.1 * zUp(file) + (x, y, z)"
};

struct MarkingSet {
    int version = 1;
    ComponentIdentity component;
    ReferenceIdentity reference;
    std::vector<Marking> markings;
    const Marking* find(MarkingType t) const;   // the first of a type, or null
    const Marking* byId(const std::string& id) const;
};

bool writeMarkings(const std::string& path, const MarkingSet& s, std::string* err);
bool readMarkings(const std::string& path, MarkingSet& s, std::string* err);

// FNV-1a over the vertex positions (printed %.6f) and the triangles: the reference's identity.
std::string referenceDigest(const TriMesh& m);
Vec3 pointPosition(const TriMesh& m, const SurfacePoint& p);
// A triangle containing both points (on its edges or inside), or -1: the triangle a path segment lies in.
int commonTriangle(const TriMesh& m, const SurfacePoint& a, const SurfacePoint& b);

struct MarkingCheck {
    Status status = Status::Unresolved;   // FAIL: refused (identity, binding or semantics); UNRESOLVED: listed open items; PASS
    std::vector<std::string> refusals;    // why the set cannot be used at all
    std::vector<std::string> open;        // what the mesher will not consume (a cut), proposals not yet edited
    int points = 0, proposals = 0, edited = 0, manual = 0;
    std::string text() const;
};
// The set against the reference it claims and the component the caller expects.
MarkingCheck validateMarkings(const MarkingSet& s, const TriMesh& reference, const ComponentIdentity& expected);

// A path made to lie on the surface segment by segment: wherever two consecutive points share no triangle, the straight
// segment between them is traced across the triangles (each crossing of a triangle edge becomes a point). For points that
// lie on one plane piece of the surface (a rebinding onto a refinement, an editor's stroke inside a facet); false when a
// segment cannot be traced (it leaves the surface), with the path unchanged.
bool densifyPath(const TriMesh& m, std::vector<SurfacePoint>& path, int* inserted);

struct RebindReport {
    int points = 0, refused = 0, inserted = 0;   // inserted: points added by densifyPath
    double maxDistance = 0;               // [mm] from a point's old position to the new reference
    std::string detail;
};
// Rebinding onto another triangulation of the same surface: every point goes to the closest point of newRef; refused (false)
// when one moves farther than tolerance. Region triangles cannot be rebound this way and are refused. On success `out`
// names the new reference (digest, counts) and keeps ids, types, provenance and settings, with a note of the rebinding.
bool rebindMarkings(const MarkingSet& in, const TriMesh& oldRef, const TriMesh& newRef, double tolerance, MarkingSet& out, RebindReport& report);

} // namespace knit::refmark
