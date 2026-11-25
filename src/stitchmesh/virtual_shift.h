#pragma once
// The shift-path's virtual helical connectivity (W19 section 3.4: relaxation
// forces are computed "considering the modified mesh after the shift operation ... on-the-fly without actually
// modifying the topology").
//
// The stored mesh keeps closed rings. The shift-path is one wale edge per strip, (ring k vertex, ring k+1 vertex), a
// connected chain on the cut guide. The virtual mesh re-reads only the face whose W- edge is the path edge (slot 0 of
// each strip): with the frozen rule (decision D-3, PathSideDown) its two path-side corners are read one ring lower,
// so the path edge that face (k, last) leaves through becomes the edge face (k + 1, 0) enters through. Positions never
// change; only corner ids are re-read. Face (0, 0) has no ring below: its W- edge collapses to the cast-on path vertex
// (the start step); face (R - 1, last) leaves through a path edge with no face after it (the end step).
//
// A gathered pole start: strip 0 is a closed round of pole-start triangles with no path edge,
// so the rule leaves it as stored; strip 1 has no path edge below, so its slot 0 becomes the start step. The validator then
// expects that closed round plus one open strand through every other face, starting in strip 1.
//
// The validator never uses the rule: it rebuilds edge adjacency from the virtual polygons and traces the yarn, then
// compares the course interlocks and boundary course edges with the stored mesh's, re-derived the same way.
#include <string>
#include <vector>

#include "model/status.h"
#include "stitchmesh/typed_mesh.h"

namespace knit {

enum class ShiftReading : char {
    PathSideDown,   // the frozen rule: slot 0's path-side corners read one ring lower
    PathSideUp,     // negative control: one ring higher (the helix would climb from the bind-off)
    None            // negative control: the stored closed rings
};
const char* toString(ShiftReading r);

struct VirtualShift {
    ShiftReading reading = ShiftReading::PathSideDown;
    std::vector<std::vector<int>> corners;   // per face: the virtual corner ids (a start step repeats one vertex)
    std::vector<char> reRead;                // per face: some corner was re-read
    std::string rule;                        // the rule in words, for the report
};

// Needs strips, slots and the shift-path; a strip without exactly one path edge is left as stored (the validator then
// reports it). False only when the record has no strips.
bool buildVirtualShift(const TypedStitchMesh& s, ShiftReading reading, VirtualShift& out, std::string* err);

struct VirtualShiftCheck {
    Status status = Status::Unresolved;
    // the path itself
    int strips = 0, pathEdges = 0, stripsWithoutPathEdge = 0, stripsWithSeveralPathEdges = 0, pathBreaks = 0, pathEdgesNotWale = 0;
    // the yarn on the virtual polygons
    int strands = 0, openStrands = 0, closedStrands = 0, longestStrand = 0;
    bool oneStrandThroughAll = false, startsInFirstStrip = false, endsInLastStrip = false;
    bool poleStart = false;                  // strip 0 is all pole-start faces: the helix starts in strip 1
    int poleRoundFaces = 0;                  // faces of strip 0 traced as one closed strand (pole start only)
    int startSteps = 0, endSteps = 0;        // faces with a collapsed W- edge; faces whose W+ edge has nothing after it
    int startFace = -1, endFace = -1;
    int virtualFacesWithoutTwoWaleEdges = 0; // excluding the start step
    int nonManifoldVirtualEdges = 0;
    // interlocks and boundaries against the stored mesh
    int interlocksStored = 0, interlocksVirtual = 0, interlocksMissing = 0, interlocksExtra = 0;
    int boundaryBottomStored = 0, boundaryBottomVirtual = 0, boundaryTopStored = 0, boundaryTopVirtual = 0;
    bool countsUnchanged = false;            // faces, kinds and roles per face as stored
    std::vector<int> strandOrder;            // faces of the longest strand, in yarn order
    std::string detail;
    std::string text() const;
};

VirtualShiftCheck validateVirtualShift(const TypedStitchMesh& s, const VirtualShift& v);

// The virtual polygons as a typed mesh of their own (the same vertices and bindings): a repeated corner is dropped with
// its collapsed edge, so the start step becomes a triangle (B, W, T), typed Unsupported with reason "helix-start-step";
// every other face keeps its kind, roles and boundary roles. No rings (the count invariant is the stored mesh's).
TypedStitchMesh virtualTypedMesh(const TypedStitchMesh& s, const VirtualShift& v);

} // namespace knit
