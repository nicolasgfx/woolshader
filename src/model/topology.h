#pragma once
// Knit topology: which stitch is drawn through which.
//
// Every stage that needs knitting structure (signature, contact classes, relaxation
// anchors) reads it from here instead of assuming one yarn per course and a wale
// index, so the same code serves the flat patch (one yarn per course) and a spiral
// piece (one yarn, many turns, same-yarn interlocks).
//
// A stitch record is one loop: a parameter range [uBegin, uEnd) on its yarn, its
// turn (= course of a flat patch) and index in that turn (= wale), and the ordered
// loops it is drawn through (its parents; 0 for the first turn, 1 for K and INC
// children, 2 for DEC). Terminal records (CAST ring, CLOSE drawstring) have
// variable arity and are kept separate.
#include <utility>
#include <vector>

#include "model/vec.h"

namespace knit {

enum class StitchOp { K, Inc, Dec };
const char* toString(StitchOp op);

// The face a loop is drawn towards. Knit: through the parent from back to front
// (the loop shows on +N, the technical face); Purl: from front to back (the same loop
// mirrored in N). expected_signature.h derives the purl crossing order from this.
enum class StitchFace : unsigned char { Knit, Purl };
const char* toString(StitchFace f);

// Lattice stitches (chart_lattice.h): chart, lattice column i and row r, and the run (one open
// yarn of consecutive kept nodes in a row). Pattern parity and parents use (chart, i, r) only;
// turn = run for yarn construction. chart = -1: not a lattice stitch.
struct LatticeCoord {
    int chart = -1, i = 0, r = 0, run = -1;
};

struct StitchRecord {
    int id = -1;
    int yarn = 0;
    double uBegin = 0, uEnd = 0;   // half-open ownership on the yarn
    int turn = 0, index = 0;       // flat patch: course, wale
    StitchOp op = StitchOp::K;
    int parents[2] = {-1, -1};     // ordered along the parent turn
    int nParents = 0;
    int children[2] = {-1, -1};    // stitches drawn through this loop (INC: 2), in order
    int nChildren = 0;
    int incSlot = 0;               // INC children: 0 / 1 (order in the turn)
    StitchFace face = StitchFace::Knit;  // K stitches only; INC / DEC are knit
    LatticeCoord lattice;              // tiled scenes
    // Local frame at the layout node (template origin): U along the course in
    // knitting order, V up the wale (towards the own head), N = U x V (front).
    Vec3 origin, U{1, 0, 0}, V{0, 1, 0}, N{0, 0, 1};
};

enum class TerminalOp { Cast, Close };

struct TerminalRecord {
    TerminalOp op = TerminalOp::Cast;
    int yarn = 0;
    double uBegin = 0, uEnd = 0;
    std::vector<int> stitches;  // CAST: the first-turn stitches worked into the ring; CLOSE: the loops gathered
    // CLOSE: the drawstring's u-range through each gathered loop (aligned with stitches)
    std::vector<std::pair<double, double>> segments;
    Vec3 centre, axis;          // pole point and outward pole axis
    double radius = 0;          // ring radius [mm]
};

// Fills children[] from parents[] (ids must equal indices).
void linkChildren(std::vector<StitchRecord>& stitches);

struct KnitTopology {
    std::vector<StitchRecord> stitches;   // id == index
    std::vector<TerminalRecord> terminals;

    bool empty() const { return stitches.empty(); }
    // (Re)builds the per-yarn u index; call after filling stitches. yarnEnd[y] is the
    // yarn's uEnd: a stitch ending there also owns u == uEnd (as YarnPath does).
    void buildIndex(const std::vector<double>& yarnEnd);
    // Stitch owning u on the yarn (half-open, the last stitch owns its end), or -1.
    int stitchAt(int yarn, double u) const;
    // The last stitch starting at or before u, or before that the yarn's first; -1: the yarn has
    // no stitch. Unlike stitchAt it answers on the padding and the tail too.
    int stitchNear(int yarn, double u) const;
    // Stitch containing [u0, u1] entirely, or -1.
    int stitchContaining(int yarn, double u0, double u1) const;
    bool isParent(int parent, int child) const;
    // Stitch owning u, or for a CLOSE drawstring segment the loop it passes through.
    int stitchOrLoopAt(int yarn, double u) const;

private:
    std::vector<std::vector<std::pair<double, int>>> byYarn_;  // sorted (uBegin, id)
    std::vector<double> yarnEnd_;
};

} // namespace knit
