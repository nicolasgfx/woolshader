#pragma once
// Spiral stitch program of one knitted piece.
//
// Knitting in the round consumes the loops of the previous turn in order. The
// program is built turn by turn from GENERATIONAL queues, so a turn never mixes
// generations (a DEC never takes a loop of the previous turn together with one of
// the current turn):
//   turn 0          n_0 stitches without parents (open start, or worked into a CAST
//                   ring: the terminal record says which);
//   turn r >= 1     input I_r = the output loops of turn r - 1, in order; ops
//                     K    pop 1, produce 1
//                     INC  pop 1, produce 2 adjacent stitches sharing that parent
//                          (atomic: both children belong to the same op)
//                     DEC  pop 2 consecutive loops of I_r, produce 1 (never on the
//                          last loop of I_r)
//                   balance |I_r| = K + INC + 2 DEC, n_r = K + 2 INC + DEC.
// The spiral (shift-path) transition needs no extra connectivity: stitch ids follow
// the yarn, so the last stitch of turn r and the first of turn r + 1 are yarn
// neighbours, and the first stitch of turn r + 1 consumes the first loop of turn r.
#include <string>
#include <vector>

#include "model/status.h"
#include "model/topology.h"

namespace knit {

struct TurnPlan {
    int input = 0;   // |I_r| (0 for turn 0)
    int k = 0, inc = 0, dec = 0;
    int output = 0;  // n_r
};

struct StitchProgram {
    std::vector<TurnPlan> turns;
    std::vector<StitchRecord> stitches;  // ids in yarn order; frames and u ranges unset
    std::vector<int> turnStart;          // first stitch id of each turn
    bool castRing = false;               // turn 0 is worked into a CAST ring
    bool closeDrawstring = false;        // the last turn is gathered by CLOSE
};

// Feasible counts n_0 .. n_{R-1} for real targets t_r (t_0 is ignored: n_0 is given):
// n_r = round(t_r) limited to |n_r - n_{r-1}| <= floor(n_{r-1} * maxShaping), so the
// shaping ops of one turn (all INC or all DEC) fit the density limit; the deficit is
// carried by the next turn's comparison with the actual count.
std::vector<int> scheduleCounts(int n0, const std::vector<double>& targets, double maxShaping = 0.5);

// The program for the given counts. Shaping positions are spread Bresenham-style
// with a per-turn phase (golden-ratio sequence) so they do not stack in one wale.
StitchProgram buildProgram(const std::vector<int>& counts, bool castRing = false, bool closeDrawstring = false);

struct GraphCheck {
    Status status = Status::Pass;
    int turns = 0, stitches = 0, inc = 0, dec = 0;
    std::string detail;  // first violation, or a summary
};
// Combinatorial checks: per-turn balance; every loop of turn r - 1 consumed exactly
// once by turn r (an INC parent by exactly its two adjacent children); DEC parents
// consecutive in the same input turn; turn-0 stitches parentless; ids in yarn order;
// with closeDrawstring, every loop of the last turn is gathered (no loop left).
GraphCheck checkProgram(const StitchProgram& p);

} // namespace knit
