#pragma once
// Settle: a local geometric post-pass of the construction, and the certified repair.
//
// The placed yarn of a piece is right topologically (local signature) but can bend tighter than the
// yarn radius and come too close to itself next to shaping and poles; under a hand it is also
// irregular on purpose. Settle validates the construction, frees the defect stitches (with a hand:
// every stitch) plus `ring` neighbour rings, and relaxes them once (KnitNumerics::relax, the NCP
// problem of knit_problem.h: curvature, contact, edge lengths to the hand's rest, band); then it
// validates the result. Controls outside the free set do not move.
//
// One round, kept whole. The relax result is stored as one state, never mixed with the state before
// it: restoring the neighbourhoods of stitches whose local signature broke put back the
// construction's interpenetrating controls next to relaxed ones (knit-perf 2026-11-07: every final
// clearance FAIL of the Rustic teddy came from that restoration, and the restored stitches failed
// again in every later round). Further rounds over the remaining defects changed almost nothing (the
// relaxed state is feasible for every pair it represents; what remains are local signature failures
// that a motion without pass-through cannot change, and UNRESOLVED certifications).
//
// What this is NOT: a certified repair. Nothing proves that a strand did not pass through another
// during settling; that is why the result goes through gate (a) again (curvature, clearance, and the
// local signature, which a pass-through would break). Certified repair from a valid state is
// certifiedRepair below: one move, certified.
#include <string>
#include <vector>

#include "knitting/layout/base_surface.h"
#include "knitting/stitch/hand_tension.h"
#include "knitting/checks/checks.h"
#include "knitting/numerics.h"
#include "model/yarn_path.h"

namespace knit {

struct SettleParams {
    int iterations = 120;       // per round
    int refresh = 10;           // contact candidate refresh / band projection interval
    int ring = 1;               // free set: defect stitches plus this many neighbour rings
    double contactMargin = 0.05;  // [mm] beyond D + c
    double bendFactor = 1.15;   // minBendRadius = bendFactor * (r + c/2)
    double lengthStiffness = 0.3;
    double bandHalf = 2.0;      // [mm] max |offset along the surface normal|
    int samplesPerSpan = 4;
    ExecPolicy exec;            // the inner validations (identical results)
};

struct SettleReport {
    int rounds = 0;             // 1 when there was anything to relax, else 0
    int freeStitches = 0;       // the relaxed stitches (defects and their rings)
    double maxMove = 0;         // max control displacement from the input [mm]
    double ms = 0;
    std::vector<int> defectsPerRound;  // the defect stitches at the start (one entry)
};

// What settle and repair move: the yarn, the parameters that judge it, and where it must stay.
struct RelaxTarget {
    CurveSet& curves;
    const PatternParams& pp;                // D, c of the validation; repair: the target diameter
    const BaseSurface* surface = nullptr;   // null: no band
    // With a binding the band uses the ownership-preserving projection (base_surface.h): a control's
    // state advances only on Owned (otherwise the band is not applied to it and its result records
    // why), and every restore of controls restores their states with them.
    SurfaceBinding* bind = nullptr;
    // A hand's rest shape (HandResult::rest, hand_tension.h): every stitch is free, and edges rest on it.
    const Controls* handRest = nullptr;
};

// Settles the curves in place (one yarn or several; topology unchanged). Returns the final report
// through `last` (gate a on the settled curves).
SettleReport settle(const RelaxTarget& target, const SettleParams& sp, KnitNumerics& numerics, ValidationReport* last = nullptr);

// ----------------------------------------------------------- certified repair (b) ----
// One certified move from a verified start diameter to the target, in a fixed order with no loop:
//   1. The start diameter D0 is read off the report at the target (startDiameter): for a fixed
//      centreline the folds pass at D <= 2 rho_lo and every nonlocal distance d passes at
//      D <= d - c - deltaD (monotonicity in D), so D0 is the smallest of those bounds. One
//      validation at D0 verifies it (UNRESOLVED counts as failure); if it fails, or D0 < D0Min, or the
//      local signature fails (a signature FAIL fails at every diameter), nothing is repaired.
//   2. One relax call: the defect stitches at the target plus `ring` rings are free, the hard
//      thresholds are D0's, the soft targets the target's, every control's displacement is capped at
//      displacementFraction (D0 + c).
//   3. The moved state is validated at the target. It is kept only when it passes gate (a) there and
//      the straight motion from the start state (valid at D0) to it (valid at the target, so at D0)
//      is certified at D0 (certifyStraightMotion, displacement < D0 + c). Otherwise the curves and
//      the binding stay as settle left them.
// So repair makes no decision on a result it reads back in the middle of an iteration: each of the
// three steps runs at most once. (Before 2026-11-08 a bisection for D0, up to six move attempts per
// stage and up to 80 growth stages each decided on the host from a validation; on the teddy every
// repair that reached the target did so in one stage, and the one that did not made moves of
// 0.000 mm for 13 stages.)
// Certificate scope of an accepted move: no centreline passes through another for all inter-yarn
// pairs and same-yarn pairs farther apart than (pi/2) delta along the yarn; no tube contact during
// the motion only if delta <= c; the shorter same-yarn interval is NOT certified (diagnostic:
// curvature bound of the midpoint state).
struct RepairParams {
    double D0Min = 0.05;       // a start diameter below this is not tried [mm]
    int ring = 1;
    double displacementFraction = 0.5;  // the move's cap: fraction of D0 + c
    SettleParams settle;       // iterations, margins, band
};

// The one move of a repair.
struct RepairMove {
    bool accepted = false;     // passed at the target and certified: the curves are the moved ones
    double delta = 0;          // max control displacement of the move [mm]
    bool noPassThrough = false, noTubeContact = false;
    double localScale = 0;     // (pi/2) delta: uncertified same-yarn arc distance [mm]
    double midRhoLo = 0;       // curvature bound of the midpoint state (diagnostic)
    int freeStitches = 0;
    std::string why;           // why it was not kept
};

struct RepairReport {
    double D0 = 0, targetD = 0;  // D0: the verified start diameter (0: none)
    bool reachedTarget = false;  // the curves pass gate (a) at the target
    bool moveTried = false;
    RepairMove move;
    int validations = 0;
    double ms = 0;
    std::string detail;
};

// `current`: gate (a) of the curves as they are, at the target diameter (settle's last report), used
// instead of validating them again; null: repair validates them first.
RepairReport certifiedRepair(const RelaxTarget& target, const RepairParams& rp, KnitNumerics& numerics, ValidationReport* final = nullptr,
                             const ValidationReport* current = nullptr);

// Stitches touched by validation defects: clearance FAIL/UNRESOLVED regions, spans
// with kappa_hi > 1 / (r + c/2), local-signature failures.
std::vector<int> defectStitches(const CurveSet& cs, const ValidationReport& r, const PatternParams& pp);

} // namespace knit
