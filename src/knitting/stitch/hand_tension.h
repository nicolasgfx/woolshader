#pragma once
// A knitter's hand: every stitch a little different. Draws one record per stitch from its stable
// identity, then says what that does to the controls and to the yarn's rest length.
#include <cstdint>
#include <vector>

#include "knitting/numerics.h"
#include "model/choices.h"
#include "model/yarn_path.h"

namespace knit {

// One stitch's share of the hand. All zero: the template.
struct StitchTension {
    float tension = 0;  // how much tighter than the template, as a fraction: < 0 is looser
    float lean = 0;     // how far the stitch's head leans along U per mm up V
    float raise = 0;    // [mm] a raised loop: how far its front stands out along N
    float spot = 0;     // a thick (> 0) or thin (< 0) spot of the yarn centred here: change of diameter, fraction
};

constexpr double kMaxThickness = 1.5;  // FinalSurfaceSpec::supportRadius bounds a ring at 1 + maxExpansion

// One part's hand. settings_mapping.h handRequest fills it.
struct HandRequest {
    HandKnit settings;
    uint64_t seed = 0;  // HandKnit::seed mixed with the part's name
    double stitchArea = 18.7;  // [mm2] turns rates per cm2 into probabilities per stitch
    double wale = 5.5;         // [mm] the width of a stitch
    double maxThickness = kMaxThickness;  // the thickest ring the surface can draw, as a factor of r; 0: all round r
    bool inUse() const { return settings.inUse(); }
};

// One record per stitch, by stitch id. A pure function of the topology and the request, so a reload
// gets the same records. Keyed on the lattice (chart, i, r) or on (turn, index): a stitch's id changes
// between diameter attempts. Shaping and the ends of runs stay as they are.
std::vector<StitchTension> drawStitchTensions(const KnitTopology& topo, const HandRequest& request);

// What knitByHand leaves for settle.
struct HandResult {
    std::vector<StitchTension> tensions;  // the records, by stitch id
    Controls rest;                        // the rest shape settle's edges rest on (restOfTensions); empty: no hand
};

// How the rest shape follows from the records, for the knit cache keys: change it when the rule changes.
inline constexpr const char* kHandRestRule = "1 - t of the stitch before its shrink (2026-11-15)";

// Draws the records, moves the controls and sets the thickness. The rest shape is taken from the
// controls with every move but the tension's shrink: the shrink already makes a tight stitch smaller,
// and the rest says how much yarn it holds. Scaling the shrunk stitch again would count the tension twice.
HandResult knitByHand(CurveSet& cs, const HandRequest& request, KnitNumerics& numerics, const ExecPolicy& exec);

// Sets every yarn's thickness from the records. knitByHand and the knit restores call it.
void setYarnThickness(CurveSet& cs, const std::vector<StitchTension>& tensions, const HandRequest& request);

// The rest shape settle measures edges on: each stitch's edges scaled by 1 - its tension, integrated
// along the yarn, so no jump appears where two stitches meet. P: the controls before the shrink.
Controls restOfTensions(const Controls& P, const KnitTopology& topo, const std::vector<StitchTension>& tensions);

// How far plain stitches stray from their mean shape, in their own U, V, N frames [mm].
struct StitchSpread {
    int stitches = 0;
    double shape = 0, height = 0;  // RMS over all controls; height: along N only
};
StitchSpread measureStitchSpread(const CurveSet& cs);

// How many raised loops still stand: a stitch counts when its front stands at least half as far from
// the plain stitches' as it was raised.
struct RaisedLoops {
    int asked = 0, standing = 0;
};
RaisedLoops measureRaisedLoops(const CurveSet& cs, const std::vector<StitchTension>& tensions);

} // namespace knit
