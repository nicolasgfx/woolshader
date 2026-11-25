#pragma once
// Small knitted pieces with known answers, shared by knit_check, the viewer and the tests:
//   presets       parameter sets with the verdict each must produce
//   shaping demos
//     "inc"       flat patch, rows 8, 8, 9, 9, 9 (one INC), one yarn per row
//     "dec"       flat patch, rows 9, 9, 8, 8, 8 (one DEC), one yarn per row
//     "cone"      spiral on a surface of revolution whose rings hold n_r * w: 30 stitches, widening
//                 to 40 (INC turns) and back to 30 (DEC turns)
//   strand demo   one straight yarn, for fibre close-ups (knit_check --strand-demo)
#include <string>
#include <vector>

#include "knitting/layout/stitch_layout.h"
#include "model/params.h"
#include "model/status.h"

namespace knit {

struct Preset {
    std::string name;
    std::string description;
    PatternParams params;
    Status expected;  // verdict this preset must produce (asserted by tests)
};

const std::vector<Preset>& presets();
const Preset* findPreset(const std::string& name);

struct ShapingDemo {
    StitchProgram prog;
    PieceLayout lay;
    CurveSet curves;
    PatternParams pp;
};

bool makeShapingDemo(const std::string& kind, const ShapingParams& sh, ShapingDemo& out, const GaugeParams& g = {});

// One straight yarn of `length` mm along +x through the origin, diameter D; its fibre front is +z
// (the flat frame's fixed N).
struct StrandDemo {
    CurveSet curves;
    PatternParams pp;
};
StrandDemo makeStrandDemo(double length = 12.0, double D = 1.16);

} // namespace knit
