#pragma once
// formatKnitReport: the text report of a knitted scene (parts, verdicts, clearance).
//
// From a Scene and its settings. Called by the viewer's report and knit_check. Cannot fail.
#include <string>
#include <vector>

#include "knitting/checks/checks.h"
#include "scene/scene.h"

namespace knit {

std::string formatKnitReport(const Scene& scene, const SceneSettings& settings, bool verbose);

// A guided part's outcome in lines for people (the report, the Knit step's details, the export's report.json): what was
// made from what, how each stage ended, generated against retained yarn, and the limits that hold for every guided knit.
std::vector<std::string> guidedPartLines(const GuidedPart& g);

// Diagnostics (-v): spans bending tighter than r + c/2 and clearance FAIL regions, grouped by op and
// template span, to stdout.
void printDefects(const CurveSet& cs, const ValidationReport& part, const PatternParams& pp);

} // namespace knit
