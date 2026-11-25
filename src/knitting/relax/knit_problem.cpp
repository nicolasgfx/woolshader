// The out-of-line part of the problem's shared types (knit_problem.h).
#include "knitting/relax/knit_problem.h"

#include <algorithm>

namespace knit {

const char* toString(RowFamily f) {
    switch (f) {
    case RowFamily::Edge: return "edge";
    case RowFamily::Curvature: return "curvature";
    case RowFamily::Contact: return "contact";
    case RowFamily::BandUpper: return "band+";
    case RowFamily::BandLower: return "band-";
    case RowFamily::Cap: return "cap";
    case RowFamily::ContactTarget: return "contact target";
    case RowFamily::CurvatureTarget: return "curvature target";
    case RowFamily::Anchor: return "anchor";
    }
    return "?";
}

double KnitProblemScore::maxConstraint() const {
    double m = 0;
    for (int f = 0; f < kRowFamilies; ++f)
        if (isHard(static_cast<RowFamily>(f))) m = std::max(m, maxViol[f]);
    return m;
}

std::string KnitProblemScore::line() const {
    char b[512];
    std::snprintf(b, sizeof b,
                  "objective %.6g, strain max %.4f rms %.4f; max violation K %.3g (%lld) C %.3g (%lld) B %.3g (%lld) M %.3g (%lld); "
                  "target shortfall C %.3g (%lld) K %.3g (%lld); degenerate %lld, not owned %lld; coverage %.4f mm%s",
                  objective, strainMax, strainRms, maxViol[1], nViol[1], maxViol[2], nViol[2], std::max(maxViol[3], maxViol[4]),
                  nViol[3] + nViol[4], maxViol[5], nViol[5], maxViol[6], nViol[6], maxViol[7], nViol[7], degenerate, notOwned, coverage,
                  covered ? "" : " (EXCEEDED: pairs may be missing)");
    return b;
}

} // namespace knit
