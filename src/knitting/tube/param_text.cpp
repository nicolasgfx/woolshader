#include "knitting/tube/param_text.h"

#include <cstdio>

#include "model/hash.h"

namespace knit {

// Listed field by field; the size checks make an added field fail the build here until listed
// (MSVC x64 layout).
static_assert(sizeof(PatternParams) == 128, "PatternParams changed: update patternParamsText");
static_assert(sizeof(SurfaceParams) == 40, "SurfaceParams changed: update surfaceParamsText");


std::string patternParamsText(const PatternParams& p) {
    char b[512];
    std::snprintf(b, sizeof b, "%d %d %a %a %a %a %a %a %a %a %d %a %a %a", p.wales, p.courses, p.waleSpacing, p.courseSpacing,
                  p.yarnDiameter, p.clearance, p.depthAmplitude, p.headHalfWidth, p.legHalfWidth, p.sampleTolerance, int(p.pattern),
                  p.transNeckY, p.transNeckZ, p.transLegY);
    // Appended only where they shape the pattern: every other knit's key stays as it was.
    std::string text = b;
    char extra[48];
    if (patternHasRibLayers(p.pattern) && p.ribLayerOffset != 0) {
        std::snprintf(extra, sizeof extra, " rib %a", p.ribLayerOffset);
        text += extra;
    }
    if (p.pattern == StitchPattern::Rib2x2 && p.ribPairSpread != 0) {
        std::snprintf(extra, sizeof extra, " spread %a", p.ribPairSpread);
        text += extra;
    }
    if (p.pattern == StitchPattern::Garter && p.garterWave != 0) {
        std::snprintf(extra, sizeof extra, " wave %a", p.garterWave);
        text += extra;
    }
    return text;
}

std::string surfaceParamsText(const SurfaceParams& p) {
    // Stored surfaces are keyed by this text. The 1 and " fit multi" are the retired preserveArea and
    // fit fields at the only values production used; they stay so that the keys do not change.
    char b[512];
    std::snprintf(b, sizeof b, "%a %d %d %a %a %a", p.ringSpacing, int(p.compression), 1, p.maxCompression, p.maxExpansion, p.smoothing);
    return std::string(b) + (p.compression ? " fit multi" : "");
}

uint64_t surfaceBuildDigest(const CurveSet& cs, const PatternParams& pp, const SurfaceParams& sp) {
    char b[32];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(curvesDigest(cs)));
    return fnv1a64(std::string("surface build v1\ncurves = ") + b + "\npattern = " + patternParamsText(pp) + "\nsurface = " +
                   surfaceParamsText(sp) + "\n");
}

uint64_t curvesDigest(const CurveSet& cs) {
    std::string s;
    char b[512];
    for (const YarnPath& y : cs.yarns) {
        s.append(reinterpret_cast<const char*>(y.curve.ctrl.data()), y.curve.ctrl.size() * sizeof(Vec3));
        s.append(reinterpret_cast<const char*>(y.thickness.data()), y.thickness.size() * sizeof(float));
        for (const StitchSpan& sp : y.stitches) {
            std::snprintf(b, sizeof b, "%d %d %a %a %d;", sp.wale, sp.course, sp.uBegin, sp.uEnd, sp.id);
            s += b;
        }
        s += "|";
    }
    return fnv1a64(s);
}

} // namespace knit
