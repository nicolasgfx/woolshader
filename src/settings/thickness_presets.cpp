// The thickness presets of the Knit step: yarn diameter and stitch size together.
#include "settings/settings_mapping.h"
#include "settings/thickness_presets.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace knit {

std::vector<ThicknessPreset> thicknessPresets(const SceneEntry&) {
    // Every scene, the teddy's spiral pieces included: the teddy's calibrated stitch (5.5 x 3.4)
    // with three drawn yarns. Dense and Medium are thicker than the stitch fits: knitted compressed
    // at the contact diameter (yarnFit, 1.16 mm for stockinette) and drawn at D, so all three knit
    // the same centrelines. Thin is the calibrated round yarn.
    return {
        {"thick-dense", "Dense", 1.86, 5.5, 3.4,
         "Chunky yarn: 1.86 mm drawn at 5.5 x 3.4 mm, knitted compressed at 1.16 mm and flattened where yarns touch. The "
         "teddy's spiral pieces pass after settle and repair (14 of 18; the ears fail); a flat swatch fails the clearance.",
         {1.80, 1.73, 1.67, 1.60}},
        {"medium", "Medium", 1.50, 5.5, 3.4,
         "1.50 mm drawn at 5.5 x 3.4 mm: the same knit as Dense, a thinner yarn. The teddy's spiral pieces pass 14 of 18.",
         {1.45, 1.40, 1.34, 1.29}},
        {"thin", "Thin", 1.16, 5.5, 3.4,
         "The calibrated round yarn: 1.16 mm at 5.5 x 3.4 mm, not compressed. The teddy's spiral pieces pass after settle and "
         "repair (14 of 18; the ears fail).",
         {1.12, 1.08, 1.04, 1.0}},
    };
}

const ThicknessPreset* matchingPreset(const SceneSettings& s, const SceneEntry& e) {
    static thread_local std::vector<ThicknessPreset> list;
    list = thicknessPresets(e);
    for (const ThicknessPreset& p : list)
        if (s.thicknessPreset == p.id && std::fabs(s.diameter - p.D) < 1e-9 && std::fabs(s.wale - p.w) < 1e-9 && std::fabs(s.course - p.h) < 1e-9)
            return &p;
    return nullptr;
}

void applyThicknessPreset(SceneSettings& s, const ThicknessPreset& p) {
    s.thicknessPreset = p.id;
    s.stitchSize = "custom";  // a preset is a stitch size of its own
    s.diameter = p.D;
    s.wale = p.w;
    s.course = p.h;
    if (!s.diameterFallback.empty()) s.diameterFallback = p.fallback;  // Draft keeps none
}

std::vector<GaugeChoice> gaugeChoices(const SceneEntry& e) {
    std::vector<GaugeChoice> list;
    for (const ThicknessPreset& p : thicknessPresets(e)) list.push_back({p.id, p.label, p.note, ""});
    if (!list.empty()) list.front().group = "Same stitch, different yarn";
    // The tightness factors of stitchSizeFromYarn (settings_mapping.cpp) against the smallest stitch.
    list.push_back({"yarn-tight", "Snug",
                    "The smallest stitch a round yarn of this diameter fits, so the yarn is never compressed. The stitch width and "
                    "row height follow the yarn diameter.",
                    "Stitch sized for the yarn"});
    list.push_back({"yarn-standard", "Regular",
                    "A roomier stitch for this yarn: 1.36 x the smallest stitch's width and 1.56 x its height. The stitch width and "
                    "row height follow the yarn diameter.",
                    ""});
    list.push_back({"yarn-loose", "Loose",
                    "A loose stitch for this yarn: 1.92 x the smallest stitch's width and 2.28 x its height. The stitch width and "
                    "row height follow the yarn diameter.",
                    ""});
    list.push_back({"custom", "Custom", "Your own yarn diameter, stitch width and row height (the fields below).", "Your own values"});
    return list;
}

std::string currentGaugeChoice(const SceneSettings& s, const SceneEntry& e) {
    if (s.stitchSize == "yarn") return s.tightness == "tight" || s.tightness == "loose" ? "yarn-" + s.tightness : "yarn-standard";
    const ThicknessPreset* p = matchingPreset(s, e);
    return p ? p->id : "custom";
}

void applyGaugeChoice(SceneSettings& s, const SceneEntry& e, const std::string& id) {
    for (const ThicknessPreset& p : thicknessPresets(e))
        if (p.id == id) {
            applyThicknessPreset(s, p);
            return;
        }
    s.thicknessPreset = "custom";
    if (id.rfind("yarn-", 0) == 0) {
        s.stitchSize = "yarn";
        s.tightness = id.substr(5);
        applyStitchSize(s);  // w and h follow D, c and the tightness
    } else {
        s.stitchSize = "custom";  // Custom: the current values stay
    }
}

std::string yarnFitNote(const SceneSettings& s) {
    const YarnFit f = yarnFit(s);
    if (!f.compressed()) return {};
    char b[400];
    std::snprintf(b, sizeof b,
                  "The stitch is too tight for a round %.2f mm yarn: knitted at a %.2f mm contact diameter (yarn compressed to %.0f %%) "
                  "and drawn at %.2f mm, %s.",
                  s.diameter, f.k * s.diameter, 100 * f.k, s.diameter,
                  s.plies > 1 ? "its round plies crossing where yarns touch" : "flattened where yarns touch");
    return b;
}

std::string thickYarnWarning(const SceneSettings& s, const SceneEntry&) {
    const YarnFit f = yarnFit(s);
    if (!f.tooTight) return {};
    // The widest yarn this stitch fits as a round tube (the stitch needs fitLimits() (D + c)).
    const FitLimits L = fitLimits();
    const double fitD = std::min(s.wale / L.w, s.course / L.h) - s.clearance;
    char least[64];
    std::snprintf(least, sizeof least, "even compressed to %.0f %%", 100 * kMinCompression);
    const char* why = !s.compressToFit ? "compression is off" : least;
    char b[512];
    std::snprintf(b, sizeof b,
                  "Yarn %.2f mm does not fit a %.2f x %.2f mm stitch (a round yarn fits up to %.2f mm; %s): parts are expected "
                  "to FAIL and knitting takes long. Use a thinner yarn or a larger stitch.",
                  s.diameter, s.wale, s.course, std::max(fitD, 0.0), why);
    return b;
}

} // namespace knit
