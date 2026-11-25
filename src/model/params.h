#pragma once
// What a knit asks for: the stitch pattern, the gauge and yarn (PatternParams), and a pattern's flat reference
// verdict (PatternReference). Settings (settings_mapping) and the knitting fixtures make PatternParams; knitting,
// fibres, gpu, scene and export read it. Knitting makes PatternReference (flatPatternReference); jobs and the
// viewer show it.
//
// Fabric-local coordinates: x along U (course direction: the way a row runs),
// y along V (wale direction: up a column of stitches), z along N = U x V,
// the designated FRONT normal (technical/knit face). The back is -N.
#include <optional>
#include <string>

#include "model/status.h"

namespace knit {

// Knit / purl patterns. The face of stitch (wale, course) is a function of
// its LATTICE coordinates only (never of a yarn or run index): see stitchFaceIsPurl.
enum class StitchPattern { Stockinette, ReverseStockinette, Garter, Rib1x1, Rib2x2 };

inline const char* toString(StitchPattern p) {
    switch (p) {
    case StitchPattern::Stockinette: return "stockinette";
    case StitchPattern::ReverseStockinette: return "reverse-stockinette";
    case StitchPattern::Garter: return "garter";
    case StitchPattern::Rib1x1: return "rib1x1";
    case StitchPattern::Rib2x2: return "rib2x2";
    }
    return "?";
}

// The pattern named s ("garter"), nullopt for an unknown name.
inline std::optional<StitchPattern> parsePattern(const std::string& s) {
    for (StitchPattern q : {StitchPattern::Stockinette, StitchPattern::ReverseStockinette, StitchPattern::Garter,
                            StitchPattern::Rib1x1, StitchPattern::Rib2x2})
        if (s == toString(q)) return q;
    return std::nullopt;
}

inline bool stitchFaceIsPurl(StitchPattern p, int wale, int course) {
    auto mod = [](int a, int m) { return ((a % m) + m) % m; };
    switch (p) {
    case StitchPattern::Stockinette: return false;
    case StitchPattern::ReverseStockinette: return true;
    case StitchPattern::Garter: return mod(course, 2) == 1;
    case StitchPattern::Rib1x1: return mod(wale, 2) == 1;
    case StitchPattern::Rib2x2: return mod(wale, 4) >= 2;
    }
    return false;
}

// Rib 1x1 and 2x2: the patterns whose wales alternate faces (PatternParams::ribLayerOffset applies).
inline bool patternHasRibLayers(StitchPattern p) { return p == StitchPattern::Rib1x1 || p == StitchPattern::Rib2x2; }

struct PatternParams {
    int wales = 6;    // stitches per course (along U)
    int courses = 6;  // number of courses / rows (along V)

    // Authoritative spacings [mm]. Densities are derived: density = 1/spacing.
    double waleSpacing = 8.0;    // stitch width along U  [mm]
    double courseSpacing = 6.5;  // row height along V    [mm]

    double yarnDiameter = 1.0;   // nominal diameter D [mm]; collision radius r = D/2
    double clearance = 0.1;      // requested minimum surface gap c [mm]

    // Template shape (normalized by the spacings) and depth.
    double depthAmplitude = 0.8;  // control depth +-A, in units of D
    double headHalfWidth = 0.375; // loop head half width at the shoulder, in stitch widths
    double legHalfWidth = 0.125;  // leg half width inside the old head, in stitch widths

    double sampleTolerance = 0.005;  // polyline sampling bound eps [mm]

    StitchPattern pattern = StitchPattern::Stockinette;
    // The garter transition (stockinette.h): stitches whose own face differs from their
    // child's flip depth between the leg (parent interlock) and the neck (child
    // interlock). Neck control j = 12 at height transNeckY with depth transNeckZ x
    // sigma_child (0: on the mid-plane), leg control j = 13 at height transLegY. The
    // plain template's values would be (0.60, 1, 0.30), which FAIL (depth gap 0.97 mm at
    // 8 x 6.5, D 1.0); (0.60, 0.4, 0.40) is the best of a measured grid search
    // (clearance 1.153 mm at 8 x 6.5 / D 1.0, 1.047 mm at 5.5 x 3.4 / D 1.0).
    double transNeckY = 0.60, transNeckZ = 0.4, transLegY = 0.40;
    // Rib layers (rib 1x1 and 2x2 only; stockinette.h): a relaxed rib is two-layered, its knit wales in
    // front of its purl wales (Spencer, Knitting Technology; Niu, Dion, Kamien 2024). Every control of a
    // knit wale moves +ribLayerOffset x D along N, of a purl wale -ribLayerOffset x D; the sinker between
    // a knit and a purl wale stays on the mid-plane. 0: one layer (the template before 2026-11).
    double ribLayerOffset = 0.0;  // [D]
    // Rib 2x2 pairs (rib 2x2 only; stockinette.h): the two wales of a pair keep the stockinette wale
    // spacing between them, and the knit pair overlaps the purl pair behind it at the pair boundaries
    // (the relaxed 2x2 rib draws in to ~0.65 of stockinette's width; Hwang et al. 2025). The pair's
    // first wale moves -ribPairSpread wale spacings along U, its second +ribPairSpread; the sinker
    // between two wales lies halfway between their moved centres. 0: evenly spaced wales.
    double ribPairSpread = 0.0;  // [wale spacings]
    // Garter ridges (garter only; stockinette.h): the relaxed garter fabric folds like an accordion, one
    // fold per row (Kaldor, James, Marschner 2008; Niu, Dion, Kamien 2024), which lets its rows be
    // ~0.70 of stockinette's (patternRestScale). Every control moves along N by
    // -garterWave x D x sigma_own x cos(pi Y): one smooth wave over two rows, a function of the
    // position in the fabric only, so every crossing keeps its depth order and gap; its slope cancels
    // the depth flip between a garter stitch's leg and neck. 0: flat rows (before 2026-11).
    double garterWave = 0.0;  // [D]

    double radius() const { return 0.5 * yarnDiameter; }
    // How far the yarn centreline reaches in front of / behind the mid-plane: the template's depth
    // A D, plus a rib's layer offset [mm].
    double centrelineDepth() const { return (depthAmplitude + (patternHasRibLayers(pattern) ? ribLayerOffset : 0.0)) * yarnDiameter; }

    bool operator==(const PatternParams&) const = default;
};

// Structured failure categories of a validation report (failCategories in knitting/checks/validate.cpp; tests
// assert these, never the wording of a detail line).
enum FailCategory : unsigned {
    kFailCurvature = 1u,    // regular / no local fold
    kFailSignature = 2u,    // crossing order of an intermeshed pair, or local signature
    kFailCrossing = 4u,     // unexpected crossings, depth gap at crossings
    kFailClearance = 8u,
    kFailGenerator = 16u,   // periodicity / interior signatures (flat generator only)
};
// The categories in a mask of FailCategory bits, "curvature+clearance", or "none". Shown by the viewer's knit
// panel for PatternReference::fail and ::unresolved.
inline std::string failCategoryText(unsigned failMask) {
    std::string s;
    const char* names[] = {"curvature", "signature", "crossing", "clearance", "generator"};
    for (int b = 0; b < 5; ++b)
        if (failMask & (1u << b)) s += (s.empty() ? "" : "+") + std::string(names[b]);
    return s.empty() ? "none" : s;
}

// The flat reference verdict of a knit / purl pattern (knitting/checks/checks.h, flatPatternReference).
// The generator (stockinette.h) says what is supported and a scene's own gates judge that scene;
// between the two, this verdict holds for the parameter tuple it was computed at (params) and nothing
// else: not for other slider values, not for a curved scene.
struct PatternReference {
    StitchPattern pattern = StitchPattern::Stockinette;
    Status verdict = Status::Unresolved;
    unsigned fail = 0, unresolved = 0;
    double minInterYarn = 0;  // [mm]
    std::string params;       // the tuple the verdict holds for, e.g. "D 1.10, 5.50 x 3.40, A 0.80, c 0.10"
};

} // namespace knit
