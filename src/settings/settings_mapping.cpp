// From settings to the numbers the stages use; see settings_mapping.h.
#include "settings/settings_mapping.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>

#include "model/choices.h"
#include "model/coordinates.h"
#include "fibres/fibre.h"
#include "settings/file_presets.h"

namespace knit {

bool fibreSettingsValid(const SceneSettings& s, std::string* err) {
    // The limits the settings are checked against are the fibre model's defaults (pieceFibres).
    const SceneSettings::Population &d = s.directed, &r = s.random;
    const FibreParams base = pieceFibres();
    auto bad = [&](const std::string& why) {
        if (err) *err = why;
        return false;
    };
    if (d.maxAngleDeg < base.liftMinDeg) return bad("directed max angle below the minimum lift " + std::to_string(base.liftMinDeg));
    if (r.maxAngleDeg < base.randomElevMinDeg)
        return bad("random max angle below the minimum start elevation " + std::to_string(base.randomElevMinDeg));
    if (!(d.lengthMin > 0 && d.lengthMin <= d.lengthMedian && d.lengthMedian <= d.lengthMax) ||
        !(r.lengthMin > 0 && r.lengthMin <= r.lengthMedian && r.lengthMedian <= r.lengthMax))
        return bad("fibre lengths need 0 < min <= median <= max");
    if (r.curlRadius < 0 || (r.curlRadius > 0 && !(r.curlLength > 0)))
        return bad("fibres.random.curl_radius and curl_length: a radius of 0 or more, and a positive length per turn");
    if (s.crimpAmplitude < 0 || (s.crimpAmplitude > 0 && !(s.crimpWavelength > 0)))
        return bad("fibres.crimp_amplitude and crimp_wavelength: an amplitude of 0 or more, and a positive wavelength");
    // The wool look's export values, which the wool profiles own (file_presets.cpp fibreKeys).
    if (!(s.tailTaper >= 0 && s.tailTaper <= 1) || !(s.tailTaperMinDiameters >= 0 && s.tailTaperMinDiameters <= 100) ||
        !(s.tailTaperTip > 0 && s.tailTaperTip <= 1))
        return bad("fibres.tail_taper 0..1, tail_taper_min_diameters 0..100, tail_taper_tip above 0 up to 1");
    if (!(s.normalStrength > 0 && s.normalStrength <= 4)) return bad("strand.normal_strength above 0 up to 4");
    if (!(s.yarnWhite.translucency >= 0 && s.yarnWhite.translucency <= 1) || !(s.yarnBrown.translucency >= 0 && s.yarnBrown.translucency <= 1))
        return bad("material.yarn_*.translucency 0..1");
    if (!(s.directedHugHeight >= 0 && s.directedHugHeight <= 1))
        return bad("fibres.directed.hug_height: 0 (off) to 1 yarn radius");
    if (s.directedFollowPlies && !(s.plyPitch > 0)) return bad("fibres.directed.follow_plies needs strand.ply_pitch > 0");
    if (!(s.randomSpreadDeg >= 5 && s.randomSpreadDeg <= 180))
        return bad("fibres.random.spread_deg: 5 to 180 degrees (180 = any direction around the yarn)");
    // Segments per fibre: at least minSegments (2), at most the GPU's polyline bound (root + anchor +
    // segments <= kMaxFibreVerts = 50, src/gpu/fibres_flat.h), so both backends accept every value.
    if (s.fibreMaxSegments < base.minSegments || s.fibreMaxSegments > 48)
        return bad("fibres.max_segments = " + std::to_string(s.fibreMaxSegments) + ": 2 to 48 segments per fibre");
    if (!(s.fibreSegmentLength >= 0.1 && s.fibreSegmentLength <= 4))
        return bad("fibres.segment_length = " + std::to_string(s.fibreSegmentLength) + ": 0.1 to 4 mm");
    return true;
}

std::vector<std::string> furExcludedParts(const SceneSettings& s) {
    std::vector<std::string> out;
    std::stringstream ss(s.pile.excludedParts);
    std::string item;
    while (std::getline(ss, item, ','))
        if (!item.empty()) out.push_back(item);
    return out;
}

void applySceneDefaults(SceneSettings& s, const SceneEntry& e) {
    // The Dense preset (thickness_presets.h): drawn 1.86, knitted compressed; the fallbacks are the
    // calibrated 1.12 / 1.08 / 1.04 / 1.0 scaled by 1.86 / 1.16 (knitted at about those diameters).
    s.diameter = 1.86;
    s.diameterFallback = {1.80, 1.73, 1.67, 1.60};
    // Knitted as spiral pieces (the teddy): its gauge (its Dense preset).
    if (e.knitsAsPieces()) {
        s.thicknessPreset = "thick-dense";
        s.diameter = e.piecesD;
        s.diameterFallback = e.piecesDFallback;
        s.wale = e.piecesW;
        s.course = e.piecesH;
    }
    if (!e.defaultYarnColourSource.empty()) s.yarnColourSource = e.defaultYarnColourSource;
    if (e.defaultFibres) s.fibres = true;
    if (!e.defaultFibrePreset.empty())
        if (const FilePreset* p = presetWithId(PresetGroup::Fibres, e.defaultFibrePreset)) applyPreset(PresetGroup::Fibres, s, *p);
    // The uv lattice's seams overlapped (knit.coverage, scene_settings.h): spiral pieces have no
    // chart borders to cover.
    if (!e.knitsAsPieces()) s.coverage = "overlap";
    applySceneSettings(s, e);
}

void applySceneSettings(SceneSettings& s, const SceneEntry& e) {
    for (const bool sizeMode : {true, false})
        for (const auto& [key, value] : e.defaultSettings) {
            if ((key == "strand.stitch_size") != sizeMode) continue;
            std::string err;
            if (!setSetting(s, key, value, &err)) std::fprintf(stderr, "scene %s: default setting %s = %s refused: %s\n", e.id.c_str(), key.c_str(), value.c_str(), err.c_str());
        }
}

bool parseTextureProfile(const std::string& id, int& plies, PlyProfileModel& model, std::string* err) {
    if (id == "none") {
        plies = 1;
        return true;
    }
    const size_t dash = id.find('-');
    const std::string n = id.substr(0, dash), m = dash == std::string::npos ? "" : id.substr(dash + 1);
    char* end = nullptr;
    const long v = std::strtol(n.c_str(), &end, 10);
    if (end == n.c_str() || *end || v < 2 || v > 8 || (m != "swept" && m != "circles")) {
        if (err) *err = "texture profile " + id + ": expected none or <n>-swept / <n>-circles with n in 2..8";
        return false;
    }
    plies = static_cast<int>(v);
    model = m == "swept" ? PlyProfileModel::Swept : PlyProfileModel::Circles;
    return true;
}

DiameterPlan partDiameters(const SceneSettings& s) {
    DiameterPlan plan;
    plan.attempt.push_back(knitDiameter(s, s.diameter));
    for (double d : s.diameterFallback)
        if (d < s.diameter) plan.attempt.push_back(knitDiameter(s, d));
    if (s.diameterAttempts > 0 && static_cast<int>(plan.attempt.size()) > s.diameterAttempts) {
        plan.notAttempted.assign(plan.attempt.begin() + s.diameterAttempts, plan.attempt.end());
        plan.attempt.resize(s.diameterAttempts);
    }
    return plan;
}

SceneMapping sceneMetres(const Vec3& centreMm) {
    SceneMapping m;
    m.sigma = 1000;
    m.centre = centreMm / 1000;
    return m;
}

FitLimits fitLimits() { return {}; }

PatternScale patternRestScale(StitchPattern p) {
    // Hwang et al. 2025, supplemental Table 2, measured patch size over stockinette's (61.77 x 59.41 mm):
    // garter-1x1 67.79 x 41.75, rib-1x1 45.63 x 73.00, rib-2x2 40.23 x 66.66 (course x wale, mm).
    switch (p) {
    case StitchPattern::Stockinette:
    case StitchPattern::ReverseStockinette: return {1.0, 1.0};
    case StitchPattern::Garter: return {1.10, 0.70};
    case StitchPattern::Rib1x1: return {0.74, 1.23};
    case StitchPattern::Rib2x2: return {0.65, 1.12};
    }
    return {};
}

StitchSize patternStitch(const SceneSettings& s, StitchPattern p) {
    const PatternScale k = patternRestScale(p);
    return {s.wale * k.w, s.course * k.h};
}

void applyPatternShape(PatternParams& pp, const SceneSettings& s, StitchPattern p) {
    const StitchSize st = patternStitch(s, p);
    pp.pattern = p;
    pp.waleSpacing = st.w;
    pp.courseSpacing = st.h;
    if (p == StitchPattern::Garter) {
        // The garter ridges and transition (knit_check flat 8 x 8 scan at 6.05 x 2.38..3.74 / D 1.164).
        pp.garterWave = kGarterWave;
        pp.transNeckY = kGarterNeckY;
        pp.transNeckZ = kGarterNeckZ;
        pp.transLegY = kGarterLegY;
        return;
    }
    if (!patternHasRibLayers(p)) return;
    // The loops keep (about) stockinette's width: in stitch widths of the narrower rib wale they are
    // wider. Rib 1x1 up to kRibMaxWiden, so that every head stays inside its stitch's sinkers (each
    // crossing then lies on its own stitch's span); rib 2x2 at stockinette's width, its pairs spread
    // to stockinette's wale spacing (a pair's two wales share a layer).
    const double widen = 1.0 / patternRestScale(p).w;
    const double loop = p == StitchPattern::Rib1x1 ? std::min(widen, kRibMaxWiden) : widen;
    pp.headHalfWidth = PatternParams{}.headHalfWidth * loop;
    pp.legHalfWidth = PatternParams{}.legHalfWidth * loop;
    pp.ribLayerOffset = p == StitchPattern::Rib1x1 ? kRibLayerOffset1x1 : kRibLayerOffset2x2;
    if (p == StitchPattern::Rib2x2) pp.ribPairSpread = 0.5 * (widen - 1.0);
}

namespace {

bool tightnessFactors(const std::string& tightness, double& mW, double& mH) {
    if (tightness == "tight") mW = mH = 1;
    else if (tightness == "standard") mW = 6.5 / 1.1 / kFitW, mH = 4.6 / 1.1 / kFitH;  // former Medium: 1.0 mm at 6.5 x 4.6
    else if (tightness == "loose") mW = 7.5 / 0.9 / kFitW, mH = 5.5 / 0.9 / kFitH;     // former Thin: 0.8 mm at 7.5 x 5.5
    else return false;
    return true;
}

bool stitchSizeFromYarn(const SceneSettings& s, double& w, double& h) {
    double mW, mH;
    if (!tightnessFactors(s.tightness, mW, mH)) return false;
    const FitLimits L = fitLimits();
    const double dc = s.diameter + s.clearance;
    auto up = [](double x) { return std::ceil(x / 0.05 - 1e-9) * 0.05; };  // 0.05 mm steps, never below x (beyond 1e-9 steps)
    w = up(L.w * mW * dc);
    h = up(L.h * mH * dc);
    return true;
}

} // namespace

void applyStitchSize(SceneSettings& s) {
    if (s.stitchSize != "yarn") return;
    double w, h;
    if (stitchSizeFromYarn(s, w, h)) s.wale = w, s.course = h;
}

YarnFit yarnFit(const SceneSettings& s) {
    YarnFit f;
    if (!(s.diameter > 0 && s.wale > 0 && s.course > 0)) return f;
    const FitLimits L = fitLimits();
    f.need = (std::min(s.wale / L.w, s.course / L.h) - s.clearance) / s.diameter;
    if (f.need >= 1) return f;
    if (!s.compressToFit) {
        f.tooTight = true;
        return f;
    }
    f.k = std::max(std::round(f.need * 1000) / 1000, kMinCompression);
    f.tooTight = f.need < kMinCompression;
    return f;
}

double knitDiameter(const SceneSettings& s, double D) { return D * yarnFit(s).k; }

bool renderFlattened(const SceneSettings& s) { return s.flatten || yarnFit(s).compressed(); }

double renderInflate(const SceneSettings& s) {
    if (!renderFlattened(s)) return 0;
    const double e = s.flatten ? s.flattenInflate : 0.0, k = yarnFit(s).k;
    return k == 1 ? e : (1 + e) / k - 1;  // not compressed: e exactly (the inflation as set)
}

bool tubeSettingsValid(const SceneSettings& s, std::string* why) {
    auto bad = [&](const std::string& m) {
        if (why) *why = m;
        return false;
    };
    for (const auto& [key, sides] : {std::pair<const char*, int>{"strand.tube_sides", s.tubeSides}, {"strand.ply_tube_sides", s.plyTubeSides}})
        if (sides < 3 || sides > 64) return bad(std::string(key) + " = " + std::to_string(sides) + " (3 to 64 sides)");
    for (const auto& [key, ring] : {std::pair<const char*, double>{"strand.tube_ring", s.tubeRing}, {"strand.ply_tube_ring", s.plyTubeRing}})
        if (!(ring >= 0.05 && ring <= 2.0)) return bad(std::string(key) + " = " + settingNumber(ring) + " mm (0.05 to 2 mm)");
    return true;
}

bool overlapCoverage(const SceneSettings& s) { return knitWorkflow(s) && s.coverage == "overlap"; }

bool pileSettingsValid(const SceneSettings& s, std::string* why) {
    const SceneSettings::Pile& p = s.pile;
    auto bad = [&](const std::string& m) {
        if (why) *why = m;
        return false;
    };
    struct Bound {
        const char* key;
        double value, lo, hi;
        const char* unit;
    };
    const Bound bounds[] = {
        {"pile.density", p.density, 1, 2000, "fibres per cm^2"},
        {"pile.length", p.length, 0.05, 10, "loop lengths"},
        {"pile.straight_length", p.straightLength, 0.1, 50, "mm"},
        {"pile.thickness", p.thickness, 0.01, 3, "mm"},
        {"pile.taper", p.taper, 0.05, 1, "of the root diameter"},
        {"pile.bend", p.bend, 0.1, 4, "quarter turns"},
        {"pile.curl", p.curl, 0, 10, "mm"},
        {"pile.gravity", p.gravity, 0, 50, "m/s^2"},
        {"pile.gravity_full_length", p.gravityFullLength, 0, 500, "mm"},
        {"pile.var.density", p.varDensity, 0, 1, ""},
        {"pile.var.density_scale", p.varDensityScale, 0.5, 500, "mm"},
        {"pile.var.length", p.varLength, 0, 1, ""},
        {"pile.var.thickness", p.varThickness, 0, 1, ""},
        {"pile.var.direction", p.varDirection, 0, 1, ""},
        {"pile.var.bend", p.varBend, 0, 1, ""},
        {"pile.var.curl", p.varCurl, 0, 1, ""},
        {"pile.fibres.density", p.fibresDensity, 0, 2000, "fibres per cm^2"},
        {"pile.fibres.length", p.fibresLength, 0.1, 50, "mm"},
        {"pile.fibres.length_spread", p.fibresLengthSpread, 0, 2, ""},
        {"pile.fibres.thickness", p.fibresThickness, 0.005, 1, "mm"},
        {"pile.fibres.taper", p.fibresTaper, 0.05, 1, "of the root diameter"},
        {"pile.fibres.lift_deg", p.fibresLiftDeg, 0, 90, "deg"},
        {"pile.fibres.direction", p.fibresDirection, 0, 1, ""},
        {"pile.islands.rate", p.islandsRate, 0, 10, "per cm^2"},
        {"pile.islands.size", p.islandsSize, 1, 100, "mm"},
        {"pile.islands.strength", p.islandsStrength, 0, 1, ""},
        {"pile.islands.matted", p.islandsMatted, 0, 100, ""},
        {"pile.islands.fluffy", p.islandsFluffy, 0, 100, ""},
        {"pile.islands.thin", p.islandsThin, 0, 100, ""},
        {"pile.irregular.clump_size", p.clumpSize, 0.5, 20, "mm"},
        {"pile.irregular.clump_strength", p.clumpStrength, 0, 1, ""},
        {"pile.irregular.patch_length", p.patchLength, 0, 1, ""},
        {"pile.irregular.patch_lift", p.patchLift, 0, 1, ""},
        {"pile.irregular.wear", p.wear, 0, 1, ""},
    };
    for (const Bound& b : bounds)
        if (!(b.value >= b.lo && b.value <= b.hi))
            return bad(std::string(b.key) + " = " + settingNumber(b.value) + ": " + settingNumber(b.lo) + " to " + settingNumber(b.hi) + (*b.unit ? " " : "") + b.unit);
    // A loop is a fibre bent round a circle: its smallest radius (curl variation) must exceed its thickest fibre's
    // radius (thickness variation), else the tube folds into itself (the exporter's curvature rule).
    const double smallestLoop = p.curl * (1 - 0.5 * p.varCurl), thickestRadius = 0.5 * p.thickness * (1 + p.varThickness);
    if (p.curl > 0 && !(smallestLoop > thickestRadius))
        return bad("pile.curl = " + settingNumber(p.curl) + " mm: its smallest loop (radius " + settingNumber(smallestLoop) +
                   " mm with pile.var.curl) is not wider than the thickest fibre's radius (" + settingNumber(thickestRadius) + " mm with pile.thickness)");
    return true;
}

} // namespace knit
