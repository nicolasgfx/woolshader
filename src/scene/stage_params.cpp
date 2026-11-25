// From settings to knitting's and fibres' parameters; see stage_params.h.
#include "scene/stage_params.h"

#include <cstdio>

#include <algorithm>

#include "knitting/checks/checks.h"
#include "model/hash.h"
#include "settings/settings_mapping.h"
#include "settings/yarn_look.h"

namespace knit {

HandRequest handRequest(const SceneSettings& s, const std::string& partName) {
    HandRequest r{s.hand, splitmix64(static_cast<uint64_t>(s.hand.seed) ^ fnv1a64(partName)), s.wale * s.course, s.wale};
    // Flattening takes a thick ring at its contacts. Without it a thick ring meets its neighbour at the
    // clearance. Explicit plies have one radius, so a thick yarn would only spread them apart.
    r.maxThickness = s.plies > 1 ? 0 : yarnFit(s).compressed() ? kMaxThickness : 1 + s.clearance / s.diameter;
    return r;
}

FinalSurfaceSpec finalSurface(const SceneSettings& s, double D) {
    FinalSurfaceSpec f;
    f.D = D;
    f.c = s.clearance;
    f.compression = renderFlattened(s);  // strand.flatten or a compressed yarn (yarnFit)
    f.inflate = renderInflate(s);
    f.plies = s.plies;
    f.compressedYarn = yarnFit(s).compressed() && !s.flatten;
    f.maxCompression = s.maxFlattening;
    return f;
}

bool yarnBodyParams(const SceneSettings& s, YarnBodyParams& p, std::string* err) {
    const SceneSettings::YarnBody& b = s.yarnBody;
    if (s.plies > 1) {
        if (err) *err = "yarn_body.enabled needs strand.plies = 1: the yarn body makes its own plies (yarn_body.plies) on the single drawn yarn";
        return false;
    }
    p = YarnBodyParams{};
    p.plies = b.plies;
    p.fibres = b.fibres;
    p.interiorFibres = b.interiorFibres;
    p.fibreRadius = 0.5e-3 * b.fibreDiameterUm;  // [mm]
    p.shellDepth = b.shellDepth;
    p.fibreAngleDeg = b.fibreAngleDeg;
    p.plyAngleDeg = b.plyAngleDeg;
    p.irregularity = b.irregularity;
    p.core = b.core == "opaque" ? BodyCore::Opaque : b.core == "translucent" ? BodyCore::Translucent : BodyCore::None;
    // Baked maps write no core (addPartMappedTubes) and the bake never sees one (bakedMapKey): under maps the core
    // setting is no input, so its room inside the shell is not checked either.
    if (yarnBodyRenderOf(s) == YarnBodyRender::Maps) p.core = BodyCore::None;
    p.seed = static_cast<uint64_t>(b.seed);
    return true;
}

SurfaceParams yarnSurfaceParams(const SceneSettings& s) {
    SurfaceParams sp;
    sp.ringSpacing = kExportModelRing;
    sp.compression = renderFlattened(s);
    if (sp.compression) {
        sp.maxCompression = s.maxFlattening;
        sp.ringSpacing = 0.025;
    }
    return sp;
}

BakeParams bakeParams(const SceneSettings& s) {
    BakeParams p;
    p.texel = 1e-3 * s.yarnBody.mapTexelUm;  // [mm]
    return p;
}

double fillingOffset(const SceneSettings& s, double D) {
    PatternParams pp;  // the centreline's depth: A D, plus a rib's layer offset
    pp.yarnDiameter = D;
    pp.depthAmplitude = s.depthAmplitude;
    applyPatternShape(pp, s, parsePattern(s.pattern).value_or(StitchPattern::Stockinette));
    return pp.centrelineDepth() + 0.5 * finalSurface(s, D).renderedDiameter() + s.clearance + s.fillingMargin;
}

FurLook furLookFromSettings(const SceneSettings& s) {
    // The settle of the approved look (the furry-teddy prototype's groom, 2026-11-19): shape stiffness 0.75 at
    // the root and half that at the tip; 0.15 of gravity at 9.81 m/s^2 (a styling value, fur.h).
    constexpr double kStiffness = 0.75, kTipStiffnessRatio = 0.5, kGravityShare = 0.15, kStandardGravity = 9.81;
    const SceneSettings::Pile& p = s.pile;
    FurLook look;
    look.loops.density = p.density / 100;  // per cm^2 -> per mm^2
    look.loops.turns = p.length;
    look.loops.straightLength = p.straightLength;
    look.loops.rootDiameter = p.thickness;
    look.loops.taper = p.taper;
    look.loops.bend = p.bend;
    look.loops.curl = p.curl;
    look.variation.density = p.varDensity;
    look.variation.densityScale = p.varDensityScale;
    look.variation.length = p.varLength;
    look.variation.thickness = p.varThickness;
    look.variation.direction = p.varDirection;
    look.variation.bend = p.varBend;
    look.variation.curl = p.varCurl;
    look.fibres.enabled = p.fibresEnabled;
    look.fibres.density = p.fibresDensity / 100;
    look.fibres.lengthMedian = p.fibresLength;
    look.fibres.lengthSpread = p.fibresLengthSpread;
    look.fibres.rootDiameter = p.fibresThickness;
    look.fibres.taper = p.fibresTaper;
    look.fibres.liftDeg = p.fibresLiftDeg;
    look.fibres.direction = p.fibresDirection;
    look.islands = {p.islandsRate, p.islandsSize, p.islandsStrength, p.islandsMatted, p.islandsFluffy, p.islandsThin};
    look.irregular = {p.clumpSize, p.clumpStrength, p.patchLength, p.patchLift, p.wear};
    look.settle.strand.stiffness = kStiffness;
    look.settle.strand.tipStiffness = std::clamp(kStiffness * kTipStiffnessRatio, 0.0, 1.0);
    look.settle.strand.gravity = kGravityShare * (p.gravity / kStandardGravity);
    look.settle.gravityFullLength = p.gravityFullLength;
    look.seed ^= static_cast<uint64_t>(static_cast<uint32_t>(p.seed)) * 0x9E3779B97F4A7C15ull;
    return look;
}

bool fibreParamsFromSettings(const SceneSettings& s, FibreParams& fp, std::string* err) {
    const SceneSettings::Population &d = s.directed, &r = s.random;
    fp = pieceFibres();
    fp.lightChecks = s.fibreChecks == "light";
    fp.density = s.fibreDensity;
    fp.islands = s.islands;
    fp.stitchArea = s.wale * s.course;
    fp.frontCos = s.sideBandCos;  // the side band (pieceFibres' -0.1 by default)
    fp.randomFraction = s.randomFraction;
    fp.curvatureWeight = s.curvatureBias;
    fp.lengthMin = d.lengthMin;
    fp.lengthMax = d.lengthMax;
    fp.lengthMean = d.lengthMedian;
    fp.lengthSpread = d.lengthSpread;
    fp.thickness = d.thickness;
    fp.thicknessSigma = d.thicknessSigma;
    // Curvature: a scale of the default bow (directed) and the random bend amplitude.
    fp.bowTheta = FibreParams{}.bowTheta * d.curvature;
    fp.bowRadial = FibreParams{}.bowRadial * d.curvature;
    fp.liftMaxDeg = d.maxAngleDeg;  // the maximum NOMINAL lift
    fp.randomBend = r.curvature;
    fp.randomElevMaxDeg = r.maxAngleDeg;
    fp.randomCurlRadius = r.curlRadius;
    fp.randomCurlLength = r.curlLength;
    fp.randomSpreadDeg = s.randomSpreadDeg;
    fp.crimpAmplitude = s.crimpAmplitude;
    fp.crimpWavelength = s.crimpWavelength;
    // The random population's own lengths / thickness (inherit where equal: bit-identical).
    if (r.lengthMedian != d.lengthMedian) fp.randomLengthMean = r.lengthMedian;
    if (r.lengthMin != d.lengthMin) fp.randomLengthMin = r.lengthMin;
    if (r.lengthMax != d.lengthMax) fp.randomLengthMax = r.lengthMax;
    if (r.lengthSpread != d.lengthSpread) fp.randomLengthSpread = r.lengthSpread;
    if (r.thickness != d.thickness) fp.randomThickness = r.thickness;
    if (r.thicknessSigma != d.thicknessSigma) fp.randomThicknessSigma = r.thicknessSigma;
    if (!fibreSettingsValid(s, err)) return false;
    fp.maxSegments = s.fibreMaxSegments;
    fp.segLenTarget = s.fibreSegmentLength;
    // The polyline's allowed distance from the smooth curve grows with the chord's sagitta (0.45 mm: unchanged).
    const double chordRatio = s.fibreSegmentLength / FibreParams{}.segLenTarget;
    fp.maxCurveDev = FibreParams{}.maxCurveDev * chordRatio * chordRatio;
    if (chordRatio > 1) {
        // Over a stitch's outer bump the yarn curves away under a long chord: its first chord starts
        // flat, and the fibre turns more per vertex. At the 0.45 mm limits most of the fibres on the
        // visible side of the yarn would be refused, so they rise more and may start flat.
        fp.liftMinDeg = std::min(fp.liftMaxDeg, fp.liftMinDeg * chordRatio);
        fp.minEmergeDeg = 0;
        fp.maxTurnDeg = std::min(60.0, FibreParams{}.maxTurnDeg * chordRatio);
    }
    fp.hugHeight = s.directedHugHeight;
    fp.drape = s.directedDrape;
    if (s.directedFollowPlies) {
        // The ply helix's slope from the radius (setUpCurve); the draw keeps a small spread around it.
        fp.plyPitch = s.plyPitch;
        fp.handedness = s.plyHandedness;
        fp.helixDeg = 0;
        fp.helixSigmaDeg = 3;
    }
    return true;
}

PatternParams flatPatternParams(const SceneSettings& s) {
    PatternParams pp;
    pp.yarnDiameter = knitDiameter(s, s.diameter);  // a compressed yarn knits at k D and is drawn at D (renderInflate)
    pp.waleSpacing = s.wale;
    pp.courseSpacing = s.course;
    pp.depthAmplitude = s.depthAmplitude;
    pp.clearance = s.clearance;
    return pp;
}

PatternReference flatPatternReference(const SceneSettings& s, StitchPattern p) {
    PatternParams pp = flatPatternParams(s);
    applyPatternShape(pp, s, p);  // the pattern's relaxed shape, as the chart knit lays it out
    return flatPatternReference(pp, p);
}

} // namespace knit
