// The fibre stage's vocabulary: the preset, the report's bookkeeping and the setup of a part.
#include "fibres/fibre.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "fibres/fibre_islands.h"
#include "fibres/strands_file.h"
#include "knitting/tube/plies.h"
#include "model/hash.h"

namespace knit {

const char* toString(FibrePopulation p) { return p == FibrePopulation::Random ? "random" : "directed"; }

const char* toString(FibreOutcome o) {
    static const char* n[] = {"accepted", "back-facing", "hidden (anchor)", "ply inner (early)", "origin inside", "occluded",
                              "occluded by self", "probability", "too short to bend", "root depth", "ply inner (late)", "yarn end",
                              "root range", "no emergence (elevation)", "root not inside", "no emergence", "strand rules",
                              "turn too sharp", "no bend", "vertex inside", "approximation", "emergence", "no clearance",
                              "collision", "collision (truncated too short)", "approaches", "hidden (tip)", "behind the fabric", "not processed",
                              "serialized root check"};
    static_assert(sizeof n / sizeof n[0] == static_cast<size_t>(FibreOutcome::Count));
    const int k = static_cast<int>(o);
    return k >= 0 && k < static_cast<int>(FibreOutcome::Count) ? n[k] : "?";
}

namespace {

// True for the outcomes decided at the anchor, before any growth.
bool rejectedAtAnchor(FibreOutcome o) {
    switch (o) {
    case FibreOutcome::BackFacing:
    case FibreOutcome::HiddenAnchor:
    case FibreOutcome::PlyInnerEarly:
    case FibreOutcome::OriginInside:
    case FibreOutcome::Occluded:
    case FibreOutcome::OccludedSelf:
    case FibreOutcome::Probability: return true;
    default: return false;
    }
}

} // namespace

FibreParams pieceFibres() {
    FibreParams fp;
    // Candidates per mm of YARN (on plies too): with the 50 % mixture the body generates about
    // 1.33 M segments.
    fp.density = 125.0;
    fp.liftMinDeg = 12;  // shallower fibres need roots longer than their tubes allow
    fp.liftMaxDeg = 24;
    fp.lengthMean = 0.8;
    fp.lengthSpread = 0.35;
    fp.lengthMin = 0.4;
    fp.lengthMax = 1.6;
    fp.length = 0.8;
    fp.thickness = 0.07;     // hair diameter: library mm are close to real mm at this gauge
    fp.thicknessSigma = 0.04;
    fp.taper = 0.5;
    fp.frontCos = -0.1;      // anchors on the yarn's sides too
    fp.lipschitzSampling = true;
    fp.sampleStepFactor = 0.5;
    fp.keepProbability = 1.0;
    fp.curvatureWeight = 0.0;
    fp.randomFraction = 0.5;
    fp.densityPerParent = true;
    return fp;
}

void FibreSet::addFibre(Fibre f, std::span<const Vec3f> pts) {
    f.firstPoint = static_cast<uint32_t>(points.size());
    f.pointCount = static_cast<int>(pts.size());
    points.insert(points.end(), pts.begin(), pts.end());
    fibres.push_back(f);
}

void FibreSet::addFibre(Fibre f, const std::vector<Vec3>& pts) {
    f.firstPoint = static_cast<uint32_t>(points.size());
    f.pointCount = static_cast<int>(pts.size());
    for (const Vec3& p : pts) points.push_back(toVec3f(p));
    fibres.push_back(f);
}

void FibreSet::appendSet(const FibreSet& other) {
    const uint32_t base = static_cast<uint32_t>(points.size());
    points.insert(points.end(), other.points.begin(), other.points.end());
    fibres.reserve(fibres.size() + other.fibres.size());
    for (Fibre f : other.fibres) {
        f.firstPoint += base;
        fibres.push_back(f);
    }
}

std::vector<Vec3> fibrePoints(const FibreSet& fs, const Fibre& f) {
    std::vector<Vec3> pts;
    pts.reserve(static_cast<size_t>(f.pointCount));
    for (const Vec3f& p : fs.pointsOf(f)) pts.push_back(toVec3(p));
    return pts;
}

double fibreLength(const FibreSet& fs, const Fibre& f) {
    const std::span<const Vec3f> pts = fs.pointsOf(f);
    double L = 0;
    for (size_t i = 1; i < pts.size(); ++i) L += norm(toVec3(pts[i]) - toVec3(pts[i - 1]));
    return L;
}

double visibleLength(const FibreSet& fs, const Fibre& f) { return fibreLength(fs, f) - f.emergeL; }

double vertexRadius(const Fibre& f, int k) { return taperedRadius(f.radius0, f.radius1, k, std::max(1, f.pointCount - 1)); }

void FibreReport::count(FibreOutcome o) {
    switch (o) {
    case FibreOutcome::Accepted: ++accepted; break;
    case FibreOutcome::BackFacing: ++backFacing; break;
    case FibreOutcome::HiddenAnchor:
    case FibreOutcome::HiddenTip: ++hiddenInPiece; break;
    case FibreOutcome::PlyInnerEarly:
    case FibreOutcome::PlyInnerLate: ++plyInner; break;
    case FibreOutcome::OriginInside: ++originInside; break;
    case FibreOutcome::Occluded: ++occluded; break;
    case FibreOutcome::OccludedSelf:
        ++occluded;
        ++occludedBySelf;
        break;
    case FibreOutcome::Probability: ++probability; break;
    case FibreOutcome::TooShortToBend: ++tooShortToBend; break;
    case FibreOutcome::RootDepth:
    case FibreOutcome::RootRange:
    case FibreOutcome::RootNotInside: ++rootRejected; break;
    case FibreOutcome::YarnEnd: ++yarnEnd; break;
    case FibreOutcome::NoEmergenceElev:
    case FibreOutcome::NoEmergence: ++noEmergence; break;
    case FibreOutcome::RuleRejected: ++ruleRejected; break;
    case FibreOutcome::TurnTooSharp: ++turnTooSharp; break;
    case FibreOutcome::NoBend: ++noBend; break;
    case FibreOutcome::VertexInside: ++vertexInside; break;
    case FibreOutcome::ApproxError: ++approxError; break;
    case FibreOutcome::EmergeRejected: ++emergeRejected; break;
    case FibreOutcome::NoClearance: ++noClearance; break;
    case FibreOutcome::Collision: ++collisionRejected; break;
    case FibreOutcome::CollisionShort:
        ++collisionRejected;
        ++truncatedTooShort;
        break;
    case FibreOutcome::Approaches: ++approaches; break;
    case FibreOutcome::BehindFabric: ++behindFabric; break;
    case FibreOutcome::SerializedRejected: ++serializedRejected; break;
    case FibreOutcome::NotProcessed:
    case FibreOutcome::Count: break;
    }
}

const std::array<int FibreReport::*, kFibreReportCounters>& fibreReportCounters() {
    static const std::array<int FibreReport::*, kFibreReportCounters> fields = {
        &FibreReport::candidates, &FibreReport::backFacing, &FibreReport::originInside, &FibreReport::occluded,
        &FibreReport::occludedBySelf, &FibreReport::probability, &FibreReport::collisionRejected, &FibreReport::truncated,
        &FibreReport::accepted, &FibreReport::hiddenInPiece, &FibreReport::tooShortToBend, &FibreReport::yarnEnd,
        &FibreReport::approxError, &FibreReport::noBend, &FibreReport::noClearance, &FibreReport::approaches,
        &FibreReport::vertexInside, &FibreReport::turnTooSharp, &FibreReport::plyInner, &FibreReport::noEmergence,
        &FibreReport::emergeRejected, &FibreReport::rootRejected, &FibreReport::serializedRejected, &FibreReport::ruleRejected,
        &FibreReport::truncatedTooShort, &FibreReport::exposed, &FibreReport::behindFabric};
    return fields;
}

void FibreReport::addCandidate(FibrePopulation population, FibreOutcome outcome, bool truncatedByCut) {
    if (pop.size() < 2) pop.resize(2);
    FibreReport& own = pop[static_cast<size_t>(population)];
    for (FibreReport* r : {this, &own}) {
        ++r->candidates;
        if (!rejectedAtAnchor(outcome)) ++r->exposed;
        r->count(outcome);
        if (outcome == FibreOutcome::Accepted && truncatedByCut) ++r->truncated;
    }
}

void FibreReport::addCounters(const FibreReport& other) {
    for (int FibreReport::*const field : fibreReportCounters()) this->*field += other.*field;
    if (pop.size() < other.pop.size()) pop.resize(other.pop.size());
    for (size_t k = 0; k < other.pop.size(); ++k) pop[k].addCounters(other.pop[k]);
}

void resumFibreLengths(FibreSet& fs) {
    fs.rep.totalLength = 0;
    for (FibreReport& p : fs.rep.pop) p.totalLength = 0;
    for (const Fibre& f : fs.fibres) {
        const double length = visibleLength(fs, f);
        fs.rep.totalLength += length;
        const size_t k = static_cast<size_t>(f.population);
        if (k < fs.rep.pop.size()) fs.rep.pop[k].totalLength += length;
    }
}

namespace {

// Depth of the buried root below the smooth surface that keeps it inside a polygonal tube of
// `sides` with rings every `ring` [mm] around strands of radius `tubeRadius` whose curvature is at
// most `kappaMax`: polygon inset + chord sagitta + 0.02 mm. The curvature is clamped to the fold
// limit 1 / r: beyond it the tube is invalid anyway (FAIL, labelled residue).
double fibreRootDepth(double tubeRadius, int sides, double ring, double kappaMax) {
    kappaMax = std::min(kappaMax, 1.0 / tubeRadius);
    return tubeRadius * (1 - std::cos(kPi / sides)) + ring * ring * kappaMax * (1 + tubeRadius * kappaMax) / 8 + 0.02;
}

} // namespace

void setupPieceFibres(FibreParams& fp, FibreFrame& fr, const PatternParams& yarnPp, const FibreSetup& s) {
    fp.rootDepth = fibreRootDepth(s.tubeRadius, s.sides, s.ring, s.kappaMax);
    fr.exposureLength = 2 * (yarnPp.centrelineDepth() + yarnPp.yarnDiameter) + std::max(fp.lengthMax, fp.randomLengthMax);
    fr.backDepth = yarnPp.centrelineDepth() + yarnPp.radius();
    fr.pieceSeed = s.pieceIndex + 1;
    fr.candidates.clear();
    if (s.plies) {
        fr.plies = s.plies;
        fr.parentCurves = s.parentCurves;
        fr.parentSurface = s.parentSurface;
        if (fp.densityPerParent && s.parentCurves) fr.candidates = plyCandidateBudget(*s.plies, *s.parentCurves, fp.density);
    }
    fr.islands.clear();
    if (fp.islands.inUse() && s.strands && s.stitches)
        fr.islands = islandWeights({*s.strands, *s.stitches, fp.islands, fp.stitchArea, splitmix64(fp.seed ^ fr.pieceSeed)});
}

// Listed field by field; the size check makes an added field fail the build here until listed.
static_assert(sizeof(FibreParams) == 480, "FibreParams changed: update fibreParamsText");
std::string fibreParamsText(const FibreParams& p) {
    // Parameters of retired growth models keep their values here: every saved fibre result's key holds them.
    constexpr double curl = 0.8, tilt = 0.7, directionalVariation = 0.35, attachLength = 0.15, attachArc = 0.3;
    constexpr double segmentLength = 0, rootCollinearCos = 0, flatten = 0;
    constexpr int segments = 8, alongYarn = 1, anchorFrame = 0;
    std::string o;
    char b[512];
    std::snprintf(b, sizeof b, "%a %a %a %a %a %a %a %a %a %d %a %a %llu ", p.density, p.length, p.thickness, curl, tilt, directionalVariation,
             p.curvatureWeight, p.keepProbability, p.frontCos, segments, attachLength, attachArc,
             static_cast<unsigned long long>(p.seed));
    o += b;
    std::snprintf(b, sizeof b, "%a %a %a %a %a %a %a %d %a %d ", p.lengthMean, p.lengthMin, p.lengthMax, segmentLength, p.taper, p.rootDepth,
             rootCollinearCos, int(p.lipschitzSampling), flatten, alongYarn);
    o += b;
    std::snprintf(b, sizeof b, "%a %a %d %a %a %a %a %a %a %a %d %d %a ", p.helixDeg, p.helixSigmaDeg, p.handedness, p.liftMinDeg, p.liftMaxDeg,
             p.bowTheta, p.bowRadial, p.lengthSpread, p.thicknessSigma, p.segLenTarget, p.minSegments, p.maxSegments, p.maxCurveDev);
    o += b;
    std::snprintf(b, sizeof b, "%a %a %a %a %a %a %d %a %a %a ", p.minTurnDeg, p.maxTurnDeg, p.attachMax, p.sampleStepFactor, p.serializedMargin,
             p.exportChordFactor, anchorFrame, p.rootTurnDeg, p.minEmergeDeg, p.plyOuterCos);
    o += b;
    std::snprintf(b, sizeof b, "%a %a %a %a %a %a %a %a %a %a %d", p.randomFraction, p.randomElevMinDeg, p.randomElevMaxDeg, p.randomBend,
             p.randomLengthMean, p.randomLengthMin, p.randomLengthMax, p.randomLengthSpread, p.randomThickness,
             p.randomThicknessSigma, int(p.densityPerParent));
    o += b;
    if (p.lightChecks) o += " light";  // only when set: the keys of full checks are unchanged
    if (p.islands.inUse()) {           // only when set: the keys of even fuzz are unchanged
        std::snprintf(b, sizeof b, " islands %a %a %a %a %a %a %a", p.islands.rate, p.islands.size, p.islands.strength, p.islands.fuzzy,
                      p.islands.matted, p.islands.worn, p.stitchArea);
        o += b;
    }
    if (p.randomCurlRadius > 0) {      // only when set: the keys of plain strays are unchanged
        std::snprintf(b, sizeof b, " curl %a %a", p.randomCurlRadius, p.randomCurlLength);
        o += b;
    }
    if (p.randomSpreadDeg < 180) {     // only when set: the keys of strays in any direction are unchanged
        std::snprintf(b, sizeof b, " spread %a", p.randomSpreadDeg);
        o += b;
    }
    if (p.crimpAmplitude > 0) {        // only when set: the keys of uncrimped fibres are unchanged
        std::snprintf(b, sizeof b, " crimp %a %a", p.crimpAmplitude, p.crimpWavelength);
        o += b;
    }
    if (p.hugHeight > 0 || p.plyPitch > 0) {  // only when set: the keys of rising fibres are unchanged
        std::snprintf(b, sizeof b, " hug %a %a", p.hugHeight, p.plyPitch);
        o += b;
    }
    if (p.drape) o += " drape";      // only when set: the keys of cut fibres are unchanged
    return o;
}

} // namespace knit
