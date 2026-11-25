#pragma once
// The fibre stage's vocabulary: what controls fibres, what a fibre is and what a run reports.
// src/fibres/ARCHITECTURE.md gives the reading order and how this relates to the literature.
#include <array>
#include <cstdint>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "fibres/fibre_record.h"
#include "model/choices.h"
#include "model/mesh.h"
#include "model/parallel.h"
#include "model/params.h"
#include "model/yarn_path.h"

namespace knit {

// The yarn surface the fibres grow from (knitting/tube/yarn_surface.h, plies.h): read through pointers
// here, so the users of these headers do not compile knitting's.
struct PlySet;
struct SurfaceModel;
struct YarnSurface;

const char* toString(FibrePopulation p);

// What became of one candidate. Both generators report these values, in this order; the GPU kernels
// return them and store them as bytes (gpu/fibres/fibres.h outcomeCode).
enum class FibreOutcome : uint8_t {
    Accepted = 0, BackFacing, HiddenAnchor, PlyInnerEarly, OriginInside, Occluded, OccludedSelf, Probability,
    TooShortToBend, RootDepth, PlyInnerLate, YarnEnd, RootRange, NoEmergenceElev, RootNotInside, NoEmergence,
    RuleRejected, TurnTooSharp, NoBend, VertexInside, ApproxError, EmergeRejected, NoClearance, Collision,
    CollisionShort, Approaches, HiddenTip, BehindFabric, NotProcessed, SerializedRejected, Count
};
const char* toString(FibreOutcome o);

struct FibreParams {
    // --- candidates
    double density = 3.0;              // candidates per mm of yarn
    bool densityPerParent = false;     // plies: `density` per mm of the PARENT yarn (FibreFrame::candidates)
    double frontCos = 0.25;            // an anchor faces out if n . N exceeds this
    double curvatureWeight = 0.5;      // keep probability: 0 uniform, 1 strongly towards high centreline curvature
    double keepProbability = 0.7;      // artistic thinning after the exposure test
    uint64_t seed = 1;
    // --- length and thickness (a truncated log-normal; see fibreTargetLength)
    double length = 0.5;               // [mm] fixed length when lengthMean is 0
    double lengthMean = 0;             // [mm] median
    double lengthMin = 0, lengthMax = 0;
    double lengthSpread = 0;           // sigma of ln(length); 0: a clamped exponential of mean lengthMean
    double thickness = 0.012;          // [mm] diameter
    double thicknessSigma = 0;         // sigma of ln(radius), the normal clamped to +-2
    double taper = 1.0;                // tip radius / root radius
    // --- directed fibres: a curve in yarn coordinates (Zhao et al. 2016 Eq. 11)
    double helixDeg = 8, helixSigmaDeg = 10;  // nominal circumferential slope p = tan(psi), psi ~ N clamped +-2.5 sigma
    int handedness = 1;                       // sense of the helix (PlyGeometry::handedness)
    double liftMinDeg = 8, liftMaxDeg = 20;   // nominal radial slope q = tan(e), e ~ U
    double bowTheta = 0.10;                   // half-sine bow around the yarn [rad], |a| in [0.5, 1] x this
    double bowRadial = 0.03;                  // half-sine bow outwards [mm], a in [-1, 1] x this
    // Crimp (directed only): a planar wave across the fibre, amplitude [mm] and wavelength [mm] along it,
    // in a plane within 45 deg of the yarn surface; 0: none (fibre_shape.h DirectedCurve).
    double crimpAmplitude = 0, crimpWavelength = 1.5;
    double hugHeight = 0;                     // [yarn radii] > 0: the rise levels off at this height (DirectedCurve); 0: linear
    double plyPitch = 0;                      // [mm] > 0: psi follows the ply helix of this pitch, helixDeg is the noise's mean
    // Directed only, light checks only: a fibre that sinks deeper than the deep-contact limit into a strand is
    // pushed out over it (it rests on the stitch) instead of being cut there; cut only if the draped fibre fails.
    bool drape = false;
    // --- random fibres: flyaways, with their own start, bend and optional curl
    double randomFraction = 0;                // P(candidate is Random), drawn before any rejection
    double randomElevMinDeg = 20, randomElevMaxDeg = 60;  // elevation of the start tangent above the tangent plane
    double randomBend = 0.3;                  // u^2 bow amplitude / span, |a| in [0.5, 1] x this
    double randomSpreadDeg = 180;             // start azimuth within this of +-T (seen from above); 180: any
    double randomCurlRadius = 0, randomCurlLength = 8;    // > 0: a helix of this radius per turn of this length [mm]
    // Own length and thickness of the random population, drawn from the same hash channels; <= 0
    // (spread, sigma: < 0) inherits the value above. populationParams applies them.
    double randomLengthMean = 0, randomLengthMin = 0, randomLengthMax = 0, randomLengthSpread = -1;
    double randomThickness = 0, randomThicknessSigma = -1;
    // --- polyline
    double segLenTarget = 0.45;               // segments = clamp(round(L / this), minSegments, maxSegments)
    int minSegments = 2, maxSegments = 3;
    double maxCurveDev = 0.04;                // [mm] smooth curve to polyline
    double minTurnDeg = 3, maxTurnDeg = 35;   // world turn at every interior vertex
    // --- root and emergence
    double rootDepth = 0;                     // [mm] root below the surface (fibreRootDepth)
    double rootTurnDeg = 20;                  // turn between the root segment and the first chord
    double minEmergeDeg = 5;                  // the first chord rises at least this far above the surface
    double plyOuterCos = 0.5;                 // plies: anchor normal . outward-from-the-yarn must exceed this
    // --- read by the GPU kernels under fibres.checks = full only
    bool lightChecks = false;                 // the GPU's choice; the CPU generator always runs the light rules
    bool lipschitzSampling = false;           // own-surface samples tested with radius + spacing / 2
    double attachMax = 0.45;                  // [mm] a fibre must clear its own yarn within this
    double sampleStepFactor = 1.0;            // Lipschitz sample step = factor * radius
    double serializedMargin = 0;              // [mm] float32 envelope of the export, added to the fibre radius
    double exportChordFactor = 0;             // [mm^2] ringSpacing^2 / 8 of the export: facet bulge = this x kappa
    // --- patches of different fuzz (fibre_islands.h)
    FibreIslands islands;
    double stitchArea = 18.7;                 // [mm2] of the knit the islands sit on
    bool operator==(const FibreParams&) const = default;
};

// The calibrated preset for knitted parts (library mm, D = 1): short surface hairs filling the gaps
// between strands.
FibreParams pieceFibres();

struct FibreReport {
    int candidates = 0, backFacing = 0, originInside = 0, occluded = 0, occludedBySelf = 0, probability = 0;
    int collisionRejected = 0, truncated = 0, accepted = 0;
    int hiddenInPiece = 0;             // anchor or tip inside another piece
    int tooShortToBend = 0, yarnEnd = 0, approxError = 0, noBend = 0, noClearance = 0, approaches = 0;
    int vertexInside = 0;              // an interior vertex or the tip inside a tube
    int turnTooSharp = 0;              // an interior turn above maxTurnDeg (noBend: below minTurnDeg)
    int plyInner = 0;                  // plies: anchor on a ply side facing into the yarn
    int noEmergence = 0;               // segment 1 does not leave the own tube, or the first chord rises < minEmergeDeg
    int emergeRejected = 0;            // front / exposure / hidden test failed at the emergence
    int rootRejected = 0;              // root outside the yarn's arc range or not inside the own tube
    int serializedRejected = 0;        // buried / single exit / clearance on the serialized tubes
    int ruleRejected = 0;              // strand rules (segment >= 4.2 r) on the polyline
    int truncatedTooShort = 0;         // a truncation would leave fewer than minSegments
    int exposed = 0;                   // passed every anchor test (entered growth)
    int behindFabric = 0;              // cut at the fabric's back would leave fewer than minSegments
    double totalLength = 0;            // sum of accepted visible lengths [mm]
    double ms = 0;
    std::vector<FibreReport> pop;      // the same per population [Directed, Random]; candidates add up
    void count(FibreOutcome o);        // the counter(s) this outcome stands for
    // One candidate's outcome, in the total and in its population's report: candidates, exposed (it passed
    // the anchor tests), the outcome's counters, and truncated when a cut shortened an accepted fibre.
    // Both generators count through this.
    void addCandidate(FibrePopulation population, FibreOutcome outcome, bool truncatedByCut);
    void addCounters(const FibreReport& other);  // lengths and timing are not added
};
// The counters in their persisted order: the one list for adding, snapshots and the cache.
constexpr int kFibreReportCounters = 27;
const std::array<int FibreReport::*, kFibreReportCounters>& fibreReportCounters();

// A run's fibres: one record per fibre (fibre_record.h) and one packed pool of float32 points.
struct FibreSet {
    std::vector<Fibre> fibres;
    std::vector<Vec3f> points;
    FibreReport rep;
    std::span<const Vec3f> pointsOf(const Fibre& f) const { return {points.data() + f.firstPoint, static_cast<size_t>(f.pointCount)}; }
    // Adds f with these points, setting its firstPoint and pointCount.
    void addFibre(Fibre f, std::span<const Vec3f> pts);
    void addFibre(Fibre f, const std::vector<Vec3>& pts);
    void appendSet(const FibreSet& other);  // its fibres and points; the report stays
};

inline Vec3 toVec3(const Vec3f& p) { return {p.x, p.y, p.z}; }
inline Vec3f toVec3f(const Vec3& p) { return {float(p.x), float(p.y), float(p.z)}; }
std::vector<Vec3> fibrePoints(const FibreSet& fs, const Fibre& f);  // in double, for the checks' geometry
double fibreLength(const FibreSet& fs, const Fibre& f);             // the whole polyline
double visibleLength(const FibreSet& fs, const Fibre& f);           // from the emergence to the tip
// Radius at vertex k: taperedRadius from radius0 towards radius1 over the polyline (strands_file.h).
double vertexRadius(const Fibre& f, int k);
// FibreReport::totalLength (and per population) re-summed serially in fibre order.
void resumFibreLengths(FibreSet& fs);

// The point of the fabric's base surface nearest a query, and the unit outward normal there.
struct FabricPoint {
    Vec3 p, n;
};

// What one generation knows of the part it grows on.
struct FibreFrame {
    // The fabric's base surface near a point, the exposure ray's length and how far the fabric's back
    // lies behind the base surface. Null: a flat patch facing +z, its ray reaching the top of the set,
    // its back the lowest point of the tubes. May be called concurrently (exec). A function rather than a
    // mesh: the strand demo passes a plane and the tests chart surfaces; the GPU takes the part's mesh
    // instead (gpu::FibreFrameGpu::outward).
    std::function<FabricPoint(const Vec3&)> fabric;
    double exposureLength = 0;         // [mm] > 0 (required with `fabric`)
    double backDepth = 0;              // [mm] > 0 (required with `fabric`)
    uint64_t pieceSeed = 0;            // mixed into the seed: every spiral has yarnId 0
    // Plies: fibres anchored on them (cs, sm), directed ones authored in the PARENT yarn's
    // coordinates: plyParent(*plies, ply, ply u) gives (parent yarn, parent arc length). All three or none.
    const PlySet* plies = nullptr;
    const CurveSet* parentCurves = nullptr;
    const SurfaceModel* parentSurface = nullptr;
    std::vector<int> candidates;       // optional: candidates per strand (index yarnId), replacing floor(density L)
    std::vector<std::vector<IslandWeight>> islands;  // per strand, per control (fibre_islands.h); empty: none
    // Optional trace for comparing generators: per candidate, in strand and candidate order, its
    // FibreOutcome in bits 0..6 and the Random population in bit 7. No effect on results. Here rather than
    // a parameter: one test reads it (gpu_fibres_properties_and_cpu_agreement).
    std::vector<uint8_t>* outcomes = nullptr;
    ExecPolicy exec;
};

// The exported tube polygons; the viewers use them for the root depth too, so the preview buries
// roots as the export does.
constexpr int kExportSides = 8;
constexpr double kExportRing = 0.4;      // [mm]
constexpr double kExportModelRing = 0.2; // [mm] the surface model the tubes come from
constexpr int kExportPlySides = 6;
constexpr double kExportPlyRing = 0.3;   // [mm] (0.25 mm brought the GPU to 11.8 of 12 GB)

// The setup every caller shares, so preview and export derive root depth, exposure, seed, ply
// parents and the candidate budget identically. The caller adds the fabric's base surface.
struct FibreSetup {
    double tubeRadius = 0;   // radius of the supporting strands (ply r_p, or the yarn r)
    int sides = kExportSides;   // tube polygon of the export (or display)
    double ring = kExportRing;  // tube ring spacing [mm]
    double kappaMax = 0;     // curvature bound of the supporting strands [1/mm]
    size_t pieceIndex = 0;   // pieceSeed = pieceIndex + 1
    // Plies (optional): fibres grow on ps.curves, directed ones authored in the parent's
    // coordinates; with fp.densityPerParent the candidates are budgeted per parent yarn.
    const PlySet* plies = nullptr;
    const CurveSet* parentCurves = nullptr;
    const SurfaceModel* parentSurface = nullptr;
    // The strands as drawn and the yarn's stitches, for the islands; null: none.
    const CurveSet* strands = nullptr;
    const KnitTopology* stitches = nullptr;
};
// Sets fp.rootDepth, fr's exposure length (2 (A D + D) + the longest fibre, of the YARN pattern
// `yarnPp`) and back depth (A D + D / 2: the stitches' depth plus a yarn radius), piece seed, ply parent
// fields, candidate budget and island weights.
void setupPieceFibres(FibreParams& fp, FibreFrame& fr, const PatternParams& yarnPp, const FibreSetup& s);

// Every field as text, doubles as exact hex floats (cache keys).
std::string fibreParamsText(const FibreParams& p);

} // namespace knit
