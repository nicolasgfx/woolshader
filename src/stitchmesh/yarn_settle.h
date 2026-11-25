#pragma once
// The yarn settle of the stitch-mesh experiment: an
// embedded yarn (yarn_embed.h) converted to production curves (yarn_curves.h) and moved by the existing knitting solver
// (KnitNumerics::relax, knitting/numerics.h) to the target diameter by radius continuation. One implementation for the
// lab (stitchmesh_lab yarn-settle, yarn-settle-fixtures) and the app's guided Knit route (scene/knit/guided_part.cpp).
//
// The pipeline, each state checked on its own so every source of change can be told apart:
//   raw        the embedding as given, checked as yarn (checkYarn)
//   converted  the yarn as uniform cubic B-splines through equal-arc points, re-embedded and checked again
//   settled    after the solve, re-embedded and checked again
// A construction that fails its check, raw or after the conversion, is not settled: the solver is never asked to repair a
// crossing order. SettleSettings::failedConstruction says what happens instead: the lab refuses it; the app's Final returns
// the converted yarn unsettled (shown and exported as a diagnostic); the lab's explicit diagnostic opt-in settles it as it
// is and marks it, so the construction can no longer pass whatever the motion does. A yarn whose strand or ports are
// broken is refused under every choice.
//
// The solve: production's NCP problem (knit_problem.h: soft edge lengths at the converted control polygon, hard contact
// and curvature rows, a band around the reference, a shape anchor to each call's start, a displacement cap) with its own
// step-wise swept pass-through test. A constrained geometric post-pass, not a physical yarn simulation. The tubes start
// interpenetrating, so the target D is reached in stages D_k = D0 (D / D0)^(k / K), K = ceil(log(D / D0) / log(growth)),
// from the start diameter D0 at which the converted curves pass the certified checks (certifiedRepair's rule,
// knitting/relax/settle.h), plus calls at D until the gate is valid there (SettleStop). After each call the curves are validated at D_k (certified curvature
// and clearance, knitting/checks/checks.h; the local signature does not apply: this yarn has no production stitch
// records) and the call's straight motion is certified at D_{k-1} when both its end states pass there
// (certifyStraightMotion). A certified chain from a verified start would prove that no centreline passed through another,
// except same-yarn pairs closer than (pi/2) delta in arc; no state has had one so far (the start is never verified).
#include <functional>
#include <ostream>
#include <string>
#include <vector>

#include "knitting/checks/checks.h"
#include "knitting/numerics.h"
#include "knitting/relax/knit_problem.h"
#include "model/mesh.h"
#include "model/status.h"
#include "model/yarn_path.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/yarn_curves.h"
#include "stitchmesh/yarn_embed.h"

namespace knit {

// What the settle does with a construction whose crossing check fails (raw, or after the conversion).
enum class FailedConstruction {
    Refuse,           // the lab's default: refused before the solve (SettleEnd::RawConstructionFailed, ConvertedConstructionFailed)
    ReturnUnsettled,  // the app's Final: converted and returned without a relax call (SettleEnd::ConstructionFailedNotSettled)
    SettleAnyway,     // the lab's diagnostic opt-in (--settle-failed-construction): settled as it is and marked (research only)
};
const char* toString(FailedConstruction f);

// The frozen settle settings (tests/data/guided/settle_settings.json; guided_settle_defaults_frozen
// compares the defaults with that file). Production's values where production has one: SettleParams (margins, bend factor, length
// stiffness, band, sampling, iterations), RepairParams (displacement fraction, smallest start diameter), PatternParams
// (clearance).
struct SettleSettings {
    double spacing = 0.2;            // conversion: arc spacing of the interpolation points [mm]
    int densePerSpan = 32;           // conversion: dense Hermite points per span
    double clearance = 0.1;          // c [mm] (PatternParams)
    double contactMargin = 0.05;     // soft contact target D + c + this [mm] (SettleParams)
    double bendFactor = 1.15;        // soft bend target 1 / (bendFactor (r + c / 2)) (SettleParams)
    double lengthStiffness = 0.3;    // (SettleParams)
    double bandHalf = 2.0;           // max |offset| from the reference [mm] (SettleParams)
    int samplesPerSpan = 4;          // (SettleParams)
    int iterations = 120, refresh = 10;   // (SettleParams)
    double displacementFraction = 0.5;    // per call: cap = this x (D_{k-1} + c) (RepairParams)
    double D0Min = 0.05;             // a start diameter below this is not tried [mm] (RepairParams)
    double growth = 1.6;             // stage diameter ratio at most
    int maxRepeats = 60;             // calls at the target diameter at most; the stall rule and the budget end a run before that (the hem needed 21 in the lab and more than 24 in the app)
    int stallCalls = 4;              // target calls in a row with neither fewer open regions nor a smaller shortfall: stalled
    double budgetMinutes = 50;       // no call starts after this much solve time (the one-hour cap of every run)
    int pinnedEnd = 3;               // controls held at each yarn end (the declared yarn ends: position and direction)
    bool band = true;
    // The solver options this settle is made for (knit_problem.h). The caller builds its numerics with them: the app in
    // scene/knit/knit_stage.cpp, the lab unless a flag says otherwise. Live rows: a quarter of the memory and half the
    // time of compact rows (torso_lower, 1,222 stitches: 4.6 against 9.8 GB, 358 against 722 s). Span-maximum rows: the
    // fixed curvature samples let a bend hide between two of them.
    ContactStorage contactStorage = ContactStorage::LiveRows;
    CurvatureRows curvatureRows = CurvatureRows::SamplesAndSpanMaximum;
    bool convertOnly = false;        // stop after the conversion and its checks (the lab's spacing probe, the app's Draft)
    FailedConstruction failedConstruction = FailedConstruction::Refuse;
    // diagnostics only (the lab): every converted control coordinate moved this many ulps up (verdict sensitivity); the
    // first head of this face mirrored through the reference before the solve (a deliberate pass-through)
    int perturbUlps = 0;
    int mirrorHeadOfFace = -1;
};

// Production's certified checks at diameter D: curvature (regular, no fold) and clearance.
struct SettleGate {
    double D = 0;
    Status regular = Status::Unresolved, fold = Status::Unresolved, clearance = Status::Unresolved;
    double rhoLo = 0, minDistance = 0;
    int regions = 0, failRegions = 0, unresolvedRegions = 0;
    ValidationReport report;
    bool valid() const { return regular == Status::Pass && fold == Status::Pass && clearance == Status::Pass; }
    std::string line() const;
};
SettleGate settleGateAt(const CurveSet& cs, double D, double c);

// One state's geometry against the converted curves: displacement at equal parameter, length, edge strain, offsets.
struct SettleDrift {
    double maxControlMove = 0, maxPointMove = 0, p95PointMove = 0;
    double length = 0, lengthRest = 0;
    double strainMax = 0, strainRms = 0;
    double meanAbsDepth = 0, maxAbsDepth = 0, meanDepth = 0;   // offset from the reference [mm]
};
SettleDrift settleDrift(const CurveSet& cur, const CurveSet& rest, const TriMesh& reference);

// One relax call of the continuation.
struct SettleStage {
    int k = 0;
    double Dprev = 0, D = 0, cap = 0;
    double ms = 0;
    std::string solver;   // the solver's call line (SettleHooks::lastSolverCall), or "relax threw: ..."
    SettleGate gate;      // at D
    bool startValid = false, certified = false;
    double delta = 0, localScale = 0;
    std::string certificate;
};

// Why the calls at the target diameter ended.
enum class SettleStop {
    GatePassed,        // the gate of a target call is valid and the construction is kept: done
    ConstructionLost,  // the gate of a target call is valid but a crossing was lost on the way: no further call restores it
    Stalled,           // SettleSettings::stallCalls target calls without progress
    MaxRepeats,        // SettleSettings::maxRepeats target calls ran
    Budget,            // SettleSettings::budgetMinutes of solve time spent
    CallRefused,       // a relax call threw (the device cap, or a device failure)
};
const char* toString(SettleStop s);

// The stop rule after one call at the target diameter D (clearance c). `progress` is its memory across those calls: the
// fewest open regions so far and, among states with that many, the smallest shortfall (the bend radius and the distance
// the gate still lacks [mm]). constructionKept is read only for a valid gate.
struct TargetProgress {
    int openRegions = 1 << 30;
    double shortfall = 1e300;
    int callsWithoutGain = 0;
};
enum class AfterTargetCall { Continue, GatePassed, ConstructionLost, Stalled };
AfterTargetCall afterTargetCall(const SettleGate& gate, double D, double c, bool constructionKept, int stallCalls, TargetProgress& progress);

// How far the pipeline went.
enum class SettleEnd {
    RawConstructionFailed,         // refused (FailedConstruction::Refuse): the raw construction check failed
    YarnDisconnected,              // refused (every choice): the strand or its ports are broken, nothing can settle it
    ConversionRefused,             // refused: a strand cannot be converted (closed)
    ConvertedConstructionFailed,   // refused (FailedConstruction::Refuse): the conversion changed the construction
    ConvertedOnly,                 // SettleSettings::convertOnly: converted and checked, not settled
    ConstructionFailedNotSettled,  // FailedConstruction::ReturnUnsettled: a construction check failed; converted, not settled
    Solved,                        // the schedule ran (see refusedCalls, budgetHit for how far)
};
const char* toString(SettleEnd e);

struct SettleResult {
    SettleEnd end = SettleEnd::RawConstructionFailed;
    std::string verdict;               // one line: "settled: target PASS, construction kept, ..." or "refused (...)"
    YarnCheck raw, converted, settled; // construction checks of the three states
    YarnConversion conversion;         // conversion.curves: the converted curves before the solve (mm, the reference's frame)
    YarnEmbedding convertedEmbedding, settledEmbedding;
    CurveSet curves;                   // the retained curves after the solve (Solved; else the converted ones)
    SettleGate gateRaw, gateFinal;     // at D: before the solve, and of the retained curves
    double D0 = 0;
    bool startVerified = false;
    std::vector<SettleStage> stages;
    int certifiedStages = 0;
    bool chainCertified = false;
    double solveMs = 0;
    SettleDrift drift;
    bool constructionKept = false;
    bool budgetHit = false;
    SettleStop stop = SettleStop::MaxRepeats;   // meaningful for SettleEnd::Solved
    int projectionChanged = 0;         // faces as required in one frame only (reference or own plane)
    int letterSwapped = 0;             // faces not as required in either frame
    int refusedCalls = 0;              // relax calls that threw (the device cap's refusal, or a device failure)
    bool rawConstructionFailed = false;        // the raw construction check failed
    bool convertedConstructionFailed = false;  // the converted construction check failed (run whenever the conversion was)
    bool conversionAddedDefects = false;       // the converted yarn has more wrongly crossed faces or fewer interlocks than the raw one
    bool settledFailedConstruction = false;    // settled although a construction check failed (FailedConstruction::SettleAnyway)
    std::string refusal;               // the first thrown call's message
    bool constructionFailed() const { return rawConstructionFailed || convertedConstructionFailed; }
};

// What the caller lets the settle report through. Both may be empty.
struct SettleHooks {
    std::function<void(int call, int calls, double D)> onCall;   // before each relax call (progress)
    std::function<std::string()> lastSolverCall;                  // the device solver's last call line (gpu::ncpNumericsCalls)
};

// What is settled: embedYarn's yarn on mesh (faces: VirtualShift::corners) and its reference (mm); name labels the log.
struct SettleInput {
    const YarnEmbedding& raw;
    const TypedStitchMesh& mesh;
    const std::vector<std::vector<int>>& faces;
    const TriMesh& reference;
    const YarnSettings& yarn;
    std::string name;
};

// The report gets the lines of the lab's settle_report.txt, the log one line per call. Throws Cancelled when the thread's
// token fires; a relax call that throws anything else (the device budget) ends the schedule and counts in refusedCalls.
SettleResult settleYarn(const SettleInput& in, const SettleSettings& settings, KnitNumerics& numerics, const SettleHooks& hooks, std::ostream& report,
                        std::ostream& log);

} // namespace knit
