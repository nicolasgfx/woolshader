#pragma once
// A scene on its way through the pipeline: its source, its parts and what has been made of them, in
// the same form for every layout. Plain records and a few questions about them (hasStitches,
// partIncomplete), none of which can fail.
//
// Filled by loadPartsToKnit and knitScene; read by the fibre stage, the export, the caches and the
// viewer.
#include <memory>
#include <string>
#include <vector>

#include "io/glb_load.h"
#include "knitting/knit/mesh_knit.h"
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"

namespace knit {

// UV-chart lattices, or one spiral from pole to pole (a sewn piece of the teddy). Each part has its own; the
// knitted result has the same form either way. Whether the scene knits sewn pieces is a question of the
// settings (knitsSewnPieces: settings/scene_registry.h).
// Guided: one imported piece knitted as a stockinette stitch mesh from the Mark step's markings (knit.route = guided,
// scene/knit/guided_part.cpp); the file's other pieces stay plain meshes.
enum class KnitLayout { UvTiling, Spiral, Guided };

// The guided route's pipeline (templates, yarn settings, the stitch-mesh build options, the settle settings): a stage
// input (derived.guided_version), so a change of any of them makes every guided knit stale. Bump with any of them.
// 2: the trial's layout fixes (knitted-over bumps, cups, the far finish), Draft without a solve, and Final skipping the
// settle of a construction that fails its check.
// 3: knots placed on the face, the relaxation holding the yarn ends only, the
// binding thread's end, the settle's stop rule, live contact rows and span-maximum curvature rows.
// 4: the border rings placed on their boundary loops.
// 5: rings turned for upright wales, border rings placed first, the
// relaxation's border spring, corner and step safeguards.
// 6: a plain bind-off strip, the relaxation's fixed time step.
constexpr int kGuidedVersion = 6;

// How a guided part ended. Settled, SettledMeshFails, NotSettled, Drafted and SettleSkipped carry a knit (its curves: the
// ones the settle retained, or the generated yarn as curves when nothing settled it); Stopped and Refused carry none (the
// yarn, when one was made, is only a labelled preview).
enum class GuidedOutcome {
    Settled,     // the target gate (curvature, clearance at D + c) passes, the construction is kept in both frames, the mesh checks pass
    SettledMeshFails,  // the yarn settled as Settled, but the stitch mesh's own checks (connectivity, geometry, coverage) do not pass
    NotSettled,  // the settle reached the target diameter, but the target gate or the construction fails or stays unresolved
    Drafted,     // Draft (settle off): the generated yarn as curves, no solve; only the checks before a settle ran
    SettleSkipped,  // Final: the yarn's crossing order failed its check (raw or converted), so it was not settled (no solve)
    Stopped,     // the solver refused a call (device memory budget) or ran out of time: not knitted (KnitStop::DeviceBudget)
    Refused,     // the markings, the stitch mesh or the yarn construction cannot be used: not knitted, the findings say why
};
const char* toString(GuidedOutcome o);
// The outcome of a guided part the solver settled, from its three checks. The yarn counts as settled when the target gate
// passes and its construction is kept; the stitch mesh's checks describe the mesh, not the yarn, so a failing mesh check
// alone gives SettledMeshFails, never NotSettled. (The part's verdict combines all three either way.)
GuidedOutcome outcomeAfterSettle(Status targetGate, Status construction, Status meshChecks);

// What a guided part was made from and how each stage ended (the Knit step's details, the export's report).
struct GuidedPart {
    std::string inputs;                // the identity: source digest, piece, its triangles, markings digest, preparation, gauge, D, c, version
    GuidedOutcome outcome = GuidedOutcome::Refused;
    std::string reason;                // one sentence for people
    std::vector<std::string> findings; // the stitch-mesh review's findings, and the checks that did not pass
    std::string route;                 // band (start curve) or pole (start point)
    std::string preparation;           // facets or loop_limit (the reference the layout used)
    int faces = 0, stitches = 0, increases = 0, decreases = 0;
    Status meshChecks = Status::Unresolved;  // the stitch mesh: connectivity, geometry, coverage of the knit region
    int facesChecked = 0, facesFailBothFrames = 0, facesFailOneFrame = 0, interlocks = 0, interlocksExpected = 0;
    double targetD = 0, reachedD = 0, startD = 0;
    bool startVerified = false;
    int calls = 0, callsCertified = 0;
    std::string settleVerdict;         // the settle's own one-line verdict (stitchmesh/yarn_settle.h)
    std::string gateLine;              // the target gate in numbers
    double generatedLengthMm = 0, retainedLengthMm = 0;
    double generatedCloserThanD = 0, retainedCloserThanD = 0;   // fraction of sampled yarn closer than D to other yarn
    int generatedTightBends = 0, retainedTightBends = 0;        // samples bent tighter than D / 2
    uint64_t stitchMeshDigest = 0;      // typedDigest of the relaxed stitch mesh
    uint64_t curvesDigest = 0;          // curvesDigest of the retained curves (part frame): two runs of the same inputs compare equal
    std::string constructionBefore;    // the construction checks before any settle: the generated yarn, then as curves
    std::string yarnShown, yarnExported;  // which geometry the viewport shows and the export writes (the report's words)
    // The generated yarn's tubes (part-local mm): shown instead of a knit when the part has none, or instead of the knit's
    // converted curves when the conversion added defects (the generated yarn is then the faithful one).
    std::shared_ptr<const TriMesh> previewTubes;
    double msMesh = 0, msYarn = 0, msSettle = 0;
    std::string report;                // the settle's report text
};

// Why a part that should have a knit has none. None: it has one, or it is too small by design. A
// stopped part's result is never cached.
enum class KnitStop {
    None,
    DeviceBudget,  // refused: a solve exceeded the device memory (machine-dependent)
    TimeCap,       // stopped at its deadline (UNRESOLVED "time cap")
    NotStarted,    // the scene's deadline passed first
};

struct ScenePart {
    int prim = -1, component = -1;
    std::string name;           // "<mesh>_<material>_c<k>"
    Vec3 centre;                // library mm (the component was re-centred here)
    // Re-centred. On the heap and immutable: the knit's surface and atlas refer to it, so it keeps its
    // address whenever the part moves (parts grow in the gap-patch pass); chart parts share it.
    std::shared_ptr<const TriMesh> mesh;
    std::unique_ptr<PartKnit> knit;  // null: not knitted (too small, error)
    double yarnD = 0;           // after the fallback chain
    std::string note;           // why not knitted, or the fallback trail
    KnitStop stop = KnitStop::None;
    double ms = 0;
    std::string profile;        // --stage-profile: one line per D attempt
    // knit.coverage = overlap (scene_settings.h): the chart this part tiles (-1: every chart of the
    // component), or a gap patch walked from `patch` (TileParams::patch); componentName: the
    // component's part name in the per-component mode (the parts of one component share it).
    int chart = -1;
    bool gapPatch = false;
    TileParams::WalkPatch patch;
    std::string componentName;
    KnitLayout layout = KnitLayout::UvTiling;
    // Spiral: `mesh` holds the chart as uv (a, b) (chart.h); rings is its domain; the poles are tested for
    // lying inside the other spiral parts of the same object.
    int rings = 0, object = 0;
    std::shared_ptr<const GuidedPart> guided;  // Guided: how it was made and how it ended (null otherwise)

    bool hasStitches() const;
};

struct CrossComponentReport {
    Status status = Status::Pass;  // sampled: polyline segments of different components
    int pairs = 0;                 // component pairs with a sample closer than D + c
    long closeSegments = 0;
    double minDistance = 1e300;
    std::vector<std::string> lines;
};

struct SceneStitchColours;
class SceneSolids;

struct Scene {
    SourceScene source;
    std::vector<ScenePart> parts;
    CrossComponentReport cross;
    // Restored from the tiled knit cache (saved_knit.h): gate (a) is rerun on the stored
    // curves (verifyKnit) before fibres or an export use it.
    bool restored = false, verified = true;
    std::vector<uint64_t> afterDigests;  // per part: the stored gate (a) report digest
    // The colour of every stitch, made once on first use (sceneStitchColours): colourwork is an
    // Export setting, so it can be switched on after the knit, and every reader shares this one.
    mutable std::shared_ptr<const SceneStitchColours> stitchColours;
    // The closed components' inside tests (scene_solids.h) that Knit, Fibres and fur read: made from
    // `source` when it loads (loadMeshScene, loadPiecesScene), before any part, and never changed after,
    // because `source` never changes. Every Scene a stage reads comes from one of those loads.
    std::shared_ptr<const SceneSolids> solids;
};

// The name each plain primitive of a scene is written under: "<mesh>_<material>", numbered for
// instances. Empty for a primitive that is not plain.
std::vector<std::string> plainPrimitiveNames(const SourceScene& source);

// The object a primitive belongs to: its source node (a bear of the teddy, an import's one object),
// else -1 - prim. The export groups an object's parts under it, and the fibres' segment budget
// holds per object.
int sceneObjectOf(const SourceScene& source, int prim);
// The name of that object in a log: its source node, else "object <n>" (sceneObjectOf).
std::string sceneObjectName(const SourceScene& source, int prim);

// Whether a part selected by the part filter `only` (empty: every part) is missing its knit: stopped
// at a time cap, never started, refused at the device memory budget, or a spiral part that was not
// knitted (every piece must be), or a guided part that was not knitted (stopped or refused). An unknitted part on UV charts is otherwise so by design (too small
// for the gauge). A knitted part whose checks FAIL is complete: its verdict says so.
bool partIncomplete(const ScenePart& r, const std::string& only);

// A part knitted on UV charts: gap patches and the clearance between parts are theirs.
bool hasChartParts(const Scene& t);

} // namespace knit
