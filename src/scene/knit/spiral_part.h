#pragma once
// One sewn piece knitted as one spiral (the teddy): the piece is charted from pole to
// pole and knitted at the scene's yarn diameters (partDiameters), or taken from the per-piece knit
// cache (spiral_cache.h), then checked by gate (a) and, when requested, by the surface gate (d).
// knitScene knits every spiral part this way, the teddy's pieces (sewn_pieces.h) and an import's closed
// spheres alike: the part carries its piece. The part iteration, time caps, cancellation, progress and
// numerics are the common ones (knit_stage.h); only the piece's own stages are here. Private to the
// Knit stage.
#include <chrono>
#include <string>
#include <vector>

#include "scene/knit/knit_cache.h"
#include "knitting/checks/checks.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/knit/pieces.h"
#include "scene/scene.h"
#include "scene/run_log.h"

namespace knit {

struct KnitRun;

// A piece's log: collected here, printed as one block when the piece finishes, and the knitting's part
// of it stored in the knit cache.
struct PieceLog {
    std::string text;
    bool stageProfile = false;  // --stage-profile
    // "     profile: <stage> <wall> ms (threads busy ...)" after a stage.
    void profile(const StageClock& clock, const std::string& label);
};

// One piece in its stages: its part and log, and the clocks its summary lines read.
struct PieceWork {
    size_t k;              // the part's index
    const ScenePart& part;
    std::string object, material;  // of the part's primitive: the chart line
    std::string stem;      // the part name: the cache file names
    PieceLog& log;
    std::chrono::steady_clock::time_point started = std::chrono::steady_clock::now();  // the chart line's time
    StageClock clock;      // --stage-profile: "load + chart" and "piece total"
    Status gateD = Status::Pass;  // the surface gate (d): the worst of its passes (Pass when not requested)
    std::string stopReason;       // why the piece could not be knitted
};

// One knitting attempt at one yarn diameter, with its checks.
struct KnitAttempt {
    ScenePieceKnit spk;
    ValidationReport r;      // gate (a)
    PatternParams pp;
    GraphCheck gc;           // the stitch program
    LayoutCheck lc;          // the stitch layout
    bool genChecks = true;   // gate (a) includes the generator's checks (not after settle or repair)
};

// The piece's kept knitting result and its record (what the knit cache stores with it).
struct PieceKnitResult {
    KnitAttempt kept;
    std::vector<KnitCacheAttempt> attempts;  // every diameter tried, with its gate (a) verdict
    bool residue = false;  // kept with its failures (no diameter passes, or the budget stopped the search)
    bool cached = false;   // from the knit cache
    double msKnit = 0, msValidate = 0;  // the last knitted attempt's times (0 from the knit cache)
};

// The spiral knitting's parameters: the scene's gauge, the spiral's fixed shaping (piece_knit.h
// defaults, in the knit cache key) and the piece's hand.
PieceKnitParams pieceKnitParams(const KnitRun& run, const PieceWork& piece);

// The piece of a spiral part, as prepareScenePiece and knitScenePiece take it.
ScenePieceRequest pieceRequest(const KnitRun& run, const PieceWork& piece);

} // namespace knit
