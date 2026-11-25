#pragma once
// knitScene: the Knit stage of every scene, whatever its layout.
//
// Loads the parts (scene_parts.h) and knits them one after another: each part is knitted by
// knitPart, gap patches are added as a second pass, and the clearance across parts is measured.
// Called by the Knit job (jobs/scene_builder.cpp) and knit_check. Returns false with the reason for
// an error of the run, and throws Cancelled when the caller's token fires. A part that fails its
// checks is a result, not an error.
#include <functional>
#include <string>
#include <vector>

#include "scene/scene.h"

namespace knit {

// A part starting or finished (KnitRunOptions::onPart).
enum class PartEvent { Started, Finished };

struct KnitRunOptions {
    bool verbose = false;
    SourceOptions source;          // how the scene loads
    // Progress for the caller (the viewer's builder), on the caller's thread, in order. onLoaded each time
    // the part list is set: once the scene is loaded, and again after gap patches add parts; onPart when
    // part `index` starts and when it is done. Neither may keep a reference past the call.
    std::function<void(const Scene&)> onLoaded;
    std::function<void(size_t index, PartEvent event, int done, int total)> onPart;
    // Inside a long part (the guided route's settle calls): the step and how far it is, on the caller's thread.
    std::function<void(const std::string& step, int done, int total)> onStep;
    bool stageProfile = false;     // per-attempt stage times (ScenePart::profile)
    // Inside a part (validation, settle): parallel on the shared executor; identical results.
    ExecPolicy inner{true};
    // Time caps [s] (0: none): the whole knitting run, and each part. A part past its deadline
    // stops at its next cancellation checkpoint (cooperative: a long device call or library
    // operation finishes first; the external capped_run.ps1 is the hard limit), is reported
    // "time cap" and the next parts still run; past the scene deadline the remaining parts are
    // "not started".
    double capS = 0, partCapS = 0;
    std::string only;          // non-empty: knit only the parts whose name contains it (knit_check --part)
    // The spiral layout's per-piece knit cache files (spiral_cache.h, knit_cache.h; knit_check
    // --save-knit / --load-knit): every knitted piece written to, or every piece taken from, the
    // directory (a cache whose inputs differ is refused: the run fails, nothing is knitted again).
    std::string saveKnitCache, loadKnitCache;
};

// Knits every part of scene e: loads the scene by its layout (a GLB or an import on its UV charts; a
// scene of sewn pieces as spirals, sewn_pieces.h), then one loop over the parts for every layout
// (time caps, cancellation, progress, the part filter, the device numerics), each part knitted by
// its layout's operation. False with the reason for an error of the run; throws
// Cancelled when the caller's token fires.
bool knitScene(const SceneEntry& entry, const SceneSettings& settings, const KnitRunOptions& options, Scene& out, std::string* err);

// Cross-component clearance (report only, sampled, does not block export): the
// gate-(a) sample polylines of different knitted parts closer than R_a + R_b + c, R the final
// surface's support radius at each part's D (D + c for round tubes at one D).
bool crossComponentClearance(const std::vector<ScenePart>& parts, const SceneSettings& settings, CrossComponentReport& out, std::string* err);

} // namespace knit
