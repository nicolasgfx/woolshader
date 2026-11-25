#pragma once
// KnitResult and FibreResult: what the user reviews, immutable once published.
//
// They compute nothing: each holds a stage's scene data, its counts and the settings it was made with.
// Written by the jobs (knit_job.cpp, fibre_job.cpp); read by the next job and the viewer, which also asks
// here whether a result is still current and whether it may be exported. Downstream stages take these
// objects, never settings, because the solver is not deterministic: knitting the same settings again
// need not give the reviewed result.
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "export/export.h"
#include "model/status.h"
#include "settings/scene_settings.h"

namespace knit {

struct Scene;
struct SceneFibres;

// A completed knit, as the user reviews it: the computational data the fibres and the export
// consume, and its settings snapshot. Immutable once published (the tiled scene's lazy
// verification is guarded by `m`).
struct KnitResult {
    int id = 0;                      // per process, increasing
    uint64_t digest = 0;             // curvesDigest over the knitted parts, in part order
    SceneSettings settings;          // what it was knitted with
    bool final = false;              // isFinalProfile(settings): settle, repair, fallback, gate (d), all of them
    // false: shown from the scene cache without its computational knit (an entry written before
    // knit.bin existed, or one that failed to restore): fibres and the export need a new knit.
    bool computational = true;
    std::string viewOnlyReason;
    std::string frozenKnitDir;      // the knit cache directory its pieces were taken from (KnitRequest::loadKnitCache), else ""
    std::shared_ptr<Scene> scene;   // the knitted scene, any layout (UV charts, spiral pieces)
    int parts = 0, pass = 0, fail = 0, unresolved = 0, plain = 0;
    Status cross = Status::Pass;         // the cross-component clearance (tiled)
    std::vector<std::string> failing;    // knitted parts whose verdict is not PASS
    mutable std::mutex m;                // the tiled scene's lazy verification (fibres, export)
};

// Completed fibres of one knit result, exactly what the export writes.
// Fur (workflow = fur): grown on the scene's base meshes, no knit (knitId 0); `base` is the scene it grew
// on, which the export writes as plain meshes with the pile. Wool: `scene` holds every part's fibres, for
// any knit layout.
// One object's fibre segments (fibres.segment_budget holds per object): every one grown, those below
// the object's threshold, and those as the export writes them at the result's late taper.
// The longest knitted yarn [mm] (the yarn maps' atlas is planned per yarn).
double longestYarnMm(const KnitResult& k);

struct ObjectSegments {
    long long grown = 0, kept = 0, exported = 0;
};

struct FibreResult {
    int id = 0;
    int knitId = 0;
    Workflow workflow = Workflow::Knit;
    uint64_t key = 0;                    // fibreResultKey(settings)
    SceneSettings settings;
    std::shared_ptr<const SceneFibres> scene;  // the fibres (fur: the pile) per part
    std::shared_ptr<const Scene> base;    // fur: the base scene
    std::vector<std::string> baseParts;   // fur: the base scene's part names, in order
    long long fibres = 0, segments = 0;
    long long displayedFibres = 0;       // the viewer draws at most kDisplaySegments (by key); the export all
    std::vector<ObjectSegments> objectSegments;  // per object, the Export step's forecast
    std::string log;
};

// The key a FibreResult of these settings carries (FibreResult::key; scene/result_keys.h).
uint64_t fibreResultKeyOf(const SceneSettings& s);
// The knitting settings changed since the knit was made (its geometry key differs).
bool knitOutdated(const KnitResult& k, const SceneSettings& current);
// The fibre settings changed since the fibres were grown.
bool fibresOutdated(const FibreResult& f, const SceneSettings& current);

// --- may it be exported?

enum class FibreChoice { LeftOut, Included };

// One rule for the Export step and the job, computed from the settings and the reviewed results. It
// cannot fail: every outcome is a verdict. Ready, Diagnostic (written, with the failed checks noted
// in report.json) or Blocked (refused, with the reasons and the next action).
ExportEligibility exportEligibility(const SceneSettings& current, const KnitResult* knit, const FibreResult* fibres,
                                    FibreChoice choice);

} // namespace knit
