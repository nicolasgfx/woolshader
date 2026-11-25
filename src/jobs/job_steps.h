#pragma once
// What the jobs of this folder share: which scene a job runs on, reporting, a result's identity and counts,
// the last knit in memory, and the progress helpers of the Knit and Fibres jobs. Every body is in
// job_steps.cpp. Private to the folder: the viewer and knit_check include jobs/jobs.h and jobs/scene_builder.h.
#include <chrono>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "jobs/display_scene.h"
#include "jobs/jobs.h"

namespace knit {

// The scene entry a job runs on, or nullptr with the reason: an unprepared import, a Knit job with knitting
// skipped, a later job on conflicting steps.
const SceneEntry* jobScene(const SceneSettings& s, BuiltScene::Kind kind, std::string* err);

// The running job's progress line (the viewer's status bar); silent without a reporter.
void reportProgress(const JobControl& control, const std::string& stage, int done, int total);

// A job's result scene, stamped with the job's kind, revision and settings snapshot.
std::shared_ptr<BuiltScene> newResultScene(BuiltScene::Kind kind, const JobControl& control, const SceneSettings& s);

// The knit result a downstream job consumes: complete, computational, current, verified.
bool usableKnit(const KnitResult* knit, const SceneSettings& current, const JobControl& control, std::string* err);

// The tiled scene becomes the knit result's computational knit, with the digest of its curves
// and the clearance between its parts (a computed knit and a restored one alike).
void adoptKnit(KnitResult& knit, const std::shared_ptr<Scene>& t);

// Counts and the failing parts of a knit result, from its part rows.
void summarize(KnitResult& k, const std::vector<PartRow>& rows);

// A part of the knit was refused at the device memory budget or stopped by the time cap (KnitStop): such a
// knit depends on this machine and moment, so it is neither kept in memory nor cached.
bool refusedParts(const Scene& t);

// Stores the scene in the persistent cache (scene_cache.h). A failed write is noted in the
// scene's report; the job's result stands.
void storeInSceneCache(const std::string& cacheRoot, const SceneSettings& s, uint64_t asset, BuiltScene& out, const std::string* knitBytes);

// A new result's id, increasing over the process.
int nextResultId();
// The last knit in memory, by its Knit-stage key (null: another key): a job that changes only later
// settings (plies, materials) shares the same knit, so nothing is knitted again.
std::shared_ptr<Scene> rememberedKnit(uint64_t key);
void rememberKnit(uint64_t key, std::shared_ptr<Scene> scene);
// The scene's parts without a knit (fur grows on them): cached by the Load stage key.
std::shared_ptr<Scene> loadBaseScene(const SceneEntry& e, const SceneSettings& s, std::string* err);

// The partial scenes of a running Knit job: first the stand-ins (each part's source surface behind the fabric),
// then each part's coarse tube as it finishes. Without a reporter nothing is copied or published.
class PartialKnitScenes {
public:
    PartialKnitScenes(const JobControl& control, std::shared_ptr<const BuiltScene> head);
    // knitScene's onLoaded: this part list's stand-ins replace the previous ones (finished parts stay); reports
    // "knitting" and publishes at once.
    void showStandIns(const SceneSettings& s, const Scene& scene);
    // A finished part's coarse tube (called from the knit's part threads).
    void showPart(size_t part, DisplayMesh tube);
    void publishNow();
    void publishIfDue();  // at most every 2 s

private:
    void fillMeshesLocked(BuiltScene& scene) const;  // the stand-ins, each replaced by its part's tube when done
    void send(std::shared_ptr<BuiltScene> scene) const;
    const JobControl& control_;
    std::shared_ptr<const BuiltScene> head_;  // the result's revision and settings
    std::mutex m_;
    KnitStandIns standIns_;                   // plain primitives, per part its stand-in, their textures
    std::vector<std::unique_ptr<DisplayMesh>> done_;
    std::chrono::steady_clock::time_point last_{};  // the clock's epoch: the first publishIfDue publishes
};

// Work progress of a Fibres job: how much of the job's work is done, not a time estimate. The generation
// reports three phases over all knitted parts (prepare, generate, finish), whatever the knit's layout. Work
// before the result exists (the display lines) stays below 1: only finish() publishes 1, when the result is
// ready to be adopted. One mutex serializes each update with its report, so the published value never
// decreases and an older value is never reported after a newer one.
class FibreProgress {
public:
    static constexpr int kScale = 1000;        // reported as done / kScale
    static constexpr int kBeforeResult = 980;  // the cap until finish()
    explicit FibreProgress(const JobControl& control);

    // Phase "fibres: preparing parts" / "fibres: generating" / "fibres: finishing parts" with done
    // of total parts.
    void phase(const std::string& name, int done, int total);
    // The job's own final stage (display lines): the value stays at its cap.
    void finishing(const std::string& what);
    // The result is ready: reports kScale of kScale.
    void finish();

    int published() const;  // the last reported value (tests)

private:
    void publishLocked(const std::string& label, double fraction);

    mutable std::mutex m_;
    const JobControl& control_;
    int published_ = 0;
    std::string label_;
};

// The Render preview's package (render_preview.cpp; defined in export_job.cpp): the Export job's writer and its
// checks of what it writes from, without the deliverable rule (a draft knit renders too); a package, marked
// diagnostic.
// `progress` (may be empty) gets the stage texts, a bake's yarns counted.
bool exportPreviewPackage(const ExportRequest& r, const SceneEntry& e, const std::function<void(const std::string&)>& progress, std::string* err);

} // namespace knit
