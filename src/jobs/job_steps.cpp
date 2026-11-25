// What the jobs share: which scene a job runs on, reporting, a result's identity, the last knit in memory,
// and the progress of the Knit and Fibres jobs.
#include "jobs/job_steps.h"
#include "jobs/display_scene.h"

#include <algorithm>
#include <atomic>
#include <cstdio>

#include "gpu/knitting.h"
#include "knitting/tube/param_text.h"
#include "model/hash.h"
#include "model/log.h"
#include "jobs/scene_cache.h"
#include "scene/knit/saved_knit.h"
#include "scene/result_keys.h"
#include "scene/source/scene_parts.h"
#include "scene/stage_params.h"
#include "settings/scene_registry.h"

namespace knit {

namespace {

std::atomic<int> g_resultIds{0};

struct RememberedScene {
    uint64_t key = 0;
    std::shared_ptr<Scene> scene;
};
std::mutex g_memoryMutex;
RememberedScene g_lastKnit, g_lastBase;

} // namespace

const ComputeCaps& computeCaps() {
    static const ComputeCaps c = gpu::detectComputeCaps();
    return c;
}

const SceneEntry* jobScene(const SceneSettings& s, BuiltScene::Kind kind, std::string* err) {
    const SceneEntry* e = sceneEntryFor(s);
    if (!e) {
        if (err) *err = s.scene == "import" ? "the imported mesh is not prepared (import it again)" : "unknown scene " + s.scene;
        return nullptr;
    }
    if (kind == BuiltScene::Kind::Knit && furWorkflow(s)) {
        if (err) *err = "the Fur workflow has no knit (workflow = fur)";
        return nullptr;
    }
    return e;
}

int nextResultId() { return ++g_resultIds; }

std::shared_ptr<Scene> rememberedKnit(uint64_t key) {
    std::lock_guard<std::mutex> lk(g_memoryMutex);
    return g_lastKnit.key == key ? g_lastKnit.scene : nullptr;
}

void rememberKnit(uint64_t key, std::shared_ptr<Scene> scene) {
    std::lock_guard<std::mutex> lk(g_memoryMutex);
    g_lastKnit = {key, std::move(scene)};
}

std::shared_ptr<Scene> loadBaseScene(const SceneEntry& e, const SceneSettings& s, std::string* err) {
    const uint64_t key = stageKey(s, Stage::Load, 0);
    {
        std::lock_guard<std::mutex> lk(g_memoryMutex);
        if (g_lastBase.scene && g_lastBase.key == key) return g_lastBase.scene;
    }
    auto t = std::make_shared<Scene>();
    if (!loadBaseParts(e, s, SourceOptions{}, *t, err)) return nullptr;
    std::lock_guard<std::mutex> lk(g_memoryMutex);
    g_lastBase = {key, t};
    return t;
}

void reportProgress(const JobControl& control, const std::string& stage, int done, int total) {
    if (control.reporter) control.reporter->progress(stage, done, total);
}

std::shared_ptr<BuiltScene> newResultScene(BuiltScene::Kind kind, const JobControl& control, const SceneSettings& s) {
    auto out = std::make_shared<BuiltScene>();
    out->kind = kind;
    out->revision = control.revision;
    out->settings = s;
    return out;
}

bool usableKnit(const KnitResult* knit, const SceneSettings& current, const JobControl& control, std::string* err) {
    auto no = [&](const std::string& m) {
        if (err) *err = m;
        return false;
    };
    if (!knit) return no("there is no knitted result: knit first");
    if (!knit->computational) return no("the knitted result was restored for viewing only (" + knit->viewOnlyReason + "): knit again");
    if (knitOutdated(*knit, current)) return no("the knitting settings changed since the knit: update the knit first");
    if (knit->scene) {
        std::lock_guard<std::mutex> lk(knit->m);
        reportProgress(control, "verifying the restored knit (gate a)", 0, 0);
        std::string why;
        if (!verifyKnit(*knit->scene, knit->settings, &why)) return no("the restored knit failed verification (" + why + "): knit again");
    }
    return true;
}

void adoptKnit(KnitResult& knit, const std::shared_ptr<Scene>& t) {
    knit.scene = t;
    for (const ScenePart& r : t->parts) {
        if (r.knit) knit.digest = fnv1a(&knit.digest, sizeof knit.digest, curvesDigest(r.knit->t.curves));
    }
    knit.cross = t->cross.status;
}

void summarize(KnitResult& k, const std::vector<PartRow>& rows) {
    for (const PartRow& r : rows) {
        ++k.parts;
        if (!r.knitted) ++k.plain;
        else if (r.verdict == Status::Pass) ++k.pass;
        else {
            r.verdict == Status::Fail ? ++k.fail : ++k.unresolved;
            k.failing.push_back(r.name + " (" + toString(r.verdict) + ")");
        }
    }
}

bool refusedParts(const Scene& t) {
    for (const ScenePart& r : t.parts)
        if (r.stop != KnitStop::None) return true;
    return false;
}

void storeInSceneCache(const std::string& cacheRoot, const SceneSettings& s, uint64_t asset, BuiltScene& out, const std::string* knitBytes) {
    std::string e2;
    if (carriesYarnMaps(out)) {  // the maps' pixels live in the bake cache: the stored scene has plain tubes
        if (!storeSceneCache(cacheRoot, s, asset, withoutYarnMaps(out), &e2, knitBytes)) out.report += "scene cache not written: " + e2 + "\n";
        return;
    }
    if (!storeSceneCache(cacheRoot, s, asset, out, &e2, knitBytes)) out.report += "scene cache not written: " + e2 + "\n";
}

// --- PartialKnitScenes

PartialKnitScenes::PartialKnitScenes(const JobControl& control, std::shared_ptr<const BuiltScene> head)
    : control_(control), head_(std::move(head)) {}

void PartialKnitScenes::showStandIns(const SceneSettings& s, const Scene& scene) {
    // The stand-ins lie behind the fabric (for the requested D, the largest tried): the knit is centred
    // on the surface, so a finished neighbour's yarn never cuts them. Computed outside the lock.
    KnitStandIns standIns = knitStandIns(scene, fillingOffset(s, s.diameter));
    {
        std::lock_guard<std::mutex> lock(m_);
        standIns_ = std::move(standIns);
    }
    reportProgress(control_, "knitting", 0, static_cast<int>(scene.parts.size()));
    publishNow();
}

void PartialKnitScenes::showPart(size_t part, DisplayMesh tube) {
    std::lock_guard<std::mutex> lk(m_);
    if (done_.size() <= part) done_.resize(part + 1);
    done_[part] = std::make_unique<DisplayMesh>(std::move(tube));
}

void PartialKnitScenes::publishNow() {
    if (!control_.reporter) return;
    auto scene = std::make_shared<BuiltScene>(*head_);
    {
        std::lock_guard<std::mutex> lk(m_);
        last_ = std::chrono::steady_clock::now();
        fillMeshesLocked(*scene);
    }
    send(std::move(scene));
}

void PartialKnitScenes::publishIfDue() {
    if (!control_.reporter) return;
    auto scene = std::make_shared<BuiltScene>(*head_);
    {
        std::lock_guard<std::mutex> lk(m_);
        const auto now = std::chrono::steady_clock::now();
        if (now - last_ < std::chrono::seconds(2)) return;
        last_ = now;
        fillMeshesLocked(*scene);
    }
    send(std::move(scene));
}

void PartialKnitScenes::fillMeshesLocked(BuiltScene& scene) const {
    const std::vector<DisplayMesh>& base = standIns_.base;
    scene.meshes = standIns_.fixed;
    scene.albedo = standIns_.albedo;
    for (size_t i = 0; i < base.size(); ++i) scene.meshes.push_back(i < done_.size() && done_[i] ? *done_[i] : base[i]);
    for (size_t i = base.size(); i < done_.size(); ++i)
        if (done_[i]) scene.meshes.push_back(*done_[i]);
}

void PartialKnitScenes::send(std::shared_ptr<BuiltScene> scene) const {
    updateSceneBounds(*scene);
    control_.reporter->partial(std::move(scene));
}

// --- FibreProgress: the fibre stage's phases as one progress bar that never runs backwards

FibreProgress::FibreProgress(const JobControl& control) : control_(control) {}

void FibreProgress::phase(const std::string& name, int done, int total) {
    std::lock_guard<std::mutex> lk(m_);
    double lo = 0, hi = 0.10;  // preparing
    if (name == "fibres: generating") lo = 0.10, hi = 0.90;
    else if (name == "fibres: finishing parts") lo = 0.90, hi = 1.0;
    const double f = total > 0 ? lo + (hi - lo) * std::clamp(double(done) / total, 0.0, 1.0) : lo;
    // The tiled phases report a part when they start it (done = the parts before it).
    char b[64];
    if (done < total) std::snprintf(b, sizeof b, " (part %d of %d)", done + 1, total);
    else std::snprintf(b, sizeof b, " (%d of %d parts)", done, total);
    publishLocked(name + b, f);
}

void FibreProgress::finishing(const std::string& what) {
    std::lock_guard<std::mutex> lk(m_);
    publishLocked("fibres: " + what, 1.0);
}

void FibreProgress::finish() {
    std::lock_guard<std::mutex> lk(m_);
    published_ = kScale;
    label_ = "fibres: done";
    reportProgress(control_, label_, kScale, kScale);
}

int FibreProgress::published() const {
    std::lock_guard<std::mutex> lk(m_);
    return published_;
}

void FibreProgress::publishLocked(const std::string& label, double fraction) {
    // Work before the result exists ends at kBeforeResult; never backwards.
    const int v = std::clamp(static_cast<int>(fraction * kBeforeResult), 0, kBeforeResult);
    published_ = std::max(published_, v);
    label_ = label;
    reportProgress(control_, label_, published_, kScale);
}

} // namespace knit
