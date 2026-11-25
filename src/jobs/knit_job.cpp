// The Knit job: the knit read from the scene cache, else knitted with partial scenes, then displayed and cached.
#include <algorithm>

#include "jobs/display_scene.h"
#include "jobs/job_steps.h"
#include "jobs/scene_cache.h"
#include "model/cancel.h"
#include "model/log.h"
#include "scene/knit/knit_stage.h"
#include "scene/knit/saved_knit.h"
#include "scene/result_keys.h"
#include "scene/source/scene_parts.h"
#include "scene/stage_params.h"
#include "settings/settings_profiles.h"

namespace knit {

// Where a knit job reads and writes the persistent scene cache (scene_cache.h).
struct KnitCacheSlot {
    bool enabled = false;  // a cache root is set and the scene's assets have a hash
    uint64_t asset = 0;    // sceneAssetHash
};

static KnitCacheSlot knitCacheSlot(const KnitRequest& r, const JobControl& control, const SceneEntry& e) {
    KnitCacheSlot cache;
    if (r.cacheRoot.empty()) return cache;
    reportProgress(control, "reading the scene cache", 0, 0);
    cache.asset = sceneAssetHash(e);
    cache.enabled = cache.asset != 0;
    return cache;
}

// A cached tiled knit's computational knit: its knit.bin, restored; without it the knit is
// for viewing only.
static void attachCachedKnit(const KnitRequest& r, const JobControl& control, const SceneEntry& e, uint64_t asset, KnitResult& knit) {
    std::string bytes, why;
    auto t = std::make_shared<Scene>();
    if (!loadSceneCacheKnit(r.cacheRoot, *r.settings, asset, bytes)) {
        knit.computational = false;
        knit.viewOnlyReason = "the cache entry predates stored knits";
        return;
    }
    reportProgress(control, "restoring the knit", 0, 0);
    if (!restoreKnit(bytes, e, *r.settings, *t, &why)) {
        knit.computational = false;
        knit.viewOnlyReason = why;
        return;
    }
    adoptKnit(knit, t);
}

// The persistent cache (scene_cache.h): a knit that was computed before is read, with its
// computational knit (tiled: knit.bin; teddy: the knit cache files) when it has one.
// nullptr on a miss.
static std::shared_ptr<BuiltScene> restoreCachedKnit(const KnitRequest& r, const JobControl& control, const SceneEntry& e, const KnitCacheSlot& cache,
                                              const std::shared_ptr<KnitResult>& knit) {
    const SceneSettings& s = *r.settings;
    std::shared_ptr<BuiltScene> c = loadSceneCache(r.cacheRoot, s, cache.asset);
    if (!c) return nullptr;
    c->kind = BuiltScene::Kind::Knit;
    c->revision = control.revision;
    c->settings = s;
    if (const double mapD = mostCommonD(c->parts); mapD > 0)
        if (std::shared_ptr<const TexturePreview> texture = texturePreview(s, mapD)) c->texture = texture;
    attachCachedKnit(r, control, e, cache.asset, *knit);
    if (knit->scene) {  // a yarn of fibres rendered as maps: the tubes under the cached maps (the stored scene has none)
        reportProgress(control, "showing the yarn maps", 0, 0);
        std::string why;
        if (!showCachedYarnMaps(*c, *knit->scene, s, &why)) c->report += "yarn maps not shown: " + why + "\n";
    }
    summarize(*knit, c->parts);
    knit->id = nextResultId();
    c->knit = knit;
    return c;
}

// The scene's knit: the last one in memory when its Knit-stage key matches, else knitted now, the
// partial scenes showing each part as it finishes. With a frozen knit directory (loadKnitCache) the
// pieces are read from it instead, and that knit is not kept in memory. nullptr on failure.
static std::shared_ptr<Scene> rememberedOrNewKnit(const SceneSettings& s, const JobControl& control, const SceneEntry& e, PartialKnitScenes& partials,
                                                  const std::string& loadKnitCache, std::string* err) {
    const uint64_t key = stageKey(s, Stage::Knit, sceneAssetHash(e));
    const bool frozen = !loadKnitCache.empty();
    if (!frozen)
        if (std::shared_ptr<Scene> remembered = rememberedKnit(key)) return remembered;
    auto t = std::make_shared<Scene>();
    KnitRunOptions ro;
    ro.capS = s.timeCapS;
    ro.loadKnitCache = loadKnitCache;
    reportProgress(control, "loading the scene", 0, 0);
    ro.onLoaded = [&](const Scene& ts) { partials.showStandIns(s, ts); };
    ro.onStep = [&](const std::string& step, int done, int total) { reportProgress(control, "guided stitch mesh: " + step, done, total); };
    ro.onPart = [&](size_t i, PartEvent event, int done, int total) {
        const ScenePart& r = t->parts[i];
        if (event == PartEvent::Started) {
            reportProgress(control, "knitting " + r.name, done, total);
            return;
        }
        reportProgress(control, "knitted " + r.name, done, total);
        if (!control.reporter || !r.knit) return;
        std::optional<DisplayMesh> tube = coarsePartTube(r, s);
        if (!tube) return;
        partials.showPart(i, std::move(*tube));
        if (done == total) partials.publishNow();
        else partials.publishIfDue();
    };
    if (!knitScene(e, s, ro, *t, err)) return nullptr;
    // A part refused at the device memory budget depends on this machine's free memory at
    // the time: such a knit is shown, but neither kept in memory nor cached on disk.
    if (!frozen && !refusedParts(*t)) rememberKnit(key, t);
    return t;
}

// The tiled knit's bytes for the scene cache (knit.bin); a knit with a refused or stopped part
// is machine- or time-dependent, so the cache is turned off for it and the report says so.
static std::string knitForCache(const JobControl& control, const Scene& t, KnitCacheSlot& cache, BuiltScene& out) {
    if (cache.enabled && refusedParts(t)) {
        cache.enabled = false;
        out.report += "scene cache not written: a part was refused at the GPU memory budget or stopped by the time cap (knit again to retry)\n";
    }
    if (!cache.enabled) return "";
    reportProgress(control, "storing the knit", 0, 0);
    return serializeKnit(t);
}

// A new knit result of these settings; Final is decided by the snapshot, also for a knit restored from the cache.
static std::shared_ptr<KnitResult> newKnitResult(const SceneSettings& s, const SceneEntry& e) {
    auto knit = std::make_shared<KnitResult>();
    knit->settings = s;
    knit->final = isFinalProfile(s, e);
    return knit;
}

// The Knit job: the knit restored from the scene cache, else knitted (any layout: knitScene)
// with partial scenes while it runs, then displayed and cached.
std::shared_ptr<BuiltScene> knitObject(const KnitRequest& r, const JobControl& control, std::string* err) {
    const SceneEntry* found = jobScene(*r.settings, BuiltScene::Kind::Knit, err);
    if (!found) return nullptr;
    const SceneEntry& e = *found;
    const auto t0 = std::chrono::steady_clock::now();
    const SceneSettings& s = *r.settings;
    std::shared_ptr<KnitResult> knit = newKnitResult(s, e);
    // A frozen knit (loadKnitCache) neither reads nor writes the scene cache: what is shown is that knit. Nor does a guided
    // knit: a restore is verified by gate (a) with production's local signature, which does not apply to its stitches, so
    // it is knitted again from its inputs (and kept in memory by its key like any knit).
    const bool frozen = !r.loadKnitCache.empty();
    KnitCacheSlot cache = frozen || knitsGuidedPiece(s) ? KnitCacheSlot{} : knitCacheSlot(r, control, e);
    if (cache.enabled)
        if (std::shared_ptr<BuiltScene> c = restoreCachedKnit(r, control, e, cache, knit)) {
            if (knit->computational || !r.knitAgain) {
                c->ms = msSince(t0);
                return c;
            }
            knit = newKnitResult(s, e);  // view only again: knit, and replace the entry
        }
    std::shared_ptr<BuiltScene> out = newResultScene(BuiltScene::Kind::Knit, control, s);
    PartialKnitScenes partials(control, std::make_shared<const BuiltScene>(*out));
    std::shared_ptr<Scene> t = rememberedOrNewKnit(s, control, e, partials, r.loadKnitCache, err);
    if (!t) return nullptr;
    knit->frozenKnitDir = r.loadKnitCache;
    reportProgress(control, "building the display", 0, 0);
    if (!knitDisplay(*out, e, *t, s, err)) return nullptr;
    adoptKnit(*knit, t);
    const std::string knitBytes = knitForCache(control, *t, cache, *out);
    updateSceneBounds(*out);
    summarize(*knit, out->parts);
    knit->id = nextResultId();
    out->ms = msSince(t0);
    if (cache.enabled) {
        reportProgress(control, "writing the scene cache", 0, 0);
        storeInSceneCache(r.cacheRoot, s, cache.asset, *out, &knitBytes);
    }
    out->knit = knit;
    return out;
}

} // namespace knit
