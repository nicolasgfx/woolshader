// The strand mesh rebuild: the shown knit's yarn tubes at a new strand tessellation, nothing knitted.
#include <algorithm>

#include "jobs/display_scene.h"
#include "scene/appearance/baked_maps.h"
#include "scene/appearance/rendered_part.h"
#include "scene/appearance/yarn_tube.h"
#include "scene/stage_params.h"
#include "settings/yarn_look.h"
#include "jobs/job_steps.h"
#include "model/cancel.h"
#include "model/log.h"
#include "scene/knit/saved_knit.h"
#include "scene/source/scene_parts.h"
#include "settings/settings_mapping.h"

namespace knit {

// The rebuild needs the shown scene of a usable knit and a valid strand tessellation.
static bool canRebuildStrandMesh(const StrandMeshRequest& r, const JobControl& control, std::string* err) {
    if (!r.shown) {
        if (err) *err = "the strand mesh rebuild needs the shown knit";
        return false;
    }
    if (!usableKnit(r.shown->knit.get(), *r.settings, control, err)) return false;
    if (std::string why; !tubeSettingsValid(*r.settings, &why)) {
        if (err) *err = why;
        return false;
    }
    return true;
}

// The shown scene, stamped with this job; its yarn length and display note are set again by
// the rebuilt tubes.
static std::shared_ptr<BuiltScene> copyShownScene(const StrandMeshRequest& r, const JobControl& control) {
    auto out = std::make_shared<BuiltScene>(*r.shown);
    out->kind = BuiltScene::Kind::Knit;
    out->revision = control.revision;
    out->partial = false;
    out->cached = false;
    out->settings = *r.settings;
    out->yarnLength = 0;
    out->displayNote.clear();
    return out;
}

// Bakes a part's yarn maps into the cache when it lacks them (the Knit step's Bake the yarn maps), a yarn at a
// time and reported; a part whose maps are cached costs a read of their plans.
static bool bakePartYarnMaps(const ScenePart& p, const SceneSettings& s, const JobControl& control, std::string* err) {
    const std::optional<RenderedPart> rp = renderedPart(p, s, err);
    if (!rp) return false;
    BakedMapSpec spec;
    if (!bakedMapSpecOf(s, p.name, *rp->strands, rp->strandPp, spec, err)) return false;
    BakedMapEntry entry;
    if (loadBakedMapSet(bakeCacheRoot(), spec, entry)) return true;
    SurfaceModel surface;
    if (!buildYarnSurface(*rp->strands, rp->strandPp, yarnSurfaceParams(s), surface, err)) return false;
    BakeProgress progress;
    progress.report = [&](int yarnsDone, int yarns) {
        reportProgress(control, "baking " + p.name + ": yarn " + std::to_string(yarnsDone) + " of " + std::to_string(yarns), yarnsDone, yarns);
    };
    return ensureBakedMapSet(bakeCacheRoot(), spec, *rp->strands, surface, BakeInterrupt::Cancellable, progress, entry, err);
}

// The tubes of a tiled knit at the new tessellation, from its retained Scene, by part name in `tubes`, as the
// settings render the yarn; a request to bake first bakes the maps the cache lacks, so the tubes come up under them.
static bool rebuildKnitTubes(const StrandMeshRequest& r, const JobControl& control, const KnitTubeLook& look, BuiltScene& tubes, std::string* err) {
    const SceneSettings& s = *r.settings;
    const Scene& t = *r.shown->knit->scene;
    const bool bake = r.bakeMaps && yarnBodyRenderOf(s) == YarnBodyRender::Maps;
    int done = 0, total = 0;
    for (const ScenePart& p : t.parts) total += p.knit != nullptr;
    for (const ScenePart& p : t.parts) {
        cancelCheckpoint();
        if (!p.knit) continue;
        reportProgress(control, "rebuilding the strand mesh", done++, total);
        if (bake && !bakePartYarnMaps(p, s, control, err)) return false;
        if (!addKnittedPartTubes(tubes, t, p, s, look, err)) return false;
    }
    return true;
}

// The scene cache entry of this tessellation, with the knit's knit.bin; a knit with refused parts
// is never cached.
static void storeRebuiltStrandMesh(const StrandMeshRequest& r, const JobControl& control, const SceneEntry& e, BuiltScene& out) {
    const Scene& t = *r.shown->knit->scene;
    if (r.cacheRoot.empty() || refusedParts(t)) return;
    const uint64_t asset = sceneAssetHash(e);
    if (!asset) return;
    reportProgress(control, "writing the scene cache", 0, 0);
    const std::string bytes = serializeKnit(t);
    storeInSceneCache(r.cacheRoot, *r.settings, asset, out, &bytes);
}

// The strand tessellation changed (strand.tube_*, a Fibres-stage setting: never a knitting input):
// the knit's yarn tubes are rebuilt from its knitted curves (the retained Scene; nothing is
// knitted) and replace the tubes of the shown scene in place; every other mesh, the rows and
// the knit result (its id) stay. The scene cache gets the entry for this tessellation.
std::shared_ptr<BuiltScene> rebuildStrandMesh(const StrandMeshRequest& r, const JobControl& control, std::string* err) {
    const SceneEntry* e = jobScene(*r.settings, BuiltScene::Kind::Knit, err);
    if (!e) return nullptr;
    const auto t0 = std::chrono::steady_clock::now();
    const SceneSettings& s = *r.settings;
    if (!canRebuildStrandMesh(r, control, err)) return nullptr;
    std::shared_ptr<BuiltScene> out = copyShownScene(r, control);
    // The rebuilt tubes by part name, then put where the old ones were.
    const KnitTubeLook look = knitTubeLook(*r.shown->knit->scene, s, mostCommonD(out->parts));
    if (look.texture) out->texture = look.texture;
    out->displayNote = look.displayNote;
    BuiltScene tubes;
    if (!rebuildKnitTubes(r, control, look, tubes, err) || !replaceYarnMeshes(*out, tubes, err)) return nullptr;
    out->yarnLength = tubes.yarnLength;
    if (!tubes.displayNote.empty()) out->displayNote = tubes.displayNote;
    out->ms = msSince(t0);
    storeRebuiltStrandMesh(r, control, *e, *out);
    out->knit = r.shown->knit;  // the same knit result (id unchanged)
    return out;
}

} // namespace knit
