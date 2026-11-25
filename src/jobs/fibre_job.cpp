// The Fibres job: grows the fibres of the reviewed knit, or hair on the meshes, and builds their lines.
#include <algorithm>
#include <cstdio>
#include <map>
#include <optional>

#include "jobs/display_scene.h"
#include "jobs/job_steps.h"
#include "model/cancel.h"
#include "model/log.h"
#include "scene/appearance/stitch_colours.h"
#include "scene/fibres/fibre_stage.h"
#include "settings/settings_mapping.h"

namespace knit {

// --- Fibres job

// How the viewer grows fibres: per-part stage times in the log, without the diagnostics (the fibres
// are identical without them: fibre_options_do_not_change_strands).
static FibreRunOptions viewerFibreOptions() {
    FibreRunOptions fo;
    fo.inner = ExecPolicy{true};
    fo.stageProfile = true;
    return fo;
}

// Adds one part's fibres to the result's counts: all of them, the displayed ones (key below the
// display threshold) and the export's segments.
static void countPartFibres(FibreResult& fr, const PartFibres& pf, float displayBelow) {
    fr.fibres += static_cast<long long>(pf.fs.fibres.size());
    for (const Fibre& f : pf.fs.fibres) fr.displayedFibres += float(f.key) < displayBelow;
    fr.segments += segmentCount(pf.fs);
}

// The fibre segments per object (the segment budget's unit, sceneObjectOf), as grown, kept and
// exported: what the Export step forecasts from. The late taper is the settings'.
static std::vector<ObjectSegments> objectSegmentsOf(const Scene& t, const SceneFibres& tf, const SceneSettings& s) {
    const StrandTaper taper{s.tailTaper, s.tailTaperMinDiameters, s.tailTaperTip};
    std::map<int, ObjectSegments> byObject;
    for (const PartFibres& pf : tf.parts) {
        cancelCheckpoint();
        const ScenePart& part = t.parts[static_cast<size_t>(pf.part)];
        ObjectSegments& o = byObject[sceneObjectOf(t.source, part.prim)];
        const KeptSegments kept = keptSegments(pf, part.centre, taper);
        o.grown += segmentCount(pf.fs);
        o.kept += kept.kept;
        o.exported += kept.exported;
    }
    std::vector<ObjectSegments> out;
    for (const auto& [object, segments] : byObject) out.push_back(segments);
    return out;
}

// One part's display lines: each fibre with its part's colour and its own (linear, as part colours;
// the viewer picks one from the current colour source): a wool fibre its stitch's colour, a hair the
// texture colour at its root. Fur hairs are drawn as one population of their own.
static void addPartFibreLines(DisplayLines& lines, const Scene& t, const PartFibres& pf, const SceneSettings& s,
                       const StitchColours& colours, float displayBelow) {
    const ScenePart& r = t.parts[pf.part];
    const Vec3 pc = objectColourOf(t, r, s);
    const bool fur = furWorkflow(s);
    const bool rootColoured = fur && !pf.fibreColour.empty();
    // Fibres on separate ply tubes grew on the plies: they find their stitch through the ply set.
    const bool stitchColoured = !fur && r.knit && !colours.key.empty();
    std::optional<RenderedPart> plied;
    if (stitchColoured && s.plies > 1) plied = renderedPart(r, s, nullptr);
    std::optional<StitchAlongYarn> along;
    if (stitchColoured && (s.plies <= 1 || plied)) along.emplace(r.knit->t.curves, plied ? plied->plies.get() : nullptr);
    for (const Fibre& f : pf.fs.fibres) {
        if (!(float(f.key) < displayBelow)) continue;  // display decimation (the export keeps them)
        // The fibre's own colour: a hair's root colour, a wool fibre's stitch colour, else its part's.
        Vec3 own = pc;
        if (rootColoured) {
            own = stitchColourLinear(pf.fibreColour[static_cast<size_t>(&f - pf.fs.fibres.data())]);
        } else if (along) {
            const int st = along->at(f.yarnId, f.s);
            own = st >= 0 && st < static_cast<int>(colours.key.size()) ? stitchColourLinear(colours.key[st]) : pc;
        }
        const uint8_t population = fur ? 2 : f.population == FibrePopulation::Random ? 1 : 0;
        addFibreLine(lines, pf.fs.pointsOf(f), r.centre, population, pc, own);
    }
}

// The display lines of a tiled knit's fibres (at most kDisplaySegments segments, the smallest
// keys first), with the result's counts.
static void addFibreLinesOfKnit(BuiltScene& out, FibreResult& fr, const Scene& t, const SceneFibres& tf, const SceneSettings& s) {
    const float qd = displayThreshold(tf);
    const SceneStitchColours& stitchColours = sceneStitchColours(t);
    for (const PartFibres& pf : tf.parts) {
        cancelCheckpoint();
        addPartFibreLines(out.fibres, t, pf, s, stitchColours.parts[pf.part], qd);
        countPartFibres(fr, pf, qd);
    }
}

// The fibre log's closing lines: where the wool fibres were generated, and this job's own stages.
static std::string fibreJobLog(const SceneSettings& s, double sceneS, double displayS) {
    char b[400];
    const ComputeCaps& caps = computeCaps();
    if (furWorkflow(s)) {
        std::snprintf(b, sizeof b, "fur: load the base meshes %.1f s, display lines %.1f s\n", sceneS, displayS);
        return b;
    }
    std::snprintf(b, sizeof b,
                  "fibres: generated on the GPU (device %s)\n"
                  "fibres: verify the knit %.1f s, display lines %.1f s\n",
                  caps.device ? caps.deviceName.c_str() : caps.whyNoDevice.c_str(), sceneS, displayS);
    return b;
}

// The scene a fibre job grows on: fur on the base scene's meshes (nothing knitted), wool on the
// job's knit, which must be usable. Null with the reason in err.
static std::shared_ptr<const Scene> fibreScene(const FibreRequest& r, const JobControl& control, const SceneEntry& e, std::string* err) {
    if (furWorkflow(*r.settings)) {
        reportProgress(control, "loading the base meshes", 0, 0);
        return loadBaseScene(e, *r.settings, err);
    }
    if (!usableKnit(r.knit.get(), *r.settings, control, err)) return nullptr;
    return r.knit->scene;
}

// The Fibres job: the fibres of the settings' profile grown on their scene (growSceneFibres) and
// their display lines; the result is adopted when the job returns.
std::shared_ptr<BuiltScene> growFibres(const FibreRequest& r, const JobControl& control, std::string* err) {
    const SceneEntry* found = jobScene(*r.settings, BuiltScene::Kind::Fibres, err);
    if (!found) return nullptr;
    const SceneEntry& e = *found;
    const auto t0 = std::chrono::steady_clock::now();
    const SceneSettings& s = *r.settings;
    const bool fur = furWorkflow(s);
    const std::shared_ptr<const Scene> scene = fibreScene(r, control, e, err);
    if (!scene) return nullptr;
    const double sceneS = secondsSince(t0);
    // Work done (FibreProgress, job_steps.h).
    FibreProgress prog(control);
    std::shared_ptr<BuiltScene> out = newResultScene(BuiltScene::Kind::Fibres, control, s);
    auto fr = std::make_shared<FibreResult>();
    fr->workflow = workflowOf(s);
    fr->key = fibreResultKeyOf(s);
    fr->settings = s;
    if (fur) {
        fr->base = scene;
        for (const ScenePart& part : scene->parts) fr->baseParts.push_back(part.name);
    } else {
        out->knit = r.knit;
        fr->knitId = r.knit->id;
    }
    reportProgress(control, fur ? "growing fur" : "generating fibres", 0, 0);
    FibreRunOptions fo = viewerFibreOptions();
    fo.onPart = [&prog](const std::string& stage, int done, int total) { prog.phase(stage, done, total); };
    std::optional<SceneFibres> grown = growSceneFibres(e, *scene, s, fo, err);
    if (!grown) return nullptr;
    auto tf = std::make_shared<SceneFibres>(std::move(*grown));
    prog.finishing("building the display lines");
    const auto tDisplay = std::chrono::steady_clock::now();
    // Fur shows its base meshes with the hairs; wool's knit is on display already.
    if (fur && !knitDisplay(*out, e, *scene, s, err)) return nullptr;
    addFibreLinesOfKnit(*out, *fr, *scene, *tf, s);
    fr->objectSegments = objectSegmentsOf(*scene, *tf, s);
    fr->log = tf->log + fibreJobLog(s, sceneS, secondsSince(tDisplay));
    fr->scene = tf;
    arrangeFibreLines(out->fibres);
    fr->id = nextResultId();
    out->fibreResult = fr;
    if (fur) updateSceneBounds(*out);
    else out->lo = {1e300, 1e300, 1e300};  // no bounds of its own: the knit's stay
    out->ms = msSince(t0);
    prog.finish();  // the result is complete: it is adopted when this job returns
    return out;
}

} // namespace knit
