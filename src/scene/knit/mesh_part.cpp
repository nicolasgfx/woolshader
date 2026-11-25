// One chart part: parameters for each attempt, then knitMesh, trying each diameter.
#include "scene/knit/knit_run.h"

#include <algorithm>
#include <cstdio>

#include "gpu/device.h"
#include "gpu/knitting.h"
#include "knitting/tube/yarn_surface.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"
#include "scene/run_log.h"
#include "scene/scene.h"

namespace knit {

namespace {

// One attempt's knitting parameters at diameter D.
MeshKnitParams attemptParams(const KnitRun& run, const ScenePart& part, double D) {
    const SceneSettings& s = run.settings;
    MeshKnitParams p;
    p.tile.pp.yarnDiameter = D;
    p.tile.pp.clearance = s.clearance;
    p.tile.pp.depthAmplitude = s.depthAmplitude;
    // The pattern's relaxed shape: its stitch (rib narrower, garter shorter than the settings'
    // stockinette stitch) and a rib's layers.
    applyPatternShape(p.tile.pp, s, run.pattern);
    p.tile.g.w = p.tile.pp.waleSpacing;
    p.tile.g.h = p.tile.pp.courseSpacing;
    p.tile.angleDeg = s.angleDeg;
    if (part.chart >= 0 || part.gapPatch) {
        // knit.coverage = overlap: this chart alone, continued beyond it (or a gap patch);
        // no boundary inset (the overhang overlaps the neighbours instead).
        p.tile.onlyChart = part.chart;
        p.tile.overhang = part.gapPatch ? 0 : std::max(0, s.overlapStitches);
        p.tile.footprintCheck = false;
        if (part.gapPatch) p.tile.patch = part.patch;
    }
    p.settle = s.settle;
    p.surface = finalSurface(s, D);
    // Gate (d) does not enter the diameter choice (gate (a) decides): it runs once, on the kept
    // attempt, after the diameter loop.
    p.surfaceGate = false;
    p.surfaceHalf = s.surfaceHalf;
    p.exec = run.options.inner;
    p.tile.solids = run.solids.otherThan(part);
    p.hand = handRequest(s, part.name);
    return p;
}

// --stage-profile: gate (d)'s passes, " [pass <g>: build <ms>, check <ms>]" each.
std::string gateDPassTimes(const PartKnit& k) {
    std::string text;
    char b[160];
    for (size_t g = 0; g < k.msGateDBuild.size(); ++g) {
        std::snprintf(b, sizeof b, " [pass %zu: build %.0f, check %.0f]", g, k.msGateDBuild[g], k.msGateDCheck[g]);
        text += b;
    }
    return text;
}

// --stage-profile: one attempt's stages (tile, gate a, settle with its validate calls, boundary,
// gate d per pass).
std::string attemptProfile(const StageClock& clock, double D, const PartKnit* k) {
    char lb[96];
    std::snprintf(lb, sizeof lb, "D %.2f attempt (%s)", D, k ? toString(k->verdict) : "not knitted");
    std::string line = "    profile: " + clock.line(lb);
    if (k) {
        char b[640];
        const ValidateStats& v = k->vSettle;
        std::snprintf(b, sizeof b, "\n      tile %.0f, gate a %.0f, settle %.0f (%d rounds; validate %lld calls %.0f ms: "
                      "curvature %.0f, sampling %.0f, signature %.0f, clearance %.0f; other %.0f), boundary %.0f, gate d %.0f",
                      k->msTile, k->msValidate, k->msSettle, k->settle.rounds, v.calls, v.ms, v.msCurvature, v.msSampling,
                      v.msSignature, v.msClearance, k->msSettle - v.ms, k->msBoundary, k->msGateD);
        line += b;
        line += gateDPassTimes(*k);
    }
    return line + " ms\n";
}

} // namespace

void knitMeshPart(const KnitRun& run, ScenePart& part) {
    const SceneSettings& s = run.settings;
    const DiameterPlan plan = partDiameters(s);
    const std::vector<double>& dList = plan.attempt;
    const std::vector<double>& notAttempted = plan.notAttempted;
    MeshKnitParams kept;  // the kept attempt's parameters (gate (d) after the loop)
    for (size_t di = 0; di < dList.size(); ++di) {
        MeshKnitParams p = attemptParams(run, part, dList[di]);
        auto k = std::make_unique<PartKnit>();
        std::string why;
        const StageClock attemptClock;
        bool knitOk = false;
        try {
            knitOk = knitMesh(*part.mesh, p, run.numerics, *k, &why);
        } catch (const gpu::DeviceBudgetExceeded& ex) {
            why = ex.what();  // refused at the device budget: the part stays unknitted, with the reason
            part.stop = KnitStop::DeviceBudget;
        }
        if (run.options.stageProfile) part.profile += attemptProfile(attemptClock, dList[di], knitOk ? k.get() : nullptr);
        if (!knitOk) {
            part.note = why;
            part.knit.reset();
            break;
        }
        part.knit = std::move(k);
        part.yarnD = dList[di];
        kept = p;
        if (part.knit->after.verdict == Status::Pass) break;
        if (di + 1 < dList.size()) part.note += "D " + std::to_string(dList[di]).substr(0, 4) + ": gate a " + toString(part.knit->after.verdict) + "; ";
    }
    if (part.knit && part.knit->after.verdict != Status::Pass && !notAttempted.empty()) {
        part.note += "diameter search stopped by budget (knit.diameter_attempts = " + std::to_string(s.diameterAttempts) + "): not attempted D";
        for (double d : notAttempted) part.note += " " + std::to_string(d).substr(0, 4);
        part.note += " mm (not a claim that they would fail); ";
    }
    if (part.knit && (s.surfaceGate || s.surfaceHalf)) {
        kept.surfaceGate = true;
        kept.gateDInVerdict = !yarnFit(s).compressed();  // compressed yarn: gate (d) is a report
        const StageClock gateClock;
        runSurfaceGate(*part.knit, kept, run.numerics);
        if (run.options.stageProfile) {
            part.profile += "    profile: " + gateClock.line("gate d (kept attempt)") + gateDPassTimes(*part.knit) + " ms\n";
        }
    }
}

} // namespace knit
