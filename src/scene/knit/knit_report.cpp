// The text report of a knitted scene, part by part.
#include "scene/knit/knit_report.h"

#include <algorithm>
#include <cstdio>
#include <map>

#include "knitting/knit/mesh_knit.h"
#include "model/log.h"
#include "scene/scene.h"
#include "settings/settings_mapping.h"

namespace knit {

namespace {

// Gate (a)'s lines: all of them, or the ones that did not pass.
void reportGateA(std::string& options, const PartKnit& k, bool verbose) {
    for (const CheckLine& l : k.after.lines)
        if (verbose || l.status != Status::Pass) appendFormat(options, "    gate a [%s] %s: %s\n", toString(l.status), l.name.c_str(), l.detail.c_str());
}

std::string noteInBrackets(const ScenePart& part) { return part.note.empty() ? "" : " (" + part.note + ")"; }

// A spiral piece has no chart, lattice, projection or settle record: its knitting line, gate (a) (the
// CPU audit) and gate (d).
void reportSpiralPart(std::string& options, const ScenePart& part, const SceneSettings& settings, bool verbose) {
    const PartKnit& k = *part.knit;
    appendFormat(options, ": %s | D %.2f%s | %zu stitches (%.0f ms)\n    %s\n", toString(k.verdict), part.yarnD, noteInBrackets(part).c_str(),
                 k.t.prog.stitches.size(), part.ms, k.detail.c_str());
    reportGateA(options, k, verbose);
    options += "    gate d " + std::string(settings.surfaceGate || settings.surfaceHalf ? toString(k.gateDStatus) : "not requested") + "\n";
}

void reportCharts(std::string& options, const PartKnit& k, bool verbose) {
    const UvAtlas* atlas = k.atlas.get();
    for (const ChartTiling& ct : k.t.charts) {
        const UvChart& ch = atlas->charts()[ct.chart];
        if (!verbose && ct.kept == 0 && ch.status == Status::Pass && ch.area < 1) continue;
        appendFormat(options, "    chart %d %s: area %.0f mm2, s1 p50/max %.2f/%.2f, s2 min/p50 %.2f/%.2f, outside [0.75,1.3] %.0f %%, "
                        "up %s %.0f deg | nodes %d -> kept %d (outside %d, ambiguous %d, boundary %d, inside other %d, run ends %d, "
                        "short runs %d) | runs %d | coverage %.2f\n",
                     ct.chart, toString(ch.status), ch.area, ch.distortion.s1[2], ch.distortion.s1[4], ch.distortion.s2[0],
                     ch.distortion.s2[2], 100 * ch.distortion.outsideFraction, ct.upFromWorld ? "world" : "uv-v",
                     ct.baseAngleDeg, ct.candidates, ct.kept, ct.outside, ct.ambiguous, ct.boundary, ct.insideOther, ct.runEnd,
                     ct.shortRun, ct.runs, ct.coverage);
    }
}

void reportMeshPart(std::string& options, const ScenePart& part, const SceneSettings& settings, bool verbose) {
    const PartKnit& k = *part.knit;
    appendFormat(options, ": %s | D %.2f%s | %zu stitches, %zu runs, %d parentless, %d uncaught (%.0f ms)\n", toString(k.verdict), part.yarnD,
                 noteInBrackets(part).c_str(), k.t.prog.stitches.size(), k.t.prog.turns.size(), k.t.lattice.parentless, k.t.lattice.uncaught, part.ms);
    reportCharts(options, k, verbose);
    int projected[5] = {0, 0, 0, 0, 0};
    for (const auto& yarn : k.t.bind.result)
        for (ProjResult result : yarn) ++projected[static_cast<int>(result)];
    appendFormat(options, "    projection: %d owned, %d budget, %d normal-turn, %d open-boundary, %d depth-band\n", projected[0], projected[1], projected[2],
                 projected[3], projected[4]);
    options += "    lattice " + std::string(toString(k.t.lattice.status)) + ": " + k.t.lattice.detail + "\n";
    options += "    layout " + std::string(toString(k.t.layout.status)) + ": " + k.t.layout.detail + "\n";
    reportGateA(options, k, verbose);
    appendFormat(options, "    gate a before settle %s, after %s (settle %d rounds, max move %.3f mm); ms: tile %.0f, validate %.0f, settle %.0f, "
                    "boundary %.0f, gate d %.0f\n",
                 toString(k.before.verdict), toString(k.after.verdict), k.settle.rounds, k.settle.maxMove, k.msTile, k.msValidate, k.msSettle,
                 k.msBoundary, k.msGateD);
    options += "    " + k.boundary.detail + "\n";
    for (size_t g = 0; g < k.gateD.size(); ++g) {
        const SurfaceCheck& sc = k.gateD[g];
        appendFormat(options, "    gate d (rendered surface, %s rings): folds %s (%d rings, min rho / support %.3f), collisions %s (%ld of %ld "
                        "points penetrating, max depth %.4f mm)%s\n",
                     g == 0 ? "0.025 mm" : "0.0125 mm", toString(sc.fold), sc.foldRings, sc.minFoldRatio < 1e299 ? sc.minFoldRatio : 0.0,
                     toString(sc.collisions), sc.penetratingPoints, sc.pointsTested, sc.maxDepth,
                     yarnFit(settings).compressed() ? "; report only (compressed yarn, not in the verdict)" : "");
    }
}

// The scene's closing lines: the counts, and the clearance across parts.
void reportTotals(std::string& options, const Scene& scene, const SceneSettings& settings) {
    int knitted = 0, pass = 0, fail = 0, unresolved = 0;
    long stitches = 0;
    for (const ScenePart& part : scene.parts) {
        if (!part.knit || part.stop == KnitStop::TimeCap || part.stop == KnitStop::NotStarted) continue;
        ++knitted;
        stitches += static_cast<long>(part.knit->t.prog.stitches.size());
        pass += part.knit->verdict == Status::Pass;
        fail += part.knit->verdict == Status::Fail;
        unresolved += part.knit->verdict == Status::Unresolved;
    }
    appendFormat(options, "parts: %zu, knitted %d (%d PASS, %d FAIL, %d UNRESOLVED), %ld stitches; pattern %s, gauge %.2f x %.2f mm, angle %.0f deg\n",
                 scene.parts.size(), knitted, pass, fail, unresolved, stitches, settings.pattern.c_str(), settings.wale, settings.course, settings.angleDeg);
    if (!hasChartParts(scene)) {
        for (const std::string& l : scene.cross.lines) options += "cross-component clearance: " + l + "\n";
        return;
    }
    char minDistance[48] = "n/a (no other component within reach)";
    if (scene.cross.minDistance < 1e299) std::snprintf(minDistance, sizeof minDistance, "%.3f mm", scene.cross.minDistance);
    appendFormat(options, "cross-component clearance (report only, sampled polylines) %s: %d component pairs, %ld close segment pairs, min distance %s\n",
                 toString(scene.cross.status), scene.cross.pairs, scene.cross.closeSegments, minDistance);
    for (const std::string& l : scene.cross.lines) options += "    " + l + "\n";
}

// A guided part: its outcome first (knitted or not), then gate (a) when it has a knit.
void reportGuidedPart(std::string& options, const ScenePart& part, bool verbose) {
    options += part.knit ? ": " + std::string(toString(part.knit->verdict)) + "\n" : ": NOT KNITTED\n";
    for (const std::string& l : guidedPartLines(*part.guided)) options += "    " + l + "\n";
    if (part.knit) reportGateA(options, *part.knit, verbose);
    if (verbose) options += part.guided->report;
}

} // namespace

std::vector<std::string> guidedPartLines(const GuidedPart& g) {
    std::vector<std::string> lines;
    char b[512];
    lines.push_back(std::string("guided stitch mesh (experimental): ") + toString(g.outcome) + ": " + g.reason);
    std::snprintf(b, sizeof b, "stitch mesh: %s route, %s reference, %d faces, %d stitches (%d INC, %d DEC); its checks %s; digest %016llx", g.route.c_str(),
                  g.preparation.c_str(), g.faces, g.stitches, g.increases, g.decreases, toString(g.meshChecks), static_cast<unsigned long long>(g.stitchMeshDigest));
    lines.push_back(b);
    if (g.calls > 0 || g.outcome == GuidedOutcome::Settled || g.outcome == GuidedOutcome::SettledMeshFails || g.outcome == GuidedOutcome::NotSettled) {
        std::snprintf(b, sizeof b, "settle: D %.3f of %.3f mm reached in %d calls (%d with a certified motion); start D0 %.4f mm %s; %s", g.reachedD, g.targetD, g.calls,
                      g.callsCertified, g.startD, g.startVerified ? "verified" : "not verified (no certified chain from the start)", g.settleVerdict.c_str());
        lines.push_back(b);
        if (!g.gateLine.empty()) lines.push_back("target gate: " + g.gateLine);
        std::snprintf(b, sizeof b, "construction after the motion: %d of %d faces as required in both frames, %d in one frame only, %d in neither; %d of %d interlocks",
                      g.facesChecked - g.facesFailBothFrames - g.facesFailOneFrame, g.facesChecked, g.facesFailOneFrame, g.facesFailBothFrames, g.interlocks,
                      g.interlocksExpected);
        lines.push_back(b);
        std::snprintf(b, sizeof b, "yarn generated -> retained: length %.0f -> %.0f mm, sampled yarn closer than D %.1f %% -> %.1f %%, bends tighter than D/2 %d -> %d samples",
                      g.generatedLengthMm, g.retainedLengthMm, 100 * g.generatedCloserThanD, 100 * g.retainedCloserThanD, g.generatedTightBends, g.retainedTightBends);
        lines.push_back(b);
    }
    if (g.outcome == GuidedOutcome::Drafted || g.outcome == GuidedOutcome::SettleSkipped) lines.push_back("settle: not run, 0 solver calls (" + g.settleVerdict + ")");
    if (!g.constructionBefore.empty()) lines.push_back("construction before any settle: " + g.constructionBefore);
    if (g.curvesDigest) {
        std::snprintf(b, sizeof b, "retained curves: digest %016llx", static_cast<unsigned long long>(g.curvesDigest));
        lines.push_back(b);
    }
    for (const std::string& f : g.findings) lines.push_back("finding: " + f);
    if (!g.yarnShown.empty()) lines.push_back("shown: " + g.yarnShown);
    if (!g.yarnExported.empty()) lines.push_back("exported: " + g.yarnExported);
    std::snprintf(b, sizeof b, "time: stitch mesh %.1f s, yarn %.1f s, %s %.1f s", g.msMesh / 1000, g.msYarn / 1000, g.calls > 0 ? "settle" : "conversion and checks",
                  g.msSettle / 1000);
    lines.push_back(b);
    lines.push_back("limits: one settle run, not an ensemble (verdicts at this setting can change under rounding); the local signature does not apply to these "
                    "templates; no motion is certified from the start");
    lines.push_back("inputs: " + g.inputs);
    return lines;
}

std::string formatKnitReport(const Scene& scene, const SceneSettings& settings, bool verbose) {
    std::string options;
    for (const ScenePart& part : scene.parts) {
        appendFormat(options, "[%s] area %.0f mm2, %zu tris", part.name.c_str(), part.mesh->area(), part.mesh->tri.size());
        if (part.layout == KnitLayout::Guided && part.guided) reportGuidedPart(options, part, verbose);
        else if (part.stop == KnitStop::TimeCap || part.stop == KnitStop::NotStarted) options += std::string(part.stop == KnitStop::TimeCap ? ": UNRESOLVED (" : ": NOT STARTED (") + part.note + ")\n";
        else if (!part.knit) options += ": NOT KNITTED (" + part.note + ") -> left out of the export\n";
        else if (part.layout == KnitLayout::Spiral) reportSpiralPart(options, part, settings, verbose);
        else reportMeshPart(options, part, settings, verbose);
        options += part.profile;
    }
    reportTotals(options, scene, settings);
    return options;
}

void printDefects(const CurveSet& cs, const ValidationReport& part, const PatternParams& pp) {
    // Spans bending tighter than r + c/2, by op and template position (span j = j-th
    // control interval of the stitch).
    const double Rc = pp.radius() + 0.5 * pp.clearance;
    std::map<std::string, int> tight;
    for (const CurvatureResult& c : part.curvature)
        for (size_t k = 0; k < c.spanKappaHi.size(); ++k) {
            if (c.spanKappaHi[k] * Rc <= 1) continue;
            int scene = cs.topo.stitchAt(0, k + 0.5);
            char key[96];
            if (scene < 0) std::snprintf(key, sizeof key, "padding");
            else {
                const StitchRecord& sr = cs.topo.stitches[scene];
                const char* role = sr.nChildren == 2 ? "INC-parent" : sr.nChildren == 0 ? "no-child"
                                   : cs.topo.stitches[sr.children[0]].op == StitchOp::Dec ? "DEC-parent" : "K-parent";
                std::snprintf(key, sizeof key, "%s/%s span %d%s", toString(sr.op), role,
                              static_cast<int>(k - sr.uBegin), c.spanKappaHi[k] * pp.radius() > 1 ? " FOLD" : "");
            }
            if (++tight[key] == 1 && c.spanKappaHi[k] * pp.radius() > 1)
                std::printf("    first %s at stitch %d (rho_lo %.3f)\n", key, scene, 1.0 / c.spanKappaHi[k]);
        }
    for (const auto& [k, n] : tight) std::printf("    tight: %-40s %d\n", k.c_str(), n);
    // Clearance FAIL regions by the stitches involved.
    std::map<std::string, int> fails;
    for (const ClearanceRegion& g : part.clearance.regions) {
        if (g.status != Status::Fail) continue;
        const ClearanceIssue& is = part.clearance.issues[g.worstIssue];
        int ta = cs.topo.stitchAt(is.yarnA, is.uA), tb = cs.topo.stitchAt(is.yarnB, is.uB);
        auto desc = [&](int scene, double u) {
            if (scene < 0) return std::string("pad");
            const StitchRecord& sr = cs.topo.stitches[scene];
            int span = static_cast<int>(u - sr.uBegin);
            return std::string(toString(sr.op)) + (sr.nChildren == 2 ? "/INCp" : "") + ":" + std::to_string(span);
        };
        int dt = ta >= 0 && tb >= 0 ? cs.topo.stitches[tb].turn - cs.topo.stitches[ta].turn : 99;
        std::string key = desc(ta, is.uA) + " vs " + desc(tb, is.uB) + " dturn " + std::to_string(dt);
        if (++fails[key] == 1)
            std::printf("    first clearance FAIL %s: stitches %d/%d d %.3f\n", key.c_str(), ta, tb, is.distance);
    }
    for (const auto& [k, n] : fails) std::printf("    clearance FAIL: %-36s %d\n", k.c_str(), n);
}

} // namespace knit
