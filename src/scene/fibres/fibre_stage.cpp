// The Fibres stage: prepare each part, generate in budget order, finish with the checks.
#include "scene/fibres/fibre_stage.h"

#include <chrono>
#include <cstdio>
#include <map>
#include <memory>
#include <stdexcept>

#include "fibres/fibre_budget.h"
#include "model/cancel.h"
#include "model/log.h"
#include "model/parallel.h"
#include "scene/fibres/fur_stage.h"
#include "scene/fibres/part_fibre_job.h"
#include "scene/run_log.h"
#include "scene/source/scene_solids.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"

namespace knit {

namespace {

// What the stages of one wool run share.
struct WoolRun {
    const Scene& scene;
    const SceneSettings& settings;
    const FibreRunOptions& options;
    SceneFibres& out;
    int knitted = 0;  // parts with stitches
    std::vector<std::unique_ptr<PartFibreJob>> jobs;  // knitted parts, in part order
    std::vector<FibreSet> sets;                       // one per job
    // Budget order, per job: what its object's scheduler run covered and that run's threshold.
    std::vector<BudgetPartResult> covered;
    std::vector<float> q;

    void progress(const char* stage, int done) const {
        if (options.onPart) options.onPart(stage, done, knitted);
    }
};

// Every knitted part's strands, support and frame. Above the bound on support residency (budget mode,
// FibreRunOptions::fibreSupportMB) a part's support waits released; evaluateWindow rebuilds it when a
// window needs the part.
bool prepareParts(WoolRun& run, const FibreParams& base, std::string* err) {
    const SceneSolids& solids = *run.scene.solids;
    const double modelRing = exportModelRing(run.settings);
    for (const ScenePart& part : run.scene.parts) {
        cancelCheckpoint();
        if (!part.hasStitches()) continue;
        run.progress("fibres: preparing parts", static_cast<int>(run.jobs.size()));
        const PartFibreSetup setup{part, static_cast<int>(&part - &run.scene.parts[0]), run.settings, base, solids, modelRing, run.options};
        std::unique_ptr<PartFibreJob> job = PartFibreJob::prepare(setup, err);
        if (!job) return false;
        if (run.settings.fibreGeneration == "budget" && run.options.fibreSupportMB > 0 && processPrivateMB() > run.options.fibreSupportMB)
            job->releaseSupport();
        run.jobs.push_back(std::move(job));
    }
    return true;
}

// One budget window's fibres per part, and the support it rebuilt and released.
struct EvaluatedWindow {
    std::vector<FibreSet> fibres;
    int prepares = 0, releases = 0;
};

// Every part of `jobs` with candidates in the window: its support rebuilt if it waits released, its fibres
// generated on the GPU, and its support released again while the process holds more than the bound
// (FibreRunOptions::fibreSupportMB; the teddy's 18 parts hold ~3 GB of surface model at the 0.025 mm rings
// together, one is rebuilt in ~0.1-0.5 s). In parallel over the parts, unless the support is bounded: then
// one part's support is built at a time.
EvaluatedWindow evaluateWindow(const WoolRun& run, const std::vector<size_t>& jobs, const KeyWindow& window) {
    const bool bounded = run.options.fibreSupportMB > 0;
    ExecPolicy policy;
    policy.parallel = run.options.inner.parallel && !bounded;
    EvaluatedWindow out;
    out.fibres.resize(jobs.size());
    parallelFor(static_cast<int64_t>(jobs.size()), 1, [&](int64_t begin, int64_t end, int64_t) {
        for (int64_t k = begin; k < end; ++k) {
            const size_t p = static_cast<size_t>(k);
            if (window.candidates[p] == 0) continue;
            PartFibreJob& job = *run.jobs[jobs[p]];
            std::string why;
            if (job.released()) {
                if (!job.rebuildSupport(&why)) throw std::runtime_error(why);
                ++out.prepares;
            }
            if (!job.generate(window.keyLo, window.keyHi, out.fibres[p], &why)) throw std::runtime_error(why);
            if (bounded && processPrivateMB() > run.options.fibreSupportMB) {
                job.releaseSupport();
                ++out.releases;
            }
        }
    }, policy);
    return out;
}

// Every object's candidates in a budget order of its own (fibres.segment_budget holds per object).
// Only those below its object's key threshold are evaluated; the export is byte-identical to full
// generation, the per-part diagnostics cover the evaluated subset.
void generateInBudgetOrder(WoolRun& run) {
    std::map<int, std::vector<size_t>> jobsOfObject;
    for (size_t i = 0; i < run.jobs.size(); ++i) jobsOfObject[sceneObjectOf(run.scene.source, run.jobs[i]->part().prim)].push_back(i);
    run.covered.resize(run.jobs.size());
    run.q.resize(run.jobs.size());
    for (const auto& [object, jobs] : jobsOfObject) {
        // 1. Every candidate's key into its bucket, per part, in parallel over the parts.
        const auto started = std::chrono::steady_clock::now();
        std::vector<std::vector<int64_t>> perBucket(jobs.size());
        ExecPolicy parts;
        parts.parallel = run.options.inner.parallel;
        parallelFor(static_cast<int64_t>(jobs.size()), 1, [&](int64_t begin, int64_t end, int64_t) {
            for (int64_t k = begin; k < end; ++k) perBucket[static_cast<size_t>(k)] = candidatesPerBucket(run.jobs[jobs[static_cast<size_t>(k)]]->candidates());
        }, parts);
        const double msCount = msSince(started);
        // 2. Windows of whole buckets until the budget is crossed or nothing is left.
        BudgetSchedule schedule(std::move(perBucket), run.settings.segmentBudget);
        int prepares = 0, releases = 0;
        double msEvaluate = 0;
        while (!schedule.finished()) {
            cancelCheckpoint();
            const KeyWindow& window = schedule.nextWindow();
            const auto evaluating = std::chrono::steady_clock::now();
            EvaluatedWindow evaluated = evaluateWindow(run, jobs, window);
            msEvaluate += msSince(evaluating);
            prepares += evaluated.prepares;
            releases += evaluated.releases;
            schedule.addWindow(std::move(evaluated.fibres));
            appendFormat(run.out.log, "fibres: %s\n", schedule.windowLine().c_str());
        }
        cancelCheckpoint();  // a cancellation during the last window publishes nothing either
        // 3. Every part's fibres in candidate order, with what they cover.
        BudgetRun budget = schedule.finish();
        appendFormat(run.out.log,
                     "fibres: %s: budget %lld: q %.6f (crossing bucket %d), %s; candidates %lld, covered %lld (key < %.6f), evaluated %lld; "
                     "%d batches; count %.0f ms, evaluate %.0f ms, total %.0f ms; support prepares %d, releases %d\n",
                     sceneObjectName(run.scene.source, run.jobs[jobs[0]]->part().prim).c_str(), budget.budget, double(budget.q),
                     budget.crossingBucket, budget.complete ? "complete" : "partial", static_cast<long long>(budget.candidates),
                     static_cast<long long>(budget.covered), double(budget.keyBound), static_cast<long long>(budget.evaluated), budget.batches,
                     msCount, msEvaluate, msSince(started), prepares, releases);
        for (size_t k = 0; k < jobs.size(); ++k) {
            run.sets[jobs[k]] = std::move(budget.parts[k].fs);
            run.covered[jobs[k]] = budget.parts[k];
            run.covered[jobs[k]].keyBound = budget.keyBound;  // the object's common coverage
            run.q[jobs[k]] = budget.q;
        }
    }
}

bool generateEveryCandidate(WoolRun& run, std::string* err) {
    for (size_t i = 0; i < run.jobs.size(); ++i) {
        cancelCheckpoint();
        run.progress("fibres: generating", static_cast<int>(i));
        if (!run.jobs[i]->generate(0.0f, 2.0f, run.sets[i], err)) return false;
    }
    return true;
}

// One part's lines of the run's log.
void logPartFibres(WoolRun& run, const ScenePart& part, const FinishedFibres& done) {
    std::string& log = run.out.log;
    const FibreReport& report = done.fibres.rep;
    const char* name = part.name.c_str();
    if (run.options.fibreCheck)
        appendFormat(log, "  %s fibres: %d candidates -> %d accepted (directed %d, random %d); angle to the yarn at emergence p50 %.1f / max %.1f deg, "
                          "along max p95 %.1f / max %.1f deg (directed; nominal lift %.0f..%.0f)\n",
                     name, report.candidates, report.accepted, done.angles.pop[0].fibres, done.angles.pop[1].fibres, done.angles.pop[0].emergence[0],
                     done.angles.pop[0].emergence[2], done.angles.pop[0].along[1], done.angles.pop[0].along[2], done.params.liftMinDeg, done.params.liftMaxDeg);
    else
        appendFormat(log, "  %s fibres: %d candidates -> %d accepted (angles not measured)\n", name, report.candidates, report.accepted);
    // Why the others were not: "hidden" is inside any of the scene's solids (FibreFrame::hidden).
    appendFormat(log, "  %s fibres: rejected %d back-facing, %d hidden inside a solid, %d origin inside, %d occluded, %d collision, %d truncated, "
                      "%d no bend, %d too sharp, %d approximation, %d no emergence, %d root, %d rules\n",
                 name, report.backFacing, report.hiddenInPiece, report.originInside, report.occluded, report.collisionRejected, report.truncated,
                 report.noBend, report.turnTooSharp, report.approxError, report.noEmergence, report.rootRejected, report.ruleRejected);
    if (run.options.fibreCheck) log += done.checks;
    if (!done.gpuLine.empty()) log += "    " + done.gpuLine + "\n";
    log += done.profile + done.digest;
}

// What the run's closing lines count.
struct FibreTotals {
    long long candidates = 0, retained = 0, segments = 0;
    double anglesSeconds = 0;
};

// Angles, checks and digest of every part; each part's support is freed as it finishes.
void finishParts(WoolRun& run, FibreTotals& totals) {
    const bool budgetMode = run.out.budgetMode;
    for (size_t i = 0; i < run.jobs.size(); ++i) {
        cancelCheckpoint();
        run.progress("fibres: finishing parts", static_cast<int>(i));
        const ScenePart& part = run.jobs[i]->part();
        FinishedFibres done = run.jobs[i]->finish(std::move(run.sets[i]));
        float coveredBelow = 2.0f;
        if (budgetMode) {
            const BudgetPartResult& covered = run.covered[i];
            coveredBelow = covered.keyBound;
            appendFormat(run.out.log, "  %s fibres: budget-ordered subset key < %.6f: %lld of %lld candidates evaluated\n", part.name.c_str(),
                         double(covered.keyBound), static_cast<long long>(covered.covered), static_cast<long long>(covered.candidates));
        }
        run.jobs[i].reset();
        logPartFibres(run, part, done);
        totals.candidates += done.fibres.rep.candidates;
        totals.retained += static_cast<long long>(done.fibres.fibres.size());
        totals.segments += segmentCount(done.fibres);
        totals.anglesSeconds += done.anglesSeconds;
        PartFibres pf{static_cast<int>(&part - &run.scene.parts[0]), part.name, std::move(done.fibres), {}};
        pf.coveredBelow = coveredBelow;
        run.out.parts.push_back(std::move(pf));
    }
}

} // namespace

std::optional<SceneFibres> growSceneFibres(const SceneEntry& entry, const Scene& scene, const SceneSettings& settings, const FibreRunOptions& options,
                                           std::string* err) {
    switch (workflowOf(settings)) {
    case Workflow::Fur: return growSceneFur(entry, scene, settings, err);
    case Workflow::Knit: return growWoolFibres(scene, settings, options, err);
    }
    return std::nullopt;
}

std::optional<SceneFibres> growWoolFibres(const Scene& scene, const SceneSettings& settings, const FibreRunOptions& options, std::string* err) {
    FibreParams base;
    if (!fibreParamsFromSettings(settings, base, err)) return std::nullopt;
    SceneFibres out;
    WoolRun run{scene, settings, options, out};
    for (const ScenePart& part : scene.parts) run.knitted += part.hasStitches();
    // The wall time of the three phases over all parts is always logged: where a slow run spends it.
    const auto preparing = std::chrono::steady_clock::now();
    if (!prepareParts(run, base, err)) return std::nullopt;
    const double prepareSeconds = secondsSince(preparing), preparedMB = processPrivateMB();

    out.budgetMode = settings.fibreGeneration == "budget";
    appendFormat(out.log, "fibres: %s generation on the GPU\n", out.budgetMode ? "budget-ordered" : "full");
    // The drape is a rule of the light checks (cutAtDeepContact); the full checks' collision cut has no such rule.
    if (base.drape && !base.lightChecks) out.log += "fibres: drape applies to light checks only; full checks unchanged\n";
    const auto generating = std::chrono::steady_clock::now();
    run.progress("fibres: generating", 0);
    run.sets.resize(run.jobs.size());
    const bool scheduled = out.budgetMode && !run.jobs.empty();
    if (scheduled) generateInBudgetOrder(run);
    else if (!generateEveryCandidate(run, err)) return std::nullopt;
    const double generateSeconds = secondsSince(generating), generatedMB = processPrivateMB();

    const auto finishing = std::chrono::steady_clock::now();
    FibreTotals totals;
    finishParts(run, totals);
    // The segment budget's threshold per object. Budget order evaluated only the candidates below the
    // scheduler's threshold; the covered subsets must reproduce it (the scheduler's exactness).
    out.log += applyObjectBudgets(scene, settings.segmentBudget, out.parts);
    for (size_t i = 0; scheduled && i < out.parts.size(); ++i)
        if (out.parts[i].keyBelow != run.q[i]) {
            if (err) *err = "fibres: the budget threshold of the evaluated fibres differs from the budget scheduler's (" + out.parts[i].name + ")";
            return std::nullopt;
        }
    appendFormat(out.log, "fibres: %lld candidates, %lld retained, %lld segments (the segment budget %lld per object)\n", totals.candidates,
                 totals.retained, totals.segments, settings.segmentBudget);
    appendFormat(out.log, "fibres: wall time: prepare %.1f s, generate %.1f s, finish %.1f s (of it measureFibreAngles %.1f s) (%d parts)\n",
                 prepareSeconds, generateSeconds, secondsSince(finishing), totals.anglesSeconds, run.knitted);
    appendFormat(out.log, "fibres: private memory after prepare %.0f MB, after generate %.0f MB, after finish %.0f MB\n", preparedMB, generatedMB,
                 processPrivateMB());
    if (const double peak = peakWorkingSetMB(); peak > 0) appendFormat(out.log, "memory: peak %.0f MB\n", peak);
    return out;
}

} // namespace knit
