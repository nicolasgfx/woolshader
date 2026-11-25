// The Knit stage: load the parts, knit them in order, gap patches, clearance, details.
#include "scene/knit/knit_stage.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <mutex>
#include <stdexcept>
#include <thread>

#include "gpu/device.h"
#include "gpu/knitting.h"
#include "gpu/surfaces.h"
#include "model/cancel.h"
#include "model/log.h"
#include "scene/knit/knit_run.h"
#include "scene/knit/gap_patches.h"
#include "scene/source/scene_parts.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"
#include "scene/scene.h"
#include "scene/run_log.h"

namespace knit {

namespace {

// --stage-profile: one stage's wall time.
void logStage(const KnitRunOptions& options, const StageClock& clock, const char* label) {
    if (options.stageProfile) std::printf("profile: %s\n", clock.line(label).c_str());
}

// The time caps of a knitting run (KnitRunOptions::capS / partCapS). The current part runs under a
// token of its own; a watchdog thread cancels it at the part's deadline (the part is timed out) or
// when the caller's token fires. Without caps there is no thread: the caller's token applies.
class PartWatchdog {
public:
    using Clock = std::chrono::steady_clock;

    PartWatchdog(double capS, double partCapS, const CancelToken* caller)
        : partCapS_(partCapS), caller_(caller), capped_(capS > 0 || partCapS > 0),
          sceneDeadline_(capS > 0 ? after(Clock::now(), capS) : Clock::time_point::max()) {
        if (capped_) thread_ = std::thread([this] { watch(); });
    }
    ~PartWatchdog() {
        watching_ = false;
        if (thread_.joinable()) thread_.join();
    }

    // Past the scene's deadline: the parts not yet started stay so.
    bool sceneExpired() const { return capped_ && Clock::now() >= sceneDeadline_; }

    // One part: watched from here to the end of the scope.
    class PartTimer {
    public:
        explicit PartTimer(PartWatchdog& watchdog) : watchdog_(watchdog) {
            std::lock_guard<std::mutex> lock(watchdog_.mutex_);
            watchdog_.partToken_ = &token_;
            watchdog_.timedOut_ = false;
            watchdog_.partDeadline_ = watchdog_.partCapS_ > 0 ? std::min(watchdog_.sceneDeadline_, after(Clock::now(), watchdog_.partCapS_)) : watchdog_.sceneDeadline_;
        }
        ~PartTimer() {
            std::lock_guard<std::mutex> lock(watchdog_.mutex_);
            watchdog_.partToken_ = nullptr;
        }
        PartTimer(const PartTimer&) = delete;
        PartTimer& operator=(const PartTimer&) = delete;

        const CancelToken* token() const { return watchdog_.capped_ ? &token_ : watchdog_.caller_; }  // the part runs under it
        bool timedOut() const { return watchdog_.timedOut_; }

    private:
        PartWatchdog& watchdog_;
        CancelToken token_;
    };

private:
    static Clock::time_point after(Clock::time_point scene, double seconds) {
        return scene + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(seconds));
    }

    void watch() {
        while (watching_) {
            {
                std::lock_guard<std::mutex> lock(mutex_);
                if (caller_ && caller_->cancelled()) {
                    if (partToken_) partToken_->cancel();
                } else if (partToken_ && Clock::now() >= partDeadline_) {
                    timedOut_ = true;
                    partToken_->cancel();
                }
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
    }

    const double partCapS_;
    const CancelToken* caller_;
    const bool capped_;
    const Clock::time_point sceneDeadline_;
    std::mutex mutex_;
    CancelToken* partToken_ = nullptr;
    Clock::time_point partDeadline_ = Clock::time_point::max();
    std::atomic<bool> timedOut_{false}, watching_{true};
    std::thread thread_;
};

// What a run refuses before anything is loaded. A relaxation needs the solver's linear algebra
// (computeNotes), and the per-piece cache files are the spiral layout's.
bool checkKnitRequest(const SceneEntry& entry, const SceneSettings& settings, const KnitRunOptions& options, std::string* err) {
    if (std::string why; knitWorkflow(settings) && (settings.settle || settings.repair) && !gpu::ncpAvailable(&why)) {
        if (err) *err = "the knitting solver is unavailable (" + why + "): settle and repair cannot run on this machine";
        return false;
    }
    if (!knitsSewnPieces(entry, settings) && (!options.saveKnitCache.empty() || !options.loadKnitCache.empty())) {
        if (err) *err = "scene " + entry.id + ": the per-piece knit cache files are the spiral pieces' (--save-knit / --load-knit); "
                        "this scene is knitted on its mesh (its knits are kept in the scene cache)";
        return false;
    }
    return true;
}

// The run's numerics. Scenes knit on the GPU with no CPU stand-in: null with the reason without a
// usable device.
std::unique_ptr<KnitNumerics> deviceNumerics(const SceneSettings& settings, const KnitRunOptions& options, std::string* err) {
    if (std::string why; !gpu::available(&why)) {
        if (err) *err = "knitting needs a usable CUDA GPU: " + why;
        return nullptr;
    }
    if (options.verbose) std::printf("GPU knitting: %s\n", gpu::deviceSummary().c_str());
    std::string why;
    // A guided knit runs on the solver options its settle is made for (SettleSettings); every other route on the
    // environment's (v5 unless a switch is set).
    gpu::NcpNumericsOptions solver;
    if (knitsGuidedPiece(settings)) {
        const SettleSettings settle = guidedSettleSettings(settings);
        solver.contactStorage = settle.contactStorage;
        solver.curvatureRows = settle.curvatureRows;
    }
    std::unique_ptr<KnitNumerics> numerics = gpu::makeKnitNumerics(&why, solver);
    if (!numerics && err) *err = "GPU knitting: " + why;
    return numerics;
}

// Part `index` by its layout's operation. False for an error of the run; a device failure throws.
bool knitPart(const KnitRun& run, size_t index, ScenePart& part, std::string* err) {
    switch (part.layout) {
    case KnitLayout::Spiral: return knitSpiralPart(run, index, part, err);
    case KnitLayout::UvTiling: knitMeshPart(run, part); return true;
    case KnitLayout::Guided: knitGuidedPart(run, part); return true;
    }
    return true;
}

// A part stopped at its deadline keeps no knit, and its note says so.
void noteTimeCap(ScenePart& part, const KnitRunOptions& options) {
    part.stop = KnitStop::TimeCap;
    part.knit.reset();
    std::string note;
    appendFormat(note, "time cap: stopped after %.0f s (part cap %.0f s, scene cap %.0f s)", part.ms / 1000, options.partCapS, options.capS);
    part.note = note + (part.note.empty() ? "" : "; " + part.note);
}

// The parts from `first` on, in order and one at a time: the solves take the device one after
// another, and parts knitting side by side only added their device memory to the running solve (man,
// 8 jobs: 4053 s against 851 s with one). False with the reason for an error
// of the run: a device failure is one, never a silent CPU fallback. The caller's Cancelled passes
// through.
bool knitParts(const KnitRun& run, size_t first, Scene& scene, PartWatchdog& watchdog, std::string* err) {
    const KnitRunOptions& options = run.options;
    const int total = static_cast<int>(scene.parts.size());
    int done = static_cast<int>(first);
    for (size_t i = first; i < scene.parts.size(); ++i) {
        cancelCheckpoint();
        ScenePart& part = scene.parts[i];
        if (!options.only.empty() && part.name.find(options.only) == std::string::npos) {
            part.note = "skipped (part filter \"" + options.only + "\")";
            ++done;
            continue;
        }
        if (watchdog.sceneExpired()) return true;  // this and the remaining parts: not started
        if (options.onPart) options.onPart(i, PartEvent::Started, done, total);
        const auto started = std::chrono::steady_clock::now();
        PartWatchdog::PartTimer timer(watchdog);
        ThreadCancelScope partScope(timer.token());
        bool timedOut = false;
        try {
            if (!knitPart(run, i, part, err)) return false;
        } catch (const Cancelled&) {
            if (!timer.timedOut()) throw;
            timedOut = true;
        } catch (const std::exception& ex) {
            if (err) *err = ex.what();
            return false;
        }
        part.ms = msSince(started);
        if (timedOut) noteTimeCap(part, options);
        if (options.verbose) {
            const char* verdict = part.knit ? toString(part.knit->verdict) : timedOut ? "time cap" : "not knitted";
            std::printf("  %s: %s (%.0f ms)\n", part.name.c_str(), verdict, part.ms);
            std::fflush(stdout);
        }
        if (options.onPart) options.onPart(i, PartEvent::Finished, ++done, total);
    }
    return true;
}

// Past the scene's deadline: the parts that never ran, filtered ones too. A part that ran has a time.
void markNotStarted(Scene& scene, double capS) {
    for (ScenePart& part : scene.parts)
        if (!part.knit && part.stop == KnitStop::None && part.ms == 0) {
            part.stop = KnitStop::NotStarted;
            part.note = "not started: the scene's time cap (" + std::to_string(static_cast<int>(capS)) + " s) passed";
        }
}

} // namespace

bool knitScene(const SceneEntry& entry, const SceneSettings& settings, const KnitRunOptions& options, Scene& out, std::string* err) {
    if (!checkKnitRequest(entry, settings, options, err)) return false;
    const StageClock loadClock;
    const std::optional<KnitSource> source = loadPartsToKnit(entry, settings, options.source, out, err);
    if (!source) return false;
    logStage(options, loadClock, "load");
    // The inside tests of the lattice nodes (UV charts); the sewn pieces touch by design.
    const SceneSolids none;
    const SceneSolids& solids = hasChartParts(out) ? *out.solids : none;
    std::vector<PlacedPiece> spirals(out.parts.size());
    for (size_t i = 0; i < out.parts.size(); ++i)
        if (out.parts[i].layout == KnitLayout::Spiral) spirals[i] = {out.parts[i].mesh.get(), out.parts[i].centre, out.parts[i].object};
    const std::optional<SpiralKnitCache> cache = spiralKnitCache(settings, source->pieces, options.saveKnitCache, options.loadKnitCache, err);
    if (!cache) return false;
    const std::unique_ptr<KnitNumerics> numerics = deviceNumerics(settings, options, err);
    if (!numerics) return false;
    if (options.onLoaded) options.onLoaded(out);

    const KnitRun run{settings, options, source->pattern, solids, source->pieces, spirals, out.source, *cache, *numerics, source->guided.get()};
    PartWatchdog watchdog(options.capS, options.partCapS, threadCancelToken());
    bool ok = knitParts(run, 0, out, watchdog, err);
    // knit.coverage = overlap, knit.gap_patches: local patches over the regions the chart parts left
    // uncovered, knitted as a second pass of independent parts.
    if (ok && overlapCoverage(settings) && settings.gapPatches && !watchdog.sceneExpired()) {
        const size_t firstPatch = out.parts.size();
        const StageClock gapClock;
        addGapPatches(out, settings, options);
        logStage(options, gapClock, "gap coverage");
        if (out.parts.size() > firstPatch) {
            if (options.onLoaded) options.onLoaded(out);
            ok = knitParts(run, firstPatch, out, watchdog, err);
        }
    }
    if (watchdog.sceneExpired()) markNotStarted(out, options.capS);
    if (!ok) return false;
    const StageClock clearanceClock;
    if (hasChartParts(out) && !crossComponentClearance(out.parts, settings, out.cross, err)) return false;  // the pieces' line says why not (sewn_pieces.h)
    logStage(options, clearanceClock, "cross-part clearance");
    return true;
}

bool crossComponentClearance(const std::vector<ScenePart>& parts, const SceneSettings& settings, CrossComponentReport& out, std::string* err) {
    out = CrossComponentReport{};
    // Sampled centreline polylines in library mm, with each part's final-surface support radius.
    // knit.coverage = overlap: the patches of one component overlap by design (one group, not checked).
    std::vector<gpu::ClearanceSegment> segments;
    std::map<std::pair<int, int>, int> groupOfComponent;
    for (size_t i = 0; i < parts.size(); ++i) {
        const ScenePart& part = parts[i];
        if (!part.knit || part.layout != KnitLayout::UvTiling) continue;  // spiral pieces touch their neighbours by design
        const double radius = finalSurface(settings, part.yarnD).supportRadius();
        const int group = overlapCoverage(settings) ? groupOfComponent.emplace(std::make_pair(part.prim, part.component), static_cast<int>(groupOfComponent.size())).first->second
                                                    : static_cast<int>(i);
        for (const Polyline& line : part.knit->after.polylines)
            for (size_t k = 0; k + 1 < line.p.size(); ++k)
                segments.push_back({line.p[k] + part.centre, line.p[k + 1] + part.centre, radius, static_cast<int>(i), group});
    }
    gpu::SceneClearance found;
    if (!gpu::checkSceneClearance({segments, settings.clearance}, found, err)) return false;
    out.minDistance = found.minDistance;
    out.pairs = static_cast<int>(found.pairs.size());
    for (const gpu::ClosePartPair& pair : found.pairs) {
        out.closeSegments += static_cast<long>(pair.segments);
        std::string line;
        appendFormat(line, "%s vs %s: %lld sample segment pairs closer than R_a + R_b + c = %.3f mm (min %.3f mm)", parts[pair.partA].name.c_str(),
                     parts[pair.partB].name.c_str(), pair.segments, pair.threshold, pair.minDistance);
        out.lines.push_back(line);
    }
    out.status = found.pairs.empty() ? Status::Pass : Status::Fail;
    return true;
}

} // namespace knit
