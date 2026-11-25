// One part's fibres: prepare, generate, finish.
#include "scene/fibres/part_fibre_job.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "gpu/fibre_generator.h"
#include "gpu/surfaces.h"
#include "knitting/tube/param_text.h"
#include "model/cancel.h"
#include "model/compute_caps.h"
#include "model/coordinates.h"
#include "model/hash.h"
#include "model/log.h"
#include "model/parallel.h"
#include "scene/source/scene_solids.h"
#include "settings/settings_mapping.h"

namespace knit {

namespace {

// A part's final fibre parameters and frame.
struct PartFrame {
    FibreParams params;
    FibreFrame frame;
    double serializedMargin = 0;
};

// The pieces' setup (root depth, exposure, serialized root check), on the final surface
// (FinalSurfaceSpec: D_r, compression; or the plies).
void setPartFrame(const ScenePart& part, int index, const RenderedPart& rendered, const FibreParams& base, PartFrame& out) {
    out.params = base;
    FibreSetup setup;
    setup.tubeRadius = rendered.strandPp.radius();
    setup.sides = rendered.sides;
    setup.ring = rendered.ring;
    setup.kappaMax = rendered.kappaMax;
    setup.pieceIndex = index;
    setup.strands = rendered.strands;
    setup.stitches = &rendered.yarn->topo;
    if (rendered.plies) {
        setup.plies = rendered.plies.get();
        setup.parentCurves = rendered.yarn;
        setup.parentSurface = rendered.parentSm.get();
    }
    setupPieceFibres(out.params, out.frame, part.knit->t.pp, setup);
    const SceneMapping map = sceneMetres(part.centre);
    Vec3 lo, hi;
    part.mesh->bounds(lo, hi);
    double maxCoord = 0;
    for (const Vec3& corner : {lo, hi}) {
        const Vec3 q = map.toExport(corner);
        maxCoord = std::max({maxCoord, std::fabs(q.x), std::fabs(q.y), std::fabs(q.z)});
    }
    out.serializedMargin = 2 * maxCoord * 1.1920929e-7 * map.sigma;
    out.params.serializedMargin = out.serializedMargin;
    out.params.exportChordFactor = rendered.ring * rendered.ring / 8;
}

} // namespace

PartFibreJob::PartFibreJob(const PartFibreSetup& setup)
    : part_(setup.part), solids_(setup.solids), options_(setup.options), modelRing_(setup.modelRing) {}

PartFibreJob::~PartFibreJob() = default;

SurfaceParams PartFibreJob::surfaceParams() const {
    SurfaceParams surface;
    surface.compression = compressed_;
    surface.maxCompression = rendered_.fs.maxCompression;
    return surface;
}

std::unique_ptr<PartFibreJob> PartFibreJob::prepare(const PartFibreSetup& setup, std::string* err) {
    std::unique_ptr<PartFibreJob> job(new PartFibreJob(setup));
    StageClock clock;
    const auto stageDone = [&](const char* label) {
        if (setup.options.stageProfile) job->profile_ += "    profile: " + clock.line(label) + "\n";
        clock = StageClock{};
    };
    std::optional<RenderedPart> rendered = renderedPart(setup.part, setup.settings, err);
    if (!rendered) return nullptr;
    job->rendered_ = std::move(*rendered);
    stageDone("fibres: renderedPart");
    job->map_ = sceneMetres(setup.part.centre);
    job->compressed_ = job->rendered_.fs.compression && !job->rendered_.plies;
    PartFrame frame;
    setPartFrame(setup.part, setup.index, job->rendered_, setup.base, frame);
    job->params_ = frame.params;
    job->frame_ = std::move(frame.frame);
    job->serializedMargin_ = frame.serializedMargin;
    stageDone("fibres: setupPieceFibres");
    // The surface model on the device (float32 rings, compressed or round). A device failure fails the
    // fibres, never a CPU stand-in.
    std::string gpuError;
    if (!gpu::buildSurface({*job->rendered_.strands, job->rendered_.strandPp, job->surfaceParams()}, job->surface_, nullptr, &gpuError)) {
        if (err) *err = setup.part.name + ": fibre surface on the GPU: " + gpuError;
        return nullptr;
    }
    stageDone("fibres: buildSurface (GPU)");
    // The index's users: the GPU audit and the fibre checks (fibreCheck).
    if (setup.options.fibreCheck) {
        job->buildTubes();
        stageDone(job->tubes_->sharedSurface() ? "fibres: ExportTubeIndex (surface model shared)" : "fibres: ExportTubeIndex");
    }
    job->source_ = candidateSource(job->surface_, job->params_, &job->frame_);  // the candidates the GPU enumerates
    return job;
}

void PartFibreJob::buildTubes() {
    PartFibreJob* const self = this;
    ExportTubeSpec spec{modelRing_, rendered_.ring, rendered_.sides, compressed_, [self](const Vec3& q) { return self->map_.toExportF(q); }, serializedMargin_,
                        rendered_.fs.maxCompression};
    tubes_ = std::make_unique<ExportTubeIndex>(*rendered_.strands, rendered_.strandPp, std::move(spec), options_.inner,
                                               surface_.yarns.empty() ? nullptr : &surface_,
                                               surfaceBuildDigest(*rendered_.strands, rendered_.strandPp, surfaceParams()));
}

void PartFibreJob::releaseSupport() {
    tubes_.reset();
    gpuInputs_.reset();
    SurfaceModel().yarns.swap(surface_.yarns);
    released_ = true;
}

bool PartFibreJob::rebuildSupport(std::string* err) {
    // What prepare built: the surface model on the device and, with the checks, the export tube index.
    // The GPU generator rebuilds its flattened inputs itself.
    if (!gpu::buildSurface({*rendered_.strands, rendered_.strandPp, surfaceParams()}, surface_, nullptr, err)) return false;
    if (options_.fibreCheck) buildTubes();
    released_ = false;
    return true;
}

bool PartFibreJob::generate(float keyLo, float keyHi, FibreSet& out, std::string* err) {
    StageClock clock;
    // GPU fibres: the frame's callbacks as device data.
    gpu::FibreFrameGpu frame = gpu::frameFrom(frame_);
    frame.outward = &part_.knit->surface->mesh();
    frame.hidden = solids_.queries();
    frame.hiddenScale = 1;
    frame.hiddenOffset = part_.centre;
    frame.quantize = true;
    frame.quantSigma = map_.sigma;
    frame.quantCentre = map_.centre;
    frame.checkpoint = [] { cancelCheckpoint(); };
    frame.keyLo = keyLo;
    frame.keyHi = keyHi;
    frame.inputs = &gpuInputs_;
    gpu::FibreRunReport report;
    if (!gpu::generateFibres({*rendered_.strands, surface_, rendered_.strandPp, params_, frame}, out, &report, nullptr, err)) return false;
    checks_ += "    " + report.line() + "\n";
    if (!gpuLine_.empty()) gpuLine_ += "\n    ";  // one line per key window (budget order)
    if (!(keyHi > 1.0f && keyLo <= 0.0f)) appendFormat(gpuLine_, "keys [%.6f, %.6f): ", double(keyLo), double(keyHi));
    gpuLine_ += report.line();
    if (options_.stageProfile) profile_ += "    profile: " + clock.line("fibres: generate (GPU)") + "\n";
    if (!options_.fibreCheck) return true;
    clock = StageClock{};
    // In parallel: on 3 plies one accepts() costs 10-20 ms
    // (teddy piece 7: 82,551 fibres, 133 s at 12 threads busy; serially the article's teddy command
    // ran past its 1200 s cap on one thread). ExportTubeIndex::accepts holds no mutable state; each
    // fibre writes its own flag, so the count is the serial one.
    const ExportTubeIndex& tubes = *tubes_;
    std::vector<char> rejects(out.fibres.size(), 0);
    parallelFor(static_cast<int64_t>(out.fibres.size()), 2048, [&rejects, &out, &tubes](int64_t begin, int64_t end, int64_t) {
        cancelCheckpoint();
        for (int64_t q = begin; q < end; ++q) rejects[static_cast<size_t>(q)] = !tubes.accepts(out, out.fibres[static_cast<size_t>(q)], 1e-4);
    }, options_.inner);
    long rejected = 0;
    for (char c : rejects) rejected += c;
    appendFormat(checks_, "    GPU fibres audit (CPU, not applied): the serialized root check would reject %ld of %zu\n", rejected, out.fibres.size());
    if (options_.stageProfile) profile_ += "    profile: " + clock.line("fibres: audit of the omitted root check (CPU)") + "\n";
    return true;
}

FinishedFibres PartFibreJob::finish(FibreSet fibres) {
    FinishedFibres out;
    out.fibres = std::move(fibres);
    out.params = params_;
    out.checks = checks_;
    out.profile = profile_;
    out.gpuLine = gpuLine_;
    if (options_.fibreCheck) {
        const StageClock clock;
        const auto started = std::chrono::steady_clock::now();
        out.angles = measureFibreAngles(*rendered_.yarn, out.fibres, options_.inner);
        out.anglesSeconds = secondsSince(started);
        if (options_.stageProfile) out.profile += "    profile: " + clock.line("fibres: measureFibreAngles") + "\n";
    }
    std::string checkText;  // fibreDigest: every check result, in order
    if (options_.fibreCheck) {
        if (!tubes_) buildTubes();  // released (bounded residency)
        const AttachmentResult attachment = rootAttachmentCheck(*tubes_, out.fibres, 1e-4 / map_.sigma, options_.inner);
        if (options_.fibreDigest) checkText += checkSnapshot(attachment) + "\n";
        if (options_.stageProfile) appendFormat(out.profile, "    profile: check attachment: %.0f ms\n", attachment.ms);
        appendFormat(out.checks, "    attachment: %d/%d roots inside the exported tube, %d/%d tips outside\n", attachment.rootInside,
                     attachment.fibres, attachment.tipOutside, attachment.fibres);
        if (rendered_.fs.referenceCheckApplicable()) {
            const FibreReferenceResult reference = fibreReferenceCheck(*tubes_, out.fibres, 1e-4, options_.inner);
            if (options_.fibreDigest) checkText += checkSnapshot(reference) + "\n";
            if (options_.stageProfile) appendFormat(out.profile, "    profile: check reference: %.0f ms\n", reference.ms);
            appendFormat(out.checks, "    reference check: %d OK, %d FAIL, %d UNRESOLVED fibres\n", reference.ok, reference.fail, reference.unresolved);
        } else {
            out.checks += "    reference check: UNRESOLVED (not applicable: compressed tubes)\n";
        }
    }
    if (options_.stageProfile) out.profile += "    profile: " + partClock_.line("fibres + checks, part total") + "\n";
    if (options_.fibreDigest)
        appendFormat(out.digest, "    digest (%s): fibres %016llx (%zu), checks %016llx\n", buildConfiguration().c_str(),
                     static_cast<unsigned long long>(fnv1a64(fibreSnapshot(out.fibres))), out.fibres.fibres.size(),
                     static_cast<unsigned long long>(fnv1a64(checkText)));
    return out;
}

double exportModelRing(const SceneSettings& settings) { return renderFlattened(settings) ? 0.025 : kExportModelRing; }

std::optional<std::string> partFibreKey(const ScenePart& part, int index, const SceneSettings& settings, const FibreParams& base, std::string* err) {
    const std::optional<RenderedPart> drawn = renderedPart(part, settings, err);
    if (!drawn) return std::nullopt;
    const RenderedPart& rendered = *drawn;
    PartFrame frame;
    setPartFrame(part, index, rendered, base, frame);
    std::string key;
    appendFormat(key, "strands = %016llx %zu\nexport tube = %a %a %d ", static_cast<unsigned long long>(curvesDigest(*rendered.strands)),
                 rendered.strands->yarns.size(), exportModelRing(settings), rendered.ring, rendered.sides);
    return key + "\nstrand params = " + patternParamsText(rendered.strandPp) + "\nfibre params = " + fibreParamsText(frame.params) +
           "\nfibre placement = legacy\n" + fibreBackendLine();
}

std::string fibreBackendLine() {
    // v3 (2026-11-14): the fabric-back test bounds far samples by a slope of 2 instead of searching each (long strays).
    return "fibre backend = gpu float32 (gpu fibres v3: none behind the fabric, far samples bounded by slope 2, no serialized root check)\n";
}

} // namespace knit
