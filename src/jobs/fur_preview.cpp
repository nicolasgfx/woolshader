// The Fur step's disc preview (fur_preview.h has what is cached and why).
#include "jobs/fur_preview.h"

#include <algorithm>
#include <chrono>
#include <cmath>

#include "fibres/fur_field.h"
#include "jobs/display_scene.h"
#include "model/colour.h"
#include "model/log.h"
#include "scene/appearance/stitch_colours.h"
#include "scene/fibres/fur_stage.h"
#include "scene/result_keys.h"
#include "scene/scene.h"
#include "scene/source/scene_parts.h"
#include "scene/stage_params.h"
#include "settings/scene_registry.h"
#include "settings/settings_mapping.h"

namespace knit {

namespace {

constexpr double kCapacityHeadroom = 1.25;  // a resampled cache holds this x the asked density: a drag upward reuses it
constexpr double kMarkShare = 0.6;          // a direction mark's length, in disc radii
// While the edits come faster than this, a request is drawn as a draft of this share of the fibres first (even:
// the lowest budget keys), and in full once the edits pause for kFullAfterMs.
constexpr int kEditingMs = 300, kFullAfterMs = 250;
constexpr float kDraftShare = 0.25f;

uint32_t rgba8(const Vec3& srgb) {
    const auto byte = [](double v) { return static_cast<uint32_t>(std::lround(std::clamp(v, 0.0, 1.0) * 255)); };
    return byte(srgb.x) | byte(srgb.y) << 8 | byte(srgb.z) << 16 | 0xFF000000u;
}

// A value in [0, 1] on a blue -> white -> red ramp (the attribute colourings).
Vec3 ramp(double t) {
    t = std::clamp(t, 0.0, 1.0);
    return t < 0.5 ? Vec3{2 * t, 2 * t, 1.0} : Vec3{1.0, 2 - 2 * t, 2 - 2 * t};
}

} // namespace

// The worker thread's own state: the loaded scene, each part's root caches and the loaded mask.
struct FurPreviewWorker::Cache {
    std::string sceneKey;
    std::shared_ptr<Scene> scene;
    std::shared_ptr<const BuiltScene> body;  // the scene's meshes for the view
    struct PartRoots {
        bool sampled = false;
        uint64_t seed = 0;
        double loopsCapacity = 0, fibresCapacity = 0;
        bool fibresOn = false;
        FurPartRoots roots;
        std::string modifiersKey;
        FurPartModifiers modifiers;
    };
    std::vector<PartRoots> parts;
    std::string maskKey;
    FurMask mask;
};

namespace {

// The disc colour of one fibre in the chosen colouring (sRGB).
Vec3 discColour(const FurDisc& d, FurDiscColour colour, const Vec3& plain, const FurLook& look) {
    switch (colour) {
    case FurDiscColour::Plain: return plain;
    case FurDiscColour::Turns: return ramp(look.loops.turns > 0 ? d.turns / look.loops.turns : 0);
    case FurDiscColour::Thickness: {
        const double r0 = 0.5 * look.loops.rootDiameter, v = std::max(1e-9, look.variation.thickness);
        return ramp(d.layer == FurLayer::Loops ? 0.5 + 0.5 * (d.rootRadius / r0 - 1) / v : 0.0);
    }
    case FurDiscColour::Tilt: return {0.5 + 0.5 * d.d.x, 0.5 + 0.5 * d.d.y, 0.5 + 0.5 * d.d.z};
    case FurDiscColour::Layer: return d.layer == FurLayer::Loops ? Vec3{0.95, 0.55, 0.2} : Vec3{0.3, 0.85, 0.95};
    }
    return plain;
}

// The preview of one request: scene, roots (cached), fields, buffers.
bool buildFurPreview(const SceneSettings& s, FurDiscColour colour, float keyBelow, FurPreviewWorker::Cache& cache, FurPreview& out,
                     std::string* err) {
    if (!pileSettingsValid(s, err)) return false;
    const SceneEntry* e = sceneEntryFor(s);
    if (!e) {
        if (err) *err = "the object is not loaded (import it again)";
        return false;
    }
    // The scene, reloaded on a new object (the Load stage's inputs).
    const std::string sceneKey = stageInputs(s, Stage::Load) + "plain=" + s.importPlainNodes;
    if (sceneKey != cache.sceneKey || !cache.scene) {
        const auto t0 = std::chrono::steady_clock::now();
        auto scene = std::make_shared<Scene>();
        if (!loadBaseParts(*e, s, SourceOptions{}, *scene, err)) return false;
        auto body = std::make_shared<BuiltScene>();
        body->kind = BuiltScene::Kind::Fibres;
        if (!knitDisplay(*body, *e, *scene, s, err)) return false;
        updateSceneBounds(*body);
        cache.body = std::move(body);
        cache.scene = std::move(scene);
        cache.sceneKey = sceneKey;
        cache.parts.assign(cache.scene->parts.size(), {});
        out.msScene = msSince(t0);
    }
    const Scene& t = *cache.scene;
    out.body = cache.body;
    const FurLook look = furLookFromSettings(s);
    std::string maskPath;
    if (!furMaskPath(s, maskPath, err)) return false;
    const std::string maskKey = maskPath + "|" + s.pile.maskHash;
    if (maskKey != cache.maskKey) {
        cache.mask = FurMask{};
        if (!maskPath.empty() && !loadFurMask(maskPath, cache.mask, err)) return false;
        cache.maskKey = maskKey;
    }
    // Each growing part: its roots (resampled when the seed changed or a density passed the capacity), its field.
    const ExecPolicy exec{true};
    struct PartDiscs {
        FurField field;
        Vec3 plain;  // the part's stored colour, sRGB (the Plain colouring)
    };
    std::vector<PartDiscs> parts;
    SceneSettings objectColour = s;
    objectColour.yarnColourSource = "object";
    for (size_t i = 0; i < t.parts.size(); ++i) {
        cancelCheckpoint();
        if (!furPartGrows(s, t.parts[i])) continue;
        FurPreviewWorker::Cache::PartRoots& pr = cache.parts[i];
        const bool fibresOn = look.fibres.enabled && look.fibres.density > 0;
        const double loopsBound = furDensityBound(look, FurLayer::Loops), fibresBound = furDensityBound(look, FurLayer::Fibres);
        if (!pr.sampled || pr.seed != look.seed || loopsBound > pr.loopsCapacity || (fibresOn && (!pr.fibresOn || fibresBound > pr.fibresCapacity))) {
            const auto t0 = std::chrono::steady_clock::now();
            const bool sameSeed = pr.sampled && pr.seed == look.seed;
            const double loopsCapacity = std::max(sameSeed ? pr.loopsCapacity : 0.0, loopsBound * kCapacityHeadroom);
            const double fibresCapacity = fibresOn ? std::max(sameSeed ? pr.fibresCapacity : 0.0, fibresBound * kCapacityHeadroom) : 0;
            FurPartRoots roots = sampleFurPartRoots(t, i, look, loopsCapacity, fibresCapacity, exec);  // a cancel leaves the cache as it was
            pr.roots = std::move(roots);
            pr.loopsCapacity = loopsCapacity;
            pr.fibresCapacity = fibresCapacity;
            pr.seed = look.seed;
            pr.fibresOn = fibresOn;
            pr.sampled = true;
            out.msRoots += msSince(t0);
        }
        FurCoverage coverage;
        if (!furPartCoverage(t, t.parts[i], cache.mask, coverage, err)) return false;
        // The irregularity fields, rebuilt only when their own inputs change (the seed, the patch scale, the amounts).
        const std::string modifiersKey = std::to_string(look.seed) + "|" + std::to_string(look.variation.density) + "|" +
                                         std::to_string(look.variation.densityScale) + "|" + std::to_string(look.irregular.clumpSize) + "|" +
                                         std::to_string(look.irregular.clumpStrength) + "|" + std::to_string(look.irregular.patchLength) + "|" +
                                         std::to_string(look.irregular.patchLift) + "|" + std::to_string(look.irregular.wear) + "|" +
                                         std::to_string(look.islands.rate) + "|" + std::to_string(look.islands.size) + "|" +
                                         std::to_string(look.islands.strength) + "|" + std::to_string(look.islands.matted) + "|" +
                                         std::to_string(look.islands.fluffy) + "|" + std::to_string(look.islands.thin);
        if (modifiersKey != pr.modifiersKey) {
            const auto t0 = std::chrono::steady_clock::now();
            FurPartModifiers mods = furPartModifiers(t, i, look);
            pr.modifiers = std::move(mods);
            pr.modifiersKey = modifiersKey;
            out.msRoots += msSince(t0);
        }
        const auto t1 = std::chrono::steady_clock::now();
        PartDiscs pd;
        FurFieldRequest request{pr.roots.loops, pr.roots.fibres, look, coverage, &pr.modifiers};
        request.keyBelow = keyBelow;
        pd.field = buildFurField(request, exec);
        out.msField += msSince(t1);
        const Vec3 linear = s.pileColourSource == "single" ? s.pileMaterial.color : partYarnColour(t, t.parts[i], objectColour);
        pd.plain = {srgbEncode(linear.x), srgbEncode(linear.y), srgbEncode(linear.z)};
        const FurField& f = pd.field;
        out.loopCounts.sampled += f.loops.sampled;
        out.loopCounts.hidden += f.loops.hidden;
        out.loopCounts.masked += f.loops.masked;
        out.loopCounts.kept += f.loops.kept;
        out.fibreCounts.sampled += f.fibres.sampled;
        out.fibreCounts.hidden += f.fibres.hidden;
        out.fibreCounts.masked += f.fibres.masked;
        out.fibreCounts.kept += f.fibres.kept;
        parts.push_back(std::move(pd));
    }
    // The buffers: each disc's slot first (serial, cheap), then every disc written in parallel into its slot. The
    // budget key rides along per vertex, so the view thins by key (fur_discs.h) without a sorted order.
    const auto t2 = std::chrono::steady_clock::now();
    struct Slots {
        std::vector<int32_t> slot;  // per disc: its index among the loop discs, or among the needles
    };
    std::vector<Slots> slots(parts.size());
    for (size_t pi = 0; pi < parts.size(); ++pi) {
        const std::vector<FurDisc>& discs = parts[pi].field.discs;
        slots[pi].slot.resize(discs.size());
        for (size_t q = 0; q < discs.size(); ++q) {
            const bool disc = discs[q].layer == FurLayer::Loops && discs[q].radius > 0;
            slots[pi].slot[q] = static_cast<int32_t>(disc ? out.loops++ : out.fibres++);
        }
    }
    const size_t nd = static_cast<size_t>(out.loops), nn = static_cast<size_t>(out.fibres);
    out.discPos.resize(nd * 9);
    out.discNormal.resize(nd * 9);
    out.discLocal.resize(nd * 9);
    out.discColour.resize(nd * 3);
    out.markPos.resize(nd * 6);
    out.markKey.resize(nd * 2);
    out.markColour.resize(nd * 2);
    out.needlePos.resize(nn * 6);
    out.needleKey.resize(nn * 2);
    out.needleColour.resize(nn * 2);
    // The equilateral triangle about the unit circle: corners 2 radii from the centre at 90, 210 and 330 deg.
    constexpr double kCornerX[3] = {0.0, -1.7320508075688772, 1.7320508075688772}, kCornerY[3] = {2.0, -1.0, -1.0};
    const auto put = [](float* to, const Vec3& v) {
        to[0] = float(v.x);
        to[1] = float(v.y);
        to[2] = float(v.z);
    };
    for (size_t pi = 0; pi < parts.size(); ++pi) {
        const std::vector<FurDisc>& discs = parts[pi].field.discs;
        const std::vector<int32_t>& slot = slots[pi].slot;
        const Vec3 plain = parts[pi].plain;
        parallelFor(static_cast<int64_t>(discs.size()), 4096, [&](int64_t b0, int64_t e0, int64_t) {
            cancelCheckpoint();
            for (int64_t q = b0; q < e0; ++q) {
                const FurDisc& d = discs[static_cast<size_t>(q)];
                const size_t k = static_cast<size_t>(slot[static_cast<size_t>(q)]);
                const uint32_t c = rgba8(discColour(d, colour, plain, look));
                const float key = d.budgetKey;
                if (d.layer == FurLayer::Loops && d.radius > 0) {
                    const Vec3 centre = d.root + d.c * d.radius, n = cross(d.d, d.c);
                    for (int v = 0; v < 3; ++v) {
                        put(&out.discPos[k * 9 + v * 3], centre + (d.d * kCornerX[v] + d.c * kCornerY[v]) * d.radius);
                        put(&out.discNormal[k * 9 + v * 3], n);
                        out.discLocal[k * 9 + v * 3] = float(kCornerX[v]);
                        out.discLocal[k * 9 + v * 3 + 1] = float(kCornerY[v]);
                        out.discLocal[k * 9 + v * 3 + 2] = key;
                        out.discColour[k * 3 + v] = c;
                    }
                    put(&out.markPos[k * 6], d.root);
                    put(&out.markPos[k * 6 + 3], d.root + d.d * (kMarkShare * d.radius));
                    out.markKey[k * 2] = out.markKey[k * 2 + 1] = key;
                    out.markColour[k * 2] = out.markColour[k * 2 + 1] = c;
                } else {
                    put(&out.needlePos[k * 6], d.root);
                    put(&out.needlePos[k * 6 + 3], d.root + d.d * (d.radius > 0 ? 0.0 : d.straightLength));
                    out.needleKey[k * 2] = out.needleKey[k * 2 + 1] = key;
                    out.needleColour[k * 2] = out.needleColour[k * 2 + 1] = c;
                }
            }
        }, exec);
    }
    out.msBuffers = msSince(t2);
    out.note = "Discs show each fibre's root, plane, size and direction before the settle; gravity and the exact bend show after Grow fur.";
    out.keyBelow = keyBelow;
    if (out.loops + out.fibres == 0) out.note = "No fibre grows here (coverage, excluded parts or density).";
    return true;
}

} // namespace

FurPreviewWorker::FurPreviewWorker() : cache_(std::make_unique<Cache>()) { worker_ = std::thread([this] { loop(); }); }

FurPreviewWorker::~FurPreviewWorker() {
    {
        std::lock_guard<std::mutex> lk(m_);
        stop_ = true;
        if (token_) token_->cancel();
    }
    cv_.notify_all();
    worker_.join();
}

uint64_t FurPreviewWorker::request(const SceneSettings& s, FurDiscColour colour) {
    const std::string text = writeSettings(s) + "|colour=" + std::to_string(static_cast<int>(colour));
    const auto now = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> lk(m_);
    if (text == lastRequest_) return generation_;
    lastRequest_ = text;
    pending_ = std::make_unique<SceneSettings>(s);
    pendingColour_ = colour;
    pendingGeneration_ = ++generation_;
    pendingDraft_ = now - lastRequestAt_ < std::chrono::milliseconds(kEditingMs);  // edits in a row: a draft first
    lastRequestAt_ = now;
    // The newest wins: a full build is cancelled for it; a running draft finishes (it is quick, and a drag that
    // cancelled every draft would show none), then the newest pending request is drawn.
    if (running_ && !runningDraft_ && token_) token_->cancel();
    cv_.notify_all();
    return generation_;
}

std::shared_ptr<const FurPreview> FurPreviewWorker::latest() const {
    std::lock_guard<std::mutex> lk(m_);
    return latest_;
}

std::string FurPreviewWorker::error() const {
    std::lock_guard<std::mutex> lk(m_);
    return error_;
}

bool FurPreviewWorker::busy() const {
    std::lock_guard<std::mutex> lk(m_);
    return running_ || pending_ != nullptr;
}

uint64_t FurPreviewWorker::requested() const {
    std::lock_guard<std::mutex> lk(m_);
    return generation_;
}

void FurPreviewWorker::loop() {
    for (;;) {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return stop_ || pending_ != nullptr; });
        if (stop_) return;
        std::unique_ptr<SceneSettings> s = std::move(pending_);
        const FurDiscColour colour = pendingColour_;
        const uint64_t generation = pendingGeneration_;
        const bool draft = pendingDraft_;
        running_ = true;
        runningDraft_ = draft;
        token_ = std::make_shared<CancelToken>();
        std::shared_ptr<CancelToken> token = token_;
        lk.unlock();
        auto preview = std::make_shared<FurPreview>();
        preview->generation = generation;
        std::string err;
        bool ok = false, cancelled = false;
        try {
            ThreadCancelScope scope(token.get());
            ok = buildFurPreview(*s, colour, draft ? kDraftShare : 2.0f, *cache_, *preview, &err);
        } catch (const Cancelled&) {
            cancelled = true;
        } catch (const std::exception& e) {
            err = e.what();
        }
        lk.lock();
        running_ = false;
        if (cancelled) continue;  // a newer request is pending; the shown preview stays
        error_ = ok ? std::string() : err;
        if (ok) latest_ = std::move(preview);
        // A draft is drawn in full once the edits pause (the same request; a newer one replaces it).
        if (ok && draft && !cv_.wait_for(lk, std::chrono::milliseconds(kFullAfterMs), [&] { return stop_ || pending_ != nullptr; })) {
            pending_ = std::move(s);
            pendingColour_ = colour;
            pendingGeneration_ = generation;
            pendingDraft_ = false;
        }
    }
}

} // namespace knit
