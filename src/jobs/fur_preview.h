#pragma once
// The Fur step's disc preview (handover rev. 3, section 12): one disc per loop fibre and one needle per fine fibre,
// built live from the Fur settings without growing a curve. A disc is the fibre's control: the circle of radius r
// through its root in the plane (d, c) (fur_field.h). It shows density, variation, curl and direction; it does not
// predict the settle (gravity, collisions) or the exact bend: those wait for Grow fur.
//
// The worker keeps what the edits rarely change and rebuilds the rest (handover section 8, live preview):
//   - the loaded scene, keyed by the Load stage's inputs (scene/import keys); reloaded on a new object;
//   - per part and layer, the candidate roots with the hidden test (fur_field.h FurRootCache), keyed by the seed
//     and a capacity: a density at or below it reuses them (the fields do not depend on the capacity), a higher one
//     resamples at 1.25 x the density (headroom for a drag upward);
//   - per request, the fields and the display buffers: no mesh query, parallel.
// Requests coalesce: the newest wins, a running build of an older one is cancelled; the last complete preview stays
// shown until the next is ready. The viewer uploads the buffers on its own thread (viewer/fur_discs.h).
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "jobs/built_scene.h"
#include "model/cancel.h"
#include "settings/scene_settings.h"

namespace knit {

// What the discs are coloured by (the Fur step's View section).
enum class FurDiscColour { Plain, Turns, Thickness, Tilt, Layer };

// One layer's candidates as the preview counts them (fur_look.h FurLayerCounts, before any settle).
struct FurPreviewCounts {
    long long sampled = 0, hidden = 0, masked = 0, kept = 0;
};

// One preview's display data, immutable once published (scene mm).
struct FurPreview {
    uint64_t generation = 0;          // the request it answers (FurPreviewWorker::request)
    std::shared_ptr<const BuiltScene> body;  // the object's meshes under the discs (made once per loaded scene)
    // Loop discs: one triangle each, circumscribed about the disc's circle; `discLocal` holds each vertex's
    // coordinates in disc radii about the centre (x, y: the circle is x^2 + y^2 <= 1, the fragment test) and the
    // fibre's budget key (a uniform random number: drawing the keys below a limit thins evenly).
    std::vector<float> discPos, discNormal, discLocal;  // 3, 3, 3 floats per vertex
    std::vector<uint32_t> discColour;                    // RGBA8 per vertex (sRGB)
    // The loops' direction marks (from the root along d) and the fine fibres (from the root along d, their length):
    // one line, 2 vertices, each, with the fibre's budget key per vertex.
    std::vector<float> markPos, needlePos;               // 3 floats per vertex
    std::vector<float> markKey, needleKey;               // 1 float per vertex
    std::vector<uint32_t> markColour, needleColour;
    long long loops = 0, fibres = 0;                     // discs and needles in the buffers
    float keyBelow = 2;                                  // < 1: a draft holding only the fibres with keys below it
    FurPreviewCounts loopCounts, fibreCounts;            // per layer
    double msScene = 0, msRoots = 0, msField = 0, msBuffers = 0;  // this preview's work [ms] (0: cached)
    std::string note;                                    // what the preview does not show, or why it is empty
};

class FurPreviewWorker {
public:
    FurPreviewWorker();
    ~FurPreviewWorker();
    FurPreviewWorker(const FurPreviewWorker&) = delete;
    FurPreviewWorker& operator=(const FurPreviewWorker&) = delete;
    // A preview of these settings in this colouring; returns the request's generation. The same request as the
    // last one (settings text and colouring) starts nothing.
    uint64_t request(const SceneSettings& s, FurDiscColour colour);
    std::shared_ptr<const FurPreview> latest() const;   // the last complete preview (any settings)
    std::string error() const;                          // of the last finished request ("" none)
    bool busy() const;
    uint64_t requested() const;                         // the newest request's generation
    struct Cache;                                       // the worker thread's own state (fur_preview.cpp)

private:
    void loop();
    mutable std::mutex m_;
    std::condition_variable cv_;
    std::unique_ptr<SceneSettings> pending_;
    FurDiscColour pendingColour_ = FurDiscColour::Plain;
    bool pendingDraft_ = false;                     // drawn as a draft first (FurPreview::keyBelow)
    std::chrono::steady_clock::time_point lastRequestAt_{};
    std::string lastRequest_;
    uint64_t generation_ = 0, pendingGeneration_ = 0;
    std::shared_ptr<const FurPreview> latest_;
    std::string error_;
    std::shared_ptr<CancelToken> token_;
    std::unique_ptr<Cache> cache_;  // the worker thread's own (scene, roots, mask)
    bool stop_ = false, running_ = false, runningDraft_ = false;
    std::thread worker_;
};

} // namespace knit
