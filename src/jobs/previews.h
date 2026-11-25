#pragma once
// What the viewer computes outside the job queue (the queue would cancel a running job): the Fibres step's
// reference preview and the flat-swatch references, each on its own worker; the ply cross-section sketch; the
// hair coverage overlay.
//
// The reference preview: a small sample the controls act on in about a second,
// instead of the whole model (which "Apply to model" computes).
//   Wool: a flat, square patch (5 wales, as many courses as fill the same length) of the current
//     pattern, yarn diameter, gauge and depth
//     (generateStockinette), rendered with the model's own strands (renderedStrands: D_r,
//     flattening, ply tubes) and its twisted-ply appearance as REAL relief on the tube (the maps'
//     profile rho(phi), ply_profile.h - the viewer's fixed-function display cannot show a normal
//     map, so the relief is geometry here); fibres from the model's generator and settings,
//     kept where they root in the central square, 60 % of the side (no patch-edge artefacts).
//   Fur: a sphere cap (radius max(10, 1.5 x longest hair) mm) under -z gravity, grown with the
//     model's fur generator and settings.
// Hairs and fibres are drawn as lit tubes with their real radii and taper (not lines), so
// thickness and taper edits show. Approximate by construction: flat (wool) / a sphere (fur), no
// hidden test, no export quantisation or root check; the note says so.
#include <condition_variable>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "jobs/built_scene.h"
#include "model/cancel.h"
#include "model/choices.h"

namespace knit {

// The viewer's colours of the two wool fibre populations, on the patch and on the model: fibres
// along the yarn yellow, stray fibres blue (display only; the export keeps the materials).
constexpr Vec3 kOrientedFibreColour{1.0, 0.85, 0.1};
constexpr Vec3 kStrayFibreColour{0.15, 0.4, 1.0};

struct FibrePreview {
    uint64_t key = 0;                 // fibrePreviewKey of the settings it was built for
    Workflow workflow = Workflow::Knit;
    std::vector<DisplayMesh> meshes;  // the patch / sphere and the fibre tubes
    Vec3 lo, hi;                      // bounds [mm]
    long long fibres = 0;
    std::string note;
    double ms = 0;
};

// The settings a preview depends on: everything but the scene, the import and the display.
uint64_t fibrePreviewKey(const SceneSettings& s);
bool buildFibrePreview(const SceneSettings& s, FibrePreview& out, std::string* err);

// The reference preview on its own worker: latest wins - a newer request cancels a running build
// and waits 150 ms for further edits.
class FibrePreviewWorker {
public:
    FibrePreviewWorker();
    ~FibrePreviewWorker();
    void request(const SceneSettings& s);
    std::shared_ptr<const FibrePreview> latest() const;  // the last completed preview (any settings)
    std::string error() const;                           // of the last request ("" none)
    bool busy() const;
    uint64_t requestedKey() const;  // the settings' key of the latest request (built, building or waiting)

private:
    void loop();
    mutable std::mutex m_;
    std::condition_variable cv_;
    std::unique_ptr<SceneSettings> pending_;
    uint64_t pendingKey_ = 0, runningKey_ = 0, doneKey_ = 0, requestedKey_ = 0;
    std::shared_ptr<const FibrePreview> latest_;
    std::string error_;
    std::shared_ptr<CancelToken> token_;
    bool stop_ = false, running_ = false;
    std::thread worker_;
};

// A twisted-ply style as its cross-section: the maps' radius profile rho(phi) around a unit yarn
// (knitting/tube/ply_profile.h), as `points` points (x, y) = rho (cos phi, sin phi).
std::vector<Vec2> plyCrossSection(int plies, PlyProfileModel model, double fill, int points);

// The viewer's coverage overlay of a fur result's scene: a copy of `b` whose base-colour textures
// are tinted red where the selected mask is below 0.5 (per texel, in the textures' own uv space),
// and whose excluded parts are drawn red as a whole. Untextured meshes show the part exclusion
// only (the mask needs colour texture coordinates). `covered` receives the share of texels kept.
std::shared_ptr<BuiltScene> furCoverageOverlay(const BuiltScene& b, const SceneSettings& s, double* covered, std::string* err);

// Flat-swatch references of the current D, w, h, A, c on a worker (cached by parameters): the Knit
// stage's pattern availability without a knit.
class FlatReferences {
public:
    FlatReferences();
    ~FlatReferences();
    // The references for these parameters, or null while they are computed (request() starts it).
    std::shared_ptr<const std::vector<PatternReference>> get(const SceneSettings& s);
    void request(const SceneSettings& s);

private:
    static std::string key(const SceneSettings& s);
    void loop();
    std::mutex m_;
    std::condition_variable cv_;
    std::map<std::string, std::shared_ptr<const std::vector<PatternReference>>> done_;
    std::string pending_, running_;
    SceneSettings pendingSettings_;
    bool stop_ = false;
    std::thread worker_;
};

} // namespace knit
