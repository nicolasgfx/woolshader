#pragma once
// The live render preview's lifecycle: one `ppt serve` process (the photon path tracer's preview mode),
// kept warm while the preview is on, and the residency contract with wool's own GPU jobs
// (gpu/device.h):
// - the worker is a ResidencyTenant: a job's PriorityLease may ask it to give the renderer's memory back.
//   It asks the renderer to suspend and drop its scene (graceful); after kGracefulReleaseMs it terminates
//   the process and returns only once the exit is confirmed. The lease then re-reads the budget.
// - while any lease lives, the worker loads, rebuilds and resumes nothing; after the last lease it loads
//   the scene again (mayRebuild).
// - the renderer checks the preview's own admission: a scene whose forecast does not fit its budget is
//   refused, the last frame stays, and the state says why.
// - a new scene replaces the old one release-first (the renderer drops the old mirror before loading).
//
// Threads: a control thread (starts and stops the process, sends what the newest request changed) and a
// reader thread (the frame pipe). request() and state() are cheap and called every viewer frame.
#include <cstdint>
#include <memory>
#include <string>

#include "io/preview_pipe.h"
#include "jobs/jobs.h"

namespace knit {

struct RenderPreviewRequest {
    bool enabled = false;
    bool paused = false;
    std::string ppt;                // ppt.exe
    std::string logPath;            // the renderer's output
    // What to render: a package's scene.json that exists already (sceneKey changes when it changes;
    // the residency drill renders the swatch project this way) ...
    std::string scenePath;
    uint64_t sceneKey = 0;
    // ... or results the worker exports into a package first, on its own thread, after any GPU job
    // (dir, name and format are the worker's). packageKey changes when the results or the materials
    // changed; 0: nothing to export.
    ExportRequest package;
    uint64_t packageKey = 0;
    // The results are about to change (a job runs, or the edits have not paused yet): load nothing, not
    // even the current package again after a GPU job; the last frame stays until the new package is there.
    bool hold = false;
    std::string packageRoot;        // the packages' folder: one subfolder per revision, the last two kept
    // The renderer's compiled textures (textures_compiled/*.ppttex, named by their content's digest), kept
    // across revisions and sessions: each new package starts with them, so a texture is compiled once.
    std::string textureCache;
    long long forecastBytes = 0;    // the preview's memory forecast, checked by the renderer against its budget
    // The camera in the package's frame (metres, Y up). cameraKey 0: the package's own camera.
    double eye[3] = {0, 0, 0}, target[3] = {0, 0, 0}, up[3] = {0, 1, 0};
    double fovDeg = 35;
    bool interacting = false;
    uint64_t cameraKey = 0;
    int width = 512, height = 320;
};

struct RenderPreviewState {
    // off, starting, waiting (a GPU job holds priority), exporting, loading, rendering, converged, paused,
    // released, refused (the preview does not fit), error
    std::string phase = "off";
    std::string message;
    std::shared_ptr<const PreviewFrame> frame;  // the newest frame (kept through refusals and releases)
    uint64_t frameNumber = 0;                   // increments with each new frame
    int spp = 0;
    double msPerSpp = 0;
    long long processBytes = 0, availableBytes = 0;
};

// The worker registers itself as a gpu::ResidencyTenant while its renderer process runs (render_preview.cpp
// holds the adapter, so this header stays free of the gpu module).
class RenderPreviewWorker {
public:
    RenderPreviewWorker();
    ~RenderPreviewWorker();
    RenderPreviewWorker(const RenderPreviewWorker&) = delete;
    RenderPreviewWorker& operator=(const RenderPreviewWorker&) = delete;

    void request(const RenderPreviewRequest& r);  // the newest wins
    RenderPreviewState state() const;

    // The last release: how long the graceful part took, whether the process had to be terminated.
    struct ReleaseReport {
        double ms = 0;
        bool terminated = false;
        long long freedBytes = 0;
    };
    ReleaseReport lastRelease() const;

private:
    struct RendererSession;  // the renderer process, its threads and what it was told (render_preview.cpp)
    std::unique_ptr<RendererSession> session_;
};

} // namespace knit
