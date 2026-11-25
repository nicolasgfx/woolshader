// The live preview's residency contract through the real path: RenderPreviewWorker starts a file-loaded
// `ppt serve`, and a wool GPU job's PriorityLease takes the memory back (the live preview's
// feasibility study, the owner's G constraints). Needs the device, a ppt.exe with `serve` (KNIT_PREVIEW_PPT) and a package
// (KNIT_PREVIEW_SCENE, default: the swatch project); SKIPs otherwise. One case runs the whole drill, in
// order, and prints its measurements.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <optional>
#include <string>
#include <thread>
#include <vector>

#include <cuda_runtime.h>

#include "gpu_common.h"
#include "jobs/render_preview.h"
#include "model/env.h"

using namespace knit;

namespace {

using Clock = std::chrono::steady_clock;
constexpr size_t MiB = size_t(1) << 20;

double msSince(Clock::time_point t0) { return std::chrono::duration<double, std::milli>(Clock::now() - t0).count(); }

std::string envOr(const char* name, const char* fallback) {
    const std::optional<std::string> v = envVar(name);
    return v && !v->empty() ? *v : fallback;
}

// Polls the worker until `done(state)` or the timeout; returns the last state.
RenderPreviewState waitFor(RenderPreviewWorker& w, double seconds, const std::function<bool(const RenderPreviewState&)>& done) {
    const auto t0 = Clock::now();
    RenderPreviewState s = w.state();
    while (!done(s) && msSince(t0) < seconds * 1000) {
        std::this_thread::sleep_for(std::chrono::milliseconds(20));
        s = w.state();
    }
    return s;
}

// The job: `bytes` of CUDA memory in 1 GiB blocks, each written (resident). False if refused.
bool runJob(size_t bytes) {
    std::vector<void*> blocks;
    bool ok = true;
    for (size_t left = bytes; left > 0 && ok;) {
        const size_t n = left < (size_t(1) << 30) ? left : (size_t(1) << 30);
        void* p = nullptr;
        ok = cudaMalloc(&p, n) == cudaSuccess && cudaMemset(p, 0x5a, n) == cudaSuccess;
        if (p) blocks.push_back(p);
        left -= n;
    }
    cudaDeviceSynchronize();
    for (void* p : blocks) cudaFree(p);
    cudaGetLastError();
    return ok;
}

} // namespace

TEST(render_preview_lease_drill) {
    REQUIRE_GPU();
    const std::string ppt = envOr("KNIT_PREVIEW_PPT", "");
    const std::string scene = envOr("KNIT_PREVIEW_SCENE", "");
    if (ppt.empty() || !std::filesystem::exists(ppt) || !std::filesystem::exists(scene)) SKIP("KNIT_PREVIEW_PPT (ppt with serve) or the scene is missing");

    RenderPreviewWorker worker;
    RenderPreviewRequest r;
    r.enabled = true;
    r.ppt = ppt;
    r.logPath = (std::filesystem::temp_directory_path() / "knit-preview-test-ppt.log").string();
    r.scenePath = scene;
    r.sceneKey = 1;

    // 1. The preview renders.
    auto t0 = Clock::now();
    worker.request(r);
    RenderPreviewState s = waitFor(worker, 180, [](const RenderPreviewState& x) { return (x.frameNumber > 0 && x.spp >= 4) || x.phase == "refused"; });
    // The swatch (39.9 M segments, ~8.6 GiB) is at the edge of a 12 GB card: when other processes hold
    // memory, the renderer's admission refuses it -- the contract working, but no drill to run.
    static std::string busy;  // Skipped keeps the pointer
    if (s.phase == "refused" && s.frameNumber == 0) {
        busy = "the card is too busy for the preview now: " + s.message;
        SKIP(busy.c_str());
    }
    std::printf("start + load -> rendering: %.1f s, %d spp, %.1f ms/spp, renderer %.0f MiB\n", msSince(t0) / 1000, s.spp, s.msPerSpp,
                s.processBytes / double(MiB));
    CHECK(s.frameNumber > 0);
    const long long held = s.processBytes;
    CHECK(held > static_cast<long long>(1000 * MiB));

    // 2. A GPU job that does not fit while the preview holds its scene: graceful release.
    const size_t jobBytes = static_cast<size_t>(held) - 1024 * MiB;
    {
        const gpu::MemoryAvailability a = gpu::memoryAvailability();
        std::printf("before the lease: budget known %d, budget %.0f MiB, usage %.0f MiB, physically free %.0f MiB, available %.0f MiB\n",
                    a.budgetKnown, a.budgetBytes / double(MiB), a.usageBytes / double(MiB), a.physicalFreeBytes / double(MiB),
                    a.availableBytes / double(MiB));
    }
    t0 = Clock::now();
    {
        gpu::PriorityLease lease(jobBytes);
        const double grantMs = msSince(t0);
        const RenderPreviewWorker::ReleaseReport rep = worker.lastRelease();
        std::printf("lease of %.0f MiB granted after %.0f ms: graceful release %.0f ms, terminated %d, freed %.0f MiB\n", jobBytes / double(MiB),
                    grantMs, rep.ms, rep.terminated, rep.freedBytes / double(MiB));
        CHECK(!rep.terminated);
        CHECK(rep.ms < gpu::kGracefulReleaseMs);
        CHECK(lease.freedBytes() >= jobBytes);
        CHECK(gpu::memoryAvailability().availableBytes >= jobBytes);
        const uint64_t frames = worker.state().frameNumber;
        CHECK(runJob(jobBytes));
        std::this_thread::sleep_for(std::chrono::milliseconds(1500));
        CHECK(worker.state().frameNumber == frames);  // nothing loads or renders while the lease lives
    }
    // 3. After the lease: the preview loads again by itself.
    t0 = Clock::now();
    const uint64_t before = worker.state().frameNumber;
    s = waitFor(worker, 180, [&](const RenderPreviewState& x) { return x.frameNumber > before && x.phase == "rendering"; });
    std::printf("after the lease -> rendering again: %.1f s\n", msSince(t0) / 1000);
    CHECK(s.frameNumber > before);

    // 4. A lease while the renderer is loading (a load cannot be interrupted): terminate after the
    // graceful timeout, return only when the exit is confirmed.
    r.sceneKey = 2;
    worker.request(r);
    waitFor(worker, 10, [](const RenderPreviewState& x) { return x.phase == "loading"; });
    std::this_thread::sleep_for(std::chrono::milliseconds(300));
    t0 = Clock::now();
    {
        gpu::PriorityLease lease(jobBytes);
        const RenderPreviewWorker::ReleaseReport rep = worker.lastRelease();
        std::printf("lease during a load granted after %.0f ms: release %.0f ms, terminated %d, freed %.0f MiB\n", msSince(t0), rep.ms,
                    rep.terminated, rep.freedBytes / double(MiB));
        CHECK(rep.terminated);
        CHECK(rep.ms < gpu::kGracefulReleaseMs + 5000);
        CHECK(gpu::memoryAvailability().availableBytes >= jobBytes);
        CHECK(runJob(jobBytes));
    }
    t0 = Clock::now();
    s = waitFor(worker, 240, [](const RenderPreviewState& x) { return x.phase == "rendering" && x.spp >= 2; });
    std::printf("after the terminating lease -> restart + load -> rendering: %.1f s\n", msSince(t0) / 1000);
    CHECK(s.phase == "rendering");

    // 5. Admission: a preview that does not fit keeps the last image and says why.
    const uint64_t framesBefore = worker.state().frameNumber;
    r.sceneKey = 3;
    r.forecastBytes = 64LL << 30;
    worker.request(r);
    s = waitFor(worker, 30, [](const RenderPreviewState& x) { return x.phase == "refused"; });
    std::printf("refused: %s\n", s.message.c_str());
    CHECK(s.phase == "refused");
    CHECK(s.frame != nullptr);  // the last image stays
    CHECK(s.frameNumber >= framesBefore);
    const uint64_t framesAtRefusal = s.frameNumber;
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    CHECK(worker.state().frameNumber == framesAtRefusal);  // nothing renders after the refusal

    // 6. Off: the process ends.
    r.enabled = false;
    worker.request(r);
    s = waitFor(worker, 15, [](const RenderPreviewState& x) { return x.phase == "off"; });
    CHECK(s.phase == "off");
}
