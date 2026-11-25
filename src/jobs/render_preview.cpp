// The live render preview's lifecycle (render_preview.h).
//
// Locks: `m` guards everything below except the stdin writes (`writeM`). A thread never holds `m` while
// it waits on the lease machinery (TenantRegistration's destructor waits for a release in flight, and
// that release needs `m` to see its acknowledgement), so the registration is made and dropped outside `m`.
#include "jobs/render_preview.h"

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <mutex>
#include <thread>

#include "gpu/device.h"
#include "io/child_process.h"
#include "jobs/job_steps.h"
#include "settings/scene_registry.h"
#include "model/cancel.h"

namespace knit {

namespace {

using Clock = std::chrono::steady_clock;

double msSince(Clock::time_point t0) { return std::chrono::duration<double, std::milli>(Clock::now() - t0).count(); }

// A JSON string literal.
std::string quotedJson(const std::string& s) {
    std::string out = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') out += '\\';
        out += c;
    }
    return out + "\"";
}

std::string vec3Json(const double v[3]) {
    char b[128];
    std::snprintf(b, sizeof b, "[%.9g,%.9g,%.9g]", v[0], v[1], v[2]);
    return b;
}

// How long a confirmed exit may take after TerminateProcess.
constexpr int kExitConfirmMs = 5000;

// Copies each compiled texture (*.ppttex) of folder `from` that folder `to` lacks. The names are the
// textures' content digests, so a file of the same name is the same texture. A copy is written under a
// temporary name first: a reader never sees half a file.
void copyCompiledTextures(const std::filesystem::path& from, const std::filesystem::path& to) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (from.empty() || to.empty() || !fs::is_directory(from, ec)) return;
    fs::create_directories(to, ec);
    for (const fs::directory_entry& f : fs::directory_iterator(from, ec)) {
        const fs::path name = f.path().filename();
        std::error_code exists;
        if (f.path().extension() != ".ppttex" || fs::exists(to / name, exists)) continue;
        const fs::path part = to / (name.string() + ".part");
        std::error_code copied, renamed;
        if (fs::copy_file(f.path(), part, fs::copy_options::overwrite_existing, copied)) fs::rename(part, to / name, renamed);
        if (copied || renamed) fs::remove(part, ec);
    }
}

// The renderer's folder of compiled textures inside a package (photon_path_tracer core/ppttex.h).
std::filesystem::path compiledTexturesOf(const std::filesystem::path& packageDir) { return packageDir / "textures_compiled"; }

} // namespace

struct RenderPreviewWorker::RendererSession {
    // The residency adapter: a job's lease calls these on the job's thread.
    struct Tenant final : gpu::ResidencyTenant {
        RendererSession* d = nullptr;
        size_t reclaimableBytes() override { return d->reclaimable(); }
        size_t pendingBytes() override { return d->pending(); }
        size_t release(int gracefulMs) override { return d->releaseNow(gracefulMs); }
        void mayRebuild() override { d->cv.notify_all(); }
    };
    Tenant tenant;
    std::mutex m;
    std::condition_variable cv;
    std::mutex writeM;

    RenderPreviewRequest wanted;
    bool stop = false;
    RenderPreviewState st;
    ReleaseReport report;

    // The process. Started and cleaned up by the control thread; terminated by release() if it must.
    std::unique_ptr<ChildProcess> child;
    std::unique_ptr<PreviewPipe> pipe;
    std::unique_ptr<gpu::TenantRegistration> registration;
    std::thread reader;
    bool connected = false;     // the renderer is warm and takes commands
    bool readerEnded = false;   // the pipe closed
    bool terminatedByRelease = false;
    Clock::time_point restartAfter = Clock::now();  // after a failure: no restart storm
    Clock::time_point refusedSince{};               // when the renderer refused the scene (retry after 10 s)

    // What the renderer was told.
    uint64_t loadedSceneKey = 0;
    bool needsLoad = true;
    int sceneRev = 0, cameraRev = 0;
    uint64_t sentCameraKey = 0;
    int sentWidth = 0, sentHeight = 0;
    bool sentPaused = false, suspendedForLease = false;
    long long nextRequestId = 1;
    std::string lastAck;
    long long lastAckId = 0, lastAckFreed = 0;

    std::thread control;

    // The timeline (preview-timeline.log beside the renderer's log): one line per step with the seconds
    // since the worker began, so a slow update shows where its time went.
    std::mutex noteM;
    std::string timelinePath;
    const Clock::time_point born = Clock::now();
    void note(const std::string& what) {
        std::lock_guard<std::mutex> lk(noteM);
        if (timelinePath.empty()) return;
        std::ofstream f(timelinePath, std::ios::app);
        char t[32];
        std::snprintf(t, sizeof t, "%9.2f s  ", msSince(born) / 1000);
        f << t << what << "\n";
    }
    void setTimeline(const std::string& logPath) {
        std::lock_guard<std::mutex> lk(noteM);
        timelinePath = logPath.empty() ? "" : (std::filesystem::path(logPath).parent_path() / "preview-timeline.log").string();
    }
    bool firstFrameNoted = true;  // the first frame after a load is on the timeline

    // The package export (its own thread: an export takes seconds and must not hold up the control loop).
    std::thread exporter;
    uint64_t exportingKey = 0, exportedKey = 0;
    std::string exportedScene, exportError;
    std::shared_ptr<CancelToken> exportToken;
    int packageRev = 0;

    void exporterLoop();

    size_t reclaimable();
    size_t pending();
    size_t releaseNow(int gracefulMs);
    // A load was sent and the renderer has not answered it yet (its memory is still growing).
    int loadingRev = 0;
    long long loadingForecast = 0;

    void send(const std::string& line) {
        std::lock_guard<std::mutex> lk(writeM);
        if (child) child->writeLine(line);
    }

    void readLoop(PreviewPipe* p) {
        if (p->waitForClient()) {
            PreviewStatus s;
            PreviewFrame f;
            for (;;) {
                const PreviewPipe::Message msg = p->read(s, f);
                if (msg == PreviewPipe::Message::Closed) break;
                std::lock_guard<std::mutex> lk(m);
                if (msg == PreviewPipe::Message::Status) {
                    if (!s.state.empty()) st.phase = s.state;
                    st.message = s.message;
                    st.spp = s.spp;
                    st.msPerSpp = s.msPerSpp;
                    if (s.processBytes) st.processBytes = s.processBytes;
                    st.availableBytes = s.availableBytes;
                    if (s.state != "starting" && !s.state.empty() && !connected) {
                        connected = true;
                        note("renderer ready");
                    }
                    const bool settled = s.state == "rendering" || s.state == "converged" || s.state == "refused" || s.state == "error";
                    if (loadingRev != 0 && settled && s.sceneRev >= loadingRev) {
                        char b[160];
                        std::snprintf(b, sizeof b, "scene %d %s after %.0f ms in the renderer", loadingRev, s.state.c_str(), s.loadMs);
                        note(b + (s.message.empty() ? std::string() : ": " + s.message));
                        loadingRev = 0;
                        firstFrameNoted = false;
                    }
                    if (!s.ack.empty()) {
                        lastAck = s.ack;
                        lastAckId = s.requestId;
                        lastAckFreed = s.freedBytes;
                    }
                } else {
                    st.frame = std::make_shared<PreviewFrame>(f);
                    if (!firstFrameNoted && f.sceneRev == sceneRev) {
                        note("first frame of scene " + std::to_string(f.sceneRev));
                        firstFrameNoted = true;
                    }
                    ++st.frameNumber;
                    st.spp = f.spp;
                }
                cv.notify_all();
            }
        }
        std::lock_guard<std::mutex> lk(m);
        readerEnded = true;
        connected = false;
        cv.notify_all();
    }

    // Whether `ppt` has the serve mode: an older ppt.exe opens its editor for an unknown subcommand, so
    // its help text is read first (once per program).
    std::string timelineFor;  // the logPath the timeline belongs to
    std::string checkedPpt;
    bool pptServes = false;
    bool hasServe(const std::string& ppt, const std::string& logPath, std::string* err) {
        if (ppt == checkedPpt) {
            if (!pptServes && err) *err = ppt + " has no `serve` mode (photon_path_tracer branch feature/live-preview-serve)";
            return pptServes;
        }
        checkedPpt = ppt;
        pptServes = false;
        if (!std::filesystem::exists(ppt)) {
            if (err) *err = "no renderer at " + ppt + " (set KNIT_PPT_EXE)";
            return false;
        }
        const std::string helpLog = logPath + ".help";
        ChildProcess help;
        ChildProcessStart what;
        what.program = ppt;
        what.arguments = {"--help"};
        what.logPath = helpLog;
        if (help.start(what, err) && help.exitedWithin(10000)) {
            std::ifstream text(helpLog);
            const std::string all((std::istreambuf_iterator<char>(text)), std::istreambuf_iterator<char>());
            pptServes = all.find("serve --frames") != std::string::npos;
        }
        if (!pptServes && err) *err = ppt + " has no `serve` mode (photon_path_tracer branch feature/live-preview-serve)";
        return pptServes;
    }

    // Called without `m`.
    bool startProcess(const RenderPreviewRequest& r, std::string* err) {
        std::error_code ec;
        const std::string logPath = r.logPath.empty() ? (std::filesystem::temp_directory_path() / "knit-preview-ppt.log").string() : r.logPath;
        std::filesystem::create_directories(std::filesystem::path(logPath).parent_path(), ec);
        if (!hasServe(r.ppt, logPath, err)) return false;
        auto newPipe = std::make_unique<PreviewPipe>();
        if (!newPipe->create(err)) return false;
        ChildProcessStart what;
        what.program = r.ppt;
        what.arguments = {"serve", "--frames", newPipe->name(), "--width", std::to_string(r.width), "--height", std::to_string(r.height)};
        what.logPath = logPath;
        what.pipeStdin = true;
        auto newChild = std::make_unique<ChildProcess>();
        if (!newChild->start(what, err)) return false;
        auto newRegistration = std::make_unique<gpu::TenantRegistration>(&tenant);
        std::lock_guard<std::mutex> lk(m);
        pipe = std::move(newPipe);
        child = std::move(newChild);
        registration = std::move(newRegistration);
        connected = false;
        readerEnded = false;
        terminatedByRelease = false;
        loadedSceneKey = 0;
        needsLoad = true;
        sentCameraKey = 0;
        sentWidth = r.width;
        sentHeight = r.height;
        sentPaused = false;
        suspendedForLease = false;
        st.phase = "starting";
        st.message = "warming the renderer";
        note("renderer started: " + r.ppt);
        reader = std::thread([this, p = pipe.get()] { readLoop(p); });
        return true;
    }

    // Ends the process (quit, then terminate) and cleans up. Called without `m`.
    void stopProcess(bool graceful) {
        registration.reset();  // no lease calls into us any more; waits for one in flight
        if (child) {
            if (graceful) {
                send("{\"cmd\":\"quit\"}");
                if (!child->exitedWithin(gpu::kGracefulReleaseMs)) child->terminate(kExitConfirmMs);
            } else {
                child->terminate(kExitConfirmMs);
            }
        }
        if (pipe) pipe->unblock();
        if (reader.joinable()) reader.join();
        std::lock_guard<std::mutex> lk(m);
        {
            std::lock_guard<std::mutex> w(writeM);
            child.reset();
        }
        pipe.reset();
        connected = false;
        loadingRev = 0;
        loadedSceneKey = 0;
        needsLoad = true;
    }

    void controlLoop() {
        std::unique_lock<std::mutex> lk(m);
        while (!stop) {
            cv.wait_for(lk, std::chrono::milliseconds(50));
            if (stop) break;
            const RenderPreviewRequest r = wanted;
            if (r.logPath != timelineFor) {
                timelineFor = r.logPath;
                setTimeline(r.logPath);
            }
            // An export for results that are no longer wanted stops.
            if (exportToken && exportingKey != 0 && (!r.enabled || r.packageKey != exportingKey)) exportToken->cancel();

            // The process ended (crash, or terminated by a release): clean up, and say so unless expected.
            if (child && child->exitedWithin(0)) {
                const bool expected = terminatedByRelease;
                const int code = child->exitCode();
                lk.unlock();
                stopProcess(false);
                lk.lock();
                if (!expected) {
                    st.phase = "error";
                    st.message = "the renderer ended (exit code " + std::to_string(code) + ")";
                    restartAfter = Clock::now() + std::chrono::seconds(5);
                } else {
                    st.phase = "released";
                    st.message = "the renderer was stopped for a GPU job";
                }
                continue;
            }
            if (!r.enabled) {
                if (child) {
                    lk.unlock();
                    stopProcess(true);
                    lk.lock();
                }
                st.phase = "off";
                st.message.clear();
                continue;
            }
            const bool jobRuns = gpu::priorityHeld();
            if (!child) {
                if (jobRuns) {
                    st.phase = "waiting";
                    st.message = "a GPU job runs; the preview starts after it";
                    continue;
                }
                if (Clock::now() < restartAfter || r.ppt.empty()) continue;
                lk.unlock();
                std::string err;
                const bool ok = startProcess(r, &err);
                lk.lock();
                if (!ok) {
                    st.phase = "error";
                    st.message = err;
                    restartAfter = Clock::now() + std::chrono::seconds(5);
                }
                continue;
            }
            if (!connected) continue;

            // A GPU job holds priority: no rendering, loading or resuming until it ends.
            if (jobRuns) {
                if (!suspendedForLease && !sentPaused) {
                    send("{\"cmd\":\"suspend\",\"request_id\":" + std::to_string(nextRequestId++) + "}");
                    suspendedForLease = true;
                    note("suspended: a GPU job holds priority");
                }
                st.message = "paused for a GPU job";
                continue;
            }
            if (r.paused != sentPaused) {
                send(r.paused ? "{\"cmd\":\"suspend\",\"request_id\":" + std::to_string(nextRequestId++) + "}" : "{\"cmd\":\"resume\"}");
                sentPaused = r.paused;
            }
            if (r.paused) continue;
            if (suspendedForLease) {
                suspendedForLease = false;
                note("the GPU job ended");
                if (!needsLoad) send("{\"cmd\":\"resume\"}");
            }
            // A refused scene (not enough memory at that moment) is asked for again every 10 s; the last
            // frame stays meanwhile.
            if (st.phase == "refused" && loadedSceneKey != 0) {
                if (refusedSince == Clock::time_point{}) refusedSince = Clock::now();
                else if (Clock::now() - refusedSince > std::chrono::seconds(10)) {
                    needsLoad = true;
                    refusedSince = Clock::time_point{};
                }
            } else {
                refusedSince = Clock::time_point{};
            }
            if (r.width != sentWidth || r.height != sentHeight) {
                send("{\"cmd\":\"resolution\",\"width\":" + std::to_string(r.width) + ",\"height\":" + std::to_string(r.height) + "}");
                sentWidth = r.width;
                sentHeight = r.height;
            }
            // The scene: the request's own package, or the one the exporter wrote for the current results.
            std::string scenePath = r.scenePath;
            uint64_t sceneKey = r.sceneKey;
            if (scenePath.empty() && r.packageKey != 0 && exportedKey == r.packageKey && !exportedScene.empty()) {
                scenePath = exportedScene;
                sceneKey = exportedKey;
            }
            if (!scenePath.empty() && !r.hold && (sceneKey != loadedSceneKey || needsLoad)) {
                send("{\"cmd\":\"load\",\"scene\":" + quotedJson(scenePath) + ",\"scene_rev\":" + std::to_string(++sceneRev) +
                     ",\"forecast_bytes\":" + std::to_string(r.forecastBytes) + "}");
                note("load scene " + std::to_string(sceneRev) + ": " + scenePath);
                loadedSceneKey = sceneKey;
                loadingRev = sceneRev;
                loadingForecast = r.forecastBytes;
                needsLoad = false;
                sentCameraKey = 0;
                st.phase = "loading";
            }
            if (r.cameraKey != 0 && r.cameraKey != sentCameraKey && loadedSceneKey != 0) {
                send("{\"cmd\":\"camera\",\"eye\":" + vec3Json(r.eye) + ",\"target\":" + vec3Json(r.target) + ",\"up\":" + vec3Json(r.up) +
                     ",\"fov_deg\":" + std::to_string(r.fovDeg) + ",\"interacting\":" + (r.interacting ? "true" : "false") +
                     ",\"camera_rev\":" + std::to_string(++cameraRev) + "}");
                sentCameraKey = r.cameraKey;
            }
        }
        lk.unlock();
        stopProcess(true);
    }
};

void RenderPreviewWorker::RendererSession::exporterLoop() {
    std::unique_lock<std::mutex> lk(m);
    while (!stop) {
        cv.wait_for(lk, std::chrono::milliseconds(50));
        if (stop) break;
        const RenderPreviewRequest r = wanted;
        if (!r.enabled || r.packageKey == 0 || r.packageKey == exportedKey || !r.scenePath.empty() || r.packageRoot.empty()) continue;
        if (gpu::priorityHeld()) continue;  // a GPU job first; the export may need the device too
        exportingKey = r.packageKey;
        auto token = std::make_shared<CancelToken>();
        exportToken = token;
        const int rev = ++packageRev;
        const std::filesystem::path dir = std::filesystem::path(r.packageRoot) / ("rev" + std::to_string(rev));
        lk.unlock();
        note("export rev" + std::to_string(rev) + " started");
        const Clock::time_point exportStart = Clock::now();
        ExportRequest e = r.package;
        e.dir = dir.string();
        e.name = "preview";
        e.format = kPackageFormat;
        std::string err;
        bool result = false;
        try {
            ThreadCancelScope scope(token.get());
            const SceneEntry* entry = e.settings ? sceneEntryFor(*e.settings) : nullptr;
            if (!entry) err = "no scene to export";
            else result = exportPreviewPackage(e, *entry, [&](const std::string& text) {
                const std::lock_guard<std::mutex> lock(m);
                st.message = text;
            }, &err);
        } catch (const Cancelled&) {
            err = "cancelled";
        } catch (const std::exception& ex) {
            err = ex.what();
        }
        char took[64];
        std::snprintf(took, sizeof took, " after %.1f s", msSince(exportStart) / 1000);
        note("export rev" + std::to_string(rev) + (result ? std::string(" written") : " failed (" + err + ")") + took);
        // The compiled textures: the earlier revisions' go into the cache (the renderer has loaded them by
        // now), the cache's into the new package.
        const std::filesystem::path root(r.packageRoot);
        for (int earlier = std::max(1, rev - 2); earlier < rev; ++earlier)
            copyCompiledTextures(compiledTexturesOf(root / ("rev" + std::to_string(earlier))), r.textureCache);
        if (result) copyCompiledTextures(r.textureCache, compiledTexturesOf(dir));
        // The packages before the previous one are not loaded any more.
        std::error_code ec;
        if (rev > 2) std::filesystem::remove_all(root / ("rev" + std::to_string(rev - 2)), ec);
        lk.lock();
        exportingKey = 0;
        if (token->cancelled()) continue;
        exportedKey = r.packageKey;  // a failed export is not retried for the same results
        if (result) {
            exportedScene = (dir / "scene.json").string();
            exportError.clear();
        } else {
            exportError = err;
            st.phase = "error";
            st.message = "the preview's export failed: " + err;
        }
    }
}

RenderPreviewWorker::RenderPreviewWorker() : session_(std::make_unique<RendererSession>()) {
    session_->tenant.d = session_.get();
    session_->control = std::thread([this] { session_->controlLoop(); });
    session_->exporter = std::thread([this] { session_->exporterLoop(); });
}

RenderPreviewWorker::~RenderPreviewWorker() {
    {
        std::lock_guard<std::mutex> lk(session_->m);
        session_->stop = true;
    }
    session_->cv.notify_all();
    {
        std::lock_guard<std::mutex> lk(session_->m);
        if (session_->exportToken) session_->exportToken->cancel();
    }
    session_->exporter.join();
    session_->control.join();
}

void RenderPreviewWorker::request(const RenderPreviewRequest& r) {
    std::lock_guard<std::mutex> lk(session_->m);
    session_->wanted = r;
}

RenderPreviewState RenderPreviewWorker::state() const {
    std::lock_guard<std::mutex> lk(session_->m);
    RenderPreviewState s = session_->st;
    if (session_->exportingKey != 0 && s.phase != "off") {
        s.phase = "exporting";
        s.message = "writing the preview's package";
    }
    return s;
}

RenderPreviewWorker::ReleaseReport RenderPreviewWorker::lastRelease() const {
    std::lock_guard<std::mutex> lk(session_->m);
    return session_->report;
}

size_t RenderPreviewWorker::RendererSession::reclaimable() {
    std::lock_guard<std::mutex> lk(m);
    if (!child || readerEnded) return 0;
    return st.processBytes > 0 ? static_cast<size_t>(st.processBytes) : 1;  // 1: alive, size not reported yet
}

size_t RenderPreviewWorker::RendererSession::pending() {
    std::lock_guard<std::mutex> lk(m);
    if (!child || readerEnded || loadingRev == 0) return 0;
    return loadingForecast > 0 ? static_cast<size_t>(loadingForecast) : gpu::kUnknownGrowth;
}

size_t RenderPreviewWorker::RendererSession::releaseNow(int gracefulMs) {
    RendererSession& d = *this;
    const Clock::time_point t0 = Clock::now();
    std::unique_lock<std::mutex> lk(d.m);
    if (!d.child || d.readerEnded) return 0;
    d.needsLoad = true;  // whatever happens, the scene is loaded again after the lease
    const long long heldBefore = d.st.processBytes;
    const long long id = d.nextRequestId++;
    lk.unlock();
    // Graceful: finish the running sample, drop the scene mirror (the Integrator stays warm).
    d.send("{\"cmd\":\"suspend\",\"request_id\":" + std::to_string(d.nextRequestId++) + "}");
    d.send("{\"cmd\":\"release_scene\",\"request_id\":" + std::to_string(id) + "}");
    lk.lock();
    const bool acked = d.cv.wait_until(lk, t0 + std::chrono::milliseconds(gracefulMs),
                                       [&] { return (d.lastAck == "release_scene" && d.lastAckId == id) || d.readerEnded; });
    if (acked && !d.readerEnded) {
        d.report = ReleaseReport{msSince(t0), false, d.lastAckFreed};
        d.note("released for a GPU job in " + std::to_string(static_cast<int>(d.report.ms)) + " ms");
        d.st.phase = "released";
        d.st.message = "the preview gave its GPU memory to a job";
        return static_cast<size_t>(d.lastAckFreed);
    }
    // Not within the timeout (a scene load cannot be interrupted): terminate, and return only when the
    // exit is confirmed -- sending the termination is not proof that the memory is free.
    d.terminatedByRelease = true;
    ChildProcess* process = d.child.get();
    lk.unlock();
    const bool exited = process->terminate(kExitConfirmMs);
    lk.lock();
    d.report = ReleaseReport{msSince(t0), true, exited ? heldBefore : 0};
    d.note("terminated for a GPU job after " + std::to_string(static_cast<int>(d.report.ms)) + " ms" + (exited ? "" : " (no exit)"));
    d.st.phase = "released";
    d.st.message = exited ? "the renderer was stopped for a GPU job" : "the renderer did not exit";
    d.cv.notify_all();
    return exited ? static_cast<size_t>(heldBefore) : 0;
}

} // namespace knit
