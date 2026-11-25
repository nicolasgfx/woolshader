// The job queue: one worker thread, the newest request wins, each finished job publishes one scene.
#include "jobs/scene_builder.h"

#include <chrono>
#include <memory>

#include "gpu/device.h"
#include "model/cancel.h"
#include "model/log.h"

namespace knit {

// The running job's reports: progress, and partial scenes (published only while that job is
// still the running one).
class SceneBuilder::Reporter : public JobReporter {
public:
    Reporter(SceneBuilder& b, int revision) : b_(b), revision_(revision) {}
    void progress(const std::string& stage, int done, int total) override {
        std::lock_guard<std::mutex> lk(b_.m_);
        if (b_.runningRevision_ != revision_) return;
        b_.progress_.stage = stage;
        b_.progress_.done = done;
        b_.progress_.total = total;
    }
    void partial(std::shared_ptr<BuiltScene> scene) override {
        scene->partial = true;
        scene->knit.reset();  // a partial scene is never a result
        std::lock_guard<std::mutex> lk(b_.m_);
        if (b_.runningRevision_ != revision_ || (b_.running_ && b_.running_->cancelled())) return;
        b_.latest_ = std::move(scene);
    }

private:
    SceneBuilder& b_;
    int revision_;
};

JobSnapshot SceneBuilder::snapshot() const {
    std::lock_guard<std::mutex> lk(m_);
    JobSnapshot j;
    j.busy = running_ != nullptr || pending_ != nullptr;
    j.cancelRequested = running_ && cancelRequested_;
    j.progress = progress_;
    j.progress.revision = runningRevision_;
    j.progress.seconds = runningRevision_ >= 0 ? secondsSince(jobStart_) : 0;
    j.latest = latest_;
    j.latestComplete = complete_;
    j.lastOutcome = outcome_;
    j.lastError = error_;
    j.published = published_;
    return j;
}

SceneBuilder::SceneBuilder(JobFn fn) : fn_(std::move(fn)) { worker_ = std::thread([this] { loop(); }); }

SceneBuilder::~SceneBuilder() {
    {
        std::lock_guard<std::mutex> lk(m_);
        stop_ = true;
        pending_.reset();
        if (running_) running_->cancel();
    }
    cv_.notify_all();
    worker_.join();
}

int SceneBuilder::queue(JobRequest r) {
    int revision = 0;
    {
        std::lock_guard<std::mutex> lk(m_);
        revision = nextRevision_++;
        pending_ = std::make_unique<QueuedJob>(QueuedJob{std::move(r), revision});
        if (running_ && runningKind_ != BuiltScene::Kind::Export) running_->cancel();
    }
    cv_.notify_all();
    return revision;
}

int SceneBuilder::startKnit(KnitRequest r) { return queue(JobRequest(std::move(r))); }
int SceneBuilder::startStrandMesh(StrandMeshRequest r) { return queue(JobRequest(std::move(r))); }
int SceneBuilder::startFibres(FibreRequest r) { return queue(JobRequest(std::move(r))); }
int SceneBuilder::startExport(ExportRequest r) { return queue(JobRequest(std::move(r))); }

void SceneBuilder::cancel() {
    std::lock_guard<std::mutex> lk(m_);
    pending_.reset();
    if (running_) {
        running_->cancel();
        cancelRequested_ = true;
    }
}

void SceneBuilder::loop() {
    for (;;) {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return stop_ || pending_; });
        if (stop_) return;
        const QueuedJob job = std::move(*pending_);
        pending_.reset();
        CancelToken token;
        running_ = &token;
        cancelRequested_ = false;
        runningKind_ = jobKind(job.request);
        runningRevision_ = job.revision;
        progress_ = JobProgress{};
        progress_.stage = "starting";
        jobStart_ = std::chrono::steady_clock::now();
        Reporter reporter(*this, job.revision);
        const JobControl control{job.revision, &reporter};
        lk.unlock();
        // A job on another scene: the meshes the previous scene left on the device go (gpu::releaseResident).
        const SceneSettings& settings = jobSettings(job.request);
        const std::string scene = settings.scene + "\n" + settings.importSource;
        if (scene != lastScene_) gpu::releaseResident();
        lastScene_ = scene;
        std::shared_ptr<BuiltScene> res;
        std::string err;
        JobOutcome oc = JobOutcome::Failed;
        try {
            ThreadCancelScope scope(&token);
            // The job has priority over the live preview for its whole run (gpu/device.h): the preview gives
            // its memory back when the job does not fit, and loads or rebuilds nothing until the job ends.
            const gpu::PriorityLease lease(gpu::kJobBudget);
            res = fn_ ? fn_(job.request, control, &err) : runJob(job.request, control, &err);
            oc = res ? JobOutcome::Done : JobOutcome::Failed;
        } catch (const Cancelled&) {
            oc = JobOutcome::Cancelled;
        } catch (const std::exception& e) {
            err = e.what();
        }
        lk.lock();
        running_ = nullptr;
        runningRevision_ = -1;
        cancelRequested_ = false;
        outcome_ = oc;
        error_ = oc == JobOutcome::Failed ? (err.empty() ? std::string("failed") : err) : std::string();
        if (oc == JobOutcome::Done) {
            latest_ = res;
            complete_ = std::move(res);
            ++published_;
        }
    }
}

BuiltScene::Kind jobKind(const JobRequest& r) {
    if (std::holds_alternative<FibreRequest>(r)) return BuiltScene::Kind::Fibres;
    if (std::holds_alternative<ExportRequest>(r)) return BuiltScene::Kind::Export;
    return BuiltScene::Kind::Knit;  // a knit or a strand mesh rebuild
}

const SceneSettings& jobSettings(const JobRequest& r) {
    if (const KnitRequest* k = std::get_if<KnitRequest>(&r)) return *k->settings;
    if (const StrandMeshRequest* p = std::get_if<StrandMeshRequest>(&r)) return *p->settings;
    if (const FibreRequest* f = std::get_if<FibreRequest>(&r)) return *f->settings;
    return *std::get<ExportRequest>(r).settings;
}

std::shared_ptr<BuiltScene> runJob(const JobRequest& r, const JobControl& control, std::string* err) {
    if (const KnitRequest* k = std::get_if<KnitRequest>(&r)) return knitObject(*k, control, err);
    if (const StrandMeshRequest* p = std::get_if<StrandMeshRequest>(&r)) return rebuildStrandMesh(*p, control, err);
    if (const FibreRequest* f = std::get_if<FibreRequest>(&r)) return growFibres(*f, control, err);
    return exportResults(std::get<ExportRequest>(r), control, err);
}

} // namespace knit
