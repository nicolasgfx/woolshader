#pragma once
// The job queue (an active object): one worker thread runs the jobs of jobs.h, GL-free (tested without a window).
//
// The queue numbers the requests (revisions):
//   - the newest pending request replaces an older pending one (latest wins);
//   - a running Knit, strand mesh or Fibres job is cancelled by a newer request; an Export job runs to
//     completion from its snapshot unless cancel() or the destructor stops it;
//   - a finished job publishes ONE immutable BuiltScene (geometry, settings, verdicts and its result object
//     together), assigned under the queue's lock; a cancelled or failed job publishes nothing
//     (JobSnapshot::latestComplete keeps the previous result; lastOutcome says what happened);
//   - a long Knit job may publish PARTIAL scenes on the way (BuiltScene::partial: the stand-ins first, then the
//     parts finished so far) and reports progress; a partial scene is never a result (no KnitResult);
//   - a Knit job whose result is in the persistent cache (scene_cache.h) is not computed.
#include <chrono>
#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <variant>

#include "jobs/jobs.h"
#include "model/cancel.h"

namespace knit {

// A request queued at the builder: exactly one of the four.
using JobRequest = std::variant<KnitRequest, StrandMeshRequest, FibreRequest, ExportRequest>;
// The kind of result a request publishes (a strand mesh rebuild: Knit) and its settings snapshot.
BuiltScene::Kind jobKind(const JobRequest& r);
const SceneSettings& jobSettings(const JobRequest& r);
// The request's job function (knitObject, rebuildStrandMesh, growFibres or exportResults).
std::shared_ptr<BuiltScene> runJob(const JobRequest& r, const JobControl& control, std::string* err);

struct JobProgress {
    int revision = -1;   // -1: idle
    std::string stage;
    int done = 0, total = 0;
    double seconds = 0;  // since the job started
};

enum class JobOutcome { None, Done, Cancelled, Failed };

// One coherent read of the builder (taken under its lock: the fields agree with each other).
struct JobSnapshot {
    bool busy = false;             // a job runs or waits
    bool cancelRequested = false;  // a cancel is pending on the running job
    JobProgress progress;          // revision -1: idle
    std::shared_ptr<const BuiltScene> latest;          // newest published, possibly partial
    std::shared_ptr<const BuiltScene> latestComplete;  // newest completed job's scene
    JobOutcome lastOutcome = JobOutcome::None;
    std::string lastError;
    int published = 0;             // completed jobs
};

class SceneBuilder {
public:
    // The job function (tests inject synthetic jobs; empty: the real ones, runJob).
    using JobFn = std::function<std::shared_ptr<BuiltScene>(const JobRequest&, const JobControl&, std::string* err)>;
    explicit SceneBuilder(JobFn fn = {});
    ~SceneBuilder();  // cancels and joins
    // Each queues its request (replacing a pending one, cancelling a running Knit, strand mesh or
    // Fibres job) and returns its revision.
    int startKnit(KnitRequest r);
    int startStrandMesh(StrandMeshRequest r);
    int startFibres(FibreRequest r);
    int startExport(ExportRequest r);
    // Cancels the running job (any kind) and drops the pending one; the worker stops at its next
    // checkpoint (JobSnapshot::cancelRequested until then). Nothing restarts by itself.
    void cancel();
    JobSnapshot snapshot() const;

private:
    class Reporter;
    struct QueuedJob {
        JobRequest request;
        int revision = 0;
    };
    int queue(JobRequest r);
    void loop();
    JobFn fn_;
    mutable std::mutex m_;
    std::condition_variable cv_;
    std::unique_ptr<QueuedJob> pending_;
    int nextRevision_ = 1;
    std::shared_ptr<const BuiltScene> latest_, complete_;
    CancelToken* running_ = nullptr;        // token of the running job (under m_)
    bool cancelRequested_ = false;
    BuiltScene::Kind runningKind_ = BuiltScene::Kind::Knit;  // queue() does not cancel a running Export
    int runningRevision_ = -1;
    JobOutcome outcome_ = JobOutcome::None;
    std::string error_;
    int published_ = 0;
    JobProgress progress_;
    std::chrono::steady_clock::time_point jobStart_;
    std::string lastScene_;  // the scene of the last job (worker thread only)
    bool stop_ = false;
    std::thread worker_;
};

} // namespace knit
