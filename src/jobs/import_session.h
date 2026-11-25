#pragma once
// The Import step's background work (GL-free, tested without a window): read and analyse a file
// (mesh_read.h), make its display preview, then prepare and store the chosen options (import_pipeline.h).
// One file at a time on one worker thread; open() of another file supersedes the current one.
// Every request (open, prepare) carries a request number; cancel() and reset() retire it AT ONCE
// (the phase changes before they return: a cancelled preparation is back at Analysed, a cancelled
// read is Cancelled, reset is Idle) and the worker publishes progress and results only while its
// request is still the current one, so a late completion never overwrites a newer state.
#include <condition_variable>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "model/cancel.h"
#include "import/import_analysis.h"
#include "import/prepare_import.h"
#include "jobs/built_scene.h"

namespace knit {

class ImportSession {
public:
    enum class Phase { Idle, Reading, Analysed, Failed, Preparing, Prepared, Cancelled };

    // All the session says, read at one moment (as SceneBuilder's JobSnapshot).
    struct Snapshot {
        Phase phase = Phase::Idle;
        std::string path, stage, error;
        double fraction = 0;
        int generation = 0;  // increases with every completed phase (the viewer polls it)
        // Analysed and later: the file, its analysis, its decoded textures' warnings and a display
        // preview (in FILE coordinates and axes: the viewer applies the unit / up-axis transform it shows).
        std::shared_ptr<const ImportedMesh> mesh;
        std::shared_ptr<const ImportAnalysis> analysis;
        std::shared_ptr<const BuiltScene> preview;
        std::vector<std::string> textureWarnings;
        // Prepared: the stored entry.
        uint64_t identity = 0;
        ImportOptions options;
        std::vector<std::string> changes;
    };

    ImportSession();
    ~ImportSession();  // cancels and joins
    void open(const std::string& path, const ImportLimits& lim = {});
    // After Analysed: prepare with o and store it in the import cache (root: importCacheRoot()).
    void prepare(const ImportOptions& o);
    void cancel();
    void reset();      // back to Idle (a closed import)
    bool busy() const;  // a task queued or still running (also one already retired)
    // Tests: called on the worker outside the lock at "start" (task taken, nothing run yet) and
    // "finish" (work done, nothing published yet), so a test can interleave deterministically.
    void setTestHook(std::function<void(const char* point)> hook);
    Snapshot snapshot() const;
    // Analysed: whether this file prepared with o is the stored import named by storedIdentity (settings
    // import.identity, hex) and that import is still in the cache, so it need not be prepared again.
    // Hashes the file's mesh (importIdentity): call it once per analysed file.
    bool storedImportMatches(const ImportOptions& o, const std::string& storedIdentity) const;
    // The settings of the prepared import (scene = import, import.* keys) on top of `base`.
    SceneSettings applyTo(const SceneSettings& base) const;

private:
    struct Task;
    void loop();
    void dropFile();  // under the lock: nothing of this file can mix with the next
    mutable std::mutex m_;
    std::condition_variable cv_;
    std::unique_ptr<Task> pending_;
    CancelToken* running_ = nullptr;
    uint64_t request_ = 0;  // the current request; a task with another number publishes nothing
    std::function<void(const char*)> hook_;
    Snapshot state_;
    ImportLimits limits_;  // of the open file: its preparation reads its textures under them
    bool stop_ = false;
    std::thread worker_;
};

// import.repairs (settings) <-> ImportOptions's repair flags.
std::string repairsText(const ImportOptions& o);
void parseRepairs(const std::string& text, ImportOptions& o);
// The ImportOptions stored in settings (import.unit_mm, import.up, import.repairs, import.surface,
// import.wrap_stitch_mm, import.plain_nodes).
ImportOptions importOptionsFrom(const SceneSettings& s);

} // namespace knit
