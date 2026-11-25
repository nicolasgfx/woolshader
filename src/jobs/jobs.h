#pragma once
// The four jobs the viewer and knit_check start: Knit, the strand mesh rebuild, Fibres and Export.
//
// Each job takes the results the user reviewed (KnitResult, FibreResult), never settings keys: the NCP
// solver is not deterministic, so knitting equal settings again need not give the reviewed knit. A job
// returns its result scene, or nullptr with the reason in err; a cancellation throws Cancelled
// (model/cancel.h checkpoints inside the long loops). scene_builder.h runs the jobs on a worker thread;
// tests and knit_check call them directly. The display data are CPU arrays: the viewer's main thread
// uploads them to GL.
#include <memory>
#include <string>

#include "export/export.h"
#include "jobs/built_scene.h"
#include "jobs/results.h"
#include "model/compute_caps.h"
#include "settings/scene_settings.h"

namespace knit {

// Where a running job reports: its progress line and its partial scenes. The queue's reporter shows them in
// its snapshot; tests record them. The folder's one interface.
class JobReporter {
public:
    virtual ~JobReporter() = default;
    virtual void progress(const std::string& stage, int done, int total) = 0;
    virtual void partial(std::shared_ptr<BuiltScene> scene) = 0;
};

// What a running job may do besides computing: report under its revision. Cancellation is the thread's
// token (model/cancel.h), installed by whoever runs the job.
struct JobControl {
    int revision = 0;
    JobReporter* reporter = nullptr;  // null: nobody listens
};

// The Knit job: the settings' object knitted (any layout), or opened from the scene cache.
struct KnitRequest {
    std::shared_ptr<const SceneSettings> settings;
    std::string cacheRoot;   // persistent scene cache (scene_cache.h); empty: none
    // The user asked to knit again a result shown for viewing only: a cache entry without a
    // restorable computational knit is then knitted again (and replaced) instead of shown again.
    bool knitAgain = false;
    // The viewer's --load-knit (scripts): every spiral piece is taken from this knit cache directory
    // (knit_check --save-knit); the scene cache and the knit kept in memory are bypassed, and a
    // directory that does not fit the settings fails the job (nothing is knitted). "": none.
    std::string loadKnitCache;
};

// The strand mesh rebuild: the yarn tubes of a shown knit at the settings' strand tessellation
// (strand.tube_*). Nothing is knitted; the result keeps the scene's other meshes and the same knit.
struct StrandMeshRequest {
    std::shared_ptr<const SceneSettings> settings;
    std::shared_ptr<const BuiltScene> shown;  // the completed Knit result on screen (with its knit)
    std::string cacheRoot;
    // A yarn of fibres rendered as maps: bake the maps the cache lacks (the Knit step's "Bake the yarn maps"),
    // then show them on the tubes; without it the tubes stay plain until the maps exist.
    bool bakeMaps = false;
};

// The Fibres job: the fibres of the reviewed knit (wool), or the hair grown on the object's
// meshes when the Knit step is skipped (knit null).
struct FibreRequest {
    std::shared_ptr<const SceneSettings> settings;
    std::shared_ptr<const KnitResult> knit;
};

// The Export job: the reviewed knit (null when the Knit step is skipped) and fibre results written
// in a format and published to a folder (publish.h). Under Fur: `fibres` carries the base meshes.
struct ExportRequest {
    std::shared_ptr<const SceneSettings> settings;
    std::shared_ptr<const KnitResult> knit;
    std::shared_ptr<const FibreResult> fibres;  // written with includeFibres; under Fur: always
    std::string dir;                            // written into (created if missing); other files stay
    std::string name;                           // stem of the OBJ / Blender files ("": the folder's name)
    std::string format = kPackageFormat;        // an ExportFormat::id
    bool includeFibres = false;
};

// Knit (knit_job.cpp): load and knit (knitScene: any layout), or read the scene cache; the display arrays and
// a KnitResult. Partial scenes go to the reporter while it runs.
std::shared_ptr<BuiltScene> knitObject(const KnitRequest& r, const JobControl& control, std::string* err);
// The strand mesh rebuild (strand_mesh_job.cpp): the shown knit's tubes at the new tessellation.
std::shared_ptr<BuiltScene> rebuildStrandMesh(const StrandMeshRequest& r, const JobControl& control, std::string* err);
// Fibres (fibre_job.cpp): a FibreResult and its display lines.
std::shared_ptr<BuiltScene> growFibres(const FibreRequest& r, const JobControl& control, std::string* err);
// Export (export_job.cpp): after exportEligibility (Blocked: refused; Diagnostic: written, and marked so in
// report.json), through the publication contract.
std::shared_ptr<BuiltScene> exportResults(const ExportRequest& r, const JobControl& control, std::string* err);

// What this machine computes with, detected once (gpu::detectComputeCaps).
const ComputeCaps& computeCaps();

} // namespace knit
