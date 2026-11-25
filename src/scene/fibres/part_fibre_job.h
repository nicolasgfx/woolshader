#pragma once
// One part's fibres in three steps, so that every part's candidates can run in one budget order:
// prepare (strands, support, frame), generate (all candidates or one key window), finish (angles,
// checks, digest). Private to the Fibres stage.
#include <memory>
#include <optional>
#include <string>

#include "model/coordinates.h"
#include "fibres/fibre_checks.h"
#include "fibres/fibre_generator.h"
#include "scene/fibres/fibre_stage.h"
#include "scene/appearance/rendered_part.h"
#include "scene/run_log.h"

namespace knit {

namespace gpu {
struct GpuFibreInputs;
}

class SceneSolids;

// What a part's fibres are prepared from.
struct PartFibreSetup {
    const ScenePart& part;
    int index;                  // in Scene::parts
    const SceneSettings& settings;
    const FibreParams& base;    // the settings' parameters, before the part's own
    const SceneSolids& solids;  // fibres hide inside them
    double modelRing;           // the exported tubes' model ring (exportModelRing)
    FibreRunOptions options;
};

// What finish hands back: the fibres, and what the run's log and totals take from them.
struct FinishedFibres {
    FibreSet fibres;
    FibreAngles angles;
    FibreParams params;
    std::string checks, profile, digest, gpuLine;
    double anglesSeconds = 0;  // wall time of measureFibreAngles
};

class PartFibreJob {
public:
    // Null with the reason. On the heap: the tube index's callback refers to the job.
    static std::unique_ptr<PartFibreJob> prepare(const PartFibreSetup& setup, std::string* err);
    ~PartFibreJob();

    // The candidates with float key in [keyLo, keyHi) (budget order); every candidate at (0, 2).
    bool generate(float keyLo, float keyHi, FibreSet& out, std::string* err);
    FinishedFibres finish(FibreSet fibres);

    // Bounded residency (budget mode): free the support while the part waits, and build it again.
    void releaseSupport();
    bool rebuildSupport(std::string* err);
    bool released() const { return released_; }

    const CandidateSource& candidates() const { return source_; }  // valid while released
    const ScenePart& part() const { return part_; }

private:
    explicit PartFibreJob(const PartFibreSetup& setup);
    SurfaceParams surfaceParams() const;
    void buildTubes();

    const ScenePart& part_;
    const SceneSolids& solids_;
    const FibreRunOptions options_;
    const double modelRing_;
    RenderedPart rendered_;
    SceneMapping map_;
    bool compressed_ = false;
    FibreParams params_;
    FibreFrame frame_;
    double serializedMargin_ = 0;  // the float32 envelope of the exported coordinates
    SurfaceModel surface_;
    std::unique_ptr<ExportTubeIndex> tubes_;
    std::shared_ptr<gpu::GpuFibreInputs> gpuInputs_;  // reused by the part's key windows
    std::string profile_, checks_, gpuLine_;          // stage lines, check lines, the generator's run reports
    StageClock partClock_;                            // since prepare
    CandidateSource source_;
    bool released_ = false;
};

// The model ring of the exported tubes (the fibres' root depth follows it): the export's.
double exportModelRing(const SceneSettings& settings);

// The inputs of a part's fibres, one "name = value" per line (saved fibres, fibre_files.h): the
// strands, the export tube, the strand and fibre parameters, the placement and the producer. Empty with
// the reason.
std::optional<std::string> partFibreKey(const ScenePart& part, int index, const SceneSettings& settings, const FibreParams& base,
                                        std::string* err);
// The key's producer line: a generator change that alters fibres without a settings change gets a new one, so
// saved fibres of the former generator are refused (stale: fibre backend).
std::string fibreBackendLine();

} // namespace knit
