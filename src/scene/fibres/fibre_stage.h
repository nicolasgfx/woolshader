#pragma once
// growSceneFibres: the Fibres stage: wool fibres on a knit, or fur on its meshes.
//
// From a Scene and the settings; makes the SceneFibres (scene_fibres.h) the export writes as it is.
// Called by the Fibres job (jobs/fibre_job.cpp) and knit_check. An empty result comes with the reason
// for a failed part, never a CPU stand-in for a device failure.
#include <functional>
#include <optional>
#include <string>

#include "model/parallel.h"
#include "scene/fibres/scene_fibres.h"
#include "scene/scene.h"

namespace knit {

// How the fibres are grown: diagnostics that never change them, progress and memory bounds.
struct FibreRunOptions {
    // The fibre diagnostics (knit_check --fiber-check), none of which changes the fibres (their strands,
    // keys and candidate ids): the attachment and reference checks per part, the measured emergence and
    // along angles (measureFibreAngles) and, on the GPU, the CPU audit of the omitted serialized root
    // check (counted, never applied). Off, as the viewer grows fibres: the light checks only.
    bool fibreCheck = false;
    bool stageProfile = false, fibreDigest = false;  // per-part profile lines, digests
    ExecPolicy inner{true};        // fibre generation and checks inside a part (identical results)
    double fibreSupportMB = 3000;  // bounded fibre-support residency in budget mode [MB private]; 0 = unbounded
    // Progress: (stage, parts done, parts total), per part.
    std::function<void(const std::string& stage, int done, int total)> onPart;
};

// The scene's fibres for its fibre profile (fibres.profile), as the viewer reviews them and the
// export writes them: wool fibres on the knitted parts, or fur on the parts' meshes (fur_stage.h).
std::optional<SceneFibres> growSceneFibres(const SceneEntry& entry, const Scene& scene, const SceneSettings& settings, const FibreRunOptions& options,
                                           std::string* err);
// The wool fibres: every knitted part's candidates in one budget order (fibres.generation =
// budget), or every candidate.
std::optional<SceneFibres> growWoolFibres(const Scene& scene, const SceneSettings& settings, const FibreRunOptions& options, std::string* err);

} // namespace knit
