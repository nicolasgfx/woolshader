#pragma once
// The named settings profiles: Draft and Final, the viewer's opening settings, and what this machine
// cannot compute.
//
// Each function changes or reads a SceneSettings; the caps come from the caller (jobs/jobs.h
// computeCaps). Called by the viewer, knit_check, the jobs and exportEligibility. None of them can fail: computeNotes turns a
// missing device or solver into a sentence for the user.
#include <string>
#include <vector>

#include "model/compute_caps.h"
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"

namespace knit {

// The viewer's interactive profile: no settle, no certified repair, no diameter fallback and
// no gate (d) (the knit, gate a and the tiling checks only); "full validation" restores them. The CLI keeps the full
// profile (teddy parity).
// The full profile tries one diameter
// (knit.diameter_attempts=1; the result is labelled residue when it fails); the preview has no
// diameter search.
void applyPreviewProfile(SceneSettings& s);
void applyFullProfile(SceneSettings& s, const SceneEntry& e);
bool isPreviewProfile(const SceneSettings& s);
// Final quality: EVERY step of the full profile enabled as applyFullProfile sets it (settle,
// repair, the scene's diameter fallback, the surface gate and its half step). A settings file that
// enables only some of them is neither Draft nor Final ("custom") and is not exportable. The one
// predicate of the UI, the knit result (computed and restored from the cache) and the export rule.
bool isFinalProfile(const SceneSettings& s, const SceneEntry& e);
// What cannot run on this machine, for the user (empty: everything can): no usable GPU, or a Final
// knit without the solver's linear algebra (Draft needs no relaxation). The pipelines refuse the job
// with the same reason; the viewer says it before a job.
std::vector<std::string> computeNotes(const SceneSettings& s, const ComputeCaps& c);
// The scene's own fibre density and islands (SceneEntry::defaultFibreDensity / defaultFibreIslands),
// over what s holds: the viewer's opening settings and a scene switch.
void applySceneFibreDefaults(SceneSettings& s, const SceneEntry& e);
// The settings the viewer opens a scene with (what knit_check --scene <id> --cache-build
// precomputes): the scene defaults in the preview profile, on the GPU when available.
SceneSettings viewerDefaults(const SceneEntry& e);

} // namespace knit
