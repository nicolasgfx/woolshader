#pragma once
// The Fibres step's reference samples: what the settings grow on a small piece, in about a second,
// instead of on the whole model (which "Apply to model" computes). The viewer draws them
// (jobs/previews.h).
//   Wool: a flat, square patch (5 wales, as many courses as fill the same length) of the current
//     pattern, yarn diameter, gauge and depth (generateStockinette), with the model's own strands
//     (renderedStrands: D_r, flattening, ply tubes), its twisted-ply appearance as a profile, and
//     fibres from the model's generator and settings, kept where they root in the central square,
//     60 % of the side (no patch-edge artefacts).
//   Fur: a sphere cap (radius max(10, 1.5 x longest hair) mm) under -z gravity, grown with the
//     model's fur generator and settings, thinned to about 4000 hairs.
// Approximate by construction: flat (wool) / a sphere (fur), no hidden test, no export rounding or
// root check.
#include <memory>
#include <optional>
#include <string>

#include "fibres/fibre.h"
#include "knitting/tube/ply_profile.h"
#include "knitting/tube/yarn_surface.h"
#include "model/mesh.h"
#include "settings/scene_settings.h"

namespace knit {

struct WoolSample {
    SurfaceModel surface;                      // the drawn yarn (the ply tubes' envelope when plied)
    std::unique_ptr<PlyShadingProfile> relief; // the twisted-ply appearance of a single tube; null: none
    int reliefPlies = 1;
    bool plyTubes = false;
    FibreSet along, stray;                     // the fibres rooted in the central square, by population
    long long candidates = 0;
    double square = 0;                         // the central square's side [mm]
};
// Empty with the reason.
std::optional<WoolSample> growWoolSample(const SceneSettings& s, std::string* err);

struct FurSample {
    TriMesh sphere;
    double radius = 0;   // [mm]
    double thinning = 1; // the share of the coat's density grown (< 1: thinned)
    FibreSet hairs;
    long long sampled = 0, stuck = 0;
};
FurSample growFurSample(const SceneSettings& s);

} // namespace knit
