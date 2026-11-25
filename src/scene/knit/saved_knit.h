#pragma once
// The computational knit of a tiled scene, kept with its scene-cache entry (scene_cache.h, file
// knit.bin), so a knit restored from the cache feeds the fibres and the export WITHOUT knitting
// again: what the user reviewed is what later stages use (the NCP solver is not deterministic,
// so knitting the same settings again need not give the same yarn).
//
// Stored per part: its identity (name, primitive, component, centre, the digest of its mesh) and,
// knitted, the curves and stitch records bit for bit (knit_cache.h's encoding), the pattern
// parameters, the kept diameter and note, the verdicts and check summaries, and the digest of
// gate (a)'s report with the validator that made it. The atlas (UV tiling) and the mesh surface are
// rebuilt from the part's mesh (deterministic; the mesh digest must match). Restored parts are
// UNVERIFIED: verifyKnit reruns gate (a) on the stored curves with the same validator (the
// knitting's numerics, or the CPU validator of a spiral piece) and requires the stored digest (a
// mismatch, a missing GPU, or any other difference refuses the cache).
#include <string>

#include "scene/scene.h"

namespace knit {

// v2: an FNV-1a 64 checksum of everything before it closes the record; a damaged byte anywhere
// (also in report text that no verification reads) refuses the cache before parsing. v3: the
// gate (a) validator per part, and the layout (spiral pieces are cached like every scene).
// A v2 record is still read: it lacks only the validator flag, and every v2 knit was validated by
// the knitting's numerics (spiral pieces were not stored in this format before v3). v4: the gap
// patches (knit.coverage = overlap), parts made after the first knitting pass that loading the scene
// does not make: each one's seed, so it is made again before the parts are compared (v2 / v3
// records have none).
constexpr int kSavedKnitVersion = 4;

std::string serializeKnit(const Scene& t);
// Loads the scene's parts (loadPartsToKnit, any layout), makes the gap patches again from their stored
// seeds (no coverage is measured), and fills every part from `bytes`.
bool restoreKnit(const std::string& bytes, const SceneEntry& e, const SceneSettings& s, Scene& out, std::string* err);
// Gate (a) on the restored curves, compared with the stored digests; fills each knit's `after`.
bool verifyKnit(Scene& t, const SceneSettings& s, std::string* err);

} // namespace knit
