#pragma once
// The keys that say whether a stage's result is current: a stage is rebuilt iff its key changes, and
// the key covers every upstream input and the source asset's hash.
#include <cstdint>
#include <string>

#include "settings/scene_settings.h"

namespace knit {

// What `stage` computes from, one "key=value" line per input (the canonical description the stage
// keys hash): the keys of the stage and of the stages before it that can change its result under
// these settings, and the values derived from them. A feature's parameters count only while it is
// in use; a label of other values never does. Settings with the same text build the same result.
std::string stageInputs(const SceneSettings& s, Stage stage);
// Stage key: hash of stageInputs, the asset hash and the generator version.
uint64_t stageKey(const SceneSettings& s, Stage stage, uint64_t assetHash);
// The Fur result's key: the Load keys (workflow among them), the asset hash, every pile.* key and the segment
// budget - never a knitting key (fur grows on the mesh, not on a knit) nor a material (export only).
uint64_t furKey(const SceneSettings& s, uint64_t assetHash);
// The geometry a knit result depends on: the Surface stage key (with the Load and Knit inputs). A
// result is current iff its key equals the current settings' key.
uint64_t geometryKey(const SceneSettings& s);
// The fibres a fibre result depends on: every Fibres-stage key and the segment budget (the budget
// order generates only what the budget keeps), whatever fibres.enabled says. Fur: furKey (the
// mesh and the pile keys; no knitting key).
uint64_t fibreResultKey(const SceneSettings& s);

} // namespace knit
