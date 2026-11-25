#pragma once
// The persistent cache of what the scene viewer shows (a BuiltScene of kind
// Knit), so a scene opens without knitting again. One entry per key under
// out/cache/scenes/<scene>-<id> (git-ignored):
//   * key: the scene, the Knit and Surface stage keys (scene_settings.h) and the asset
//     hash, a cache format version and a generator version (bumped when the pipeline's
//     output changes for the same settings). Material and display settings are not in it:
//     the viewer applies them when drawing;
//   * content: display meshes, part rows, report, bounds, yarn
//     length; not the texture preview (it is re-read from the map set cache). Tiled scenes
//     also keep the computational knit (knit.bin, saved_knit.h), so a knit opened from
//     the cache feeds fibres and the export without knitting again;
//   * written through publish.h (staging + completion marker with checksums); a hit is a
//     verified publication of the key, so a corrupt or partial entry is a miss.
#include <cstdint>
#include <memory>
#include <string>

#include "jobs/built_scene.h"

namespace knit {

std::string sceneCacheRoot();  // <repo>/out/cache/scenes
std::string sceneCacheKey(const SceneSettings& s, uint64_t assetHash);

// The binary form (versioned, little-endian); read fails on any truncation or mismatch.
std::string serializeBuiltScene(const BuiltScene& b);
bool deserializeBuiltScene(const std::string& bytes, BuiltScene& b);

// knitBytes (optional): the tiled knit (serializeKnit), stored as knit.bin in the entry.
bool storeSceneCache(const std::string& root, const SceneSettings& s, uint64_t assetHash, const BuiltScene& b, std::string* err,
                     const std::string* knitBytes = nullptr);
// nullptr on a miss.
std::shared_ptr<BuiltScene> loadSceneCache(const std::string& root, const SceneSettings& s, uint64_t assetHash);
// The entry's knit.bin (false: no such file in a verified entry).
bool loadSceneCacheKnit(const std::string& root, const SceneSettings& s, uint64_t assetHash, std::string& bytes);
// The directory of the entry (whether or not it exists).
std::string sceneCacheEntryDir(const std::string& root, const SceneSettings& s, uint64_t assetHash);

} // namespace knit
