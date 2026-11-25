#pragma once
// Saved fibres (knit_check --save-fibres / --load-fibres): a wool fibre result as one file per
// knitted part, so a later run on the same knit exports without growing the fibres again. Loading
// gives the same SceneFibres as growSceneFibres (fibre_stage.h); nothing downstream can tell
// them apart. The files are fibre-cache v6, read and written here; older versions are refused.
//
// A file is refused, never used, if its part's strands, fibre parameters, tube, placement or
// producer differ from the current ones, or if its fibres do not reproduce their stored digest. The
// stored check results are not read (the shared result has no place for them): diagnostics such as
// --fiber-check need grown fibres. A budget-ordered file holds only the candidates below its
// coverage; a budget that needs more is refused, not resumed: grow the fibres again.
#include <cstdint>
#include <optional>
#include <string>

#include "settings/scene_registry.h"
#include "settings/scene_settings.h"
#include "scene/scene.h"
#include "scene/fibres/scene_fibres.h"

namespace knit {

// What identifies one part's file besides its fibres.
struct FibreFileIdentity {
    std::string key;           // the fibres' inputs, one "name = value" per line
    float keyBound = 2.0f;     // every candidate with float key below it was evaluated (2: all)
    uint64_t fibreDigest = 0;  // fibreSnapshotChunkDigest of the fibres
};
// One part's file. What else a file holds (the candidate source of a budget run, the check results
// and their digest) is skipped and written empty.
struct FibreCacheEntry {
    FibreFileIdentity identity;
    FibreSet fs;
};
bool writeFibreCache(const std::string& path, const FibreFileIdentity& identity, const FibreSet& fs, std::string* err);
std::optional<FibreCacheEntry> readFibreCache(const std::string& path, std::string* err);
// The name of the first line of `want` that the stored key `have` does not hold ("" if it holds them
// all): value by value, numbers to rounding, and a stored line may continue past `want`'s.
std::string staleFibreKeyField(const std::string& have, const std::string& want);

// Writes every part of f to <dir>/<part name>.fibrecache.
bool saveSceneFibres(const std::string& dir, const Scene& t, const SceneSettings& s, const SceneFibres& f, std::string* err);
// Reads every knitted part of t from dir and decides the segment budget's threshold as a growth does.
bool loadSceneFibres(const std::string& dir, const Scene& t, const SceneSettings& s, SceneFibres& out, std::string* err);

} // namespace knit
