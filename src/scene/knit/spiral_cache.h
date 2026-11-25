#pragma once
// The per-piece knit cache of a run (knit_check --save-knit / --load-knit): files keyed by every
// knitting input, under the producer's subdirectory (<dir>/gpu-f32/ncp). Private to the Knit stage.
#include <optional>
#include <string>

#include "settings/settings_mapping.h"
#include "scene/source/sewn_pieces.h"

namespace knit {

struct KnitRun;
struct PieceKnitResult;
struct PieceWork;

struct SpiralKnitCache {
    std::string saveDir, loadDir, key;  // no directories: no cache, and no key
};
// The directories with the producer's subdirectory, and the key. Empty with the reason if an input
// cannot be hashed.
std::optional<SpiralKnitCache> spiralKnitCache(const SceneSettings& settings, const SpiralSource& src, const std::string& saveDir, const std::string& loadDir,
                                               std::string* err);
// The piece's knitting result from the cache, refused (an error of the run, with the reason) if any
// knitting input changed; gate (a) rerun on the cached curves must reproduce the knitting run's report.
bool loadPieceFromCache(const KnitRun& run, PieceWork& piece, const DiameterPlan& plan, PieceKnitResult& out, std::string* err);
// The knitting result into the knit cache, if gate (a) of the final curves reproduces it (the cache
// reruns it). `knitLog` is the knitting's part of the piece's log.
void savePieceToCache(const KnitRun& run, PieceWork& piece, const PieceKnitResult& out, const std::string& knitLog, double msKnit);

} // namespace knit
