#pragma once
// What every part of one knitting run shares, and the three operations a part knits by: knitPart
// (knit_stage.cpp) chooses one by the part's layout. Private to the Knit stage.
#include "stitchmesh/yarn_settle.h"
#include <string>

#include "knitting/knit/mesh_knit.h"
#include "knitting/knit/piece_knit.h"
#include "scene/knit/knit_stage.h"
#include "scene/source/guided_source.h"
#include "scene/source/scene_solids.h"
#include "scene/source/sewn_pieces.h"
#include "scene/knit/spiral_cache.h"

namespace knit {

struct KnitRun {
    const SceneSettings& settings;
    const KnitRunOptions& options;
    StitchPattern pattern;
    const SceneSolids& solids;     // empty when every part is a spiral
    const SpiralSource& pieces;    // the teddy's manifest (its knit cache key); empty for other scenes
    const std::vector<PlacedPiece>& spirals;  // per part, in part order: the spiral parts' pieces (the pole tests)
    const SourceScene& source;     // the parts' primitives (the names in a piece's log)
    const SpiralKnitCache& cache;  // no directories: no cache
    KnitNumerics& numerics;        // the device's
    const GuidedInputs* guided;    // the guided route's file, markings and piece (null on the other routes)
};

// One part on its UV charts (mesh_part.cpp): every diameter of partDiameters until gate (a) passes, the
// last knitted attempt kept, gate (d) once on it. Fills the part's knit, yarnD and note; a part that
// fits no stitch keeps the reason in its note, one refused at the device budget is stopped
// (KnitStop::DeviceBudget). Any other device failure throws, and the run reports it.
void knitMeshPart(const KnitRun& run, ScenePart& part);

// Spiral part `index` (spiral_part.cpp): knitted at the scene's diameters until gate (a) passes
// (keeping the requested diameter's result, labelled residue, when none passes), or taken from the knit
// cache; GPU curves audited by gate (a) on the CPU; the knit cache entry written; gate (d) when
// requested. Sets the part's knit, yarnD and note; its log goes to stdout as one block. A piece that
// cannot be knitted stays unknitted with the reason in its note (partIncomplete). False for an error of
// the run (a refused knit cache, an invalid final surface), with the reason.
bool knitSpiralPart(const KnitRun& run, size_t index, ScenePart& part, std::string* err);

// The guided part (guided_part.cpp): its stitch mesh built from the markings, its yarn embedded, converted and settled
// to the scene's diameter on the device (stitchmesh/yarn_settle.h), with stitch records derived from the stitch mesh.
// Settled and not settled parts keep the retained curves as their knit; a refused construction or a stopped solve leaves
// the part unknitted (a stopped one KnitStop::DeviceBudget or TimeCap) with the generated yarn as a labelled preview.
// Sets part.guided in every case. Throws Cancelled; a device failure other than the budget throws.
void knitGuidedPart(const KnitRun& run, ScenePart& part);
// The guided route's settle settings for a scene: the stitch-mesh defaults, Final's rule for a failed construction, and
// Draft converting without a solve. Their solver options decide the run's numerics (knit_stage.cpp) and identity line.
SettleSettings guidedSettleSettings(const SceneSettings& s);
// The guided part's identity line (GuidedPart::inputs): source and piece digests, the markings file, preparation, gauge,
// D, c, the guided version and the solver identity with its options, always in full.
std::string guidedInputsText(const SceneSettings& s, const GuidedInputs& in, double D);

} // namespace knit
