#pragma once
// The Fur stage (workflow = fur): a loop pile and fine fibres grown on a scene's parts themselves, with no
// knitting. The parts are the Scene's components loaded without a knit (loadBaseParts); the result is a
// SceneFibres, so the viewer and the export use the wool path's strands, keys and writer. Per part:
//   - the root caches of both layers (fur_field.h sampleFurRoots): the hidden test drops a root inside another
//     eligible solid; parts listed in pile.excluded_parts grow nothing;
//   - the field (buildFurField): one disc per kept root, with the coverage (the whole surface, or a mask image
//     in the colour texture's uv space, white = full density, only when chosen in the settings);
//   - the segment budget: one global key threshold over every part's discs of both layers BEFORE the settle
//     (a fibre's segment count is fixed by its disc), the export's own threshold rule;
//   - the curves (fur_loops.h growLoops) with the body: every non-skipped primitive (collision, MeshQuery's
//     local side, the stuck test);
//   - one population (Directed); colour per part, one colour, or each fibre the texture colour at its root
//     (material.pile_colour_source, an Export setting).
#include <optional>
#include <string>
#include <vector>

#include "fibres/fur.h"
#include "fibres/fur_field.h"
#include "scene/scene.h"
#include "scene/fibres/scene_fibres.h"

namespace knit {

struct FurMask {
    std::string path;
    int width = 0, height = 0;
    std::vector<uint8_t> value;  // per texel, row 0 = the top of the image
    uint64_t hash = 0;           // of the file bytes
};
bool loadFurMask(const std::string& path, FurMask& out, std::string* err);
// The mask at colour uv (v up, repeat wrap, nearest texel), in [0, 1].
double sampleFurMask(const FurMask& m, const Vec2& uv);

// The mask file the settings select ("" for the whole surface); err: no mask loaded, or a mask file that
// changed since it was loaded (pile.mask_hash: Reload it).
bool furMaskPath(const SceneSettings& s, std::string& path, std::string* err);

// Whether the part grows fur (not listed in pile.excluded_parts).
bool furPartGrows(const SceneSettings& s, const ScenePart& part);
// The coverage on the part's triangles (null: the whole surface, no mask loaded); false with the reason when the
// mask needs colour texture coordinates the part lacks. The function refers to `mask` and the scene: keep both.
bool furPartCoverage(const Scene& t, const ScenePart& part, const FurMask& mask, FurCoverage& out, std::string* err);
// One part's root caches, both layers (the fine fibres' only when the layer is on), drawn at the capacities
// (>= the look's densities; a larger capacity gives the same fields: fur_field.h) with the hidden test against
// every eligible solid of the scene. The stage builds them per run; the viewer's disc preview keeps them.
// One part's irregularity fields (fur_field.h buildFurModifiers), on its mesh in scene coordinates.
FurPartModifiers furPartModifiers(const Scene& t, size_t part, const FurLook& look);
struct FurPartRoots {
    int part = -1;
    FurRootCache loops, fibres;
};
FurPartRoots sampleFurPartRoots(const Scene& t, size_t part, const FurLook& look, double loopsCapacity, double fibresCapacity,
                                const ExecPolicy& exec);

// A fibre as a polyline relative to `centre` (points[0], the buried root, is the strand's first vertex; the
// visible fibre starts at points[1]).
void addHairFibre(FibreSet& fs, const Hair& h, const Vec3& centre);

std::optional<SceneFibres> growSceneFur(const SceneEntry& e, const Scene& t, const SceneSettings& s, std::string* err);

} // namespace knit
