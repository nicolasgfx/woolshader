#pragma once
// Overlap coverage (knit.coverage = overlap, scene_settings.h): a component is knitted as several
// independent parts (one per chart, continued beyond it, and gap patches). Its coverage is that
// of all its parts' yarn together: the parts' rendered tubes are joined into one tube (each tube
// vertex keeps its binding foot on the component's mesh, which the parts share) and measured with
// measureCoverage (knitting/knit/coverage.h).
#include <string>
#include <vector>

#include "knitting/knit/coverage.h"
#include "settings/scene_settings.h"
#include "scene/scene.h"

namespace knit {

struct ComponentGroup {
    int prim = -1, component = -1;
    std::vector<int> parts;            // indices into Scene::parts (same prim and component)
};
std::vector<ComponentGroup> componentGroups(const Scene& scene);

struct GroupTube {
    TubeMesh tube;                     // pos and tri of every knitted part's tube (component frame)
    std::vector<int> feet;             // per tube vertex: its binding foot triangle (-1: none)
    int knitted = 0;                   // parts contributing
    std::string note;                  // parts skipped (and why)
};
// The yarn tubes of the group's knitted parts at the coverage tube (0.6 mm rings, 6 sides, whatever
// the tessellation), with separate plies the yarn around them.
GroupTube groupTube(const Scene& scene, const ComponentGroup& g, const SceneSettings& settings, const ExecPolicy& exec);

// The group's coverage on its mesh, at CoverageParams' defaults.
CoverageReport groupCoverage(const Scene& scene, const ComponentGroup& g, const GroupTube& own, const ExecPolicy& exec);

struct KnitRunOptions;
// knit.coverage = overlap, knit.gap_patches (knitScene's second pass): per component, the coverage of
// its knitted chart parts together; every uncovered region at least one course wide and one stitch in
// area becomes a gap patch part (TileParams::patch) seeded at its centre, radius half its box diagonal
// plus (knit.overlap_stitches + 1) max(w, h), wale along the world up. Appended to out.parts, not
// knitted here.
void addGapPatches(Scene& out, const SceneSettings& settings, const KnitRunOptions& options);

} // namespace knit
