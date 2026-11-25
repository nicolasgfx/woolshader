// A component knitted as several parts: their tubes joined and measured, and patches for what stays open.
#include "scene/knit/gap_patches.h"

#include "scene/appearance/yarn_tube.h"
#include "scene/stage_params.h"

#include <algorithm>
#include <map>

#include "model/cancel.h"
#include "scene/knit/knit_stage.h"
#include "scene/scene.h"
#include "settings/settings_mapping.h"

namespace knit {

namespace {

// The tube gap coverage measures. Fixed, not strand.tube_*: those are Fibres rows (the yarn's
// tessellation), and which gap patches exist is part of the knit, whose key holds Load and Knit rows
// only. The values are those of Reduced, the default tessellation, so default settings knit as before.
constexpr double kCoverageTubeRing = 0.6;  // [mm]
constexpr int kCoverageTubeSides = 6;

} // namespace

std::vector<ComponentGroup> componentGroups(const Scene& scene) {
    std::vector<ComponentGroup> out;
    std::map<std::pair<int, int>, size_t> at;
    for (size_t i = 0; i < scene.parts.size(); ++i) {
        const ScenePart& part = scene.parts[i];
        if (part.layout != KnitLayout::UvTiling) continue;  // a spiral covers its whole piece: no gaps to patch
        const std::pair<int, int> key{part.prim, part.component};
        auto it = at.find(key);
        if (it == at.end()) {
            it = at.emplace(key, out.size()).first;
            out.push_back({part.prim, part.component, {}});
        }
        out[it->second].parts.push_back(static_cast<int>(i));
    }
    return out;
}

GroupTube groupTube(const Scene& scene, const ComponentGroup& g, const SceneSettings& settings, const ExecPolicy& exec) {
    GroupTube out;
    for (int pi : g.parts) {
        cancelCheckpoint();
        const ScenePart& part = scene.parts[pi];
        if (!part.knit) continue;
        // The yarn's drawn tube, with or without separate plies: the plies lie inside it, so it is
        // what the fabric covers.
        const PatternParams& pp = part.knit->t.pp;
        const FinalSurfaceSpec fs = finalSurface(settings, pp.yarnDiameter);
        std::string entry;
        if (!fs.valid(&entry)) {
            out.note += part.name + ": " + entry + "; ";
            continue;
        }
        TubeStyle style;
        style.modelRing = SurfaceParams{}.ringSpacing;  // the model's own spacing
        style.ring = kCoverageTubeRing;
        style.sides = kCoverageTubeSides;
        style.compression = fs.compression;
        style.maxCompression = fs.maxCompression;
        TubeMesh tube;
        if (!buildYarnTube(part.knit->t.curves, fs.rendered(pp), style, tube, &entry)) {
            out.note += part.name + ": " + entry + "; ";
            continue;
        }
        const std::vector<int> feet = tubeVertexFeet(tube, part.knit->t.curves, part.knit->t.bind.state);
        // Parts of one component share its mesh and centre: positions and feet need no change.
        const uint32_t base = static_cast<uint32_t>(out.tube.pos.size());
        out.tube.pos.insert(out.tube.pos.end(), tube.pos.begin(), tube.pos.end());
        for (uint32_t v : tube.tri) out.tube.tri.push_back(base + v);
        out.feet.insert(out.feet.end(), feet.begin(), feet.end());
        ++out.knitted;
    }
    return out;
}

CoverageReport groupCoverage(const Scene& scene, const ComponentGroup& g, const GroupTube& own, const ExecPolicy& exec) {
    CoverageInput in;
    in.mesh = scene.parts[g.parts.front()].mesh.get();
    in.own.tube = &own.tube;
    in.own.vertexFoot = own.feet;
    return measureCoverage(in, CoverageParams{}, exec);
}

void addGapPatches(Scene& out, const SceneSettings& settings, const KnitRunOptions& options) {
    // The chart knit's own stitch (the pattern's relaxed size; checkKnitSettings accepted the pattern).
    const StitchSize st = patternStitch(settings, parsePattern(settings.pattern).value_or(StitchPattern::Stockinette));
    const double w = st.w, h = st.h;
    const int maxPatches = 400;  // per component
    for (const ComponentGroup& g : componentGroups(out)) {
        cancelCheckpoint();
        // A copy, not a reference into out.parts: the push_back below can move every part, and a
        // reference read after it was the null mesh of the explicit-ply crash (UvAtlas on *part.mesh).
        const ScenePart& head = out.parts[g.parts.front()];
        if (head.chart < 0) continue;  // not an overlap component
        struct ComponentOf {
            int prim, component;
            std::string componentName;
            Vec3 centre;
            std::shared_ptr<const TriMesh> mesh;
        };
        const ComponentOf first{head.prim, head.component, head.componentName, head.centre, head.mesh};
        const GroupTube gt = groupTube(out, g, settings, options.inner);
        if (gt.knitted == 0) continue;  // nothing measured: that is no hole
        const CoverageReport c = groupCoverage(out, g, gt, options.inner);
        int n = 0;
        for (const CoverageRegion& reg : c.holes.items) {
            // Openings of the fabric itself (narrower than a course) are not patched.
            if (n >= maxPatches || reg.centreTri < 0 || reg.width < h || reg.area < w * h) continue;
            ScenePart pr;
            pr.prim = first.prim;
            pr.component = first.component;
            pr.componentName = first.componentName;
            pr.name = first.componentName + "_gap" + std::to_string(n);
            pr.centre = first.centre;
            pr.mesh = first.mesh;
            pr.gapPatch = true;
            pr.patch.tri = reg.centreTri;
            for (int k = 0; k < 3; ++k) pr.patch.bary[k] = reg.centreBary[k];
            pr.patch.wale = TileParams{}.up;
            pr.patch.radius = 0.5 * norm(reg.hi - reg.lo) + (std::max(0, settings.overlapStitches) + 1) * std::max(w, h);
            out.parts.push_back(std::move(pr));
            ++n;
        }
        if (options.verbose)
            std::printf("  %s: coverage after the chart parts %.2f %% uncovered (%d regions, max width %.1f mm): %d gap patches\n",
                        first.componentName.c_str(), 100 * c.uncoveredFrac, c.holes.count, c.holes.maxWidth, n);
    }
}

} // namespace knit
