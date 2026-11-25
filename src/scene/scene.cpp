// Questions about a scene: what is knitted and what is missing.
#include "scene/scene.h"

#include <set>

namespace knit {

bool partIncomplete(const ScenePart& part, const std::string& only) {
    if (part.knit || (!only.empty() && part.name.find(only) == std::string::npos)) return false;
    return part.stop != KnitStop::None || part.layout == KnitLayout::Spiral || part.layout == KnitLayout::Guided;
}

const char* toString(GuidedOutcome o) {
    switch (o) {
    case GuidedOutcome::Settled: return "settled";
    case GuidedOutcome::SettledMeshFails: return "settled (the stitch mesh checks do not pass)";
    case GuidedOutcome::NotSettled: return "not settled";
    case GuidedOutcome::Drafted: return "draft (not settled)";
    case GuidedOutcome::SettleSkipped: return "not settled (settle skipped)";
    case GuidedOutcome::Stopped: return "stopped";
    case GuidedOutcome::Refused: return "refused";
    }
    return "?";
}

GuidedOutcome outcomeAfterSettle(Status targetGate, Status construction, Status meshChecks) {
    if (targetGate != Status::Pass || construction != Status::Pass) return GuidedOutcome::NotSettled;
    return meshChecks == Status::Pass ? GuidedOutcome::Settled : GuidedOutcome::SettledMeshFails;
}

bool hasChartParts(const Scene& scene) {
    for (const ScenePart& part : scene.parts)
        if (part.layout == KnitLayout::UvTiling) return true;
    return false;
}

bool ScenePart::hasStitches() const { return knit && !knit->t.prog.stitches.empty(); }

int sceneObjectOf(const SourceScene& source, int prim) {
    const int node = source.prims[static_cast<size_t>(prim)].nodeIndex;
    return node >= 0 ? node : -1 - prim;
}

std::string sceneObjectName(const SourceScene& source, int prim) {
    const SourcePrimitive& g = source.prims[static_cast<size_t>(prim)];
    return g.node.empty() ? "object " + std::to_string(sceneObjectOf(source, prim)) : g.node;
}

std::vector<std::string> plainPrimitiveNames(const SourceScene& source) {
    std::vector<std::string> names(source.prims.size());
    std::set<std::string> seen;
    for (size_t i = 0; i < source.prims.size(); ++i) {
        const SourcePrimitive& g = source.prims[i];
        if (g.cls != SourcePrimitive::Class::Plain) continue;
        std::string name = g.mesh + "_" + g.material;
        for (int n = 2; seen.count(name); ++n) name = g.mesh + "_" + g.material + "_" + std::to_string(n);  // instances
        seen.insert(name);
        names[i] = name;
    }
    return names;
}

} // namespace knit
