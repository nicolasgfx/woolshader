// The Fibres stage's result and what is read from it; see scene_fibres.h.
#include "scene/fibres/scene_fibres.h"

#include <map>

#include "fibres/fibre_budget.h"
#include "model/log.h"
#include "settings/settings_mapping.h"

namespace knit {

PartStrands partStrands(const PartFibres& fibres, const Vec3& centreMm, FibrePopulation population, const StrandTaper& taper) {
    PartStrands out;
    if (taper.tail > 0) {
        TaperedStrands t = fibresToTaperedStrands(fibres.fs, sceneMetres(centreMm), int(population), taper);
        out.strands = std::move(t.strands);
        out.tapered = t.tapered;
        out.untapered = t.untapered;
    } else {
        out.strands = fibresToStrands(fibres.fs, sceneMetres(centreMm), int(population));
    }
    for (const Fibre& f : fibres.fs.fibres)
        if (f.population == population) out.keys.push_back(float(f.key));
    return out;
}

KeptSegments keptSegments(const PartFibres& fibres, const Vec3& centreMm, const StrandTaper& taper) {
    KeptSegments out;
    for (const Fibre& f : fibres.fs.fibres)
        if (f.pointCount >= 2 && float(f.key) < fibres.keyBelow) out.kept += f.segments();
    out.exported = taperedSegmentCount(fibres.fs, sceneMetres(centreMm), taper, fibres.keyBelow);
    return out;
}

TaperReport partTaperReport(const PartFibres& fibres, const Vec3& centreMm, const StrandTaper& taper) {
    return measureTaper(fibres.fs, sceneMetres(centreMm), taper, fibres.keyBelow);
}

void addKeySegments(const FibreSet& fibres, std::vector<std::pair<float, int>>& out) {
    for (const Fibre& f : fibres.fibres) out.push_back({float(f.key), f.segments()});
}

long long segmentCount(const FibreSet& fibres) {
    long long segments = 0;
    for (const Fibre& f : fibres.fibres) segments += f.segments();
    return segments;
}

std::string applyObjectBudgets(const Scene& scene, long long budget, std::vector<PartFibres>& parts) {
    std::map<int, std::vector<PartFibres*>> partsOfObject;
    for (PartFibres& pf : parts) partsOfObject[sceneObjectOf(scene.source, scene.parts[static_cast<size_t>(pf.part)].prim)].push_back(&pf);
    std::string log;
    for (const auto& [object, members] : partsOfObject) {
        std::vector<std::pair<float, int>> keySegments;
        for (const PartFibres* pf : members) addKeySegments(pf->fs, keySegments);
        const float q = budgetThreshold(keySegments, budget);
        long long kept = 0;
        for (PartFibres* pf : members) {
            pf->keyBelow = q;
            for (const Fibre& f : pf->fs.fibres) kept += float(f.key) < q ? f.segments() : 0;
        }
        appendFormat(log, "fibres: %s: the segment budget %lld keeps key < %.6f: %lld segments in %zu parts\n",
                     sceneObjectName(scene.source, scene.parts[static_cast<size_t>(members[0]->part)].prim).c_str(), budget, double(q), kept, members.size());
    }
    return log;
}

} // namespace knit
