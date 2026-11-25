// The questions asked of the reviewed results; see results.h.
#include "jobs/results.h"

#include <algorithm>
#include <cstdio>

#include "model/yarn_path.h"
#include "scene/result_keys.h"
#include "scene/scene.h"
#include "settings/scene_registry.h"
#include "settings/settings_profiles.h"

namespace knit {

uint64_t fibreResultKeyOf(const SceneSettings& s) { return fibreResultKey(s); }

double longestYarnMm(const KnitResult& k) {
    double longest = 0;
    if (!k.scene) return 0;
    for (const ScenePart& part : k.scene->parts) {
        if (!part.knit) continue;
        for (const YarnPath& y : part.knit->t.curves.yarns) longest = std::max(longest, y.arcLengthAt(y.curve.uMax()));
    }
    return longest;
}

bool knitOutdated(const KnitResult& k, const SceneSettings& current) { return geometryKey(k.settings) != geometryKey(current); }

bool fibresOutdated(const FibreResult& f, const SceneSettings& current) { return f.key != fibreResultKey(current); }

ExportEligibility exportEligibility(const SceneSettings& current, const KnitResult* knit, const FibreResult* fibres, FibreChoice choice) {
    const bool includeFibres = choice == FibreChoice::Included;
    ExportEligibility e;
    bool blocked = false;
    auto block = [&](const std::string& why, const std::string& fix) {
        blocked = true;
        e.reasons.push_back(why);
        if (e.fix.empty()) e.fix = fix;
    };
    if (furWorkflow(current)) {
        // Fur: no knit is involved; the export writes the fur result's base meshes and its pile.
        if (!fibres) block("There is no fur yet.", "Grow fur");
        else if (fibres->workflow != Workflow::Fur || !fibres->base) block("The fibres shown are wool fibres, not fur.", "Grow fur");
        else if (includeFibres && fibresOutdated(*fibres, current)) block("Fur settings changed since the fur was grown.", "Update fur");
        if (!blocked) e.state = ExportEligibility::State::Ready;
        return e;
    }
    if (!knit) {
        block("There is no knitted result yet.", "Knit at Final quality");
        return e;
    }
    if (!knit->computational) block("The knitted result was restored for viewing only (" + knit->viewOnlyReason + ").", "Knit again");
    if (!knit->final)
        block(isPreviewProfile(knit->settings) ? "The knitted result is a Draft: it skipped settle and the surface check."
                                               : "The knitted result used a custom profile that skips some Final validation steps.",
              "Knit at Final quality");
    if (knitOutdated(*knit, current)) block("Knitting settings changed since the knitted result was made.", "Update the knit");
    if (current.yarnBody.enabled) {
        // The export would fail on these fibres: said here, with what to change, before it is started.
        const YarnBodyCheck body = checkYarnBody(current);
        if (!body.accepted) {
            block("The yarn's fibres are refused: " + body.reason + ".", "Correct the yarn's fibres");
            e.reasons.push_back(body.remedy + " (Fibres step, Yarn fibres)");
        }
        const YarnBodyCheck maps = checkYarnMaps(current, longestYarnMm(*knit));
        if (!maps.accepted) {
            block("The yarn's maps are refused: " + maps.reason + ".", "Correct the yarn's fibres");
            e.reasons.push_back(maps.remedy + " (Fibres step, Yarn fibres)");
        }
    }
    if (includeFibres) {
        if (!fibres) block("Fibres are included but none have been generated for this knit.", "Generate fibres");
        else if (fibres->knitId != knit->id) block("The fibres belong to an earlier knitted result.", "Generate fibres");
        else if (fibresOutdated(*fibres, current)) block("Fibre settings changed since the fibres were generated.", "Update fibres");
    }
    if (blocked) return e;  // e.state is Blocked
    if (!knit->failing.empty()) {
        char b[160];
        std::snprintf(b, sizeof b, "%zu of %d knitted parts did not pass validation.", knit->failing.size(), knit->parts - knit->plain);
        e.reasons.push_back(b);
    }
    if (knit->cross != Status::Pass)
        e.reasons.push_back(std::string("The clearance between parts is ") + toString(knit->cross) + ".");
    e.state = e.reasons.empty() ? ExportEligibility::State::Ready : ExportEligibility::State::Diagnostic;
    return e;
}

} // namespace knit
