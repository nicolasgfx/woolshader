// The strands of a knitted part as they are drawn: the yarn at its final diameter, or its plies.
#include "scene/appearance/rendered_part.h"

#include <algorithm>

#include "gpu/knitting.h"
#include "gpu/surfaces.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"
#include "knitting/knit/mesh_knit.h"
#include "scene/scene.h"

namespace knit {

std::optional<RenderedPart> renderedPart(const ScenePart& part, const SceneSettings& settings, std::string* err) {
    if (!part.knit) {
        if (err) *err = part.name + ": not knitted";
        return std::nullopt;
    }
    const PartKnit& k = *part.knit;
    double kappaMax = 0;
    for (const CurvatureResult& c : k.after.curvature) kappaMax = std::max(kappaMax, std::min(c.kappaHi, 1e3));
    std::optional<RenderedPart> out = renderedStrands(k.t.curves, k.t.pp, kappaMax, settings, err);
    if (!out && err) *err = part.name + ": " + *err;
    return out;
}

std::optional<RenderedPart> renderedStrands(const CurveSet& yarn, const PatternParams& pp, double kappaMaxYarn, const SceneSettings& settings,
                                            std::string* err) {
    const SceneSettings& s = settings;
    RenderedPart out;
    out.yarn = &yarn;
    out.fs = finalSurface(s, pp.yarnDiameter);
    if (!out.fs.valid(err)) return std::nullopt;
    out.ppR = out.fs.rendered(pp);
    if (!tubeSettingsValid(s, err)) return std::nullopt;
    if (s.plies > 1) {
        // Plies inside the drawn yarn (plies.h), authored on its round model: the validated
        // envelope, or for a yarn compressed to fit its full diameter, past that envelope.
        PlyGeometry g;
        g.count = s.plies;
        g.pitch = s.plyPitch;
        g.fill = s.plyFill;
        g.handedness = s.plyHandedness;
        SurfaceParams spY;
        spY.compression = false;
        out.parentSm = std::make_unique<SurfaceModel>();
        out.plies = std::make_unique<PlySet>();
        std::string gpuError;
        if (!gpu::buildSurface({yarn, out.ppR, spY}, *out.parentSm, nullptr, &gpuError) ||
            !gpu::buildPlies({yarn, *out.parentSm, g, out.ppR.radius()}, *out.plies, nullptr, &gpuError)) {
            if (err) *err = "plies on the GPU: " + gpuError;
            return std::nullopt;
        }
        const std::unique_ptr<KnitNumerics> numerics = gpu::makeKnitNumerics(&gpuError);
        if (!numerics) {
            if (err) *err = "plies on the GPU: " + gpuError;
            return std::nullopt;
        }
        out.strands = &out.plies->curves;
        out.strandPp = pp;
        out.strandPp.yarnDiameter = 2 * out.plies->layout.rp;
        out.ring = s.plyTubeRing;  // the export's ply tube (strand.ply_tube_*)
        out.sides = s.plyTubeSides;
        const Allowances al = Allowances::forExtent(out.plies->curves.extent());
        for (const YarnPath& yp : out.plies->curves.yarns) {
            const CurvatureResult cr = numerics->curvature(yp, out.plies->layout.rp, al, ExecPolicy{});
            out.plyFolds += cr.noFold == Status::Fail;
            out.plyFoldsUnresolved += cr.noFold == Status::Unresolved;
            for (double kh : cr.spanKappaHi)
                if (kh < 1e6) out.kappaMax = std::max(out.kappaMax, kh);
        }
    } else {
        out.strands = &yarn;
        out.strandPp = out.ppR;
        out.ring = s.tubeRing;  // the export's yarn tube (strand.tube_*)
        out.sides = s.tubeSides;
        out.kappaMax = kappaMaxYarn;
    }
    return out;
}

} // namespace knit
