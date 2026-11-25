// A part's yarn as a tube mesh, on the device; see yarn_tube.h.
#include "scene/appearance/yarn_tube.h"

#include "gpu/surfaces.h"

namespace knit {

SurfaceParams surfaceParamsOf(const TubeStyle& style) {
    SurfaceParams sp;
    sp.compression = style.compression;
    sp.maxCompression = style.maxCompression;
    sp.ringSpacing = style.modelRing;
    return sp;
}

bool buildYarnTube(const CurveSet& cs, const PatternParams& pp, const TubeStyle& style, TubeMesh& out, std::string* err) {
    const SurfaceParams sp = surfaceParamsOf(style);
    TubeParams tp;
    tp.ringSpacing = style.ring;
    tp.sides = style.sides;
    tp.uvRepeat = style.uvRepeat;
    tp.uvWrap = style.uvWrap;
    tp.uvPhaseJitter = style.uvPhaseJitter;
    tp.uvPitchJitter = style.uvPitchJitter;
    return gpu::buildTubeMesh({cs, pp, sp}, tp, out, nullptr, err);
}

bool buildYarnSurface(const CurveSet& cs, const PatternParams& pp, const TubeStyle& style, SurfaceModel& out, std::string* err) {
    return buildYarnSurface(cs, pp, surfaceParamsOf(style), out, err);
}

bool buildYarnSurface(const CurveSet& cs, const PatternParams& pp, const SurfaceParams& sp, SurfaceModel& out, std::string* err) {
    return gpu::buildSurface({cs, pp, sp}, out, nullptr, err);
}

} // namespace knit
