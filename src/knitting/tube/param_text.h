#pragma once
// Parameters and curves as exact text and digests: the inputs of cache keys. Every field as text,
// doubles as exact hex floats; a static_assert on each struct's size makes a new field a build error
// until it is listed.
#include <cstdint>
#include <string>

#include "model/params.h"
#include "knitting/tube/yarn_surface.h"
#include "model/yarn_path.h"

namespace knit {

std::string patternParamsText(const PatternParams& p);
std::string surfaceParamsText(const SurfaceParams& p);
uint64_t curvesDigest(const CurveSet& cs);   // control points and stitch spans of every yarn
// The digest of a buildSurface's complete inputs (the curves, every PatternParams and SurfaceParams
// field): equal digests -> the same model.
uint64_t surfaceBuildDigest(const CurveSet& cs, const PatternParams& pp, const SurfaceParams& sp);

} // namespace knit
