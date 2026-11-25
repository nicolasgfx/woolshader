#pragma once
// renderedPart: the strands a knitted part is drawn and exported with: the yarn, or its plies.
//
// From a ScenePart and the settings. Called by the fibre stage, the export, the display and the
// fibre preview. Empty with the reason when the strands cannot be built.
#include <memory>
#include <optional>
#include <string>

#include "knitting/tube/plies.h"
#include "knitting/tube/yarn_surface.h"
#include "scene/scene.h"

namespace knit {

// A knitted part's rendered strands: the yarn at the final surface's D_r, or the
// explicit plies (inside the drawn yarn), with the export ring / sides and the
// curvature bound for the fibre root depth.
struct RenderedPart {
    FinalSurfaceSpec fs;
    PatternParams ppR;                 // the yarn at D_r
    const CurveSet* yarn = nullptr;    // the knitted yarn (fibres are authored on it)
    const CurveSet* strands = nullptr; // the exported strands: the yarn, or the plies
    PatternParams strandPp;            // their diameter
    std::unique_ptr<PlySet> plies;
    std::unique_ptr<SurfaceModel> parentSm;
    double ring = 0;
    int sides = 0;
    double kappaMax = 0;
    int plyFolds = 0, plyFoldsUnresolved = 0;
};
// Moving one keeps its pointers: `strands` points into the heap PlySet or at the knitted yarn.
//
// A knitted part's strands. The reason starts with the part's name.
std::optional<RenderedPart> renderedPart(const ScenePart& part, const SceneSettings& settings, std::string* err);
// The same for any knitted yarn (the reference preview's flat patch): `yarn` must outlive the result;
// kappaMaxYarn: the single tube's curvature bound (the knit's gate (a) result, or measured). The reason
// names no part: the caller prefixes it.
std::optional<RenderedPart> renderedStrands(const CurveSet& yarn, const PatternParams& pp, double kappaMaxYarn, const SceneSettings& settings,
                                            std::string* err);

} // namespace knit
