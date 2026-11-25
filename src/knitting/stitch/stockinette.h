#pragma once
// Stockinette (jersey) patch generator.
//
// One continuous open uniform cubic B-spline per course. A stitch template of
// kControlsPerStitch normalized controls (X in stitch widths, Y in course
// heights, Z in units of the depth amplitude) is replicated with
//   x = (X + i - cx) * waleSpacing,  y = (Y + r - cy) * courseSpacing,  z = Z * A * D.
// The template's absolute depth sequence is a CHOICE that realizes the expected
// signature (expected_signature.h); it is not implied by topology.
//
// Interior patch only: course 0 has no course below (cast-on omitted), the top
// course's heads are uncaught (bind-off omitted), and row ends are open (edge
// turns omitted). Odd courses are traversed in -U, as in flat knitting.
#include <vector>

#include "model/params.h"
#include "model/topology.h"
#include "model/yarn_path.h"

namespace knit {

inline constexpr int kControlsPerStitch = 16;

struct TemplateControl { double X, Y, Z; };

// Normalized controls j = 0..15; j = 0 is the sinker at X = -0.5, j = 8 the crown.
std::vector<TemplateControl> stitchTemplate(const PatternParams& pp);
// The head-side X of the template (stitch_layout.cpp's capHeadSide): the widest head control, X_10 = X_11 =
// headHalfWidth + 0.02.
double templateHeadSideX(const PatternParams& pp);

// Patterns: the depth signs of a stitch from the faces involved, sigma = +1
// knit, -1 purl:
//   j 0-3, 13-15 (sinker, foot, leg bottom, leg; Y <= 0.30): the stitch's interlock with
//                its PARENT (its legs over the parent's head)       -> sigma_own
//   j 4-12       (neck, head sides, crown; Y >= 0.60): passed through by its CHILD
//                                                                   -> sigma_child
//   j 0          the sinker shared with the left neighbour: depth 0 where the faces
//                of the two differ (the yarn passes between front and back there).
//   rib 1x1 / 2x2 with pp.ribLayerOffset != 0: every control except a transition sinker moves by
//                sigma_own x ribLayerOffset x D along N (knit wales in front, purl wales behind).
// All knit = stitchTemplate(pp) exactly (the signs multiply by +1.0).
struct StitchFaces {
    StitchFace own = StitchFace::Knit, child = StitchFace::Knit, left = StitchFace::Knit;
    bool operator==(const StitchFaces&) const = default;
};
std::vector<TemplateControl> stitchTemplate(const PatternParams& pp, const StitchFaces& f);
// Faces of the flat patch's stitch (wale i, course r) for pp.pattern. The patch is a
// WINDOW of the infinite pattern: the virtual neighbours outside it (the padding
// controls' stitches i = -1 and i = wales, and the top course's child row) take the
// pattern's faces, so every stitch's template depends on its lattice neighbourhood only.
StitchFaces flatStitchFaces(const PatternParams& pp, int wale, int course);

// Template control index -> local position for stitch (wale i, course r).
Vec3 templateToLocal(const TemplateControl& tc, int wale, int course, const PatternParams& pp);

CurveSet generateStockinette(const PatternParams& pp);

} // namespace knit
