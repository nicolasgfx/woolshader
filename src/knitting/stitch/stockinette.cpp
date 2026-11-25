#include "knitting/stitch/stockinette.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace knit {

double templateHeadSideX(const PatternParams& pp) { return pp.headHalfWidth + 0.02; }  // right[2].X = right[3].X below

std::vector<TemplateControl> stitchTemplate(const PatternParams& pp) {
    const double Xh = pp.headHalfWidth, Xl = pp.legHalfWidth;
    // Right half, crown (j=8) to sinker (j=16 == next stitch's j=0).
    // Depth: crown region B, head side above the depth transition B, head side
    // below it F, neck/leg F, leg bottom and foot/sinker B. The yarn changes depth
    // where it passes through the neighbouring course's loop.
    const TemplateControl right[9] = {
        {0.0, 1.40, -1.0},            // j=8  crown
        {0.64 * Xh, 1.38, -1.0},      // j=9  head top
        {Xh + 0.02, 1.24, -1.0},      // j=10 head side, above the transition
        {Xh + 0.02, 0.95, +1.0},      // j=11 head side, below the transition
        {0.80 * Xh, 0.60, +1.0},      // j=12 neck (upper)
        {Xl, 0.30, +1.0},             // j=13 leg
        {Xl, -0.05, -1.0},            // j=14 leg bottom
        {0.5 * (Xl + 0.5), -0.30, -1.0},  // j=15 foot
        {0.5, -0.36, -1.0},           // j=16 sinker
    };
    std::vector<TemplateControl> t(kControlsPerStitch);
    for (int j = 8; j < 16; ++j) t[j] = right[j - 8];
    for (int j = 0; j < 8; ++j) {  // mirror: j <-> 16 - j
        const TemplateControl& m = right[16 - j - 8];
        t[j] = {-m.X, m.Y, m.Z};
    }
    return t;
}

std::vector<TemplateControl> stitchTemplate(const PatternParams& pp, const StitchFaces& f) {
    std::vector<TemplateControl> t = stitchTemplate(pp);
    const double own = f.own == StitchFace::Purl ? -1.0 : 1.0, child = f.child == StitchFace::Purl ? -1.0 : 1.0;
    for (int j = 0; j < kControlsPerStitch; ++j) t[j].Z *= (j >= 4 && j <= 12) ? child : own;
    if (f.left != f.own) t[0].Z = 0.0;
    if (f.child != f.own) {
        // Both halves (j and 16 - j mirror each other).
        for (int j : {4, 12}) {
            t[j].Y = pp.transNeckY;
            t[j].Z *= pp.transNeckZ;
        }
        for (int j : {3, 13}) t[j].Y = pp.transLegY;
    }
    if (pp.pattern == StitchPattern::Garter && pp.garterWave != 0) {
        // cos(pi (course + Y)) = sigma_own cos(pi Y) in garter (knit rows even): the same wave for
        // every stitch, continuous across stitches and rows. Z is in units of A D.
        const double a = pp.garterWave / pp.depthAmplitude;
        for (TemplateControl& q : t) q.Z -= a * own * std::cos(3.14159265358979323846 * q.Y);
    }
    if (pp.pattern == StitchPattern::Rib2x2 && pp.ribPairSpread != 0) {
        // A pair's first wale (its left neighbour has the other face) moves -spread, its second
        // +spread. The sinker j = 0 lies halfway between this wale and its left neighbour: inside a
        // pair (first -spread, second +spread) at 0, between pairs (second +spread, first -spread) at 0.
        const double shift = (f.left != f.own ? -1.0 : 1.0) * pp.ribPairSpread;
        for (int j = 1; j < kControlsPerStitch; ++j) t[j].X += shift;
    }
    if (patternHasRibLayers(pp.pattern) && pp.ribLayerOffset != 0) {
        // Z is in units of A x D (templateToLocal): the layer offset in D becomes ribLayerOffset / A.
        const double layer = own * pp.ribLayerOffset / pp.depthAmplitude;
        for (int j = 0; j < kControlsPerStitch; ++j)
            if (j != 0 || f.left == f.own) t[j].Z += layer;
    }
    return t;
}

StitchFaces flatStitchFaces(const PatternParams& pp, int wale, int course) {
    auto face = [&](int i, int r) { return stitchFaceIsPurl(pp.pattern, i, r) ? StitchFace::Purl : StitchFace::Knit; };
    return {face(wale, course), face(wale, course + 1), face(wale - 1, course)};
}

Vec3 templateToLocal(const TemplateControl& tc, int wale, int course, const PatternParams& pp) {
    const double cx = 0.5 * (pp.wales - 1), cy = 0.5 * (pp.courses - 1) + 0.52;
    return {(tc.X + wale - cx) * pp.waleSpacing, (tc.Y + course - cy) * pp.courseSpacing,
            tc.Z * pp.depthAmplitude * pp.yarnDiameter};
}

CurveSet generateStockinette(const PatternParams& pp) {
    CurveSet cs;
    cs.kind = CurveKind::FlatPatch;
    cs.wales = pp.wales;
    cs.courses = pp.courses;
    cs.waleSpacing = pp.waleSpacing;
    cs.courseSpacing = pp.courseSpacing;
    const int K = kControlsPerStitch, n = pp.wales;
    // One template per face combination (all knit: stitchTemplate(pp) itself).
    std::vector<std::pair<StitchFaces, std::vector<TemplateControl>>> tpls;
    auto tplFor = [&](int i, int r) -> const std::vector<TemplateControl>& {
        const StitchFaces f = flatStitchFaces(pp, i, r);
        for (const auto& [k, t] : tpls)
            if (k == f) return t;
        tpls.push_back({f, stitchTemplate(pp, f)});
        return tpls.back().second;
    };
    for (int r = 0; r < pp.courses; ++r) {
        YarnPath y;
        y.id = r;
        y.course = r;
        y.name = "course " + std::to_string(r);
        y.reversed = (r % 2) == 1;
        // Global control index g in [-1, nK+1]: stitch i = floor(g/K), j = g mod K.
        // The knot at u = g (before reversal) is centred on control g, so stitch i
        // occupies u in [iK, (i+1)K] with its sinkers at the stitch boundaries.
        for (int g = -1; g <= n * K + 1; ++g) {
            int i = (g >= 0) ? g / K : -1;
            int j = g - i * K;
            y.curve.ctrl.push_back(templateToLocal(tplFor(i, r)[j], i, r, pp));
        }
        if (y.reversed) std::reverse(y.curve.ctrl.begin(), y.curve.ctrl.end());
        for (int i = 0; i < n; ++i) {
            StitchSpan st;
            st.wale = i;
            st.course = r;
            st.id = r * n + i;
            if (!y.reversed) { st.uBegin = i * K; st.uEnd = (i + 1) * K; }
            else { st.uBegin = (n - 1 - i) * K; st.uEnd = (n - i) * K; }
            y.stitches.push_back(st);
            // Topology: stitch (i, r) is drawn through (i, r - 1). The frame is the
            // fabric frame at the template origin (U = +x is the template's X axis).
            StitchRecord rec;
            rec.id = st.id;
            rec.yarn = r;
            rec.uBegin = st.uBegin;
            rec.uEnd = st.uEnd;
            rec.turn = r;
            rec.index = i;
            rec.face = flatStitchFaces(pp, i, r).own;
            if (r > 0) { rec.parents[0] = (r - 1) * n + i; rec.nParents = 1; }
            rec.origin = templateToLocal({0, 0, 0}, i, r, pp);
            cs.topo.stitches.push_back(rec);
        }
        std::sort(y.stitches.begin(), y.stitches.end(),
                  [](const StitchSpan& a, const StitchSpan& b) { return a.uBegin < b.uBegin; });
        y.finalize();
        cs.yarns.push_back(std::move(y));
    }
    std::vector<double> ends;
    for (const YarnPath& y : cs.yarns) ends.push_back(y.uEnd());
    cs.topo.buildIndex(ends);
    return cs;
}

} // namespace knit
