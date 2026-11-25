#pragma once
// The yarn-settle continuation: the embedded yarn
// (piecewise cubic Hermite pieces, yarn_embed.h) in production's curve representation (uniform cubic B-splines,
// model/yarn_path.h), so the existing knitting solver (knitting/numerics.h, KnitNumerics::relax) can move it, and back.
//
// Conversion. Not exact: the Hermite strand is C1 inside a piece and G1 at its ports, a uniform cubic B-spline is C2. Each
// open strand is evaluated densely (densePerSpan points per Hermite span), resampled at equal arc length (n intervals of
// `spacing` = L / n, n = round(L / requested spacing)), and replaced by the uniform cubic B-spline that interpolates those
// points at u = 0, 1, ..., n with natural ends (interpolateBSpline, knitting/tube/plies.h). So u = arc length / spacing on
// the input. Every piece keeps its interval [u0, u1] (the arc positions of its two ports) and every Hermite knot its
// parameter: the pieces, ports, strands, faces and templates of the embedding survive a motion of the controls, and
// embeddingFromCurves re-evaluates them on the moved curve. The approximation is measured (YarnConversion), never assumed.
//
// A closed strand is not converted (production's curves are open): yarnToCurves refuses it with the reason.
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/yarn_path.h"
#include "stitchmesh/yarn_embed.h"

namespace knit {

struct YarnConversion {
    CurveSet curves;                          // one yarn per strand, CurveKind::Surface, no stitch records (topo empty)
    std::vector<double> spacing;              // per strand: the arc spacing of its interpolation points [mm]
    std::vector<int> strandOfPiece;           // per piece
    std::vector<double> u0, u1;               // per piece: its interval on its strand's curve
    std::vector<std::vector<double>> knotU;   // per piece: the parameter of each of its Hermite knots
    // measured on the input
    double hermiteLength = 0, curveLength = 0;  // dense Hermite polyline; the B-spline (arc-length table) [mm]
    double maxDeviation = 0, p99Deviation = 0;  // dense Hermite points to the nearest point of the curve (16 samples per span) [mm]
    double maxEqualParameter = 0;               // dense Hermite points to the curve at their own parameter [mm]
    int densePoints = 0;
    std::string text() const;
};

// Every strand of e must be open. st: the embedding's settings (its tangent rule). False with the reason otherwise.
bool yarnToCurves(const YarnEmbedding& e, const YarnSettings& st, double spacing, int densePerSpan, YarnConversion& out, std::string* err);

// The embedding on curves of the same structure as cv.curves (for example after a relaxation moved the controls): every
// piece gets as many samples as before, uniform in u over its interval; its knots the curve points at their parameters;
// tangents the curve's derivative there; depth the signed offset along the reference normal [D]. Strands are reassembled
// as embedYarn does (a shared port sample once). Ports, faces, templates and the log stay.
YarnEmbedding embeddingFromCurves(const YarnEmbedding& e, const YarnConversion& cv, const CurveSet& curves, const TriMesh& reference, double D);

} // namespace knit
