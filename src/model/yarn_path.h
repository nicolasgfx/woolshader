#pragma once
// A yarn path: one continuous physical strand (here: one open course of the patch).
//
// Two parameterizations, kept distinct:
//   u  native B-spline parameter in [0, numSpans]  (knots at integers; material coordinate)
//   s  arc length in mm, s in [0, L]
// Arc length is an ESTIMATE (adaptive Gauss-Kronrod); where window membership
// matters, use lengthBounds(), which encloses the true length.
#include <string>
#include <vector>

#include "model/bspline.h"
#include "model/parallel.h"
#include "model/topology.h"

namespace knit {

// One stitch (one repeat of the structure) as a parameter range on a yarn path.
struct StitchSpan {
    int wale = 0;
    int course = 0;
    double uBegin = 0, uEnd = 0;  // uBegin < uEnd; ownership is half-open [uBegin, uEnd)
    int id = -1;                  // stitch id in CurveSet::topo (-1: none recorded)
};

// An enclosure of a true arc length [mm]: chords below, Bezier control polygons above.
struct ArcLengthBounds {
    double lo = 0;
    double hi = 0;
};

class YarnPath {
public:
    int id = 0;            // stable identifier (= course index for the stockinette patch)
    std::string name;
    int course = 0;
    bool reversed = false;  // traversal runs in -U (flat knitting alternates rows)
    BSpline curve;
    std::vector<StitchSpan> stitches;
    // Radius factor per control, for the drawn yarn only; empty: 1 everywhere. The centreline's gates
    // never read it.
    std::vector<float> thickness;

    double thicknessAt(double u) const;
    double largestThickness() const;

    // Build the arc-length table; tolS is the target absolute length tolerance [mm].
    // exec: the spans' integrations and bounds in parallel; the sums are accumulated serially
    // (identical tables).
    void finalize(double tolS = 1e-9, const ExecPolicy& exec = {});

    double uBegin() const { return 0.0; }
    double uEnd() const { return curve.uMax(); }
    double length() const { return tables_.length; }
    double lengthErrorEstimate() const { return tables_.lengthErr; }

    double arcLengthAt(double u) const;
    double paramAtArcLength(double s) const;

    Vec3 positionAtParam(double u) const { return curve.eval(u); }
    Vec3 start() const { return curve.eval(uBegin()); }
    Vec3 end() const { return curve.eval(uEnd()); }

    // Enclosure of the true arc length between u0 and u1:
    // lo = sum of chords (8 per span), hi = sum of Bezier control-polygon lengths.
    ArcLengthBounds lengthBounds(double u0, double u1) const;

    // Stitch whose half-open range contains [u0,u1] entirely, or -1.
    int stitchContaining(double u0, double u1) const;
    // The topology id of entry spanIndex of `stitches` (-1: none).
    int stitchId(int spanIndex) const { return spanIndex < 0 ? -1 : stitches[spanIndex].id; }

    // What finalize() derives: saved by the knit cache, written by the GPU.
    struct Tables {
        std::vector<double> tu, ts;  // table nodes (u, s)
        // lengthBounds' prefix sums (compensated) of the per-span chord-sum lower and
        // control-polygon upper bounds.
        std::vector<double> preLo, preHi;
        double length = 0, lengthErr = 0;
        // stitchContaining: stitches sorted by uBegin and non-overlapping at finalize() time (and
        // still the same count) are found by a binary search; otherwise by the linear scan, which is
        // correct for any spans. No knitting path makes unsorted or empty spans today; a set built by
        // hand, or edited after finalize(), takes the scan.
        bool stitchesSorted = false;
        size_t sortedCount = 0;
    };
    const Tables& tables() const { return tables_; }
    void setTables(Tables t) { tables_ = std::move(t); }

private:
    Tables tables_;
    ArcLengthBounds spanBounds(int k, double t0, double t1) const;
};

// Polyline sampling with a rigorous chord-error bound.
//
// On a cubic span, C'' is linear in the parameter, so max|C''| over the span is
// attained at an end. Linear interpolation over a sub-interval of length h has
// error <= h^2/8 * max|C''|. Each span is split into n equal sub-intervals with
// n = ceil(sqrt(M / (8 eps))), so every polyline point is within eps of the
// curve point at the same parameter (and vice versa). Distances between two
// sampled curves are therefore exact to +-2 eps (plus floating-point allowance).
struct Polyline {
    int yarnId = 0;
    std::vector<Vec3> p;      // vertices (on the curve)
    std::vector<double> u;    // native parameter of each vertex (material coordinate)
    std::vector<double> s;    // arc-length ESTIMATE of each vertex (display / t mapping)
    // Per segment k = [k, k+1]: enclosure of the true arc length of that curve piece.
    std::vector<double> segLo;  // chord length          (lower bound)
    std::vector<double> segHi;  // control-polygon length (upper bound)
    std::vector<double> cumLo, cumHi;  // prefix sums, size = segments + 1
    std::vector<double> segEps;        // chord-error bound of each segment [mm]
    double epsBound = 0;  // achieved bound (max over spans) [mm]

    size_t segments() const { return p.empty() ? 0 : p.size() - 1; }
};

// maxPerSpan caps the work; if the cap prevents reaching eps, epsBound reports the
// bound actually achieved.
// exec: spans and the per-vertex arc lengths in parallel, appended in span order; the prefix sums
// stay serial (identical polyline).
Polyline samplePath(const YarnPath& y, double eps, int maxPerSpan = 512, const ExecPolicy& exec = {});

// FlatPatch: the planar stockinette patch (one yarn per course, fabric frame with N
// = +z, translation-periodic). Only there do the front-projection checks, generator
// periodicity, cell reuse, plates and patch metrics apply. Surface: a knitted piece
// on a curved base surface; structure comes from `topo` only.
enum class CurveKind { FlatPatch, Surface };

// A set of yarn paths with the knitting structure they realize. Downstream stages
// (tubes, UVs, plies, fibres) consume a CurveSet and never depend on how it was made.
struct CurveSet {
    std::vector<YarnPath> yarns;
    CurveKind kind = CurveKind::FlatPatch;
    KnitTopology topo;  // filled whenever a yarn has stitches (a set built by hand: fillFlatTopology)
    int wales = 0, courses = 0;
    double waleSpacing = 0, courseSpacing = 0;
    double extent() const;  // characteristic size L_ref [mm] for allowances
};

// A set built by hand: fills cs.topo from the StitchSpans, stitch (wale i, course r) drawn through
// (i, r - 1), records ordered by (course, wale), and sets each span's id.
void fillFlatTopology(CurveSet& cs);

} // namespace knit
