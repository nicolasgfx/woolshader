#pragma once
// Knitting one scene piece: a spiral from pole to pole on the piece's chart, with CAST at the start
// pole and a drawstring gathering the last turn at the end pole.
//
// PieceSurface extends the chart (valid on a in [aMin, aMax]) over the two pole caps:
// for a < aMin the point is the closest mesh point to the chord from pole A to the
// chart point (aMin, b), at fraction a / aMin (and symmetrically at pole B). At the
// cap boundary the chord ends ON the surface, so the extension is continuous there,
// and b keeps its meaning around the pole.
//
// Rows are placed PER MERIDIAN: row rho lies at arc length s_start + rho h' k(b) from
// the start pole along meridian b, where k(b) scales each meridian so that all reach
// the end at the same last row (pieces that are not rotationally symmetric about
// their pole axis, e.g. the cupped ears, would otherwise get rows far apart on one
// side and crowded on the other). Turn r's target count is the length of the row at
// rho = r + 0.5 divided by w.
//   start (CAST)  turn 0 = castCount loops from s_start = castCount w / 2 pi off the start
//                 pole; their lower parts (Y < 0.6) are pulled up into a ring (no legs into
//                 the pole). Turn 1 is drawn through them as through ordinary loops.
//   end (CLOSE)   the last turn's heads are mapped into virtual children and the yarn
//                 continues as a drawstring through all of them (alternately back to
//                 front / front to back, see placeSpiral). Threaded: the virtual children
//                 one row further, endGap (and the rows' remainder) off the end pole.
//                 Gathered, mirroring the start as a knitter closes a piece with a
//                 drawstring pulled tight: the last turn ends s_start off the end pole with
//                 castCount loops (where the halving limit lets the decreases get there),
//                 the row step h' stretched so that nothing is left over, and the virtual
//                 children on a ring of radius n (2 D + c) / 2 pi round the end pole (n
//                 loops, each passed twice by the drawstring: the tightest ring the yarn
//                 fits in). The remaining opening is measured (poleOpening).
#include <memory>
#include <string>
#include <vector>

#include "knitting/layout/chart.h"
#include "knitting/layout/stitch_layout.h"
#include "knitting/numerics.h"
#include "knitting/knit/pieces.h"
#include "knitting/stitch/hand_tension.h"

namespace knit {

class PieceSurface : public BaseSurface {
public:
    // query (optional): the device's batched at / project for the placement (PieceSurface semantics),
    // which also evaluates the arc-length tables; null: this surface. It must outlive the surface.
    // at and project here stay CPU, double precision: the band of settle and the reference the
    // placement tests compare with.
    explicit PieceSurface(const Chart& c, const BaseSurface* query = nullptr);
    SurfaceSample at(double a, double b) const override;
    SurfaceSample project(const Vec3& q) const override;
    const BaseSurface& placementQuery() const { return query_ ? *query_ : *this; }
    const Chart& chart() const { return c_; }
    // Mean (over b) meridian length.
    double totalArc() const { return arc_.back(); }
    // Per meridian b (interpolated between 64 tabulated meridians): the inverse of the arc
    // length from pole A, and the meridian's total length.
    double aAtArc(double s, double b) const;
    double totalArc(double b) const;
    Vec3 poleA() const { return c_.poleA(); }
    Vec3 poleB() const { return c_.poleB(); }

private:
    // PieceSurface::at for each (a, b), in order (non-virtual: used by the constructor).
    void atBatchDirect(const std::vector<Vec2>& ab, std::vector<SurfaceSample>& out) const;

    const Chart& c_;
    const BaseSurface* query_;
    std::vector<double> aTab_, arc_;  // arc length table over a in [0, 1]
    std::vector<std::vector<double>> arcB_;  // per meridian j (b = (j + 0.5) / 64)
};

// How the drawstring closes the last turn. Threaded: through the loops where they sit, the virtual row endGap off the
// end pole, which leaves an open ring (6 to 10 mm); for an end inside another piece, never seen. Gathered: the last
// turn mirrors the cast ring and the drawstring pulls its loops into the tightest ring the yarn fits in; for an end that
// shows. Gathering crowds the last turn's loops: gate (a) finds extra crossings there.
enum class EndClosure { Threaded, Gathered };

struct PieceKnitParams {
    GaugeParams g;
    ShapingParams sh;
    int castCount = 6;    // loops of the first turn, and of a gathered end's last turn
    double endGap = 12.0; // [mm] a threaded end: the virtual row after the last turn this far off the end pole (plus the rows' remainder)
    EndClosure end = EndClosure::Threaded;  // knitScenePiece gathers an end whose pole shows
    double walePull = 0.5;  // turnFractions: wale springs vs course springs
    HandRequest hand;       // a knitter's hand (off by default): knitScenePiece applies it last
};

struct PieceKnit {
    StitchProgram prog;
    PieceLayout lay;
    CurveSet curves;
    bool startAtA = true;
    double dir = 1;             // knitting direction in b
    std::vector<int> counts;
    double startArc = 0, endArc = 0;  // arc length of the first / last row from the start pole
    double rowStep = 0;               // the row spacing: h, or h' for a gathered end
    double closeRing = 0;             // a gathered end: the drawstring's ring round the end pole (0: threaded)
    std::string detail;
};

// Where the spiral starts.
struct PieceStart {
    const PieceSurface& surface;
    bool atPoleA = true;
};
// The yarn and the gauge: pp carries D, c, A (w, h come from kp.g).
struct PieceYarn {
    const PatternParams& pp;
    const PieceKnitParams& kp;
};

// Builds program, layout and yarn for the piece, on the surface's placementQuery (measureRows,
// planTurns, placeTurns, closeLastTurn); numerics stores the placed curves. Returns false with a
// message if the piece is too small for the gauge.
bool knitPiece(const PieceStart& start, const PieceYarn& yarn, KnitNumerics& numerics, PieceKnit& out, std::string* err = nullptr);

// Hole radius at a pole: certified lower bound of the distance from the pole to the
// yarn's centreline (per span: distance to a ball around the Bezier control points,
// subdivided to 1e-3 mm), minus the yarn radius. It counts yarn anywhere, the
// drawstring's tail under the surface too.
double poleHoleRadius(const CurveSet& cs, const Vec3& pole, double yarnRadius);
// The opening seen from outside at a pole: over the yarn within one diameter of the
// surface's level at the pole (along its normal), the smallest distance from the pole
// across the surface, minus the yarn radius. Sampled every 0.02 of a span: a measurement,
// not a bound. About 1e300 when no yarn is that near the surface.
double poleOpening(const CurveSet& cs, const Vec3& pole, const Vec3& normal, double yarnDiameter);

// Winding number of a closed triangle mesh around p (about 1 inside, 0 outside).
double windingNumber(const TriMesh& m, const Vec3& p);

// One scene piece, knitted: chart, surface, pole choice (the end pole inside another piece of the same
// object where possible), the knit and the hand.
struct ScenePieceKnit {
    std::unique_ptr<Chart> chart;
    std::unique_ptr<BaseSurface> query;  // the placement's queries (null: the knit cache's rebuild)
    std::unique_ptr<PieceSurface> surface;
    PieceKnit knit;
    HandResult hand;  // the hand's records by stitch id and its rest shape, for settle (empty: no hand)
    bool poleAInside = false, poleBInside = false;
};
// A piece where it stands: its mesh in library mm about `centre`, and the object it belongs to.
struct PlacedPiece {
    const TriMesh* mesh = nullptr;  // null: no piece (a part on UV charts)
    Vec3 centre;
    int object = 0;
};
// Which piece, in memory: its mesh holds the chart as uv (a, b), `rings` the chart's domain. `others`:
// every piece of the scene, itself included (the pole tests read those of its object).
struct ScenePieceRequest {
    PlacedPiece piece;
    int rings = 0;
    const std::vector<PlacedPiece>& others;
};
// The piece before knitting: mesh, chart, surface and the pole tests (knitScenePiece's first
// step). placement: makes the placement's queries; null only for the knit cache's rebuild of a
// cached piece, which places nothing and keeps the CPU arc-length tables.
bool prepareScenePiece(const ScenePieceRequest& request, KnitNumerics* placement, ScenePieceKnit& out, std::string* err = nullptr);
bool knitScenePiece(const ScenePieceRequest& request, const PieceYarn& yarn, KnitNumerics& numerics, ScenePieceKnit& out,
                    std::string* err = nullptr);

} // namespace knit
