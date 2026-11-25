#include "knitting/knit/piece_knit.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "model/cancel.h"

namespace knit {

PieceSurface::PieceSurface(const Chart& c, const BaseSurface* query) : c_(c), query_(query) {
    // Meridian arc lengths over 64 meridians, 800 steps in a; their mean. The surface points in one
    // batch (a-major, as the loops below read them).
    const int NA = 800, NB = 64;
    aTab_.resize(NA + 1);
    arc_.assign(NA + 1, 0.0);
    arcB_.assign(NB, std::vector<double>(NA + 1, 0.0));
    std::vector<Vec2> ab;
    ab.reserve(size_t(NA + 1) * NB);
    for (int i = 0; i <= NA; ++i)
        for (int j = 0; j < NB; ++j) ab.push_back({double(i) / NA, (j + 0.5) / NB});
    std::vector<SurfaceSample> pts;
    if (query) query->atBatch(ab, pts);
    else PieceSurface::atBatchDirect(ab, pts);
    std::vector<Vec3> prev(NB);
    for (int i = 0; i <= NA; ++i) {
        const double a = double(i) / NA;
        aTab_[i] = a;
        double acc = 0;
        for (int j = 0; j < NB; ++j) {
            Vec3 p = pts[size_t(i) * NB + j].p;
            if (i > 0) {
                const double d = norm(p - prev[j]);
                acc += d;
                arcB_[j][i] = arcB_[j][i - 1] + d;
            }
            prev[j] = p;
        }
        if (i > 0) arc_[i] = arc_[i - 1] + acc / NB;
    }
}

namespace {

// Meridian blend: the two tabulated meridians around b and the weight of the second.
void meridianPair(double b, int nb, int& j0, int& j1, double& t) {
    double x = (b - std::floor(b)) * nb - 0.5;
    if (x < 0) x += nb;
    j0 = static_cast<int>(std::floor(x)) % nb;
    j1 = (j0 + 1) % nb;
    t = x - std::floor(x);
}

double invert(const std::vector<double>& aTab, const std::vector<double>& arc, double s) {
    if (s <= 0) return 0;
    if (s >= arc.back()) return 1;
    auto it = std::upper_bound(arc.begin(), arc.end(), s);
    const size_t i = static_cast<size_t>(it - arc.begin()) - 1;
    const double t = (s - arc[i]) / std::max(1e-300, arc[i + 1] - arc[i]);
    return aTab[i] + t * (aTab[i + 1] - aTab[i]);
}

} // namespace

double PieceSurface::aAtArc(double s, double b) const {
    int j0, j1;
    double t;
    meridianPair(b, static_cast<int>(arcB_.size()), j0, j1, t);
    return (1 - t) * invert(aTab_, arcB_[j0], s) + t * invert(aTab_, arcB_[j1], s);
}

double PieceSurface::totalArc(double b) const {
    int j0, j1;
    double t;
    meridianPair(b, static_cast<int>(arcB_.size()), j0, j1, t);
    return (1 - t) * arcB_[j0].back() + t * arcB_[j1].back();
}

SurfaceSample PieceSurface::at(double a, double b) const {
    const double a0 = c_.aMin(), a1 = c_.aMax();
    if (a >= a0 && a <= a1) {
        ChartPoint p = c_.eval(a, b);
        return {p.p, p.n};
    }
    // Pole cap: the chord from the pole to the cap boundary point, projected.
    const bool capA = a < a0;
    const Vec3 pole = capA ? c_.poleA() : c_.poleB();
    const Vec3 edge = c_.eval(capA ? a0 : a1, b).p;
    const double t = capA ? std::max(0.0, a) / a0 : std::max(0.0, 1.0 - a) / (1.0 - a1);
    ChartPoint q = c_.closestPoint(pole + (edge - pole) * t);
    return {q.p, q.n};
}

void PieceSurface::atBatchDirect(const std::vector<Vec2>& ab, std::vector<SurfaceSample>& out) const {
    out.resize(ab.size());
    for (size_t i = 0; i < ab.size(); ++i) out[i] = PieceSurface::at(ab[i].x, ab[i].y);
}

SurfaceSample PieceSurface::project(const Vec3& q) const {
    ChartPoint p = c_.closestPoint(q);
    return {p.p, p.n};
}

namespace {

// Where the rows of a piece lie. Per meridian: row rho at s_start + rho h k(b) from the start pole,
// with k(b) = (S(b) - s_start - endFromPole) / (S - s_start - endFromPole), so every meridian reaches
// its end gap at the last row. Q answers the surface queries of the placement in batches.
struct PieceRows {
    const PieceSurface& s;
    const BaseSurface& Q;
    bool startAtA;
    double startArc, endFromPole, h;
    int count = 0;  // R, the rows knitted

    double aAt(double rho, double b) const {
        const double Sb = s.totalArc(b);
        const double k = (Sb - startArc - endFromPole) / (s.totalArc() - startArc - endFromPole);
        const double sFromStart = startArc + rho * h * k;
        return s.aAtArc(startAtA ? sFromStart : Sb - sFromStart, b);
    }

    // a at arc length `distance` from the end pole along meridian b.
    double aOffEndPole(double distance, double b) const {
        const double Sb = s.totalArc(b);
        const double sFromStart = Sb - distance;
        return s.aAtArc(startAtA ? sFromStart : Sb - sFromStart, b);
    }

    // Length of the row at each rho (256 chords around b, closed), the rows in one batch.
    std::vector<double> lengths(const std::vector<double>& rhos) const {
        const int n = 256;
        std::vector<Vec2> ab;
        ab.reserve(rhos.size() * n);
        for (double rho : rhos)
            for (int j = 0; j < n; ++j) {
                const double b = double(j) / n;
                ab.push_back({aAt(rho, b), b});
            }
        std::vector<SurfaceSample> pts;
        Q.atBatch(ab, pts);
        std::vector<double> out(rhos.size());
        for (size_t q = 0; q < rhos.size(); ++q) {
            const SurfaceSample* P = pts.data() + q * n;
            double L = 0;
            Vec3 first = P[0].p, prev = first;
            for (int j = 1; j <= n; ++j) {
                Vec3 p = j == n ? first : P[j].p;
                L += norm(p - prev);
                prev = p;
            }
            out[q] = L;
        }
        return out;
    }
};

// Where the stitches of a turn lie along it.
struct TurnFractions {
    static constexpr int kPerRow = 512;   // table entries per integer row
    static constexpr int kPerTurn = 32;   // cumulative length samples per turn
    int rows = 0;
    double dir = 1;
    std::vector<std::vector<double>> bOfF;  // per integer row, b at equal arc-length fractions
    std::vector<std::vector<double>> cumLength;  // per turn, cumulative row length along the helix

    // b at fraction f of the row at rho: the tables of the two rows around it, linear between them
    // (the helix rises one row per turn).
    double bAt(double rho, double f) const {
        const double fl = std::floor(f);
        const double fr = f - fl;  // the table covers one turn; whole turns add +-1 in b
        const double x = fr * kPerRow;
        const int q = std::min(kPerRow - 1, static_cast<int>(x));
        const double u = x - q;
        auto rowB = [&](int r) {
            r = std::clamp(r, 0, rows + 2);
            return bOfF[r][q] + u * (bOfF[r][q + 1] - bOfF[r][q]);
        };
        const double rr = std::clamp(rho, 0.0, rows + 2.0);
        const int r0 = static_cast<int>(std::floor(rr));
        const double t = rr - r0;
        return (1 - t) * rowB(r0) + t * rowB(r0 + 1) + dir * fl;
    }

    // Equal arc length along the helix within a turn: the inverse of the cumulative length, so the
    // stitch at fraction f of the turn's stitches goes where that fraction of the helical length is
    // reached.
    double warped(int r, double f) const {
        const std::vector<double>& cum = cumLength[std::clamp(r, 0, rows - 1)];
        const double target = f * cum.back();
        auto it = std::lower_bound(cum.begin(), cum.end(), target);
        const int i = std::clamp(static_cast<int>(it - cum.begin()), 1, kPerTurn);
        const double t = (target - cum[i - 1]) / std::max(1e-300, cum[i] - cum[i - 1]);
        return (i - 1 + t) / kPerTurn;
    }
};

// The stitch counts of the piece from the length of each row at its middle (rho = r + 0.5), and the
// program for them. A gathered end mirrors the start: the last turn has castCount loops.
void measureRows(const PieceRows& rows, const PieceKnitParams& kp, bool gathered, PieceKnit& out) {
    const int R = rows.count;
    std::vector<double> targets(R);
    targets[0] = kp.castCount;
    std::vector<double> rhos;
    for (int r = 1; r < R; ++r) rhos.push_back(r + 0.5);
    const std::vector<double> L = rows.lengths(rhos);
    for (int r = 1; r < R; ++r) targets[r] = L[r - 1] / kp.g.w;
    if (gathered) targets[R - 1] = kp.castCount;
    out.counts = scheduleCounts(kp.castCount, targets);
    out.prog = buildProgram(out.counts, true, true);  // a CAST ring, a CLOSE drawstring
}

// The knitting direction in b: c x v must point outwards (v along the rows, c along b).
double knittingDirection(const PieceRows& rows) {
    const double am = rows.aAt(0.5 * rows.count, 0.25), da = 1e-3;
    std::vector<SurfaceSample> q;
    rows.Q.atBatch({{am, 0.25}, {am + da, 0.25}, {am - da, 0.25}, {am, 0.25 + 1e-3}, {am, 0.25 - 1e-3}}, q);
    Vec3 n = q[0].n;
    Vec3 dpda = q[1].p - q[2].p;
    Vec3 dpdb = q[3].p - q[4].p;
    const double rowDir = rows.aAt(1.0, 0.25) >= rows.aAt(0.0, 0.25) ? 1.0 : -1.0;
    return dot(cross(dpdb, dpda * rowDir), n) > 0 ? 1.0 : -1.0;
}

// Arc-length fractions along each row: a table of b(f) per integer row (b in the knitting direction).
std::vector<std::vector<double>> rowFractionTables(const PieceRows& rows, double dir) {
    const int R = rows.count, NF = TurnFractions::kPerRow;
    std::vector<Vec2> ab;
    ab.reserve(size_t(R + 3) * (NF + 1));
    for (int r = 0; r < R + 3; ++r)
        for (int j = 0; j <= NF; ++j) {
            const double b = dir * double(j) / NF;
            ab.push_back({rows.aAt(r, b), b});
        }
    std::vector<SurfaceSample> rowPts;  // every table row's points in one batch (row-major)
    rows.Q.atBatch(ab, rowPts);
    cancelCheckpoint();
    std::vector<std::vector<double>> bOfF(R + 3);
    for (int r = 0; r < R + 3; ++r) {
        std::vector<double> bs(NF + 1), cum(NF + 1, 0.0);
        Vec3 prev;
        for (int j = 0; j <= NF; ++j) {
            const double b = dir * double(j) / NF;
            bs[j] = b;
            Vec3 p = rowPts[size_t(r) * (NF + 1) + j].p;
            if (j > 0) cum[j] = cum[j - 1] + norm(p - prev);
            prev = p;
        }
        std::vector<double>& tab = bOfF[r];
        tab.resize(NF + 1);
        for (int q = 0; q <= NF; ++q) {
            const double target = cum.back() * q / NF;
            auto it = std::lower_bound(cum.begin(), cum.end(), target);
            size_t i = std::min<size_t>(NF, static_cast<size_t>(it - cum.begin()));
            if (i == 0) { tab[q] = bs[0]; continue; }
            const double t = (target - cum[i - 1]) / std::max(1e-300, cum[i] - cum[i - 1]);
            tab[q] = bs[i - 1] + t * (bs[i] - bs[i - 1]);
        }
    }
    return bOfF;
}

// The cumulative row length within each turn: the turn runs from row r to row r + 1, whose lengths
// differ a lot near a pole.
std::vector<std::vector<double>> turnLengthTables(const PieceRows& rows) {
    const int R = rows.count, NW = TurnFractions::kPerTurn;
    std::vector<double> rhos;
    for (int r = 0; r < R; ++r)
        for (int q = 1; q <= NW; ++q) rhos.push_back(r + (q - 0.5) / NW);
    const std::vector<double> L = rows.lengths(rhos);
    std::vector<std::vector<double>> tables(R);
    for (int r = 0; r < R; ++r) {
        std::vector<double> cum(NW + 1, 0.0);
        for (int q = 1; q <= NW; ++q) cum[q] = cum[q - 1] + L[size_t(r) * NW + q - 1];
        tables[r] = cum;
    }
    return tables;
}

TurnFractions planTurns(const PieceRows& rows, double dir) {
    TurnFractions turns;
    turns.rows = rows.count;
    turns.dir = dir;
    turns.bOfF = rowFractionTables(rows, dir);
    turns.cumLength = turnLengthTables(rows);
    return turns;
}

// One stitch at row coordinate rho and turn fraction f: b from the turn's arc-length table, a from the
// row at that b.
HelixSpot spotOnTurn(const PieceRows& rows, const TurnFractions& turns, double rho, double f) {
    HelixSpot spot;
    spot.f = f;
    spot.b = turns.bAt(rho, f);
    spot.a = rows.aAt(rho, spot.b);
    return spot;
}

// Where the stitches sit on the chart: each turn's relaxed fractions (turnFractions) warped to equal arc
// length along the helix, b from the row tables and a from the rows; the virtual row on the drawstring's
// ring round the end pole (closeRing > 0) or one row on.
HelixNodes turnsOnChart(const PieceRows& rows, const TurnFractions& turns, const StitchProgram& prog, double walePull, double closeRing) {
    std::vector<double> f = turnFractions(prog, walePull);
    for (const StitchRecord& st : prog.stitches) {
        const double fl = std::floor(f[st.id]);
        f[st.id] = fl + turns.warped(st.turn, f[st.id] - fl);
    }
    HelixNodes out;
    out.aSign = rows.aAt(1.0, 0.0) >= rows.aAt(0.0, 0.0) ? 1.0 : -1.0;
    for (const StitchRecord& st : prog.stitches) out.stitches.push_back(spotOnTurn(rows, turns, st.turn + f[st.id], f[st.id]));
    // The virtual yarn neighbours one stitch before the first and after the last stitch.
    const int n0 = prog.turns.front().output, nl = prog.turns.back().output;
    const StitchRecord &first = prog.stitches.front(), &last = prog.stitches.back();
    const double fb = f[first.id] - 1.0 / n0, fa = f[last.id] + 1.0 / nl;
    out.before = spotOnTurn(rows, turns, first.turn + fb, fb);
    out.after = spotOnTurn(rows, turns, last.turn + fa, fa);
    // The virtual row: one spot per stitch of the last turn T - 1, at rho = T + f.
    const int T = static_cast<int>(prog.turns.size());
    for (int k = prog.turnStart[T - 1]; k < prog.turnStart[T - 1] + prog.turns[T - 1].output; ++k) {
        HelixSpot v;
        v.f = f[k];
        v.b = turns.bAt(T + v.f, v.f);
        v.a = closeRing > 0 ? rows.aOffEndPole(closeRing, v.b) : rows.aAt(T + v.f, v.b);
        out.virtualRow.push_back(v);
    }
    return out;
}

// Layout, virtual row and yarn on the helix; numerics stores the curves.
void placeTurns(const PieceRows& rows, const TurnFractions& turns, const PieceYarn& yarn, KnitNumerics& numerics, PieceKnit& out) {
    const PieceKnitParams& kp = yarn.kp;
    out.lay = helicalLayout(rows.Q, out.prog, kp.g, turnsOnChart(rows, turns, out.prog, kp.walePull, out.closeRing));
    out.curves = placeSpiral(rows.Q, out.lay, out.prog, yarn.pp, kp.sh);
    numerics.store(out.curves, controlsOf(out.curves), {});
}

// The terminal records: CAST at the start pole, and at the end pole the drawstring placeSpiral
// built or, without one, the open last turn. endRadius: the gathered ring's, or the end gap.
void closeLastTurn(const PieceSurface& s, bool startAtA, double endRadius, PieceKnit& out) {
    const Vec3 startPole = startAtA ? s.poleA() : s.poleB(), endPole = startAtA ? s.poleB() : s.poleA();
    TerminalRecord cast;
    cast.op = TerminalOp::Cast;
    for (int k = 0; k < out.prog.turns.front().output; ++k) cast.stitches.push_back(k);
    cast.centre = startPole;
    cast.radius = out.startArc;
    std::vector<TerminalRecord> terms{cast};
    for (TerminalRecord& tr : out.curves.topo.terminals)
        if (tr.op == TerminalOp::Close) {
            tr.centre = endPole;
            tr.radius = endRadius;
            terms.push_back(tr);
        }
    if (terms.size() == 1) {
        TerminalRecord open;
        open.op = TerminalOp::Close;
        for (int k = out.prog.turnStart.back(); k < static_cast<int>(out.prog.stitches.size()); ++k) open.stitches.push_back(k);
        open.centre = endPole;
        open.radius = endRadius;
        terms.push_back(open);
    }
    out.curves.topo.terminals = terms;
}

} // namespace

bool knitPiece(const PieceStart& start, const PieceYarn& yarn, KnitNumerics& numerics, PieceKnit& out, std::string* err) {
    const PieceSurface& s = start.surface;
    const bool startAtA = start.atPoleA;
    const PieceKnitParams& kp = yarn.kp;
    out = PieceKnit{};
    out.startAtA = startAtA;
    PatternParams pp = yarn.pp;
    pp.waleSpacing = kp.g.w;
    pp.courseSpacing = kp.g.h;
    const double S = s.totalArc(), w = kp.g.w, h = kp.g.h;
    const bool gathered = kp.end == EndClosure::Gathered;
    out.startArc = kp.castCount * w / (2 * kPi);
    // Rows h apart from startArc off the start pole. Threaded: the last row h + endGap off the end pole, whatever
    // is left over adds to that gap. Gathered: the turns end startArc off the end pole too, R turns of h' each, h'
    // the nearest step to h that leaves nothing over (the last turn rises from row R - 1 to R).
    int R = 0;
    double endFromPole = 0;
    if (gathered) {
        const double span = S - 2 * out.startArc;
        R = static_cast<int>(std::lround(span / h));
        out.rowStep = span / std::max(R, 1);
        endFromPole = out.startArc + out.rowStep;  // row R - 1
    } else {
        endFromPole = h + kp.endGap;
        R = static_cast<int>(std::floor((S - out.startArc - endFromPole) / h)) + 1;
        out.rowStep = h;
    }
    if (R < 3) {
        if (err) *err = "piece too small for the gauge";
        return false;
    }
    out.endArc = out.startArc + (R - 1) * out.rowStep;
    const PieceRows rows{s, s.placementQuery(), startAtA, out.startArc, endFromPole, out.rowStep, R};
    measureRows(rows, kp, gathered, out);
    // The gathered ring: each of the last turn's n loops is passed twice by the drawstring, so the ring holds 2 n yarn
    // widths and their clearances.
    if (gathered) out.closeRing = out.counts.back() * (2 * pp.yarnDiameter + pp.clearance) / (2 * kPi);
    out.dir = knittingDirection(rows);
    const TurnFractions turns = planTurns(rows, out.dir);
    cancelCheckpoint();
    placeTurns(rows, turns, {pp, kp}, numerics, out);
    closeLastTurn(s, startAtA, gathered ? out.closeRing : kp.endGap, out);
    char buf[256];
    if (gathered)
        std::snprintf(buf, sizeof buf, "%d turns, start %s, rows from %.1f to %.1f mm of %.1f mm meridian (step %.3f mm), end gathered to a %.2f mm ring, dir %+.0f",
                      R, startAtA ? "A" : "B", out.startArc, out.endArc, S, out.rowStep, out.closeRing, out.dir);
    else
        std::snprintf(buf, sizeof buf, "%d turns, start %s, rows from %.1f to %.1f mm of %.1f mm meridian, dir %+.0f", R,
                      startAtA ? "A" : "B", out.startArc, out.endArc, S, out.dir);
    out.detail = buf;
    return true;
}

namespace {

// Lower bound of the distance from p to a cubic Bezier piece (bounding sphere of the
// control points, subdivided), and an upper bound from sampled points.
void bezierDistance(const CubicBezier& b, const Vec3& p, double& lb, double& ub, int depth) {
    Vec3 c = (b.b[0] + b.b[1] + b.b[2] + b.b[3]) * 0.25;
    double R = 0;
    for (const Vec3& q : b.b) R = std::max(R, norm(q - c));
    ub = std::min({ub, norm(b.b[0] - p), norm(b.b[3] - p)});
    const double l = norm(p - c) - R;
    if (l >= ub) return;  // cannot improve the minimum
    if (R < 1e-3 || depth > 30) {
        lb = std::min(lb, l);
        return;
    }
    const BezierHalves halves = b.split(0.5);
    bezierDistance(halves.left, p, lb, ub, depth + 1);
    bezierDistance(halves.right, p, lb, ub, depth + 1);
}

} // namespace

double poleHoleRadius(const CurveSet& cs, const Vec3& pole, double yarnRadius) {
    double lb = 1e300, ub = 1e300;
    for (const YarnPath& y : cs.yarns)
        for (int k = 0; k < y.curve.numSpans(); ++k) bezierDistance(y.curve.span(k), pole, lb, ub, 0);
    return std::min(lb, ub) - yarnRadius;
}

double poleOpening(const CurveSet& cs, const Vec3& pole, const Vec3& normal, double yarnDiameter) {
    double nearest = 1e300;
    for (const YarnPath& y : cs.yarns)
        for (double u = 0; u <= y.uEnd(); u += 0.02) {
            const Vec3 d = y.curve.eval(u) - pole;
            const double depth = dot(d, normal);
            if (std::fabs(depth) > yarnDiameter) continue;
            nearest = std::min(nearest, norm(d - normal * depth));
        }
    return nearest - 0.5 * yarnDiameter;
}

double windingNumber(const TriMesh& m, const Vec3& p) {
    // Sum of signed solid angles (Van Oosterom & Strackee), divided by 4 pi.
    double sum = 0;
    for (const auto& t : m.tri) {
        Vec3 a = m.pos[t[0]] - p, b = m.pos[t[1]] - p, c = m.pos[t[2]] - p;
        double la = norm(a), lb = norm(b), lc = norm(c);
        double num = dot(a, cross(b, c));
        double den = la * lb * lc + dot(a, b) * lc + dot(b, c) * la + dot(c, a) * lb;
        sum += 2 * std::atan2(num, den);
    }
    return sum / (4 * kPi);
}

// Whether a point of `piece` (in its own frame) lies inside another piece of its object. The winding number is taken in
// the other piece's frame, so no mesh is copied.
bool insideAnotherPiece(const PlacedPiece& piece, const std::vector<PlacedPiece>& others, const Vec3& q) {
    for (const PlacedPiece& other : others)
        if (other.mesh && other.mesh != piece.mesh && other.object == piece.object && windingNumber(*other.mesh, q + piece.centre - other.centre) > 0.5)
            return true;
    return false;
}

bool prepareScenePiece(const ScenePieceRequest& request, KnitNumerics* placement, ScenePieceKnit& out, std::string* err) {
    if (!request.piece.mesh) {
        if (err) *err = "a spiral part without its piece";
        return false;
    }
    out.chart = std::make_unique<Chart>(*request.piece.mesh, request.rings);
    cancelCheckpoint();  // between the phases of a piece
    out.query = placement ? placement->pieceQuery(*out.chart) : nullptr;
    out.surface = std::make_unique<PieceSurface>(*out.chart, out.query.get());
    cancelCheckpoint();
    out.poleAInside = insideAnotherPiece(request.piece, request.others, out.chart->poleA());
    out.poleBInside = insideAnotherPiece(request.piece, request.others, out.chart->poleB());
    return true;
}

bool knitScenePiece(const ScenePieceRequest& request, const PieceYarn& yarn, KnitNumerics& numerics, ScenePieceKnit& out,
                    std::string* err) {
    if (!prepareScenePiece(request, &numerics, out, err)) return false;
    cancelCheckpoint();
    const bool startAtA = !(out.poleAInside && !out.poleBInside);
    // The end is gathered where it shows; inside another piece it stays threaded, the way a piece is sewn on at its
    // open joint end (Igarashi, Igarashi, Suzuki 2008, 3.2).
    PieceKnitParams kp = yarn.kp;
    kp.end = (startAtA ? out.poleBInside : out.poleAInside) ? EndClosure::Threaded : EndClosure::Gathered;
    if (!knitPiece({*out.surface, startAtA}, {yarn.pp, kp}, numerics, out.knit, err)) return false;
    // A knitter's hand, last, as knitMesh applies it after the lattice's placement; settle reads the records.
    out.hand = {};
    if (kp.hand.inUse()) out.hand = knitByHand(out.knit.curves, kp.hand, numerics, {});
    return true;
}

} // namespace knit
