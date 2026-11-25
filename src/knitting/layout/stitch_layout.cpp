// Stitch layout and yarn placement on a base surface: nodes, frames, control points (stitch_layout.h).
#include "knitting/layout/stitch_layout.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "knitting/stitch/stockinette.h"
#include "model/cancel.h"

namespace knit {

namespace {

Vec3 tangential(const Vec3& d, const Vec3& n) { return d - n * dot(d, n); }

} // namespace

std::vector<double> turnFractions(const StitchProgram& prog, double wale, int iterations) {
    std::vector<double> f(prog.stitches.size(), 0.0);
    auto wrap = [](double d) { return d - std::floor(d + 0.5); };  // to [-0.5, 0.5)
    for (size_t r = 0; r < prog.turns.size(); ++r) {
        const int s0 = prog.turnStart[r], n = prog.turns[r].output;
        if (r == 0) {
            for (int j = 0; j < n; ++j) f[s0 + j] = double(j) / n;
            continue;
        }
        // Targets from the parents: K above its parent, DEC above the midpoint of its
        // two parents, INC slots half a child stitch to either side of their parent.
        std::vector<double> t(n);
        for (int j = 0; j < n; ++j) {
            const StitchRecord& st = prog.stitches[s0 + j];
            double fp = f[st.parents[0]];
            if (st.nParents == 2) fp += 0.5 * wrap(f[st.parents[1]] - fp);
            if (st.op == StitchOp::Inc) fp += (st.incSlot == 0 ? -0.5 : 0.5) / n;
            t[j] = double(j) / n + wrap(fp - double(j) / n);
        }
        // Minimize sum (f_{j+1} - f_j - 1/n)^2 + wale * sum (f_j - t_j)^2, periodic
        // (f_n = f_0 + 1), by Gauss-Seidel from the targets.
        std::vector<double> x = t;
        for (int it = 0; it < iterations; ++it)
            for (int j = 0; j < n; ++j) {
                const double prev = j > 0 ? x[j - 1] : x[n - 1] - 1.0;
                const double next = j + 1 < n ? x[j + 1] : x[0] + 1.0;
                const double d = 1.0 / n;
                // d/dx_j: (x - prev - d) - (next - x - d) + wale (x - t) = 0
                x[j] = (prev + d + next - d + wale * t[j]) / (2.0 + wale);
            }
        for (int j = 0; j < n; ++j) f[s0 + j] = x[j];
    }
    return f;
}

namespace {

LayoutNode nodeAt(const HelixSpot& spot) {
    LayoutNode nd;
    nd.f = spot.f;
    nd.b = spot.b;
    nd.a = spot.a;
    return nd;
}

// The virtual row after the last turn on the surface, in one batch, with frames: course along the
// ring of virtual children, wale away from each one's parent in the last turn.
void addVirtualRow(const BaseSurface& s, const StitchProgram& prog, const std::vector<HelixSpot>& row, PieceLayout& lay) {
    const int T = static_cast<int>(prog.turns.size());
    const int s0 = prog.turnStart[T - 1], n = prog.turns[T - 1].output;
    lay.virtualNodes.assign(n, {});
    lay.virtualChildOf.assign(lay.nodes.size(), -1);
    for (int j = 0; j < n; ++j) {
        lay.virtualNodes[j] = nodeAt(row[j]);
        lay.virtualChildOf[s0 + j] = j;
    }
    {
        std::vector<Vec2> ab;
        for (const LayoutNode& v : lay.virtualNodes) ab.push_back({v.a, v.b});
        std::vector<SurfaceSample> q;
        s.atBatch(ab, q);
        for (int j = 0; j < n; ++j) {
            lay.virtualNodes[j].p = q[j].p;
            lay.virtualNodes[j].n = q[j].n;
        }
    }
    for (int j = 0; j < n; ++j) {
        LayoutNode& v = lay.virtualNodes[j];
        const Vec3 prev = lay.virtualNodes[(j + n - 1) % n].p, next = lay.virtualNodes[(j + 1) % n].p;
        v.c = normalized(tangential(next - prev, v.n));
        v.wL = norm(v.p - prev);
        v.wR = norm(next - v.p);
        v.w = 0.5 * (v.wL + v.wR);
        const LayoutNode& par = lay.nodes[s0 + j];
        v.v = normalized(tangential(v.p - par.p, v.n));
        v.h = norm(v.p - par.p);
    }
}

// One stitch of the plain helix: row rho at a = (rho + 1) rowHeight / length, b = f.
HelixSpot evenSpot(double rho, double f, double rowHeight, double length) { return {(rho + 1.0) * rowHeight / length, f, f}; }

} // namespace

HelixNodes evenHelix(const StitchProgram& prog, double rowHeight, double length, double walePull) {
    const std::vector<double> f = turnFractions(prog, walePull);
    HelixNodes out;
    out.aSign = (1.0 + 1.0) * rowHeight / length >= (0.0 + 1.0) * rowHeight / length ? 1.0 : -1.0;
    for (const StitchRecord& st : prog.stitches) out.stitches.push_back(evenSpot(st.turn + f[st.id], f[st.id], rowHeight, length));
    const int n0 = prog.turns.front().output, nl = prog.turns.back().output;
    const StitchRecord &first = prog.stitches.front(), &last = prog.stitches.back();
    const double fb = f[first.id] - 1.0 / n0, fa = f[last.id] + 1.0 / nl;
    out.before = evenSpot(first.turn + fb, fb, rowHeight, length);
    out.after = evenSpot(last.turn + fa, fa, rowHeight, length);
    return out;
}

PieceLayout helicalLayout(const BaseSurface& s, const StitchProgram& prog, const GaugeParams& g, const HelixNodes& nodes) {
    PieceLayout lay;
    lay.helical = true;
    lay.aSign = nodes.aSign;
    for (const HelixSpot& spot : nodes.stitches) lay.nodes.push_back(nodeAt(spot));
    lay.before = nodeAt(nodes.before);
    lay.after = nodeAt(nodes.after);
    {
        std::vector<Vec2> ab;
        ab.reserve(lay.nodes.size() + 2);
        for (const LayoutNode& nd : lay.nodes) ab.push_back({nd.a, nd.b});
        ab.push_back({lay.before.a, lay.before.b});
        ab.push_back({lay.after.a, lay.after.b});
        std::vector<SurfaceSample> q;
        s.atBatch(ab, q);
        const size_t N = lay.nodes.size();
        for (size_t k = 0; k < N; ++k) {
            lay.nodes[k].p = q[k].p;
            lay.nodes[k].n = q[k].n;
        }
        lay.before.p = q[N].p;
        lay.before.n = q[N].n;
        lay.after.p = q[N + 1].p;
        lay.after.n = q[N + 1].n;
    }
    computeFrames(s, prog, g, lay);
    if (nodes.virtualRow.empty()) return lay;
    cancelCheckpoint();  // between the two batches
    addVirtualRow(s, prog, nodes.virtualRow, lay);
    return lay;
}

LayoutCheck checkLayout(const StitchProgram& prog, const PieceLayout& lay, const GaugeParams& g, double lo, double hi) {
    LayoutCheck c;
    double worstDev = 0;
    auto ratio = [&](double v, double target, int k) {
        const double q = v / target;
        if (q < lo || q > hi) ++c.outOfBounds;
        const double dev = std::max(q, 1.0 / std::max(q, 1e-300));
        if (dev > worstDev) { worstDev = dev; c.worst = k; }
        return q;
    };
    for (const StitchRecord& st : prog.stitches) {
        const LayoutNode& nd = lay.nodes[st.id];
        ++c.cells;
        const double wq = ratio(0.5 * (nd.wL + nd.wR), g.w, st.id);
        c.minW = std::min(c.minW, wq);
        c.maxW = std::max(c.maxW, wq);
        if (st.nParents > 0) {
            const double hq = ratio(nd.h, g.h, st.id);
            c.minH = std::min(c.minH, hq);
            c.maxH = std::max(c.maxH, hq);
            const double ang = std::acos(std::clamp(dot(nd.c, nd.v), -1.0, 1.0)) * 180 / kPi;
            c.minAngle = std::min(c.minAngle, ang);
            c.maxAngle = std::max(c.maxAngle, ang);
            if (ang < 60 || ang > 120) ++c.skewed;
            if (dot(cross(nd.c, nd.v), nd.n) <= 0) ++c.inverted;
        }
        if (st.nChildren == 2) {
            const LayoutNode &c0 = lay.nodes[st.children[0]], &c1 = lay.nodes[st.children[1]];
            ratio(norm(c1.p - c0.p) + 0.5 * (c0.wL + c1.wR), 2 * g.w, st.id);
        }
    }
    c.status = c.outOfBounds || c.skewed || c.inverted ? Status::Fail : Status::Pass;
    char buf[300];
    std::snprintf(buf, sizeof buf,
                  "%d cells: width %.2f..%.2f w, height %.2f..%.2f h, course/wale angle %.0f..%.0f deg; %d out of [%.2f, "
                  "%.2f], %d skewed, %d inverted",
                  c.cells, c.minW, c.maxW, c.minH, c.maxH, c.minAngle, c.maxAngle, c.outOfBounds, lo, hi, c.skewed,
                  c.inverted);
    c.detail = buf;
    return c;
}

PieceLayout planarLayout(const StitchProgram& prog, const GaugeParams& g) {
    // Row 0 centred and uniform; later rows relaxed as in turnFractions, with free row
    // ends instead of a periodic turn (x in stitch widths).
    std::vector<double> x(prog.stitches.size(), 0.0);
    for (size_t r = 0; r < prog.turns.size(); ++r) {
        const int s0 = prog.turnStart[r], n = prog.turns[r].output;
        std::vector<double> t(n);
        for (int j = 0; j < n; ++j) {
            const StitchRecord& st = prog.stitches[s0 + j];
            if (st.nParents == 0) { t[j] = j - 0.5 * (n - 1); continue; }
            double xp = x[st.parents[0]];
            if (st.nParents == 2) xp = 0.5 * (xp + x[st.parents[1]]);
            if (st.op == StitchOp::Inc) xp += st.incSlot == 0 ? -0.5 : 0.5;
            t[j] = xp;
        }
        std::vector<double> v = t;
        const double wale = r == 0 ? 1e9 : 0.5;
        for (int it = 0; it < 400; ++it)
            for (int j = 0; j < n; ++j) {
                double num = wale * t[j], den = wale;
                if (j > 0) { num += v[j - 1] + 1.0; den += 1; }
                if (j + 1 < n) { num += v[j + 1] - 1.0; den += 1; }
                v[j] = num / den;
            }
        for (int j = 0; j < n; ++j) x[s0 + j] = v[j];
    }
    PieceLayout lay;
    for (const StitchRecord& st : prog.stitches) {
        LayoutNode nd;
        nd.p = {x[st.id] * g.w, st.turn * g.h, 0};
        nd.n = {0, 0, 1};
        lay.nodes.push_back(nd);
    }
    PlaneSurface plane(1, 1);
    computeFrames(plane, prog, g, lay);
    return lay;
}

namespace {

// Course direction and widths from the neighbours in the same turn where they exist (yarn neighbours
// of a spiral, including across the transition), else one-sided.
void courseFrames(const StitchProgram& prog, const GaugeParams& g, PieceLayout& lay) {
    const int N = static_cast<int>(lay.nodes.size());
    auto sameRow = [&](int a, int b) { return b >= 0 && b < N && prog.stitches[a].turn == prog.stitches[b].turn; };
    const bool helix = lay.helical, explicitRows = !lay.rowPrev.empty();
    for (int k = 0; k < N; ++k) {
        LayoutNode& nd = lay.nodes[k];
        Vec3 prev, next;
        double span = 2;
        if (helix) {
            prev = k > 0 ? lay.nodes[k - 1].p : lay.before.p;
            next = k + 1 < N ? lay.nodes[k + 1].p : lay.after.p;
        } else if (explicitRows) {
            const int ip = lay.rowPrev[k], in = lay.rowNext[k];
            prev = ip >= 0 ? lay.nodes[ip].p : nd.p;
            next = in >= 0 ? lay.nodes[in].p : nd.p;
            span = (ip >= 0 ? 1 : 0) + (in >= 0 ? 1 : 0);
        } else {
            const bool hp = sameRow(k, k - 1), hn = sameRow(k, k + 1);
            prev = hp ? lay.nodes[k - 1].p : nd.p;
            next = hn ? lay.nodes[k + 1].p : nd.p;
            span = (hp ? 1 : 0) + (hn ? 1 : 0);
        }
        nd.c = normalized(tangential(next - prev, nd.n));
        nd.w = span > 0 ? norm(next - prev) / span : g.w;
        // Half widths: to the boundary shared with each neighbour (one-sided rows use the other
        // side's distance).
        nd.wL = norm(nd.p - prev);
        nd.wR = norm(next - nd.p);
        if (nd.wL == 0) nd.wL = nd.wR > 0 ? nd.wR : g.w;
        if (nd.wR == 0) nd.wR = nd.wL;
    }
}

// Two surface points per parentless stitch without a wale hint, around it along the chart's
// a-direction, evaluated in one batch (in stitch order, up before down). at[k]: the stitch's first
// sample, -1 for the others.
std::vector<SurfaceSample> parentlessWaleSamples(const BaseSurface& s, const StitchProgram& prog, const PieceLayout& lay, std::vector<int>& at) {
    const int N = static_cast<int>(lay.nodes.size());
    at.assign(N, -1);
    std::vector<Vec2> ab;
    const double da = 1e-4;
    for (int k = 0; k < N; ++k) {
        if (prog.stitches[k].nParents != 0 || !lay.waleHint.empty()) continue;
        const LayoutNode& nd = lay.nodes[k];
        at[k] = static_cast<int>(ab.size());
        ab.push_back({nd.a + lay.aSign * da, nd.b});
        ab.push_back({nd.a - lay.aSign * da, nd.b});
    }
    std::vector<SurfaceSample> samples;
    if (!ab.empty()) s.atBatch(ab, samples);
    return samples;
}

// Wale direction and height from the parent(s) (INC children: their slot); a parentless stitch takes
// its wale hint or the chart's a-direction, and the row height.
void waleFrames(const BaseSurface& s, const StitchProgram& prog, const GaugeParams& g, PieceLayout& lay) {
    const int N = static_cast<int>(lay.nodes.size());
    std::vector<int> sampleAt;
    const std::vector<SurfaceSample> updn = parentlessWaleSamples(s, prog, lay, sampleAt);
    for (int k = 0; k < N; ++k) {
        LayoutNode& nd = lay.nodes[k];
        const StitchRecord& st = prog.stitches[k];
        if (st.nParents == 0) {
            const Vec3 d = lay.waleHint.empty() ? tangential(updn[sampleAt[k]].p - updn[sampleAt[k] + 1].p, nd.n) : tangential(lay.waleHint[k], nd.n);
            nd.v = norm(d) > 0 ? normalized(d) : normalized(cross(nd.n, nd.c));
            nd.h = g.h;
            continue;
        }
        Vec3 base = st.op == StitchOp::Dec ? (lay.nodes[st.parents[0]].p + lay.nodes[st.parents[1]].p) * 0.5
                                           : lay.nodes[st.parents[0]].p;
        if (st.op == StitchOp::Inc) {
            const LayoutNode& par = lay.nodes[st.parents[0]];
            base = base + par.c * ((st.incSlot == 0 ? -0.5 : 0.5) * par.w);
        }
        nd.v = normalized(tangential(nd.p - base, nd.n));
        nd.h = norm(nd.p - base);
    }
}

void copyFrame(LayoutNode& dst, const LayoutNode& src) {
    dst.c = src.c;
    dst.v = src.v;
    dst.w = src.w;
    dst.wL = src.wL;
    dst.wR = src.wR;
    dst.h = src.h;
    if (dst.n.x == 0 && dst.n.y == 0 && dst.n.z == 0) dst.n = src.n;
}

} // namespace

void computeFrames(const BaseSurface& s, const StitchProgram& prog, const GaugeParams& g, PieceLayout& lay) {
    courseFrames(prog, g, lay);
    waleFrames(s, prog, g, lay);
    copyFrame(lay.before, lay.nodes.front());
    copyFrame(lay.after, lay.nodes.back());
}

namespace {

Vec3 frameMap(const LayoutNode& nd, double X, double Y) {
    const double wl = nd.wL > 0 ? nd.wL : nd.w, wr = nd.wR > 0 ? nd.wR : nd.w;
    return nd.p + nd.c * (X * (X < 0 ? wl : wr)) + nd.v * (Y * nd.h);
}

// Weight of the children's frame for a control at template height Y (0 up to the legs, 1 from the
// head side below the depth transition, smooth in between).
double headWeight(double Y, double y0, double y1) {
    const double t = std::clamp((Y - y0) / (y1 - y0), 0.0, 1.0);
    return t * t * (3 - 2 * t);  // smoothstep: own frame at y0, children's from y1
}

// What every control of a placement reads.
struct Placement {
    const PieceLayout& lay;
    const StitchProgram& prog;
    const PatternParams& pp;
    const ShapingParams& sh;
};

// A head's X in child frame cn: on a child narrower than the nominal stitch width, the head sides
// keep the nominal ABSOLUTE gap to the shared boundary (X = +-0.5) instead of scaling it down with
// the child (no effect at or above nominal width). Only the head sides move (narrowing the whole
// half-head was measured worse: more head-head clearance failures, same shoulder curvature).
double capHeadSide(const LayoutNode& cn, double Xc, const PatternParams& pp) {
    if (std::fabs(Xc) <= 0.3) return Xc;
    const double half = Xc < 0 ? (cn.wL > 0 ? cn.wL : cn.w) : (cn.wR > 0 ? cn.wR : cn.w);
    const double side = templateHeadSideX(pp);  // the template's head-side X (stitchTemplate)
    const double lim = 0.5 - (0.5 - side) * pp.waleSpacing / half;
    return Xc < 0 ? -std::min(-Xc, lim) : std::min(Xc, lim);
}

// Where stitch k's head point (X, Y) maps in the frame of the stitch(es) drawn through it.
Vec3 headInChildFrames(const Placement& in, int k, double X, double Y, bool decParent) {
    const StitchRecord& st = in.prog.stitches[k];
    const PatternParams& pp = in.pp;
    if (st.nChildren == 2) {  // INC parent: the head stretched over both children
        const LayoutNode &c0 = in.lay.nodes[st.children[0]], &c1 = in.lay.nodes[st.children[1]];
        return X < 0 ? frameMap(c0, capHeadSide(c0, 2 * X + 0.5, pp), Y - 1) : frameMap(c1, capHeadSide(c1, 2 * X - 0.5, pp), Y - 1);
    }
    const StitchRecord& ch = in.prog.stitches[st.children[0]];
    const LayoutNode& cn = in.lay.nodes[ch.id];
    double Xc = X;
    if (decParent) Xc += ch.parents[1] == k ? in.sh.decShift : -in.sh.decShift;  // k2tog: the later parent on top
    return frameMap(cn, capHeadSide(cn, Xc, pp), Y - 1);
}

// A control's query: the point p0 to project and the depth Z (template units).
struct ControlQuery {
    Vec3 p0;
    double Z = 0;
};

// Template control tc of stitch k, through the frames of stitch_layout.h.
ControlQuery controlQuery(const Placement& in, int k, const TemplateControl& tc) {
    const PieceLayout& lay = in.lay;
    const StitchProgram& prog = in.prog;
    const ShapingParams& sh = in.sh;
    const double X = tc.X, Y = tc.Y;
    double Z = tc.Z;
    const LayoutNode& nd = lay.nodes[k];
    const StitchRecord& st = prog.stitches[k];
    Vec3 p0 = frameMap(nd, X, Y);
    const bool decParent = st.nChildren == 1 && prog.stitches[st.children[0]].op == StitchOp::Dec;
    // Weight of the children's frame. DEC parents cross each other on the way to their shared child:
    // they keep their own frame longer (blend from the neck to the head side) and are layered in
    // depth before that (legs to neck), so the crossing happens fully layered.
    const int virt = st.nChildren == 0 && k < static_cast<int>(lay.virtualChildOf.size()) ? lay.virtualChildOf[k] : -1;
    const double gamma = st.nChildren == 0 && virt < 0 ? 0.0
                         : decParent ? headWeight(Y, sh.decBlend0, sh.decBlend1)
                                     : headWeight(Y, sh.headBlend0, sh.headBlend1);
    // CAST ring: the first turn's lower part (weights from the template Y, so 0 there) is pulled up
    // into a ring at the head base.
    if (prog.castRing && st.turn == 0 && Y < 0.6) p0 = frameMap(nd, X, 0.6 + (Y - 0.6) * 0.3);
    if (gamma > 0) {
        Vec3 q;
        if (virt >= 0) {  // last turn: head into the virtual row
            const LayoutNode& cn = lay.virtualNodes[virt];
            q = frameMap(cn, capHeadSide(cn, X, in.pp), Y - 1);
        } else {
            q = headInChildFrames(in, k, X, Y, decParent);
        }
        p0 = p0 * (1 - gamma) + q * gamma;
    }
    if (decParent) {
        const bool later = prog.stitches[st.children[0]].parents[1] == k;
        Z += headWeight(Y, sh.decDepth0, sh.decDepth1) * (later ? sh.decDepth : -sh.decDepth);
    }
    if (st.op == StitchOp::Dec && Y <= 0.6) Z *= 1 + sh.decLegDepth;
    return {p0, Z};
}

// A control from its projected foot.
Vec3 controlAt(const SurfaceSample& q, double Z, const PatternParams& pp) { return q.p + q.n * (Z * pp.depthAmplitude * pp.yarnDiameter); }

// One template per face combination, made on first use.
class FaceTemplates {
public:
    explicit FaceTemplates(const PatternParams& pp) : pp_(pp), plain_(stitchTemplate(pp)) {}
    const std::vector<TemplateControl>& plain() const { return plain_; }
    const std::vector<TemplateControl>& of(const StitchFaces& f) {
        for (const auto& [key, t] : byFaces_)
            if (key == f) return t;
        byFaces_.push_back({f, stitchTemplate(pp_, f)});
        return byFaces_.back().second;
    }

private:
    const PatternParams& pp_;
    std::vector<TemplateControl> plain_;
    std::vector<std::pair<StitchFaces, std::vector<TemplateControl>>> byFaces_;
};

// Consecutive stitch ids on one yarn.
struct StitchRange {
    int first, count;
};

// How one yarn's controls are queried: they run from g = -1 to gEnd, and those beyond the stitches'
// ends use the frames of the virtual neighbours `before` and `after`.
struct YarnControls {
    StitchRange stitches;
    int gEnd;
    const LayoutNode& before;
    const LayoutNode& after;
    const std::vector<StitchFaces>* faces;  // null: the plain template
};

struct StitchControl {
    ControlQuery q;
    int stitch;  // whose node state starts its continuation
};

// Every control's query along one yarn. Run-end padding has no wale transition: it uses the end
// stitch's faces with left = own.
std::vector<StitchControl> stitchControlQueries(const Placement& in, const YarnControls& yarn, FaceTemplates& templates) {
    const int K = kControlsPerStitch, s0 = yarn.stitches.first, n = yarn.stitches.count;
    std::vector<StitchControl> out;
    for (int g = -1; g <= yarn.gEnd; ++g) {
        const int i = g >= 0 ? g / K : -1;
        const int j = g - i * K;
        const int stitch = s0 + std::clamp(i, 0, n - 1);
        const bool padding = i < 0 || i >= n;
        StitchFaces f;
        if (yarn.faces) {
            f = (*yarn.faces)[stitch];
            if (padding) f.left = f.own;
        }
        const TemplateControl& tc = yarn.faces ? templates.of(f)[j] : templates.plain()[j];
        if (padding) out.push_back({{frameMap(i < 0 ? yarn.before : yarn.after, tc.X, tc.Y), tc.Z}, stitch});
        else out.push_back({controlQuery(in, s0 + i, tc), s0 + i});
    }
    return out;
}

void fillRecord(StitchRecord& rec, const StitchRecord& src, int yarn, const StitchSpan& span, const LayoutNode& nd) {
    rec = src;
    rec.yarn = yarn;
    rec.uBegin = span.uBegin;
    rec.uEnd = span.uEnd;
    rec.origin = nd.p;
    rec.N = nd.n;
    rec.U = normalized(tangential(nd.c, nd.n));
    rec.V = cross(rec.N, rec.U);
}

// Each stitch's span on yarn index `yarn` and its topology record with the frame at its node.
void recordStitchSpans(const Placement& in, const StitchRange& stitches, int yarn, YarnPath& y, CurveSet& cs) {
    const int K = kControlsPerStitch;
    for (int i = 0; i < stitches.count; ++i) {
        const StitchRecord& src = in.prog.stitches[stitches.first + i];
        StitchSpan sp{src.index, src.turn, double(i * K), double((i + 1) * K), src.id};
        y.stitches.push_back(sp);
        fillRecord(cs.topo.stitches[src.id], src, yarn, sp, in.lay.nodes[src.id]);
    }
}

// Sizes the curve set's records and sets its stitch counts.
void startCurveSet(const StitchProgram& prog, CurveKind kind, CurveSet& cs) {
    cs.kind = kind;
    cs.topo.stitches.resize(prog.stitches.size());
    int maxN = 0;
    for (const TurnPlan& t : prog.turns) maxN = std::max(maxN, t.output);
    cs.wales = maxN;
    cs.courses = static_cast<int>(prog.turns.size());
}

constexpr int kStrandControlsPerLoop = 5;

// The CLOSE strand's queries: through every loop of the last turn, alternately back to front and
// front to back, at height closeY of the virtual child frame; then one past the last loop's side.
void addDrawstringQueries(const PieceLayout& lay, const ShapingParams& sh, int loops, std::vector<StitchControl>& out) {
    const double Xs[kStrandControlsPerLoop] = {-0.4, -0.2, 0.0, 0.2, 0.4};
    const double zs[kStrandControlsPerLoop] = {-1, -1, 0, 1, 1};
    for (int k = 0; k < loops; ++k) {
        const LayoutNode& vn = lay.virtualNodes[k];
        const double sign = k % 2 == 0 ? 1.0 : -1.0;
        for (int q = 0; q < kStrandControlsPerLoop; ++q) out.push_back({{frameMap(vn, Xs[q], sh.closeY), sign * zs[q] * sh.closeDepth}, 0});
    }
    const LayoutNode& vn = lay.virtualNodes[loops - 1];
    const double sign = (loops - 1) % 2 == 0 ? 1.0 : -1.0;
    out.push_back({{frameMap(vn, 0.65, sh.closeY), sign * sh.closeDepth}, 0});
}

// Appends the strand's controls (from query `first` on) and its tail to the yarn; returns the CLOSE
// record.
TerminalRecord addCloseDrawstring(const Placement& in, const std::vector<StitchControl>& queries, const std::vector<SurfaceSample>& feet,
                                  size_t first, YarnPath& y) {
    const PatternParams& pp = in.pp;
    const int s0 = in.prog.turnStart.back(), m = in.prog.turns.back().output;
    TerminalRecord close;
    close.op = TerminalOp::Close;
    close.yarn = 0;
    size_t next = first;
    const int c0 = static_cast<int>(y.curve.ctrl.size());  // first strand control index
    for (int k = 0; k < m; ++k) {
        for (int q = 0; q < kStrandControlsPerLoop; ++q, ++next) y.curve.ctrl.push_back(controlAt(feet[next], queries[next].q.Z, pp));
        // Control c sits at the knot u = c - 1: this loop's controls span u in [c0 + 5k - 1,
        // c0 + 5k + 3]; half a knot of margin on each side.
        close.stitches.push_back(s0 + k);
        close.segments.push_back({c0 + 5.0 * k - 1.5, c0 + 5.0 * k + 3.5});
    }
    // One more strand control past the last loop's side (the curve ends short of its last controls),
    // then two padding controls: the tail goes into the piece through the gathered ring (towards the
    // centroid of the last loops, then along -n), as a drawstring's tail is hidden inside. Straight
    // down at the last loop it passed under the previous turn's heads converging at the pole.
    y.curve.ctrl.push_back(controlAt(feet[next], queries[next].q.Z, pp));
    close.segments.back().second += 1.0;
    const size_t e = y.curve.ctrl.size();
    Vec3 hub{0, 0, 0}, nHub{0, 0, 0};
    for (int k = 0; k < m; ++k) {
        hub = hub + in.lay.virtualNodes[k].p;
        nHub = nHub + in.lay.virtualNodes[k].n;
    }
    hub = hub * (1.0 / m);
    const Vec3 nIn = -normalized(nHub);
    y.curve.ctrl.push_back(lerp(y.curve.ctrl[e - 1], hub, 0.5) + nIn * (1.0 * pp.yarnDiameter));
    y.curve.ctrl.push_back(hub + nIn * (2.0 * pp.yarnDiameter));
    close.uBegin = c0 - 1.5;
    close.uEnd = c0 + 5.0 * m - 0.5;
    return close;
}

} // namespace

void finalizeCurves(CurveSet& cs) {
    for (YarnPath& y : cs.yarns) y.finalize();
}

CurveSet placeSpiral(const BaseSurface& s, const PieceLayout& lay, const StitchProgram& prog, const PatternParams& pp,
                     const ShapingParams& sh) {
    const Placement in{lay, prog, pp, sh};
    const int K = kControlsPerStitch, n = static_cast<int>(prog.stitches.size());
    FaceTemplates templates(pp);
    CurveSet cs;
    startCurveSet(prog, CurveKind::Surface, cs);
    YarnPath y;
    y.id = 0;
    y.name = "spiral";
    const bool close = prog.closeDrawstring && !lay.virtualNodes.empty();
    const StitchRange stitches{0, n};
    // Every projected control's query first (the stitches, then the CLOSE strand); the projections
    // in one batch (BaseSurface::projectBatch), in the same order.
    std::vector<StitchControl> queries = stitchControlQueries(in, {stitches, close ? n * K : n * K + 1, lay.before, lay.after, nullptr}, templates);
    const size_t nStitchQueries = queries.size();
    if (close) addDrawstringQueries(lay, sh, prog.turns.back().output, queries);
    std::vector<Vec3> qp(queries.size());
    for (size_t q = 0; q < queries.size(); ++q) qp[q] = queries[q].q.p0;
    std::vector<SurfaceSample> feet;
    s.projectBatch(qp, feet);
    for (size_t q = 0; q < nStitchQueries; ++q) y.curve.ctrl.push_back(controlAt(feet[q], queries[q].q.Z, pp));
    TerminalRecord closeRec;
    if (close) closeRec = addCloseDrawstring(in, queries, feet, nStitchQueries, y);
    recordStitchSpans(in, stitches, 0, y, cs);
    if (close) {
        closeRec.uEnd = std::min(closeRec.uEnd, y.uEnd());
        cs.topo.terminals.push_back(closeRec);
    }
    cs.topo.buildIndex({y.uEnd()});
    cs.yarns.push_back(std::move(y));
    return cs;
}

namespace {

// Where every control of every turn lands. With a binding: the continuations, and for those not
// Owned the unrestricted projection as a position fallback. Without: the projections.
struct RowFeet {
    std::vector<ProjState> next;
    std::vector<ProjResult> result;
    std::vector<SurfaceSample> feet;
    SurfaceBinding* bind = nullptr;
    size_t consumed = 0;  // queries already turned into controls
};

RowFeet projectRowControls(const BaseSurface& s, const PieceLayout& lay, const std::vector<std::vector<StitchControl>>& rowQueries, SurfaceBinding* bind) {
    RowFeet out;
    out.bind = bind;
    std::vector<ProjState> starts;
    std::vector<Vec3> qp;
    for (const std::vector<StitchControl>& row : rowQueries)
        for (const StitchControl& x : row) {
            if (bind) starts.push_back(lay.nodeState[x.stitch]);
            qp.push_back(x.q.p0);
        }
    if (!bind) {
        s.projectBatch(qp, out.feet);
        return out;
    }
    bind->surface->continueBatch(starts, qp, bind->cp, out.next, out.result);
    std::vector<Vec3> fallback;
    for (size_t q = 0; q < qp.size(); ++q)
        if (out.result[q] != ProjResult::Owned) fallback.push_back(qp[q]);
    std::vector<SurfaceSample> fb;
    if (!fallback.empty()) s.projectBatch(fallback, fb);
    out.feet.resize(qp.size());
    size_t at = 0;
    for (size_t q = 0; q < qp.size(); ++q)
        if (out.result[q] != ProjResult::Owned) out.feet[q] = fb[at++];
    return out;
}

// Turn r's controls; with a binding each records its state, which advances on Owned only.
void appendRowControls(const Placement& in, const std::vector<StitchControl>& row, RowFeet& at, int r, YarnPath& y) {
    const PatternParams& pp = in.pp;
    SurfaceBinding* bind = at.bind;
    for (const StitchControl& x : row) {
        const size_t q = at.consumed++;
        const int c = static_cast<int>(y.curve.ctrl.size());
        if (!bind) {
            y.curve.ctrl.push_back(controlAt(at.feet[q], x.q.Z, pp));
            continue;
        }
        const ProjResult res = at.result[q];
        const bool owned = res == ProjResult::Owned;
        bind->state[r][c] = owned ? at.next[q] : in.lay.nodeState[x.stitch];
        bind->result[r][c] = res;
        y.curve.ctrl.push_back(owned ? at.next[q].foot + at.next[q].n * (x.q.Z * pp.depthAmplitude * pp.yarnDiameter)
                                     : controlAt(at.feet[q], x.q.Z, pp));
    }
}

} // namespace

CurveSet placeRows(const BaseSurface& s, const PieceLayout& lay, const StitchProgram& prog, const PatternParams& pp,
                   const RowPlacement& rows) {
    const Placement in{lay, prog, pp, rows.sh};
    SurfaceBinding* bind = rows.bind;
    FaceTemplates templates(pp);
    const int K = kControlsPerStitch, T = static_cast<int>(prog.turns.size());
    CurveSet cs;
    startCurveSet(prog, rows.kind, cs);
    if (bind) {
        bind->state.assign(prog.turns.size(), {});
        bind->result.assign(prog.turns.size(), {});
    }
    // Every control's query of every turn first; then the continuations or projections in one
    // batch, in order. Run-end padding uses the turn's end frames shifted by one width.
    std::vector<std::vector<StitchControl>> rowQueries(prog.turns.size());
    for (int r = 0; r < T; ++r) {
        const int s0 = prog.turnStart[r], n = prog.turns[r].output;
        LayoutNode pre = lay.nodes[s0], post = lay.nodes[s0 + n - 1];
        pre.p = pre.p - pre.c * pre.w;
        post.p = post.p + post.c * post.w;
        rowQueries[r] = stitchControlQueries(in, {{s0, n}, n * K + 1, pre, post, rows.faces}, templates);
    }
    RowFeet at = projectRowControls(s, lay, rowQueries, bind);
    std::vector<double> ends;
    for (int r = 0; r < T; ++r) {
        const int s0 = prog.turnStart[r], n = prog.turns[r].output;
        YarnPath y;
        y.id = r;
        y.course = r;
        y.name = "row " + std::to_string(r);
        if (bind) {
            bind->state[r].assign(n * K + 3, ProjState{});
            bind->result[r].assign(n * K + 3, ProjResult::Budget);
        }
        appendRowControls(in, rowQueries[r], at, r, y);
        recordStitchSpans(in, {s0, n}, r, y, cs);
        ends.push_back(y.uEnd());
        cs.yarns.push_back(std::move(y));
    }
    cs.topo.buildIndex(ends);
    return cs;
}

} // namespace knit
