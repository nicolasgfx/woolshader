// Settle, the local geometric post-pass of the construction, and the certified repair (settle.h).
#include "knitting/relax/settle.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <set>
#include <vector>

#include "model/cancel.h"
#include "model/log.h"

namespace knit {

namespace {

using Clock = std::chrono::steady_clock;
using Ctrl = Controls;
using Mask = ControlMask;

using MoveParams = RelaxMoveParams;

// Free controls of a stitch set grown by `ring` neighbour rings (yarn neighbours,
// parents, children).
Mask freeMask(const CurveSet& cs, const KnitTopology& topo, const std::vector<int>& seed, int ring, int* count) {
    std::set<int> freeSt(seed.begin(), seed.end());
    for (int r = 0; r < ring; ++r) {
        std::set<int> grown = freeSt;
        for (int t : freeSt) {
            const StitchRecord& st = topo.stitches[t];
            for (int d : {-1, 1}) {
                int u = t + d;
                if (u >= 0 && u < static_cast<int>(topo.stitches.size()) && topo.stitches[u].yarn == st.yarn) grown.insert(u);
            }
            for (int q = 0; q < st.nParents; ++q) grown.insert(st.parents[q]);
            for (int q = 0; q < st.nChildren; ++q) grown.insert(st.children[q]);
        }
        freeSt.swap(grown);
    }
    if (count) *count = static_cast<int>(freeSt.size());
    Mask m(cs.yarns.size());
    for (size_t y = 0; y < cs.yarns.size(); ++y) m[y].assign(cs.yarns[y].curve.ctrl.size(), 0);
    for (int t : freeSt) {
        const StitchRecord& st = topo.stitches[t];
        const int c0 = static_cast<int>(std::floor(st.uBegin)) + 1, c1 = static_cast<int>(std::floor(st.uEnd));
        for (int c = std::max(1, c0); c <= c1 && c + 1 < static_cast<int>(m[st.yarn].size()); ++c) m[st.yarn][c] = 1;
    }
    // Terminal segments (CLOSE drawstring) of free last-turn stitches move with them,
    // together with the tail behind the strand (all but the yarn's last control).
    for (const TerminalRecord& tr : topo.terminals) {
        if (tr.op != TerminalOp::Close || tr.uEnd <= tr.uBegin) continue;
        bool any = false;
        for (int t : tr.stitches) any = any || freeSt.count(t);
        if (!any) continue;
        for (int c = static_cast<int>(tr.uBegin) + 1; c + 1 < static_cast<int>(m[tr.yarn].size()); ++c)
            m[tr.yarn][c] = 1;
    }
    return m;
}

MoveParams moveParams(const PatternParams& pp, const SettleParams& sp) {
    MoveParams mp;
    const double Dc = pp.yarnDiameter + pp.clearance;
    mp.thr = Dc + sp.contactMargin;
    mp.kappaMax = 1.0 / (sp.bendFactor * (pp.radius() + 0.5 * pp.clearance));
    mp.window = 0.5 * kPi * Dc * 1.3;
    mp.iterations = sp.iterations;
    mp.refresh = sp.refresh;
    mp.samplesPerSpan = sp.samplesPerSpan;
    mp.lengthStiffness = sp.lengthStiffness;
    mp.bandHalf = sp.bandHalf;
    mp.exec = sp.exec;
    return mp;
}

} // namespace

std::vector<int> defectStitches(const CurveSet& cs, const ValidationReport& r, const PatternParams& pp) {
    const KnitTopology& topo = cs.topo;
    std::set<int> s;
    auto addRange = [&](int y, double u0, double u1) {
        for (double u = std::floor(u0); u <= u1 + 1e-9; u += 1.0) {
            int t = topo.stitchOrLoopAt(y, std::clamp(u, 0.0, cs.yarns[y].uEnd()));
            if (t < 0 && !topo.stitches.empty())  // padding or tail: blame the nearest stitch
                t = topo.stitchAt(y, std::clamp(u, topo.stitches.front().uBegin + 0.5, topo.stitches.back().uEnd - 0.5));
            if (t >= 0) s.insert(t);
        }
    };
    for (const ClearanceRegion& g : r.clearance.regions) {
        if (g.status == Status::Pass) continue;
        addRange(g.yarnA, g.uA0, g.uA1);
        addRange(g.yarnB, g.uB0, g.uB1);
    }
    const double Rc = pp.radius() + 0.5 * pp.clearance;
    for (size_t y = 0; y < r.curvature.size(); ++y)
        for (size_t k = 0; k < r.curvature[y].spanKappaHi.size(); ++k)
            if (r.curvature[y].spanKappaHi[k] * Rc > 1) {
                int t = topo.stitchOrLoopAt(static_cast<int>(y), k + 0.5);
                if (t < 0)  // padding: blame the nearest stitch
                    t = topo.stitchAt(static_cast<int>(y), std::min(k + 0.5, topo.stitches.back().uEnd - 0.5));
                if (t >= 0) s.insert(t);
            }
    for (const LocalPairResult& lp : r.localSignature.pairs)
        if (lp.status != Status::Pass) {
            s.insert(lp.child);
            for (int q = 0; q < topo.stitches[lp.child].nParents; ++q) s.insert(topo.stitches[lp.child].parents[q]);
        }
    return std::vector<int>(s.begin(), s.end());
}

namespace {

std::vector<int> everyStitch(const KnitTopology& topo) {
    std::vector<int> all(topo.stitches.size());
    for (size_t t = 0; t < all.size(); ++t) all[t] = static_cast<int>(t);
    return all;
}

} // namespace

SettleReport settle(const RelaxTarget& target, const SettleParams& sp, KnitNumerics& num, ValidationReport* last) {
    auto t0 = Clock::now();
    CurveSet& cs = target.curves;
    const PatternParams& pp = target.pp;
    SettleReport rep;
    const KnitTopology topo = cs.topo;  // copy: cs is modified
    Ctrl P = controlsOf(cs);
    const Ctrl P0 = P;
    const Ctrl Prest = target.handRest ? *target.handRest : P0;
    MoveParams mp = moveParams(pp, sp);
    ValidationReport rpt = num.validate(cs, pp, {.generatorChecks = false, .exec = sp.exec});
    const std::vector<int> defects = target.handRest ? everyStitch(topo) : defectStitches(cs, rpt, pp);
    rep.defectsPerRound.push_back(static_cast<int>(defects.size()));
    if (!defects.empty()) {
        cancelCheckpoint();
        rep.rounds = 1;
        const Mask freeCtrl = freeMask(cs, topo, defects, sp.ring, &rep.freeStitches);
        mp.call = {"settle", 0, -1, pp};
        num.relax({P, Prest, freeCtrl, cs, target.surface, mp, target.bind});
        num.store(cs, P, sp.exec);
        rpt = num.validate(cs, pp, {.generatorChecks = false, .exec = sp.exec});
    }
    for (size_t y = 0; y < P.size(); ++y)
        for (size_t c = 0; c < P[y].size(); ++c) rep.maxMove = std::max(rep.maxMove, norm(P[y][c] - P0[y][c]));
    if (last) *last = rpt;
    rep.ms = msSince(t0);
    return rep;
}

// ------------------------------------------------------- certified repair ----

namespace {

PatternParams atDiameter(const PatternParams& pp, double D) {
    PatternParams q = pp;
    q.yarnDiameter = D;
    return q;
}

bool fullyValid(const ValidationReport& r) { return r.verdict == Status::Pass; }

// The largest diameter at which the report at the target says the fixed centreline can pass: folds
// need D <= 2 rho_lo, every nonlocal distance d (the minima and each region that did not pass)
// needs d >= D + c + deltaD. Monotonicity in D makes every smaller D pass as well;
// the caller verifies the value with one validation, so a bound the report does not capture shows
// up as a failed verification, never as an unverified start.
double startDiameter(const ValidationReport& atTarget, const PatternParams& target) {
    const ClearanceResult& cl = atTarget.clearance;
    const double allowance = target.clearance + cl.deltaD;
    double D = target.yarnDiameter;
    for (const CurvatureResult& c : atTarget.curvature) D = std::min(D, 2.0 * c.rhoLo);
    if (cl.minInterYarn < 1e299) D = std::min(D, cl.minInterYarn - allowance);
    if (cl.minSelfNonlocal < 1e299) D = std::min(D, cl.minSelfNonlocal - allowance);
    for (const ClearanceRegion& g : cl.regions)
        if (g.status != Status::Pass) D = std::min(D, g.worstDistance - allowance);
    return D * (1 - 1e-6);  // strictly inside every bound despite rounding
}

// Diagnostic only: the curvature bound of the motion's midpoint state.
double midpointCurvatureRadius(const RelaxTarget& yarn, KnitNumerics& num, const Ctrl& Pi, const Ctrl& P, const PatternParams& pD,
                               const ExecPolicy& exec) {
    Ctrl Pm = Pi;
    for (size_t y = 0; y < Pm.size(); ++y)
        for (size_t c = 0; c < Pm[y].size(); ++c) Pm[y][c] = (Pi[y][c] + P[y][c]) * 0.5;
    CurveSet mid = yarn.curves;
    num.store(mid, Pm, exec);
    double rhoMid = 1e300;
    for (const YarnPath& y : mid.yarns)
        rhoMid = std::min(rhoMid, num.curvature(y, pD.radius(), Allowances::forExtent(mid.extent()), exec).rhoLo);
    return rhoMid;
}

std::string repairSummary(const RepairReport& rep) {
    const RepairMove& m = rep.move;
    char buf[320];
    std::snprintf(buf, sizeof buf,
                  "D0 %.3f, target %.3f mm: move %s (delta %.3f mm%s, free %d stitches, mid rho_lo %.3f)%s%s; "
                  "same-yarn pairs closer than %.3f mm of arc not certified; %d validations",
                  rep.D0, rep.targetD, m.accepted ? "kept" : "not kept", m.delta, m.noTubeContact ? ", no tube contact" : "",
                  m.freeStitches, m.midRhoLo, m.why.empty() ? "" : " - ", m.why.c_str(), m.localScale, rep.validations);
    return buf;
}

} // namespace

RepairReport certifiedRepair(const RelaxTarget& yarn, const RepairParams& rp, KnitNumerics& num, ValidationReport* final, const ValidationReport* current) {
    auto t0 = Clock::now();
    CurveSet& cs = yarn.curves;
    const PatternParams& target = yarn.pp;
    const ExecPolicy& exec = rp.settle.exec;
    RepairReport rep;
    rep.targetD = target.yarnDiameter;
    auto validateAt = [&](const CurveSet& curves, const PatternParams& pp) {
        ++rep.validations;
        return num.validate(curves, pp, {.generatorChecks = false, .exec = exec});
    };
    auto finish = [&](const ValidationReport& r) {
        if (final) *final = r;
        rep.reachedTarget = fullyValid(r);
        rep.ms = msSince(t0);
        return rep;
    };
    ValidationReport atTarget;
    if (current) atTarget = *current;
    else atTarget = validateAt(cs, target);
    if (fullyValid(atTarget)) {
        rep.D0 = target.yarnDiameter;
        rep.detail = "valid at the target diameter; no repair needed";
        return finish(atTarget);
    }
    // A local signature FAIL does not depend on the diameter (only its UNRESOLVED depth signs do):
    // no diameter can pass, so no start diameter exists and none is validated.
    if (atTarget.localSignature.status == Status::Fail) {
        rep.detail = "no verified start diameter: the local signature fails at the target diameter (" + std::to_string(atTarget.localSignature.fail) +
                     " stitches), and a failing local signature fails at every diameter (threading or local folds are construction errors, not repairable)";
        return finish(atTarget);
    }

    // 1. The start diameter, from the report, verified once.
    const double D0 = startDiameter(atTarget, target);
    char buf[160];
    if (!(D0 >= rp.D0Min)) {
        std::snprintf(buf, sizeof buf, "no verified start diameter: the report allows at most D = %.3f mm (below %.3f)", D0, rp.D0Min);
        rep.detail = buf;
        return finish(atTarget);
    }
    cancelCheckpoint();
    const PatternParams p0 = atDiameter(target, D0);
    const ValidationReport r0 = validateAt(cs, p0);
    if (!fullyValid(r0)) {
        std::string failing;
        for (const CheckLine& l : r0.lines)
            if (l.status != Status::Pass) failing += " [" + l.name + ": " + toString(l.status) + ", " + l.detail + "]";
        std::snprintf(buf, sizeof buf, "no verified start diameter: the construction fails at the report's D0 = %.3f mm:", D0);
        rep.detail = buf + failing;
        return finish(atTarget);
    }
    rep.D0 = D0;

    // 2. One move: hard thresholds at D0, soft targets at the target, displacement capped.
    cancelCheckpoint();
    rep.moveTried = true;
    RepairMove& m = rep.move;
    MoveParams mp = moveParams(target, rp.settle);
    mp.maxDisplacement = rp.displacementFraction * (D0 + target.clearance);
    mp.call = {"repair", -1, 0, p0};
    const KnitTopology topo = cs.topo;
    const Mask freeCtrl = freeMask(cs, topo, defectStitches(cs, atTarget, target), rp.ring, &m.freeStitches);
    const Ctrl Prest = controlsOf(cs);
    const Ctrl Pi = Prest;
    Ctrl P = Pi;
    SurfaceBinding trialBind;
    if (yarn.bind) trialBind = *yarn.bind;
    num.relax({P, Prest, freeCtrl, cs, yarn.surface, mp, yarn.bind ? &trialBind : nullptr});
    CurveSet trial = cs;
    num.store(trial, P, exec);

    // 3. Kept only when valid at the target and certified at D0.
    const ValidationReport rT = validateAt(trial, target);
    const MotionCertificate mc = certifyStraightMotion(cs, r0, trial, rT, p0);
    m.delta = mc.maxDisplacement;
    m.noPassThrough = mc.noPassThrough;
    m.noTubeContact = mc.noTubeContact;
    m.localScale = mc.localScale;
    m.midRhoLo = midpointCurvatureRadius(yarn, num, Pi, P, p0, exec);
    if (!fullyValid(rT)) m.why = std::string("gate (a) ") + toString(rT.verdict) + " at the target";
    else if (!mc.noPassThrough) m.why = "the straight motion is not certified (" + mc.detail + ")";
    m.accepted = m.why.empty();
    if (!m.accepted) {
        rep.detail = repairSummary(rep);
        return finish(atTarget);
    }
    cs = std::move(trial);
    if (yarn.bind) *yarn.bind = std::move(trialBind);
    rep.detail = repairSummary(rep);
    return finish(rT);
}

} // namespace knit
