// The yarn settle: conversion, checks and the radius continuation on the knitting solver (yarn_settle.h).
#include "stitchmesh/yarn_settle.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <exception>
#include <memory>

#include "knitting/layout/base_surface.h"
#include "knitting/relax/settle.h"
#include "model/mesh_query.h"

namespace knit {

namespace {

using Clock = std::chrono::steady_clock;

double msSince(Clock::time_point t0) { return std::chrono::duration<double, std::milli>(Clock::now() - t0).count(); }

void putf(std::ostream& o, const char* fmt, ...) {
    char b[4096];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(b, sizeof b, fmt, ap);
    va_end(ap);
    o << b;
}

// The strand and its ports: one open strand, every port shared by exactly one piece in and one out. A yarn without this
// is broken, not wrongly crossed: no solve can mend it, and it cannot become stitch records.
bool yarnConnected(const YarnCheck& c) {
    return c.portsUnmatched == 0 && c.portsOverused == 0 && c.openStrands == 1 && c.closedStrands == 0 && c.facesWithoutTemplate == 0;
}

int facesWronglyCrossed(const YarnCheck& c) { return c.facesChecked - c.facesCrossingsOk; }

void putDefects(std::ostream& rep, const YarnCheck& c) {
    for (size_t i = 0; i < c.defects.size() && i < 6; ++i) putf(rep, "  defect: %s\n", c.defects[i].c_str());
}

PatternParams paramsAt(double D, double c) {
    PatternParams pp;
    pp.yarnDiameter = D;
    pp.clearance = c;
    return pp;
}

// certifiedRepair's start diameter (knitting/relax/settle.cpp): folds pass at D <= 2 rho_lo, every nonlocal distance d at
// D <= d - c - deltaD.
double startDiameter(const SettleGate& atTarget, double target, double c) {
    const ClearanceResult& cl = atTarget.report.clearance;
    const double allowance = c + cl.deltaD;
    double D = target;
    for (const CurvatureResult& cr : atTarget.report.curvature) D = std::min(D, 2.0 * cr.rhoLo);
    if (cl.minInterYarn < 1e299) D = std::min(D, cl.minInterYarn - allowance);
    if (cl.minSelfNonlocal < 1e299) D = std::min(D, cl.minSelfNonlocal - allowance);
    for (const ClearanceRegion& g : cl.regions)
        if (g.status != Status::Pass) D = std::min(D, g.worstDistance - allowance);
    return D * (1 - 1e-6);
}

// The diagnostic pass-through: the first head of `face` mirrored through the reference (its offset negated). The certificate
// of that motion must refuse it; the line saying so goes to the report.
void mirrorHead(const YarnEmbedding& e, const YarnConversion& cv, const TriMesh& reference, int face, double D0Min, double c, CurveSet& cur,
                std::ostream& report) {
    TriMesh rr = reference;
    const MeshQuery q(rr);
    for (size_t i = 0; i < e.pieces.size(); ++i) {
        const YarnPiece& p = e.pieces[i];
        if (p.face != face || p.kind != YarnPieceKind::Head) continue;
        BSpline& b = cur.yarns[cv.strandOfPiece[i]].curve;
        const int c0 = static_cast<int>(std::ceil(cv.u0[i])) + 2, c1 = static_cast<int>(std::floor(cv.u1[i]));
        for (int k = c0; k <= c1; ++k) {
            const MeshQuery::Hit hit = q.closest(b.ctrl[k]);
            b.ctrl[k] = hit.p * 2.0 - b.ctrl[k];
        }
        putf(report, "negative control: controls %d..%d of face %d's head %d mirrored through the reference before the solve\n", c0, c1, face, p.index);
        break;
    }
    for (YarnPath& y : cur.yarns) y.finalize(1e-6);
    const SettleGate a0 = settleGateAt(cv.curves, D0Min, c), a1 = settleGateAt(cur, D0Min, c);
    const MotionCertificate mc = certifyStraightMotion(cv.curves, a0.report, cur, a1.report, paramsAt(D0Min, c));
    putf(report, "negative control: the mirroring as a straight motion at D %.3f: %s -> %s\n", D0Min, mc.detail.c_str(), mc.noPassThrough ? "certified (WRONG)" : "not certified");
}

} // namespace

std::string SettleGate::line() const {
    char b[300];
    std::snprintf(b, sizeof b, "D %.3f: regular %s, no fold (rho_lo %.4f vs r %.4f) %s, clearance (min nonlocal %.4f vs D + c) %s (%d FAIL, %d UNRESOLVED regions)", D,
                  toString(regular), rhoLo, 0.5 * D, toString(fold), minDistance, toString(clearance), failRegions, unresolvedRegions);
    return b;
}

SettleGate settleGateAt(const CurveSet& cs, double D, double c) {
    SettleGate g;
    g.D = D;
    ExecPolicy ex;
    ex.parallel = true;
    g.report = validate(cs, paramsAt(D, c), {.generatorChecks = false, .exec = ex});
    g.regular = Status::Pass;
    g.fold = Status::Pass;
    g.rhoLo = 1e300;
    for (const CurvatureResult& cr : g.report.curvature) {
        g.regular = combine(g.regular, cr.regular);
        g.fold = combine(g.fold, cr.noFold);
        g.rhoLo = std::min(g.rhoLo, cr.rhoLo);
    }
    g.clearance = g.report.clearance.status;
    g.minDistance = std::min(g.report.clearance.minSelfNonlocal, g.report.clearance.minInterYarn);
    for (const ClearanceRegion& r : g.report.clearance.regions) {
        ++g.regions;
        if (r.status == Status::Fail) ++g.failRegions;
        else if (r.status == Status::Unresolved) ++g.unresolvedRegions;
    }
    return g;
}

SettleDrift settleDrift(const CurveSet& cur, const CurveSet& rest, const TriMesh& reference) {
    SettleDrift d;
    TriMesh r = reference;
    if (r.vnormal.size() != r.pos.size()) r.computeNormals();
    const MeshQuery q(r);
    std::vector<double> moves;
    double sumAbs = 0, sum = 0, s2 = 0;
    long long n = 0, ne = 0;
    for (size_t y = 0; y < cur.yarns.size(); ++y) {
        const BSpline& A = cur.yarns[y].curve;
        const BSpline& B = rest.yarns[y].curve;
        d.length += cur.yarns[y].length();
        d.lengthRest += rest.yarns[y].length();
        for (size_t c = 0; c < A.ctrl.size(); ++c) d.maxControlMove = std::max(d.maxControlMove, norm(A.ctrl[c] - B.ctrl[c]));
        for (size_t c = 0; c + 1 < A.ctrl.size(); ++c) {
            const double l0 = norm(B.ctrl[c + 1] - B.ctrl[c]), l1 = norm(A.ctrl[c + 1] - A.ctrl[c]);
            if (l0 <= 0) continue;
            const double e = (l1 - l0) / l0;
            d.strainMax = std::max(d.strainMax, std::fabs(e));
            s2 += e * e;
            ++ne;
        }
        for (int k = 0; k <= 4 * A.numSpans(); ++k) {
            const double u = 0.25 * k;
            const Vec3 x = A.eval(u);
            moves.push_back(norm(x - B.eval(u)));
            const MeshQuery::Hit h = q.closest(x);
            if (h.tri < 0) continue;
            const std::array<int, 3>& t = r.tri[h.tri];
            const TrianglePoint tp = closestOnTriangle(x, r.pos[t[0]], r.pos[t[1]], r.pos[t[2]]);
            Vec3 nn{0, 0, 0};
            for (int m = 0; m < 3; ++m) nn = nn + r.vnormal[t[m]] * tp.bary[m];
            const double z = dot(x - tp.p, normalized(nn));
            sumAbs += std::fabs(z);
            sum += z;
            d.maxAbsDepth = std::max(d.maxAbsDepth, std::fabs(z));
            ++n;
        }
    }
    if (n) {
        d.meanAbsDepth = sumAbs / n;
        d.meanDepth = sum / n;
    }
    if (ne) d.strainRms = std::sqrt(s2 / ne);
    if (!moves.empty()) {
        std::sort(moves.begin(), moves.end());
        d.maxPointMove = moves.back();
        d.p95PointMove = moves[static_cast<size_t>(0.95 * (moves.size() - 1))];
    }
    return d;
}

const char* toString(FailedConstruction f) {
    switch (f) {
    case FailedConstruction::Refuse: return "refused";
    case FailedConstruction::ReturnUnsettled: return "returned unsettled (no relax call)";
    case FailedConstruction::SettleAnyway: return "settled anyway (diagnostic opt-in)";
    }
    return "?";
}

const char* toString(SettleStop s) {
    switch (s) {
    case SettleStop::GatePassed: return "gate passed";
    case SettleStop::ConstructionLost: return "construction lost";
    case SettleStop::Stalled: return "stalled";
    case SettleStop::MaxRepeats: return "repeat limit";
    case SettleStop::Budget: return "budget";
    case SettleStop::CallRefused: return "call refused";
    }
    return "";
}

const char* toString(SettleEnd e) {
    switch (e) {
    case SettleEnd::RawConstructionFailed: return "refused: raw construction check failed";
    case SettleEnd::YarnDisconnected: return "refused: strand or ports broken";
    case SettleEnd::ConversionRefused: return "refused: conversion";
    case SettleEnd::ConvertedConstructionFailed: return "refused: converted construction check failed";
    case SettleEnd::ConvertedOnly: return "converted only";
    case SettleEnd::ConstructionFailedNotSettled: return "not settled: construction check failed";
    case SettleEnd::Solved: return "solved";
    }
    return "?";
}

AfterTargetCall afterTargetCall(const SettleGate& gate, double D, double c, bool constructionKept, int stallCalls, TargetProgress& progress) {
    if (gate.valid()) return constructionKept ? AfterTargetCall::GatePassed : AfterTargetCall::ConstructionLost;
    const int open = gate.failRegions + gate.unresolvedRegions;
    const double lacking = std::max(0.0, 0.5 * D - gate.rhoLo) + std::max(0.0, D + c - gate.minDistance);
    // a gain: fewer open regions, or as many with a shortfall smaller by more than a micrometre
    if (open < progress.openRegions || (open == progress.openRegions && lacking < progress.shortfall - 1e-3)) {
        progress.openRegions = open;
        progress.shortfall = lacking;
        progress.callsWithoutGain = 0;
        return AfterTargetCall::Continue;
    }
    return ++progress.callsWithoutGain >= stallCalls ? AfterTargetCall::Stalled : AfterTargetCall::Continue;
}

namespace {
// The yarn of a state re-embedded and checked: what the settle's stop rule and its final verdict both judge.
struct Construction {
    YarnEmbedding embedding;
    YarnCheck check;
    int letterSwapped = 0;      // faces not as required in either frame
    int projectionChanged = 0;  // faces as required in one frame only
    bool kept = false;
};
} // namespace

SettleResult settleYarn(const SettleInput& in, const SettleSettings& ss, KnitNumerics& num, const SettleHooks& hooks, std::ostream& rep, std::ostream& log) {
    const TypedStitchMesh& s = in.mesh;
    const std::vector<std::vector<int>>& faces = in.faces;
    const TriMesh& ref = in.reference;
    const YarnSettings& st = in.yarn;
    const YarnEmbedding& e = in.raw;
    SettleResult R;
    R.raw = checkYarn(s, faces, ref, e, st);
    putf(rep, "raw: %s", R.raw.text().c_str());
    putf(rep, "templates: %s\n", st.templateSet.c_str());
    putf(rep, "a failed construction: %s\n", toString(ss.failedConstruction));
    R.rawConstructionFailed = R.raw.status != Status::Pass;
    if (R.rawConstructionFailed && ss.failedConstruction == FailedConstruction::Refuse) {
        R.end = SettleEnd::RawConstructionFailed;
        R.verdict = "refused (raw construction check FAIL: the solver is not asked to repair a construction)";
        putf(rep, "%s\n", R.verdict.c_str());
        putDefects(rep, R.raw);
        return R;
    }
    if (!yarnConnected(R.raw)) {
        R.end = SettleEnd::YarnDisconnected;
        R.verdict = "refused (the yarn's strand or ports are broken: " + std::to_string(R.raw.openStrands) + " open and " + std::to_string(R.raw.closedStrands) +
                    " closed strands, " + std::to_string(R.raw.portsUnmatched) + " unmatched and " + std::to_string(R.raw.portsOverused) + " over-used ports)";
        putf(rep, "%s\n", R.verdict.c_str());
        putDefects(rep, R.raw);
        return R;
    }
    if (R.rawConstructionFailed) {
        putf(rep, "raw construction check FAIL: %d of %d faces wrongly crossed\n", facesWronglyCrossed(R.raw), R.raw.facesChecked);
        putDefects(rep, R.raw);
    }
    // conversion
    std::string err;
    if (!yarnToCurves(e, st, ss.spacing, ss.densePerSpan, R.conversion, &err)) {
        R.end = SettleEnd::ConversionRefused;
        R.verdict = "refused (" + err + ")";
        putf(rep, "%s\n", R.verdict.c_str());
        return R;
    }
    const YarnConversion& cv = R.conversion;
    R.convertedEmbedding = embeddingFromCurves(e, cv, cv.curves, ref, st.diameter);
    R.curves = cv.curves;
    R.converted = checkYarn(s, faces, ref, R.convertedEmbedding, st);
    putf(rep, "%sconverted: %s", cv.text().c_str(), R.converted.text().c_str());
    R.convertedConstructionFailed = R.converted.status != Status::Pass;
    R.conversionAddedDefects = R.convertedConstructionFailed &&
                               (!R.rawConstructionFailed || facesWronglyCrossed(R.converted) > facesWronglyCrossed(R.raw) ||
                                R.converted.interlocksVerified < R.raw.interlocksVerified || !yarnConnected(R.converted));
    if (R.conversionAddedDefects) {
        putf(rep, "the conversion added defects: wrongly crossed faces %d -> %d, interlocks %d -> %d of %d\n", facesWronglyCrossed(R.raw), facesWronglyCrossed(R.converted),
             R.raw.interlocksVerified, R.converted.interlocksVerified, R.converted.interlocksExpected);
        putDefects(rep, R.converted);
    }
    if (R.convertedConstructionFailed && !R.rawConstructionFailed && ss.failedConstruction == FailedConstruction::Refuse) {
        R.end = SettleEnd::ConvertedConstructionFailed;
        R.verdict = "refused (the conversion changed the construction)";
        putf(rep, "%s\n", R.verdict.c_str());
        return R;
    }
    if (ss.convertOnly) {
        R.end = SettleEnd::ConvertedOnly;
        R.verdict = std::string("converted only, not settled") + (R.constructionFailed() ? " (construction check FAIL)" : "");
        putf(rep, "%s\n", R.verdict.c_str());
        return R;
    }
    if (R.constructionFailed() && ss.failedConstruction == FailedConstruction::ReturnUnsettled) {
        R.end = SettleEnd::ConstructionFailedNotSettled;
        R.verdict = std::string("not settled: the ") + (R.rawConstructionFailed ? "raw" : "converted") + " construction check FAILS (" +
                    std::to_string(facesWronglyCrossed(R.rawConstructionFailed ? R.raw : R.converted)) +
                    " faces wrongly crossed), so the settle was skipped (the solver is not asked to repair a crossing order)";
        putf(rep, "%s\n", R.verdict.c_str());
        return R;
    }
    if (R.constructionFailed()) {   // FailedConstruction::SettleAnyway
        R.settledFailedConstruction = true;
        putf(rep, "diagnostic settle of a failed construction (opt-in): settled as it is and marked; the construction cannot pass\n");
    }
    R.end = SettleEnd::Solved;
    const double D = st.diameter, c = ss.clearance;
    CurveSet& cur = R.curves;
    if (ss.perturbUlps > 0) {  // the sensitivity control: the same start moved by the last bits only
        for (YarnPath& y : cur.yarns) {
            for (Vec3& q : y.curve.ctrl)
                for (int k = 0; k < ss.perturbUlps; ++k) {
                    q.x = std::nextafter(q.x, 1e300);
                    q.y = std::nextafter(q.y, 1e300);
                    q.z = std::nextafter(q.z, 1e300);
                }
        }
        putf(rep, "sensitivity control: every converted control coordinate moved %d ulp(s) up before the solve\n", ss.perturbUlps);
    }
    if (ss.mirrorHeadOfFace >= 0) mirrorHead(e, cv, ref, ss.mirrorHeadOfFace, ss.D0Min, c, cur, rep);
    R.gateRaw = settleGateAt(cur, D, c);
    putf(rep, "production checks before the solve (local signature: not applicable, no stitch records): %s\n", R.gateRaw.line().c_str());
    // start diameter
    R.D0 = startDiameter(R.gateRaw, D, c);
    double Dprev = std::max(R.D0, ss.D0Min);
    SettleGate g0 = settleGateAt(cur, Dprev, c);
    R.startVerified = R.D0 >= ss.D0Min && g0.valid();
    putf(rep, "start diameter D0 %.4f mm (%s): %s\n", R.D0, R.startVerified ? "verified" : "NOT verified: the chain cannot be certified from the start", g0.line().c_str());
    // the stages
    const int K = Dprev >= D ? 0 : static_cast<int>(std::ceil(std::log(D / Dprev) / std::log(ss.growth)));
    std::vector<double> schedule;
    for (int k = 1; k <= K; ++k) schedule.push_back(k == K ? D : Dprev * std::pow(D / Dprev, static_cast<double>(k) / K));
    for (int k = 0; k < ss.maxRepeats; ++k) schedule.push_back(D);
    std::unique_ptr<MeshSurface> surface = ss.band ? num.meshSurface(ref) : nullptr;
    ControlMask freeCtrl(cur.yarns.size());
    for (size_t y = 0; y < cur.yarns.size(); ++y) {
        const int n = static_cast<int>(cur.yarns[y].curve.ctrl.size());
        freeCtrl[y].assign(n, 1);
        for (int k = 0; k < std::min(ss.pinnedEnd, n); ++k) freeCtrl[y][k] = freeCtrl[y][n - 1 - k] = 0;
    }
    const Controls Prest = controlsOf(cv.curves);
    bool chain = R.startVerified;
    SettleGate prevGate = g0;
    // The stop rule at the target diameter (SettleStop). The construction of a state: its yarn re-embedded and checked.
    auto constructionOf = [&](const CurveSet& state) {
        Construction made;
        made.embedding = embeddingFromCurves(e, cv, state, ref, D);
        made.check = checkYarn(s, faces, ref, made.embedding, st);
        for (const YarnFaceCheck& fc : made.check.faces) {
            if (!fc.crossingsOk && !fc.ownPlaneOk) ++made.letterSwapped;
            if (fc.crossingsOk != fc.ownPlaneOk) ++made.projectionChanged;
        }
        const YarnCheck& sc = made.check;
        made.kept = sc.portsUnmatched == 0 && sc.portsOverused == 0 && sc.openStrands == 1 && sc.closedStrands == 0 && sc.maxPortGap < 1e-6 && sc.maxTangentJump < 1.0 &&
                    made.letterSwapped == 0 && !R.constructionFailed();   // a construction wrong before the motion is never "kept"
        return made;
    };
    TargetProgress progress;
    R.stop = SettleStop::MaxRepeats;
    const Clock::time_point tSolve = Clock::now();
    for (size_t k = 0; k < schedule.size(); ++k) {
        if (msSince(tSolve) > ss.budgetMinutes * 60000.0) {
            R.budgetHit = true;
            R.stop = SettleStop::Budget;
            putf(rep, "budget: %.0f min of solve time spent; calls %zu..%zu not started (UNRESOLVED: budget)\n", ss.budgetMinutes, k + 1, schedule.size());
            break;
        }
        SettleStage row;
        row.k = static_cast<int>(k) + 1;
        row.Dprev = Dprev;
        row.D = schedule[k];
        if (hooks.onCall) hooks.onCall(row.k, static_cast<int>(schedule.size()), row.D);
        const PatternParams pk = paramsAt(row.D, c);
        RelaxMoveParams mp;
        mp.thr = row.D + c + ss.contactMargin;
        mp.kappaMax = 1.0 / (ss.bendFactor * (0.5 * row.D + 0.5 * c));
        mp.window = 0.5 * kPi * (row.D + c) * 1.3;
        mp.iterations = ss.iterations;
        mp.refresh = ss.refresh;
        mp.samplesPerSpan = ss.samplesPerSpan;
        mp.lengthStiffness = ss.lengthStiffness;
        mp.bandHalf = ss.bandHalf;
        mp.maxDisplacement = ss.displacementFraction * (row.Dprev + c);
        mp.exec.parallel = true;
        mp.call = {"repair", 0, row.k, pk};
        row.cap = mp.maxDisplacement;
        Controls P = controlsOf(cur);
        const CurveSet before = cur;
        const Clock::time_point t0c = Clock::now();
        try {
            num.relax({P, Prest, freeCtrl, cur, surface.get(), mp, nullptr});
        } catch (const std::exception& ex) {  // the device cap's refusal (DeviceBudgetExceeded) or a device failure; Cancelled passes
            row.solver = std::string("relax threw: ") + ex.what();
            if (R.refusal.empty()) R.refusal = ex.what();
            ++R.refusedCalls;
            R.stop = SettleStop::CallRefused;
            R.stages.push_back(row);
            putf(rep, "stage %d: %s\n", row.k, row.solver.c_str());
            break;
        }
        ExecPolicy ex;
        ex.parallel = true;
        num.store(cur, P, ex);
        row.ms = msSince(t0c);
        row.solver = hooks.lastSolverCall ? hooks.lastSolverCall() : "";
        row.gate = settleGateAt(cur, row.D, c);
        // the motion's certificate at D_{k-1}: both end states must pass there
        const SettleGate atPrev = row.D == row.Dprev ? row.gate : settleGateAt(cur, row.Dprev, c);
        row.startValid = prevGate.valid();
        const MotionCertificate mc = certifyStraightMotion(before, prevGate.report, cur, atPrev.report, paramsAt(row.Dprev, c));
        row.certified = mc.noPassThrough;
        row.delta = mc.maxDisplacement;
        row.localScale = mc.localScale;
        row.certificate = mc.detail;
        chain = chain && row.certified;
        putf(rep, "stage %d: D %.4f (from %.4f), cap %.4f mm, %.0f ms\n  solver: %s\n  after: %s\n  motion: %s -> %s\n", row.k, row.D, row.Dprev, row.cap, row.ms, row.solver.c_str(),
             row.gate.line().c_str(), mc.detail.c_str(), row.certified ? "certified (no centreline passed through another)" : "NOT certified");
        log << "  " << in.name << " stage " << row.k << ": D " << row.D << ", " << (row.gate.valid() ? "valid" : "not valid") << ", " << (row.certified ? "certified" : "not certified")
            << ", " << static_cast<int>(row.ms) << " ms" << std::endl;
        R.stages.push_back(row);
        prevGate = row.gate;
        Dprev = row.D;
        if (row.D != D) continue;
        // at the target: done when the gate is valid and the construction is still the generated one; a valid gate on a
        // changed construction ends the run as well (no further call can restore a crossing)
        const bool kept = row.gate.valid() && constructionOf(cur).kept;
        const AfterTargetCall after = afterTargetCall(row.gate, D, c, kept, ss.stallCalls, progress);
        if (after == AfterTargetCall::Continue) continue;
        R.stop = after == AfterTargetCall::GatePassed ? SettleStop::GatePassed : after == AfterTargetCall::ConstructionLost ? SettleStop::ConstructionLost : SettleStop::Stalled;
        if (after == AfterTargetCall::Stalled)
            putf(rep, "stop after call %d: %d calls at the target without fewer open regions or a smaller shortfall (UNRESOLVED: stalled)\n", row.k, progress.callsWithoutGain);
        else
            putf(rep, "stop after call %d: the target gate is valid and the construction is %s\n", row.k, kept ? "kept" : "NOT kept (a crossing was lost on the way)");
        break;
    }
    R.solveMs = msSince(tSolve);
    for (const SettleStage& r : R.stages) R.certifiedStages += r.certified;
    R.chainCertified = chain && !R.stages.empty() && (R.stop == SettleStop::GatePassed || R.stop == SettleStop::MaxRepeats);   // every call that ran, and a complete run
    // the settled state, judged independently of the solver's status
    R.gateFinal = R.stages.empty() ? R.gateRaw : R.stages.back().gate.D == D ? R.stages.back().gate : settleGateAt(cur, D, c);  // a refused call's gate is empty (D 0)
    const Construction settled = constructionOf(cur);
    R.settledEmbedding = settled.embedding;
    R.settled = settled.check;
    R.drift = settleDrift(cur, cv.curves, ref);
    // The construction after the motion: the strand, its ends and ports, and per face the required crossings, judged in the
    // reference frame (the construction's) and in the yarn's own best-fit plane. A face failing in both frames counts as a
    // changed construction; a face failing in one frame only is listed for inspection.
    R.letterSwapped = settled.letterSwapped;
    R.projectionChanged = settled.projectionChanged;
    R.constructionKept = settled.kept;
    for (const YarnFaceCheck& fc : R.settled.faces)
        if (!fc.crossingsOk || !fc.ownPlaneOk)
            putf(rep, "face %d (%s): expected %s; reference frame %s (%s); own plane (tilt %.1f deg) %s (%s)\n", fc.face, fc.templ.c_str(), fc.expected.c_str(), fc.found.c_str(),
                 fc.crossingsOk ? "as required" : "NOT as required", fc.ownPlaneTilt, fc.foundOwnPlane.c_str(), fc.ownPlaneOk ? "as required" : "NOT as required");
    const SettleDrift d0 = settleDrift(cv.curves, cv.curves, ref);
    putf(rep, "settled: %s", R.settled.text().c_str());
    putf(rep, "production checks at the target after the solve: %s\n", R.gateFinal.line().c_str());
    putf(rep,
         "drift from the converted curves: controls max %.4f mm; curve points at equal parameter max %.4f mm, p95 %.4f mm; length %.2f -> %.2f mm (%+.3f %%); control "
         "edge strain max %.4f, rms %.4f; offset from the reference: mean |z| %.4f -> %.4f mm, mean z %+.4f -> %+.4f mm, max |z| %.4f -> %.4f mm\n",
         R.drift.maxControlMove, R.drift.maxPointMove, R.drift.p95PointMove, R.drift.lengthRest, R.drift.length, 100.0 * (R.drift.length / R.drift.lengthRest - 1), R.drift.strainMax,
         R.drift.strainRms, d0.meanAbsDepth, R.drift.meanAbsDepth, d0.meanDepth, R.drift.meanDepth, d0.maxAbsDepth, R.drift.maxAbsDepth);
    putf(rep, "motion chain: %d of %zu calls certified, start %s: %s\n", R.certifiedStages, R.stages.size(), R.startVerified ? "verified" : "not verified",
         R.chainCertified ? "CERTIFIED from the converted state to the final state" : "NOT certified end to end");
    R.verdict = std::string(R.settledFailedConstruction ? "diagnostic settle of a failed construction (construction FAIL before the solve): " : "") +
                "settled: target " + (R.gateFinal.valid() ? "PASS" : "FAIL") + ", construction " + (R.constructionKept ? "kept" : "UNRESOLVED") +
                (R.letterSwapped ? " (" + std::to_string(R.letterSwapped) + " faces not as required in either frame: a pass-through is possible there)" : std::string("")) +
                (R.projectionChanged ? " (" + std::to_string(R.projectionChanged) + " faces as required in one frame only)" : std::string("")) + ", motion " +
                (R.chainCertified ? "certified" : "not certified") + (R.budgetHit ? ", budget hit" : "") + (R.stop == SettleStop::Stalled ? ", stalled" : "") +
                (R.refusedCalls ? ", REFUSED: " + std::to_string(R.refusedCalls) + " relax call(s) threw, the schedule did not complete" : std::string(""));
    putf(rep, "verdict: %s\n", R.verdict.c_str());
    return R;
}

} // namespace knit
