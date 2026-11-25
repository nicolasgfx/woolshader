// Gate (a): runs the checks of checks.h on yarn curves, combines them into a verdict and a report,
// certifies a straight motion between two valid states, and checks the flat generator and the
// patterns on flat swatches, and writes a report as exact text. The signature checks are in
// signature.cpp.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <sstream>

#include "knitting/checks/checks.h"
#include "knitting/stitch/stockinette.h"
#include "model/log.h"

namespace knit {

namespace {
using Clock = std::chrono::steady_clock;
double localParam(const YarnPath& y, int stitchIdx, double u) {
    const StitchSpan& st = y.stitches[stitchIdx];
    return y.reversed ? st.uEnd - u : u - st.uBegin;
}
} // namespace

namespace {

// Periodicity on coefficients: stitch (i, r) controls equal the controls of the first stitch
// (i0, r0) of the same face class (StitchFaces) translated by ((i - i0) w, (r - r0) h, 0).
// Stockinette has one class, reference (0, 0).
void checkPeriodicity(const CurveSet& cs, const PatternParams& pp, const Allowances& al, GeneratorChecks& out) {
    const int K = kControlsPerStitch, n = cs.wales;
    auto ctrlOf = [&](const YarnPath& y, int g) -> const Vec3& {
        // g = global (unreversed) index in [-1, nK+1]; stored index g+1, reversed if needed.
        int idx = g + 1;
        if (y.reversed) idx = static_cast<int>(y.curve.ctrl.size()) - 1 - idx;
        return y.curve.ctrl[idx];
    };
    std::vector<std::pair<StitchFaces, std::pair<int, int>>> classRef;  // faces -> (yarn, wale) of the reference
    auto refOf = [&](const StitchFaces& f) -> std::pair<int, int> {
        for (const auto& [k, v] : classRef)
            if (k == f) return v;
        for (size_t y = 0; y < cs.yarns.size(); ++y)
            for (int i = 0; i < n; ++i)
                if (flatStitchFaces(pp, i, cs.yarns[y].course) == f) {
                    classRef.push_back({f, {static_cast<int>(y), i}});
                    return classRef.back().second;
                }
        return {-1, -1};  // a virtual stitch's class with no real stitch: not compared
    };
    for (const YarnPath& y : cs.yarns) {
        for (int g = -1; g <= n * K + 1; ++g) {
            int i = g >= 0 ? g / K : -1;
            int j = g - i * K;
            const auto [ry, i0] = refOf(flatStitchFaces(pp, i, y.course));
            if (ry < 0) continue;
            const YarnPath& yr = cs.yarns[ry];
            Vec3 ref = ctrlOf(yr, i0 * K + j) + Vec3{(i - i0) * pp.waleSpacing, 0, 0} +
                       Vec3{0, (y.course - yr.course) * pp.courseSpacing, 0};
            out.maxPeriodicityDeviation = std::max(out.maxPeriodicityDeviation, norm(ctrlOf(y, g) - ref));
        }
    }
    // Rounding of (X + i - cx) * w accumulates a few ulps of the extent.
    const double tol = 64 * al.etaX;
    out.periodicity = out.maxPeriodicityDeviation <= tol ? Status::Pass : Status::Fail;
}

// One crossing of a pair: its parameters local to the two stitches, and the newer yarn's depth.
struct PairCrossing {
    double ln, lo, depth;
};

// The pair's crossings sorted along the newer yarn.
std::vector<PairCrossing> pairCrossings(const CurveSet& cs, const CrossingResult& cr, const PairSignature& ps) {
    std::vector<PairCrossing> v;
    for (int k : ps.crossingIdx) {
        const Crossing& x = cr.crossings[k];
        bool newerA = cs.yarns[x.yarnA].course == ps.olderCourse + 1;
        const YarnPath& yn = cs.yarns[newerA ? x.yarnA : x.yarnB];
        const YarnPath& yo = cs.yarns[newerA ? x.yarnB : x.yarnA];
        double un = newerA ? x.uA : x.uB, uo = newerA ? x.uB : x.uA;
        double dn = newerA ? x.depth.mid() : -x.depth.mid();
        v.push_back({localParam(yn, newerA ? x.stitchA : x.stitchB, un), localParam(yo, newerA ? x.stitchB : x.stitchA, uo), dn});
    }
    std::sort(v.begin(), v.end(), [](const PairCrossing& a, const PairCrossing& b) { return a.ln < b.ln; });
    return v;
}

// Interior signatures: every interior pair reproduces the reference pair of its class: the faces of
// the 3 x 3 stitches around the child (one class for stockinette).
void checkInteriorSignatures(const CurveSet& cs, const PatternParams& pp, const ValidationReport& r, GeneratorChecks& out) {
    const CrossingResult& cr = r.crossings;
    Status st = Status::Pass;
    std::vector<std::pair<int, const PairSignature*>> refs;
    auto pairClass = [&](const PairSignature& ps) {
        int key = 0;
        for (int dr = -1; dr <= 1; ++dr)
            for (int dw = -1; dw <= 1; ++dw)
                key = 2 * key + stitchFaceIsPurl(pp.pattern, ps.wale + dw, ps.olderCourse + 1 + dr);
        return key;
    };
    auto interior = [&](const PairSignature& ps) {
        return ps.wale > 0 && ps.wale + 1 < cs.wales && ps.olderCourse > 0 && ps.olderCourse + 2 < cs.courses;
    };
    int compared = 0;
    double maxDu = 0, maxDz = 0;
    for (const PairSignature& ps : r.signature.pairs) {
        if (!interior(ps) || ps.status != Status::Pass) continue;
        const int cls = pairClass(ps);
        const PairSignature* ref = nullptr;
        for (const auto& [c, candidate] : refs)
            if (c == cls) ref = candidate;
        if (!ref) { refs.push_back({cls, &ps}); continue; }
        std::vector<PairCrossing> a = pairCrossings(cs, cr, *ref), b = pairCrossings(cs, cr, ps);
        if (a.size() != b.size()) { st = Status::Fail; continue; }
        for (size_t q = 0; q < a.size(); ++q) {
            maxDu = std::max({maxDu, std::fabs(a[q].ln - b[q].ln), std::fabs(a[q].lo - b[q].lo)});
            maxDz = std::max(maxDz, std::fabs(a[q].depth - b[q].depth));
        }
        ++compared;
    }
    if (maxDu > 1e-7 || maxDz > 1e-7 * r.allowances.Lref) st = Status::Fail;
    if (refs.empty()) st = Status::Unresolved;  // patch too small to have interior pairs
    out.interiorSignatures = st;
    char buf[256];
    std::snprintf(buf, sizeof buf, "controls max dev %.2e mm; %d interior pairs vs reference: max du %.1e, max dz %.1e mm",
                  out.maxPeriodicityDeviation, compared, maxDu, maxDz);
    out.detail = buf;
}

GeneratorChecks generatorChecks(const CurveSet& cs, const PatternParams& pp, const ValidationReport& r) {
    GeneratorChecks out;
    checkPeriodicity(cs, pp, r.allowances, out);
    checkInteriorSignatures(cs, pp, r, out);
    return out;
}

// --- the checks in the order validate runs them

void runCurvature(ValidationReport& r, const CurveSet& cs, const PatternParams& pp, const ExecPolicy& exec) {
    auto t0 = Clock::now();
    r.curvature.clear();
    for (const YarnPath& y : cs.yarns) r.curvature.push_back(checkCurvature(y, pp.radius(), r.allowances, exec));
    r.times.curvatureMs = msSince(t0);
}

void runSampling(ValidationReport& r, const CurveSet& cs, const PatternParams& pp, const ExecPolicy& exec) {
    auto t0 = Clock::now();
    r.polylines.clear();
    for (const YarnPath& y : cs.yarns) r.polylines.push_back(samplePath(y, pp.sampleTolerance, 512, exec));
    r.times.samplingMs = msSince(t0);
}

void runCrossings(ValidationReport& r, const CurveSet& cs, const PatternParams& pp) {
    auto t0 = Clock::now();
    r.crossings = findCrossings(cs, r.allowances);
    r.signature = checkSignature(cs, r.crossings, pp);
    r.times.crossingsMs = msSince(t0);
}

void runClearance(ValidationReport& r, const CurveSet& cs, const PatternParams& pp, const ExecPolicy& exec) {
    auto t0 = Clock::now();
    r.clearance = checkClearance({cs, r.polylines, r.curvature}, pp, r.allowances, exec);
    r.times.clearanceMs = msSince(t0);
}

} // namespace

void summarize(ValidationReport& r, const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt) {
    r.lines.clear();
    Status reg = Status::Pass, fold = Status::Pass;
    double speedLo = 1e300, rhoLo = 1e300, rhoEst = 1e300;
    for (const CurvatureResult& c : r.curvature) {
        reg = combine(reg, c.regular);
        fold = combine(fold, c.noFold);
        speedLo = std::min(speedLo, c.speedLo);
        rhoLo = std::min(rhoLo, c.rhoLo);
        rhoEst = std::min(rhoEst, c.rhoEst);
    }
    char buf[512];
    std::snprintf(buf, sizeof buf, "min speed lower bound %.4g mm/param", speedLo);
    r.lines.push_back({"regular (|C'| > 0)", reg, buf});
    std::snprintf(buf, sizeof buf, "rho_lo %.4f mm (bound), rho_est %.4f mm (sampled), r = %.4f mm", rhoLo, rhoEst,
                  pp.radius());
    r.lines.push_back({"no local fold (rho >= r)", fold, buf});
    if (cs.kind == CurveKind::FlatPatch) {
        r.lines.push_back({"crossing signature", r.signature.status, r.signature.detail});
        std::snprintf(buf, sizeof buf, "min |depth gap| at crossings %.4f mm vs D+c = %.4f mm",
                      r.signature.minDepthGapLo, pp.yarnDiameter + pp.clearance);
        r.lines.push_back({"depth gap at crossings (necessary)", r.signature.depthGapVsClearance, buf});
    } else {
        r.lines.push_back({"local signature", r.localSignature.status, r.localSignature.detail});
    }
    int fails = 0, unres = 0, pen = 0;
    for (const ClearanceRegion& g : r.clearance.regions) {
        fails += g.status == Status::Fail;
        unres += g.status == Status::Unresolved;
        pen += g.penetrating;
    }
    std::snprintf(buf, sizeof buf,
                  "min inter-yarn %.4f mm, min self nonlocal %.4f mm (-1: none beyond pi*rho), threshold D+c = %.4f mm, "
                  "max allowance +-%.2e mm; %d FAIL regions (%d penetrating), %d UNRESOLVED",
                  r.clearance.minInterYarn > 1e299 ? -1.0 : r.clearance.minInterYarn,
                  r.clearance.minSelfNonlocal > 1e299 ? -1.0 : r.clearance.minSelfNonlocal, pp.yarnDiameter + pp.clearance,
                  r.clearance.deltaD, fails, pen, unres);
    r.lines.push_back({"clearance", r.clearance.status, buf});
    r.generatorChecksApplied = false;
    if (opt.generatorChecks && cs.kind == CurveKind::FlatPatch) {
        r.generator = generatorChecks(cs, pp, r);
        r.generatorChecksApplied = true;
        r.lines.push_back({"generator: control periodicity", r.generator.periodicity, r.generator.detail});
        r.lines.push_back({"generator: interior signatures identical", r.generator.interiorSignatures, ""});
    }
    r.verdict = Status::Pass;
    for (const CheckLine& l : r.lines) r.verdict = combine(r.verdict, l.status);
}

ValidateStats& validateStats() {
    thread_local ValidateStats s;
    return s;
}

ValidationReport validate(const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt) {
    ValidateStats& st = validateStats();
    const auto t0 = Clock::now();
    ValidationReport r;
    r.allowances = Allowances::forExtent(cs.extent());
    runCurvature(r, cs, pp, opt.exec);
    const auto t1 = Clock::now();
    runSampling(r, cs, pp, opt.exec);
    const auto t2 = Clock::now();
    // A global front projection only means something for the flat patch; on a
    // surface the signature is checked per stitch in its own frame.
    if (cs.kind == CurveKind::FlatPatch) runCrossings(r, cs, pp);
    else r.localSignature = checkLocalSignature(cs, pp, opt.exec, opt.attributionWindow);
    const auto t3 = Clock::now();
    runClearance(r, cs, pp, opt.exec);
    const auto t4 = Clock::now();
    summarize(r, cs, pp, opt);
    auto d = [](Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
    ++st.calls;
    st.msCurvature += d(t0, t1);
    st.msSampling += d(t1, t2);
    st.msSignature += d(t2, t3);
    st.msClearance += d(t3, t4);
    st.ms += msSince(t0);
    return r;
}

std::string formatReport(const ValidationReport& r, const CurveSet& cs, const PatternParams& pp, bool verbose) {
    std::ostringstream os;
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "patch %dx%d  wale spacing %.3f mm (%.4f wales/mm, %.2f per 10cm)  course spacing %.3f mm (%.4f "
                  "courses/mm, %.2f per 10cm)\n  D %.3f mm  c %.3f mm  A %.3f D  eps %.4f mm  head %.3f  leg %.3f\n",
                  pp.wales, pp.courses, pp.waleSpacing, 1.0 / pp.waleSpacing, 100.0 / pp.waleSpacing, pp.courseSpacing,
                  1.0 / pp.courseSpacing, 100.0 / pp.courseSpacing, pp.yarnDiameter, pp.clearance, pp.depthAmplitude,
                  pp.sampleTolerance, pp.headHalfWidth, pp.legHalfWidth);
    os << buf;
    const char* verdict = r.verdict == Status::Pass ? "VALID (verified within tolerance)"
                          : r.verdict == Status::Fail ? "INVALID" : "UNRESOLVED";
    os << "  verdict: " << verdict << "\n";
    for (const CheckLine& l : r.lines) {
        std::snprintf(buf, sizeof buf, "  [%-10s] %-40s %s\n", toString(l.status), l.name.c_str(), l.detail.c_str());
        os << buf;
    }
    // Loop (stitch) length of an interior stitch, and yarn lengths.
    if (!cs.yarns.empty() && !cs.yarns[cs.yarns.size() / 2].stitches.empty()) {  // a yarn without stitch records has no interior stitch
        const YarnPath& y = cs.yarns[cs.yarns.size() / 2];
        const StitchSpan& st = y.stitches[y.stitches.size() / 2];
        const ArcLengthBounds loop = y.lengthBounds(st.uBegin, st.uEnd);
        double est = y.arcLengthAt(st.uEnd) - y.arcLengthAt(st.uBegin);
        std::snprintf(buf, sizeof buf, "  loop length (interior stitch) %.4f mm  [enclosure %.4f .. %.4f]; yarn %d length %.4f mm (est. err %.1e)\n",
                      est, loop.lo, loop.hi, y.id, y.length(), y.lengthErrorEstimate());
        os << buf;
    }
    std::snprintf(buf, sizeof buf,
                  "  timings: curvature %.1f ms, sampling %.1f ms, crossings %.1f ms (%ld piece pairs), clearance %.1f ms (%ld segment pairs)\n",
                  r.times.curvatureMs, r.times.samplingMs, r.times.crossingsMs, r.crossings.pairsVisited,
                  r.times.clearanceMs, r.clearance.segmentPairs);
    os << buf;
    if (verbose) {
        for (const UnresolvedItem& u : r.crossings.unresolved) {
            std::snprintf(buf, sizeof buf, "    unresolved crossing search: yarn %d u[%.4f,%.4f] vs yarn %d u[%.4f,%.4f]: %s\n",
                          u.yarnA, u.uA0, u.uA1, u.yarnB, u.uB0, u.uB1, u.reason.c_str());
            os << buf;
        }
        for (int k : r.signature.unexpected) {
            const Crossing& x = r.crossings.crossings[k];
            std::snprintf(buf, sizeof buf, "    unexpected crossing: yarn %d u=%.4f / yarn %d u=%.4f at (%.3f, %.3f) %s\n",
                          x.yarnA, x.uA, x.yarnB, x.uB, x.pA.x, x.pA.y, x.note.c_str());
            os << buf;
        }
        for (const PairSignature& ps : r.signature.pairs) {
            if (ps.status == Status::Pass) continue;
            std::snprintf(buf, sizeof buf, "    pair wale %d courses %d/%d: %s %s\n", ps.wale, ps.olderCourse,
                          ps.olderCourse + 1, toString(ps.status), ps.detail.c_str());
            os << buf;
        }
        int shown = 0;
        for (const ClearanceRegion& g : r.clearance.regions) {
            if (++shown > 30) { os << "    ...\n"; break; }
            std::snprintf(buf, sizeof buf, "    clearance %s: yarn %d u[%.3f,%.3f] vs yarn %d u[%.3f,%.3f] worst d %.4f (thr %.4f)%s [%s]\n",
                          toString(g.status), g.yarnA, g.uA0, g.uA1, g.yarnB, g.uB0, g.uB1, g.worstDistance,
                          g.threshold, g.penetrating ? " PENETRATING" : "",
                          g.worstIssue >= 0 ? r.clearance.issues[g.worstIssue].reason.c_str() : "");
            os << buf;
        }
    }
    return os.str();
}

MotionCertificate certifyStraightMotion(const CurveSet& a, const ValidationReport& ra, const CurveSet& b,
                                        const ValidationReport& rb, const PatternParams& pp) {
    MotionCertificate m;
    auto ok = [](const ValidationReport& r) {
        Status fold = Status::Pass;
        for (const CheckLine& l : r.lines)
            if (l.name.rfind("no local fold", 0) == 0 || l.name.rfind("regular", 0) == 0) fold = combine(fold, l.status);
        return r.clearance.status == Status::Pass && fold == Status::Pass;
    };
    m.endStatesValid = ok(ra) && ok(rb);
    if (a.yarns.size() != b.yarns.size()) { m.detail = "different curve structure"; return m; }
    for (size_t y = 0; y < a.yarns.size(); ++y) {
        const auto &P = a.yarns[y].curve.ctrl, &Q = b.yarns[y].curve.ctrl;
        if (P.size() != Q.size()) { m.detail = "different curve structure"; return m; }
        for (size_t k = 0; k < P.size(); ++k) m.maxDisplacement = std::max(m.maxDisplacement, norm(Q[k] - P[k]));
    }
    // Rounding allowance on the displacement (positions are exact to etaX each).
    const double delta = m.maxDisplacement + 2 * Allowances::forExtent(std::max(a.extent(), b.extent())).etaX;
    const double g = pp.yarnDiameter + pp.clearance;  // certified nonlocal distance at both ends
    m.noPassThrough = m.endStatesValid && delta < g;
    m.noTubeContact = m.endStatesValid && delta <= pp.clearance;
    m.localScale = 0.5 * kPi * delta;
    char buf[200];
    std::snprintf(buf, sizeof buf, "max control displacement %.4f mm vs D + c = %.3f mm (pass-through) and c = %.3f mm (tube contact)%s",
                  m.maxDisplacement, g, pp.clearance, m.endStatesValid ? "" : "; an end state is not valid");
    m.detail = buf;
    return m;
}

namespace {

// Categories (FailCategory, model/params.h) with a FAIL line, and with an UNRESOLVED one.
unsigned failCategories(const ValidationReport& r, unsigned* unresolved) {
    unsigned f = 0, u = 0;
    auto add = [&](Status s, unsigned cat) {
        if (s == Status::Fail) f |= cat;
        else if (s == Status::Unresolved) u |= cat;
    };
    for (const CurvatureResult& c : r.curvature) {
        add(c.regular, kFailCurvature);
        add(c.noFold, kFailCurvature);
    }
    if (!r.signature.pairs.empty() || !r.signature.unexpected.empty()) {
        Status pairs = Status::Pass;
        for (const PairSignature& ps : r.signature.pairs) pairs = combine(pairs, ps.status);
        add(pairs, kFailSignature);
        add(r.signature.overlapWrongOrder || r.signature.layerWrongOrder ? Status::Fail : Status::Pass, kFailSignature);
        if (!r.signature.unexpected.empty()) add(Status::Fail, kFailCrossing);
        add(r.signature.depthGapVsClearance, kFailCrossing);
    }
    add(r.localSignature.status, kFailSignature);
    add(r.clearance.status, kFailClearance);
    if (r.generatorChecksApplied) {
        add(r.generator.periodicity, kFailGenerator);
        add(r.generator.interiorSignatures, kFailGenerator);
    }
    if (unresolved) *unresolved = u;
    return f;
}

} // namespace

PatternReference flatPatternReference(const PatternParams& base, StitchPattern p, int n) {
    PatternParams pp = base;
    pp.pattern = p;
    pp.wales = pp.courses = n;
    const CurveSet cs = generateStockinette(pp);
    const ValidationReport r = validate(cs, pp);
    PatternReference out;
    out.pattern = p;
    out.verdict = r.verdict;
    out.fail = failCategories(r, &out.unresolved);
    out.minInterYarn = r.clearance.minInterYarn;
    char b[160];
    std::snprintf(b, sizeof b, "D %.2f, %.2f x %.2f mm, A %.2f, c %.2f", pp.yarnDiameter, pp.waleSpacing, pp.courseSpacing,
                  pp.depthAmplitude, pp.clearance);
    out.params = b;
    return out;
}

namespace {

// One formatted piece of text, cut at 511 characters: stored digests hold that cut.
void appendFormatted(std::string& o, const char* format, ...) {
    char b[512];
    va_list values;
    va_start(values, format);
    std::vsnprintf(b, sizeof b, format, values);
    va_end(values);
    o += b;
}

} // namespace

std::string validationSnapshot(const ValidationReport& r) {
    std::string o;
    if (r.gpuDiagnostic) appendFormatted(o, "producer gpu float32 diagnostic\n");  // absent for the CPU (digests unchanged)
    appendFormatted(o, "verdict %s generator %d\n", toString(r.verdict), int(r.generatorChecksApplied));
    for (const CheckLine& l : r.lines) appendFormatted(o, "line %s %s %s\n", l.name.c_str(), toString(l.status), l.detail.c_str());
    for (const CurvatureResult& c : r.curvature) {
        appendFormatted(o, "curv %d %a %a %a %a %a %a %a %s %s %d\n", c.yarnId, c.speedLo, c.kappaHi, c.rhoLo, c.rhoEst, c.uRhoEst,
            c.kappaWitnessLo, c.uWitness, toString(c.regular), toString(c.noFold), c.leaves);
        for (double k : c.spanKappaHi) appendFormatted(o, " %a", k);
        o += "\n";
    }
    for (const Polyline& p : r.polylines) appendFormatted(o, "poly %zu %a\n", p.p.size(), p.epsBound);
    const ClearanceResult& c = r.clearance;
    appendFormatted(o, "clr %s %a %a %a %ld %ld %zu %zu\n", toString(c.status), c.minInterYarn, c.minSelfNonlocal, c.deltaD, c.segmentPairs,
        c.schurCertified, c.issues.size(), c.regions.size());
    for (const ClearanceIssue& is : c.issues)
        appendFormatted(o, "is %d %d %a %a %a %a %a %a %a %a %a %a %s %d %s\n", is.yarnA, is.yarnB, is.uA, is.uB, is.pA.x, is.pA.y, is.pA.z, is.pB.x,
            is.pB.y, is.pB.z, is.distance, is.threshold, toString(is.status), int(is.penetrating), is.reason.c_str());
    for (const ClearanceRegion& g : c.regions)
        appendFormatted(o, "rg %d %d %a %a %a %a %a %a %s %d %d\n", g.yarnA, g.yarnB, g.uA0, g.uA1, g.uB0, g.uB1, g.worstDistance, g.threshold,
            toString(g.status), int(g.penetrating), g.worstIssue);
    const LocalSignatureResult& s = r.localSignature;
    appendFormatted(o, "sig %s %d %d %d %d %s\n", toString(s.status), s.pass, s.fail, s.unresolved, s.neighbourhoodSelf, s.detail.c_str());
    for (const LocalPairResult& lp : s.pairs)
        appendFormatted(o, "lp %d %d %s %d %d %d %d %d %s\n", lp.child, lp.parent, toString(lp.status), lp.crossings, lp.unexpected, lp.overlap,
            lp.terminal, lp.neighbourhoodSelf, lp.detail.c_str());
    appendFormatted(o, "crossings %zu %zu sig %s %zu\n", r.crossings.crossings.size(), r.crossings.unresolved.size(), toString(r.signature.status),
        r.signature.pairs.size());
    return o;
}

} // namespace knit
