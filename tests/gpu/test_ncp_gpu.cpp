// ncp-gpu: the device problem against the independent CPU evaluator
// (knit_problem.h) on compact fixtures built from production types (B-spline yarns, the cylinder
// spiral's knitting, a bound planar mesh), and the solver on fixtures with known answers:
//   rows      every CPU row (value and gradient per variable) has its device row; extra device
//             contact rows only beyond the coverage pad (the float boxes' guard);
//   strand    a jittered strand between fixed ends returns to its rest line (objective -> 0);
//   pair      two parallel strands closer than the hard contact are pushed apart to it;
//   infeasible crossing strands whose only free control is capped: UNRESOLVED, never CONVERGED;
//   swept     a trial that moves one strand through another is caught by the swept test even
//             though both end states are separated;
//   coverage  a strongly jittered patch: after the solve no pair (brute force over all segment
//             pairs) is closer than the hard contact minus the allowance;
//   scaling   cylinder spirals of increasing size (total problem, not a free mask): time, steps.
// SKIP without a device.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gpu/device/runtime.h"
#include "gpu/port/ncp_call_file.h"
#include "gpu/solver/ncp.h"
#include "gpu/surface/placement.h"
#include "gpu_common.h"
#include "../cpu_numerics.h"
#include "knitting/checks/checks.h"
#include "../ncp_reference.h"
#include "knitting/layout/stitch_layout.h"
#include "model/mesh.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/knit/fixtures.h"
#include "knitting/relax/settle.h"
#include "knitting/stitch/stitch_program.h"
#include "knitting/stitch/stockinette.h"
#include "knitting/layout/uv_atlas.h"
#include "knitting/knit/mesh_knit.h"
#include "gpu/solver/ncp.h"
#include "gpu/device.h"
#include "gpu/knitting.h"
#include "model/cancel.h"
#include "model/env.h"

using namespace knit;

namespace {

// A device for the row tests: no solve follows.
gpu::NcpCall rowsCall(const CurveSet& cs, Controls& P, const Controls& rest, const ControlMask& freeCtrl, const BaseSurface* surface, SurfaceBinding* bind,
                      const RelaxMoveParams& mp, const KnitProblemParams& kp) {
    static const NcpSolveParams solve;
    static const gpu::NcpOptions options;
    return {cs, P, rest, freeCtrl, surface, bind, mp, kp, solve, options};
}

} // namespace

namespace {

#define REQUIRE_NCP()                                                     \
    do {                                                                  \
        REQUIRE_GPU();                                                    \
        std::string why_;                                                 \
        if (!gpu::ncpAvailable(&why_)) SKIP("the ncp solver cannot run"); \
    } while (0)

YarnPath makeYarn(std::vector<Vec3> ctrl, int id) {
    YarnPath y;
    y.id = id;
    y.course = id;
    y.name = "test " + std::to_string(id);
    y.curve.ctrl = std::move(ctrl);
    y.stitches.push_back({0, id, 0.0, y.curve.uMax()});
    y.finalize();
    return y;
}
CurveSet makeSet(std::vector<YarnPath> ys) {
    CurveSet cs;
    cs.yarns = std::move(ys);
    cs.wales = 1;
    cs.courses = static_cast<int>(cs.yarns.size());
    return cs;
}
std::vector<Vec3> lineCtrl(Vec3 a, Vec3 b, int n) {
    std::vector<Vec3> c;
    const Vec3 step = (b - a) / n;
    for (int i = -1; i <= n + 1; ++i) c.push_back(a + step * i);
    return c;
}

struct CylinderPiece {
    PatternParams pp;
    GaugeParams g;
    StitchProgram prog;
    PieceLayout lay;
    CurveSet cs;
    std::unique_ptr<ConeSurface> surf;
};
CylinderPiece cylinderPiece(int n, int turns) {  // tests.cpp cylinderPiece
    CylinderPiece c;
    c.pp = findPreset("tightest")->params;
    c.pp.waleSpacing = c.g.w;
    c.pp.courseSpacing = c.g.h;
    const double R = n * c.g.w / (2 * kPi), H = (turns + 2) * c.g.h;
    c.surf = std::make_unique<ConeSurface>(R, R, H);
    c.prog = buildProgram(std::vector<int>(turns, n));
    c.lay = helicalLayout(*c.surf, c.prog, c.g, evenHelix(c.prog, c.g.h, H));
    c.cs = placeSpiral(*c.surf, c.lay, c.prog, c.pp);
    finalizeCurves(c.cs);
    return c;
}

// The cylinder spiral with stitches [first, first + n) free, jittered (tests.cpp ncpFixture).
struct Fixture {
    CylinderPiece c;
    ControlMask freeCtrl;
    RelaxMoveParams mp;
    Controls P;
};
Fixture cylinderFixture(int n, int turns, double D, double cap, int first, int stitches, double jitter) {
    Fixture f;
    f.c = cylinderPiece(n, turns);
    PatternParams pp = f.c.pp;
    pp.yarnDiameter = D;
    const double Dc = pp.yarnDiameter + pp.clearance;
    f.mp.thr = Dc + 0.05;
    f.mp.kappaMax = 1.0 / (1.15 * (pp.radius() + 0.5 * pp.clearance));
    f.mp.window = 0.5 * kPi * Dc * 1.3;
    f.mp.maxDisplacement = cap;
    const YarnPath& y = f.c.cs.yarns[0];
    f.freeCtrl.assign(1, std::vector<char>(y.curve.ctrl.size(), 0));
    for (int t = first; t < first + stitches && t < static_cast<int>(f.c.cs.topo.stitches.size()); ++t) {
        const StitchRecord& st = f.c.cs.topo.stitches[t];
        for (int c = static_cast<int>(st.uBegin) + 1; c <= static_cast<int>(st.uEnd) && c + 1 < static_cast<int>(y.curve.ctrl.size()); ++c)
            f.freeCtrl[0][c] = 1;
    }
    f.P = {y.curve.ctrl};
    for (size_t c = 0; c < f.P[0].size(); ++c)
        if (f.freeCtrl[0][c]) {
            const double a = 0.7 * c, b = 1.3 * c;
            f.P[0][c] += Vec3{std::sin(a), std::cos(b), std::sin(a + b)} * jitter;
        }
    return f;
}

// Device rows against KnitProblem::linearize at the same controls.
struct RowCompare {
    long long cpuRows = 0, matched = 0, missing = 0, extra = 0, extraInsidePad = 0, gradientsCompared = 0, contactGradientsCompared = 0;
    double maxValueDiff[kRowFamilies] = {}, maxGradDiff[kRowFamilies] = {};
};
RowCompare compareRows(const KnitProblem& pr, const Controls& P, const std::vector<gpu::NcpDevice::DebugRow>& dev, double padMinusGuard) {
    RowCompare rc;
    std::vector<Row> cpu;
    pr.linearize(P, cpu);
    std::map<uint64_t, const gpu::NcpDevice::DebugRow*> byKey;
    for (const auto& d : dev) byKey[d.key] = &d;
    std::map<uint64_t, bool> seen;
    for (const Row& r : cpu) {
        ++rc.cpuRows;
        auto it = byKey.find(r.key);
        if (it == byKey.end()) {
            ++rc.missing;
            continue;
        }
        ++rc.matched;
        seen[r.key] = true;
        const auto& d = *it->second;
        const int f = static_cast<int>(r.family);
        rc.maxValueDiff[f] = std::max(rc.maxValueDiff[f], std::fabs(d.value - r.value));
        // The device evaluates a contact row's gradient only where the row has a term in the merit.
        if (!d.hasGradient) continue;
        ++rc.gradientsCompared;
        if (r.family == RowFamily::Contact) ++rc.contactGradientsCompared;
        for (int e = 0; e < r.n; ++e) {
            int k = -1;
            for (int j = 0; j < d.n; ++j)
                if (d.var[j] == r.e[e].var) k = j;
            const double gx = k >= 0 ? d.g[k][0] : 0, gy = k >= 0 ? d.g[k][1] : 0, gz = k >= 0 ? d.g[k][2] : 0;
            const double diff = std::max({std::fabs(gx - r.e[e].g.x), std::fabs(gy - r.e[e].g.y), std::fabs(gz - r.e[e].g.z)});
            rc.maxGradDiff[f] = std::max(rc.maxGradDiff[f], diff);
        }
    }
    for (const auto& d : dev)
        if (!seen.count(d.key)) {
            ++rc.extra;
            // a device contact row the CPU did not emit must lie beyond the CPU's pad (thr + 2 delta)
            if (d.family == static_cast<int>(RowFamily::Contact)) {
                const double dist = pr.params().hardContact > 0 ? pr.params().hardContact - d.value * pr.D() : pr.moveParams().thr - d.value * pr.D();
                if (dist < padMinusGuard) ++rc.extraInsidePad;
            } else {
                ++rc.extraInsidePad;  // only contact rows may differ
            }
        }
    return rc;
}

void reportRows(const char* what, const RowCompare& rc) {
    std::printf("    %s: %lld CPU rows, %lld matched (%lld with gradients), %lld missing, %lld extra (%lld inside the pad)\n", what, rc.cpuRows, rc.matched,
                rc.gradientsCompared, rc.missing, rc.extra, rc.extraInsidePad);
    for (int f = 0; f < kRowFamilies; ++f)
        if (rc.maxValueDiff[f] > 0 || rc.maxGradDiff[f] > 0)
            std::printf("      %-17s max |value diff| %.3g, max |gradient diff| %.3g\n", toString(static_cast<RowFamily>(f)), rc.maxValueDiff[f],
                        rc.maxGradDiff[f]);
}

// Tolerances: double rows on both sides, the band included (the device evaluates the surface in
// double on the CPU's mesh: closest point, continuation with the CPU's constants).
void checkRows(const RowCompare& rc) {
    CHECK(rc.missing == 0);
    CHECK(rc.extraInsidePad == 0);
    for (int f = 0; f < kRowFamilies; ++f) {
        CHECK(rc.maxValueDiff[f] <= 1e-9);
        CHECK(rc.maxGradDiff[f] <= 1e-9);
    }
}

} // namespace

TEST(ncp_gpu_rows_match_cpu_cylinder) {
    REQUIRE_NCP();
    // The cylinder band (analytic), the cap, curvature and contact targets (hard thresholds), anchors.
    // A thick yarn on the cylinder's gauge: contact rows closer than thr, with a term in the merit.
    Fixture f = cylinderFixture(24, 5, 1.3, 0.5, 30, 6, 0.05);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.04;
    kp.hardKappa = f.mp.kappaMax * 1.1;
    kp.wA = 0.01;
    KnitProblem pr(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp);
    pr.refresh(f.P, nullptr);
    std::string err;
    gpu::NcpDevice dev(rowsCall(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp));
    dev.refresh();
    gpu::NcpEvalSummary s;
    s = dev.linearize(gpu::ControlSet::Current, 10);
    std::vector<gpu::NcpDevice::DebugRow> rows;
    rows = dev.debugRows();
    const double pad = f.mp.thr + 2 * kp.delta;
    RowCompare rc = compareRows(pr, f.P, rows, pad - 1e-6);
    reportRows("at the refresh", rc);
    checkRows(rc);
    std::printf("      contact rows with gradients compared: %lld\n", rc.contactGradientsCompared);
    CHECK(rc.contactGradientsCompared > 0);  // contact rows with a term in the merit
    CHECK(rc.matched > 1000);
    // Moved within the coverage (no refresh): the band and every row follow the current controls.
    Controls Q = f.P;
    for (int v = 0; v < pr.vars(); ++v) {
        const auto [y, c] = pr.ctrlOf(v);
        Q[y][c] += Vec3{0.03 * std::cos(v), 0.02 * std::sin(2.0 * v), 0.025};
    }
    dev.setTrial(Q);
    s = dev.linearize(gpu::ControlSet::Trial, 10);
    dev.accept();
    rows = dev.debugRows();
    rc = compareRows(pr, Q, rows, pad - 2 * 0.05 - 1e-6);  // an extra row (beyond the pad at P) may have come up to 2 x 0.05 closer
    reportRows("moved by <= 0.05 mm", rc);
    checkRows(rc);
}

TEST(ncp_gpu_rows_match_cpu_bound_mesh) {
    REQUIRE_NCP();
    // A strand over a bound planar mesh (MeshSurface continuation on the device, float).
    TriMesh m;
    const int N = 10;
    for (int j = 0; j <= N; ++j)
        for (int i = 0; i <= N; ++i) m.pos.push_back({2.0 * i, 2.0 * j, 0.03 * std::sin(0.9 * i + 0.4 * j)});
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            const int a = j * (N + 1) + i;
            m.tri.push_back({a, a + 1, a + N + 2});
            m.tri.push_back({a, a + N + 2, a + N + 1});
        }
    m.computeNormals();
    MeshSurface ms(m);
    std::vector<Vec3> ctrl = lineCtrl({1, 10, 0.25}, {19, 10, 0.25}, 9);
    for (size_t c = 0; c < ctrl.size(); ++c) ctrl[c].z += 0.05 * std::sin(1.7 * c);
    CurveSet cs = makeSet({makeYarn(ctrl, 0)});
    ControlMask freeCtrl(1, std::vector<char>(ctrl.size(), 0));
    for (size_t c = 1; c + 1 < ctrl.size(); ++c) freeCtrl[0][c] = 1;
    SurfaceBinding bind;
    bind.surface = &ms;
    bind.state.assign(1, std::vector<ProjState>(ctrl.size()));
    bind.result.assign(1, std::vector<ProjResult>(ctrl.size(), ProjResult::Owned));
    for (size_t c = 0; c < ctrl.size(); ++c) {
        const Vec3 q = ctrl[c];
        const int i = std::clamp(static_cast<int>(q.x / 2), 0, N - 1), j = std::clamp(static_cast<int>(q.y / 2), 0, N - 1);
        const double b[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
        bind.state[0][c] = ms.stateAt(2 * (j * N + i) + 1, b, q);
    }
    RelaxMoveParams mp;
    mp.thr = 0.5;
    mp.kappaMax = 10;
    mp.window = 100;
    mp.bandHalf = 0.2;
    Controls P{ctrl};
    SurfaceBinding bc = bind, bg = bind;
    KnitProblem pr(cs, P, P, freeCtrl, &ms, &bc, mp);
    pr.refresh(P, &bc);
    std::string err;
    gpu::NcpDevice dev(rowsCall(cs, P, P, freeCtrl, &ms, &bg, mp, {}));
    dev.refresh();
    gpu::NcpEvalSummary s;
    s = dev.linearize(gpu::ControlSet::Current, 10);
    std::vector<gpu::NcpDevice::DebugRow> rows;
    rows = dev.debugRows();
    RowCompare rc = compareRows(pr, P, rows, mp.thr + 0.2 - 1e-6);
    reportRows("bound mesh", rc);
    checkRows(rc);
    CHECK(s.notOwned == 0);
}

TEST(ncp_gpu_ownership_loss_pins_the_control) {
    REQUIRE_NCP();
    // A strand on a bound planar mesh running towards its open boundary (x = 20), rest lengths 1.4x
    // the current ones: the stretched edges push the free end past the boundary, where the
    // continuation is not Owned (OpenBoundary). A trial that loses a binding is rejected; the
    // controls that lose it are pinned and the step is retried at the same damping, so the rest of
    // the strand still relaxes (before the pins every such step was shortened towards the boundary
    // with growing damping until the solver stalled; sphere, GEOMETRY 8.21). After the solve every
    // binding is still Owned, at least one control was pinned and the objective decreased.
    TriMesh m;
    const int N = 10;
    for (int j = 0; j <= N; ++j)
        for (int i = 0; i <= N; ++i) m.pos.push_back({2.0 * i, 2.0 * j, 0});
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            const int a = j * (N + 1) + i;
            m.tri.push_back({a, a + 1, a + N + 2});
            m.tri.push_back({a, a + N + 2, a + N + 1});
        }
    m.computeNormals();
    MeshSurface ms(m);
    // every control (lineCtrl's end controls at x = 7, 19) starts inside the mesh, Owned
    const std::vector<Vec3> ctrl = lineCtrl({8, 10, 0.05}, {18, 10, 0.05}, 10);
    std::vector<Vec3> rest = ctrl;
    for (size_t c = 0; c < rest.size(); ++c) rest[c] = ctrl[0] + (ctrl[c] - ctrl[0]) * 1.4;
    CurveSet cs = makeSet({makeYarn(ctrl, 0)});
    ControlMask freeCtrl(1, std::vector<char>(ctrl.size(), 0));
    for (size_t c = 3; c < ctrl.size(); ++c) freeCtrl[0][c] = 1;
    SurfaceBinding bind;
    bind.surface = &ms;
    bind.state.assign(1, std::vector<ProjState>(ctrl.size()));
    bind.result.assign(1, std::vector<ProjResult>(ctrl.size(), ProjResult::Owned));
    for (size_t c = 0; c < ctrl.size(); ++c) {
        const Vec3 q = ctrl[c];
        const int i = std::clamp(static_cast<int>(q.x / 2), 0, N - 1), j = std::clamp(static_cast<int>(q.y / 2), 0, N - 1);
        const double b[3] = {1.0 / 3, 1.0 / 3, 1.0 / 3};
        bind.state[0][c] = ms.stateAt(2 * (j * N + i) + 1, b, q);
    }
    RelaxMoveParams mp;
    mp.thr = 0.5;
    mp.kappaMax = 10;
    mp.window = 100;
    mp.bandHalf = 0.2;
    Controls P{ctrl}, Prest{rest};
    auto objective = [&](const Controls& X, const SurfaceBinding& b) {
        SurfaceBinding bc = b;
        KnitProblem pr(cs, X, Prest, freeCtrl, &ms, &bc, mp);
        pr.refresh(X, &bc);
        return pr.scoreGeometry(X, &bc, 1e-3).objective;
    };
    const double obj0 = objective(P, bind);
    long long notOwned0 = 0;
    {
        SurfaceBinding bc = bind;
        KnitProblem pr(cs, P, Prest, freeCtrl, &ms, &bc, mp);
        pr.refresh(P, &bc);
        for (const ProjResult r : bc.result[0]) notOwned0 += r != ProjResult::Owned;
    }
    CHECK(notOwned0 == 0);
    SurfaceBinding bg = bind;
    Controls Q = P;
    NcpSolveReport rep;
    gpu::NcpStats st;
    std::string err;
    NcpSolveParams sp;
    sp.maxOuter = 10;
    CHECK(gpu::solveNcp({cs, Q, Prest, freeCtrl, &ms, &bg, mp, {}, sp, {}}, rep, &st, &err));
    long long notOwned = 0;
    for (const ProjResult r : bg.result[0]) notOwned += r != ProjResult::Owned;
    const double obj1 = objective(Q, bg);
    std::printf("    %s; %lld LM steps, %lld ownership rejections, %lld controls pinned; objective %.6g -> %.6g; end x %.3f; not Owned %lld\n",
                rep.status.c_str(), st.lmSteps, st.ownershipRejected, st.ownershipPinned, obj0, obj1, Q[0].back().x, notOwned);
    CHECK(err.empty());
    CHECK(st.ownershipPinned > 0);
    CHECK(notOwned == 0);
    // KKT with the pinned controls held (the rest lengths cannot be reached: the end is pinned)
    CHECK(rep.converged);
    CHECK(obj1 < obj0);
}

TEST(ncp_gpu_strand_returns_to_rest_line) {
    REQUIRE_NCP();
    // Fixed ends, rest lengths 0.98 of the even spacing: every edge is stretched, so the unique
    // minimum is the straight evenly spaced line (tension stiffens the lateral directions), with
    // objective sum (wL / 2) ((l - rest) / D)^2 over the edges at l = spacing. The free middle
    // starts jittered. (At exactly the rest length the valley is quartic: not a sharp test.)
    const std::vector<Vec3> line = lineCtrl({0, 0, 0}, {20, 0, 0}, 16);
    std::vector<Vec3> rest = line;
    for (size_t c = 0; c < rest.size(); ++c) rest[c] = line[c] * 0.98;
    std::vector<Vec3> ctrl = line;
    ControlMask freeCtrl(1, std::vector<char>(ctrl.size(), 0));
    for (size_t c = 3; c + 3 < ctrl.size(); ++c) {
        freeCtrl[0][c] = 1;
        ctrl[c] += Vec3{0.1 * std::sin(1.3 * c), 0.15 * std::cos(0.7 * c), 0.1 * std::sin(0.4 * c)};
    }
    CurveSet cs = makeSet({makeYarn(ctrl, 0)});
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 2;
    mp.window = 100;
    Controls P{ctrl}, Prest{rest};
    for (const char* lin : {"matrix-free pcg"}) {
        Controls Q = P;
        gpu::NcpOptions o;
        NcpSolveReport rep;
        gpu::NcpStats st;
        std::string err;
        if (!gpu::solveNcp({cs, Q, Prest, freeCtrl, nullptr, nullptr, mp, {}, {}, o}, rep, &st, &err)) std::printf("    error: %s\n", err.c_str());
        CHECK(err.empty());
        double dmax = 0;
        for (size_t c = 0; c < Q[0].size(); ++c) dmax = std::max(dmax, norm(Q[0][c] - line[c]));
        int edges = 0;  // edges with a free end
        for (size_t c = 0; c + 1 < Q[0].size(); ++c) edges += freeCtrl[0][c] || freeCtrl[0][c + 1];
        const double r = (1.25 - 0.98 * 1.25) / mp.thr, expected = 0.5 * r * r * edges;
        std::printf("    %s: %s, objective %.9g (expected %.9g), max distance to the line %.3g mm, %lld LM steps\n",
                    lin, rep.status.c_str(), rep.score.objective, expected, dmax, st.lmSteps);
        CHECK(rep.converged);
        CHECK(std::fabs(rep.score.objective - expected) <= 1e-6 * expected);
        CHECK(dmax < 1e-3);
    }
}

TEST(ncp_gpu_step_budget) {
    REQUIRE_NCP();
    // NcpSolveParams::maxLmSteps on the strand above: a budget of 3 ends the call UNRESOLVED: step
    // budget after at most one more inner iteration (<= 12 trial steps); a budget the call does not
    // reach changes nothing (identical controls and status to no budget).
    const std::vector<Vec3> line = lineCtrl({0, 0, 0}, {20, 0, 0}, 16);
    std::vector<Vec3> rest = line;
    for (size_t c = 0; c < rest.size(); ++c) rest[c] = line[c] * 0.98;
    std::vector<Vec3> ctrl = line;
    ControlMask freeCtrl(1, std::vector<char>(ctrl.size(), 0));
    for (size_t c = 3; c + 3 < ctrl.size(); ++c) {
        freeCtrl[0][c] = 1;
        ctrl[c] += Vec3{0.1 * std::sin(1.3 * c), 0.15 * std::cos(0.7 * c), 0.1 * std::sin(0.4 * c)};
    }
    CurveSet cs = makeSet({makeYarn(ctrl, 0)});
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 2;
    mp.window = 100;
    const Controls P{ctrl}, Prest{rest};
    for (const char* lin : {"matrix-free pcg"}) {
        gpu::NcpOptions o;
        auto run = [&](long long budget, Controls& Q, NcpSolveReport& rep, gpu::NcpStats& st) {
            Q = P;
            NcpSolveParams sp;
            sp.maxLmSteps = budget;
            std::string err;
            if (!gpu::solveNcp({cs, Q, Prest, freeCtrl, nullptr, nullptr, mp, {}, sp, o}, rep, &st, &err)) std::printf("    error: %s\n", err.c_str());
            CHECK(err.empty());
        };
        Controls Q0, Q3, Qbig;
        NcpSolveReport r0, r3, rbig;
        gpu::NcpStats s0, s3, sbig;
        run(0, Q0, r0, s0);
        run(3, Q3, r3, s3);
        run(s0.lmSteps + 100, Qbig, rbig, sbig);
        std::printf("    %s: no budget %s (%lld LM steps); budget 3: %s (%lld LM steps); budget %lld: %s\n",
                    lin, r0.status.c_str(), s0.lmSteps, r3.status.c_str(), s3.lmSteps,
                    s0.lmSteps + 100, rbig.status.c_str());
        CHECK(r0.converged && !r0.budgetHit && s0.lmSteps > 3);
        CHECK(r3.budgetHit && !r3.converged && !r3.stopped);
        CHECK(r3.status.rfind("UNRESOLVED: step budget", 0) == 0);
        CHECK(s3.lmSteps >= 3 && s3.lmSteps < 3 + 12);
        CHECK(rbig.status == r0.status && sbig.lmSteps == s0.lmSteps && !rbig.budgetHit);
        bool same = true;
        for (size_t c = 0; c < Q0[0].size(); ++c) same = same && Q0[0][c].x == Qbig[0][c].x && Q0[0][c].y == Qbig[0][c].y && Q0[0][c].z == Qbig[0][c].z;
        CHECK(same);
    }
}

TEST(ncp_gpu_pair_pushed_to_hard_contact) {
    REQUIRE_NCP();
    // Strand b bulges towards the straight strand a: 0.6 apart in the middle, 3 mm at the fixed
    // ends (outside the contact). Hard contact 1.0 (no targets); both free in the middle.
    std::vector<Vec3> a = lineCtrl({0, 0, 0}, {16, 0, 0}, 16), b = a;
    for (size_t c = 0; c < b.size(); ++c) {
        const double x = (static_cast<double>(c) - 0.5 * (b.size() - 1)) / (0.5 * (b.size() - 1));
        b[c].y = 0.6 + 2.4 * x * x;
    }
    CurveSet cs = makeSet({makeYarn(a, 0), makeYarn(b, 1)});
    ControlMask freeCtrl(2, std::vector<char>(a.size(), 0));
    for (size_t c = 6; c + 6 < a.size(); ++c) freeCtrl[0][c] = freeCtrl[1][c] = 1;
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 1.5;
    mp.window = 100;
    Controls P{a, b};
    NcpSolveReport rep;
    gpu::NcpStats st;
    std::string err;
    NcpSolveParams sp;
    if (!gpu::solveNcp({cs, P, P, freeCtrl, nullptr, nullptr, mp, {}, sp, {}}, rep, &st, &err)) std::printf("    error: %s\n", err.c_str());
    CHECK(err.empty());
    // brute force: min distance between the two sample polylines (4 samples per span)
    KnitProblem pr(cs, P, P, freeCtrl, nullptr, nullptr, mp);
    pr.refresh(P, nullptr);
    const KnitProblemScore sc = pr.scoreGeometry(P, nullptr, 1e-3);
    std::printf("    %s; contact violation %.3g, curvature violation %.3g, objective %.4g\n", rep.status.c_str(), sc.maxViol[2], sc.maxViol[1],
                sc.objective);
    CHECK(rep.converged);
    CHECK(sc.maxConstraint() <= 1e-3);
}

TEST(ncp_gpu_infeasible_is_unresolved) {
    REQUIRE_NCP();
    // Two strands crossing at distance 0.2 < hard contact 1.0; only one control of one strand is
    // free and its cap is 0.05 mm: no feasible state exists. Never CONVERGED.
    std::vector<Vec3> a = lineCtrl({0, 0, 0}, {16, 0, 0}, 16), b = lineCtrl({8, -8, 0.2}, {8, 8, 0.2}, 16);
    CurveSet cs = makeSet({makeYarn(a, 0), makeYarn(b, 1)});
    ControlMask freeCtrl{std::vector<char>(a.size(), 0), std::vector<char>(b.size(), 0)};
    freeCtrl[0][9] = 1;
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 1.5;
    mp.window = 100;
    mp.maxDisplacement = 0.05;
    Controls P{a, b};
    NcpSolveParams sp;
    sp.maxOuter = 30;
    NcpSolveReport rep;
    gpu::NcpStats st;
    std::string err;
    CHECK(gpu::solveNcp({cs, P, P, freeCtrl, nullptr, nullptr, mp, {}, sp, {}}, rep, &st, &err));
    std::printf("    %s; max violation %.3g\n", rep.status.c_str(), rep.score.maxConstraint());
    CHECK(!rep.converged);
    CHECK(rep.status.rfind("UNRESOLVED", 0) == 0);
    CHECK(rep.score.maxConstraint() > 0.1);
    // the cap holds exactly (projection)
    CHECK(norm(P[0][9] - a[9]) <= 0.05 + 1e-12);
}

TEST(ncp_gpu_infeasible_stops_early) {
    REQUIRE_NCP();
    // The fixture of ncp_gpu_infeasible_is_unresolved with the infeasible rule (window 3): the call
    // ends before the iteration cap, labelled "UNRESOLVED: infeasible", never CONVERGED or STOPPED,
    // and the cap still holds.
    std::vector<Vec3> a = lineCtrl({0, 0, 0}, {16, 0, 0}, 16), b = lineCtrl({8, -8, 0.2}, {8, 8, 0.2}, 16);
    CurveSet cs = makeSet({makeYarn(a, 0), makeYarn(b, 1)});
    ControlMask freeCtrl{std::vector<char>(a.size(), 0), std::vector<char>(b.size(), 0)};
    freeCtrl[0][9] = 1;
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 1.5;
    mp.window = 100;
    mp.maxDisplacement = 0.05;
    NcpSolveParams sp;
    sp.maxOuter = 30;
    sp.infeasibleWindow = 3;
    sp.progressTol = 1e-2;  // as the pipeline: must not fire on an infeasible state
    Controls P{a, b};
    NcpSolveReport rep;
    gpu::NcpStats st;
    std::string err;
    CHECK(gpu::solveNcp({cs, P, P, freeCtrl, nullptr, nullptr, mp, {}, sp, {}}, rep, &st, &err));
    std::printf("    %s after %d outer iterations; max violation %.3g\n", rep.status.c_str(), rep.outer, rep.score.maxConstraint());
    CHECK(!rep.converged && !rep.stopped);
    CHECK(rep.status.rfind("UNRESOLVED: infeasible", 0) == 0);
    CHECK(rep.outer < sp.maxOuter);
    CHECK(rep.rho >= sp.infeasibleRho);
    CHECK(rep.score.maxConstraint() > 0.1);
    CHECK(norm(P[0][9] - a[9]) <= 0.05 + 1e-12);
}

TEST(ncp_gpu_progress_stop_is_feasible_and_labelled) {
    REQUIRE_NCP();
    // The fixture of ncp_gpu_pair_pushed_to_hard_contact with the progress rule. A coarse tolerance
    // (0.5) forces the stop: it happens only at a feasible state (independent score), is labelled
    // STOPPED and not CONVERGED. The pipeline's 1e-2 also ends feasible (STOPPED or CONVERGED).
    std::vector<Vec3> a = lineCtrl({0, 0, 0}, {16, 0, 0}, 16), b = a;
    for (size_t c = 0; c < b.size(); ++c) {
        const double x = (static_cast<double>(c) - 0.5 * (b.size() - 1)) / (0.5 * (b.size() - 1));
        b[c].y = 0.6 + 2.4 * x * x;
    }
    CurveSet cs = makeSet({makeYarn(a, 0), makeYarn(b, 1)});
    ControlMask freeCtrl(2, std::vector<char>(a.size(), 0));
    for (size_t c = 6; c + 6 < a.size(); ++c) freeCtrl[0][c] = freeCtrl[1][c] = 1;
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 1.5;
    mp.window = 100;
    for (double tol : {0.5, 1e-2}) {
        Controls P{a, b};
        NcpSolveParams sp;
        sp.progressTol = tol;
        sp.infeasibleWindow = 3;
        NcpSolveReport rep;
        gpu::NcpStats st;
        std::string err;
        CHECK(gpu::solveNcp({cs, P, P, freeCtrl, nullptr, nullptr, mp, {}, sp, {}}, rep, &st, &err));
        KnitProblem pr(cs, P, P, freeCtrl, nullptr, nullptr, mp);
        pr.refresh(P, nullptr);
        const KnitProblemScore sc = pr.scoreGeometry(P, nullptr, 1e-3);
        std::printf("    tol %g: %s after %d outer iterations; independent max violation %.3g\n", tol, rep.status.c_str(), rep.outer,
                    sc.maxConstraint());
        CHECK(rep.converged != rep.stopped);
        CHECK(rep.status.rfind(rep.stopped ? "STOPPED: feasible" : "CONVERGED", 0) == 0);
        CHECK(sc.maxConstraint() <= sp.tolViol);
        if (tol == 0.5) CHECK(rep.stopped);
    }
}

TEST(ncp_gpu_swept_catches_pass_through) {
    REQUIRE_NCP();
    // Strand b 0.5 above strand a moves to 0.5 below it: both end states are 0.5 apart, the
    // straight motion passes through. The swept test must count the pair; a motion that stays
    // above must not.
    std::vector<Vec3> a = lineCtrl({0, 0, 0}, {16, 0, 0}, 16), b = lineCtrl({8, -8, 0.5}, {8, 8, 0.5}, 16);
    CurveSet cs = makeSet({makeYarn(a, 0), makeYarn(b, 1)});
    ControlMask freeCtrl{std::vector<char>(a.size(), 0), std::vector<char>(b.size(), 1)};
    freeCtrl[1][0] = freeCtrl[1][b.size() - 1] = 0;
    RelaxMoveParams mp;
    mp.thr = 1.0;
    mp.kappaMax = 1.5;
    mp.window = 100;
    Controls P{a, b};
    KnitProblemParams kp;
    kp.delta = 1.0;  // the pair stays covered
    std::string err;
    gpu::NcpDevice dev(rowsCall(cs, P, P, freeCtrl, nullptr, nullptr, mp, kp));
    dev.refresh();
    auto movedBy = [&](double dz) {
        Controls Q = P;
        for (size_t c = 1; c + 1 < b.size(); ++c) Q[1][c].z += dz;
        return Q;
    };
    dev.setTrial(movedBy(-1.0));
    gpu::SweptResult swept = dev.swept(1e-3);
    std::printf("    through: %lld touching, %lld unresolved\n", swept.touching, swept.unresolved);
    CHECK(swept.touching + swept.unresolved > 0);
    dev.setTrial(movedBy(0.3));
    swept = dev.swept(1e-3);
    std::printf("    away: %lld touching, %lld unresolved\n", swept.touching, swept.unresolved);
    CHECK(swept.touching == 0 && swept.unresolved == 0);
}

TEST(ncp_gpu_coverage_brute_force) {
    REQUIRE_NCP();
    // A strongly jittered cylinder patch: contacts enter and leave during the solve. Afterwards
    // every pair of sample segments (brute force, the same-yarn arc window excluded) is at least
    // the hard contact minus the allowance apart, and the refreshes kept the coverage.
    Fixture f = cylinderFixture(24, 5, 0.8, 0, 26, 14, 0.15);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.05;
    kp.wA = 0.01;
    NcpSolveReport rep;
    gpu::NcpStats st;
    std::string err;
    Controls P = f.P;
    CHECK(gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, {}, {}}, rep, &st, &err));
    std::printf("    %s; %lld refreshes, %lld LM steps, max violation %.3g\n", rep.status.c_str(), st.refreshes, st.lmSteps, rep.score.maxConstraint());
    CHECK(rep.score.maxConstraint() <= 1e-3);
    // brute force over the samples of the problem (4 per span)
    const YarnPath& y = f.c.cs.yarns[0];
    const int sps = f.mp.samplesPerSpan, nSp = static_cast<int>(P[0].size()) - 3;
    std::vector<Vec3> pts;
    std::vector<double> s;
    for (int k = 0; k < nSp; ++k)
        for (int j = 0; j < sps; ++j) {
            const double t = double(j) / sps, u = 1 - t;
            const double w[4] = {u * u * u / 6, (3 * t * t * t - 6 * t * t + 4) / 6, (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6, t * t * t / 6};
            pts.push_back(P[0][k] * w[0] + P[0][k + 1] * w[1] + P[0][k + 2] * w[2] + P[0][k + 3] * w[3]);
            s.push_back(y.arcLengthAt(k + t));
        }
    double dmin = 1e300;
    for (size_t i = 0; i + 1 < pts.size(); ++i) {
        bool act = false;
        const int k = static_cast<int>(i) / sps;
        for (int q = 0; q < 4; ++q) act = act || f.freeCtrl[0][k + q];
        if (!act) continue;
        for (size_t j = 0; j + 1 < pts.size(); ++j) {
            if (j == i) continue;
            const double gap = j > i ? s[j] - s[i + 1] : s[i] - s[j + 1];
            if (gap < f.mp.window) continue;
            dmin = std::min(dmin, closestSegmentPoints(pts[i], pts[i + 1], pts[j], pts[j + 1]).distance);
        }
    }
    std::printf("    brute-force min distance %.5f mm (hard contact %.5f)\n", dmin, kp.hardContact);
    CHECK(dmin >= kp.hardContact - 1e-3 * f.mp.thr);
}

TEST(ncp_gpu_coverage_every_trial) {
    REQUIRE_NCP();
    // Coverage at every state the solver evaluates, not only the final one: for every trial Q
    // (accepted, rejected by merit / ownership / the swept test, and the shortened ones in between)
    // and the midpoint of its straight path P -> Q, every pair of active sample segments closer than
    // thr (the same-yarn arc window excluded) was within the broad phase's pad thr + 2 delta at the
    // controls of the last broad phase Xref, i.e. it was a candidate pair with rows. The iteration runs
    // on the device, so the check does too (NcpOptions::auditCoverage: brute force after every trial).
    // A thick yarn on the cylinder's gauge (crowded: many pairs below thr) and a strong jitter (long
    // motions: coverage-limited steps, refreshes).
    Fixture f = cylinderFixture(24, 5, 1.3, 0, 26, 14, 0.4);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.05;
    kp.wA = 0.01;
    kp.delta = 0.15;  // a small coverage radius: more refreshes and coverage-limited steps
    NcpSolveParams sp;
    sp.maxOuter = 8;  // test time; the trials of the first outer iterations move the most
    gpu::NcpOptions opt;
    opt.auditCoverage = true;
    NcpSolveReport rep;
    gpu::NcpStats st;
    std::string err;
    Controls P = f.P;
    CHECK(gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, sp, opt}, rep, &st, &err));
    std::printf("    %s; %lld LM steps (%lld rejected), %lld refreshes, %lld coverage-limited steps; %lld pair states closer than thr checked, %lld uncovered\n",
                rep.status.c_str(), st.lmSteps, st.rejected, st.refreshes, rep.coverageLimited, st.coverageChecks, st.coverageViolations);
    CHECK(st.coverageViolations == 0);
    CHECK(st.coverageChecks > 0);
    CHECK(st.refreshes > 1);  // the fixture must move beyond the first broad phase
    CHECK(st.rejected > 0);   // and evaluate rejected trials
}

TEST(ncp_gpu_scaling_cylinders) {
    REQUIRE_NCP();
    // Total problem growing with the spiral (every stitch but the first and last turn free).
    for (int turns : {6, 12, 24, 48}) {
        const int n = 24;
        Fixture f = cylinderFixture(n, turns, 0.8, 0, n, n * (turns - 2), 0.25);
        KnitProblemParams kp;
        kp.hardContact = f.mp.thr - 0.05;
        kp.wA = 0.01;
        for (const char* lin : {"matrix-free pcg"}) {
            gpu::NcpOptions o;
            NcpSolveReport rep;
            gpu::NcpStats st;
            std::string err;
            Controls P = f.P;
            const auto t0 = std::chrono::steady_clock::now();
            CHECK(gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, {}, o}, rep, &st, &err));
            const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
            std::printf("    %2d turns %6lld unknowns %-5s: %8.1f ms, %s, %lld LM steps, %lld refreshes, contact rows %lld (capacity %lld), peak %.0f MB\n", turns,
                        st.unknowns, lin, ms, rep.status.c_str(), st.lmSteps, st.refreshes,
                        st.maxContactRows, st.capacityContactRows, st.devPeakMB);
            CHECK(rep.score.maxConstraint() <= 1e-3);
            CHECK(st.lmSteps > 0);
        }
    }
}

TEST(ncp_gpu_trial_merit_matches_accepted_rows) {
    REQUIRE_NCP();
    // A trial's contact rows are reduced, not stored, and accept() evaluates
    // them again for the new current state. Both roads must give one merit, or the acceptance test
    // would compare a state with a number that belongs to another.
    // A yarn far thicker than the stitches give room for: contact rows and their targets are active.
    Fixture f = cylinderFixture(24, 5, 3.0, 0.5, 30, 6, 0.05);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.04;
    kp.hardKappa = f.mp.kappaMax * 1.1;
    kp.wA = 0.01;
    KnitProblem pr(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp);
    std::string err;
    gpu::NcpDevice dev(rowsCall(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp));
    dev.refresh();
    const double rho = 10;
    gpu::NcpEvalSummary current, trial, accepted;
    current = dev.linearize(gpu::ControlSet::Current, rho);
    Controls Q = f.P;
    for (int v = 0; v < pr.vars(); ++v) {
        const auto [y, c] = pr.ctrlOf(v);
        Q[y][c] += Vec3{0.03 * std::cos(v), 0.02 * std::sin(2.0 * v), 0.025};
    }
    dev.setTrial(Q);
    trial = dev.linearize(gpu::ControlSet::Trial, rho);
    dev.accept();
    accepted = dev.linearize(gpu::ControlSet::Current, rho);
    std::printf("    merit: before %.12g, trial %.12g, accepted %.12g\n", current.merit, trial.merit, accepted.merit);
    CHECK(trial.merit != current.merit);  // the fixture's motion reached the rows
    // Rows this test exists for: without active contacts the reduction is never exercised.
    CHECK(accepted.maxViol[static_cast<int>(RowFamily::Contact)] > 0);
    CHECK(accepted.maxViol[static_cast<int>(RowFamily::ContactTarget)] > 0);
    CHECK_NEAR(trial.merit, accepted.merit, 1e-9 * accepted.merit);
    CHECK_NEAR(trial.objective, accepted.objective, 1e-9 * std::max(accepted.objective, 1e-12));
    for (int family = 0; family < kRowFamilies; ++family) {
        CHECK_NEAR(trial.maxViol[family], accepted.maxViol[family], 1e-12);
        CHECK(trial.nViol[family] == accepted.nViol[family]);
    }
}

// ---------------------------------------------------------------- session ----
// The solver keeps no device state between relax calls that could go stale: every call
// uploads the controls, the free mask, the rest geometry, the parameters and the binding from
// the host (settle's and repair's rollbacks restore exactly those). The tests below check the
// consequences through settle / certifiedRepair with the ncp numerics.

namespace {

TriMesh planeGrid(int nx, int ny, double c) {
    TriMesh m;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) m.pos.push_back({i * c, j * c, 0});
    auto id = [&](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int a = id(i, j), b = id(i + 1, j), d = id(i + 1, j + 1), e = id(i, j + 1);
            auto uv = [&](int v) { return Vec2{(v % (nx + 1)) * c / 100, (v / (nx + 1)) * c / 100}; };
            m.tri.push_back({a, b, d});
            m.uv.push_back({uv(a), uv(b), uv(d)});
            m.tri.push_back({a, d, e});
            m.uv.push_back({uv(a), uv(d), uv(e)});
        }
    m.computeNormals();
    return m;
}

MeshKnitParams denseTiled() {
    PatternParams pp;
    pp.waleSpacing = 5.5;
    pp.courseSpacing = 3.4;
    pp.yarnDiameter = 1.16;  // the flat dense gauge FAILS clearance: settle has defects to work on
    MeshKnitParams p;
    p.tile.pp = pp;
    p.tile.g.w = pp.waleSpacing;
    p.tile.g.h = pp.courseSpacing;
    p.settle = false;
    return p;
}

} // namespace

TEST(ncp_settle_keeps_the_relaxed_state_whole) {
    REQUIRE_NCP();
    // tests_s6 settle_keeps_the_relaxed_state_whole with the device numerics: a round whose relax
    // result breaks a stitch (mirrored through the surface, its states poisoned) is kept as the relax
    // returned it, controls and projection states alike; nothing is restored from the state before
    // the round, and the validation reports the broken stitch.
    const TriMesh m = planeGrid(60, 40, 1.0);
    MeshKnitParams p = denseTiled();
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k, &err));
    std::unique_ptr<KnitNumerics> ncp = gpu::makeKnitNumerics(&err, gpu::NcpNumericsOptions{});
    CHECK(ncp != nullptr);
    const int y = static_cast<int>(k.t.prog.turns.size() / 2), local = k.t.prog.turns[y].output / 2;
    const int c0 = local * kControlsPerStitch + 1, c1 = c0 + kControlsPerStitch;
    knit_test::SabotagedRound sabotaged(*ncp, y, c0, c1);
    SettleParams sp;
    ValidationReport last;
    const SettleReport rep = settle({k.t.curves, k.t.pp, k.surface.get(), &k.t.bind}, sp, sabotaged, &last);
    std::printf("    %d round(s), %d free stitches; local signature %d FAIL; %s\n", rep.rounds, rep.freeStitches, last.localSignature.fail,
                gpu::knitNumericsStats(*ncp).c_str());
    CHECK(rep.rounds == 1);
    CHECK(last.localSignature.fail > 0);
    CHECK(sabotaged.relaxed.size() == k.t.curves.yarns.size());
    long long differ = 0;
    for (size_t yy = 0; yy < k.t.curves.yarns.size(); ++yy)
        for (size_t c = 0; c < k.t.curves.yarns[yy].curve.ctrl.size(); ++c) {
            const Vec3 &q = k.t.curves.yarns[yy].curve.ctrl[c], &r = sabotaged.relaxed[yy][c];
            differ += !(q.x == r.x && q.y == r.y && q.z == r.z);
            differ += !(k.t.bind.state[yy][c] == sabotaged.relaxedBinding.state[yy][c] && k.t.bind.result[yy][c] == sabotaged.relaxedBinding.result[yy][c]);
        }
    CHECK(differ == 0);
}

TEST(ncp_memory_cap_refuses_the_call) {
    REQUIRE_NCP();
    // A settle round whose solve would exceed the memory cap (1 MB here; the device budget in
    // production) throws DeviceBudgetExceeded with the reason, leaves controls and projection states
    // as they were, and leaves no cap behind: the same round with the default cap runs.
    const TriMesh m = planeGrid(60, 40, 1.0);
    MeshKnitParams p = denseTiled();
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k, &err));
    gpu::NcpNumericsOptions tiny;
    tiny.memoryCapMB = 1;
    std::unique_ptr<KnitNumerics> capped = gpu::makeKnitNumerics(&err, tiny);
    std::unique_ptr<KnitNumerics> normal = gpu::makeKnitNumerics(&err, gpu::NcpNumericsOptions{});
    CHECK(capped && normal);
    const SurfaceBinding b0 = k.t.bind;
    std::vector<std::vector<Vec3>> c0;
    for (const YarnPath& yp : k.t.curves.yarns) c0.push_back(yp.curve.ctrl);
    SettleParams sp;
    std::string what;
    try {
        settle({k.t.curves, k.t.pp, k.surface.get(), &k.t.bind}, sp, *capped);
    } catch (const gpu::DeviceBudgetExceeded& e) {
        what = e.what();
    }
    std::printf("    refused: %s\n", what.c_str());
    CHECK(what.find("cap of 1 MB") != std::string::npos);
    for (size_t yy = 0; yy < c0.size(); ++yy)
        for (size_t c = 0; c < c0[yy].size(); ++c) {
            const Vec3 &a = c0[yy][c], &q = k.t.curves.yarns[yy].curve.ctrl[c];
            CHECK(a.x == q.x && a.y == q.y && a.z == q.z);
            CHECK(k.t.bind.state[yy][c] == b0.state[yy][c] && k.t.bind.result[yy][c] == b0.result[yy][c]);
        }
    const SettleReport rep = settle({k.t.curves, k.t.pp, k.surface.get(), &k.t.bind}, sp, *normal);
    std::printf("    uncapped: %d round(s); %s\n", rep.rounds, gpu::knitNumericsStats(*normal).c_str());
    CHECK(rep.rounds == 1);
}

TEST(ncp_repair_rejected_attempt_keeps_projection_state) {
    REQUIRE_NCP();
    // tests_s6 repair_rejected_attempt_keeps_projection_state.
    const TriMesh m = planeGrid(60, 40, 1.0);
    MeshKnitParams p = denseTiled();
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k, &err));
    std::unique_ptr<KnitNumerics> ncp = gpu::makeKnitNumerics(&err, gpu::NcpNumericsOptions{});
    CHECK(ncp != nullptr);
    const SurfaceBinding b0 = k.t.bind;
    std::vector<std::vector<Vec3>> c0;
    for (const YarnPath& yp : k.t.curves.yarns) c0.push_back(yp.curve.ctrl);
    knit_test::RejectEveryAttempt rejecting(*ncp);
    const RepairReport rr = certifiedRepair({k.t.curves, k.t.pp, k.surface.get(), &k.t.bind}, RepairParams{}, rejecting);
    std::printf("    %s\n", rr.detail.c_str());
    std::printf("    %d attempts rejected; %s\n", rejecting.attempts, gpu::knitNumericsStats(*ncp).c_str());
    CHECK(rejecting.attempts > 0);
    for (size_t yy = 0; yy < c0.size(); ++yy)
        for (size_t c = 0; c < c0[yy].size(); ++c) {
            const Vec3 &a = c0[yy][c], &q = k.t.curves.yarns[yy].curve.ctrl[c];
            CHECK(a.x == q.x && a.y == q.y && a.z == q.z);
            CHECK(k.t.bind.state[yy][c] == b0.state[yy][c] && k.t.bind.result[yy][c] == b0.result[yy][c]);
        }
}

TEST(ncp_binding_only_change_same_controls) {
    REQUIRE_NCP();
    // Two relax calls on one ncp numerics object with the same controls but different projection
    // states: the second result equals a fresh object's result with the second states (bitwise),
    // so nothing of the first call's binding (or its multipliers) leaks into the second.
    const TriMesh m = planeGrid(40, 30, 1.0);
    MeshKnitParams p = denseTiled();
    PartKnit k;
    std::string err;
    CHECK(knitMesh(m, p, knit_test::cpuNumerics(), k, &err));
    Controls P;
    for (const YarnPath& yp : k.t.curves.yarns) P.push_back(yp.curve.ctrl);
    ControlMask freeCtrl;
    for (const auto& yc : P) freeCtrl.push_back(std::vector<char>(yc.size(), 1));
    for (auto& f : freeCtrl) f.front() = f.back() = 0;
    RelaxMoveParams mp;
    const PatternParams& pp = k.t.pp;
    const double Dc = pp.yarnDiameter + pp.clearance;
    mp.thr = Dc + 0.05;
    mp.kappaMax = 1.0 / (1.15 * (pp.radius() + 0.5 * pp.clearance));
    mp.window = 0.5 * kPi * Dc * 1.3;
    SurfaceBinding bA = k.t.bind, bB = k.t.bind;
    // bB: every 7th control without a valid state (its continuation is not Owned: no band row)
    int k7 = 0;
    for (auto& yy : bB.state)
        for (ProjState& s : yy)
            if (++k7 % 7 == 0) s.tri = -1;
    std::unique_ptr<KnitNumerics> n1 = gpu::makeKnitNumerics(&err, gpu::NcpNumericsOptions{});
    std::unique_ptr<KnitNumerics> n2 = gpu::makeKnitNumerics(&err, gpu::NcpNumericsOptions{});
    Controls Q1 = P, Q2 = P, Q3 = P;
    SurfaceBinding s1 = bA, s2 = bB, s3 = bB;
    n1->relax({Q1, P, freeCtrl, k.t.curves, k.surface.get(), mp, &s1});
    n1->relax({Q2, P, freeCtrl, k.t.curves, k.surface.get(), mp, &s2});  // same object, binding B
    n2->relax({Q3, P, freeCtrl, k.t.curves, k.surface.get(), mp, &s3});  // fresh object, binding B
    long long diff = 0;
    for (size_t yy = 0; yy < Q2.size(); ++yy)
        for (size_t c = 0; c < Q2[yy].size(); ++c) diff += !(Q2[yy][c].x == Q3[yy][c].x && Q2[yy][c].y == Q3[yy][c].y && Q2[yy][c].z == Q3[yy][c].z);
    std::printf("    second call vs fresh object: %lld controls differ\n", diff);
    CHECK(diff == 0);
}

TEST(ncp_gpu_capacity_fails_explicitly) {
    REQUIRE_NCP();
    // Every int-indexed total (candidates, rows, row entries) is checked in 64 bits when the call is
    // sized: with the index limit lowered below the fixture's totals the setup fails with "capacity
    // exceeded" (no wrapped count, no truncated list); above them it gives the rows of an unlimited
    // device.
    Fixture f = cylinderFixture(24, 5, 0.8, 0.5, 30, 6, 0.05);
    int failures = 0;
    for (long long limit : {100LL, 1000LL, 4000LL}) {
        gpu::NcpDevice::setIndexLimitForTests(limit);
        try {
            gpu::NcpDevice dev(rowsCall(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, {}));
            dev.refresh();
            std::printf("    limit %lld: ok\n", limit);
        } catch (const gpu::DeviceError& error) {
            const std::string message = error.what();
            std::printf("    limit %lld: %s\n", limit, message.c_str());
            CHECK(message.find("capacity exceeded") != std::string::npos);
            ++failures;
        }
    }
    gpu::NcpDevice::setIndexLimitForTests(0x7fffffffLL);
    CHECK(failures >= 2);
    gpu::NcpDevice::setIndexLimitForTests(10000000);
    gpu::NcpDevice big(rowsCall(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, {}));
    gpu::NcpDevice::setIndexLimitForTests(0x7fffffffLL);
    gpu::NcpDevice ref(rowsCall(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, {}));
    big.refresh();
    ref.refresh();
    CHECK(big.contactRows() == ref.contactRows());
}

TEST(ncp_gpu_contact_capacity_overflow_is_explicit) {
    REQUIRE_NCP();
    // The contact capacity is sized from the starting state with headroom. A refresh that would need
    // more rows is not truncated: a single refresh reports the overflow, and a solve whose refresh
    // overflows and may not grow its capacity (NcpOptions::maxRegrows 0) ends UNRESOLVED ("capacity
    // exceeded") at the last accepted state. The capacity here is one row short of the starting state,
    // so the solve's first refresh overflows: nothing is accepted and the controls come back as they went
    // in (the same check follows every refresh). Allowed to grow (the default), the same solve grows the
    // capacity once and ends exactly as the solve with the default headroom.
    Fixture f = cylinderFixture(24, 5, 1.3, 0, 26, 14, 0.4);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.05;
    kp.wA = 0.01;
    kp.delta = 0.15;
    {
        gpu::NcpDevice dev(rowsCall(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp));
        dev.setContactCapacityForTests(10);
        std::string what;
        try {
            dev.refresh();
        } catch (const gpu::DeviceError& e) {
            what = e.what();
        }
        std::printf("    refresh with 10 rows of capacity: %s\n", what.c_str());
        CHECK(what.find("capacity exceeded") != std::string::npos);
    }
    gpu::NcpOptions noRegrow;
    noRegrow.maxRegrows = 0;
    NcpSolveParams budget;
    budget.maxLmSteps = 40;  // the comparison needs a few steps, not the whole solve
    NcpSolveReport rep, grown, wide;
    gpu::NcpStats st, stGrown, stWide;
    std::string err;
    Controls P = f.P, Pgrown = f.P, Pwide = f.P;
    gpu::NcpDevice::setContactHeadroomForTests(1.0, -1);
    const bool ok = gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, {}, noRegrow}, rep, &st, &err);
    const bool okGrown = gpu::solveNcp({f.c.cs, Pgrown, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, budget, {}}, grown, &stGrown, &err);
    gpu::NcpDevice::setContactHeadroomForTests(-1, -1);
    const bool okWide = gpu::solveNcp({f.c.cs, Pwide, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, budget, {}}, wide, &stWide, &err);
    std::printf("    grown: %s, %lld regrows, %lld LM steps; default headroom: %s, %lld LM steps\n", grown.status.c_str(), stGrown.regrows, stGrown.lmSteps,
                wide.status.c_str(), stWide.lmSteps);
    CHECK(okGrown && okWide);
    CHECK(stGrown.regrows == 1 && stWide.regrows == 0);
    CHECK(Pgrown == Pwide && grown.status == wide.status && stGrown.lmSteps == stWide.lmSteps && stGrown.pcgIterations == stWide.pcgIterations &&
          grown.refreshes == wide.refreshes && stGrown.maxContactRows == stWide.maxContactRows);
    std::printf("    %s; contact rows %lld at the start, max %lld, capacity %lld, %lld LM steps\n", rep.status.c_str(), st.initialContactRows,
                st.maxContactRows, st.capacityContactRows, st.lmSteps);
    CHECK(ok);
    CHECK(rep.status.find("capacity exceeded (contact") != std::string::npos);
    CHECK(!rep.converged && !rep.stopped);
    double moved = 0;
    bool finite = true;
    for (size_t c = 0; c < P[0].size(); ++c) {
        finite = finite && std::isfinite(P[0][c].x) && std::isfinite(P[0][c].y) && std::isfinite(P[0][c].z);
        moved = std::max(moved, norm(P[0][c] - f.P[0][c]));
    }
    CHECK(finite && moved == 0);
    CHECK(st.lmSteps == 0);
}

TEST(ncp_gpu_contact_capacity_overflow_later_keeps_the_accepted_state) {
    REQUIRE_NCP();
    // The overflow at a later refresh (the test above overflows at the first). The free stitches' rest
    // state is shrunk towards their centre, so the yarn tightens and comes closer to its neighbours: the
    // contact rows grow while it relaxes (the reference solve shows it). With the capacity exactly the
    // starting rows (headroom 1, no minimum) a refresh after accepted steps overflows, and the solve ends
    // UNRESOLVED ("capacity exceeded") at the last accepted state when it may not grow the capacity:
    // moved, finite, with a lower objective than the start, and the same state in a second run. Allowed
    // to grow (the default), the solve grows the capacity at that refresh, continues where it stopped,
    // and ends bit for bit as the reference solve (same controls, status, steps, iterations, refreshes).
    Fixture f = cylinderFixture(24, 5, 1.3, 0, 26, 14, 0);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.05;
    kp.wA = 0.01;
    kp.delta = 0.15;
    Vec3 centre{0, 0, 0};
    int count = 0;
    for (size_t k = 0; k < f.P[0].size(); ++k)
        if (f.freeCtrl[0][k]) {
            centre = centre + f.P[0][k];
            ++count;
        }
    centre = centre * (1.0 / count);
    Controls rest = f.P;
    for (size_t k = 0; k < rest[0].size(); ++k)
        if (f.freeCtrl[0][k]) rest[0][k] = centre + (rest[0][k] - centre) * 0.8;
    NcpSolveParams budget;
    budget.maxLmSteps = 300;  // bounds the reference solve; the overflow comes first
    gpu::NcpOptions noRegrow;
    noRegrow.maxRegrows = 0;
    gpu::NcpOptions regrow;
    const gpu::NcpOptions* options = &noRegrow;
    auto solve = [&](Controls& P, NcpSolveReport& rep, gpu::NcpStats& st) {
        std::string err;
        P = f.P;
        const bool ok = gpu::solveNcp({f.c.cs, P, rest, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, budget, *options}, rep, &st, &err);
        if (!ok) std::printf("    %s\n", err.c_str());
        return ok;
    };
    auto objective = [&](const Controls& P) {
        Controls Q = P;  // rowsCall takes the controls it may write
        gpu::NcpDevice dev(rowsCall(f.c.cs, Q, rest, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp));
        dev.refresh();
        return dev.linearize(gpu::ControlSet::Current, 10).objective;
    };
    Controls Pref;
    NcpSolveReport ref;
    gpu::NcpStats rs;
    CHECK(solve(Pref, ref, rs));
    std::printf("    reference: %s; contact rows %lld at the start, max %lld\n", ref.status.c_str(), rs.initialContactRows, rs.maxContactRows);
    CHECK(rs.maxContactRows > rs.initialContactRows);  // the fixture exercises the later overflow
    gpu::NcpDevice::setContactHeadroomForTests(1.0, 0);
    Controls P1, P2, P3;
    NcpSolveReport rep1, rep2, rep3;
    gpu::NcpStats st1, st2, st3;
    const bool ok1 = solve(P1, rep1, st1), ok2 = solve(P2, rep2, st2);
    options = &regrow;
    const bool ok3 = solve(P3, rep3, st3);
    gpu::NcpDevice::setContactHeadroomForTests(-1, -1);
    std::printf("    grown: %s after %lld LM steps, %lld refreshes, %lld regrows (capacity %lld)\n", rep3.status.c_str(), st3.lmSteps, rep3.refreshes, st3.regrows,
                st3.capacityContactRows);
    CHECK(ok3 && st3.regrows >= 1);
    CHECK(P3 == Pref && rep3.status == ref.status && st3.lmSteps == rs.lmSteps && st3.pcgIterations == rs.pcgIterations && rep3.refreshes == ref.refreshes &&
          st3.maxContactRows == rs.maxContactRows);
    std::printf("    capacity %lld: %s after %lld LM steps, %lld refreshes\n", st1.capacityContactRows, rep1.status.c_str(), st1.lmSteps, rep1.refreshes);
    CHECK(ok1 && ok2);
    CHECK(st1.capacityContactRows == rs.initialContactRows);
    CHECK(rep1.status.find("capacity exceeded (contact") != std::string::npos);
    CHECK(!rep1.converged && !rep1.stopped && !rep1.budgetHit);
    CHECK(st1.lmSteps > 0 && rep1.refreshes > 1);  // not the first refresh: steps were accepted before it
    double moved = 0;
    bool finite = true;
    for (size_t c = 0; c < P1[0].size(); ++c) {
        finite = finite && std::isfinite(P1[0][c].x) && std::isfinite(P1[0][c].y) && std::isfinite(P1[0][c].z);
        moved = std::max(moved, norm(P1[0][c] - f.P[0][c]));
    }
    const double e0 = objective(f.P), e1 = objective(P1);
    std::printf("    moved %.4f mm; objective %.6g at the start, %.6g at the result\n", moved, e0, e1);
    CHECK(finite && moved > 0);
    CHECK(e1 < e0);
    CHECK(P1 == P2 && rep2.status == rep1.status && st2.lmSteps == st1.lmSteps);
}

TEST(ncp_gpu_rows_match_cpu_hybrid_contacts) {
    REQUIRE_NCP();
    // The hybrid contact model (knit_problem.h ContactModel): nearly parallel candidates as two
    // point-segment rows, the others as one segment row; the device classifies like the CPU (same
    // keys, type bits 57..58) and evaluates the same values and gradients.
    for (auto model : {KnitProblemParams::ContactModel::Hybrid, KnitProblemParams::ContactModel::PointSegment, KnitProblemParams::ContactModel::Point}) {
        Fixture f = cylinderFixture(24, 5, 1.3, 0.5, 30, 6, 0.05);  // thick: contacts with a term in the merit
        KnitProblemParams kp;
        kp.hardContact = f.mp.thr - 0.04;
        kp.contactModel = model;
        KnitProblem pr(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp);
        pr.refresh(f.P, nullptr);
        std::string err;
        gpu::NcpDevice dev(rowsCall(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp));
        dev.refresh();
        gpu::NcpEvalSummary sm;
        sm = dev.linearize(gpu::ControlSet::Current, 10);
        std::vector<gpu::NcpDevice::DebugRow> rows;
        rows = dev.debugRows();
        long long types[4] = {0, 0, 0, 0};
        for (const auto& d : rows)
            if (d.family == static_cast<int>(RowFamily::Contact)) ++types[(d.key >> 57) & 3];
        const RowCompare rc = compareRows(pr, f.P, rows, f.mp.thr + 2 * kp.delta - 1e-6);
        std::printf("    model %d: contact rows by type: segment %lld, point %lld, point-segment %lld\n", static_cast<int>(model), types[0], types[1],
                    types[2]);
        reportRows("rows", rc);
        checkRows(rc);
        std::printf("      contact rows with gradients compared: %lld\n", rc.contactGradientsCompared);
        CHECK(rc.contactGradientsCompared > 0);
        if (model == KnitProblemParams::ContactModel::Hybrid) CHECK(types[0] > 0 && types[2] > 0);
    }
}

TEST(ncp_early_stop_rules) {
    // ncpEarlyStop: the progress rule stops only feasible states and uses the
    // magnitude of the change, so a rise of the objective beyond the tolerance continues; the
    // infeasible rule needs rho >= infeasibleRho and a violation that did not halve over the window.
    NcpSolveParams sp;
    using E = NcpEarlyStop;
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = true, .objPrev = 1.0, .obj = 1.0, .rho = 1e4, .viol = 0, .violBack = 1.0}) == E::None);  // both rules off by default
    sp.progressTol = 1e-2;
    sp.infeasibleWindow = 3;
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = true, .objPrev = 1.0, .obj = 0.995, .rho = 10, .viol = 0, .violBack = -1}) == E::Progress);  // fell by 0.5 %
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = true, .objPrev = 1.0, .obj = 1.005, .rho = 10, .viol = 0, .violBack = -1}) == E::Progress);  // rose by 0.5 %
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = true, .objPrev = 1.0, .obj = 0.95, .rho = 10, .viol = 0, .violBack = -1}) == E::None);       // fell by 5 %: still progressing
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = true, .objPrev = 1.0, .obj = 1.2, .rho = 10, .viol = 0, .violBack = -1}) == E::None);        // rose by 20 %: continue (the corrected case)
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = false, .objPrev = 1.0, .obj = 1.0, .rho = 10, .viol = 0.5, .violBack = -1}) == E::None);     // not feasible: no progress stop
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = true, .objPrev = -1, .obj = 1.0, .rho = 10, .viol = 0, .violBack = -1}) == E::None);         // no previous outer iteration yet
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = false, .objPrev = 1.0, .obj = 1.0, .rho = 1e3, .viol = 0.6, .violBack = 1.0}) == E::Infeasible);  // violation 1.0 -> 0.6 at rho 1e3
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = false, .objPrev = 1.0, .obj = 1.0, .rho = 1e3, .viol = 0.4, .violBack = 1.0}) == E::None);        // halved: still converging
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = false, .objPrev = 1.0, .obj = 1.0, .rho = 640, .viol = 0.9, .violBack = 1.0}) == E::None);        // rho below infeasibleRho
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = false, .objPrev = 1.0, .obj = 1.0, .rho = 1e4, .viol = 0.9, .violBack = -1}) == E::None);         // window not filled yet
    CHECK(gpu::ncpEarlyStop(sp, {.feasible = true, .objPrev = 1.0, .obj = 1.0, .rho = 1e4, .viol = 0, .violBack = 1.0}) == E::Progress);       // feasible: never Infeasible
}

// --- contact memory: the compact contact storage against
// the independent CPU rows, at a linearization held fixed.
namespace {

// The merit, objective, gradient, diagonal blocks and undamped product H p of the model at the rows'
// controls, from the independent CPU rows (KnitProblem::linearize with activeBelow 0: exactly the rows
// with a term in the merit), assembled as ncp.h states the model: objective rows (edges wL, anchors wA,
// targets wT) (w / 2) r^2; hard rows (rho / 2) max(0, g + lambda / rho)^2, with lambda = max(0, rho g) in
// the second pass (one multiplier update from zero at the same state); Gauss-Newton blocks w g g^T, and a
// stretched edge's kk (I - u u^T) between its controls.
struct HostOperator {
    double merit = 0, objective = 0;
    std::vector<double> grad, hblock, product;
    std::vector<double> gradTerms;  // per component: the sum of the terms' magnitudes |wres g| (the sum's rounding scale)
};
HostOperator hostOperator(const std::vector<Row>& rows, int vars, const KnitProblemParams& kp, double D, double rho, bool updatedMultipliers,
                          const std::vector<double>& p) {
    HostOperator h;
    h.grad.assign(size_t(3 * vars), 0.0);
    h.hblock.assign(size_t(6 * vars), 0.0);
    h.product.assign(size_t(3 * vars), 0.0);
    h.gradTerms.assign(size_t(3 * vars), 0.0);
    for (const Row& r : rows) {
        double w = 0, wres = 0, m = 0;
        bool objective = true;
        switch (r.family) {
        case RowFamily::Edge: w = kp.wL; break;
        case RowFamily::Anchor: w = kp.wA; break;
        case RowFamily::ContactTarget:
        case RowFamily::CurvatureTarget: w = kp.wT; break;
        default: objective = false; break;
        }
        if (objective) {
            wres = w * r.value;
            m = 0.5 * w * r.value * r.value;
            h.objective += m;
        } else {
            const double lam = updatedMultipliers ? std::max(0.0, rho * r.value) : 0.0;
            const double phi = r.value + lam / rho;
            if (phi > 0) {
                w = rho;
                wres = rho * phi;
                m = 0.5 * rho * phi * phi;
            }
        }
        h.merit += m;
        double jp = 0;
        for (int e = 0; e < r.n; ++e) {
            const Vec3& g = r.e[e].g;
            const int a = r.e[e].var;
            h.grad[3 * a] += wres * g.x;
            h.grad[3 * a + 1] += wres * g.y;
            h.grad[3 * a + 2] += wres * g.z;
            h.gradTerms[3 * a] += std::fabs(wres * g.x);
            h.gradTerms[3 * a + 1] += std::fabs(wres * g.y);
            h.gradTerms[3 * a + 2] += std::fabs(wres * g.z);
            double* b = &h.hblock[6 * a];
            b[0] += w * g.x * g.x; b[1] += w * g.x * g.y; b[2] += w * g.x * g.z;
            b[3] += w * g.y * g.y; b[4] += w * g.y * g.z; b[5] += w * g.z * g.z;
            jp += g.x * p[3 * a] + g.y * p[3 * a + 1] + g.z * p[3 * a + 2];
        }
        for (int e = 0; e < r.n; ++e) {
            const int a = r.e[e].var;
            h.product[3 * a] += r.e[e].g.x * w * jp;
            h.product[3 * a + 1] += r.e[e].g.y * w * jp;
            h.product[3 * a + 2] += r.e[e].g.z * w * jp;
        }
        if (r.family != RowFamily::Edge || r.n == 0) continue;
        const double l = r.aux, gl = norm(r.e[0].g);
        if (!(r.value > 0 && l > 1e-12 && gl > 0)) continue;
        const Vec3 u = r.e[0].g / gl;
        const double stiff = kp.wL * r.value / (l * D);
        const int a0 = r.e[0].var, a1 = r.n > 1 ? r.e[1].var : -1;
        Vec3 dp{p[3 * a0], p[3 * a0 + 1], p[3 * a0 + 2]};
        if (a1 >= 0) dp = dp - Vec3{p[3 * a1], p[3 * a1 + 1], p[3 * a1 + 2]};
        const Vec3 s = (dp - u * dot(u, dp)) * stiff;
        for (int k = 0; k < (a1 >= 0 ? 2 : 1); ++k) {
            const int a = k == 0 ? a0 : a1;
            const double sign = k == 0 ? 1.0 : -1.0;
            h.product[3 * a] += sign * s.x;
            h.product[3 * a + 1] += sign * s.y;
            h.product[3 * a + 2] += sign * s.z;
            double* b = &h.hblock[6 * a];
            b[0] += stiff * (1 - u.x * u.x); b[1] -= stiff * u.x * u.y; b[2] -= stiff * u.x * u.z;
            b[3] += stiff * (1 - u.y * u.y); b[4] -= stiff * u.y * u.z; b[5] += stiff * (1 - u.z * u.z);
        }
    }
    return h;
}

// The largest difference of a device quantity from the host's, relative to the host's largest magnitude.
double relativeDifference(const std::vector<double>& device, const std::vector<double>& host) {
    double d = 0, scale = 1e-300;
    for (size_t i = 0; i < host.size(); ++i) {
        d = std::max(d, std::fabs(device[i] - host[i]));
        scale = std::max(scale, std::fabs(host[i]));
    }
    return d / scale;
}
double relativeDifference(double device, double host) { return std::fabs(device - host) / std::max(std::fabs(host), 1e-300); }

// The device's frozen passes against the host's from the CPU rows; returns the largest relative difference.
double compareFrozen(const std::vector<gpu::NcpDevice::FrozenPass>& passes, const std::vector<Row>& rows, int vars, const KnitProblemParams& kp, double D,
                     double rho, const std::vector<double>& p) {
    double worst = 0;
    for (size_t k = 0; k < passes.size(); ++k) {
        const HostOperator h = hostOperator(rows, vars, kp, D, rho, k == 1, p);
        const gpu::NcpDevice::FrozenPass& f = passes[k];
        const double dm = relativeDifference(f.eval.merit, h.merit), dob = relativeDifference(f.eval.objective, h.objective);
        const double dg = relativeDifference(f.grad, h.grad), dh = relativeDifference(f.hblock, h.hblock), dq = relativeDifference(f.product, h.product);
        std::printf("    pass %zu (lambda %s): merit %.10g vs host %.10g; relative differences: merit %.2e, objective %.2e, gradient %.2e, blocks %.2e, "
                    "product %.2e; %lld contact rows, %lld active rows\n",
                    k, k == 0 ? "0" : "max(0, rho g)", f.eval.merit, h.merit, dm, dob, dg, dh, dq, f.contactRows, f.activeRows);
        worst = std::max({worst, dm, dob, dg, dh, dq});
        // Where the gradient differs most, against the magnitude of the terms summed there.
        size_t at = 0;
        double termRel = 0;
        for (size_t i = 0; i < h.grad.size(); ++i) {
            if (std::fabs(f.grad[i] - h.grad[i]) > std::fabs(f.grad[at] - h.grad[at])) at = i;
            termRel = std::max(termRel, std::fabs(f.grad[i] - h.grad[i]) / std::max(h.gradTerms[i], 1e-300));
        }
        std::printf("      gradient: largest difference %.3e at variable %zu (device %.12g, host %.12g, terms' magnitude %.6g); largest difference relative "
                    "to its terms' magnitude %.2e\n",
                    std::fabs(f.grad[at] - h.grad[at]), at / 3, f.grad[at], h.grad[at], h.gradTerms[at], termRel);
    }
    return worst;
}

// The device's contact row keys against the CPU's (the same refresh state): every CPU row present; a
// device row the CPU did not emit must lie beyond the CPU's pad less the float boxes' guard.
struct KeyCompare {
    long long cpu = 0, device = 0, missing = 0, extra = 0, extraInsidePad = 0, duplicates = 0;
};
KeyCompare compareKeys(const KnitProblem& pr, const Controls& P, std::vector<uint64_t> device, double padMinusGuard) {
    KeyCompare k;
    std::vector<uint64_t> cpu;
    const auto& pairs = pr.contactRowPairs();
    const auto& types = pr.contactRowTypes();
    cpu.reserve(pairs.size());
    for (size_t q = 0; q < pairs.size(); ++q)
        cpu.push_back((2ull << 60) ^ (static_cast<uint64_t>(pairs[q].first) << 30) ^ static_cast<uint64_t>(pairs[q].second) ^ (static_cast<uint64_t>(types[q]) << 57));
    std::sort(cpu.begin(), cpu.end());
    std::sort(device.begin(), device.end());
    k.cpu = static_cast<long long>(cpu.size());
    k.device = static_cast<long long>(device.size());
    const auto uniqueEnd = std::unique(device.begin(), device.end());
    k.duplicates = static_cast<long long>(device.end() - uniqueEnd);
    device.erase(uniqueEnd, device.end());
    std::vector<uint64_t> missing, extra;
    std::set_difference(cpu.begin(), cpu.end(), device.begin(), device.end(), std::back_inserter(missing));
    std::set_difference(device.begin(), device.end(), cpu.begin(), cpu.end(), std::back_inserter(extra));
    k.missing = static_cast<long long>(missing.size());
    k.extra = static_cast<long long>(extra.size());
    for (uint64_t key : extra) {
        const size_t j = key & ((1ull << 30) - 1), i = (key >> 30) & ((1ull << 27) - 1);
        const int type = static_cast<int>((key >> 57) & 3);
        if (pr.contactRowDistance(P, i, j, type) < padMinusGuard) ++k.extraInsidePad;
    }
    std::printf("    contact keys: %lld CPU, %lld device (%lld duplicates), %lld missing, %lld extra (%lld inside the pad)\n", k.cpu, k.device, k.duplicates,
                k.missing, k.extra, k.extraInsidePad);
    return k;
}

std::vector<double> frozenDirection(int vars) {  // knit_numerics.cpp replayNcpCall's direction
    std::vector<double> p(size_t(3 * vars));
    for (size_t i = 0; i < p.size(); ++i) p[i] = std::sin(0.37 * static_cast<double>(i) + 1.0);
    return p;
}

} // namespace

TEST(ncp_gpu_frozen_operator_matches_cpu_reference) {
    REQUIRE_NCP();
    // The compact contact storage (gradients formed per active slot, merits summed from the values,
    // multipliers carried by key in the active entries' memory) against the independent CPU rows at one
    // linearization held fixed: merit, objective, gradient, diagonal blocks and the product H p, before
    // and after one multiplier update and a second refresh. The thick cylinder: contacts with terms.
    Fixture f = cylinderFixture(24, 5, 1.3, 0.5, 30, 6, 0.05);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.04;
    kp.hardKappa = f.mp.kappaMax * 1.1;
    kp.wA = 0.01;
    const double rho = 10;
    KnitProblem pr(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp);
    pr.refresh(f.P, nullptr);
    std::vector<Row> rows;
    pr.linearize(f.P, rows, 0);
    gpu::NcpDevice dev(rowsCall(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp));
    const std::vector<double> p = frozenDirection(dev.vars());
    const std::vector<gpu::NcpDevice::FrozenPass> passes = dev.frozenLinearizations(rho, p);
    const double worst = compareFrozen(passes, rows, pr.vars(), kp, pr.D(), rho, p);
    CHECK(passes.size() == 2 && passes[0].activeRows > 0);
    CHECK(passes[1].eval.merit > passes[0].eval.merit);  // violated rows: the multipliers moved (pass 1 exercises them)
    CHECK(worst <= 1e-9);
    const KeyCompare k = compareKeys(pr, f.P, dev.contactKeys(), f.mp.thr + 2 * kp.delta - 1e-6);
    CHECK(k.missing == 0 && k.extraInsidePad == 0 && k.duplicates == 0 && k.cpu > 1000);
}

TEST(ncp_gpu_dropped_rows_are_detected) {
    REQUIRE_NCP();
    // The omission control: 50 contact rows overwritten with a copy of their predecessor keep the row
    // count, and the key comparison against the CPU names all 50 as missing (a count check would pass).
    Fixture f = cylinderFixture(24, 5, 1.3, 0.5, 30, 6, 0.05);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.04;
    KnitProblem pr(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp);
    pr.refresh(f.P, nullptr);
    gpu::NcpDevice dev(rowsCall(f.c.cs, f.P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp));
    dev.refresh();
    const double pad = f.mp.thr + 2 * kp.delta - 1e-6;
    const long long rows = dev.contactRows();
    const KeyCompare intact = compareKeys(pr, f.P, dev.contactKeys(), pad);
    dev.overwriteContactRowsForTests(100, 50);
    const KeyCompare broken = compareKeys(pr, f.P, dev.contactKeys(), pad);
    CHECK(intact.missing == 0 && intact.extraInsidePad == 0 && rows > 200);
    CHECK(dev.contactRows() == rows);
    CHECK(broken.missing == 50 && broken.duplicates == 50 && broken.extra == intact.extra);
}

TEST(ncp_gpu_frozen_operator_matches_cpu_reference_on_call) {
    REQUIRE_NCP();
    // The same comparison on a captured relax call (KNIT_NCP_FROZEN_CALL: a .ncpcall file, e.g. the real
    // arm's target stage): the complete contact key set, and the frozen operator from the CPU rows with a
    // term in the merit. A full-size check within the host memory bound, where the legacy layout's
    // allocation does not fit.
    const char* path = std::getenv("KNIT_NCP_FROZEN_CALL");
    if (!path || !*path) SKIP("KNIT_NCP_FROZEN_CALL not set");
    gpu::NcpCallRecord rec;
    std::string err;
    CHECK(gpu::readNcpCall(path, rec, &err));
    if (!err.empty()) std::printf("    %s\n", err.c_str());
    rec.move.call.phase = rec.phase.c_str();
    std::unique_ptr<MeshSurface> surface;
    if (rec.hasMesh) surface = gpu::makeMeshSurface(rec.mesh, nullptr, &err);
    const double rho = 10;
    // The device first (its active rows size the CPU's row list), released before the CPU rows are built.
    std::vector<gpu::NcpDevice::FrozenPass> passes;
    std::vector<uint64_t> deviceKeys;
    std::vector<gpu::NcpDevice::DebugRow> deviceGeneral;
    std::vector<double> p;
    {
        Controls P = rec.P;
        gpu::NcpDevice dev(rowsCall(rec.curves, P, rec.Prest, rec.freeCtrl, surface.get(), nullptr, rec.move, rec.problem));
        p = frozenDirection(dev.vars());
        passes = dev.frozenLinearizations(rho, p);
        deviceKeys = dev.contactKeys();
        deviceGeneral = dev.debugRows(gpu::NcpDevice::kDebugGeneralRowsOnly);
        std::printf("    device memory report: %s\n", dev.stats().memoryLine().c_str());
    }
    const auto t0 = std::chrono::steady_clock::now();
    KnitProblem pr(rec.curves, rec.P, rec.Prest, rec.freeCtrl, surface.get(), nullptr, rec.move, rec.problem);
    pr.refresh(rec.P, nullptr);
    std::vector<Row> rows;
    rows.reserve(size_t(passes[0].activeRows) * 2 + 100000);  // targets are separate CPU rows: at most one per active contact row
    pr.linearize(rec.P, rows, 0);
    const double cpuMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    std::printf("    CPU reference: %zu contact rows, %zu rows with a term, %.0f ms\n", pr.contactRowPairs().size(), rows.size(), cpuMs);
    std::printf("    the whole operator from the CPU rows (information: includes the general rows' CPU / device differences):\n");
    compareFrozen(passes, rows, pr.vars(), rec.problem, pr.D(), rho, p);
    // The contact path alone: the CPU's curvature rows replaced by the device's (general rows, whose kernel the
    // contact storage does not touch; compared row by row below), everything else from the CPU.
    std::vector<Row> mixed;
    mixed.reserve(rows.size());
    std::map<uint64_t, const Row*> cpuCurvature;
    for (const Row& r : rows) {
        if (r.family == RowFamily::Curvature || r.family == RowFamily::CurvatureTarget) cpuCurvature[r.key] = &r;
        else mixed.push_back(r);
    }
    double curvValue = 0, curvGrad = 0;
    long long curvMatched = 0;
    for (const auto& d : deviceGeneral) {
        const RowFamily f = static_cast<RowFamily>(d.family);
        if (f != RowFamily::Curvature && f != RowFamily::CurvatureTarget) continue;
        if (f == RowFamily::Curvature && !(d.value > 0)) continue;  // the CPU's activeBelow 0
        Row r;
        r.family = f;
        r.value = d.value;
        r.key = d.key;
        r.n = d.n;
        for (int e = 0; e < d.n; ++e) r.e[e] = {d.var[e], Vec3{d.g[e][0], d.g[e][1], d.g[e][2]}};
        mixed.push_back(r);
        const auto it = cpuCurvature.find(d.key);
        if (it == cpuCurvature.end()) continue;
        ++curvMatched;
        curvValue = std::max(curvValue, std::fabs(it->second->value - d.value));
        for (int e = 0; e < d.n; ++e)
            curvGrad = std::max({curvGrad, std::fabs(it->second->e[e].g.x - d.g[e][0]), std::fabs(it->second->e[e].g.y - d.g[e][1]), std::fabs(it->second->e[e].g.z - d.g[e][2])});
    }
    std::printf("    curvature rows with a term, CPU against device: %zu CPU, %lld matched, max |value diff| %.3g, max |gradient diff| %.3g\n", cpuCurvature.size(),
                curvMatched, curvValue, curvGrad);
    std::printf("    the contact path (device curvature rows):\n");
    const double worst = compareFrozen(passes, mixed, pr.vars(), rec.problem, pr.D(), rho, p);
    // Where the gradient differs most: the rows touching that variable, CPU against device by key.
    const HostOperator h0 = hostOperator(rows, pr.vars(), rec.problem, pr.D(), rho, false, p);
    size_t at = 0;
    for (size_t i = 0; i < h0.grad.size(); ++i)
        if (std::fabs(passes[0].grad[i] - h0.grad[i]) > std::fabs(passes[0].grad[at] - h0.grad[at])) at = i;
    const int worstVar = static_cast<int>(at / 3);
    std::vector<gpu::NcpDevice::DebugRow> touching;
    {
        Controls P = rec.P;
        gpu::NcpDevice dev(rowsCall(rec.curves, P, rec.Prest, rec.freeCtrl, surface.get(), nullptr, rec.move, rec.problem));
        dev.frozenLinearizations(rho, p);
        touching = dev.debugRows(worstVar);
    }
    std::map<uint64_t, const gpu::NcpDevice::DebugRow*> byKey;
    for (const auto& d : touching) byKey[d.key] = &d;
    const auto [wy, wc] = pr.ctrlOf(worstVar);
    std::printf("    variable %d (yarn %d control %d): %zu device rows touch it%s", worstVar, wy, wc, touching.size(), "\n");
    for (const Row& r : rows) {
        int e = -1;
        for (int k = 0; k < r.n; ++k)
            if (r.e[k].var == worstVar) e = k;
        if (e < 0) continue;
        const auto it = byKey.find(r.key);
        double dg = -1, dv = 0;
        if (it != byKey.end()) {
            const auto& d = *it->second;
            dv = d.value;
            for (int k = 0; k < d.n; ++k)
                if (d.var[k] == worstVar)
                    dg = std::max({std::fabs(d.g[k][0] - r.e[e].g.x), std::fabs(d.g[k][1] - r.e[e].g.y), std::fabs(d.g[k][2] - r.e[e].g.z)});
        }
        if (dg > 1e-9 || it == byKey.end() || std::fabs(dv - r.value) > 1e-9)
            std::printf("      %-16s key %016llx: value CPU %.12g device %.12g; gradient difference %.3e%s%s", toString(r.family), static_cast<unsigned long long>(r.key), r.value,
                        dv, dg, it == byKey.end() ? " (no device row)" : "", "\n");
    }
    CHECK(worst <= 1e-9);
    const KeyCompare k = compareKeys(pr, rec.P, deviceKeys, rec.move.thr + 2 * rec.problem.delta - 1e-6);
    CHECK(k.missing == 0 && k.extraInsidePad == 0 && k.duplicates == 0);
}

TEST(ncp_gpu_budget_refusal_cleans_up) {
    REQUIRE_NCP();
    // The memory cap of a solve (AllocationCap; the port sets min(caller's cap, granted budget)):
    //   below the setup's need: refused (false, capHit), and every tracked byte is given back;
    //   enough for the setup but not for the capacity a refresh needs later: the regrow is refused the
    //   same way (never a truncated row list), and the bytes are given back;
    //   the same solve without the cap then runs.
    Fixture f = cylinderFixture(24, 5, 1.3, 0, 26, 14, 0.4);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.05;
    kp.wA = 0.01;
    kp.delta = 0.15;
    NcpSolveParams budget;
    budget.maxLmSteps = 40;
    const size_t before = gpu::memStats().current;
    // The setup's exact bytes at a capacity one row short of the starting rows (the regrow then needs more).
    gpu::NcpDevice::setContactHeadroomForTests(1.0, -1);
    size_t setupBytes = 0;
    {
        Controls Q = f.P;
        gpu::NcpDevice dev(rowsCall(f.c.cs, Q, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp));
        setupBytes = gpu::memStats().current - before;
    }
    CHECK(gpu::memStats().current == before);
    auto cappedSolve = [&](size_t capBytes, std::string& err) {
        Controls P = f.P;
        NcpSolveReport rep;
        gpu::NcpStats st;
        gpu::AllocationCap cap(capBytes);
        const bool ok = gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, budget, {}}, rep, &st, &err);
        return ok || !gpu::AllocationCap::capHit() ? std::string(ok ? "solved" : "failed") : std::string("refused");
    };
    std::string err1, err2;
    const std::string small = cappedSolve(before + setupBytes / 2, err1);
    std::printf("    cap at half the setup: %s (%s); tracked %zu -> %zu bytes\n", small.c_str(), err1.c_str(), before, gpu::memStats().current);
    CHECK(small == "refused" && gpu::memStats().current == before);
    const std::string regrow = cappedSolve(before + setupBytes, err2);
    std::printf("    cap at the setup's bytes (the regrow needs more): %s (%s); tracked %zu bytes\n", regrow.c_str(), err2.c_str(), gpu::memStats().current);
    CHECK(regrow == "refused" && gpu::memStats().current == before);
    gpu::NcpDevice::setContactHeadroomForTests(-1, -1);
    std::string err3;
    Controls P = f.P;
    NcpSolveReport rep;
    gpu::NcpStats st;
    CHECK(gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, budget, {}}, rep, &st, &err3));
    CHECK(gpu::memStats().current == before);
}

// --- live contact storage (ContactStorage::LiveRows):
// the experimental storage against the compact one (v5) and the independent CPU rows. Each test names the
// defect it fails on.
namespace {

const gpu::NcpOptions& liveOptions() {
    static gpu::NcpOptions options = [] {
        gpu::NcpOptions o;
        o.contactStorage = ContactStorage::LiveRows;
        return o;
    }();
    return options;
}

gpu::NcpCall liveRowsCall(const CurveSet& cs, Controls& P, const Controls& rest, const ControlMask& freeCtrl, const BaseSurface* surface, const RelaxMoveParams& mp,
                          const KnitProblemParams& kp) {
    static const NcpSolveParams solve;
    return {cs, P, rest, freeCtrl, surface, nullptr, mp, kp, solve, liveOptions()};
}

// The thick jittered cylinder of the frozen-operator tests: contacts with terms, near-parallel rows.
struct ThickCylinder {
    Fixture f = cylinderFixture(24, 5, 1.3, 0.5, 30, 6, 0.05);
    KnitProblemParams kp;
    ThickCylinder() {
        kp.hardContact = f.mp.thr - 0.04;
        kp.hardKappa = f.mp.kappaMax * 1.1;
        kp.wA = 0.01;
        kp.contactModel = KnitProblemParams::ContactModel::Hybrid;  // the port's model (both row kinds)
    }
    // ncp_refresh.cu nearCutoff: the float pad thr + 2 delta + 1e-3, plus max(0, hardC - thr).
    double nearCutoff() const { return static_cast<double>(static_cast<float>(f.mp.thr + 2 * kp.delta) + 1e-3f) + std::max(0.0, kp.hardContact - f.mp.thr); }
};

// Two antiparallel strands (one span per sample, segments 0.5 mm) at a distance just inside the near
// cutoff, the second shifted half a segment along the axis: each overlapping segment pair is near, while a
// near-parallel row's own distance (from its first point to the other segment) reaches past the cutoff.
struct AntiparallelStrands {
    CurveSet cs;
    Controls P;
    ControlMask freeCtrl;
    RelaxMoveParams mp;
    KnitProblemParams kp;
    double cutoff = 0, gap = 0;
};
AntiparallelStrands antiparallelStrands() {
    AntiparallelStrands s;
    s.mp.thr = 1.2;
    s.mp.kappaMax = 1;
    s.mp.window = 0.5;
    s.mp.samplesPerSpan = 1;
    s.kp.contactModel = KnitProblemParams::ContactModel::Hybrid;  // the port's model: near-parallel rows point - segment
    s.cutoff = static_cast<double>(static_cast<float>(s.mp.thr + 2 * s.kp.delta) + 1e-3f);  // hardC = thr here
    s.gap = s.cutoff - 0.002;
    s.cs = makeSet({makeYarn(lineCtrl({-3, 0, 0}, {3, 0, 0}, 12), 0), makeYarn(lineCtrl({3.25, s.gap, 0}, {-2.75, s.gap, 0}, 12), 1)});
    finalizeCurves(s.cs);
    s.P = {s.cs.yarns[0].curve.ctrl, s.cs.yarns[1].curve.ctrl};
    s.freeCtrl.resize(2);
    for (int y = 0; y < 2; ++y) s.freeCtrl[y].assign(s.P[y].size(), 1);
    return s;
}

// The CPU's contact rows split by the full segment distance at P (contactRowDistance with type 0) against
// the device's row keys: near rows the device lacks, far rows it keeps, device rows the CPU has not.
struct ResidentCompare {
    long long near = 0, far = 0, nearMissing = 0, farKept = 0, extra = 0;
};
ResidentCompare compareResident(const KnitProblem& pr, const Controls& P, std::vector<uint64_t> device, double cutoff) {
    ResidentCompare r;
    std::sort(device.begin(), device.end());
    const auto& pairs = pr.contactRowPairs();
    const auto& types = pr.contactRowTypes();
    std::vector<uint64_t> cpu;
    for (size_t q = 0; q < pairs.size(); ++q) {
        const uint64_t key = (2ull << 60) ^ (static_cast<uint64_t>(pairs[q].first) << 30) ^ static_cast<uint64_t>(pairs[q].second) ^ (static_cast<uint64_t>(types[q]) << 57);
        cpu.push_back(key);
        const double segments = pr.contactRowDistance(P, pairs[q].first, pairs[q].second, 0);
        const bool present = std::binary_search(device.begin(), device.end(), key);
        if (segments < cutoff - 1e-6) {
            ++r.near;
            if (!present) ++r.nearMissing;
        } else if (segments > cutoff + 1e-6) {
            ++r.far;
            if (present) ++r.farKept;
        }
    }
    std::sort(cpu.begin(), cpu.end());
    for (uint64_t key : device)
        if (!std::binary_search(cpu.begin(), cpu.end(), key)) ++r.extra;
    std::printf("    rows by segment distance at P: %lld near (%lld missing on the device), %lld far (%lld kept); %lld device rows the CPU lacks\n", r.near, r.nearMissing,
                r.far, r.farKept, r.extra);
    return r;
}

} // namespace

TEST(ncp_gpu_live_rows_match_compact_frozen) {
    REQUIRE_NCP();
    // At a linearization held fixed (and after one multiplier update and a refresh that carries the
    // multipliers by key), the live storage gives the compact storage's active rows and entries exactly and
    // its merit, objective, gradient, diagonal blocks and operator product to rounding (the merit sums
    // group fewer rows; the warp assembly and gather add each variable's entries in a fixed tree instead
    // of one sequence), with fewer rows; both agree with the CPU rows. Fails when a near row is left out,
    // a multiplier is carried to the wrong row, or an entry is lost or doubled in the warp kernels.
    ThickCylinder t;
    const double rho = 10;
    Controls Pc = t.f.P, Pl = t.f.P;
    gpu::NcpDevice compact(rowsCall(t.f.c.cs, Pc, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), nullptr, t.f.mp, t.kp));
    gpu::NcpDevice live(liveRowsCall(t.f.c.cs, Pl, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), t.f.mp, t.kp));
    const std::vector<double> p = frozenDirection(compact.vars());
    const std::vector<gpu::NcpDevice::FrozenPass> a = compact.frozenLinearizations(rho, p), b = live.frozenLinearizations(rho, p);
    for (size_t k = 0; k < 2; ++k) {
        std::printf("    pass %zu: merit %.17g / %.17g; rows %lld / %lld; active rows %lld / %lld\n", k, a[k].eval.merit, b[k].eval.merit, a[k].contactRows, b[k].contactRows,
                    a[k].activeRows, b[k].activeRows);
        CHECK(relativeDifference(b[k].eval.merit, a[k].eval.merit) <= 1e-14 && relativeDifference(b[k].eval.objective, a[k].eval.objective) <= 1e-14);
        const double dg = relativeDifference(b[k].grad, a[k].grad), dh = relativeDifference(b[k].hblock, a[k].hblock), dq = relativeDifference(b[k].product, a[k].product);
        std::printf("    pass %zu against compact: gradient %.2e, blocks %.2e, product %.2e (relative)\n", k, dg, dh, dq);
        CHECK(dg <= 1e-13 && dh <= 1e-13 && dq <= 1e-13);
        CHECK(a[k].activeRows == b[k].activeRows && a[k].activeEntries == b[k].activeEntries);
        CHECK(b[k].contactRows < a[k].contactRows);
    }
    CHECK(b[1].eval.merit > b[0].eval.merit);  // the multipliers moved: pass 1 exercises the carry-over
    KnitProblem pr(t.f.c.cs, t.f.P, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), nullptr, t.f.mp, t.kp);
    pr.refresh(t.f.P, nullptr);
    std::vector<Row> rows;
    pr.linearize(t.f.P, rows, 0);
    CHECK(compareFrozen(b, rows, pr.vars(), t.kp, pr.D(), rho, p) <= 1e-9);
}

TEST(ncp_gpu_live_rows_keep_every_near_pair) {
    REQUIRE_NCP();
    // The cull uses the two full sample segments: every CPU row whose segments are within the near cutoff
    // is a device row, every row beyond it is left out, and the device has no row the CPU lacks. A cull by
    // the row's own distance (point - segment for near-parallel pairs, which can exceed the segments'
    // distance) leaves out near rows and fails here: the antiparallel strands are built to have such rows.
    {
        ThickCylinder t;
        KnitProblem pr(t.f.c.cs, t.f.P, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), nullptr, t.f.mp, t.kp);
        pr.refresh(t.f.P, nullptr);
        Controls Pl = t.f.P;
        gpu::NcpDevice live(liveRowsCall(t.f.c.cs, Pl, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), t.f.mp, t.kp));
        live.refresh();
        const ResidentCompare r = compareResident(pr, t.f.P, live.contactKeys(), t.nearCutoff());
        CHECK(r.near > 1000 && r.far > 1000);  // both kinds occur on this fixture
        CHECK(r.nearMissing == 0 && r.farKept == 0 && r.extra == 0);
    }
    AntiparallelStrands s = antiparallelStrands();
    KnitProblem pr(s.cs, s.P, s.P, s.freeCtrl, nullptr, nullptr, s.mp, s.kp);
    pr.refresh(s.P, nullptr);
    Controls Pl = s.P;
    gpu::NcpDevice live(liveRowsCall(s.cs, Pl, s.P, s.freeCtrl, nullptr, s.mp, s.kp));
    live.refresh();
    long long rowFarSegmentsNear = 0;
    const auto& pairs = pr.contactRowPairs();
    const auto& types = pr.contactRowTypes();
    for (size_t q = 0; q < pairs.size(); ++q)
        if (pr.contactRowDistance(s.P, pairs[q].first, pairs[q].second, types[q]) > s.cutoff &&
            pr.contactRowDistance(s.P, pairs[q].first, pairs[q].second, 0) < s.cutoff - 1e-6)
            ++rowFarSegmentsNear;
    std::printf("    antiparallel strands: %lld rows beyond the cutoff by their own distance but near by their segments\n", rowFarSegmentsNear);
    const ResidentCompare r = compareResident(pr, s.P, live.contactKeys(), s.cutoff);
    CHECK(rowFarSegmentsNear > 10);
    CHECK(r.nearMissing == 0 && r.farKept == 0 && r.extra == 0);
}

TEST(ncp_gpu_live_swept_test_matches_compact) {
    REQUIRE_NCP();
    // The swept test runs over the row list. (1) The antiparallel strands, just inside the near cutoff,
    // moved towards each other by just under delta each (the coverage limit), tested at the largest eps
    // the cull covers (max(hardC, thr) + guard = cutoff - 2 delta): touching and unresolved counts equal the
    // compact storage's. No pair touches: the conservative advancement's first step already reaches within
    // eps / 2 of a pair's start distance, and a pair that starts near the cutoff closes by at most 2 delta,
    // so a pair the cull leaves out can never touch (which ncp_gpu_live_rows_keep_every_near_pair relies on
    // and checks key by key). (2) A strand moved through another, clear of it at both ends: both storages
    // catch it. Fails when the live storage loses a row of the row list the swept test reads.
    AntiparallelStrands s = antiparallelStrands();
    Controls Pc = s.P, Pl = s.P;
    gpu::NcpDevice compact(rowsCall(s.cs, Pc, s.P, s.freeCtrl, nullptr, nullptr, s.mp, s.kp));
    gpu::NcpDevice live(liveRowsCall(s.cs, Pl, s.P, s.freeCtrl, nullptr, s.mp, s.kp));
    compact.refresh();
    live.refresh();
    const double eps = s.cutoff - 2 * s.kp.delta;
    Controls Q1 = s.P;
    for (Vec3& c : Q1[0]) c.y += 0.999 * s.kp.delta;
    for (Vec3& c : Q1[1]) c.y -= 0.999 * s.kp.delta;
    compact.setTrial(Q1);
    live.setTrial(Q1);
    const gpu::SweptResult tc = compact.swept(eps), tl = live.swept(eps);
    std::printf("    antiparallel strands moved together by 2 x 0.999 delta: touching %lld / %lld, unresolved %lld / %lld\n", tc.touching, tl.touching, tc.unresolved,
                tl.unresolved);
    CHECK(tc.touching == tl.touching && tc.unresolved == tl.unresolved);
    // (2) Two crossing strands 0.12 mm apart; the trial moves the second 0.24 mm through the first.
    const double gap = 0.12;
    CurveSet cs = makeSet({makeYarn(lineCtrl({-3, 0, 0}, {3, 0, 0}, 12), 0), makeYarn(lineCtrl({0, -3, gap}, {0, 3, gap}, 12), 1)});
    finalizeCurves(cs);
    RelaxMoveParams mp;
    mp.thr = 0.1;
    mp.kappaMax = 1;
    mp.window = 0.5;
    KnitProblemParams kp;
    ControlMask freeCtrl(2);
    Controls P{cs.yarns[0].curve.ctrl, cs.yarns[1].curve.ctrl};
    for (int y = 0; y < 2; ++y) freeCtrl[y].assign(P[y].size(), 1);
    Controls Pc2 = P, Pl2 = P;
    gpu::NcpDevice compact2(rowsCall(cs, Pc2, P, freeCtrl, nullptr, nullptr, mp, kp));
    gpu::NcpDevice live2(liveRowsCall(cs, Pl2, P, freeCtrl, nullptr, mp, kp));
    compact2.refresh();
    live2.refresh();
    Controls Q = P;
    for (Vec3& c : Q[1]) c.z -= 2 * gap;
    compact2.setTrial(Q);
    live2.setTrial(Q);
    const gpu::SweptResult sc = compact2.swept(1e-3 * mp.thr), sl = live2.swept(1e-3 * mp.thr);
    std::printf("    strand moved through another: touching %lld / %lld\n", sc.touching, sl.touching);
    CHECK(sl.touching > 0 && sl.touching == sc.touching);
}

TEST(ncp_gpu_live_trial_gives_new_contacts_their_place) {
    REQUIRE_NCP();
    // A row without a live place enters a trial with multiplier 0; once the trial is accepted, the rows that
    // gained a term are live. After a trial that pushes contacts together (within the coverage limit),
    // the accepted state's live rows hold every compact row with a term, with the same value and gradient.
    // Fails when the trial reduction skips the rows that are not live, or the accept does not renumber.
    ThickCylinder t;
    const double rho = 10;
    Controls Pc = t.f.P, Pl = t.f.P;
    gpu::NcpDevice compact(rowsCall(t.f.c.cs, Pc, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), nullptr, t.f.mp, t.kp));
    gpu::NcpDevice live(liveRowsCall(t.f.c.cs, Pl, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), t.f.mp, t.kp));
    compact.refresh();
    live.refresh();
    const long long liveBefore = live.debugRows().size();
    compact.linearize(gpu::ControlSet::Current, rho);
    live.linearize(gpu::ControlSet::Current, rho);
    // Every free control 0.9 delta towards the cylinder's axis (inside the coverage limit, as every solver
    // step is): the rows compress.
    Controls Q = t.f.P;
    for (size_t c = 0; c < Q[0].size(); ++c) {
        if (!t.f.freeCtrl[0][c]) continue;
        Vec3 radial{Q[0][c].x, Q[0][c].y, 0};
        Q[0][c] -= normalized(radial) * (0.9 * t.kp.delta);
    }
    compact.setTrial(Q);
    live.setTrial(Q);
    const gpu::NcpEvalSummary tc = compact.linearize(gpu::ControlSet::Trial, rho), tl = live.linearize(gpu::ControlSet::Trial, rho);
    std::printf("    trial merit %.17g / %.17g (rows without a live place enter with multiplier 0)\n", tc.merit, tl.merit);
    CHECK(std::fabs(tc.merit - tl.merit) <= 1e-12 * std::fabs(tc.merit));
    compact.accept();
    live.accept();
    std::map<uint64_t, gpu::NcpDevice::DebugRow> liveRows;
    for (const gpu::NcpDevice::DebugRow& d : live.debugRows())
        if (d.family == static_cast<int>(RowFamily::Contact)) liveRows[d.key] = d;
    long long withTerm = 0, missing = 0, different = 0;
    for (const gpu::NcpDevice::DebugRow& d : compact.debugRows()) {
        if (d.family != static_cast<int>(RowFamily::Contact) || !d.hasGradient) continue;
        ++withTerm;
        const auto it = liveRows.find(d.key);
        if (it == liveRows.end()) {
            ++missing;
            continue;
        }
        bool same = it->second.value == d.value && it->second.n == d.n;
        for (int e = 0; e < d.n && same; ++e)
            for (int k = 0; k < 3; ++k) same = same && it->second.g[e][k] == d.g[e][k];
        if (!same) ++different;
    }
    std::printf("    live rows %lld before, %zu after the accepted trial; compact rows with a term %lld: %lld missing, %lld different\n", liveBefore, liveRows.size(),
                withTerm, missing, different);
    CHECK(withTerm > static_cast<long long>(liveBefore) / 2 && missing == 0 && different == 0);
}

TEST(ncp_gpu_live_multiplier_follows_the_row_key) {
    REQUIRE_NCP();
    // v5's carry-over rule: a multiplier carries over only when the same (i, j, type) is a row again.
    // After one multiplier update (frozen pass 1), the yarn moves far enough that pairs leave the broad
    // phase and near-parallel pairs change type; after the refresh both storages give the same merit and
    // the same live / term rows. Fails when the live storage keeps a multiplier v5 drops (or the reverse).
    ThickCylinder t;
    const double rho = 10;
    Controls Pc = t.f.P, Pl = t.f.P;
    gpu::NcpDevice compact(rowsCall(t.f.c.cs, Pc, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), nullptr, t.f.mp, t.kp));
    gpu::NcpDevice live(liveRowsCall(t.f.c.cs, Pl, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), t.f.mp, t.kp));
    const std::vector<double> p = frozenDirection(compact.vars());
    compact.frozenLinearizations(rho, p);  // ends with multipliers updated once and carried by a refresh
    live.frozenLinearizations(rho, p);
    // A twist: every free control turned 0.35 mm about the axis and lifted 0.2 mm (beyond delta: the next
    // refresh rebuilds the rows from scratch, types change where segments turn).
    Controls Q = t.f.P;
    for (size_t c = 0; c < Q[0].size(); ++c) {
        if (!t.f.freeCtrl[0][c]) continue;
        const Vec3 radial = normalized(Vec3{Q[0][c].x, Q[0][c].y, 0});
        Q[0][c] += Vec3{-radial.y, radial.x, 0} * (0.35 * std::sin(0.4 * c)) + Vec3{0, 0, 0.2 * std::cos(0.3 * c)};
    }
    for (gpu::NcpDevice* d : {&compact, &live}) {
        d->setTrial(Q);
        d->linearize(gpu::ControlSet::Trial, rho);
        d->accept();
        d->refresh();
    }
    const gpu::NcpEvalSummary ec = compact.linearize(gpu::ControlSet::Current, rho), el = live.linearize(gpu::ControlSet::Current, rho);
    std::map<uint64_t, double> termRows;
    for (const gpu::NcpDevice::DebugRow& d : compact.debugRows())
        if (d.family == static_cast<int>(RowFamily::Contact) && d.hasGradient) termRows[d.key] = d.value;
    long long missing = 0;
    for (const auto& [key, value] : termRows) {
        (void)value;
        bool found = false;
        for (const gpu::NcpDevice::DebugRow& d : live.debugRows())
            if (d.key == key) found = true;
        if (!found) ++missing;
        if (missing > 10) break;
    }
    std::printf("    after the twist and a refresh: merit %.17g / %.17g, max contact violation %.6g / %.6g; %zu compact rows with a term, %lld not live\n", ec.merit,
                el.merit, ec.maxViol[static_cast<int>(RowFamily::Contact)], el.maxViol[static_cast<int>(RowFamily::Contact)], termRows.size(), missing);
    CHECK(std::fabs(ec.merit - el.merit) <= 1e-12 * std::fabs(ec.merit));
    CHECK(ec.maxViol[static_cast<int>(RowFamily::Contact)] == el.maxViol[static_cast<int>(RowFamily::Contact)]);
    CHECK(!termRows.empty() && missing == 0);
}

TEST(ncp_gpu_live_dropped_row_is_detected) {
    REQUIRE_NCP();
    // The omission control on the live storage: 50 resident rows overwritten with a copy of their
    // predecessor keep the row count, and the comparison against the CPU's near rows names all 50.
    ThickCylinder t;
    KnitProblem pr(t.f.c.cs, t.f.P, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), nullptr, t.f.mp, t.kp);
    pr.refresh(t.f.P, nullptr);
    Controls Pl = t.f.P;
    gpu::NcpDevice live(liveRowsCall(t.f.c.cs, Pl, t.f.P, t.f.freeCtrl, t.f.c.surf.get(), t.f.mp, t.kp));
    live.refresh();
    const ResidentCompare intact = compareResident(pr, t.f.P, live.contactKeys(), t.nearCutoff());
    live.overwriteContactRowsForTests(100, 50);
    const ResidentCompare broken = compareResident(pr, t.f.P, live.contactKeys(), t.nearCutoff());
    CHECK(intact.nearMissing == 0);
    CHECK(broken.nearMissing == 50);
}

TEST(ncp_gpu_live_capacity_overflow_is_explicit) {
    REQUIRE_NCP();
    // A live list smaller than the rows with a term: the solve ends UNRESOLVED "capacity exceeded (live
    // contact rows ...)" at the last accepted state (here the start: nothing moved), without a regrow.
    Fixture f = cylinderFixture(24, 5, 1.3, 0, 26, 14, 0.4);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.05;
    kp.wA = 0.01;
    NcpSolveParams budget;
    budget.maxLmSteps = 40;
    gpu::NcpDevice::setLiveHeadroomForTests(0.0, 1);  // a live list of one row
    Controls P = f.P;
    NcpSolveReport rep;
    gpu::NcpStats st;
    std::string err;
    const bool ok = gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, budget, liveOptions()}, rep, &st, &err);
    gpu::NcpDevice::setLiveHeadroomForTests(-1, -1);
    std::printf("    %s; %lld LM steps, %lld regrows\n", rep.status.c_str(), st.lmSteps, st.regrows);
    CHECK(ok && !rep.converged && !rep.stopped);
    CHECK(rep.status.find("capacity exceeded (live contact rows") != std::string::npos);
    bool unmoved = true;
    for (size_t c = 0; c < P[0].size(); ++c) unmoved = unmoved && P[0][c].x == f.P[0][c].x && P[0][c].y == f.P[0][c].y && P[0][c].z == f.P[0][c].z;
    CHECK(st.lmSteps == 0 && st.regrows == 0 && unmoved);
}

TEST(ncp_gpu_live_refusal_and_cancel_clean_up) {
    REQUIRE_NCP();
    // The live storage under the memory cap and the cancel token: a cap below the setup's need is a refusal
    // that gives back every tracked byte; a cancelled solve throws Cancelled and gives back every byte; the
    // same solve then runs and stays within the cap.
    Fixture f = cylinderFixture(24, 5, 1.3, 0, 26, 14, 0.4);
    KnitProblemParams kp;
    kp.hardContact = f.mp.thr - 0.05;
    kp.wA = 0.01;
    NcpSolveParams budget;
    budget.maxLmSteps = 40;
    const size_t before = gpu::memStats().current;
    size_t setupBytes = 0;
    {
        Controls Q = f.P;
        gpu::NcpDevice dev(liveRowsCall(f.c.cs, Q, f.P, f.freeCtrl, f.c.surf.get(), f.mp, kp));
        setupBytes = gpu::memStats().current - before;
    }
    std::string err;
    {
        Controls P = f.P;
        NcpSolveReport rep;
        gpu::NcpStats st;
        gpu::AllocationCap cap(before + setupBytes / 2);
        const bool ok = gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, budget, liveOptions()}, rep, &st, &err);
        std::printf("    cap at half the setup: %s (%s); tracked %zu -> %zu bytes\n", ok ? "solved" : "refused", err.c_str(), before, gpu::memStats().current);
        CHECK(!ok && gpu::AllocationCap::capHit());
    }
    CHECK(gpu::memStats().current == before);
    bool cancelled = false;
    gpu::AllocationCap uncapped(0);  // no cap from here on; clears the refusal's flag for the tests that follow
    {
        CancelToken token;
        token.cancel();
        ThreadCancelScope scope(&token);
        Controls P = f.P;
        NcpSolveReport rep;
        gpu::NcpStats st;
        try {
            gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, budget, liveOptions()}, rep, &st, &err);
        } catch (const Cancelled&) {
            cancelled = true;
        }
    }
    std::printf("    cancelled before the launch: %s; tracked %zu bytes\n", cancelled ? "Cancelled thrown" : "NOT cancelled", gpu::memStats().current);
    CHECK(cancelled && gpu::memStats().current == before && gpu::abortFlagValue() == 0);
    Controls P = f.P;
    NcpSolveReport rep;
    gpu::NcpStats st;
    CHECK(gpu::solveNcp({f.c.cs, P, f.P, f.freeCtrl, f.c.surf.get(), nullptr, f.mp, kp, budget, liveOptions()}, rep, &st, &err));
    std::printf("    the same solve after both: %s, %lld LM steps, %lld up / %lld down bytes\n", rep.status.c_str(), st.lmSteps, st.uploadBytes, st.downloadBytes);
    CHECK(gpu::memStats().current == before && !gpu::AllocationCap::capHit());
}

TEST(ncp_gpu_live_identity_names_the_storage) {
    // The cache identity: CompactRows with the sampled curvature keeps the v5 text (every pin unchanged); LiveRows and
    // the span-maximum curvature rows are other solvers, and so is their combination.
    const ContactStorage compact = ContactStorage::CompactRows, live = ContactStorage::LiveRows;
    const CurvatureRows samples = CurvatureRows::Samples, spanMax = CurvatureRows::SamplesAndSpanMaximum;
    const std::string v5 = gpu::ncpSolverIdentity(compact, samples);
    CHECK(v5.rfind("ncp-gpu v5 (", 0) == 0 && v5.back() == ')');
    CHECK(gpu::ncpSolverIdentity(live, samples) != v5);
    CHECK(gpu::ncpSolverIdentity(live, samples).find("live contact rows") != std::string::npos);
    CHECK(gpu::ncpSolverIdentity(compact, spanMax).find("span-maximum curvature rows") != std::string::npos);
    CHECK(gpu::ncpSolverIdentity(live, spanMax) != gpu::ncpSolverIdentity(live, samples));
    CHECK(gpu::ncpSolverIdentity(live, spanMax) != gpu::ncpSolverIdentity(compact, spanMax));
    if (!envVar("KNIT_NCP_CONTACT_STORAGE") && !envVar("KNIT_NCP_CURVATURE_ROWS")) CHECK(gpu::ncpSolverIdentity() == v5);
}

// The smallest curvature radius of a uniform cubic B-spline yarn over every span, at 257 parameters per span.
double denseMinimumRadius(const std::vector<Vec3>& P) {
    double rMin = 1e300;
    for (size_t k = 0; k + 3 < P.size(); ++k)
        for (int j = 0; j <= 256; ++j) {
            const double t = j / 256.0, s = 1 - t;
            const double b1[4] = {-0.5 * s * s, 1.5 * t * t - 2 * t, -1.5 * t * t + t + 0.5, 0.5 * t * t};
            const double b2[4] = {s, 3 * t - 2, 1 - 3 * t, t};
            Vec3 d1{0, 0, 0}, d2{0, 0, 0};
            for (int i = 0; i < 4; ++i) {
                d1 += P[k + i] * b1[i];
                d2 += P[k + i] * b2[i];
            }
            const double L = norm(d1), C = norm(cross(d1, d2));
            if (C > 0) rMin = std::min(rMin, L * L * L / C);
        }
    return rMin;
}

TEST(ncp_gpu_span_maximum_row_opens_a_bend_between_samples) {
    REQUIRE_NCP();
    // tube-144's last call left a bend of radius 0.613 mm between two curvature samples of radius 0.724 mm, below
    // the hard radius 0.649 mm (the control polygon doubles back by 0.042 mm there; report section 9). These are
    // fourteen of its controls (relative to control 11515; the bend is in span 5), the middle eight free, the rest
    // lengths those of the edges as they are. With the fixed samples only, the start is feasible and optimal, so the
    // bend stays. With the span-maximum row, the solve must open it to the hard radius.
    const std::vector<Vec3> ctrl = {
        {0.543980, -0.987178, 0.385932},   {0.474342, -0.728644, 0.485289},   {0.359667, -0.466455, 0.485642},
        {0.229567, -0.230862, 0.380160},   {0.102986, -0.070556, 0.199330},   {0.000000, 0.000000, 0.000000},
        {-0.133444, 0.014034, -0.168117},  {-0.107258, 0.011250, -0.135155},  {-0.256339, 0.025001, -0.321782},
        {-0.410978, -0.063873, -0.646375}, {-0.496274, 0.017659, -1.019338},  {-0.405872, 0.163863, -1.352856},
        {-0.186502, 0.276722, -1.592665},  {0.058953, 0.342711, -1.823829}};
    const double hardRadius = 0.649;
    CHECK(denseMinimumRadius(ctrl) < 0.62);
    ControlMask freeCtrl(1, std::vector<char>(ctrl.size(), 0));
    for (size_t c = 3; c + 3 < ctrl.size(); ++c) freeCtrl[0][c] = 1;
    CurveSet cs = makeSet({makeYarn(ctrl, 0)});
    RelaxMoveParams mp;
    mp.thr = 0.05;  // no contact: one short yarn, the window spans it
    mp.kappaMax = 1 / hardRadius;
    mp.window = 100;
    KnitProblemParams kp;
    kp.curvatureSamples = 9;  // as knit_numerics.cpp
    const Controls rest{ctrl};
    for (CurvatureRows rows : {CurvatureRows::Samples, CurvatureRows::SamplesAndSpanMaximum}) {
        Controls P{ctrl};
        gpu::NcpOptions o;
        o.curvatureRows = rows;
        NcpSolveReport rep;
        gpu::NcpStats st;
        std::string err;
        CHECK(gpu::solveNcp({cs, P, rest, freeCtrl, nullptr, nullptr, mp, kp, {}, o}, rep, &st, &err));
        const double r = denseMinimumRadius(P[0]);
        const bool spanMax = rows == CurvatureRows::SamplesAndSpanMaximum;
        std::printf("    %s: %s, %lld LM steps, smallest radius %.4f mm (hard %.3f), objective %.3g\n", spanMax ? "samples and span maximum" : "samples only",
                    rep.status.c_str(), st.lmSteps, r, hardRadius, rep.score.objective);
        CHECK(rep.converged);
        if (spanMax) {
            CHECK(r >= hardRadius / (1 + 1e-3) - 1e-6);  // the hard rows' allowance tolViol = 1e-3
        } else {
            CHECK(r < 0.62);  // the defect this row closes: feasible at every sample, the bend still there
        }
    }
}
