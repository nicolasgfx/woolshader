// The tests' CPU knitting numerics (cpu_numerics.h): the former CPU relaxation (projections on the
// free controls), the boundary checks, the lattice classification and gate (d), kept only to drive the
// orchestration on machines without a device and as the reference of the device's batches.
#include "cpu_numerics.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <set>

#include "knitting/knit/piece_knit.h"
#include "model/cancel.h"
#include "model/spatial_hash.h"

namespace knit_test {

using namespace knit;

namespace {

using Clock = std::chrono::steady_clock;
using Ctrl = std::vector<std::vector<Vec3>>;
using Mask = std::vector<std::vector<char>>;
using MoveParams = RelaxMoveParams;

double msSince(Clock::time_point t0) { return std::chrono::duration<double, std::milli>(Clock::now() - t0).count(); }

void basisWeights(double t, double w[4]) {
    const double s = 1 - t;
    w[0] = s * s * s / 6;
    w[1] = (3 * t * t * t - 6 * t * t + 4) / 6;
    w[2] = (-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6;
    w[3] = t * t * t / 6;
}

struct Sample {
    int y, k;      // yarn, span (controls k..k+3)
    double w[4];
    double s;      // arc length (input curve) for the same-yarn window
};

Vec3 evalSample(const Ctrl& P, const Sample& sm) {
    const auto& c = P[sm.y];
    return c[sm.k] * sm.w[0] + c[sm.k + 1] * sm.w[1] + c[sm.k + 2] * sm.w[2] + c[sm.k + 3] * sm.w[3];
}

// Move the curve point of sample sm by delta, distributed over its FREE controls
// (minimum-norm PBD correction).
void applyToSample(Ctrl& P, const Mask& freeCtrl, const Sample& sm, const Vec3& delta) {
    double sw = 0;
    for (int i = 0; i < 4; ++i)
        if (freeCtrl[sm.y][sm.k + i]) sw += sm.w[i] * sm.w[i];
    if (sw <= 0) return;
    for (int i = 0; i < 4; ++i)
        if (freeCtrl[sm.y][sm.k + i]) P[sm.y][sm.k + i] += delta * (sm.w[i] / sw);
}


// Gauss-Seidel projections on the free controls of P: bend limiter per span, soft
// control-edge lengths (rest from Prest), contacts against all nearby spans, band
// around the base surface, and (optionally) a cap on each control's displacement.
void relaxFree(Ctrl& P, const Ctrl& Prest, const Mask& freeCtrl, const CurveSet& cs, const BaseSurface* surface,
               const MoveParams& mp, SurfaceBinding* bind = nullptr) {
    const int Y = static_cast<int>(P.size());
    const Ctrl Pstart = P;
    std::vector<Sample> samples;
    std::vector<char> active;
    for (int y = 0; y < Y; ++y) {
        const YarnPath& yp = cs.yarns[y];
        for (int k = 0; k < yp.curve.numSpans(); ++k) {
            bool act = false;
            for (int i = 0; i < 4; ++i) act = act || freeCtrl[y][k + i];
            for (int j = 0; j < mp.samplesPerSpan; ++j) {
                Sample sm;
                sm.y = y;
                sm.k = k;
                basisWeights(double(j) / mp.samplesPerSpan, sm.w);
                sm.s = yp.arcLengthAt(k + double(j) / mp.samplesPerSpan);
                samples.push_back(sm);
                active.push_back(act);
            }
        }
    }
    auto segOk = [&](size_t i) { return i + 1 < samples.size() && samples[i + 1].y == samples[i].y; };
    std::vector<std::vector<double>> rest(Y);
    for (int y = 0; y < Y; ++y)
        for (size_t c = 0; c + 1 < P[y].size(); ++c) rest[y].push_back(norm(Prest[y][c + 1] - Prest[y][c]));
    std::vector<std::pair<size_t, size_t>> cand;
    const double thr = mp.thr;
    for (int it = 0; it < mp.iterations; ++it) {
        cancelCheckpoint();
        if (it % mp.refresh == 0) {
            // The refresh is independent per sample, per segment and per control (each writes its
            // own entries): in parallel with mp.exec, results in the serial order.
            std::vector<Vec3> pos(samples.size());
            parallelFor(static_cast<int64_t>(samples.size()), 4096, [&](int64_t b, int64_t e, int64_t) {
                for (int64_t i = b; i < e; ++i) pos[i] = evalSample(P, samples[i]);
            }, mp.exec);
            SpatialHash3 h(thr + 1.0);
            std::vector<size_t> segId;
            auto box = [&](size_t i, Vec3& lo, Vec3& hi) {
                lo = {std::min(pos[i].x, pos[i + 1].x), std::min(pos[i].y, pos[i + 1].y), std::min(pos[i].z, pos[i + 1].z)};
                hi = {std::max(pos[i].x, pos[i + 1].x), std::max(pos[i].y, pos[i + 1].y), std::max(pos[i].z, pos[i + 1].z)};
            };
            for (size_t i = 0; i < samples.size(); ++i)
                if (segOk(i)) {
                    Vec3 lo, hi;
                    box(i, lo, hi);
                    h.insert(lo, hi);
                    segId.push_back(i);
                }
            cand.clear();
            const Vec3 pad{thr, thr, thr};
            const int64_t qChunk = 1024;
            std::vector<std::vector<std::pair<size_t, size_t>>> candPart(static_cast<size_t>((int64_t(samples.size()) + qChunk - 1) / qChunk));
            parallelFor(static_cast<int64_t>(samples.size()), qChunk, [&](int64_t b, int64_t e, int64_t ci) {
                for (size_t i = static_cast<size_t>(b); i < static_cast<size_t>(e); ++i) {
                    if (!segOk(i) || !active[i]) continue;
                    Vec3 lo, hi;
                    box(i, lo, hi);
                    h.query(lo - pad, hi + pad, [&](int id) {
                        const size_t j = segId[id];
                        if (j == i || (active[j] && j < i)) return;  // each active pair once
                        if (samples[j].y == samples[i].y && std::fabs(samples[j].s - samples[i].s) < mp.window) return;
                        candPart[ci].push_back({i, j});
                    });
                }
            }, mp.exec);
            for (const auto& part : candPart) cand.insert(cand.end(), part.begin(), part.end());
            if (bind && bind->surface)
                for (int y = 0; y < Y; ++y)
                    parallelFor(static_cast<int64_t>(P[y].size()), 256, [&](int64_t b, int64_t e, int64_t) {
                        for (size_t c = static_cast<size_t>(b); c < static_cast<size_t>(e); ++c) {
                            if (!freeCtrl[y][c]) continue;
                            ProjState next;
                            const ProjResult r = bind->surface->continueFrom(bind->state[y][c], P[y][c], next, bind->cp);
                            bind->result[y][c] = r;
                            if (r != ProjResult::Owned) continue;  // no band for an unresolved control
                            bind->state[y][c] = next;
                            double z = dot(P[y][c] - next.foot, next.n);
                            if (std::fabs(z) > mp.bandHalf) P[y][c] -= next.n * (z - (z > 0 ? mp.bandHalf : -mp.bandHalf));
                        }
                    }, mp.exec);
            else if (surface)
                for (int y = 0; y < Y; ++y)
                    parallelFor(static_cast<int64_t>(P[y].size()), 256, [&](int64_t b, int64_t e, int64_t) {
                        for (size_t c = static_cast<size_t>(b); c < static_cast<size_t>(e); ++c) {
                            if (!freeCtrl[y][c]) continue;
                            SurfaceSample q = surface->project(P[y][c]);
                            double z = dot(P[y][c] - q.p, q.n);
                            if (std::fabs(z) > mp.bandHalf) P[y][c] -= q.n * (z - (z > 0 ? mp.bandHalf : -mp.bandHalf));
                        }
                    }, mp.exec);
        }
        // Bend limiter per span (sampled curvature, 5 points).
        for (int y = 0; y < Y; ++y) {
            const int nSp = static_cast<int>(P[y].size()) - 3;
            for (int k = 0; k < nSp; ++k) {
                if (!freeCtrl[y][k + 1] && !freeCtrl[y][k + 2]) continue;
                CubicBezier bz = BSpline{{P[y][k], P[y][k + 1], P[y][k + 2], P[y][k + 3]}}.span(0);
                double kap = 0;
                for (int q = 0; q <= 4; ++q) {
                    Vec3 d1 = bz.d1(q / 4.0), d2 = bz.d2(q / 4.0);
                    double l1 = norm(d1);
                    if (l1 > 1e-12) kap = std::max(kap, norm(cross(d1, d2)) / (l1 * l1 * l1));
                }
                if (kap <= mp.kappaMax) continue;
                const double lam = 0.5 * (1 - mp.kappaMax / kap);
                for (int c : {k + 1, k + 2}) {
                    if (!freeCtrl[y][c] || c <= 0 || c + 1 >= static_cast<int>(P[y].size())) continue;
                    P[y][c] += ((P[y][c - 1] + P[y][c + 1]) * 0.5 - P[y][c]) * lam;
                }
            }
        }
        // Edge lengths (soft).
        for (int y = 0; y < Y; ++y)
            for (size_t c = 0; c + 1 < P[y].size(); ++c) {
                const bool fa = freeCtrl[y][c], fb = freeCtrl[y][c + 1];
                if (!fa && !fb) continue;
                Vec3 d = P[y][c + 1] - P[y][c];
                double l = norm(d);
                if (l <= 1e-12) continue;
                Vec3 corr = d * ((l - rest[y][c]) / l * mp.lengthStiffness);
                if (fa && fb) { P[y][c] += corr * 0.5; P[y][c + 1] -= corr * 0.5; }
                else if (fa) P[y][c] += corr;
                else P[y][c + 1] -= corr;
            }
        // Contacts.
        for (const auto& [i, j] : cand) {
            Vec3 a0 = evalSample(P, samples[i]), a1 = evalSample(P, samples[i + 1]);
            Vec3 b0 = evalSample(P, samples[j]), b1 = evalSample(P, samples[j + 1]);
            const ClosestSegmentPoints closestPair = closestSegmentPoints(a0, a1, b0, b1);
            const double d = closestPair.distance, s = closestPair.s, t = closestPair.t;
            if (d >= thr) continue;
            Vec3 pa = lerp(a0, a1, s), pb = lerp(b0, b1, t);
            Vec3 n = d > 1e-9 ? (pa - pb) / d : Vec3{0, 0, 1};
            const double corr = thr - d;
            const bool ma = active[i], mb = active[j];
            const double fa = ma && mb ? 0.5 : ma ? 1.0 : 0.0, fb = ma && mb ? 0.5 : mb ? 1.0 : 0.0;
            if (fa > 0) {
                applyToSample(P, freeCtrl, samples[i], n * (corr * fa * (1 - s)));
                applyToSample(P, freeCtrl, samples[i + 1], n * (corr * fa * s));
            }
            if (fb > 0) {
                applyToSample(P, freeCtrl, samples[j], n * (-corr * fb * (1 - t)));
                applyToSample(P, freeCtrl, samples[j + 1], n * (-corr * fb * t));
            }
        }
        // Displacement cap (keeps the straight-motion certificate's delta bounded).
        if (mp.maxDisplacement > 0)
            for (int y = 0; y < Y; ++y)
                for (size_t c = 0; c < P[y].size(); ++c) {
                    if (!freeCtrl[y][c]) continue;
                    Vec3 d = P[y][c] - Pstart[y][c];
                    double l = norm(d);
                    if (l > mp.maxDisplacement) P[y][c] = Pstart[y][c] + d * (mp.maxDisplacement / l);
                }
    }
}

void storeCpu(CurveSet& cs, const Ctrl& P, const ExecPolicy& exec = {}) {
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        cs.yarns[y].curve.ctrl = P[y];
        cs.yarns[y].finalize(1e-9, exec);
    }
}

} // namespace

BoundaryCheck checkBoundary(const CurveSet& cs, const SurfaceBinding& bind, const MeshSurface& surf, double Rs, const BoundaryCheckSteps& steps) {
    const double step = steps.step, minStep = steps.minStep;
    const ExecPolicy& exec = steps.exec;
    BoundaryCheck out;
    const TriMesh& m = surf.mesh();
    std::vector<std::array<Vec3, 2>> segs;
    for (size_t t = 0; t < m.tri.size(); ++t)
        for (int e = 0; e < 3; ++e)
            if (surf.isBoundaryEdge(static_cast<int>(t), e)) segs.push_back({m.pos[m.tri[t][e]], m.pos[m.tri[t][(e + 1) % 3]]});
    const double reach = Rs + 4 * step;  // query half-size: beyond it the distance is at least reach
    SpatialHash3 hash(std::max(reach, 1e-6));
    for (const auto& s : segs) {
        const Vec3 lo{std::min(s[0].x, s[1].x), std::min(s[0].y, s[1].y), std::min(s[0].z, s[1].z)};
        const Vec3 hi{std::max(s[0].x, s[1].x), std::max(s[0].y, s[1].y), std::max(s[0].z, s[1].z)};
        hash.insert(lo, hi);
    }
    auto dist = [&](const Vec3& x) {
        double d = reach;
        const Vec3 r{reach, reach, reach};
        hash.query(x - r, x + r, [&](int id) {
            const Vec3 a = segs[id][0], ab = segs[id][1] - segs[id][0];
            const double L2 = dot(ab, ab);
            const double t = L2 > 0 ? std::clamp(dot(x - a, ab) / L2, 0.0, 1.0) : 0.0;
            d = std::min(d, norm(x - (a + ab * t)));
        });
        return d;
    };
    for (size_t y = 0; y < cs.yarns.size(); ++y) {
        cancelCheckpoint();  // UX P3: per yarn (B1 / B2 samples)
        const YarnPath& yp = cs.yarns[y];
        // B1: arc-length samples; an interval [u0, u1] with enclosed arc length L <= hi is
        // certified if (d0 + d1 - hi) / 2 >= Rs (1-Lipschitz distance along the curve).
        std::vector<double> us;
        const double L = yp.length();
        const int n = std::max(1, static_cast<int>(std::ceil(L / step)));
        for (int k = 0; k <= n; ++k) us.push_back(yp.paramAtArcLength(L * k / n));
        // Per sample, then per initial interval (its refinement is independent of the others),
        // then per B2 sample, each in parallel into its own counters and minimum: sums and the
        // minimum do not depend on the order, so the result is the serial one.
        struct Acc {
            double minDistance = 1e300;
            int fail = 0, unresolved = 0;
            long b2 = 0;
            int b2Unresolved[5] = {0, 0, 0, 0, 0};
        };
        const int64_t chunk = exec.chunk > 0 ? exec.chunk : 64;
        auto blocks = [&](int64_t count) { return std::vector<Acc>(static_cast<size_t>((count + chunk - 1) / chunk)); };
        auto merge = [&](const std::vector<Acc>& v) {
            for (const Acc& a : v) {
                out.minDistance = std::min(out.minDistance, a.minDistance);
                out.b1Fail += a.fail;
                out.b1Unresolved += a.unresolved;
                out.b2Samples += a.b2;
                for (int k = 0; k < 5; ++k) out.b2Unresolved[k] += a.b2Unresolved[k];
            }
        };
        std::vector<Acc> sampleAcc = blocks(n + 1);
        parallelFor(n + 1, chunk, [&](int64_t b, int64_t e, int64_t ci) {
            for (int64_t k = b; k < e; ++k) {
                const double d = dist(yp.positionAtParam(us[k]));
                sampleAcc[ci].minDistance = std::min(sampleAcc[ci].minDistance, d);
                if (d < Rs) ++sampleAcc[ci].fail;
            }
        }, exec);
        merge(sampleAcc);
        std::vector<Acc> intervalAcc = blocks(n);
        parallelFor(n, chunk, [&](int64_t b, int64_t e, int64_t ci) {
            Acc& a = intervalAcc[ci];
            std::vector<std::pair<double, double>> work;
            for (int64_t k = b; k < e; ++k) {
                work.push_back({us[k], us[k + 1]});
                while (!work.empty()) {
                    const auto [u0, u1] = work.back();
                    work.pop_back();
                    const double d0 = dist(yp.positionAtParam(u0)), d1 = dist(yp.positionAtParam(u1));
                    if (d0 < Rs || d1 < Rs) continue;  // counted as FAIL samples above (or below, when refined)
                    const double hi = yp.lengthBounds(u0, u1).hi;
                    if (0.5 * (d0 + d1 - hi) >= Rs) continue;
                    if (hi < minStep) {
                        ++a.unresolved;
                        continue;
                    }
                    const double um = 0.5 * (u0 + u1), dm = dist(yp.positionAtParam(um));
                    a.minDistance = std::min(a.minDistance, dm);
                    if (dm < Rs) {
                        ++a.fail;
                        continue;
                    }
                    work.push_back({u0, um});
                    work.push_back({um, u1});
                }
            }
        }, exec);
        merge(intervalAcc);
        // B2 (sampled): continued projection of each sample from the nearest control's state.
        if (!bind.surface || y >= bind.state.size()) continue;
        std::vector<Acc> b2Acc = blocks(n + 1);
        parallelFor(n + 1, chunk, [&](int64_t b, int64_t e, int64_t ci) {
            for (int64_t k = b; k < e; ++k) {
                const double u = us[k];
                const int c = std::clamp(static_cast<int>(std::lround(u)) + 1, 0, static_cast<int>(bind.state[y].size()) - 1);
                ProjState next;
                ContinuationParams cp = bind.cp;
                const ProjResult r = surf.continueFrom(bind.state[y][c], yp.positionAtParam(u), next, cp);
                ++b2Acc[ci].b2;
                if (r != ProjResult::Owned) ++b2Acc[ci].b2Unresolved[static_cast<int>(r)];
            }
        }, exec);
        merge(b2Acc);
    }
    out.b1 = out.b1Fail ? Status::Fail : out.b1Unresolved ? Status::Unresolved : Status::Pass;
    const int overhang = out.b2Unresolved[static_cast<int>(ProjResult::OpenBoundary)];
    int other = 0;
    for (int k = 1; k < 5; ++k) other += out.b2Unresolved[k];
    out.b2 = overhang ? Status::Fail : other ? Status::Unresolved : Status::Pass;
    char b[320];
    std::snprintf(b, sizeof b,
                  "B1 (certified) %s: %zu open-boundary edges, min centreline distance %.4f mm vs R_s %.4f (%d FAIL, %d UNRESOLVED); "
                  "B2 (sampled, step %.3f mm) %s: %ld samples, %d overhang, %d budget, %d normal-turn, %d depth-band",
                  toString(out.b1), segs.size(), out.minDistance, Rs, out.b1Fail, out.b1Unresolved, step, toString(out.b2),
                  out.b2Samples, overhang, out.b2Unresolved[1], out.b2Unresolved[2], out.b2Unresolved[4]);
    out.detail = b;
    return out;
}

ValidationReport CpuNumerics::validate(const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt) {
    return knit::validate(cs, pp, opt);
}

CurvatureResult CpuNumerics::curvature(const YarnPath& y, double radius, const Allowances& al, const ExecPolicy& exec) {
    return checkCurvature(y, radius, al, exec);
}

void CpuNumerics::relax(const RelaxRequest& r) { relaxFree(r.P, r.Prest, r.freeCtrl, r.cs, r.surface, r.mp, r.bind); }

void CpuNumerics::store(CurveSet& cs, const Controls& P, const ExecPolicy& exec) { storeCpu(cs, P, exec); }

std::unique_ptr<MeshSurface> CpuNumerics::meshSurface(const TriMesh& mesh) { return std::make_unique<MeshSurface>(mesh); }

std::unique_ptr<BaseSurface> CpuNumerics::pieceQuery(const Chart& chart) { return std::make_unique<PieceSurface>(chart); }

// Each node on the chart: exactly one triangle, outside the solids, and (with the footprint test)
// clear of the chart's border. The run-end tests stay -1: tileMesh evaluates them when needed.
void CpuNumerics::classifyLattice(const UvAtlas& atlas, const LatticeChartQuery& q, const LatticeSolids& solids,
                                  std::vector<LatticeNodeClass>& out) {
    out.assign(static_cast<size_t>(q.ni) * q.nr, LatticeNodeClass{});
    const TriMesh& m = atlas.mesh();
    for (int r = 0; r < q.nr; ++r) {
        cancelCheckpoint();
        for (int i = 0; i < q.ni; ++i) {
            LatticeNodeClass& c = out[static_cast<size_t>(r) * q.ni + i];
            const double x = q.xmin + i * q.w, y = q.ymin + r * q.h;
            const std::vector<int> hits = atlas.locate(q.chart, q.c2 * x + q.v2 * y, c.bary);
            if (hits.empty()) continue;  // code 1: outside
            if (hits.size() > 1) {
                c.code = 2;
                continue;
            }
            c.tri = hits[0];
            const auto& tv = m.tri[c.tri];
            c.p = m.pos[tv[0]] * c.bary[0] + m.pos[tv[1]] * c.bary[1] + m.pos[tv[2]] * c.bary[2];
            if (solids.contains(c.p)) c.code = 3;
            else if (q.footprintCheck && !footprintClear(atlas, q, {x, y, 0.5 * q.w, 0.5 * q.w})) c.code = 4;
            else c.code = 0;
        }
    }
}

SurfaceCheck CpuNumerics::checkRenderedSurface(const CurveSet& cs, const PatternParams& ppR, const SurfaceParams& sp) {
    const SurfaceModel sm = buildSurface(cs, ppR, sp, exec_);
    SurfaceCheck check = checkSurface(cs, sm, {}, exec_);
    check.compressionClamped = sm.compression.clamped;
    return check;
}

BoundaryCheck CpuNumerics::checkBoundary(const CurveSet& cs, const SurfaceBinding& bind, const MeshSurface& surf, double Rs) {
    BoundaryCheckSteps steps;
    steps.exec = exec_;
    return knit_test::checkBoundary(cs, bind, surf, Rs, steps);
}

void SabotagedRound::relax(const RelaxRequest& r) {
    const bool firstRound = r.mp.call.round == 0;
    if (firstRound)
        for (int c = c0_; c < c1_; ++c) {
            controlsBefore.push_back(r.cs.yarns[yarn_].curve.ctrl[c]);
            statesBefore.push_back(r.bind->state[yarn_][c]);
            resultsBefore.push_back(r.bind->result[yarn_][c]);
        }
    inner_.relax(r);
    if (!firstRound) return;
    for (int c = c0_; c < c1_; ++c) {
        Vec3& q = r.P[yarn_][c];
        const ProjState& s = r.bind->state[yarn_][c];
        q = q - s.n * (2 * dot(q - s.foot, s.n));
        r.bind->state[yarn_][c].tri = -7;
    }
    relaxed = r.P;
    relaxedBinding = *r.bind;
}

void RejectEveryAttempt::relax(const RelaxRequest& r) {
    ++attempts;
    inner_.relax(r);
    rejectNext_ = true;
}

ValidationReport RejectEveryAttempt::validate(const CurveSet& cs, const PatternParams& pp, const ValidateOptions& opt) {
    ValidationReport report = inner_.validate(cs, pp, opt);
    if (rejectNext_) report.verdict = Status::Fail;
    rejectNext_ = false;
    return report;
}

knit::KnitNumerics& cpuNumerics() {
    static CpuNumerics n;
    return n;
}

} // namespace knit_test
