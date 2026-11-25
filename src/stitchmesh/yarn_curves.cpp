// The embedded yarn as production curves and back, and the threading counts (yarn_curves.h).
#include "stitchmesh/yarn_curves.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "knitting/tube/plies.h"
#include "model/mesh_query.h"

namespace knit {

namespace {

// Distance from q to the segment [a, b].
double pointSegment(const Vec3& q, const Vec3& a, const Vec3& b) {
    const Vec3 d = b - a;
    const double dd = dot(d, d);
    const double t = dd > 0 ? std::clamp(dot(q - a, d) / dd, 0.0, 1.0) : 0.0;
    return norm(q - (a + d * t));
}

// The signed offset of x along the reference normal at its closest point, in units of D.
double depthOf(const MeshQuery& q, const TriMesh& ref, const Vec3& x, double D) {
    const MeshQuery::Hit h = q.closest(x);
    if (h.tri < 0) return 0;
    const auto& t = ref.tri[h.tri];
    const TrianglePoint tp = closestOnTriangle(x, ref.pos[t[0]], ref.pos[t[1]], ref.pos[t[2]]);
    Vec3 n{0, 0, 0};
    for (int k = 0; k < 3; ++k) n = n + ref.vnormal[t[k]] * tp.bary[k];
    if (norm(n) == 0) n = ref.faceNormal(h.tri);
    return dot(x - tp.p, normalized(n)) / D;
}

} // namespace

std::string YarnConversion::text() const {
    char b[400];
    std::snprintf(b, sizeof b,
                  "conversion: %zu strand(s) as uniform cubic B-splines, spacing %.4f mm (%zu controls); length Hermite %.2f mm, B-spline %.2f mm (%+.4f %%); "
                  "deviation of %d dense Hermite points: nearest max %.4f mm, p99 %.4f mm; at equal parameter max %.4f mm\n",
                  curves.yarns.size(), spacing.empty() ? 0.0 : spacing[0], curves.yarns.empty() ? size_t{0} : curves.yarns[0].curve.ctrl.size(), hermiteLength, curveLength,
                  hermiteLength > 0 ? 100.0 * (curveLength / hermiteLength - 1) : 0.0, densePoints, maxDeviation, p99Deviation, maxEqualParameter);
    return b;
}

bool yarnToCurves(const YarnEmbedding& e, const YarnSettings& st, double spacing, int densePerSpan, YarnConversion& out, std::string* err) {
    out = YarnConversion{};
    out.curves.kind = CurveKind::Surface;
    out.strandOfPiece.assign(e.pieces.size(), -1);
    out.u0.assign(e.pieces.size(), 0);
    out.u1.assign(e.pieces.size(), 0);
    out.knotU.assign(e.pieces.size(), {});
    std::vector<double> deviations;
    for (size_t si = 0; si < e.strands.size(); ++si) {
        const YarnStrand& sd = e.strands[si];
        if (sd.closed) {
            if (err) *err = "strand " + std::to_string(si) + " is closed; production curves are open (closed strands are not converted)";
            return false;
        }
        // the dense polyline of the strand, each piece's first and last dense index, each knot's dense index
        std::vector<Vec3> X;
        std::vector<int> first(sd.pieces.size()), last(sd.pieces.size());
        std::vector<std::vector<int>> knotIndex(sd.pieces.size());
        for (size_t k = 0; k < sd.pieces.size(); ++k) {
            const YarnPiece& p = e.pieces[sd.pieces[k]];
            const std::vector<Vec3> d = sampleYarnHermite(p.controls, p.tangents, densePerSpan, st.spanTangents);
            const int base = X.empty() ? 0 : static_cast<int>(X.size()) - 1;   // a shared port once
            first[k] = base;
            for (size_t m = X.empty() ? 0 : 1; m < d.size(); ++m) X.push_back(d[m]);
            last[k] = static_cast<int>(X.size()) - 1;
            for (size_t j = 0; j < p.controls.size(); ++j) knotIndex[k].push_back(base + static_cast<int>(j) * densePerSpan);
        }
        std::vector<double> s(X.size(), 0);
        for (size_t k = 1; k < X.size(); ++k) s[k] = s[k - 1] + norm(X[k] - X[k - 1]);
        const double L = s.back();
        const int n = std::max(3, static_cast<int>(std::lround(L / spacing)));
        const double h = L / n;
        // n + 1 points at equal arc length on the dense polyline
        std::vector<Vec3> P;
        P.reserve(n + 1);
        size_t seg = 0;
        for (int i = 0; i <= n; ++i) {
            const double target = i == n ? L : i * h;
            while (seg + 2 < X.size() && s[seg + 1] < target) ++seg;
            const double len = s[seg + 1] - s[seg];
            const double t = len > 0 ? std::clamp((target - s[seg]) / len, 0.0, 1.0) : 0.0;
            P.push_back(X[seg] + (X[seg + 1] - X[seg]) * t);
        }
        YarnPath y;
        y.id = static_cast<int>(si);
        y.name = "strand " + std::to_string(si);
        y.curve = interpolateBSpline(P);
        y.finalize(1e-6);
        out.hermiteLength += L;
        out.curveLength += y.length();
        out.spacing.push_back(h);
        for (size_t k = 0; k < sd.pieces.size(); ++k) {
            const int pi = sd.pieces[k];
            out.strandOfPiece[pi] = static_cast<int>(si);
            out.u0[pi] = s[first[k]] / h;
            out.u1[pi] = s[last[k]] / h;
            for (int idx : knotIndex[k]) out.knotU[pi].push_back(s[idx] / h);
        }
        // the approximation: each dense point against the curve near its own parameter
        const int per = 16;
        for (size_t k = 0; k < X.size(); ++k) {
            const double u = std::min(s[k] / h, y.curve.uMax());
            out.maxEqualParameter = std::max(out.maxEqualParameter, norm(y.curve.eval(u) - X[k]));
            const int span = std::min(static_cast<int>(u), y.curve.numSpans() - 1);
            double best = 1e300;
            for (int sp = std::max(0, span - 2); sp <= std::min(y.curve.numSpans() - 1, span + 2); ++sp) {
                Vec3 a = y.curve.eval(sp);
                for (int q = 1; q <= per; ++q) {
                    const Vec3 b2 = y.curve.eval(sp + static_cast<double>(q) / per);
                    best = std::min(best, pointSegment(X[k], a, b2));
                    a = b2;
                }
            }
            deviations.push_back(best);
            out.maxDeviation = std::max(out.maxDeviation, best);
        }
        out.densePoints += static_cast<int>(X.size());
        out.curves.yarns.push_back(std::move(y));
    }
    if (!deviations.empty()) {
        std::sort(deviations.begin(), deviations.end());
        out.p99Deviation = deviations[static_cast<size_t>(0.99 * (deviations.size() - 1))];
    }
    return true;
}

YarnEmbedding embeddingFromCurves(const YarnEmbedding& e, const YarnConversion& cv, const CurveSet& curves, const TriMesh& reference, double D) {
    YarnEmbedding out = e;
    TriMesh ref = reference;
    if (ref.vnormal.size() != ref.pos.size()) ref.computeNormals();
    const MeshQuery q(ref);
    for (size_t i = 0; i < out.pieces.size(); ++i) {
        YarnPiece& p = out.pieces[i];
        const int y = cv.strandOfPiece[i];
        if (y < 0) continue;
        const BSpline& c = curves.yarns[y].curve;
        const size_t ns = p.samples.size();
        for (size_t k = 0; k < ns; ++k) {
            const double t = ns > 1 ? static_cast<double>(k) / (ns - 1) : 0.0;
            const double u = k + 1 == ns ? cv.u1[i] : cv.u0[i] + (cv.u1[i] - cv.u0[i]) * t;
            p.samples[k] = c.eval(u);
        }
        for (size_t j = 0; j < p.controls.size(); ++j) {
            const double u = j == 0 ? cv.u0[i] : j + 1 == p.controls.size() ? cv.u1[i] : cv.knotU[i][j];
            p.controls[j] = c.eval(u);
            p.tangents[j] = c.d1(u);
            p.depth[j] = depthOf(q, ref, p.controls[j], D);
        }
    }
    for (YarnStrand& sd : out.strands) {
        sd.samples.clear();
        sd.sampleFace.clear();
        sd.samplePiece.clear();
        for (int i : sd.pieces) {
            const YarnPiece& p = out.pieces[i];
            for (size_t k = sd.samples.empty() ? 0 : 1; k < p.samples.size(); ++k) {
                sd.samples.push_back(p.samples[k]);
                sd.sampleFace.push_back(p.face);
                sd.samplePiece.push_back(i);
            }
        }
        sd.arc.assign(sd.samples.size(), 0);
        for (size_t k = 1; k < sd.samples.size(); ++k) sd.arc[k] = sd.arc[k - 1] + norm(sd.samples[k] - sd.samples[k - 1]);
    }
    return out;
}

} // namespace knit
