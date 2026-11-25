#include "stitchmesh/coverage.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>

#include "model/mesh_query.h"
#include "stitchmesh/typed_mesh.h"

namespace knit {

namespace {

// A uniform grid over triangle boxes, for the normal-segment test.
struct TriangleGrid {
    Vec3 lo, hi;
    double cell = 1;
    int nx = 1, ny = 1, nz = 1;
    std::vector<std::vector<int>> cells;

    int index(int x, int y, int z) const { return (z * ny + y) * nx + x; }
    void clampCell(const Vec3& p, int& x, int& y, int& z) const {
        x = std::clamp(static_cast<int>((p.x - lo.x) / cell), 0, nx - 1);
        y = std::clamp(static_cast<int>((p.y - lo.y) / cell), 0, ny - 1);
        z = std::clamp(static_cast<int>((p.z - lo.z) / cell), 0, nz - 1);
    }
};

TriangleGrid buildGrid(const TriMesh& m, double cell) {
    TriangleGrid g;
    m.bounds(g.lo, g.hi);
    g.lo = g.lo - Vec3{cell, cell, cell};
    g.hi = g.hi + Vec3{cell, cell, cell};
    g.cell = cell;
    g.nx = std::max(1, static_cast<int>((g.hi.x - g.lo.x) / cell) + 1);
    g.ny = std::max(1, static_cast<int>((g.hi.y - g.lo.y) / cell) + 1);
    g.nz = std::max(1, static_cast<int>((g.hi.z - g.lo.z) / cell) + 1);
    g.cells.assign(static_cast<size_t>(g.nx) * g.ny * g.nz, {});
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t) {
        Vec3 a = m.pos[m.tri[t][0]], b = a;
        for (int k = 1; k < 3; ++k) {
            const Vec3& p = m.pos[m.tri[t][k]];
            a = Vec3{std::min(a.x, p.x), std::min(a.y, p.y), std::min(a.z, p.z)};
            b = Vec3{std::max(b.x, p.x), std::max(b.y, p.y), std::max(b.z, p.z)};
        }
        int x0, y0, z0, x1, y1, z1;
        g.clampCell(a, x0, y0, z0);
        g.clampCell(b, x1, y1, z1);
        for (int z = z0; z <= z1; ++z)
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) g.cells[g.index(x, y, z)].push_back(t);
    }
    return g;
}

// Moller-Trumbore on the segment p + s d, s in [-1, 1]: true with s when it crosses triangle t (edges included).
bool segmentHits(const TriMesh& m, int t, const Vec3& p, const Vec3& d, double& s) {
    const Vec3 a = m.pos[m.tri[t][0]], b = m.pos[m.tri[t][1]], c = m.pos[m.tri[t][2]];
    const Vec3 e1 = b - a, e2 = c - a;
    const Vec3 h = cross(d, e2);
    const double det = dot(e1, h);
    if (std::fabs(det) < 1e-14) return false;
    const double inv = 1.0 / det;
    const Vec3 q = p - a;
    const double u = dot(q, h) * inv;
    if (u < -1e-9 || u > 1 + 1e-9) return false;
    const Vec3 r = cross(q, e1);
    const double v = dot(d, r) * inv;
    if (v < -1e-9 || u + v > 1 + 1e-9) return false;
    s = dot(e2, r) * inv;
    return s >= -1 && s <= 1;
}

double pointPolyline(const Vec3& p, const std::vector<Vec3>& loop) {
    double best = 1e300;
    for (size_t i = 0; i < loop.size(); ++i) best = std::min(best, closestSegmentPoints(p, p, loop[i], loop[(i + 1) % loop.size()]).distance);
    return best;
}

// Points every `step` along a closed polyline (its corners included).
std::vector<Vec3> densify(const std::vector<Vec3>& loop, double step) {
    std::vector<Vec3> out;
    for (size_t i = 0; i < loop.size(); ++i) {
        const Vec3 a = loop[i], b = loop[(i + 1) % loop.size()];
        const int n = std::max(1, static_cast<int>(std::ceil(norm(b - a) / step)));
        for (int j = 0; j < n; ++j) out.push_back(a + (b - a) * (static_cast<double>(j) / n));
    }
    return out;
}

double directedDistance(const std::vector<Vec3>& from, const std::vector<Vec3>& to) {
    double worst = 0;
    for (const Vec3& p : from) worst = std::max(worst, pointPolyline(p, to));
    return worst;
}

} // namespace

std::vector<std::vector<int>> polygonBorderLoops(const PolygonMesh& m) {
    std::map<std::pair<int, int>, int> count;
    for (int f = 0; f < m.faces(); ++f)
        for (int k = 0; k < m.sides(f); ++k) {
            const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % m.sides(f));
            count[{std::min(a, b), std::max(a, b)}] += 1;
        }
    std::map<int, int> nextOf;
    for (int f = 0; f < m.faces(); ++f)
        for (int k = 0; k < m.sides(f); ++k) {
            const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % m.sides(f));
            if (count[{std::min(a, b), std::max(a, b)}] == 1) nextOf[a] = b;
        }
    std::vector<std::vector<int>> loops;
    std::set<int> seen;
    for (const auto& [v0, unused] : nextOf) {
        (void) unused;
        if (seen.count(v0)) continue;
        std::vector<int> loop;
        for (int v = v0; !seen.count(v) && nextOf.count(v); v = nextOf[v]) {
            seen.insert(v);
            loop.push_back(v);
        }
        loops.push_back(loop);
    }
    return loops;
}

CoverageReport measureCoverage(const TriMesh& reference, const std::vector<std::vector<int>>& refLoops, const PolygonMesh& out, const CoverageOptions& o) {
    CoverageReport rep;
    // the output surface fanned
    TriMesh fan;
    fan.pos = out.pos;
    std::vector<int> faceOfTri;
    for (int f = 0; f < out.faces(); ++f) {
        std::vector<std::array<int, 3>> tris;
        earClipFace(out, f, tris);
        for (const std::array<int, 3>& t : tris) {
            fan.tri.push_back(t);
            faceOfTri.push_back(f);
        }
    }
    if (fan.tri.empty()) {
        rep.status = Status::Fail;
        rep.detail = "no output faces";
        return rep;
    }
    fan.computeNormals();
    const MeshQuery outQuery(fan);
    const MeshQuery refQuery(reference);
    double meanEdge = 0;
    for (const auto& t : fan.tri) meanEdge += norm(fan.pos[t[1]] - fan.pos[t[0]]);
    meanEdge /= fan.tri.size();
    const TriangleGrid grid = buildGrid(fan, std::max(2 * meanEdge, 4 * o.coverDistance));
    // reference loops as polylines
    std::vector<std::vector<Vec3>> refPolys;
    for (const std::vector<int>& l : refLoops) {
        std::vector<Vec3> p;
        for (int v : l) p.push_back(reference.pos[v]);
        refPolys.push_back(p);
    }
    // samples
    std::vector<double> vertexArea(reference.pos.size(), 0);
    for (int t = 0; t < static_cast<int>(reference.tri.size()); ++t) {
        const auto& tr = reference.tri[t];
        const double a = 0.5 * norm(cross(reference.pos[tr[1]] - reference.pos[tr[0]], reference.pos[tr[2]] - reference.pos[tr[0]]));
        rep.referenceArea += a;
        CoverageSample s;
        s.p = (reference.pos[tr[0]] + reference.pos[tr[1]] + reference.pos[tr[2]]) / 3.0;
        s.weight = 0.5 * a;
        rep.perSample.push_back(s);
        for (int k = 0; k < 3; ++k) vertexArea[tr[k]] += a / 3.0;
    }
    const size_t centroids = rep.perSample.size();
    for (size_t v = 0; v < reference.pos.size(); ++v) {
        CoverageSample s;
        s.p = reference.pos[v];
        s.weight = 0.5 * vertexArea[v];
        s.vertex = true;
        rep.perSample.push_back(s);
    }
    std::vector<int> stamp(fan.tri.size(), -1);
    for (size_t i = 0; i < rep.perSample.size(); ++i) {
        CoverageSample& s = rep.perSample[i];
        const Vec3 n = s.vertex ? reference.vnormal[i - centroids] : reference.faceNormal(static_cast<int>(i));
        for (const std::vector<Vec3>& poly : refPolys)
            if (pointPolyline(s.p, poly) < o.boundaryBand) s.excluded = true;
        const MeshQuery::Hit hit = outQuery.closest(s.p, 4 * o.coverDistance);
        s.distance = hit.tri >= 0 ? hit.dist : 1e300;
        s.face = hit.tri >= 0 && hit.dist <= o.coverDistance ? faceOfTri[hit.tri] : -1;
        // crossings of the normal segment
        const Vec3 d = n * o.coverDistance;
        int x0, y0, z0, x1, y1, z1;
        const Vec3 a = s.p - d, b = s.p + d;
        grid.clampCell(Vec3{std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}, x0, y0, z0);
        grid.clampCell(Vec3{std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}, x1, y1, z1);
        std::vector<double> hits;
        for (int z = z0; z <= z1; ++z)
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x)
                    for (int t : grid.cells[grid.index(x, y, z)]) {
                        if (stamp[t] == static_cast<int>(i)) continue;
                        stamp[t] = static_cast<int>(i);
                        double par = 0;
                        if (segmentHits(fan, t, s.p, d, par)) hits.push_back(par);
                    }
        std::sort(hits.begin(), hits.end());
        int distinct = 0;
        for (size_t h = 0; h < hits.size(); ++h)
            if (h == 0 || hits[h] - hits[h - 1] > 1e-6) ++distinct;
        s.crossings = distinct;
        ++rep.samples;
        if (s.excluded) {
            ++rep.excluded;
            continue;
        }
        rep.consideredArea += s.weight;
        if (s.face < 0) {
            ++rep.uncovered;
            rep.uncoveredArea += s.weight;
        }
        if (s.crossings > 1) ++rep.coveredTwice;
    }
    rep.uncoveredMeasured = rep.consideredArea > 0;
    rep.uncoveredFraction = rep.uncoveredMeasured ? rep.uncoveredArea / rep.consideredArea : 0;
    // areas
    for (int t = 0; t < static_cast<int>(fan.tri.size()); ++t) {
        const auto& tr = fan.tri[t];
        const Vec3 c = cross(fan.pos[tr[1]] - fan.pos[tr[0]], fan.pos[tr[2]] - fan.pos[tr[0]]);
        const double a = 0.5 * norm(c);
        rep.outputArea += a;
        const MeshQuery::Hit hit = refQuery.closest((fan.pos[tr[0]] + fan.pos[tr[1]] + fan.pos[tr[2]]) / 3.0);
        if (hit.tri >= 0 && a > 0) rep.projectedArea += a * std::fabs(dot(c / (2 * a), reference.faceNormal(hit.tri)));
    }
    rep.projectedRatio = rep.referenceArea > 0 ? rep.projectedArea / rep.referenceArea : 0;
    // boundary correspondence
    const std::vector<std::vector<int>> outLoops = polygonBorderLoops(out);
    rep.outputLoops = static_cast<int>(outLoops.size());
    rep.referenceLoops = static_cast<int>(refLoops.size());
    std::vector<std::vector<Vec3>> refDense;
    for (const std::vector<Vec3>& p : refPolys) refDense.push_back(densify(p, 0.25));
    std::set<int> matchedRef;
    for (int i = 0; i < rep.outputLoops; ++i) {
        std::vector<Vec3> poly;
        for (int v : outLoops[i]) poly.push_back(out.pos[v]);
        const std::vector<Vec3> dense = densify(poly, 0.25);
        LoopCorrespondence lc;
        lc.outputLoop = i;
        lc.outputVertices = static_cast<int>(outLoops[i].size());
        double best = 1e300;
        for (int j = 0; j < rep.referenceLoops; ++j) {
            const double d = directedDistance(dense, refPolys[j]);
            if (d < best) {
                best = d;
                lc.referenceLoop = j;
            }
        }
        if (lc.referenceLoop >= 0) {
            lc.outputToReference = best;
            lc.referenceToOutput = directedDistance(refDense[lc.referenceLoop], poly);
            lc.hausdorff = std::max(lc.outputToReference, lc.referenceToOutput);
            matchedRef.insert(lc.referenceLoop);
        }
        rep.worstHausdorff = std::max(rep.worstHausdorff, lc.hausdorff);
        rep.loops.push_back(lc);
    }
    std::string fails;
    auto add = [&](bool bad, const std::string& what) {
        if (!bad) return;
        if (!fails.empty()) fails += "; ";
        fails += what;
    };
    char b[200];
    std::snprintf(b, sizeof b, "uncovered %.4f %% > %.2f %%", 100 * rep.uncoveredFraction, 100 * o.uncoveredLimit);
    add(rep.uncoveredMeasured && rep.uncoveredFraction > o.uncoveredLimit, b);
    add(rep.coveredTwice > 0, "samples covered twice " + std::to_string(rep.coveredTwice));
    std::snprintf(b, sizeof b, "projected area ratio %.4f outside [%.2f, %.2f]", rep.projectedRatio, o.areaLow, o.areaHigh);
    add(rep.projectedRatio < o.areaLow || rep.projectedRatio > o.areaHigh, b);
    add(rep.outputLoops != rep.referenceLoops || static_cast<int>(matchedRef.size()) != rep.referenceLoops, "boundary loops do not correspond one to one");
    std::snprintf(b, sizeof b, "boundary Hausdorff %.3f mm > %.3f", rep.worstHausdorff, o.hausdorffLimit);
    add(rep.worstHausdorff > o.hausdorffLimit, b);
    rep.status = !fails.empty() ? Status::Fail : rep.uncoveredMeasured ? Status::Pass : Status::Unresolved;
    rep.detail = fails;
    if (!rep.uncoveredMeasured)
        rep.detail += std::string(fails.empty() ? "" : "; ") + "uncovered share not measurable: every sample lies within the boundary band";
    return rep;
}

std::string CoverageReport::text() const {
    char b[900];
    std::snprintf(b, sizeof b,
                  "coverage %s: samples %d (excluded in the boundary band %d), uncovered %d = %.4f %% of %.1f mm^2 considered; covered twice %d; output area %.1f, "
                  "projected %.1f, reference %.1f mm^2, projected ratio %.5f; boundary loops output %d / reference %d, worst Hausdorff %.4f mm%s%s\n",
                  toString(status), samples, excluded, uncovered, 100 * uncoveredFraction, consideredArea, coveredTwice, outputArea, projectedArea, referenceArea,
                  projectedRatio, outputLoops, referenceLoops, worstHausdorff, detail.empty() ? "" : ": ", detail.c_str());
    return b;
}

} // namespace knit
