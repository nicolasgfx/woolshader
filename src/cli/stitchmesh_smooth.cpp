// The smooth-reference commands of stitchmesh_lab:
//
//   prep-compare      an imported piece's facets against its Loop-limit preparation: the change of shape both ways (dense
//                     samples), area, the parametric map's displacement, the levels' convergence, and (--alt-mesh) the
//                     same piece triangulated with the other quad diagonals, to show how triangulation enters Loop's surface
//   smooth-fixture    the surface treatment on a minimal case: a coarse 16-sided prism, its Loop limit at levels 1-4, a
//                     ring tube relaxed on each (and on the prism itself, and on a fine true cylinder): stop reason, the
//                     independent stationarity on the surface walked, the distance between successive levels' results
//   shaping-fixtures  minimal increase / decrease neighbourhoods (and the boundary-shaping variants) on surfaces of
//                     revolution whose ring circumferences fit the counts, relaxed with the CP3 model; each shaping cell
//                     measured with the frozen pentagon bounds and with the sub-quad assessment (typed_mesh.h), plus
//                     deliberately defective cells that the assessment must reject
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "stitchmesh_lab_shared.h"
#include "jobs/marking_session.h"
#include "model/mesh_query.h"
#include "stitchmesh/construction_fixtures.h"
#include "stitchmesh/coverage.h"
#include "stitchmesh/prepared_reference.h"
#include "stitchmesh/quantiles.h"
#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"

namespace fs = std::filesystem;
using namespace knit;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct Opts {
    std::map<std::string, std::string> kv;
    std::string get(const std::string& k, const std::string& d = "") const {
        const auto it = kv.find(k);
        return it == kv.end() ? d : it->second;
    }
    double number(const std::string& k, double d) const {
        const auto it = kv.find(k);
        return it == kv.end() ? d : std::atof(it->second.c_str());
    }
};

Opts parse(int argc, char** argv) {
    Opts o;
    for (int i = 2; i + 1 < argc; i += 2) o.kv[argv[i]] = argv[i + 1];
    return o;
}

void putf(std::ostream& o, const char* fmt, ...) {
    char b[2048];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(b, sizeof b, fmt, ap);
    va_end(ap);
    o << b;
}

void writeObj(const std::string& path, const TriMesh& m) {
    std::ofstream f(path, std::ios::binary);
    for (const Vec3& q : m.pos) putf(f, "v %.9g %.9g %.9g\n", q.x, q.y, q.z);
    for (const auto& t : m.tri) putf(f, "f %d %d %d\n", t[0] + 1, t[1] + 1, t[2] + 1);
}

// Dense samples of a mesh: its vertices, edge midpoints and triangle centroids.
std::vector<Vec3> samples(const TriMesh& m) {
    std::vector<Vec3> s = m.pos;
    std::set<std::pair<int, int>> seen;
    for (const auto& t : m.tri) {
        s.push_back((m.pos[t[0]] + m.pos[t[1]] + m.pos[t[2]]) / 3.0);
        for (int k = 0; k < 3; ++k) {
            const std::pair<int, int> e{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])};
            if (seen.insert(e).second) s.push_back(0.5 * (m.pos[e.first] + m.pos[e.second]));
        }
    }
    return s;
}

// Distances from a's samples to b's surface; signed along b's normal at the closest point when `signedOut`.
Quantiles distances(const TriMesh& a, const TriMesh& b, double* meanSigned = nullptr) {
    const MeshQuery q(b);
    std::vector<double> d;
    double sum = 0;
    for (const Vec3& x : samples(a)) {
        const MeshQuery::Hit h = q.closest(x);
        d.push_back(norm(h.p - x));
        if (h.tri >= 0) sum += dot(x - h.p, b.faceNormal(h.tri));
    }
    if (meanSigned) *meanSigned = d.empty() ? 0 : sum / d.size();
    return quantiles(d);
}

double quantile99(const TriMesh& a, const TriMesh& b) {
    const MeshQuery q(b);
    std::vector<double> d;
    for (const Vec3& x : samples(a)) d.push_back(norm(q.closest(x).p - x));
    std::sort(d.begin(), d.end());
    return d.empty() ? 0 : d[static_cast<size_t>(0.99 * (d.size() - 1))];
}

// The largest dihedral angle between neighbouring triangles [deg] and the 95th percentile.
void dihedral(const TriMesh& m, double& maxDeg, double& p95Deg) {
    std::map<std::pair<int, int>, std::vector<int>> by;
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
        for (int k = 0; k < 3; ++k) by[{std::min(m.tri[t][k], m.tri[t][(k + 1) % 3]), std::max(m.tri[t][k], m.tri[t][(k + 1) % 3])}].push_back(t);
    std::vector<double> a;
    for (const auto& [e, ts] : by)
        if (ts.size() == 2) a.push_back(std::acos(std::max(-1.0, std::min(1.0, dot(m.faceNormal(ts[0]), m.faceNormal(ts[1]))))) * 180.0 / kPi);
    const Quantiles q = quantiles(a);
    maxDeg = q.max;
    p95Deg = q.p95;
}

double signedVolume(const TriMesh& m) {
    double v = 0;
    for (const auto& t : m.tri) v += dot(m.pos[t[0]], cross(m.pos[t[1]], m.pos[t[2]])) / 6.0;
    return v;
}

int borderEdges(const TriMesh& m) {
    std::map<std::pair<int, int>, int> use;
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) use[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}] += 1;
    int n = 0;
    for (const auto& [e, c] : use) n += c == 1;
    return n;
}

// ---------------------------------------------------------------------------------------------------------- prep-compare
int commandPrepCompare(const Opts& a) {
    MarkedFile f;
    std::string err;
    if (!readMarkedFile(a.get("--mesh"), f, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const int pi = static_cast<int>(a.number("--piece", -1));
    if (pi < 0 || pi >= static_cast<int>(f.pieces->pieces.size())) {
        std::fprintf(stderr, "prep-compare: no piece %d\n", pi);
        return 1;
    }
    const double mm = a.number("--unit-mm", 1), h = a.number("--h", 3.4);
    const std::string out = a.get("--out");
    fs::create_directories(out);
    TriMesh piece = f.pieces->pieces[pi].mesh;
    piece.transform(mm, Vec3{0, 0, 0});
    piece.computeNormals();
    const int L = a.kv.count("--levels") ? static_cast<int>(a.number("--levels", 2)) : refinementLevels(piece, 0.5 * h);
    std::ofstream rep(out + "/prep_compare.txt", std::ios::binary);
    std::ofstream csv(out + "/prep_compare.csv", std::ios::binary);
    csv << "what,levels,max_mm,p95_mm,p99_mm,p50_mm,mean_signed_mm\n";
    double dmax = 0, d95 = 0;
    dihedral(piece, dmax, d95);
    putf(rep, "piece %d (%s): %zu triangles, %d border edges, area %.1f mm2, volume %.1f mm3 (closed: %s); creases of the import: max %.1f deg, p95 %.1f deg\n", pi,
         f.mesh ? "imported" : "?", piece.tri.size(), borderEdges(piece), piece.area(), signedVolume(piece), borderEdges(piece) == 0 ? "yes" : "no", dmax, d95);
    putf(rep, "the prepared reference's automatic levels (mean edge <= h / 2 = %.2f mm): %d\n", 0.5 * h, L);
    RefinedMesh prevLoop;
    for (int lv = std::max(0, L - 1); lv <= L + 1; ++lv) {
        const RefinedMesh F = refineMidpoint(piece, lv), S = refineLoop(piece, lv);
        double ms = 0;
        const Quantiles sf = distances(S.mesh, F.mesh, &ms), fs2 = distances(F.mesh, S.mesh);
        dihedral(S.mesh, dmax, d95);
        putf(csv, "loop_to_facets,%d,%.4f,%.4f,%.4f,%.4f,%.4f\n", lv, sf.max, sf.p95, quantile99(S.mesh, F.mesh), sf.p50, ms);
        putf(csv, "facets_to_loop,%d,%.4f,%.4f,%.4f,%.4f,\n", lv, fs2.max, fs2.p95, quantile99(F.mesh, S.mesh), fs2.p50);
        putf(rep, "level %d: %zu triangles; Loop limit -> facets max %.3f p99 %.3f p95 %.3f median %.3f mm (mean signed %.3f mm, + outside); facets -> Loop max %.3f p95 %.3f mm; "
                  "area %.1f mm2 (facets %.1f, %+.2f %%), volume %.1f mm3 (%+.2f %%); Loop creases max %.2f deg p95 %.2f deg\n",
             lv, S.mesh.tri.size(), sf.max, quantile99(S.mesh, F.mesh), sf.p95, sf.p50, ms, fs2.max, fs2.p95, S.mesh.area(), F.mesh.area(),
             100 * (S.mesh.area() / F.mesh.area() - 1), signedVolume(S.mesh), 100 * (signedVolume(S.mesh) / signedVolume(F.mesh) - 1), dmax, d95);
        // the parametric map: the same refined vertex (same triangle and weights) on either preparation
        double mapMax = 0;
        for (size_t v = 0; v < S.mesh.pos.size(); ++v) mapMax = std::max(mapMax, norm(S.mesh.pos[v] - F.mesh.pos[v]));
        putf(rep, "         a point of the file keeps its triangle and weights: it moves by up to %.3f mm between the preparations\n", mapMax);
        if (lv > std::max(0, L - 1)) {
            const Quantiles conv = distances(S.mesh, prevLoop.mesh), back = distances(prevLoop.mesh, S.mesh);
            putf(rep, "         level %d against level %d (the limit surface's approximation): max %.4f mm, p95 %.4f mm (back %.4f / %.4f)\n", lv, lv - 1, conv.max, conv.p95,
                 back.max, back.p95);
            putf(csv, "loop_level_change,%d,%.4f,%.4f,,%.4f,\n", lv, std::max(conv.max, back.max), std::max(conv.p95, back.p95), conv.p50);
        }
        if (lv == L) {
            writeObj(out + "/facets_L.obj", F.mesh);
            writeObj(out + "/loop_L.obj", S.mesh);
            const MeshQuery q(F.mesh);
            std::ofstream dv(out + "/loop_L_displacement.csv", std::ios::binary);
            dv << "vertex,distance_to_facets_mm,signed_mm\n";
            for (size_t v = 0; v < S.mesh.pos.size(); ++v) {
                const MeshQuery::Hit hh = q.closest(S.mesh.pos[v]);
                putf(dv, "%zu,%.5f,%.5f\n", v, norm(hh.p - S.mesh.pos[v]), hh.tri >= 0 ? dot(S.mesh.pos[v] - hh.p, F.mesh.faceNormal(hh.tri)) : 0.0);
            }
        }
        prevLoop = S;
    }
    // triangulation: the same polygons with the other diagonals (a file written by the driver)
    if (!a.get("--alt-mesh").empty()) {
        MarkedFile g;
        if (!readMarkedFile(a.get("--alt-mesh"), g, &err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        // the alternative file's piece with the same vertices (by box and count)
        int best = -1;
        Vec3 lo0, hi0;
        piece.bounds(lo0, hi0);
        double bestD = 1e300;
        for (int k = 0; k < static_cast<int>(g.pieces->pieces.size()); ++k) {
            TriMesh m = g.pieces->pieces[k].mesh;
            if (m.tri.size() != piece.tri.size()) continue;
            m.transform(mm, Vec3{0, 0, 0});
            Vec3 lo, hi;
            m.bounds(lo, hi);
            const double d = norm(lo - lo0) + norm(hi - hi0);
            if (d < bestD) {
                bestD = d;
                best = k;
            }
        }
        if (best < 0) putf(rep, "alternative triangulation: no piece with %zu triangles\n", piece.tri.size());
        else {
            TriMesh alt = g.pieces->pieces[best].mesh;
            alt.transform(mm, Vec3{0, 0, 0});
            const RefinedMesh S = refineLoop(piece, L), T = refineLoop(alt, L);
            const Quantiles st = distances(S.mesh, T.mesh), ts = distances(T.mesh, S.mesh);
            const RefinedMesh FA = refineMidpoint(alt, L);
            const Quantiles ff = distances(FA.mesh, refineMidpoint(piece, L).mesh);
            putf(rep, "alternative triangulation (the other diagonal of every quad; piece %d of the alternative file, box offset %.2g mm): its facets lie within %.4f mm of the "
                      "import's (the same polygons, so ~0 when planar quads); its Loop limit against the import's: max %.3f p95 %.3f mm (back %.3f / %.3f)\n",
                 best, bestD, ff.max, st.max, st.p95, ts.max, ts.p95);
            putf(csv, "loop_alt_triangulation,%d,%.4f,%.4f,,%.4f,\n", L, std::max(st.max, ts.max), std::max(st.p95, ts.p95), st.p50);
            writeObj(out + "/loop_L_alt.obj", T.mesh);
        }
    }
    std::printf("prep-compare: written to %s\n", out.c_str());
    return 0;
}

// --------------------------------------------------------------------------------------------------------- smooth-fixture
// A prism of `sides` sides, radius r, height Z, `rings` + 1 rings of corner vertices only; quads split along one diagonal
// (alternating per row when `alternate`).
TriMesh coarsePrism(int sides, double r, double Z, int rings, bool alternate) {
    TriMesh m;
    for (int k = 0; k <= rings; ++k)
        for (int j = 0; j < sides; ++j) m.pos.push_back({r * std::cos(2 * kPi * j / sides), r * std::sin(2 * kPi * j / sides), Z * k / rings});
    auto id = [&](int k, int j) { return k * sides + ((j % sides) + sides) % sides; };
    for (int k = 0; k < rings; ++k)
        for (int j = 0; j < sides; ++j) {
            if (!alternate || (k + j) % 2 == 0) {
                m.tri.push_back({id(k, j), id(k, j + 1), id(k + 1, j + 1)});
                m.tri.push_back({id(k, j), id(k + 1, j + 1), id(k + 1, j)});
            } else {
                m.tri.push_back({id(k, j), id(k, j + 1), id(k + 1, j)});
                m.tri.push_back({id(k, j + 1), id(k + 1, j + 1), id(k + 1, j)});
            }
        }
    m.computeNormals();
    return m;
}

TriMesh fineCylinder(double r, double Z, int around, double step) {
    TriMesh m;
    const int nz = std::max(1, static_cast<int>(std::lround(Z / step)));
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j < around; ++j) m.pos.push_back({r * std::cos(2 * kPi * j / around), r * std::sin(2 * kPi * j / around), Z * k / nz});
    auto id = [&](int k, int j) { return k * around + ((j % around) + around) % around; };
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < around; ++j) {
            m.tri.push_back({id(k, j), id(k, j + 1), id(k + 1, j + 1)});
            m.tri.push_back({id(k, j), id(k + 1, j + 1), id(k + 1, j)});
        }
    m.computeNormals();
    return m;
}

struct TubeRun {
    RelaxProduct p;
    WalkResidual wr;
    std::vector<Vec3> start;
    double radiusSpread = 0, radiusMean = 0, maxDeviation = 0;
};

// A ring tube (around x rows, gauge w x h) projected onto `surface` and relaxed there (CP3 settings), its border on the
// surface's border loops; the walk residual measured on the same surface.
TubeRun relaxTube(const TriMesh& surface, int around, int rows, double w, double h, double radius, double z0, std::vector<RelaxFace>& faces, TypedStitchMesh& tube) {
    TubeRun t;
    TriMesh unusedRef;
    std::vector<DeclaredBoundary> decl;
    std::string err;
    ringTubeFixture(std::vector<int>(rows + 1, around), std::vector<std::vector<int>>(rows), radius, h, tube, unusedRef, decl, &err);
    const MeshQuery q(surface);
    for (Vec3& p : tube.mesh.pos) {
        p.z += z0;
        p = q.closest(p).p;
    }
    VirtualShift vs;
    buildVirtualShift(tube, ShiftReading::PathSideDown, vs, &err);
    faces.clear();
    for (int f = 0; f < tube.mesh.faces(); ++f) faces.push_back(RelaxFace{vs.corners[f], tube.role[f]});
    std::vector<char> border(tube.mesh.pos.size(), 0);
    for (const std::vector<int>& l : polygonBorderLoops(tube.mesh))
        for (int v : l) border[v] = 1;
    const std::vector<std::vector<int>> loops = triangleBoundaryLoops(surface);
    MeshRelaxRequest r;
    r.w = w;
    r.h = h;
    t.start = tube.mesh.pos;
    relaxStitchMesh(surface, loops, faces, tube.mesh.pos, border, r, t.p, &err);
    t.wr = walkStationarity(surface, faces, t.p.pos, border, loops, r, nullptr, 1e-3);
    std::vector<double> rad;
    for (const Vec3& p : t.p.pos) rad.push_back(std::sqrt(p.x * p.x + p.y * p.y));
    const Quantiles qr = quantiles(rad);
    t.radiusSpread = qr.max - qr.min;
    double s = 0;
    for (double x : rad) s += x;
    t.radiusMean = rad.empty() ? 0 : s / rad.size();
    for (const Vec3& p : t.p.pos) t.maxDeviation = std::max(t.maxDeviation, norm(q.closest(p).p - p));
    return t;
}

int commandSmoothFixture(const Opts& a) {
    const std::string out = a.get("--out");
    fs::create_directories(out);
    const int sides = 16, around = 17, rows = 8;
    const double radius = 15, w = 5.5, h = 3.4, Z = rows * h;
    std::ofstream rep(out + "/smooth_fixture.txt", std::ios::binary);
    std::ofstream csv(out + "/smooth_fixture.csv", std::ios::binary);
    csv << "surface,triangulation,levels,triangles,crease_max_deg,surface_radius_spread_mm,surface_radius_mean_mm,status,stop,steps,rms_force_end,walk_rms,walk_max,walk_on_edge,"
           "tube_radius_spread_mm,max_vertex_from_previous_level_mm,max_vertex_from_true_cylinder_mm\n";
    putf(rep, "coarse prism: %d sides, radius %.1f mm (corners), height %.1f mm, 2 rows of quads; tube %d around x %d rows, gauge %.1f x %.1f (CP3 settings)\n", sides, radius, Z,
         around, rows, w, h);
    for (int alt = 0; alt < 2; ++alt) {
        const TriMesh coarse = coarsePrism(sides, radius, Z, 2, alt == 1);
        std::vector<Vec3> previous;
        std::vector<Vec3> cylinderResult;
        std::vector<std::pair<std::string, TriMesh>> surfaces;
        // the reference cylinder: the Loop limit's mean radius at the mid-height ring (filled in below)
        for (int lv = 1; lv <= 4; ++lv) surfaces.push_back({"loop_limit", refineLoop(coarse, lv).mesh});
        surfaces.insert(surfaces.begin(), {"prism_facets", refineMidpoint(coarse, 4).mesh});
        double meanR = 0;
        {
            const TriMesh& s4 = surfaces.back().second;
            int n = 0;
            for (const Vec3& p : s4.pos)
                if (std::fabs(p.z - 0.5 * Z) < 0.3 * Z) {
                    meanR += std::sqrt(p.x * p.x + p.y * p.y);
                    ++n;
                }
            meanR /= std::max(1, n);
        }
        surfaces.push_back({"true_cylinder", fineCylinder(meanR, Z, 384, 0.425)});
        // the true cylinder first, so every run is compared with it
        std::vector<RelaxFace> faces;
        TypedStitchMesh tube;
        const TubeRun ref = relaxTube(surfaces.back().second, around, rows, w, h, meanR, 0, faces, tube);
        cylinderResult = ref.p.pos;
        for (size_t si = 0; si < surfaces.size(); ++si) {
            const std::string& name = surfaces[si].first;
            const TriMesh& s = surfaces[si].second;
            const int lv = name == "loop_limit" ? static_cast<int>(si) : name == "prism_facets" ? 4 : 0;
            double cmax = 0, c95 = 0;
            dihedral(s, cmax, c95);
            std::vector<double> rad;
            for (const Vec3& p : s.pos)
                if (std::fabs(p.z - 0.5 * Z) < 0.3 * Z) rad.push_back(std::sqrt(p.x * p.x + p.y * p.y));
            const Quantiles qr = quantiles(rad);
            const TubeRun t = relaxTube(s, around, rows, w, h, meanR, 0, faces, tube);
            double fromPrev = -1, fromCyl = 0;
            if (name == "loop_limit" && !previous.empty()) {
                fromPrev = 0;
                for (size_t v = 0; v < t.p.pos.size(); ++v) fromPrev = std::max(fromPrev, norm(t.p.pos[v] - previous[v]));
            }
            for (size_t v = 0; v < t.p.pos.size(); ++v) fromCyl = std::max(fromCyl, norm(t.p.pos[v] - cylinderResult[v]));
            if (name == "loop_limit") previous = t.p.pos;
            putf(csv, "%s,%s,%d,%zu,%.3f,%.4f,%.4f,%s,%s,%d,%.4e,%.4e,%.4e,%d,%.4f,%.4f,%.4f\n", name.c_str(), alt ? "alternating" : "uniform", lv, s.tri.size(), cmax,
                 qr.max - qr.min, qr.count ? 0.5 * (qr.max + qr.min) : 0.0, toString(t.p.status), t.p.stopReason.c_str(), t.p.steps, t.p.rmsForceEnd, t.wr.rms, t.wr.max,
                 t.wr.onEdge, t.radiusSpread, fromPrev, fromCyl);
            putf(rep, "%-14s %-11s level %d: %6zu triangles, creases max %6.2f deg, surface radius spread %.3f mm; relaxation %s (%s, %d steps), walk residual rms %.2e max %.2e "
                      "(%d on an edge); tube radius spread %.3f mm; %s%.4f mm from the true-cylinder result\n",
                 name.c_str(), alt ? "alternating" : "uniform", lv, s.tri.size(), cmax, qr.max - qr.min, toString(t.p.status), t.p.stopReason.c_str(), t.p.steps, t.wr.rms,
                 t.wr.max, t.wr.onEdge, t.radiusSpread, fromPrev >= 0 ? ("previous level " + std::to_string(fromPrev).substr(0, 6) + " mm, ").c_str() : "", fromCyl);
            if (alt == 0 && name == "loop_limit" && lv == 3) writeObj(out + "/loop3_uniform.obj", s);
        }
        if (alt == 0) writeObj(out + "/coarse_prism.obj", coarse);
    }
    std::printf("smooth-fixture: written to %s\n", out.c_str());
    return 0;
}

// ------------------------------------------------------------------------------------------------------- shaping-fixtures
TriMesh lathe(const std::vector<double>& rho, const std::vector<double>& z, int around, int sub) { return revolutionThroughRings(rho, z, around, sub); }

struct ShapingCase {
    std::string name, purpose;
    std::vector<int> counts;
    std::vector<std::vector<int>> shaping;
    bool boundaryShaping = false;
    bool fitted = true;   // ring circumferences = counts x w (false: one cylinder for all rings)
};

int commandShapingFixtures(const Opts& a) {
    const std::string out = a.get("--out");
    fs::create_directories(out);
    const double w = 5.5, h = 3.4;
    std::vector<ShapingCase> cases = {
        {"K_tube", "regular knit only: the K bounds' baseline", {24, 24, 24, 24, 24}, {{}, {}, {}, {}}},
        {"INC_isolated", "one kp in an interior strip (24 -> 25)", {24, 24, 25, 25, 25}, {{}, {12}, {}, {}}},
        {"DEC_isolated", "one d12k in an interior strip (25 -> 24)", {25, 25, 24, 24, 24}, {{}, {12}, {}, {}}},
        {"INC4", "four kp in one strip (24 -> 28), spaced 6 apart", {24, 24, 28, 28, 28}, {{}, {3, 9, 15, 21}, {}, {}}},
        {"DEC4", "four d12k in one strip (28 -> 24), spaced 6 apart", {28, 28, 24, 24, 24}, {{}, {2, 8, 14, 20}, {}, {}}},
        // added after the arm's results located its failing shaping cells in densely shaped strips (pole start, narrow finish)
        {"INC_every_other", "dense: five kp on every other stitch (12 -> 17)", {12, 12, 17, 17, 17}, {{}, {1, 3, 5, 7, 9}, {}, {}}},
        {"INC_adjacent_pairs", "dense: two pairs of neighbouring kp (16 -> 20; the adjacent-shaping variant's pattern)", {16, 16, 20, 20, 20}, {{}, {2, 3, 9, 10}, {}, {}}},
        {"DEC_every_other", "dense: five d12k on every other stitch (17 -> 12)", {17, 17, 12, 12, 12}, {{}, {1, 3, 5, 7, 9}, {}, {}}},
        {"INC_isolated_cylinder", "one kp, every ring on one cylinder (circumference 24.5 w): the reference does not fit the counts", {24, 24, 25, 25, 25}, {{}, {12}, {}, {}}, false, false},
        {"kp_cast_on", "kp + cast-on (boundary-shaping variant): 24 -> 25 in the first strip", {24, 25, 25, 25}, {{12}, {}, {}}, true},
        {"d12k_bind_off", "d12k + bind-off (boundary-shaping variant): 25 -> 24 in the last strip", {25, 25, 25, 24}, {{}, {}, {12}}, true},
    };
    std::ofstream csv(out + "/shaping_fixtures.csv", std::ios::binary);
    csv << "case,face,kind,old_within,old_outside,min_course_w,max_course_w,min_wale_h,max_wale_h,min_corner,max_corner,min_sub_area_wh,folded,new_interfaces,new_sub_quads,new_within,"
           "new_outside,course_course_angle,old_wale_slant_min,relax_status,walk_rms\n";
    std::ofstream rep(out + "/shaping_fixtures.txt", std::ios::binary);
    putf(rep, "gauge %.1f x %.1f; CP3 model and settings; rings on surfaces of revolution with circumference counts x w and rows h apart along the surface (unless named)\n", w, h);
    TypedStitchMesh incRelaxed;
    TriMesh incReference;
    for (const ShapingCase& c : cases) {
        TypedStitchMesh m;
        TriMesh unusedRef;
        std::vector<DeclaredBoundary> decl;
        std::string err;
        if (!ringTubeFixture(c.counts, c.shaping, 20, h, m, unusedRef, decl, &err, c.boundaryShaping)) {
            putf(rep, "%s: not built: %s\n", c.name.c_str(), err.c_str());
            continue;
        }
        // the rings' radii and heights: circumference n w, consecutive rings h apart along the surface
        std::vector<double> rho, z;
        double mean = 0;
        for (int n : c.counts) mean += n;
        mean /= c.counts.size();
        for (size_t k = 0; k < c.counts.size(); ++k) {
            rho.push_back((c.fitted ? c.counts[k] : mean) * w / (2 * kPi));
            z.push_back(k == 0 ? 0.0 : z.back() + std::sqrt(std::max(0.0, h * h - (rho[k] - rho[k - 1]) * (rho[k] - rho[k - 1]))));
        }
        const TriMesh ref = lathe(rho, z, 192, 4);
        const MeshQuery q(ref);
        for (size_t k = 0; k < m.rings.size(); ++k)
            for (int v : m.rings[k]) {
                const double phi = std::atan2(m.mesh.pos[v].y, m.mesh.pos[v].x);
                m.mesh.pos[v] = q.closest(Vec3{rho[k] * std::cos(phi), rho[k] * std::sin(phi), z[k]}).p;
            }
        VirtualShift vs;
        buildVirtualShift(m, ShiftReading::PathSideDown, vs, &err);
        std::vector<RelaxFace> faces;
        for (int f = 0; f < m.mesh.faces(); ++f) faces.push_back(RelaxFace{vs.corners[f], m.role[f]});
        std::vector<char> border(m.mesh.pos.size(), 0);
        for (const std::vector<int>& l : polygonBorderLoops(m.mesh))
            for (int v : l) border[v] = 1;
        const std::vector<std::vector<int>> loops = triangleBoundaryLoops(ref);
        MeshRelaxRequest r;
        r.w = w;
        r.h = h;
        RelaxProduct p;
        relaxStitchMesh(ref, loops, faces, m.mesh.pos, border, r, p, &err);
        const WalkResidual wr = walkStationarity(ref, faces, p.pos, border, loops, r, nullptr, 1e-3);
        TypedStitchMesh relaxed = m;
        relaxed.mesh.pos = p.pos;
        const std::string dir = out + "/" + c.name;
        fs::create_directories(dir);
        lab::writePolygonObj(dir + "/relaxed.obj", relaxed);
        lab::writePolygonObj(dir + "/start.obj", m);
        writeObj(dir + "/reference.obj", ref);
        const std::vector<CellMeasure> old = measureCells(relaxed, w, h);
        const std::vector<ShapingCellMeasure> sc = measureShapingCells(relaxed, w, h);
        int oldIn = 0, newIn = 0, kIn = 0, kN = 0;
        for (const CellMeasure& cm : old)
            if (cm.kind == StitchFaceKind::Regular) {
                ++kN;
                kIn += cm.withinBounds;
            }
        for (const ShapingCellMeasure& s : sc) {
            const CellMeasure& o = old[s.face];
            oldIn += o.withinBounds;
            newIn += s.within;
            putf(csv, "%s,%d,%s,%d,\"%s\",%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%.4f,%d,%d,%d,%d,\"%s\",%.2f,%.4f,%s,%.3e\n", c.name.c_str(), s.face, toString(s.kind), o.withinBounds ? 1 : 0,
                 o.outside.c_str(), s.minCourse, s.maxCourse, s.minWale, s.maxWale, s.minCorner, s.maxCorner, s.minSubArea, s.foldedSubQuads, s.interfacesWithin ? 1 : 0,
                 s.subQuadsWithin ? 1 : 0, s.within ? 1 : 0, s.outside.c_str(), o.courseCourseAngle, o.minWale, toString(p.status), wr.rms);
        }
        putf(rep, "%-22s %s: relaxation %s (%s, %d steps), walk residual rms %.2e; K within %d/%d; shaping cells %zu: frozen pentagon bounds %d within, sub-quad assessment %d within\n",
             c.name.c_str(), c.purpose.c_str(), toString(p.status), p.stopReason.c_str(), p.steps, wr.rms, kIn, kN, sc.size(), oldIn, newIn);
        for (const ShapingCellMeasure& s : sc) {
            const CellMeasure& o = old[s.face];
            putf(rep, "    face %d %s: frozen [%s] course-course %.1f deg, wale/slant %.3f-%.3f | sub-quads: course/w %.3f-%.3f, wale/h %.3f-%.3f, corners %.1f-%.1f deg, area >= %.3f wh, folded %d [%s]\n",
                 s.face, toString(s.kind), o.withinBounds ? "within" : o.outside.c_str(), o.courseCourseAngle, o.minWale, o.maxWale, s.minCourse, s.maxCourse, s.minWale, s.maxWale,
                 s.minCorner, s.maxCorner, s.minSubArea, s.foldedSubQuads, s.within ? "within" : s.outside.c_str());
        }
        if (c.name == "INC_isolated") {
            incRelaxed = relaxed;
            incReference = ref;
        }
    }
    // deliberately defective cells, made from the relaxed isolated increase: each must be rejected by the assessment
    int kp = -1;
    for (int f = 0; f < incRelaxed.mesh.faces(); ++f)
        if (incRelaxed.kind[f] == StitchFaceKind::Increase) kp = f;
    if (kp >= 0) {
        const PolygonMesh& pm = incRelaxed.mesh;
        const int n = pm.sides(kp);
        int b0 = -1, b1 = -1, mid = -1, t0 = -1;
        for (int k = 0; k < n; ++k) {
            const StitchEdgeRole before = incRelaxed.role[kp][(k + n - 1) % n], after = incRelaxed.role[kp][k];
            if (after == StitchEdgeRole::Bottom) {
                b0 = pm.vertexAt(kp, k);
                b1 = pm.vertexAt(kp, (k + 1) % n);
            }
            if (before == StitchEdgeRole::Top && after == StitchEdgeRole::Top) mid = pm.vertexAt(kp, k);
            if (before == StitchEdgeRole::Top && after == StitchEdgeRole::Wale) t0 = pm.vertexAt(kp, k);
        }
        struct Defect {
            std::string name, what;
            int vertex;
            Vec3 to;
        };
        const Vec3 bm = 0.5 * (pm.pos[b0] + pm.pos[b1]);
        const Vec3 up = pm.pos[mid] - bm;
        std::vector<Defect> defects = {
            {"kp_folded", "the two top edges' shared vertex moved below the bottom edge (0.3 h): both sub-quads fold", mid, bm - up * (0.3 * h / std::max(1e-9, norm(up)))},
            {"kp_collapsed_top", "a top edge shortened to 0.3 w (its outer vertex moved towards the shared one)", t0, pm.pos[mid] + (pm.pos[t0] - pm.pos[mid]) * (0.3 * w / std::max(1e-9, norm(pm.pos[t0] - pm.pos[mid])))},
            {"kp_sheared", "the shared top vertex pushed 1.2 w sideways along the row: one sub-quad strongly sheared", mid, pm.pos[mid] + (pm.pos[b1] - pm.pos[b0]) * (1.2 * w / std::max(1e-9, norm(pm.pos[b1] - pm.pos[b0])))},
        };
        // and the frozen bounds' own ideal: the trapezoid (bottom w, top 2 w, height h) placed in the cell's frame
        const Vec3 ex = (pm.pos[b1] - pm.pos[b0]) / std::max(1e-9, norm(pm.pos[b1] - pm.pos[b0]));
        Vec3 ey = up - ex * dot(up, ex);
        ey = ey / std::max(1e-9, norm(ey));
        putf(rep, "\ndefective cells (made from INC_isolated's relaxed kp, face %d), and the trapezoid the frozen bounds assume:\n", kp);
        for (int d = 0; d <= static_cast<int>(defects.size()); ++d) {
            TypedStitchMesh x = incRelaxed;
            std::string name, what;
            if (d < static_cast<int>(defects.size())) {
                x.mesh.pos[defects[d].vertex] = defects[d].to;
                name = defects[d].name;
                what = defects[d].what;
            } else {
                // trapezoid: b0 at the origin of the cell frame, b1 at w; top corners at (-w/2, h), (w/2, h), (3w/2, h)
                const Vec3 o = pm.pos[b0];
                int t1 = -1;
                for (int k = 0; k < n; ++k)
                    if (incRelaxed.role[kp][(k + n - 1) % n] == StitchEdgeRole::Wale && incRelaxed.role[kp][k] == StitchEdgeRole::Top) t1 = pm.vertexAt(kp, k);
                x.mesh.pos[b1] = o + ex * w;
                x.mesh.pos[t1] = o + ex * (1.5 * w) + ey * h;
                x.mesh.pos[mid] = o + ex * (0.5 * w) + ey * h;
                x.mesh.pos[t0] = o + ex * (-0.5 * w) + ey * h;
                name = "kp_trapezoid";
                what = "the frozen bounds' ideal trapezoid (bottom w, top 2 w, height h; positive control for both measures)";
            }
            const std::vector<CellMeasure> old = measureCells(x, w, h);
            for (const ShapingCellMeasure& s : measureShapingCells(x, w, h)) {
                if (s.face != kp) continue;
                const CellMeasure& o = old[s.face];
                putf(rep, "    %-18s %s: frozen bounds %s; assessment %s\n", name.c_str(), what.c_str(), o.withinBounds ? "within" : ("REJECT (" + o.outside + ")").c_str(),
                     s.within ? "within" : ("REJECT (" + s.outside + ")").c_str());
                putf(csv, "%s,%d,%s,%d,\"%s\",%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%.4f,%d,%d,%d,%d,\"%s\",%.2f,%.4f,,\n", name.c_str(), s.face, toString(s.kind), o.withinBounds ? 1 : 0,
                     o.outside.c_str(), s.minCourse, s.maxCourse, s.minWale, s.maxWale, s.minCorner, s.maxCorner, s.minSubArea, s.foldedSubQuads, s.interfacesWithin ? 1 : 0,
                     s.subQuadsWithin ? 1 : 0, s.within ? 1 : 0, s.outside.c_str(), o.courseCourseAngle, o.minWale);
            }
            lab::writePolygonObj(out + "/" + name + ".obj", x);
        }
    }
    std::printf("shaping-fixtures: written to %s\n", out.c_str());
    return 0;
}

} // namespace

// 2 for a command this file does not have.
int smoothCommand(const std::string& cmd, int argc, char** argv) {
    const Opts a = parse(argc, argv);
    if (cmd == "prep-compare") return commandPrepCompare(a);
    if (cmd == "smooth-fixture") return commandSmoothFixture(a);
    if (cmd == "shaping-fixtures") return commandShapingFixtures(a);
    return 2;
}
