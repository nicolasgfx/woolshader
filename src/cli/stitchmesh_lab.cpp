// stitchmesh_lab: the command line of the stockinette stitch-mesh experiment.
// Checkpoint 1, input analysis only, no stitch mesh:
//   stitchmesh_lab surface --out <dir> [--jitter <mm> --seed <n>] [--segments 192 --spacing 1.0]
//       the reference surface tube-taper-v1 (reference.obj, reference_vertices.csv, reference.json, components.json,
//       guides.json)
//   stitchmesh_lab analyze --dir <dir> --case D|I --w <mm> --h <mm> --out <dir>
//       Method A: the time field from the case's cast-on boundary, its error against the exact meridian arc length,
//       the rings at h' with their contours, lengths, cut crossings and counts, and production's count scheduler on
//       the same mid-strip lengths (time_field.csv, rings.csv, ring_points.csv, mid_strips.csv, analyze_report.txt)
//   stitchmesh_lab field --dir <dir> --mode solved|guided --spacing <mm> --seed <n> --out <dir> [--sweeps 20]
//       Method B: RIM's 2-RoSy orientation field (field.csv, singularities.csv, field_report.txt)
// Checkpoint 2, typed meshes before relaxation:
//   stitchmesh_lab fixtures --out <dir>
//       the construction fixtures and negative controls (construction_fixtures.h), each checked; exit 1 when one is not
//       as expected (fixtures.csv, fixtures_report.txt, one folder per case)
//   stitchmesh_lab build --dir <ref dir> --case D|I --method rows|w18 --w <mm> --h <mm> --out <dir>
//       Method A (guided rows) or B (W18 pipeline): the typed mesh, its check, the virtual shift, coverage and cells
//       (candidate.obj, candidate_bindings.csv, connectivity.csv, face_types.csv, rows.csv, loops.csv, virtual_shift.csv,
//       coverage.csv, cells_before.csv, sections.csv, shaping.csv (A), report.json, build_report.txt)
//   stitchmesh_lab load --case-dir <dir> --dir <ref dir>
//       re-read the typed mesh's files, re-run the check and the virtual shift, compare the digest (load_report.txt)
// Checkpoint 3, relaxation (relax_mesh.h, Yuksel et al. 2012 section 6.1):
//   stitchmesh_lab relax-fixtures --out <dir>
//       the solver on hand-built tubes: rest state, perturbation, virtual helix, kp / d12k, the stall rule
//   stitchmesh_lab relax --case-dir <CP2 run> --dir <ref dir> --w <mm> --h <mm> --out <dir>
//       Method A on its virtual faces, Method B on its stored faces; before/after checks, coverage and cells
// Checkpoint 4, the teddy's prepared components for Method A (analysis only, no stitch mesh):
//   stitchmesh_lab component --manifest <pieces/manifest.txt> --out <dir> [--gauge-piece <file> --stitches 65 --w 5.5 --h 3.4]
//       per piece: validity, poles, hidden regions (inside another piece), the time field from each pole, its level
//       contours and events, extra maxima, ring counts, the proposed cut guide and the proposed direction
// Surface markings, on an imported OBJ or PLY as the Mark step reads it:
//   stitchmesh_lab pieces --mesh <file> [--out <dir>]          the file's pieces (pieces.csv)
//   stitchmesh_lab snap --mesh <file> --points <txt> --out <file.markings> [--piece <i>]
//       markings from 3D points (scripts, proposals)
//   stitchmesh_lab marked --mesh <file> --markings <file> --unit-mm <mm> --w <mm> --h <mm> [--piece <i>] [--depth built|relaxed] --out <dir>
//       Method A on each marked piece: review.json, stitch_mesh.obj, rings.obj, traces.obj, regions.csv, log.txt
// Every file is written here; the components do no IO.
#include <algorithm>
#include <cstdarg>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <chrono>

#include "io/marking_file.h"
#include "io/mesh_read.h"
#include "jobs/marking_session.h"
#include "knitting/knit/pieces.h"
#include "knitting/stitch/stitch_program.h"
#include "model/mesh.h"
#include "model/mesh_check.h"
#include "model/hash.h"
#include "model/mesh_query.h"
#include "model/surface_path.h"
#include "model/surface_pieces.h"
#include "stitchmesh/construction_fixtures.h"
#include "stitchmesh/coverage.h"
#include "stitchmesh/directions.h"
#include "stitchmesh/guided_rows.h"
#include "stitchmesh/labels.h"
#include "stitchmesh/marked_build.h"
#include "stitchmesh/orientation_field.h"
#include "stitchmesh/reference_surface.h"
#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/remesh_stage.h"
#include "stitchmesh/stitch_mesh.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"
#include "stitchmesh/revolution_surface.h"
#include "stitchmesh/row_rings.h"
#include "stitchmesh/time_field.h"
#include "stitchmesh_lab_shared.h"

using namespace knit;
namespace fs = std::filesystem;

namespace {

struct Args {
    std::map<std::string, std::string> values;
    std::string get(const std::string& key, const std::string& fallback = "") const {
        const auto it = values.find(key);
        return it == values.end() ? fallback : it->second;
    }
    double number(const std::string& key, double fallback) const {
        const auto it = values.find(key);
        return it == values.end() ? fallback : std::stod(it->second);
    }
};

Args parseArgs(int argc, char** argv) {
    Args a;
    for (int i = 2; i + 1 < argc; i += 2) a.values[argv[i]] = argv[i + 1];
    return a;
}

// printf into a stream (C varargs through vsnprintf; the lab's lines are short).
void putf(std::ostream& o, const char* fmt, ...) {
    char buf[4096];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    o << buf;
}

std::string hex64(uint64_t v) {
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", (unsigned long long) v);
    return b;
}

double quantile(std::vector<double> v, double q) {
    if (v.empty()) return 0;
    std::sort(v.begin(), v.end());
    const size_t i = (size_t) std::min((double) v.size() - 1, std::floor(q * (double) (v.size() - 1) + 0.5));
    return v[i];
}

// Boundary loops of a triangle mesh, each as its vertices in walking order.
std::vector<std::vector<int>> boundaryLoops(const TriMesh& m) {
    std::map<std::pair<int, int>, int> uses;
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) uses[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}] += 1;
    std::map<int, int> nextOf;   // directed as in the triangle
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) {
            const int a = t[k], b = t[(k + 1) % 3];
            if (uses[{std::min(a, b), std::max(a, b)}] == 1) nextOf[a] = b;
        }
    std::vector<std::vector<int>> loops;
    std::set<int> seen;
    for (const auto& [v0, unused] : nextOf) {
        (void) unused;
        if (seen.count(v0)) continue;
        std::vector<int> loop;
        int v = v0;
        while (!seen.count(v) && nextOf.count(v)) {
            seen.insert(v);
            loop.push_back(v);
            v = nextOf[v];
        }
        loops.push_back(loop);
    }
    return loops;
}

double loopLength(const TriMesh& m, const std::vector<int>& loop) {
    double L = 0;
    for (size_t i = 0; i < loop.size(); ++i) L += norm(m.pos[loop[(i + 1) % loop.size()]] - m.pos[loop[i]]);
    return L;
}

double meanZ(const TriMesh& m, const std::vector<int>& loop) {
    double z = 0;
    for (int v : loop) z += m.pos[v].z;
    return loop.empty() ? 0 : z / loop.size();
}

bool writeObj(const std::string& path, const TriMesh& m) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    f << "# stitchmesh_lab reference surface, mm, +z up\n";
    for (const Vec3& p : m.pos) putf(f, "v %.17g %.17g %.17g\n", p.x, p.y, p.z);
    for (const auto& t : m.tri) putf(f, "f %d %d %d\n", t[0] + 1, t[1] + 1, t[2] + 1);
    return (bool) f;
}

// The tube recipe as frozen (frozen_settings.json: tube-taper-v1).
RevolutionProfile tubeProfile() { return RevolutionProfile{}; }

int commandSurface(const Args& a) {
    const std::string out = a.get("--out");
    if (out.empty()) {
        std::fprintf(stderr, "surface: --out <dir> is required\n");
        return 2;
    }
    fs::create_directories(out);
    RevolutionMeshParams p;
    p.segments = (int) a.number("--segments", 192);
    p.ringSpacing = a.number("--spacing", 1.0);
    p.jitter = a.number("--jitter", 0);
    p.jitterSeed = (uint64_t) a.number("--seed", 0);
    const RevolutionProfile prof = tubeProfile();
    RevolutionMesh rm;
    std::string err;
    if (!makeRevolutionMesh(prof, p, rm, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const MeshCheck check = checkMesh(rm.mesh);
    const std::vector<std::vector<int>> loops = boundaryLoops(rm.mesh);
    if (!writeObj(out + "/reference.obj", rm.mesh)) {
        std::fprintf(stderr, "cannot write %s/reference.obj\n", out.c_str());
        return 1;
    }
    {
        std::ofstream f(out + "/reference_vertices.csv");
        f << "v,arc,z,phi\n";
        char b[160];
        for (size_t i = 0; i < rm.arc.size(); ++i) {
            std::snprintf(b, sizeof b, "%zu,%.17g,%.17g,%.17g\n", i, rm.arc[i], rm.height[i], rm.phi[i]);
            f << b;
        }
    }
    double maxSurfaceError = 0;   // every vertex on the exact surface: | |xy| - r(z) |
    for (const Vec3& q : rm.mesh.pos) maxSurfaceError = std::max(maxSurfaceError, std::fabs(std::hypot(q.x, q.y) - prof.radius(q.z)));
    double minEdge = 1e300, maxEdge = 0, minAngle = 180, maxAngle = 0;
    for (const auto& t : rm.mesh.tri)
        for (int k = 0; k < 3; ++k) {
            const Vec3 A = rm.mesh.pos[t[k]], B = rm.mesh.pos[t[(k + 1) % 3]], C = rm.mesh.pos[t[(k + 2) % 3]];
            minEdge = std::min(minEdge, norm(B - A));
            maxEdge = std::max(maxEdge, norm(B - A));
            const double ang = std::acos(std::clamp(dot(normalized(B - A), normalized(C - A)), -1.0, 1.0)) * 180.0 / kPi;
            minAngle = std::min(minAngle, ang);
            maxAngle = std::max(maxAngle, ang);
        }
    const ComponentCheck& c0 = check.components.front();
    std::ofstream j(out + "/reference.json", std::ios::binary);
    putf(j, "{\n  \"recipe\": \"tube-taper-v1\",\n  \"units\": \"mm\",\n  \"axis\": \"+z\",\n");
    putf(j, "  \"profile\": {\"r0\": %.17g, \"r1\": %.17g, \"z0\": %.17g, \"z1\": %.17g, \"z2\": %.17g},\n", prof.r0, prof.r1, prof.z0, prof.z1, prof.z2);
    putf(j, "  \"tessellation\": {\"segments\": %d, \"ring_spacing\": %.17g, \"rings\": %d, \"stagger\": \"ring i turned by i*pi/segments\"},\n", p.segments, p.ringSpacing, rm.rings);
    putf(j, "  \"variant\": {\"jitter_mm\": %.17g, \"seed\": %llu},\n", p.jitter, (unsigned long long) p.jitterSeed);
    putf(j, "  \"vertices\": %zu,\n  \"triangles\": %zu,\n  \"digest\": \"%s\",\n", rm.mesh.pos.size(), rm.mesh.tri.size(), hex64(meshDigest(rm.mesh)).c_str());
    putf(j, "  \"meridian_length\": %.17g,\n  \"exact_area\": %.17g,\n  \"mesh_area\": %.17g,\n", rm.meridianLength, rm.exactArea, rm.mesh.area());
    putf(j, "  \"max_vertex_surface_error\": %.6g,\n  \"edge_min_max\": [%.6g, %.6g],\n  \"angle_min_max\": [%.6g, %.6g],\n", maxSurfaceError, minEdge, maxEdge, minAngle, maxAngle);
    putf(j, "  \"check\": {\"components\": %zu, \"non_manifold_edges\": %d, \"inconsistent_edges\": %d, \"boundary_edges\": %d, \"euler\": %d, \"consistent\": %s, \"self_intersection\": \"%s\"},\n",
                 check.components.size(), c0.nonManifoldEdges, c0.inconsistentEdges, c0.boundaryEdges, c0.euler, c0.consistent ? "true" : "false", toString(c0.selfIntersection));
    putf(j, "  \"boundary_loops\": [");
    for (size_t k = 0; k < loops.size(); ++k)
        putf(j, "%s{\"vertices\": %zu, \"mean_z\": %.6f, \"length\": %.6f}", k ? ", " : "", loops[k].size(), meanZ(rm.mesh, loops[k]), loopLength(rm.mesh, loops[k]));
    putf(j, "]\n}\n");
    
    // the one prepared component and its guides (the cut guide as a half-plane, see guides.json)
    std::ofstream cj(out + "/components.json", std::ios::binary);
    putf(cj, "{\n  \"components\": [\n    {\"id\": \"tube\", \"source\": \"procedural (tube-taper-v1)\", \"triangles\": \"all\",\n"
                     "     \"boundaries\": [{\"id\": \"z0\", \"loop_mean_z\": 0, \"role_case_D\": \"cast-on\", \"role_case_I\": \"bind-off\"},\n"
                     "                    {\"id\": \"z150\", \"loop_mean_z\": 150, \"role_case_D\": \"bind-off\", \"role_case_I\": \"cast-on\"}],\n"
                     "     \"closures\": [], \"side_boundaries\": [], \"cuts\": [\"cut-phi0\"]}\n  ],\n  \"connections\": []\n}\n");
    
    std::ofstream gj(out + "/guides.json", std::ios::binary);
    putf(gj, "{\n  \"cut-phi0\": {\"kind\": \"half-plane\", \"origin\": [0, 0, 0], \"normal\": [0, 1, 0], \"side\": [1, 0, 0],\n"
                     "    \"meaning\": \"row origin and shift-path location (Method A: vertex 0 of every ring); the reference mesh is NOT cut: "
                     "the pairing across the guide is the mesh's own connectivity, so no seam is introduced\",\n"
                     "    \"limit\": \"a half-plane guide fits a surface of revolution; teddy components need a surface polyline guide (CP4)\"},\n"
                     "  \"growth\": {\"D\": \"away from the z = 0 boundary\", \"I\": \"away from the z = 150 boundary\"},\n"
                     "  \"row_axis_method_B\": {\"kind\": \"ring-axis\", \"axis\": [0, 0, 1]},\n"
                     "  \"shaping\": \"none (shaping follows the count schedule)\"\n}\n");
    
    std::printf("surface: %zu vertices, %zu triangles, digest %s, %zu boundary loops, max vertex surface error %.3g mm\n",
                rm.mesh.pos.size(), rm.mesh.tri.size(), hex64(meshDigest(rm.mesh)).c_str(), loops.size(), maxSurfaceError);
    return 0;
}

bool readReference(const std::string& dir, TriMesh& m, std::vector<double>& arc, std::string& digestSaved, double& meridian, std::string* err) {
    if (!readObj(dir + "/reference.obj", m, err)) return false;
    m.computeNormals();
    std::ifstream f(dir + "/reference_vertices.csv");
    std::string line;
    std::getline(f, line);
    while (std::getline(f, line)) {
        const size_t comma = line.find(',');   // v,arc,z,phi: the arc is the second column
        if (comma != std::string::npos) arc.push_back(std::strtod(line.c_str() + comma + 1, nullptr));
    }
    std::ifstream jf(dir + "/reference.json");
    const std::string json((std::istreambuf_iterator<char>(jf)), std::istreambuf_iterator<char>());
    const size_t d = json.find("\"digest\": \"");
    digestSaved = d == std::string::npos ? "" : json.substr(d + 11, 16);
    const size_t ml = json.find("\"meridian_length\": ");
    meridian = ml == std::string::npos ? 0 : std::stod(json.substr(ml + 19));
    if (arc.size() != m.pos.size()) {
        if (err) *err = "reference_vertices.csv does not match reference.obj";
        return false;
    }
    return true;
}

int commandAnalyze(const Args& a) {
    const std::string dir = a.get("--dir"), out = a.get("--out"), cs = a.get("--case", "D");
    const double w = a.number("--w", 5.5), h = a.number("--h", 5.5);
    if (dir.empty() || out.empty() || (cs != "D" && cs != "I")) {
        std::fprintf(stderr, "analyze: --dir, --out and --case D|I are required\n");
        return 2;
    }
    fs::create_directories(out);
    TriMesh m;
    std::vector<double> arc;
    std::string digestSaved, err;
    double T = 0;
    if (!readReference(dir, m, arc, digestSaved, T, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const std::string digestRead = hex64(meshDigest(m));
    const RevolutionProfile prof = tubeProfile();
    MeridianTable table;
    table.build(prof, 600001);
    // the case's cast-on boundary: the loop at z = 0 (D) or z = 150 (I), found from the mesh, not from the recipe
    std::vector<std::vector<int>> loops = boundaryLoops(m);
    std::sort(loops.begin(), loops.end(), [&](const std::vector<int>& x, const std::vector<int>& y) { return meanZ(m, x) < meanZ(m, y); });
    if (loops.size() != 2) {
        std::fprintf(stderr, "analyze: expected 2 boundary loops, found %zu\n", loops.size());
        return 1;
    }
    const std::vector<int>& castOn = cs == "D" ? loops.front() : loops.back();
    const std::vector<int>& bindOff = cs == "D" ? loops.back() : loops.front();
    std::vector<double> t;
    TimeFieldReport tr;
    if (!geodesicDistance(m, castOn, t, tr, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    // errors against the exact distance (the meridian arc length from the cast-on boundary: meridians are geodesics
    // that meet the parallels at right angles)
    std::vector<double> absErr;
    double maxErr = 0, sumErr = 0, maxRel = 0;
    int worstV = 0;
    {
        std::ofstream f(out + "/time_field.csv");
        f << "v,t,t_exact,error\n";
        char b[160];
        for (size_t v = 0; v < m.pos.size(); ++v) {
            const double ex = cs == "D" ? arc[v] : T - arc[v];
            const double e = t[v] - ex;
            absErr.push_back(std::fabs(e));
            sumErr += e;
            if (std::fabs(e) > maxErr) {
                maxErr = std::fabs(e);
                worstV = (int) v;
            }
            if (ex > 1) maxRel = std::max(maxRel, std::fabs(e) / ex);
            std::snprintf(b, sizeof b, "%zu,%.12g,%.12g,%.6g\n", v, t[v], ex, e);
            f << b;
        }
    }
    double tBindMin = 1e300, tBindMax = -1e300, tBindMean = 0;
    for (int v : bindOff) {
        tBindMin = std::min(tBindMin, t[v]);
        tBindMax = std::max(tBindMax, t[v]);
        tBindMean += t[v];
    }
    tBindMean /= bindOff.size();
    const double Tf = tBindMin;   // the field's extent: the last level that is inside the mesh everywhere
    const int R = std::max(1, (int) std::lround(Tf / h));
    const double hp = Tf / R;
    const Vec3 cutOrigin{0, 0, 0}, cutNormal{0, 1, 0}, cutSide{1, 0, 0};
    struct RingRow {
        int k = 0;
        double level = 0;
        int contours = 0, closedContours = 0, points = 0, cutCrossings = 0;
        double L = 0, Lexact = 0, cutArc = -1, cutZ = 0, arcMin = 0, arcMax = 0, arcMean = 0;
        int n = 0;
        bool boundary = false;
    };
    std::vector<RingRow> rows;
    std::ofstream pts(out + "/ring_points.csv");
    pts << "k,contour,i,x,y,z,a,b,f,tri\n";
    auto exactArcAt = [&](const RingPoint& rp) {   // exact arc length of a crossing point: along its edge, linearly
        const double sa = cs == "D" ? arc[rp.a] : T - arc[rp.a], sb = cs == "D" ? arc[rp.b] : T - arc[rp.b];
        return sa + rp.f * (sb - sa);
    };
    for (int k = 0; k <= R; ++k) {
        RingRow row;
        row.k = k;
        row.level = k * hp;
        if (k == 0 || k == R) {   // the boundary loops themselves
            const std::vector<int>& loop = k == 0 ? castOn : bindOff;
            row.boundary = true;
            row.contours = row.closedContours = 1;
            row.points = (int) loop.size();
            row.L = loopLength(m, loop);
            double s0 = 0;
            for (int v : loop) s0 += cs == "D" ? arc[v] : T - arc[v];
            row.arcMean = row.arcMin = row.arcMax = s0 / loop.size();
            row.cutArc = row.arcMean;
            row.cutCrossings = 1;
            for (size_t i = 0; i < loop.size(); ++i) {
                const Vec3& p = m.pos[loop[i]];
                char b[200];
                std::snprintf(b, sizeof b, "%d,0,%zu,%.9g,%.9g,%.9g,%d,%d,0,-1\n", k, i, p.x, p.y, p.z, loop[i], loop[i]);
                pts << b;
            }
        } else {
            const RingLevel lv = extractLevel(m, t, row.level);
            row.contours = (int) lv.contours.size();
            for (size_t c = 0; c < lv.contours.size(); ++c) {
                const Contour& con = lv.contours[c];
                row.closedContours += con.closed ? 1 : 0;
                for (size_t i = 0; i < con.points.size(); ++i) {
                    const RingPoint& rp = con.points[i];
                    char b[200];
                    std::snprintf(b, sizeof b, "%d,%zu,%zu,%.9g,%.9g,%.9g,%d,%d,%.9g,%d\n", k, c, i, rp.p.x, rp.p.y, rp.p.z, rp.a, rp.b, rp.f, rp.tri);
                    pts << b;
                }
            }
            if (lv.contours.size() == 1) {
                const Contour& con = lv.contours.front();
                row.points = (int) con.points.size();
                row.L = con.length;
                double smin = 1e300, smax = -1e300, ssum = 0;
                for (const RingPoint& rp : con.points) {
                    const double s = exactArcAt(rp);
                    smin = std::min(smin, s);
                    smax = std::max(smax, s);
                    ssum += s;
                }
                row.arcMin = smin;
                row.arcMax = smax;
                row.arcMean = ssum / con.points.size();
                const double where = cutCrossing(con, cutOrigin, cutNormal, cutSide, &row.cutCrossings);
                if (where >= 0) {
                    const size_t i = (size_t) where;
                    const double f = where - (double) i;
                    const RingPoint& p0 = con.points[i];
                    const RingPoint& p1 = con.points[(i + 1) % con.points.size()];
                    row.cutArc = exactArcAt(p0) + f * (exactArcAt(p1) - exactArcAt(p0));
                    row.cutZ = p0.p.z + f * (p1.p.z - p0.p.z);
                }
            }
        }
        const double zExact = table.heightAt(cs == "D" ? row.arcMean : T - row.arcMean);
        row.Lexact = 2.0 * kPi * prof.radius(zExact);
        row.n = row.L > 0 ? (int) std::lround(row.L / w) : 0;
        rows.push_back(row);
    }
    // the |dn| <= floor(n/2) limit (checked, not assumed)
    int limitBinds = 0;
    for (int k = 0; k < R; ++k)
        if (std::abs(rows[k + 1].n - rows[k].n) > rows[k].n / 2) ++limitBinds;
    {
        std::ofstream f(out + "/rings.csv");
        f << "k,level,boundary,contours,closed_contours,points,L,L_exact,L_over_L_exact,n,mean_course_over_w,cut_crossings,cut_arc,cut_z,spacing_over_hp,arc_min,arc_max,arc_spread,strip_dn\n";
        char b[512];
        for (size_t k = 0; k < rows.size(); ++k) {
            const RingRow& r = rows[k];
            const double spacing = k > 0 && r.cutArc >= 0 && rows[k - 1].cutArc >= 0 ? (r.cutArc - rows[k - 1].cutArc) / hp : 0;
            const int dn = k + 1 < rows.size() ? rows[k + 1].n - r.n : 0;
            std::snprintf(b, sizeof b, "%d,%.9g,%d,%d,%d,%d,%.9g,%.9g,%.9g,%d,%.9g,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.6g,%d\n", r.k, r.level, r.boundary ? 1 : 0, r.contours,
                          r.closedContours, r.points, r.L, r.Lexact, r.Lexact > 0 ? r.L / r.Lexact : 0, r.n, r.n > 0 ? r.L / (r.n * w) : 0, r.cutCrossings, r.cutArc, r.cutZ, spacing,
                          r.arcMin, r.arcMax, r.arcMax - r.arcMin, dn);
            f << b;
        }
    }
    // production's count scheduler (src/knitting/stitch/stitch_program.cpp scheduleCounts, max shaping 0.5) on the
    // lengths at the middle of every strip (production's measureRows takes each row's length at its middle)
    std::vector<double> targets;
    std::vector<double> midL;
    int midEvents = 0;
    for (int k = 0; k < R; ++k) {
        const RingLevel lv = extractLevel(m, t, (k + 0.5) * hp);
        const double L = lv.contours.size() == 1 ? lv.contours.front().length : 0;
        if (lv.contours.size() != 1) ++midEvents;
        midL.push_back(L);
        targets.push_back(L / w);
    }
    const std::vector<int> prodCounts = scheduleCounts(rows.front().n, targets, 0.5);
    {
        std::ofstream f(out + "/mid_strips.csv");
        f << "strip,mid_level,mid_L,target,production_n\n";
        char b[200];
        for (int k = 0; k < R; ++k) {
            std::snprintf(b, sizeof b, "%d,%.9g,%.9g,%.9g,%d\n", k, (k + 0.5) * hp, midL[k], targets[k], k < (int) prodCounts.size() ? prodCounts[k] : -1);
            f << b;
        }
    }
    // summary
    int events = 0, ringsSingleClosed = 0, cutOk = 0, spacingOut = 0, meanOut = 0;
    double spacingMin = 1e9, spacingMax = 0, worstMean = 0, worstLratio = 0, worstSpread = 0;
    for (size_t k = 0; k < rows.size(); ++k) {
        const RingRow& r = rows[k];
        if (r.contours == 1 && r.closedContours == 1) ++ringsSingleClosed; else ++events;
        if (r.cutCrossings == 1) ++cutOk;
        if (k > 0 && r.cutArc >= 0 && rows[k - 1].cutArc >= 0) {
            const double sp = (r.cutArc - rows[k - 1].cutArc) / hp;
            spacingMin = std::min(spacingMin, sp);
            spacingMax = std::max(spacingMax, sp);
            if (sp < 0.97 || sp > 1.03) ++spacingOut;
        }
        if (r.n > 0) {
            const double mc = r.L / (r.n * w);
            worstMean = std::max(worstMean, std::fabs(mc - 1));
            if (mc < 0.95 || mc > 1.05) ++meanOut;
        }
        if (r.Lexact > 0) worstLratio = std::max(worstLratio, std::fabs(r.L / r.Lexact - 1));
        if (!r.boundary) worstSpread = std::max(worstSpread, r.arcMax - r.arcMin);
    }
    int shaping = 0, prodShaping = 0, maxDn = 0;
    for (int k = 0; k < R; ++k) {
        shaping += std::abs(rows[k + 1].n - rows[k].n);
        maxDn = std::max(maxDn, std::abs(rows[k + 1].n - rows[k].n));
    }
    for (size_t k = 1; k < prodCounts.size(); ++k) prodShaping += std::abs(prodCounts[k] - prodCounts[k - 1]);
    std::ofstream rp(out + "/analyze_report.txt", std::ios::binary);
    putf(rp, "analyze: case %s, gauge w %.4g h %.4g\n", cs.c_str(), w, h);
    putf(rp, "reference digest read %s, saved %s: %s\n", digestRead.c_str(), digestSaved.c_str(), digestRead == digestSaved ? "equal" : "DIFFERENT");
    putf(rp, "time field: %d vertices accepted of %zu, %d triangle updates, %d edge updates, %.1f ms\n", tr.accepted, m.pos.size(), tr.triangleUpdates, tr.edgeUpdates, tr.ms);
    putf(rp, "time field error vs exact arc length [mm]: mean %+.4g, p50 %.4g, p95 %.4g, max %.4g at vertex %d (z %.3f); max relative (t > 1 mm) %.4g\n",
                 sumErr / m.pos.size(), quantile(absErr, 0.5), quantile(absErr, 0.95), maxErr, worstV, m.pos[worstV].z, maxRel);
    putf(rp, "field on the bind-off loop: min %.6f, mean %.6f, max %.6f (exact %.6f); extent T_field = %.6f\n", tBindMin, tBindMean, tBindMax, T, Tf);
    putf(rp, "rows R = %d, h' = %.6f (%+.3f %% of h)\n", R, hp, 100 * (hp / h - 1));
    putf(rp, "rings: %zu (incl. the 2 boundary loops); one closed contour: %d; contour events: %d; cut crossed exactly once: %d\n", rows.size(), ringsSingleClosed, events, cutOk);
    putf(rp, "ring spacing along the cut / h': min %.5f, max %.5f; outside [0.97, 1.03]: %d\n", spacingMin, spacingMax, spacingOut);
    putf(rp, "ring length / exact parallel: worst |ratio - 1| %.5f; worst exact-arc spread within one ring %.5f mm\n", worstLratio, worstSpread);
    putf(rp, "counts: n from %d to %d, shaping events %d, max per strip %d, |dn| limit binds %d; mean course / w worst |ratio - 1| %.5f, outside [0.95, 1.05]: %d\n",
                 rows.front().n, rows.back().n, shaping, maxDn, limitBinds, worstMean, meanOut);
    putf(rp, "production scheduleCounts on mid-strip lengths: %zu rows, from %d to %d, shaping events %d; mid-strip contour events %d\n",
                 prodCounts.size(), prodCounts.empty() ? 0 : prodCounts.front(), prodCounts.empty() ? 0 : prodCounts.back(), prodShaping, midEvents);
    
    std::printf("analyze %s: R %d h' %.4f, field max error %.4g mm, events %d, shaping %d\n", cs.c_str(), R, hp, maxErr, events, shaping);
    return 0;
}

int commandField(const Args& a) {
    const std::string dir = a.get("--dir"), out = a.get("--out"), mode = a.get("--mode", "solved");
    if (dir.empty() || out.empty() || (mode != "solved" && mode != "guided")) {
        std::fprintf(stderr, "field: --dir, --out and --mode solved|guided are required\n");
        return 2;
    }
    fs::create_directories(out);
    TriMesh m;
    std::vector<double> arc;
    std::string digestSaved, err;
    double T = 0;
    if (!readReference(dir, m, arc, digestSaved, T, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    OrientationFieldRequest r;
    r.spacing = a.number("--spacing", 11);
    r.seed = (uint64_t) a.number("--seed", 1);
    r.sweepsPerLevel = (int) a.number("--sweeps", 20);
    const std::vector<Vec3> ring = ringAxisGuide(m, Vec3{0, 0, 1});
    if (mode == "guided") r.guide = ring;
    OrientationField of;
    if (!computeOrientationField(m, r, of, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    std::vector<double> dev90, dev180;
    {
        std::ofstream f(out + "/field.csv");
        f << "v,qx,qy,qz,dev_lattice_deg,dev_line_deg\n";
        char b[200];
        for (size_t v = 0; v < m.pos.size(); ++v) {
            const double d90 = directionDeviation(of.q[v], ring[v], m.vnormal[v], 90);
            const double d180 = directionDeviation(of.q[v], ring[v], m.vnormal[v], 180);
            dev90.push_back(d90);
            dev180.push_back(d180);
            std::snprintf(b, sizeof b, "%zu,%.9g,%.9g,%.9g,%.6g,%.6g\n", v, of.q[v].x, of.q[v].y, of.q[v].z, d90, d180);
            f << b;
        }
    }
    {
        std::ofstream f(out + "/singularities.csv");
        f << "i,x,y,z\n";
        for (size_t i = 0; i < of.singularities.size(); ++i)
            f << i << "," << of.singularities[i].x << "," << of.singularities[i].y << "," << of.singularities[i].z << "\n";
    }
    int over5 = 0;
    for (double d : dev90) over5 += d > 5 ? 1 : 0;
    uint64_t fieldHash = 1469598103934665603ull;
    for (const Vec3& q : of.q)
        for (double c : {q.x, q.y, q.z}) {
            uint64_t bits;
            std::memcpy(&bits, &c, 8);
            fieldHash = (fieldHash ^ bits) * 1099511628211ull;
        }
    std::ofstream rp(out + "/field_report.txt", std::ios::binary);
    putf(rp, "field: mode %s, rosy 2, intrinsic, spacing %.4g mm, seed %llu, sweeps per level %d\n", mode.c_str(), r.spacing, (unsigned long long) r.seed, r.sweepsPerLevel);
    putf(rp, "RIM finest vertices %d (input %zu), levels", of.finestVertices, m.pos.size());
    for (int s : of.levelSizes) putf(rp, " %d", s);
    putf(rp, "\ntime %.1f ms; orientation singularities %zu; field digest %s\n", of.ms, of.singularities.size(), hex64(fieldHash).c_str());
    putf(rp, "deviation from the rows-around-z direction, modulo 90 deg (the lattice): mean %.4g, p50 %.4g, p95 %.4g, max %.4g deg; vertices over 5 deg: %d of %zu\n",
                 [&] { double s = 0; for (double d : dev90) s += d; return s / dev90.size(); }(), quantile(dev90, 0.5), quantile(dev90, 0.95), quantile(dev90, 1.0), over5, dev90.size());
    putf(rp, "deviation modulo 180 deg (which axis q is): p50 %.4g, max %.4g deg\n", quantile(dev180, 0.5), quantile(dev180, 1.0));
    
    std::printf("field %s: %zu singularities, lattice deviation p95 %.3g deg, max %.3g deg, %.0f ms\n", mode.c_str(), of.singularities.size(), quantile(dev90, 0.95), quantile(dev90, 1.0), of.ms);
    return 0;
}


// "a,b,c" as three integers.
bool threeInts(const std::string& line, int& a, int& b, int& c) {
    const char* p = line.c_str();
    char* e = nullptr;
    a = (int) std::strtol(p, &e, 10);
    if (e == p || *e != ',') return false;
    p = e + 1;
    b = (int) std::strtol(p, &e, 10);
    if (e == p || *e != ',') return false;
    p = e + 1;
    c = (int) std::strtol(p, &e, 10);
    return e != p;
}

// ---------------------------------------------------------------------------------------------------------------------
// Checkpoint 2: typed meshes, fixtures, reload
// ---------------------------------------------------------------------------------------------------------------------

const char* kindGroup(StitchFaceKind k) {
    switch (k) {
    case StitchFaceKind::Regular: return "K";
    case StitchFaceKind::Increase: return "INC_kp";
    case StitchFaceKind::Decrease: return "DEC_d12k";
    case StitchFaceKind::SpecialQuad: return "special_quad";
    case StitchFaceKind::PoleStart: return "pole_start";
    case StitchFaceKind::Unsupported: return "unresolved";
    }
    return "?";
}

const char* roleWord(StitchEdgeRole r) {
    switch (r) {
    case StitchEdgeRole::Wale: return "W";
    case StitchEdgeRole::Bottom: return "B";
    case StitchEdgeRole::Top: return "T";
    }
    return "?";
}

std::string constructionOf(const TypedStitchMesh& s, int f) {
    std::string c;
    switch (s.kind[f]) {
    case StitchFaceKind::Regular: c = "K"; break;
    case StitchFaceKind::Increase: c = "INC kp"; break;
    case StitchFaceKind::Decrease: c = "DEC d12k"; break;
    case StitchFaceKind::SpecialQuad: c = "special quad"; break;
    case StitchFaceKind::PoleStart: c = "pole start"; break;
    case StitchFaceKind::Unsupported: return "unresolved";
    }
    bool castOn = false, bindOff = false;
    for (EdgeBoundary b : s.boundary[f]) {
        castOn = castOn || b == EdgeBoundary::CastOn;
        bindOff = bindOff || b == EdgeBoundary::BindOff;
    }
    if (castOn) c += " + cast-on";
    if (bindOff) c += " + bind-off";
    return c;
}

struct EdgeNeighbours {
    std::map<std::pair<int, int>, std::vector<std::pair<int, int>>> uses;   // edge -> (face, corner)
    std::pair<int, int> across(const PolygonMesh& m, int f, int k) const {
        const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % m.sides(f));
        const auto it = uses.find({std::min(a, b), std::max(a, b)});
        if (it == uses.end() || it->second.size() != 2) return {-1, -1};
        return it->second[0].first == f && it->second[0].second == k ? it->second[1] : it->second[0];
    }
};

EdgeNeighbours neighboursOf(const PolygonMesh& m) {
    EdgeNeighbours n;
    for (int f = 0; f < m.faces(); ++f)
        for (int k = 0; k < m.sides(f); ++k) {
            const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % m.sides(f));
            n.uses[{std::min(a, b), std::max(a, b)}].push_back({f, k});
        }
    return n;
}

// The typed mesh's files; positions bit for bit (%.17g), faces in stored order with a group line
// wherever the kind changes.
void writeTypedMesh(const std::string& out, const TypedStitchMesh& s, const TypedMeshCheck& c, const VirtualShift* vs, const VirtualShiftCheck* vc) {
    const PolygonMesh& m = s.mesh;
    const EdgeNeighbours nb = neighboursOf(m);
    {
        std::ofstream f(out + "/candidate.obj", std::ios::binary);
        f << "# stitchmesh_lab typed stitch mesh, mm; groups by kind: K, INC_kp, DEC_d12k, special_quad, unresolved\n";
        for (const Vec3& p : m.pos) putf(f, "v %.17g %.17g %.17g\n", p.x, p.y, p.z);
        const char* group = "";
        for (int fc = 0; fc < m.faces(); ++fc) {
            if (std::strcmp(group, kindGroup(s.kind[fc])) != 0) {
                group = kindGroup(s.kind[fc]);
                f << "g " << group << "\n";
            }
            f << "f";
            for (int k = 0; k < m.sides(fc); ++k) f << " " << m.vertexAt(fc, k) + 1;
            f << "\n";
        }
    }
    {
        std::ofstream f(out + "/candidate_bindings.csv", std::ios::binary);
        f << "v,tri,b0,b1,b2,distance\n";
        for (size_t v = 0; v < s.binding.size(); ++v) {
            const VertexBinding& b = s.binding[v];
            putf(f, "%zu,%d,%.17g,%.17g,%.17g,%.6g\n", v, b.tri, b.bary[0], b.bary[1], b.bary[2], b.distance);
        }
    }
    {
        std::ofstream f(out + "/connectivity.csv", std::ios::binary);
        f << "face,corner,a,b,role,boundary,neighbour_face,neighbour_corner\n";
        for (int fc = 0; fc < m.faces(); ++fc)
            for (int k = 0; k < m.sides(fc); ++k) {
                const std::pair<int, int> o = nb.across(m, fc, k);
                putf(f, "%d,%d,%d,%d,%s,%s,%d,%d\n", fc, k, m.vertexAt(fc, k), m.vertexAt(fc, (k + 1) % m.sides(fc)), roleWord(s.role[fc][k]), toString(s.boundary[fc][k]),
                     o.first, o.second);
            }
    }
    std::vector<int> posInStrand(m.faces(), -1);
    for (size_t i = 0; i < c.strands.strands.size(); ++i)
        for (size_t j = 0; j < c.strands.strands[i].size(); ++j) posInStrand[c.strands.strands[i][j]] = (int) j;
    {
        std::ofstream f(out + "/face_types.csv", std::ios::binary);
        f << "face,sides,kind,construction,strip,slot,strand,strand_pos,cx,cy,cz,reason\n";
        for (int fc = 0; fc < m.faces(); ++fc) {
            const Vec3 ce = m.faceCentroid(fc);
            putf(f, "%d,%d,%s,%s,%d,%d,%d,%d,%.9g,%.9g,%.9g,%s\n", fc, m.sides(fc), kindGroup(s.kind[fc]), constructionOf(s, fc).c_str(),
                 s.strip.empty() ? -1 : s.strip[fc], s.slot.empty() ? -1 : s.slot[fc], c.strands.strandOfFace.empty() ? -1 : c.strands.strandOfFace[fc], posInStrand[fc],
                 ce.x, ce.y, ce.z, fc < (int) s.reason.size() ? s.reason[fc].c_str() : "");
        }
    }
    {
        std::ofstream f(out + "/rows.csv", std::ios::binary);
        f << "strand,closed,faces,face_list\n";
        for (size_t i = 0; i < c.strands.strands.size(); ++i) {
            putf(f, "%zu,%d,%zu,", i, (int) c.strands.closed[i], c.strands.strands[i].size());
            for (size_t j = 0; j < c.strands.strands[i].size(); ++j) f << (j ? " " : "") << c.strands.strands[i][j];
            f << "\n";
        }
    }
    {
        std::ofstream f(out + "/loops.csv", std::ios::binary);
        f << "a,b,kind,offering_face,receiving_face\n";
        for (const auto& [key, list] : nb.uses) {
            const int f0 = list[0].first, k0 = list[0].second;
            if (s.role[f0][k0] == StitchEdgeRole::Wale && (list.size() < 2 || s.role[list[1].first][list[1].second] == StitchEdgeRole::Wale)) continue;
            std::string kind;
            int offering = -1, receiving = -1;
            if (list.size() == 1) {
                const bool top = s.role[f0][k0] == StitchEdgeRole::Top;
                kind = std::string("border-") + (top ? "top-" : "bottom-") + toString(s.boundary[f0][k0]);
                (top ? offering : receiving) = f0;
            } else if (list.size() == 2) {
                const int f1 = list[1].first, k1 = list[1].second;
                const StitchEdgeRole r0 = s.role[f0][k0], r1 = s.role[f1][k1];
                if (r0 == StitchEdgeRole::Wale || r1 == StitchEdgeRole::Wale) kind = "course-against-wale";
                else if (r0 == r1) kind = r0 == StitchEdgeRole::Top ? "both-top" : "both-bottom";
                else {
                    kind = "interlock";
                    offering = r0 == StitchEdgeRole::Top ? f0 : f1;
                    receiving = r0 == StitchEdgeRole::Top ? f1 : f0;
                }
            } else {
                kind = "non-manifold";
            }
            putf(f, "%d,%d,%s,%d,%d\n", key.first, key.second, kind.c_str(), offering, receiving);
        }
    }
    {
        std::ofstream f(out + "/ring_vertices.csv", std::ios::binary);
        f << "ring,j,v\n";
        for (size_t k = 0; k < s.rings.size(); ++k)
            for (size_t j = 0; j < s.rings[k].size(); ++j) putf(f, "%zu,%zu,%d\n", k, j, s.rings[k][j]);
    }
    {
        std::ofstream f(out + "/shift_path.csv", std::ios::binary);
        f << "strip,bottom,top\n";
        for (const ShiftPathEdge& e : s.shiftPath) putf(f, "%d,%d,%d\n", e.strip, e.bottom, e.top);
    }
    if (vs && vc) {
        std::vector<int> order(m.faces(), -1);
        for (size_t i = 0; i < vc->strandOrder.size(); ++i) order[vc->strandOrder[i]] = (int) i;
        std::ofstream f(out + "/virtual_shift.csv", std::ios::binary);
        putf(f, "# rule: %s\n", vs->rule.c_str());
        f << "face,strip,slot,re_read,helix_order,virtual_corners,stored_corners\n";
        for (int fc = 0; fc < m.faces(); ++fc) {
            putf(f, "%d,%d,%d,%d,%d,", fc, s.strip.empty() ? -1 : s.strip[fc], s.slot.empty() ? -1 : s.slot[fc], (int) vs->reRead[fc], order[fc]);
            for (size_t k = 0; k < vs->corners[fc].size(); ++k) f << (k ? " " : "") << vs->corners[fc][k];
            f << ",";
            for (int k = 0; k < m.sides(fc); ++k) f << (k ? " " : "") << m.vertexAt(fc, k);
            f << "\n";
        }
    }
}

void writeCells(const std::string& out, const std::vector<CellMeasure>& cells) {
    std::ofstream f(out + "/cells_before.csv", std::ios::binary);
    f << "face,kind,width_ratio,height_ratio,min_course,max_course,min_wale,max_wale,min_angle,max_angle,course_course_angle,area_over_wh,within,outside\n";
    for (const CellMeasure& c : cells)
        putf(f, "%d,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%s\n", c.face, kindGroup(c.kind), c.widthRatio, c.heightRatio, c.minCourse, c.maxCourse, c.minWale,
             c.maxWale, c.minAngle, c.maxAngle, c.courseCourseAngle, c.area, c.withinBounds ? 1 : 0, c.outside.c_str());
}

struct CellSummary {
    int faces = 0, within = 0;
    std::map<std::string, std::pair<int, int>> byKind;   // kind -> (faces, within)
    std::vector<double> width, height, minAngle, maxAngle;
};

CellSummary summarize(const std::vector<CellMeasure>& cells) {
    CellSummary s;
    for (const CellMeasure& c : cells) {
        ++s.faces;
        s.within += c.withinBounds ? 1 : 0;
        auto& bk = s.byKind[kindGroup(c.kind)];
        bk.first += 1;
        bk.second += c.withinBounds ? 1 : 0;
        if (c.kind == StitchFaceKind::Regular) {
            s.width.push_back(c.widthRatio);
            s.height.push_back(c.heightRatio);
            s.minAngle.push_back(c.minAngle);
            s.maxAngle.push_back(c.maxAngle);
        }
    }
    return s;
}

void writeVerdicts(std::ostream& o, const TypedMeshCheck& c, const CoverageReport* cov, const VirtualShiftCheck* vc, const std::string& shiftNote,
                   const CellSummary& cs) {
    o << c.text();
    if (vc) o << vc->text(); else o << "virtual shift: " << shiftNote << "\n";
    if (cov) o << cov->text();
    putf(o, "cells (frozen bounds, reported, not a CP2 verdict): %d of %d faces within bounds", cs.within, cs.faces);
    for (const auto& [k, n] : cs.byKind) putf(o, "; %s %d/%d", k.c_str(), n.second, n.first);
    o << "\n";
    if (!cs.width.empty())
        putf(o, "  K width ratio p5 %.3f p50 %.3f p95 %.3f [min %.3f max %.3f]; height ratio p5 %.3f p50 %.3f p95 %.3f [min %.3f max %.3f]; angles min %.1f max %.1f deg\n",
             quantile(cs.width, 0.05), quantile(cs.width, 0.5), quantile(cs.width, 0.95), quantile(cs.width, 0), quantile(cs.width, 1), quantile(cs.height, 0.05),
             quantile(cs.height, 0.5), quantile(cs.height, 0.95), quantile(cs.height, 0), quantile(cs.height, 1), quantile(cs.minAngle, 0), quantile(cs.maxAngle, 1));
}

// The component verdict: PASS when connectivity and coverage pass with no unresolved face.
Status componentVerdict(const TypedMeshCheck& c, const CoverageReport* cov) {
    if (c.connectivity == Status::Fail || (cov && cov->status == Status::Fail)) return Status::Fail;
    if (c.connectivity != Status::Pass || c.unsupported > 0 || (cov && cov->status != Status::Pass)) return Status::Unresolved;
    return Status::Pass;
}

int commandFixtures(const Args& a) {
    const std::string out = a.get("--out");
    if (out.empty()) {
        std::fprintf(stderr, "fixtures: --out <dir> is required\n");
        return 2;
    }
    fs::create_directories(out);
    // --set arm: the imported arm's pole-start and adjacent-shaping fixtures (arm/), apart from the CP2 set
    const std::vector<FixtureCase> cases = a.get("--set") == "arm" ? armConstructionFixtureCases() : constructionFixtureCases();
    std::ofstream sum(out + "/fixtures.csv", std::ios::binary);
    sum << "name,negative,entry,expected_connectivity,connectivity,expected_shift,shift,as_expected,faces,K,INC,DEC,special,unsupported,interlocks,cast_on,bind_off,"
           "strands,closed,open,yarn_ends,virtual_strands,virtual_longest,detail\n";
    std::ofstream rp(out + "/fixtures_report.txt", std::ios::binary);
    int unexpected = 0;
    for (const FixtureCase& fc : cases) {
        const std::string dir = out + "/" + fc.name;
        fs::create_directories(dir);
        TypedCheckOptions o;
        o.w = 5.5;
        o.h = 5.5;
        o.boundaryTolerance = 1e-6;
        o.bindingTolerance = 0.01;
        o.faceDeviation = 1e9;   // fixtures are topological: their faces are flat facets of coarse references
        o.degenerateArea = 0.01;
        const TypedMeshCheck c = checkTypedMesh(fc.mesh, fc.reference, fc.declared, o);
        VirtualShift vs;
        VirtualShiftCheck vc;
        std::string err;
        const bool shift = fc.runShift && buildVirtualShift(fc.mesh, fc.reading, vs, &err);
        if (shift) vc = validateVirtualShift(fc.mesh, vs);
        writeTypedMesh(dir, fc.mesh, c, shift ? &vs : nullptr, shift ? &vc : nullptr);
        writeCells(dir, measureCells(fc.mesh, 5.5, 5.5));
        const Status shiftStatus = shift ? vc.status : Status::Unresolved;
        const bool ok = c.connectivity == fc.expectConnectivity && shiftStatus == fc.expectShift;
        unexpected += ok ? 0 : 1;
        {
            std::ofstream cf(dir + "/check.txt", std::ios::binary);
            cf << fc.name << ": " << fc.purpose << "\n" << c.text();
            if (shift) cf << "rule: " << vs.rule << "\n" << vc.text();
        }
        putf(rp, "%s %s [%s]: connectivity %s (expected %s), shift %s (expected %s) -> %s\n  %s\n", fc.negative ? "NEGATIVE" : "fixture ", fc.name.c_str(), fc.entry.c_str(),
             toString(c.connectivity), toString(fc.expectConnectivity), shift ? toString(vc.status) : "not run", shift ? toString(fc.expectShift) : "not run",
             ok ? "as expected" : "NOT AS EXPECTED", fc.purpose.c_str());
        rp << "  " << c.detail << "\n";
        if (shift) rp << "  " << vc.text();
        std::string detail = c.detail + (shift ? " | shift: " + vc.detail : "");
        std::replace(detail.begin(), detail.end(), ',', ';');
        std::string entry = fc.entry;
        std::replace(entry.begin(), entry.end(), ',', ';');
        putf(sum, "%s,%d,%s,%s,%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%zu,%d,%d,%d,%d,%d,%s\n", fc.name.c_str(), fc.negative ? 1 : 0, entry.c_str(),
             toString(fc.expectConnectivity), toString(c.connectivity), shift ? toString(fc.expectShift) : "not run", shift ? toString(vc.status) : "not run", ok ? 1 : 0,
             c.faces, c.regular, c.increase, c.decrease, c.special, c.unsupported, c.interlocks, c.castOnEdges, c.bindOffEdges, c.strands.strands.size(), c.closedStrands,
             c.openStrands, c.strands.yarnEnds, shift ? vc.strands : -1, shift ? vc.longestStrand : -1, detail.c_str());
        // the special-quad neighbourhood: every edge of every special quad, the face across it and the yarn through it
        if (fc.name == "special-quad-neighbourhood") {
            std::ofstream nf(dir + "/special_quad_neighbourhood.csv", std::ios::binary);
            nf << "special_quad,corner,role,neighbour,neighbour_kind,neighbour_role,strand,strand_closed,strand_faces,pos_in_strand,prev_face,next_face\n";
            const EdgeNeighbours nb = neighboursOf(fc.mesh.mesh);
            for (int f = 0; f < fc.mesh.mesh.faces(); ++f) {
                if (fc.mesh.kind[f] != StitchFaceKind::SpecialQuad) continue;
                const int st = c.strands.strandOfFace[f];
                const std::vector<int>& faces = c.strands.strands[st];
                const int pos = (int) (std::find(faces.begin(), faces.end(), f) - faces.begin());
                const int n = (int) faces.size();
                const int prev = c.strands.closed[st] || pos > 0 ? faces[(pos + n - 1) % n] : -1;
                const int next = c.strands.closed[st] || pos + 1 < n ? faces[(pos + 1) % n] : -1;
                for (int k = 0; k < 4; ++k) {
                    const std::pair<int, int> o2 = nb.across(fc.mesh.mesh, f, k);
                    putf(nf, "%d,%d,%s,%d,%s,%s,%d,%d,%d,%d,%d,%d\n", f, k, roleWord(fc.mesh.role[f][k]), o2.first, o2.first >= 0 ? kindGroup(fc.mesh.kind[o2.first]) : "-",
                         o2.first >= 0 ? roleWord(fc.mesh.role[o2.first][o2.second]) : "-", st, (int) c.strands.closed[st], n, pos, prev, next);
                }
            }
        }
    }
    putf(rp, "\n%zu cases, %d not as expected\n", cases.size(), unexpected);
    std::printf("fixtures: %zu cases, %d not as expected\n", cases.size(), unexpected);
    return unexpected == 0 ? 0 : 1;
}

double wrapAngle(double a) {
    while (a > kPi) a -= 2 * kPi;
    while (a < -kPi) a += 2 * kPi;
    return a;
}

// Method B: rows are the closed strands, ordered along the growth; the shift-path takes in each row the wale edge
// nearest the cut guide (azimuth 0) and must form a connected chain whose neighbouring faces are K. When it does,
// strips, slots and rings are assigned from it, so the same virtual shift and its validator run; otherwise the reason.
bool methodBShiftPath(TypedStitchMesh& s, const StrandTrace& t, bool growthUp, std::string& why) {
    const PolygonMesh& m = s.mesh;
    for (size_t i = 0; i < t.strands.size(); ++i)
        if (!t.closed[i]) {
            why = "not constructible: strand " + std::to_string(i) + " is open";
            return false;
        }
    std::vector<int> order(t.strands.size());
    std::vector<double> height(t.strands.size(), 0);
    for (size_t i = 0; i < t.strands.size(); ++i) {
        order[i] = (int) i;
        for (int f : t.strands[i]) height[i] += m.faceCentroid(f).z;
        height[i] /= t.strands[i].size();
    }
    std::sort(order.begin(), order.end(), [&](int x, int y) { return growthUp ? height[x] < height[y] : height[x] > height[y]; });
    const int R = (int) order.size();
    std::vector<ShiftPathEdge> path(R);
    std::vector<int> pathFace(R, -1);   // per row: the face whose W+ edge is the path edge
    for (int r = 0; r < R; ++r) {
        double best = 1e300;
        for (int f : t.strands[order[r]]) {
            const int ex = exitWaleCorner(s.role[f]);
            if (ex < 0) continue;
            const int a = m.vertexAt(f, ex), b = m.vertexAt(f, (ex + 1) % m.sides(f));
            const Vec3 mid = 0.5 * (m.pos[a] + m.pos[b]);
            const double d = std::fabs(wrapAngle(std::atan2(mid.y, mid.x)));
            if (d < best) {
                best = d;
                pathFace[r] = f;
                path[r] = ShiftPathEdge{r, a, b};
            }
        }
        if (pathFace[r] < 0) {
            why = "not constructible: row " + std::to_string(r) + " has no face with its own course direction";
            return false;
        }
    }
    int breaks = 0;
    for (int r = 1; r < R; ++r)
        if (path[r].bottom != path[r - 1].top) ++breaks;
    if (breaks > 0) {
        why = "not constructible: the per-row wale edges nearest the cut are not a connected chain (" + std::to_string(breaks) + " of " + std::to_string(R - 1) +
              " links broken); the rows stay independent closed rows";
        return false;
    }
    // strips and slots from the path; rings from each row's bottom edges and the last row's top edges
    s.strip.assign(m.faces(), -1);
    s.slot.assign(m.faces(), -1);
    s.rings.assign(R + 1, {});
    for (int r = 0; r < R; ++r) {
        const std::vector<int>& faces = t.strands[order[r]];
        const int n = (int) faces.size();
        const int at = (int) (std::find(faces.begin(), faces.end(), pathFace[r]) - faces.begin());
        // the strand's yarn order: does the face after pathFace follow it through its W+ edge?
        const int next = faces[(at + 1) % n];
        const int ex = exitWaleCorner(s.role[pathFace[r]]);
        const int a = m.vertexAt(pathFace[r], ex), b = m.vertexAt(pathFace[r], (ex + 1) % m.sides(pathFace[r]));
        bool forward = false;
        for (int k = 0; k < m.sides(next); ++k)
            if ((m.vertexAt(next, k) == b && m.vertexAt(next, (k + 1) % m.sides(next)) == a)) forward = true;
        for (int i = 0; i < n; ++i) {
            const int f = forward ? faces[(at + 1 + i) % n] : faces[((at - 1 - i) % n + n) % n];
            s.strip[f] = r;
            s.slot[f] = i;
            for (int k = 0; k < m.sides(f); ++k)
                if (s.role[f][k] == StitchEdgeRole::Bottom) s.rings[r].push_back(m.vertexAt(f, k));
            if (r == R - 1) {
                std::vector<int> tops;
                for (int k = 0; k < m.sides(f); ++k)
                    if (s.role[f][k] == StitchEdgeRole::Top) tops.push_back(m.vertexAt(f, (k + 1) % m.sides(f)));
                std::reverse(tops.begin(), tops.end());
                for (int v : tops) s.rings[R].push_back(v);
            }
        }
    }
    s.shiftPath = path;
    why = "constructed: " + std::to_string(R) + " closed rows, a connected chain of " + std::to_string(R) + " wale edges nearest azimuth 0";
    return true;
}

// Stitches crossing horizontal sections: wale edges whose ends straddle z = c, against 2 pi r(c) / w.
void writeSections(const std::string& out, const TypedStitchMesh& s, double w) {
    const PolygonMesh& m = s.mesh;
    std::set<std::pair<int, int>> wale;
    for (int f = 0; f < m.faces(); ++f)
        for (int k = 0; k < m.sides(f); ++k)
            if (s.role[f][k] == StitchEdgeRole::Wale) {
                const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % m.sides(f));
                wale.insert({std::min(a, b), std::max(a, b)});
            }
    const RevolutionProfile prof = tubeProfile();
    std::ofstream f(out + "/sections.csv", std::ios::binary);
    f << "z,wale_crossings,ideal_2pi_r_over_w\n";
    for (int zi = 1; zi < 150; ++zi) {
        const double z = zi + 0.5;
        int n = 0;
        for (const auto& e : wale) {
            const double za = m.pos[e.first].z, zb = m.pos[e.second].z;
            if ((za < z) != (zb < z)) ++n;
        }
        putf(f, "%.1f,%d,%.6f\n", z, n, 2 * kPi * prof.radius(z) / w);
    }
}

void bindAndProject(TypedStitchMesh& s, const TriMesh& m, double& maxMove) {
    const MeshQuery q(m);
    s.binding.clear();
    maxMove = 0;
    for (Vec3& p : s.mesh.pos) {
        const MeshQuery::Hit h = q.closest(p);
        VertexBinding b;
        b.tri = h.tri;
        if (h.tri >= 0) {
            const TrianglePoint tp = closestOnTriangle(p, m.pos[m.tri[h.tri][0]], m.pos[m.tri[h.tri][1]], m.pos[m.tri[h.tri][2]]);
            for (int k = 0; k < 3; ++k) b.bary[k] = tp.bary[k];
            maxMove = std::max(maxMove, norm(tp.p - p));
            p = tp.p;
            b.distance = 0;
        }
        s.binding.push_back(b);
    }
}

int commandBuild(const Args& a) {
    const std::string dir = a.get("--dir"), out = a.get("--out"), cs = a.get("--case", "D"), method = a.get("--method", "rows");
    const double w = a.number("--w", 5.5), h = a.number("--h", 5.5);
    if (dir.empty() || out.empty() || (cs != "D" && cs != "I") || (method != "rows" && method != "w18")) {
        std::fprintf(stderr, "build: --dir, --out, --case D|I and --method rows|w18 are required\n");
        return 2;
    }
    fs::create_directories(out);
    TriMesh m;
    std::vector<double> arc;
    std::string digestSaved, err;
    double T = 0;
    if (!readReference(dir, m, arc, digestSaved, T, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    std::vector<std::vector<int>> loops = boundaryLoops(m);
    std::sort(loops.begin(), loops.end(), [&](const std::vector<int>& x, const std::vector<int>& y) { return meanZ(m, x) < meanZ(m, y); });
    if (loops.size() != 2) {
        std::fprintf(stderr, "build: expected 2 boundary loops, found %zu\n", loops.size());
        return 1;
    }
    const std::vector<int>& castOn = cs == "D" ? loops.front() : loops.back();
    const std::vector<int>& bindOff = cs == "D" ? loops.back() : loops.front();
    std::vector<DeclaredBoundary> declared(2);
    declared[0].id = cs == "D" ? "z0" : "z150";
    declared[0].role = EdgeBoundary::CastOn;
    for (int v : castOn) declared[0].loop.push_back(m.pos[v]);
    declared[1].id = cs == "D" ? "z150" : "z0";
    declared[1].role = EdgeBoundary::BindOff;
    for (int v : bindOff) declared[1].loop.push_back(m.pos[v]);
    const double boundaryTol = 0.25 * h;   // the frozen boundary-correspondence bound
    std::ofstream rp(out + "/build_report.txt", std::ios::binary);
    putf(rp, "build: method %s (%s), case %s, gauge w %.4g h %.4g, reference digest %s (saved %s)\n", method.c_str(), method == "rows" ? "A, guided rows" : "B, W18 pipeline",
         cs.c_str(), w, h, hex64(meshDigest(m)).c_str(), digestSaved.c_str());
    const auto t0 = std::chrono::steady_clock::now();
    TypedStitchMesh s;
    GuidedRowsProduct gp;
    std::string shiftNote;
    double projected = 0;
    if (method == "rows") {
        std::vector<double> t;
        TimeFieldReport tr;
        if (!geodesicDistance(m, castOn, t, tr, &err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        GuidedRowsRequest req;
        req.w = w;
        req.h = h;
        if (!buildGuidedRows(m, t, castOn, bindOff, req, gp, &err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        putf(rp, "guided rows: %s; %s; contour events %d, |dn| limit binds %d\n", toString(gp.status), gp.detail.c_str(), gp.contourEvents, gp.limitBinds);
        for (const std::string& l : gp.log) rp << "  " << l << "\n";
        if (gp.status == Status::Fail) {
            putf(rp, "no mesh: component verdict FAIL (contour events shown above)\n");
            return 0;
        }
        s = gp.mesh;
        std::ofstream sf(out + "/shaping.csv", std::ios::binary);
        sf << "strip,slot,face,kind,events,index,slots_in_strip,sigma,target_fraction,target_slot,moved,placed,z,phi_deg,ring_below_n,ring_above_n,reason\n";
        for (const ShapingPlacement& p : gp.shaping) {
            const Vec3 c = s.mesh.faceCentroid(p.face);
            std::string reason = p.reason;
            std::replace(reason.begin(), reason.end(), ',', ';');
            putf(sf, "%d,%d,%d,%s,%d,%d,%d,%.6f,%.6f,%d,%d,%d,%.4f,%.3f,%d,%d,%s\n", p.strip, p.slot, p.face, kindGroup(p.kind), p.events, p.index, p.slotsInStrip, p.sigma,
                 p.targetFraction, p.targetSlot, p.moved, p.placed ? 1 : 0, c.z, std::atan2(c.y, c.x) * 180 / kPi, gp.rings[p.strip].n, gp.rings[p.strip + 1].n, reason.c_str());
        }
        // production's program (src/knitting/stitch/stitch_program.cpp buildProgram) on the same counts: where it puts
        // each shaping stitch in its turn (turn k + 1 consumes ring k's loops, as Method A's strip k does)
        {
            std::vector<int> counts;
            for (const RingSample& r : gp.rings) counts.push_back(r.n);
            const StitchProgram prog = buildProgram(counts);
            const GraphCheck gc = checkProgram(prog);
            std::ofstream pf(out + "/production_shaping.csv", std::ios::binary);
            pf << "strip,turn,op,index,turn_size,fraction,beside_path\n";
            for (const StitchRecord& st : prog.stitches) {
                if (st.op == StitchOp::K || (st.op == StitchOp::Inc && st.incSlot == 1)) continue;
                const int n = prog.turns[st.turn].output;
                const int index = st.id - prog.turnStart[st.turn];
                putf(pf, "%d,%d,%s,%d,%d,%.6f,%d\n", st.turn - 1, st.turn, st.op == StitchOp::Dec ? "DEC" : "INC", index, n, (index + 0.5) / n,
                     (index == 0 || index >= n - (st.op == StitchOp::Inc ? 2 : 1)) ? 1 : 0);
            }
            putf(rp, "production buildProgram on the same counts: graph check %s (%s)\n", toString(gc.status), gc.detail.c_str());
        }
        std::ofstream rf(out + "/ring_counts.csv", std::ios::binary);
        rf << "k,level,length,n,boundary,cut_parameter\n";
        for (const RingSample& r : gp.rings) putf(rf, "%d,%.9g,%.9g,%d,%d,%.9g\n", r.k, r.level, r.length, r.n, r.boundary ? 1 : 0, r.cutParameter);
    } else {
        ReferenceSurface surf;
        if (!prepareReferenceSurface(m, surf, &err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        RemeshRequest rq;
        rq.w = w;
        rq.h = h;
        rq.seed = 1;
        rq.maxAttempts = 3;
        rq.attemptCapS = 120;
        rq.axisGuide = globalAxisGuide("ring-z", m);
        RemeshProduct rmp;
        if (!remeshStage(surf, rq, rmp, &err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        rp << "--- remesh\n" << rmp.report.text() << rmp.kept.structural.text() << rmp.kept.quality.text();
        if (rmp.kept.mesh.faces() == 0) {
            rp << "no coarse mesh: component verdict FAIL\n";
            return 0;
        }
        {
            std::ofstream cf(out + "/coarse.obj", std::ios::binary);
            for (const Vec3& p : rmp.kept.mesh.pos) putf(cf, "v %.17g %.17g %.17g\n", p.x, p.y, p.z);
            for (int f = 0; f < rmp.kept.mesh.faces(); ++f) {
                cf << "f";
                for (int k = 0; k < rmp.kept.mesh.sides(f); ++k) cf << " " << rmp.kept.mesh.vertexAt(f, k) + 1;
                cf << "\n";
            }
        }
        LabelOptions lo;
        lo.up = cs == "D" ? Vec3{0, 0, 1} : Vec3{0, 0, -1};
        LabelProduct lp;
        if (!labelEdges(rmp.kept.mesh, lo, lp, &err)) {
            rp << "labels refused: " << err << "\ncomponent verdict FAIL\n";
            std::printf("build w18 %s: labels refused: %s\n", cs.c_str(), err.c_str());
            return 0;
        }
        rp << "--- labels\n" << lp.report.text();
        DirectionOptions dopt;
        dopt.up = lo.up;
        DirectionProduct dp;
        if (!assignDirections(lp.mesh, lp.labels, dopt, dp, &err)) {
            rp << "directions refused: " << err << "\ncomponent verdict FAIL\n";
            std::printf("build w18 %s: directions refused: %s\n", cs.c_str(), err.c_str());
            return 0;
        }
        rp << "--- directions\n" << dp.report.text();
        StitchMeshProduct sp;
        if (!buildStitchMesh(lp.mesh, dp.directed.top, sp, &err)) {
            rp << "subdivision refused: " << err << "\ncomponent verdict FAIL\n";
            return 0;
        }
        rp << "--- subdivision (ported validator)\n" << sp.report.text() << sp.check.text();
        {
            std::ofstream lf(out + "/coarse_labels.csv", std::ios::binary);
            lf << "a,b,role\n";
            for (const auto& [e, r] : lp.labels.role) putf(lf, "%d,%d,%s\n", e.first, e.second, r == EdgeRole::Course ? "course" : "wale");
            std::ofstream lm(out + "/coarse_labelled.obj", std::ios::binary);
            for (const Vec3& p : lp.mesh.pos) putf(lm, "v %.17g %.17g %.17g\n", p.x, p.y, p.z);
            for (int f = 0; f < lp.mesh.faces(); ++f) {
                lm << "f";
                for (int k = 0; k < lp.mesh.sides(f); ++k) lm << " " << lp.mesh.vertexAt(f, k) + 1;
                lm << "\n";
            }
        }
        {
            std::ofstream of(out + "/coarse_of_face.csv", std::ios::binary);
            of << "face,coarse_face,coarse_sides\n";
            for (int f = 0; f < sp.stitch.mesh.faces(); ++f) {
                const int cf = sp.stitch.coarseFace[f];
                putf(of, "%d,%d,%d\n", f, cf, cf >= 0 ? lp.mesh.sides(cf) : 0);
            }
        }
        s = typedFromStitchMesh(sp.stitch, declared, boundaryTol);
        bindAndProject(s, m, projected);
        putf(rp, "stitch vertices projected onto the reference (closest point): largest move %.4f mm\n", projected);
    }
    const double buildMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    TypedCheckOptions o;
    o.w = w;
    o.h = h;
    o.boundaryTolerance = boundaryTol;
    o.bindingTolerance = 0.001 * w;
    o.faceDeviation = 0.1 * w;
    o.requireStrips = method == "rows";
    TypedMeshCheck c = checkTypedMesh(s, m, declared, o);
    VirtualShift vs;
    VirtualShiftCheck vc;
    bool shift = false;
    if (method == "rows") {
        shift = buildVirtualShift(s, ShiftReading::PathSideDown, vs, &err);
    } else {
        std::string why;
        if (methodBShiftPath(s, c.strands, cs == "D", why)) {
            shift = buildVirtualShift(s, ShiftReading::PathSideDown, vs, &err);
            c = checkTypedMesh(s, m, declared, o);   // now with strips: the invariant per row is checked too
        }
        shiftNote = why;
        putf(rp, "method B shift-path: %s\n", why.c_str());
    }
    if (shift) vc = validateVirtualShift(s, vs);
    CoverageOptions co;
    co.coverDistance = 0.1 * w;
    co.boundaryBand = 0.5 * h;
    co.hausdorffLimit = 0.25 * h;
    const auto t1 = std::chrono::steady_clock::now();
    const CoverageReport cov = measureCoverage(m, {castOn, bindOff}, s.mesh, co);
    const double covMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
    const std::vector<CellMeasure> cells = measureCells(s, w, h);
    const CellSummary csum = summarize(cells);
    writeTypedMesh(out, s, c, shift ? &vs : nullptr, shift ? &vc : nullptr);
    writeCells(out, cells);
    writeSections(out, s, w);
    {
        std::ofstream f(out + "/coverage.csv", std::ios::binary);
        f << "i,vertex,x,y,z,weight,excluded,distance,face,crossings\n";
        for (size_t i = 0; i < cov.perSample.size(); ++i) {
            const CoverageSample& smp = cov.perSample[i];
            putf(f, "%zu,%d,%.6g,%.6g,%.6g,%.6g,%d,%.6g,%d,%d\n", i, smp.vertex ? 1 : 0, smp.p.x, smp.p.y, smp.p.z, smp.weight, smp.excluded ? 1 : 0,
                 smp.distance > 1e299 ? -1.0 : smp.distance, smp.face, smp.crossings);
        }
        std::ofstream lf(out + "/boundary_loops.csv", std::ios::binary);
        lf << "output_loop,reference_loop,output_vertices,output_to_reference,reference_to_output,hausdorff\n";
        for (const LoopCorrespondence& l : cov.loops)
            putf(lf, "%d,%d,%d,%.6f,%.6f,%.6f\n", l.outputLoop, l.referenceLoop, l.outputVertices, l.outputToReference, l.referenceToOutput, l.hausdorff);
    }
    const uint64_t digest = typedDigest(s);
    const Status verdict = componentVerdict(c, &cov);
    rp << "--- verdicts (component tube)\n";
    writeVerdicts(rp, c, &cov, shift ? &vc : nullptr, shiftNote, csum);
    putf(rp, "typed digest %s; build %.0f ms, coverage %.0f ms\n", hex64(digest).c_str(), buildMs, covMs);
    putf(rp, "component verdict (connectivity + coverage, no unresolved face): %s; geometry %s (reported apart); virtual shift %s\n", toString(verdict),
         toString(c.geometry), shift ? toString(vc.status) : "not built");
    {
        std::ofstream j(out + "/report.json", std::ios::binary);
        putf(j, "{\n  \"method\": \"%s\", \"case\": \"%s\", \"w\": %.6g, \"h\": %.6g,\n", method.c_str(), cs.c_str(), w, h);
        putf(j, "  \"typed_digest\": \"%s\",\n  \"reference_digest\": \"%s\",\n", hex64(digest).c_str(), hex64(meshDigest(m)).c_str());
        putf(j, "  \"verdict\": \"%s\", \"connectivity\": \"%s\", \"geometry\": \"%s\", \"coverage\": \"%s\", \"virtual_shift\": \"%s\",\n", toString(verdict),
             toString(c.connectivity), toString(c.geometry), toString(cov.status), shift ? toString(vc.status) : "not built");
        putf(j, "  \"faces\": %d, \"K\": %d, \"INC\": %d, \"DEC\": %d, \"special\": %d, \"unsupported\": %d,\n", c.faces, c.regular, c.increase, c.decrease, c.special,
             c.unsupported);
        putf(j, "  \"interlocks\": %d, \"both_top\": %d, \"both_bottom\": %d, \"course_against_wale\": %d, \"wale_disagree\": %d, \"cast_on\": %d, \"bind_off\": %d, "
                "\"border_course_without_role\": %d, \"border_wale\": %d,\n",
             c.interlocks, c.bothTop, c.bothBottom, c.courseAgainstWale, c.waleDisagree, c.castOnEdges, c.bindOffEdges, c.borderCourseWithoutRole, c.borderWaleEdges);
        putf(j, "  \"strands\": %zu, \"closed_strands\": %d, \"open_strands\": %d, \"yarn_ends\": %d, \"strips\": %zu, \"strips_broken\": %d,\n", c.strands.strands.size(),
             c.closedStrands, c.openStrands, c.strands.yarnEnds, c.strips.size(), c.stripsBroken);
        putf(j, "  \"inverted\": %d, \"degenerate\": %d, \"intersecting_pairs\": %d, \"max_binding_distance\": %.6g, \"max_face_deviation\": %.6g,\n", c.invertedFaces,
             c.degenerateFaces, c.intersectingPairs, c.maxBindingDistance, c.maxFaceDeviation);
        putf(j, "  \"uncovered_fraction\": %.6g, \"covered_twice\": %d, \"projected_area_ratio\": %.6g, \"worst_hausdorff\": %.6g, \"output_loops\": %d,\n",
             cov.uncoveredFraction, cov.coveredTwice, cov.projectedRatio, cov.worstHausdorff, cov.outputLoops);
        putf(j, "  \"cells_within\": %d, \"cells\": %d,\n", csum.within, csum.faces);
        putf(j, "  \"shift_note\": \"%s\", \"virtual_strands\": %d, \"virtual_longest\": %d,\n", shiftNote.c_str(), shift ? vc.strands : -1, shift ? vc.longestStrand : -1);
        putf(j, "  \"rows\": %d, \"h_prime\": %.6g, \"projected_move\": %.6g, \"build_ms\": %.1f, \"coverage_ms\": %.1f\n}\n", method == "rows" ? gp.rows : (int) c.strands.strands.size(),
             gp.hPrime, projected, buildMs, covMs);
    }
    std::printf("build %s %s w%.3g h%.3g: verdict %s (connectivity %s, coverage %s), geometry %s, shift %s; faces %d (K %d INC %d DEC %d SQ %d U %d); %.0f ms\n",
                method.c_str(), cs.c_str(), w, h, toString(verdict), toString(c.connectivity), toString(cov.status), toString(c.geometry), shift ? toString(vc.status) : "-",
                c.faces, c.regular, c.increase, c.decrease, c.special, c.unsupported, buildMs);
    return 0;
}


// ---------------------------------------------------------------------------------------------------------------------
// Checkpoint 3: relaxation
// ---------------------------------------------------------------------------------------------------------------------

// The cell statistics after relaxation: |width ratio - 1| and |height ratio - 1| of the K cells, pooled;
// the median over the cylinder regions (centroid z < 60 or z > 120) and the p95 over all K cells.
struct PostBounds {
    int kCells = 0, cylinderCells = 0, cellsWithin = 0, cellsMeasured = 0;
    double medianCylinder = 0, p95All = 0;
    bool met = false;      // median <= 0.05 and p95 <= 0.15
    bool allCells = false; // every measured cell within its frozen bounds
};

PostBounds postBounds(const std::vector<CellMeasure>& cells, const TypedStitchMesh& s, bool cylinderZones) {
    PostBounds b;
    std::vector<double> cyl, all;
    for (const CellMeasure& c : cells) {
        if (c.face < (int) s.reason.size() && s.reason[c.face] == "helix-start-step") continue;   // the virtual start triangle
        ++b.cellsMeasured;
        b.cellsWithin += c.withinBounds ? 1 : 0;
        if (c.kind != StitchFaceKind::Regular) continue;
        ++b.kCells;
        const double ew = std::fabs(c.widthRatio - 1), eh = std::fabs(c.heightRatio - 1);
        all.push_back(ew);
        all.push_back(eh);
        const double z = s.mesh.faceCentroid(c.face).z;
        if (!cylinderZones || z < 60 || z > 120) {
            ++b.cylinderCells;
            cyl.push_back(ew);
            cyl.push_back(eh);
        }
    }
    b.medianCylinder = quantile(cyl, 0.5);
    b.p95All = quantile(all, 0.95);
    b.met = b.medianCylinder <= 0.05 && b.p95All <= 0.15;
    b.allCells = b.cellsWithin == b.cellsMeasured;
    return b;
}

void writeObjPolygons(const std::string& path, const TypedStitchMesh& s) {
    std::ofstream f(path, std::ios::binary);
    f << "# stitchmesh_lab, mm\n";
    for (const Vec3& p : s.mesh.pos) putf(f, "v %.17g %.17g %.17g\n", p.x, p.y, p.z);
    const char* group = "";
    for (int fc = 0; fc < s.mesh.faces(); ++fc) {
        if (std::strcmp(group, kindGroup(s.kind[fc])) != 0) {
            group = kindGroup(s.kind[fc]);
            f << "g " << group << "\n";
        }
        f << "f";
        for (int k = 0; k < s.mesh.sides(fc); ++k) f << " " << s.mesh.vertexAt(fc, k) + 1;
        f << "\n";
    }
}

void writeCellsTo(const std::string& path, const std::vector<CellMeasure>& cells) {
    std::ofstream f(path, std::ios::binary);
    f << "face,kind,width_ratio,height_ratio,min_course,max_course,min_wale,max_wale,min_angle,max_angle,course_course_angle,area_over_wh,within,outside\n";
    for (const CellMeasure& c : cells)
        putf(f, "%d,%s,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%d,%s\n", c.face, kindGroup(c.kind), c.widthRatio, c.heightRatio, c.minCourse, c.maxCourse, c.minWale,
             c.maxWale, c.minAngle, c.maxAngle, c.courseCourseAngle, c.area, c.withinBounds ? 1 : 0, c.outside.c_str());
}

void writeHistory(const std::string& path, const RelaxProduct& p) {
    std::ofstream f(path, std::ios::binary);
    f << "step,rms_force,max_force,max_move,cg_iterations\n";
    for (const RelaxStep& s : p.history) putf(f, "%d,%.9g,%.9g,%.9g,%d\n", s.step, s.rmsForce, s.maxForce, s.maxMove, s.cgIterations);
}

std::string relaxText(const RelaxProduct& p) {
    char b[700];
    std::snprintf(b, sizeof b,
                  "relaxation %s: %s; elements: course %d, wale %d, diagonal %d springs, %d shear corners in %d sub-quads (%d faces without), %d wale struts "
                  "(%d junctions without); %d border vertices (%d snapped, largest snap %.4f mm); walk clamps %d; conflicting edge roles %d\n",
                  toString(p.status), p.detail.c_str(), p.courseSprings, p.waleSprings, p.diagonalSprings, p.shearCorners, p.subQuads, p.facesWithoutSubQuads, p.struts,
                  p.waleJunctionsWithoutStrut, p.borderVertices, p.borderVerticesSnapped, p.maxBorderSnap, p.walkClamps, p.conflictingEdgeRoles);
    return b;
}

std::vector<char> borderVerticesOf(const PolygonMesh& m) {
    std::vector<char> border(m.pos.size(), 0);
    for (const std::vector<int>& loop : polygonBorderLoops(m))
        for (int v : loop) border[v] = 1;
    return border;
}

MeshRelaxRequest relaxRequest(const Args& a, double w, double h) {
    MeshRelaxRequest r;
    r.w = w;
    r.h = h;
    r.timestep = a.number("--timestep", r.timestep);
    r.maxSteps = (int) a.number("--max-steps", r.maxSteps);
    return r;
}

std::string requestText(const MeshRelaxRequest& r) {
    char b[400];
    std::snprintf(b, sizeof b,
                  "w %.4g h %.4g k_stretch %.3g k_shear %.3g k_wale %.3g timestep %.4g max_step %.4g max_steps %d cap_s %.4g force_tolerance %.3g stall %d x %.3g "
                  "cg %d / %.1e",
                  r.w, r.h, r.kStretch, r.kShear, r.kWale, r.timestep, r.maxStep, r.maxSteps, r.capSeconds, r.forceTolerance, r.stallWindow, r.stallRatio,
                  r.cgMaxIterations, r.cgTolerance);
    return b;
}

// One relaxation with its before/after measurements on the stored faces and, when given, the virtual faces.
struct RelaxCaseResult {
    RelaxProduct relax;
    TypedMeshCheck storedBefore, storedAfter, virtualBefore, virtualAfter;
    CoverageReport coverBefore, coverAfter, coverStoredBefore, coverStoredAfter;   // cover*: solver faces; coverStored*: the stored tiling
    PostBounds boundsStoredBefore, boundsStoredAfter, boundsVirtualBefore, boundsVirtualAfter;
    uint64_t topologyBefore = 0, topologyAfter = 0;
    bool hasVirtual = false;
    TypedStitchMesh storedRelaxed, virtualRelaxed;
};

uint64_t topologyDigest(const TypedStitchMesh& s) {
    TypedStitchMesh t = s;
    for (Vec3& p : t.mesh.pos) p = Vec3{0, 0, 0};
    return typedDigest(t);
}

RelaxCaseResult relaxCase(const TypedStitchMesh& stored, const VirtualShift* vs, const TriMesh& ref, const std::vector<std::vector<int>>& refLoops,
                          const std::vector<DeclaredBoundary>& declared, const TypedCheckOptions& o, const CoverageOptions& co, const MeshRelaxRequest& rq,
                          bool cylinderZones, bool measureCoverage_) {
    RelaxCaseResult res;
    res.hasVirtual = vs != nullptr;
    std::vector<RelaxFace> faces;
    TypedStitchMesh virt;
    if (vs) {
        for (int f = 0; f < stored.mesh.faces(); ++f) faces.push_back(RelaxFace{vs->corners[f], stored.role[f]});
        virt = virtualTypedMesh(stored, *vs);
    } else {
        faces = relaxFacesOf(stored.mesh, stored.role);
    }
    TypedCheckOptions ov = o;
    ov.requireStrips = false;
    res.storedBefore = checkTypedMesh(stored, ref, declared, o);
    res.boundsStoredBefore = postBounds(measureCells(stored, o.w, o.h), stored, cylinderZones);
    if (vs) {
        res.virtualBefore = checkTypedMesh(virt, ref, declared, ov);
        res.boundsVirtualBefore = postBounds(measureCells(virt, o.w, o.h), virt, cylinderZones);
    }
    if (measureCoverage_) {
        res.coverBefore = measureCoverage(ref, refLoops, vs ? virt.mesh : stored.mesh, co);
        res.coverStoredBefore = vs ? measureCoverage(ref, refLoops, stored.mesh, co) : res.coverBefore;
    }
    std::string err;
    relaxStitchMesh(ref, refLoops, faces, stored.mesh.pos, borderVerticesOf(stored.mesh), rq, res.relax, &err);
    res.storedRelaxed = stored;
    res.storedRelaxed.mesh.pos = res.relax.pos;
    res.storedRelaxed.binding = res.relax.binding;
    res.storedAfter = checkTypedMesh(res.storedRelaxed, ref, declared, o);
    res.boundsStoredAfter = postBounds(measureCells(res.storedRelaxed, o.w, o.h), res.storedRelaxed, cylinderZones);
    res.topologyBefore = topologyDigest(stored);
    res.topologyAfter = topologyDigest(res.storedRelaxed);
    if (vs) {
        res.virtualRelaxed = virtualTypedMesh(res.storedRelaxed, *vs);
        res.virtualAfter = checkTypedMesh(res.virtualRelaxed, ref, declared, ov);
        res.boundsVirtualAfter = postBounds(measureCells(res.virtualRelaxed, o.w, o.h), res.virtualRelaxed, cylinderZones);
    }
    if (measureCoverage_) {
        res.coverAfter = measureCoverage(ref, refLoops, vs ? res.virtualRelaxed.mesh : res.storedRelaxed.mesh, co);
        res.coverStoredAfter = vs ? measureCoverage(ref, refLoops, res.storedRelaxed.mesh, co) : res.coverAfter;
    }
    return res;
}

std::string boundsText(const char* what, const PostBounds& b) {
    char t[300];
    std::snprintf(t, sizeof t, "%s: cells within frozen bounds %d / %d; K cells %d (cylinder %d): median |ratio - 1| cylinder %.4f (<= 0.05), p95 all %.4f (<= 0.15) -> %s\n",
                  what, b.cellsWithin, b.cellsMeasured, b.kCells, b.cylinderCells, b.medianCylinder, b.p95All, b.met ? "met" : "NOT met");
    return t;
}

// The CP3 component verdict: the shared criteria after relaxation. Cells and geometry on the solver's faces (Method A's
// virtual polygons); coverage on the stored tiling, as in CP2. Method A's virtual faces have the helix's end step on the
// bind-off border by construction, one row high, which the frozen boundary bound (0.25 h) does not admit: their coverage
// is reported on its own line, not hidden and not folded into this verdict.
Status relaxVerdict(const RelaxCaseResult& r) {
    const TypedMeshCheck& g = r.hasVirtual ? r.virtualAfter : r.storedAfter;
    const PostBounds& b = r.hasVirtual ? r.boundsVirtualAfter : r.boundsStoredAfter;
    if (r.relax.status == Status::Fail || r.storedAfter.connectivity == Status::Fail || g.geometry == Status::Fail || r.coverStoredAfter.status == Status::Fail ||
        r.topologyBefore != r.topologyAfter)
        return Status::Fail;
    if (r.storedAfter.connectivity != Status::Pass || !b.met || !b.allCells || g.geometry != Status::Pass || r.coverStoredAfter.status != Status::Pass ||
        r.relax.status != Status::Pass)
        return Status::Unresolved;
    return Status::Pass;
}

void writeRelaxReport(std::ostream& o, const RelaxCaseResult& r, const MeshRelaxRequest& rq) {
    o << "settings: " << requestText(rq) << "\n";
    o << relaxText(r.relax);
    putf(o, "solver faces: %s\n", r.hasVirtual ? "Method A's virtual polygons (shift-path read one ring lower at slot 0)" : "the stored faces (closed rows and strands as built)");
    putf(o, "topology digest before %s, after %s: %s\n", hex64(r.topologyBefore).c_str(), hex64(r.topologyAfter).c_str(), r.topologyBefore == r.topologyAfter ? "unchanged" : "CHANGED");
    o << "--- stored faces before\n" << r.storedBefore.text() << boundsText("stored before", r.boundsStoredBefore);
    o << "--- stored faces after\n" << r.storedAfter.text() << boundsText("stored after", r.boundsStoredAfter);
    if (r.hasVirtual) {
        o << "--- virtual faces before (diagnostic connectivity: the start triangle is listed unsupported, the end step as a border wale edge)\n"
          << r.virtualBefore.text() << boundsText("virtual before", r.boundsVirtualBefore);
        o << "--- virtual faces after\n" << r.virtualAfter.text() << boundsText("virtual after", r.boundsVirtualAfter);
    }
    o << "--- coverage of the stored tiling before\n" << r.coverStoredBefore.text() << "--- coverage of the stored tiling after\n" << r.coverStoredAfter.text();
    if (r.hasVirtual)
        o << "--- coverage of the virtual faces before (the helix end step lies on the bind-off border)\n" << r.coverBefore.text()
          << "--- coverage of the virtual faces after\n" << r.coverAfter.text();
    putf(o, "CP3 verdict (converged; connectivity; geometry and every cell within bounds on the solver faces; median/p95 met; coverage of the stored tiling): %s\n",
         toString(relaxVerdict(r)));
}

// Loads a typed mesh written by build (the same files the reload reads).
bool readTypedFiles(const std::string& dir, TypedStitchMesh& s, bool& strips, std::string& savedDigest, std::string& cs) {
    s = TypedStitchMesh{};
    {
        std::ifstream f(dir + "/candidate.obj");
        if (!f) return false;
        std::string line;
        while (std::getline(f, line)) {
            if (line.size() > 2 && line[0] == 'v' && line[1] == ' ') {
                Vec3 p;
                char* e = nullptr;
                p.x = std::strtod(line.c_str() + 2, &e);
                p.y = std::strtod(e, &e);
                p.z = std::strtod(e, &e);
                s.mesh.pos.push_back(p);
            } else if (line.size() > 2 && line[0] == 'f' && line[1] == ' ') {
                std::istringstream is(line.substr(2));
                std::vector<int> c;
                int v;
                while (is >> v) c.push_back(v - 1);
                s.mesh.addFace(c);
            }
        }
    }
    const int nF = s.mesh.faces();
    s.role.assign(nF, {});
    s.boundary.assign(nF, {});
    s.kind.assign(nF, StitchFaceKind::Unsupported);
    s.reason.assign(nF, "");
    {
        std::ifstream f(dir + "/connectivity.csv");
        std::string line;
        std::getline(f, line);
        while (std::getline(f, line)) {
            std::vector<std::string> col;
            std::stringstream ss(line);
            std::string x;
            while (std::getline(ss, x, ',')) col.push_back(x);
            if (col.size() < 6) continue;
            const int fc = std::stoi(col[0]);
            s.role[fc].push_back(col[4] == "W" ? StitchEdgeRole::Wale : col[4] == "B" ? StitchEdgeRole::Bottom : StitchEdgeRole::Top);
            s.boundary[fc].push_back(col[5] == "cast-on" ? EdgeBoundary::CastOn : col[5] == "bind-off" ? EdgeBoundary::BindOff : col[5] == "side" ? EdgeBoundary::Side
                                                                                                                                                   : EdgeBoundary::Interior);
        }
    }
    strips = false;
    {
        std::ifstream f(dir + "/face_types.csv");
        std::string line;
        std::getline(f, line);
        std::vector<int> strip(nF, -1), slot(nF, -1);
        while (std::getline(f, line)) {
            std::vector<std::string> col;
            std::stringstream ss(line);
            std::string x;
            while (std::getline(ss, x, ',')) col.push_back(x);
            if (col.size() < 11) continue;
            const int fc = std::stoi(col[0]);
            const std::string& k = col[2];
            s.kind[fc] = k == "K" ? StitchFaceKind::Regular : k == "INC_kp" ? StitchFaceKind::Increase : k == "DEC_d12k" ? StitchFaceKind::Decrease
                       : k == "special_quad" ? StitchFaceKind::SpecialQuad : k == "pole_start" ? StitchFaceKind::PoleStart : StitchFaceKind::Unsupported;
            strip[fc] = std::stoi(col[4]);
            slot[fc] = std::stoi(col[5]);
            s.reason[fc] = col.size() > 11 ? col[11] : "";
            if (strip[fc] >= 0) strips = true;
        }
        if (strips) {
            s.strip = strip;
            s.slot = slot;
        }
    }
    {
        std::ifstream f(dir + "/ring_vertices.csv");
        std::string line;
        std::getline(f, line);
        while (std::getline(f, line)) {
            int k = 0, j = 0, v = 0;
            if (!threeInts(line, k, j, v)) continue;
            if ((int) s.rings.size() <= k) s.rings.resize(k + 1);
            s.rings[k].push_back(v);
        }
        std::ifstream g(dir + "/shift_path.csv");
        std::getline(g, line);
        while (std::getline(g, line)) {
            ShiftPathEdge e;
            if (threeInts(line, e.strip, e.bottom, e.top)) s.shiftPath.push_back(e);
        }
    }
    std::ifstream jf(dir + "/report.json");
    const std::string json((std::istreambuf_iterator<char>(jf)), std::istreambuf_iterator<char>());
    const size_t d = json.find("\"typed_digest\": \"");
    savedDigest = d == std::string::npos ? "" : json.substr(d + 17, 16);
    const size_t cpos = json.find("\"case\": \"");
    cs = cpos == std::string::npos ? "D" : json.substr(cpos + 9, 1);
    return true;
}

std::vector<DeclaredBoundary> declaredOf(const TriMesh& m, const std::string& cs, std::vector<int>& castOn, std::vector<int>& bindOff) {
    std::vector<std::vector<int>> loops = boundaryLoops(m);
    std::sort(loops.begin(), loops.end(), [&](const std::vector<int>& x, const std::vector<int>& y) { return meanZ(m, x) < meanZ(m, y); });
    castOn = cs == "D" ? loops.front() : loops.back();
    bindOff = cs == "D" ? loops.back() : loops.front();
    std::vector<DeclaredBoundary> declared(2);
    declared[0].role = EdgeBoundary::CastOn;
    declared[1].role = EdgeBoundary::BindOff;
    for (int v : castOn) declared[0].loop.push_back(m.pos[v]);
    for (int v : bindOff) declared[1].loop.push_back(m.pos[v]);
    return declared;
}

// stitchmesh_lab relax --case-dir <CP2 run> --dir <reference dir> --w --h --out <dir>
int commandRelax(const Args& a) {
    const std::string dir = a.get("--case-dir"), refDir = a.get("--dir"), out = a.get("--out");
    const double w = a.number("--w", 5.5), h = a.number("--h", 5.5);
    if (dir.empty() || refDir.empty() || out.empty()) {
        std::fprintf(stderr, "relax: --case-dir, --dir and --out are required\n");
        return 2;
    }
    fs::create_directories(out);
    TypedStitchMesh s;
    bool strips = false;
    std::string saved, cs, err, digestSaved;
    if (!readTypedFiles(dir, s, strips, saved, cs)) {
        std::fprintf(stderr, "relax: cannot read %s\n", dir.c_str());
        return 1;
    }
    TriMesh m;
    std::vector<double> arc;
    double T = 0;
    if (!readReference(refDir, m, arc, digestSaved, T, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    std::vector<int> castOn, bindOff;
    const std::vector<DeclaredBoundary> declared = declaredOf(m, cs, castOn, bindOff);
    // bindings of the CP2 positions (A: on the surface; B: projected in CP2)
    {
        const MeshQuery q(m);
        for (const Vec3& p : s.mesh.pos) {
            const MeshQuery::Hit hit = q.closest(p);
            VertexBinding b;
            b.tri = hit.tri;
            b.distance = hit.dist;
            s.binding.push_back(b);
        }
    }
    const std::string loaded = hex64(typedDigest(s));
    TypedCheckOptions o;
    o.w = w;
    o.h = h;
    o.boundaryTolerance = 0.25 * h;
    o.bindingTolerance = 0.001 * w;
    o.faceDeviation = 0.1 * w;
    o.requireStrips = strips;
    CoverageOptions co;
    co.coverDistance = 0.1 * w;
    co.boundaryBand = 0.5 * h;
    co.hausdorffLimit = 0.25 * h;
    VirtualShift vs;
    const bool useVirtual = strips && !s.shiftPath.empty() && buildVirtualShift(s, ShiftReading::PathSideDown, vs, &err) && validateVirtualShift(s, vs).status == Status::Pass;
    const MeshRelaxRequest rq = relaxRequest(a, w, h);
    const RelaxCaseResult r = relaxCase(s, useVirtual ? &vs : nullptr, m, {castOn, bindOff}, declared, o, co, rq, true, true);
    writeObjPolygons(out + "/relaxed.obj", r.storedRelaxed);
    writeCellsTo(out + "/cells_before.csv", measureCells(s, w, h));
    writeCellsTo(out + "/cells_after.csv", measureCells(r.storedRelaxed, w, h));
    if (useVirtual) {
        writeObjPolygons(out + "/virtual_before.obj", virtualTypedMesh(s, vs));
        writeObjPolygons(out + "/virtual_after.obj", r.virtualRelaxed);
        writeCellsTo(out + "/cells_virtual_before.csv", measureCells(virtualTypedMesh(s, vs), w, h));
        writeCellsTo(out + "/cells_virtual_after.csv", measureCells(r.virtualRelaxed, w, h));
        std::ofstream vf(out + "/virtual_faces.csv", std::ios::binary);
        vf << "face,strip,slot,re_read,kind,reason,virtual_corners\n";
        const TypedStitchMesh vt = virtualTypedMesh(s, vs);
        for (int f = 0; f < vt.mesh.faces(); ++f) {
            putf(vf, "%d,%d,%d,%d,%s,%s,", f, s.strip[f], s.slot[f], (int) vs.reRead[f], kindGroup(vt.kind[f]), vt.reason[f].c_str());
            for (int k = 0; k < vt.mesh.sides(f); ++k) vf << (k ? " " : "") << vt.mesh.vertexAt(f, k);
            vf << "\n";
        }
    }
    writeHistory(out + "/relax_history.csv", r.relax);
    {
        std::ofstream bf(out + "/relaxed_bindings.csv", std::ios::binary);
        bf << "v,tri,b0,b1,b2\n";
        for (size_t v = 0; v < r.relax.binding.size(); ++v) {
            const VertexBinding& b = r.relax.binding[v];
            putf(bf, "%zu,%d,%.17g,%.17g,%.17g\n", v, b.tri, b.bary[0], b.bary[1], b.bary[2]);
        }
        std::ofstream cf(out + "/coverage_after.csv", std::ios::binary);
        cf << "i,vertex,x,y,z,weight,excluded,distance,face,crossings\n";
        for (size_t i = 0; i < r.coverStoredAfter.perSample.size(); ++i) {
            const CoverageSample& smp = r.coverStoredAfter.perSample[i];
            putf(cf, "%zu,%d,%.6g,%.6g,%.6g,%.6g,%d,%.6g,%d,%d\n", i, smp.vertex ? 1 : 0, smp.p.x, smp.p.y, smp.p.z, smp.weight, smp.excluded ? 1 : 0,
                 smp.distance > 1e299 ? -1.0 : smp.distance, smp.face, smp.crossings);
        }
    }
    const Status verdict = relaxVerdict(r);
    {
        std::ofstream rp(out + "/relax_report.txt", std::ios::binary);
        putf(rp, "relax: %s, case %s, gauge w %.4g h %.4g; CP2 typed digest saved %s, loaded %s (%s)\n", dir.c_str(), cs.c_str(), w, h, saved.c_str(), loaded.c_str(),
             saved == loaded ? "equal" : "DIFFERENT");
        writeRelaxReport(rp, r, rq);
    }
    {
        const TypedMeshCheck& g = useVirtual ? r.virtualAfter : r.storedAfter;
        const TypedMeshCheck& g0 = useVirtual ? r.virtualBefore : r.storedBefore;
        const PostBounds& b0 = useVirtual ? r.boundsVirtualBefore : r.boundsStoredBefore;
        const PostBounds& b1 = useVirtual ? r.boundsVirtualAfter : r.boundsStoredAfter;
        std::ofstream j(out + "/relax_report.json", std::ios::binary);
        putf(j, "{\n  \"case_dir\": \"%s\", \"case\": \"%s\", \"w\": %.6g, \"h\": %.6g, \"solver_faces\": \"%s\",\n", fs::path(dir).filename().string().c_str(), cs.c_str(), w,
             h, useVirtual ? "virtual" : "stored");
        putf(j, "  \"cp2_digest_equal\": %s, \"topology_unchanged\": %s,\n", saved == loaded ? "true" : "false", r.topologyBefore == r.topologyAfter ? "true" : "false");
        putf(j, "  \"verdict\": \"%s\", \"stop\": \"%s\", \"steps\": %d, \"ms\": %.1f, \"rms_start\": %.6g, \"rms_end\": %.6g,\n", toString(verdict), r.relax.stopReason.c_str(),
             r.relax.steps, r.relax.ms, r.relax.rmsForceStart, r.relax.rmsForceEnd);
        putf(j, "  \"connectivity_after\": \"%s\", \"geometry_before\": \"%s\", \"geometry_after\": \"%s\", \"stored_geometry_after\": \"%s\",\n", toString(r.storedAfter.connectivity),
             toString(g0.geometry), toString(g.geometry), toString(r.storedAfter.geometry));
        putf(j, "  \"inverted_after\": %d, \"intersecting_after\": %d, \"stored_inverted_after\": %d, \"stored_intersecting_after\": %d,\n", g.invertedFaces, g.intersectingPairs,
             r.storedAfter.invertedFaces, r.storedAfter.intersectingPairs);
        putf(j, "  \"coverage_before\": \"%s\", \"coverage_after\": \"%s\", \"uncovered_after\": %.6g, \"twice_after\": %d, \"area_ratio_after\": %.6g, \"hausdorff_after\": %.6g,\n",
             toString(r.coverStoredBefore.status), toString(r.coverStoredAfter.status), r.coverStoredAfter.uncoveredFraction, r.coverStoredAfter.coveredTwice,
             r.coverStoredAfter.projectedRatio, r.coverStoredAfter.worstHausdorff);
        putf(j, "  \"virtual_coverage_before\": \"%s\", \"virtual_coverage_after\": \"%s\", \"virtual_hausdorff_before\": %.6g, \"virtual_hausdorff_after\": %.6g, "
                "\"virtual_uncovered_after\": %.6g, \"virtual_twice_after\": %d,\n",
             useVirtual ? toString(r.coverBefore.status) : "-", useVirtual ? toString(r.coverAfter.status) : "-", useVirtual ? r.coverBefore.worstHausdorff : 0.0,
             useVirtual ? r.coverAfter.worstHausdorff : 0.0, useVirtual ? r.coverAfter.uncoveredFraction : 0.0, useVirtual ? r.coverAfter.coveredTwice : 0);
        putf(j, "  \"cells_within_before\": %d, \"cells_within_after\": %d, \"cells_measured\": %d, \"median_cyl_before\": %.6g, \"median_cyl_after\": %.6g, "
                "\"p95_before\": %.6g, \"p95_after\": %.6g, \"bounds_met_after\": %s,\n",
             b0.cellsWithin, b1.cellsWithin, b1.cellsMeasured, b0.medianCylinder, b1.medianCylinder, b0.p95All, b1.p95All, b1.met ? "true" : "false");
        putf(j, "  \"stored_cells_within_before\": %d, \"stored_cells_within_after\": %d, \"stored_cells\": %d, \"stored_p95_after\": %.6g,\n", r.boundsStoredBefore.cellsWithin,
             r.boundsStoredAfter.cellsWithin, r.boundsStoredAfter.cellsMeasured, r.boundsStoredAfter.p95All);
        putf(j, "  \"struts\": %d, \"sub_quads\": %d, \"faces_without_sub_quads\": %d, \"border_snapped\": %d, \"max_border_snap\": %.6g, \"walk_clamps\": %d\n}\n", r.relax.struts,
             r.relax.subQuads, r.relax.facesWithoutSubQuads, r.relax.borderVerticesSnapped, r.relax.maxBorderSnap, r.relax.walkClamps);
    }
    std::printf("relax %s: %s, %s after %d steps (%.0f ms), RMS %.3g -> %.3g; geometry %s; coverage %s; cells within %d -> %d of %d; p95 %.3f -> %.3f\n",
                fs::path(dir).filename().string().c_str(), toString(verdict), r.relax.stopReason.c_str(), r.relax.steps, r.relax.ms, r.relax.rmsForceStart, r.relax.rmsForceEnd,
                toString((useVirtual ? r.virtualAfter : r.storedAfter).geometry), toString(r.coverStoredAfter.status),
                (useVirtual ? r.boundsVirtualBefore : r.boundsStoredBefore).cellsWithin, (useVirtual ? r.boundsVirtualAfter : r.boundsStoredAfter).cellsWithin,
                (useVirtual ? r.boundsVirtualAfter : r.boundsStoredAfter).cellsMeasured, (useVirtual ? r.boundsVirtualBefore : r.boundsStoredBefore).p95All,
                (useVirtual ? r.boundsVirtualAfter : r.boundsStoredAfter).p95All);
    return 0;
}

// The relaxation fixtures: the solver checked on hand-built tubes before any surface case.
int commandRelaxFixtures(const Args& a) {
    const std::string out = a.get("--out");
    if (out.empty()) {
        std::fprintf(stderr, "relax-fixtures: --out <dir> is required\n");
        return 2;
    }
    fs::create_directories(out);
    std::ofstream sum(out + "/relax_fixtures.csv", std::ios::binary);
    sum << "name,solver_faces,stop,steps,ms,rms_start,rms_end,cells_within_before,cells_within_after,cells,p95_before,p95_after,geometry_after,max_ring_z_spread,"
           "mean_rise_per_stitch,expected_rise\n";
    std::ofstream rp(out + "/relax_fixtures_report.txt", std::ios::binary);
    struct Spec {
        std::string name;
        std::vector<int> counts;
        std::vector<std::vector<int>> shaping;
        bool virtualFaces;
        double perturb;
        double timestep;
    };
    const double w = 5.5, h = 5.5;
    std::vector<Spec> specs = {
        {"tube16x6-rest-stored", std::vector<int>(7, 16), std::vector<std::vector<int>>(6), false, 0, 1.0},
        {"tube16x6-perturbed-stored", std::vector<int>(7, 16), std::vector<std::vector<int>>(6), false, 0.15, 1.0},
        {"tube16x6-virtual-helix", std::vector<int>(7, 16), std::vector<std::vector<int>>(6), true, 0, 1.0},
        {"kp-16-17-stored", {16, 16, 16, 17, 17, 17, 17}, {{}, {}, {8}, {}, {}, {}}, false, 0, 1.0},
        {"kp-16-17-virtual", {16, 16, 16, 17, 17, 17, 17}, {{}, {}, {8}, {}, {}, {}}, true, 0, 1.0},
        {"d12k-17-16-stored", {17, 17, 17, 16, 16, 16, 16}, {{}, {}, {8}, {}, {}, {}}, false, 0, 1.0},
        {"d12k-17-16-virtual", {17, 17, 17, 16, 16, 16, 16}, {{}, {}, {8}, {}, {}, {}}, true, 0, 1.0},
        {"tube16x6-stall-control", std::vector<int>(7, 16), std::vector<std::vector<int>>(6), false, 0.15, 1e-6},
    };
    for (const Spec& sp : specs) {
        const std::string dir = out + "/" + sp.name;
        fs::create_directories(dir);
        // the radius whose chord between 16 vertices is exactly w: ring 0 is at rest
        const double radius = w / (2 * std::sin(kPi / 16));
        TypedStitchMesh s;
        TriMesh unusedRef;
        std::vector<DeclaredBoundary> declared;
        std::string err;
        ringTubeFixture(sp.counts, sp.shaping, radius, h, s, unusedRef, declared, &err);
        const int R = (int) sp.counts.size() - 1;
        const TriMesh ref = cylinderReference(radius, 0, R * h, h / 8);
        TriMesh refN = ref;
        refN.computeNormals();
        const std::vector<std::vector<int>> loops = triangleBoundaryLoops(refN);
        if (sp.perturb > 0) {   // deterministic tangential perturbation of the interior vertices
            for (size_t v = 0; v < s.mesh.pos.size(); ++v) {
                const Vec3 p = s.mesh.pos[v];
                if (p.z < 1e-9 || p.z > R * h - 1e-9) continue;
                const double u = hash01(7, (int) v, 0, 0) - 0.5, t = hash01(7, (int) v, 1, 0) - 0.5;
                const double phi = std::atan2(p.y, p.x) + 2 * sp.perturb * w * u / radius;
                s.mesh.pos[v] = Vec3{radius * std::cos(phi), radius * std::sin(phi), p.z + 2 * sp.perturb * h * t};
            }
        }
        VirtualShift vs;
        if (sp.virtualFaces) buildVirtualShift(s, ShiftReading::PathSideDown, vs, &err);
        TypedCheckOptions o;
        o.w = w;
        o.h = h;
        o.boundaryTolerance = 0.01;
        o.bindingTolerance = 0.01;
        o.faceDeviation = 0.55;
        o.requireStrips = true;
        CoverageOptions co;
        MeshRelaxRequest rq;
        rq.w = w;
        rq.h = h;
        rq.timestep = sp.timestep;
        // fixture checks against the cylinder ending at the border rings
        std::vector<DeclaredBoundary> decl2 = declared;
        const RelaxCaseResult r = relaxCase(s, sp.virtualFaces ? &vs : nullptr, refN, loops, decl2, o, co, rq, false, false);
        writeObjPolygons(dir + "/before.obj", s);
        writeObjPolygons(dir + "/after.obj", r.storedRelaxed);
        if (sp.virtualFaces) {
            writeObjPolygons(dir + "/virtual_before.obj", virtualTypedMesh(s, vs));
            writeObjPolygons(dir + "/virtual_after.obj", r.virtualRelaxed);
        }
        writeHistory(dir + "/relax_history.csv", r.relax);
        // the helix: per interior ring, the mean z rise from vertex j to j + 1 for j = 1 .. n - 2 and the ring's z spread
        double spread = 0, rise = 0;
        int nr = 0;
        for (int k = 1; k < R; ++k) {
            const std::vector<int>& ring = s.rings[k];
            double lo = 1e300, hi = -1e300;
            for (int v : ring) {
                lo = std::min(lo, r.relax.pos[v].z);
                hi = std::max(hi, r.relax.pos[v].z);
            }
            spread = std::max(spread, hi - lo);
            for (size_t j = 1; j + 1 < ring.size(); ++j) {   // vertex 0 is the path vertex: virtually the next round's corner
                rise += r.relax.pos[ring[j + 1]].z - r.relax.pos[ring[j]].z;
                ++nr;
            }
        }
        rise /= std::max(1, nr);
        const PostBounds& b0 = sp.virtualFaces ? r.boundsVirtualBefore : r.boundsStoredBefore;
        const PostBounds& b1 = sp.virtualFaces ? r.boundsVirtualAfter : r.boundsStoredAfter;
        const TypedMeshCheck& g = sp.virtualFaces ? r.virtualAfter : r.storedAfter;
        putf(sum, "%s,%s,%s,%d,%.1f,%.6g,%.6g,%d,%d,%d,%.5f,%.5f,%s,%.5f,%.5f,%.5f\n", sp.name.c_str(), sp.virtualFaces ? "virtual" : "stored", r.relax.stopReason.c_str(),
             r.relax.steps, r.relax.ms, r.relax.rmsForceStart, r.relax.rmsForceEnd, b0.cellsWithin, b1.cellsWithin, b1.cellsMeasured, b0.p95All, b1.p95All,
             toString(g.geometry), spread, rise, h / 16);
        rp << "=== " << sp.name << "\n";
        writeRelaxReport(rp, r, rq);
        putf(rp, "interior rings: largest z spread %.4f mm, mean rise per stitch %.5f mm (one row per round would be h / n = %.5f)\n", spread, rise, h / 16);
        std::printf("%s: %s after %d steps, RMS %.3g -> %.3g, cells %d -> %d of %d, p95 %.4f -> %.4f, geometry %s, rise %.4f\n", sp.name.c_str(), r.relax.stopReason.c_str(),
                    r.relax.steps, r.relax.rmsForceStart, r.relax.rmsForceEnd, b0.cellsWithin, b1.cellsWithin, b1.cellsMeasured, b0.p95All, b1.p95All,
                    toString(g.geometry), rise);
    }
    return 0;
}


// ---------------------------------------------------------------------------------------------------------------------
// Checkpoint 4: the teddy's prepared components, analysed for Method A (no stitch mesh)
// ---------------------------------------------------------------------------------------------------------------------

struct DirectionAnalysis {
    std::string name;            // "A->B" or "B->A"
    int startVertex = -1, endVertex = -1;
    double T = 0, tMax = 0;      // field at the end pole; its largest value
    int rows = 0;
    double hPrime = 0;
    int levels = 0, singleClosed = 0, events = 0, firstEvent = -1;
    int extraMaxima = 0;         // local maxima of t farther than 2 h from the end pole
    int smallRings = 0;          // single closed levels with n < 6 (pole construction needed)
    int limitBinds = 0;          // |dn| > floor(n / 2)
    int maxN = 0, firstN = 0, lastN = 0;
    double hiddenFirst = -1, hiddenLast = -1;   // t range of the vertices inside another piece
    int pathVertices = 0;
    bool pathReachesStart = false;
    double pathLength = 0;
};

std::vector<std::vector<int>> vertexNeighbours(const TriMesh& m) {
    std::vector<std::set<int>> s(m.pos.size());
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) {
            s[t[k]].insert(t[(k + 1) % 3]);
            s[t[k]].insert(t[(k + 2) % 3]);
        }
    std::vector<std::vector<int>> out(m.pos.size());
    for (size_t v = 0; v < s.size(); ++v) out[v].assign(s[v].begin(), s[v].end());
    return out;
}

int nearestVertex(const TriMesh& m, const Vec3& p, double* dist) {
    int best = 0;
    double bd = 1e300;
    for (size_t v = 0; v < m.pos.size(); ++v) {
        const double d = norm(m.pos[v] - p);
        if (d < bd) {
            bd = d;
            best = (int) v;
        }
    }
    if (dist) *dist = bd;
    return best;
}

DirectionAnalysis analyseDirection(const TriMesh& m, const std::vector<std::vector<int>>& nb, const std::vector<char>& inside, int vs, int ve, const std::string& name,
                                   double w, double h, const std::string& dir, std::vector<double>& tOut) {
    DirectionAnalysis d;
    d.name = name;
    d.startVertex = vs;
    d.endVertex = ve;
    std::vector<double> t;
    TimeFieldReport tr;
    std::string err;
    geodesicDistance(m, {vs}, t, tr, &err);
    tOut = t;
    d.T = t[ve];
    for (double x : t) d.tMax = std::max(d.tMax, x);
    // local maxima of the field (a strict maximum over the one-ring) away from the end pole: a second summit
    std::ofstream mf(dir + "/maxima_" + (name == "A->B" ? "AB" : "BA") + ".csv", std::ios::binary);
    mf << "v,x,y,z,t,distance_to_end_pole\n";
    for (size_t v = 0; v < m.pos.size(); ++v) {
        bool top = true;
        for (int u : nb[v])
            if (t[u] >= t[v]) top = false;
        if (!top) continue;
        const double de = norm(m.pos[v] - m.pos[ve]);
        putf(mf, "%zu,%.6f,%.6f,%.6f,%.6f,%.6f\n", v, m.pos[v].x, m.pos[v].y, m.pos[v].z, t[v], de);
        if (de > 2 * h) ++d.extraMaxima;
    }
    for (size_t v = 0; v < m.pos.size(); ++v)
        if (inside[v]) {
            d.hiddenFirst = d.hiddenFirst < 0 ? t[v] : std::min(d.hiddenFirst, t[v]);
            d.hiddenLast = std::max(d.hiddenLast, t[v]);
        }
    // the levels t = k h' between the poles (the rows of Method A)
    d.rows = std::max(1, (int) std::lround(d.T / h));
    d.hPrime = d.T / d.rows;
    std::ofstream lf(dir + "/levels_" + (name == "A->B" ? "AB" : "BA") + ".csv", std::ios::binary);
    lf << "k,level,contours,closed,length_total,length_largest,n,dn,event,small_ring,limit_bind,hidden_fraction\n";
    int prevN = -1;
    for (int k = 1; k < d.rows; ++k) {
        const double level = k * d.hPrime;
        const RingLevel lv = extractLevel(m, t, level);
        ++d.levels;
        int closed = 0;
        double total = 0, largest = 0;
        for (const Contour& c : lv.contours) {
            closed += c.closed ? 1 : 0;
            total += c.length;
            largest = std::max(largest, c.length);
        }
        const bool single = lv.contours.size() == 1 && closed == 1;
        const int n = single ? (int) std::lround(total / w) : -1;
        int hiddenPts = 0, pts = 0;
        for (const Contour& c : lv.contours)
            for (const RingPoint& rp : c.points) {
                ++pts;
                if (inside[rp.a] && inside[rp.b]) ++hiddenPts;
            }
        bool event = !single, small = false, bind = false;
        int dn = 0;
        if (single) {
            ++d.singleClosed;
            if (d.firstN == 0) d.firstN = n;
            d.lastN = n;
            d.maxN = std::max(d.maxN, n);
            if (n < 6) {
                small = true;
                ++d.smallRings;
            }
            dn = prevN > 0 ? n - prevN : 0;
            if (prevN > 0 && std::abs(n - prevN) > prevN / 2) {
                bind = true;
                ++d.limitBinds;
            }
            prevN = n;
        } else {
            ++d.events;
            if (d.firstEvent < 0) d.firstEvent = k;
        }
        putf(lf, "%d,%.6f,%zu,%d,%.6f,%.6f,%d,%d,%d,%d,%d,%.4f\n", k, level, lv.contours.size(), closed, total, largest, n, dn, event ? 1 : 0,
             small ? 1 : 0, bind ? 1 : 0, pts ? (double) hiddenPts / pts : 0.0);
    }
    // the proposed cut guide: the discrete steepest descent of t from the end pole to the start pole along mesh edges
    std::ofstream pf(dir + "/cutpath_" + (name == "A->B" ? "AB" : "BA") + ".csv", std::ios::binary);
    pf << "i,v,x,y,z,t\n";
    int v = ve;
    for (int i = 0; i < (int) m.pos.size(); ++i) {
        putf(pf, "%d,%d,%.6f,%.6f,%.6f,%.6f\n", i, v, m.pos[v].x, m.pos[v].y, m.pos[v].z, t[v]);
        ++d.pathVertices;
        if (v == vs) {
            d.pathReachesStart = true;
            break;
        }
        int next = -1;
        for (int u : nb[v])
            if (t[u] < t[v] && (next < 0 || t[u] < t[next])) next = u;
        if (next < 0) break;
        d.pathLength += norm(m.pos[next] - m.pos[v]);
        v = next;
    }
    return d;
}

// stitchmesh_lab component --manifest <manifest.txt> --out <dir> [--gauge-piece <file> --stitches 65 --w 5.5 --h 3.4]
int commandComponent(const Args& a) {
    const std::string manPath = a.get("--manifest"), out = a.get("--out");
    const std::string gaugeFile = a.get("--gauge-piece", "piece_bear_002_00.obj");
    const double stitches = a.number("--stitches", 65), w = a.number("--w", 5.5), h = a.number("--h", 3.4);
    if (manPath.empty() || out.empty()) {
        std::fprintf(stderr, "component: --manifest and --out are required\n");
        return 2;
    }
    fs::create_directories(out);
    Manifest man;
    std::string err;
    if (!readManifest(manPath, man, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    Gauge g;
    if (!computeGauge(man, GaugeRequest{gaugeFile, stitches, w}, g, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const int np = (int) man.pieces.size();
    std::vector<LoadedPiece> pieces(np);
    std::vector<TriMesh> common(np);   // sigma x scene coordinates: one frame for every piece of the scene
    for (int i = 0; i < np; ++i) {
        if (!loadPiece(man, i, g.sigma, pieces[i], &err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        pieces[i].mesh.computeNormals();
        common[i] = pieces[i].mesh;
        for (Vec3& p : common[i].pos) p = p + pieces[i].map.centre * g.sigma;
        common[i].computeNormals();
    }
    std::ofstream sum(out + "/pieces.csv", std::ios::binary);
    sum << "piece,object,area_mm2,vertices,triangles,components,closed,euler,consistent,self_intersection,pole_A_hidden,pole_B_hidden,hidden_area_fraction,"
           "T_AB,rows_AB,events_AB,first_event_AB,extra_maxima_AB,small_rings_AB,limit_binds_AB,max_n_AB,first_n_AB,last_n_AB,"
           "T_BA,rows_BA,events_BA,first_event_BA,extra_maxima_BA,small_rings_BA,limit_binds_BA,max_n_BA,first_n_BA,last_n_BA,proposed,reason\n";
    for (int i = 0; i < np; ++i) {
        const LoadedPiece& P = pieces[i];
        const TriMesh& m = P.mesh;
        const std::string name = fs::path(P.entry.file).stem().string();
        const std::string dir = out + "/" + name;
        fs::create_directories(dir);
        const MeshCheck mc = checkMesh(m);
        const ComponentCheck& c0 = mc.components.front();
        // inside another piece of the same object (contact and overlap are reported, never treated as a join)
        std::vector<char> inside(m.pos.size(), 0);
        std::vector<int> insideOf(m.pos.size(), -1);   // the first piece of the same object containing the vertex
        for (int j = 0; j < np; ++j) {
            if (j == i || man.pieces[j].object != P.entry.object) continue;
            const MeshQuery q(common[j]);
            for (size_t v = 0; v < m.pos.size(); ++v)
                if (q.inside(common[i].pos[v])) {
                    inside[v] = 1;
                    if (insideOf[v] < 0) insideOf[v] = j;
                }
        }
        std::vector<double> vArea(m.pos.size(), 0);
        double area = 0, hiddenArea = 0;
        for (const auto& t : m.tri) {
            const double ar = 0.5 * norm(cross(m.pos[t[1]] - m.pos[t[0]], m.pos[t[2]] - m.pos[t[0]]));
            area += ar;
            for (int k = 0; k < 3; ++k) vArea[t[k]] += ar / 3;
        }
        for (size_t v = 0; v < m.pos.size(); ++v)
            if (inside[v]) hiddenArea += vArea[v];
        const Vec3 pA = (P.entry.poleA - P.map.centre) * g.sigma, pB = (P.entry.poleB - P.map.centre) * g.sigma;
        double dA = 0, dB = 0;
        const int vA = nearestVertex(m, pA, &dA), vB = nearestVertex(m, pB, &dB);
        const std::vector<std::vector<int>> nb = vertexNeighbours(m);
        std::vector<double> tAB, tBA;
        const DirectionAnalysis ab = analyseDirection(m, nb, inside, vA, vB, "A->B", w, h, dir, tAB);
        const DirectionAnalysis ba = analyseDirection(m, nb, inside, vB, vA, "B->A", w, h, dir, tBA);
        // the proposed direction: fewer contour events, then fewer extra maxima, then the end pole hidden inside another
        // piece (a gathered or threaded closure that does not show), then A -> B
        std::string proposed, reason;
        const int sAB = ab.events * 1000 + ab.extraMaxima, sBA = ba.events * 1000 + ba.extraMaxima;
        if (sAB != sBA) {
            proposed = sAB < sBA ? "A->B" : "B->A";
            reason = "fewer contour events / extra maxima";
        } else if (inside[vA] != inside[vB]) {
            proposed = inside[vB] ? "A->B" : "B->A";
            reason = "equal contour behaviour; the end pole lies inside another piece";
        } else {
            proposed = "A->B";
            reason = "equal contour behaviour and pole visibility; manifest order";
        }
        {
            std::ofstream mo(dir + "/reference.obj", std::ios::binary);
            for (const Vec3& p : m.pos) putf(mo, "v %.9g %.9g %.9g\n", p.x, p.y, p.z);
            for (const auto& t : m.tri) putf(mo, "f %d %d %d\n", t[0] + 1, t[1] + 1, t[2] + 1);
            std::ofstream vf(dir + "/vertices.csv", std::ios::binary);
            vf << "v,x,y,z,nx,ny,nz,t_AB,t_BA,inside_other_piece,inside_piece,area\n";
            for (size_t v = 0; v < m.pos.size(); ++v)
                putf(vf, "%zu,%.6f,%.6f,%.6f,%.5f,%.5f,%.5f,%.6f,%.6f,%d,%s,%.6f\n", v, m.pos[v].x, m.pos[v].y, m.pos[v].z, m.vnormal[v].x, m.vnormal[v].y, m.vnormal[v].z,
                     tAB[v], tBA[v], (int) inside[v], insideOf[v] < 0 ? "" : fs::path(man.pieces[insideOf[v]].file).stem().string().c_str(), vArea[v]);
            std::ofstream jf(dir + "/component.json", std::ios::binary);
            putf(jf, "{\n  \"piece\": \"%s\", \"object\": \"%s\", \"material\": \"%s\", \"source\": \"model (pieces manifest)\",\n", name.c_str(), P.entry.object.c_str(),
                 P.entry.material.c_str());
            putf(jf, "  \"sigma_mm_per_scene_unit\": %.6f, \"area_mm2\": %.3f, \"vertices\": %zu, \"triangles\": %zu,\n", g.sigma, area, m.pos.size(), m.tri.size());
            putf(jf, "  \"check\": {\"components\": %zu, \"closed\": %s, \"euler\": %d, \"consistent\": %s, \"self_intersection\": \"%s\"},\n", mc.components.size(),
                 c0.closed ? "true" : "false", c0.euler, c0.consistent ? "true" : "false", toString(c0.selfIntersection));
            putf(jf, "  \"poles\": {\"A\": {\"vertex\": %d, \"offset_mm\": %.4f, \"hidden\": %s}, \"B\": {\"vertex\": %d, \"offset_mm\": %.4f, \"hidden\": %s}},\n", vA, dA,
                 inside[vA] ? "true" : "false", vB, dB, inside[vB] ? "true" : "false");
            putf(jf, "  \"hidden_area_fraction\": %.5f,\n", area > 0 ? hiddenArea / area : 0);
            for (const DirectionAnalysis* d : {&ab, &ba})
                putf(jf, "  \"%s\": {\"T\": %.4f, \"t_max\": %.4f, \"rows\": %d, \"h_prime\": %.5f, \"levels\": %d, \"single_closed\": %d, \"events\": %d, \"first_event\": %d, "
                         "\"extra_maxima\": %d, \"small_rings\": %d, \"limit_binds\": %d, \"max_n\": %d, \"first_n\": %d, \"last_n\": %d, \"hidden_t\": [%.3f, %.3f], "
                         "\"cut_path_vertices\": %d, \"cut_path_reaches_start\": %s, \"cut_path_length\": %.3f},\n",
                     d->name == "A->B" ? "A_to_B" : "B_to_A", d->T, d->tMax, d->rows, d->hPrime, d->levels, d->singleClosed, d->events, d->firstEvent, d->extraMaxima,
                     d->smallRings, d->limitBinds, d->maxN, d->firstN, d->lastN, d->hiddenFirst, d->hiddenLast, d->pathVertices, d->pathReachesStart ? "true" : "false",
                     d->pathLength);
            putf(jf, "  \"proposed_direction\": \"%s\", \"reason\": \"%s\"\n}\n", proposed.c_str(), reason.c_str());
        }
        putf(sum, "%s,%s,%.2f,%zu,%zu,%zu,%d,%d,%d,%s,%d,%d,%.4f,", name.c_str(), P.entry.object.c_str(), area, m.pos.size(), m.tri.size(), mc.components.size(),
             c0.closed ? 1 : 0, c0.euler, c0.consistent ? 1 : 0, toString(c0.selfIntersection), (int) inside[vA], (int) inside[vB], area > 0 ? hiddenArea / area : 0);
        for (const DirectionAnalysis* d : {&ab, &ba})
            putf(sum, "%.3f,%d,%d,%d,%d,%d,%d,%d,%d,%d,", d->T, d->rows, d->events, d->firstEvent, d->extraMaxima, d->smallRings, d->limitBinds, d->maxN, d->firstN, d->lastN);
        putf(sum, "%s,%s\n", proposed.c_str(), reason.c_str());
        std::printf("%s: closed %d euler %d, hidden %.1f %%, A->B events %d maxima %d rows %d n<=%d; B->A events %d maxima %d; proposed %s (%s)\n", name.c_str(),
                    c0.closed ? 1 : 0, c0.euler, 100 * (area > 0 ? hiddenArea / area : 0), ab.events, ab.extraMaxima, ab.rows, ab.maxN, ba.events, ba.extraMaxima,
                    proposed.c_str(), reason.c_str());
    }
    std::ofstream gj(out + "/gauge.json", std::ios::binary);
    putf(gj, "{\"gauge_piece\": \"%s\", \"stitches\": %.0f, \"w\": %.4g, \"h\": %.4g, \"sigma\": %.6f, \"ring_length_scene\": %.6f, \"note\": \"the scene's calibration "
             "(shared scale, not a layout)\"}\n",
         gaugeFile.c_str(), stitches, w, h, g.sigma, g.ringLength);
    return 0;
}

// Reload: re-read candidate.obj, connectivity.csv, face_types.csv, ring_vertices.csv and shift_path.csv, re-run the
// check and the virtual shift, and compare the digest with report.json.
int commandLoad(const Args& a) {
    const std::string dir = a.get("--case-dir"), refDir = a.get("--dir");
    const double w = a.number("--w", 5.5), h = a.number("--h", 5.5);
    if (dir.empty() || refDir.empty()) {
        std::fprintf(stderr, "load: --case-dir and --dir are required\n");
        return 2;
    }
    TypedStitchMesh s;
    bool strips = false;
    std::string saved, cs;
    if (!readTypedFiles(dir, s, strips, saved, cs)) {
        std::fprintf(stderr, "load: cannot read %s\n", dir.c_str());
        return 1;
    }
    const std::string again = hex64(typedDigest(s));
    TriMesh m;
    std::vector<double> arc;
    std::string digestSaved, err;
    double T = 0;
    if (!readReference(refDir, m, arc, digestSaved, T, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    std::vector<int> castOn, bindOff;
    const std::vector<DeclaredBoundary> declared = declaredOf(m, cs, castOn, bindOff);
    TypedCheckOptions o;
    o.w = w;
    o.h = h;
    o.boundaryTolerance = 0.25 * h;
    o.bindingTolerance = 1e9;   // bindings are not reloaded; the vertices are compared through the digest
    o.faceDeviation = 0.1 * w;
    o.requireStrips = strips;
    s.binding.assign(s.mesh.pos.size(), VertexBinding{0, {1, 0, 0}, 0});
    const TypedMeshCheck c = checkTypedMesh(s, m, declared, o);
    VirtualShift vs;
    VirtualShiftCheck vc;
    const bool shift = strips && !s.shiftPath.empty() && buildVirtualShift(s, ShiftReading::PathSideDown, vs, &err);
    if (shift) vc = validateVirtualShift(s, vs);
    std::ofstream rp(dir + "/load_report.txt", std::ios::binary);
    putf(rp, "load: typed digest saved %s, reloaded %s: %s\n", saved.c_str(), again.c_str(), saved == again ? "equal" : "DIFFERENT");
    rp << c.text();
    if (shift) rp << vc.text();
    std::printf("load %s: digest %s (%s), connectivity %s, shift %s\n", dir.c_str(), again.c_str(), saved == again ? "equal" : "DIFFERENT", toString(c.connectivity),
                shift ? toString(vc.status) : "-");
    return saved == again ? 0 : 1;
}

// --- Surface markings: the Mark step's input and Method A on a marked piece, without a window.

// The file as the Mark step reads it (jobs/marking_session.h: the same reader, pieces and digest).
bool readMarkedFile(const std::string& path, MarkedFile& f, std::string& err) { return knit::readMarkedFile(path, f, &err); }

std::string nodeName(const MarkedFile& f, int node) {
    return node >= 0 && node < static_cast<int>(f.mesh->nodes.size()) ? f.mesh->nodes[node].name : "(no name)";
}

std::string jsonText(const std::string& s) {
    std::string o = "\"";
    for (char c : s) {
        if (c == '"' || c == '\\') o += '\\';
        o += c;
    }
    return o + "\"";
}

const char* levelWord(FindingLevel l) {
    switch (l) {
    case FindingLevel::Ok: return "ok";
    case FindingLevel::Warning: return "warning";
    case FindingLevel::Problem: return "problem";
    case FindingLevel::Incomplete: return "incomplete";
    }
    return "problem";
}
const char* causeWord(FindingCause c) {
    switch (c) {
    case FindingCause::Marking: return "marking";
    case FindingCause::AlgorithmLimitation: return "algorithm limitation";
    case FindingCause::Input: return "input";
    case FindingCause::SolverLimitation: return "solver limitation";
    }
    return "input";
}

void writePolylinesObj(const std::string& path, const std::vector<std::vector<Vec3>>& lines) {
    std::ofstream f(path, std::ios::binary);
    f << "# stitchmesh_lab marked: polylines in file units\n";
    int base = 1;
    for (const std::vector<Vec3>& l : lines) {
        for (const Vec3& p : l) putf(f, "v %.9g %.9g %.9g\n", p.x, p.y, p.z);
        if (l.size() >= 2) {
            f << "l";
            for (size_t i = 0; i < l.size(); ++i) f << " " << base + (int) i;
            f << "\n";
        }
        base += (int) l.size();
    }
}

//   stitchmesh_lab pieces --mesh <obj>  [--out <dir>]
// The pieces of an imported file (the Mark step's list): index, component, ordinal, triangles, area, bounds.
int commandPieces(const Args& a) {
    MarkedFile f;
    std::string err;
    if (!readMarkedFile(a.get("--mesh"), f, err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    std::ostringstream o;
    o << "piece,component,ordinal,triangles,area,lo_x,lo_y,lo_z,hi_x,hi_y,hi_z\n";
    for (size_t i = 0; i < f.pieces->pieces.size(); ++i) {
        const SurfacePiece& p = f.pieces->pieces[i];
        Vec3 lo, hi;
        p.mesh.bounds(lo, hi);
        putf(o, "%zu,%s,%d,%zu,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", i, nodeName(f, p.node).c_str(), p.ordinal, p.mesh.tri.size(), p.mesh.area(), lo.x, lo.y,
             lo.z, hi.x, hi.y, hi.z);
    }
    std::printf("%s: %zu triangles, digest %s, %zu pieces\n%s", a.get("--mesh").c_str(), f.mesh->tri.size(), hex64(f.digest).c_str(), f.pieces->pieces.size(),
                o.str().c_str());
    if (!a.get("--out").empty()) {
        fs::create_directories(a.get("--out"));
        std::ofstream(a.get("--out") + "/pieces.csv", std::ios::binary) << o.str();
    }
    return 0;
}

//   stitchmesh_lab snap --mesh <obj> --points <txt> --out <file.markings> [--piece <index>]
// Markings from 3D points (file units): each point goes to the closest surface point (of one piece when --piece is
// given). The points file holds blocks
//   marking <meaning> <shape> [proposal]   meaning and shape as in the markings file; "proposal": offered to the user
//   note <text>  /  reason <text>  /  unsuitable <text>
//   p <x> <y> <z>
//   end
// so a script (or the knitting agent's proposal) can place markings without clicking.
int commandSnap(const Args& a) {
    MarkedFile f;
    std::string err;
    if (!readMarkedFile(a.get("--mesh"), f, err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const int only = (int) a.number("--piece", -1);
    std::ifstream in(a.get("--points"));
    if (!in) {
        std::fprintf(stderr, "cannot read %s\n", a.get("--points").c_str());
        return 1;
    }
    MarkingSet set;
    set.source = a.get("--mesh");
    set.digest = f.digest;
    SurfaceMarking m;
    bool open = false;
    std::string line;
    int n = 0;
    while (std::getline(in, line)) {
        ++n;
        std::istringstream ls(line);
        std::string key;
        if (!(ls >> key) || key[0] == '#') continue;
        std::string rest;
        std::getline(ls, rest);
        const size_t first = rest.find_first_not_of(' ');
        rest = first == std::string::npos ? "" : rest.substr(first);
        if (key == "marking") {
            std::istringstream ms(rest);
            std::string meaning, shape, proposal;
            ms >> meaning >> shape >> proposal;
            m = SurfaceMarking{};
            m.id = set.nextId++;
            if (!parseMeaning(meaning, m.meaning) || !parseShape(shape, m.shape)) {
                std::fprintf(stderr, "line %d: unknown meaning or shape\n", n);
                return 1;
            }
            if (proposal == "proposal") {
                m.origin = MarkingOrigin::Proposal;
                m.proposal = ProposalState::Offered;
                m.proposalId = m.id;
            }
            open = true;
        } else if (key == "note") {
            m.note = rest;
        } else if (key == "reason") {
            m.reason = rest;
        } else if (key == "unsuitable") {
            m.unsuitable = rest;
        } else if (key == "p" && open) {
            std::istringstream ps(rest);
            Vec3 p;
            if (!(ps >> p.x >> p.y >> p.z)) {
                std::fprintf(stderr, "line %d: a point is three numbers\n", n);
                return 1;
            }
            // the closest point of the file (or of one piece), as a file triangle
            double best = 1e300;
            SurfacePoint sp;
            for (int pi = 0; pi < (int) f.pieces->pieces.size(); ++pi) {
                if (only >= 0 && pi != only) continue;
                const SurfacePiece& piece = f.pieces->pieces[pi];
                const SurfacePoint local = closestSurfacePoint(piece.mesh, p);
                const double d = norm(surfacePosition(piece.mesh, local) - p);
                if (d < best) {
                    best = d;
                    sp = local;
                    sp.tri = piece.fileTriangle[local.tri];
                }
            }
            m.points.push_back(sp);
        } else if (key == "end" && open) {
            set.markings.push_back(m);
            open = false;
        } else {
            std::fprintf(stderr, "line %d: unexpected \"%s\"\n", n, key.c_str());
            return 1;
        }
    }
    if (!saveMarkingFile(a.get("--out"), set, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    std::printf("snap: %zu markings written to %s\n", set.markings.size(), a.get("--out").c_str());
    return 0;
}

//   stitchmesh_lab marked --mesh <obj> --markings <file> --unit-mm <mm per file unit> --w <mm> --h <mm>
//                         [--piece <index>] [--depth built|relaxed] --out <dir>
// Method A on every marked piece (or one): per piece review.json (findings, the four measure groups), stitch_mesh.obj,
// rings.obj and traces.obj (polylines), regions.csv (per piece triangle: its file triangle and region) and log.txt;
// all geometry in file units, as the Mark step shows it.
std::vector<std::string> splitList(const std::string& text) {
    std::vector<std::string> out;
    std::stringstream ss(text);
    std::string x;
    while (std::getline(ss, x, ','))
        if (!x.empty()) out.push_back(x);
    return out;
}

// The integrated build's evaluation (marked_build.h MarkedEvaluation), in mm: before/ and after/ typed files, the solver's
// faces, cells, coverage, the relaxation history and summary.csv (an analysis script
// outside this repository reads its keys).
constexpr double kPi = 3.14159265358979323846;

struct EvalStats {
    int measured = 0, within = 0, k = 0, kWithin = 0, pent = 0, pentWithin = 0, pole = 0, poleWithin = 0;
    double median = 0, p95 = 0, wLo = 0, wHi = 0, hLo = 0, hHi = 0, aLo = 0, aHi = 0;
};

EvalStats evalStats(const TypedStitchMesh& s, double w, double h) {
    EvalStats e;
    std::vector<double> dev, wr, hr, an;
    for (const CellMeasure& c : measureCells(s, w, h)) {
        if (c.face < (int) s.reason.size() && s.reason[c.face] == "helix-start-step") continue;
        if (c.kind == StitchFaceKind::Unsupported) continue;
        ++e.measured;
        e.within += c.withinBounds;
        if (c.kind == StitchFaceKind::Regular) {
            ++e.k;
            e.kWithin += c.withinBounds;
            dev.push_back(std::fabs(c.widthRatio - 1));
            dev.push_back(std::fabs(c.heightRatio - 1));
            wr.push_back(c.widthRatio);
            hr.push_back(c.heightRatio);
            an.push_back(c.minAngle);
            an.push_back(c.maxAngle);
        } else if (c.kind == StitchFaceKind::PoleStart) {
            ++e.pole;
            e.poleWithin += c.withinBounds;
        } else {
            ++e.pent;
            e.pentWithin += c.withinBounds;
        }
    }
    e.median = quantile(dev, 0.5);
    e.p95 = quantile(dev, 0.95);
    e.wLo = quantile(wr, 0);
    e.wHi = quantile(wr, 1);
    e.hLo = quantile(hr, 0);
    e.hHi = quantile(hr, 1);
    e.aLo = quantile(an, 0);
    e.aHi = quantile(an, 1);
    return e;
}

void writeShapingCells(const std::string& path, const std::vector<ShapingCellMeasure>& cells) {
    std::ofstream f(path, std::ios::binary);
    f << "face,kind,min_course_w,max_course_w,min_wale_h,max_wale_h,sub_quads,min_corner,max_corner,min_sub_area_wh,folded_sub_quads,interfaces_within,sub_quads_within,within,outside\n";
    for (const ShapingCellMeasure& c : cells)
        putf(f, "%d,%s,%.5f,%.5f,%.5f,%.5f,%d,%.3f,%.3f,%.5f,%d,%d,%d,%d,\"%s\"\n", c.face, toString(c.kind), c.minCourse, c.maxCourse, c.minWale, c.maxWale, c.subQuads,
             c.minCorner, c.maxCorner, c.minSubArea, c.foldedSubQuads, c.interfacesWithin ? 1 : 0, c.subQuadsWithin ? 1 : 0, c.within ? 1 : 0, c.outside.c_str());
}

void writeEvaluation(const std::string& dir, const MarkedBuildProduct& p, const MarkedBuildRequest& rq) {
    const MarkedEvaluation& ev = p.evaluation;
    if (ev.before.mesh.faces() == 0) return;
    fs::create_directories(dir + "/before");
    fs::create_directories(dir + "/after");
    const TypedStitchMesh vb = virtualTypedMesh(ev.before, ev.shift), va = virtualTypedMesh(ev.after, ev.shift);
    writeTypedMesh(dir + "/before", ev.before, ev.storedBefore, &ev.shift, &ev.shiftBefore);
    writeTypedMesh(dir + "/after", ev.after, ev.storedAfter, &ev.shift, &ev.shiftAfter);
    writeObjPolygons(dir + "/before/virtual.obj", vb);
    writeObjPolygons(dir + "/after/virtual.obj", va);
    writeCellsTo(dir + "/before/cells_virtual.csv", measureCells(vb, rq.w, rq.h));
    writeCellsTo(dir + "/after/cells_virtual.csv", measureCells(va, rq.w, rq.h));
    writeCellsTo(dir + "/after/cells_stored.csv", measureCells(ev.after, rq.w, rq.h));
    {
        std::ofstream f(dir + "/knit_region.obj", std::ios::binary);
        for (const Vec3& q : ev.knitRegion.pos) putf(f, "v %.17g %.17g %.17g\n", q.x, q.y, q.z);
        for (const auto& t : ev.knitRegion.tri) putf(f, "f %d %d %d\n", t[0] + 1, t[1] + 1, t[2] + 1);
    }
    {
        std::ofstream f(dir + "/reference.obj", std::ios::binary);
        for (const Vec3& q : ev.reference.pos) putf(f, "v %.17g %.17g %.17g\n", q.x, q.y, q.z);
        for (const auto& t : ev.reference.tri) putf(f, "f %d %d %d\n", t[0] + 1, t[1] + 1, t[2] + 1);
    }
    {
        // the markings' constraints as the solver used them (row-start path, finish contour), one polyline each
        std::vector<std::vector<Vec3>> lines;
        for (const GuideLine& g : ev.guides.lines) {
            lines.push_back(g.points);
            if (g.closed && !g.points.empty()) lines.back().push_back(g.points.front());
        }
        writePolylinesObj(dir + "/guides.obj", lines);
    }
    if (ev.relaxed) {
        std::ofstream hf(dir + "/after/relax_history.csv", std::ios::binary);
        hf << "step,rms_force,max_force,max_force_vertex,max_move,cg_iterations,max_force_tri,max_force_edge_distance\n";
        for (const RelaxStep& st : ev.relax.history)
            putf(hf, "%d,%.9g,%.9g,%d,%.9g,%d,%d,%.6g\n", st.step, st.rmsForce, st.maxForce, st.maxForceVertex, st.maxMove, st.cgIterations, st.maxForceTri,
                 st.maxForceEdgeDistance);
    }
    std::ofstream r(dir + "/evaluation.txt", std::ios::binary);
    r << "route " << ev.route << ", reference refined x" << ev.refineLevels << ", " << ev.reference.tri.size() << " triangles\n";
    r << "--- stored before\n" << ev.storedBefore.text() << "--- virtual before\n" << ev.virtualBefore.text() << ev.shiftBefore.text();
    r << "--- coverage before, stored\n" << ev.coverStoredBefore.text() << "--- coverage before, virtual\n" << ev.coverVirtualBefore.text();
    if (ev.relaxed) r << "--- relaxation\n" << relaxText(ev.relax);
    r << "--- stored after\n" << ev.storedAfter.text() << "--- virtual after\n" << ev.virtualAfter.text() << ev.shiftAfter.text();
    r << "--- coverage after, stored\n" << ev.coverStoredAfter.text() << "--- coverage after, virtual\n" << ev.coverVirtualAfter.text();
    std::ofstream f(dir + "/summary.csv", std::ios::binary);
    const EvalStats b0 = evalStats(vb, rq.w, rq.h), b1 = evalStats(va, rq.w, rq.h);
    const LayoutMeasures& m = p.review.measures;
    f << "key,before,after\n";
    putf(f, "route,%s,\nrefine,%d,\nfaces,%d,%d\nunsupported,%d,%d\n", ev.route.c_str(), ev.refineLevels, ev.storedBefore.faces, ev.storedAfter.faces, ev.storedBefore.unsupported,
         ev.storedAfter.unsupported);
    putf(f, "connectivity_stored,%s,%s\nconnectivity_virtual,%s,%s\ngeometry_virtual,%s,%s\ngeometry_stored,%s,%s\nvirtual_shift,%s,%s\n", toString(ev.storedBefore.connectivity),
         toString(ev.storedAfter.connectivity), toString(ev.virtualBefore.connectivity), toString(ev.virtualAfter.connectivity), toString(ev.virtualBefore.geometry),
         toString(ev.virtualAfter.geometry), toString(ev.storedBefore.geometry), toString(ev.storedAfter.geometry), toString(ev.shiftBefore.status), toString(ev.shiftAfter.status));
    putf(f, "coverage_stored,%s,%s\ncoverage_virtual,%s,%s\n", toString(ev.coverStoredBefore.status), toString(ev.coverStoredAfter.status), toString(ev.coverVirtualBefore.status),
         toString(ev.coverVirtualAfter.status));
    putf(f, "uncovered_fraction_stored,%.5f,%.5f\nuncovered_fraction_virtual,%.5f,%.5f\ncovered_twice_virtual,%d,%d\nprojected_ratio_virtual,%.5f,%.5f\n",
         ev.coverStoredBefore.uncoveredFraction, ev.coverStoredAfter.uncoveredFraction, ev.coverVirtualBefore.uncoveredFraction, ev.coverVirtualAfter.uncoveredFraction,
         ev.coverVirtualBefore.coveredTwice, ev.coverVirtualAfter.coveredTwice, ev.coverVirtualBefore.projectedRatio, ev.coverVirtualAfter.projectedRatio);
    putf(f, "hausdorff_stored,%.4f,%.4f\nhausdorff_virtual,%.4f,%.4f\n", ev.coverStoredBefore.worstHausdorff, ev.coverStoredAfter.worstHausdorff, ev.coverVirtualBefore.worstHausdorff,
         ev.coverVirtualAfter.worstHausdorff);
    putf(f, "cells_measured,%d,%d\ncells_within,%d,%d\nk_cells,%d,%d\nk_within,%d,%d\npentagons,%d,%d\npentagons_within,%d,%d\npole_cells,%d,%d\npole_within,%d,%d\n", b0.measured,
         b1.measured, b0.within, b1.within, b0.k, b1.k, b0.kWithin, b1.kWithin, b0.pent, b1.pent, b0.pentWithin, b1.pentWithin, b0.pole, b1.pole, b0.poleWithin, b1.poleWithin);
    putf(f, "median_k,%.5f,%.5f\np95_k,%.5f,%.5f\nwidth_lo,%.4f,%.4f\nwidth_hi,%.4f,%.4f\nheight_lo,%.4f,%.4f\nheight_hi,%.4f,%.4f\nangle_lo,%.2f,%.2f\nangle_hi,%.2f,%.2f\n", b0.median,
         b1.median, b0.p95, b1.p95, b0.wLo, b1.wLo, b0.wHi, b1.wHi, b0.hLo, b1.hLo, b0.hHi, b1.hHi, b0.aLo, b1.aLo, b0.aHi, b1.aHi);
    putf(f, "max_face_deviation_stored,%.4f,%.4f\nmax_face_deviation_virtual,%.4f,%.4f\ninverted_virtual,%d,%d\nintersecting_pairs_virtual,%d,%d\nfaces_over_deviation_stored,%d,%d\n",
         ev.storedBefore.maxFaceDeviation, ev.storedAfter.maxFaceDeviation, ev.virtualBefore.maxFaceDeviation, ev.virtualAfter.maxFaceDeviation, ev.virtualBefore.invertedFaces,
         ev.virtualAfter.invertedFaces, ev.virtualBefore.intersectingPairs, ev.virtualAfter.intersectingPairs, ev.storedBefore.facesOverDeviation, ev.storedAfter.facesOverDeviation);
    putf(f, "relax_status,,%s\nrelax_stop,,%s\nrelax_steps,,%d\nrelax_ms,,%.1f\nrms_force,%.3e,%.3e\npole_drift_mm,,%.4f\nrow_origin_mean_mm,,%.4f\nrow_origin_max_mm,,%.4f\n",
         ev.relaxed ? toString(ev.relax.status) : "not run", ev.relax.stopReason.c_str(), ev.relax.steps, ev.relax.ms, ev.relax.rmsForceStart, ev.relax.rmsForceEnd, ev.poleMoved,
         ev.rowOriginMean, ev.rowOriginMax);
    putf(f, "rows,%d,\nfinish_level,%.4f,\nstitch_vertices,%zu,\nknit_area_mm2,%.2f,\nomitted_mm2,%.2f,\nomitted_hidden_mm2,%.2f,\nomitted_visible_mm2,%.2f,\nseconds,%.3f,\n", m.rows,
         ev.finishLevel, ev.before.mesh.pos.size(), ev.knitRegion.area(), m.omittedAreaMm2, m.omittedHiddenMm2, m.omittedVisibleMm2, m.seconds);
    if (!ev.before.rings.empty() && ev.before.rings[0].size() == 1) putf(f, "ring1,%zu,\n", ev.before.rings[1].size());
    putf(f, "preparation,%s,\nrelax_preparation,%s,\nreference_digest,%016llx,\nchange_max_mm,%.4f,\nchange_p95_mm,%.4f,\nchange_back_max_mm,%.4f,\nchange_back_p95_mm,%.4f,\n"
         "area_original_mm2,%.2f,\narea_reference_mm2,%.2f,\n", ev.preparation.c_str(), ev.relaxPreparation.c_str(), (unsigned long long) ev.referenceDigest, ev.changeMax,
         ev.changeP95, ev.changeBackMax, ev.changeBackP95, ev.areaOriginal, ev.areaReference);
    if (rq.options.relaxOnOther)
        putf(f, "other_guides,%s,\n", rq.options.otherGuides == MarkedBuildOptions::OtherGuides::Carried ? "carried"
                                       : rq.options.otherGuides == MarkedBuildOptions::OtherGuides::FinishPoint ? "finish_point" : "matched_level");
    if (ev.relaxPreparation != "facets")
        putf(f, "max_face_deviation_vs_original,,%.4f\nfaces_over_deviation_vs_original,,%d\n", ev.storedAfterVsOriginal.maxFaceDeviation, ev.storedAfterVsOriginal.facesOverDeviation);
    {
        double control = 0, trace = 0;
        for (const MarkedEvaluation::MarkingShift& ms : ev.markingShift) {
            control = std::max(control, ms.controlMax);
            trace = std::max(trace, ms.traceHausdorff);
        }
        putf(f, "marking_control_moved_mm,%.4f,\nmarking_trace_hausdorff_mm,%.4f,\n", control, trace);
    }
    {
        // the shaping-cell assessment beside the frozen pentagon bounds
        const ShapingSummary sb = summarizeShaping(measureShapingCells(vb, rq.w, rq.h)), sa = summarizeShaping(measureShapingCells(va, rq.w, rq.h));
        putf(f, "shaping_cells,%d,%d\nshaping_subquad_within,%d,%d\nshaping_interfaces_within,%d,%d\nshaping_assessed_within,%d,%d\nshaping_inverted_subquads,%d,%d\n", sb.cells,
             sa.cells, sb.subQuadWithin, sa.subQuadWithin, sb.interfacesWithin, sa.interfacesWithin, sb.within, sa.within, sb.invertedSubQuads, sa.invertedSubQuads);
        writeShapingCells(dir + "/after/cells_shaping.csv", measureShapingCells(va, rq.w, rq.h));
    }
    if (ev.relaxed) {
        // stationarity of the final state under the actual faceted constraint (whatever filter the solver used)
        MeshRelaxRequest r;
        r.w = rq.w;
        r.h = rq.h;
        const bool onRegion = !ev.relaxLoops.empty();
        const WalkResidual wr = walkStationarity(onRegion ? ev.knitRegion : ev.reference, ev.solverFaces, ev.after.mesh.pos, ev.border, ev.relaxLoops, r, &ev.guides, 1e-3);
        putf(f, "walk_rms,,%.4e\nwalk_max,,%.4e\nwalk_max_vertex,,%d\nwalk_on_edge,,%d\nwalk_on_vertex,,%d\n", wr.rms, wr.max, wr.maxVertex, wr.onEdge, wr.onVertex);
        std::ofstream wf(dir + "/after/walk_residual.csv", std::ios::binary);
        wf << "vertex,walk_residual\n";
        for (size_t v = 0; v < wr.perVertex.size(); ++v) putf(wf, "%zu,%.6e\n", v, wr.perVertex[v]);
    }
}

// The creased-surface fixture of solver note R1: a ring tube relaxed on a 16-sided prism (creases of 22.5 degrees), once with
// the facet filter (CP3) and once with the smooth filter normals (R1); each result's stationarity under the faceted constraint,
// and each result restarted under the other filter (a feasible perturbation test).
int commandRelaxCrease(const Args& a) {
    const std::string out = a.get("--out");
    fs::create_directories(out);
    const int sides = (int) a.number("--sides", 16), around = (int) a.number("--stitches", 17), rows = (int) a.number("--rows", 8);
    const double radius = a.number("--radius", 15), w = 5.5, h = 3.4, Z = rows * h;
    // the prism, rings 0.85 mm apart, vertices of each ring on the polygon's corners and edge midpoints
    TriMesh prism;
    const int per = sides * 4, nz = (int) std::lround(Z / 0.85);
    for (int k = 0; k <= nz; ++k)
        for (int j = 0; j < per; ++j) {
            const int side = j / 4;
            const double f = (j % 4) / 4.0;
            const double a0 = 2 * kPi * side / sides, a1 = 2 * kPi * (side + 1) / sides;
            const Vec3 c0{radius * std::cos(a0), radius * std::sin(a0), 0}, c1{radius * std::cos(a1), radius * std::sin(a1), 0};
            Vec3 q = c0 * (1 - f) + c1 * f;
            q.z = Z * k / nz;
            prism.pos.push_back(q);
        }
    auto id = [&](int k, int j) { return k * per + ((j % per) + per) % per; };
    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < per; ++j) {
            prism.tri.push_back({id(k, j), id(k, j + 1), id(k + 1, j + 1)});
            prism.tri.push_back({id(k, j), id(k + 1, j + 1), id(k + 1, j)});
        }
    prism.computeNormals();
    std::vector<std::vector<int>> loops = triangleBoundaryLoops(prism);
    TypedStitchMesh tube;
    TriMesh unusedRef;
    std::vector<DeclaredBoundary> decl;
    std::string err;
    ringTubeFixture(std::vector<int>(rows + 1, around), std::vector<std::vector<int>>(rows), radius, h, tube, unusedRef, decl, &err);
    {
        const MeshQuery q(prism);
        for (Vec3& p : tube.mesh.pos) p = q.closest(p).p;
    }
    VirtualShift vs;
    buildVirtualShift(tube, ShiftReading::PathSideDown, vs, &err);
    std::vector<RelaxFace> faces;
    for (int f = 0; f < tube.mesh.faces(); ++f) faces.push_back(RelaxFace{vs.corners[f], tube.role[f]});
    std::vector<char> border(tube.mesh.pos.size(), 0);
    for (const std::vector<int>& l : polygonBorderLoops(tube.mesh))
        for (int v : l) border[v] = 1;
    std::ofstream rep(out + "/crease_report.txt", std::ios::binary);
    putf(rep, "prism: %d sides (crease %.2f deg), radius %.1f mm, height %.1f mm, %zu triangles; tube %d around x %d rows, gauge %.1f x %.1f\n", sides, 360.0 / sides, radius, Z,
         prism.tri.size(), around, rows, w, h);
    std::ofstream csv(out + "/crease.csv", std::ios::binary);
    csv << "run,filter,start,status,stop,steps,rms_force_end,walk_rms,walk_max,walk_on_edge,max_move_from_facet_mm\n";
    std::vector<Vec3> facetPos, smoothPos;
    bool kinks = false;
    auto run = [&](const char* name, bool smooth, const std::vector<Vec3>& start, const char* startName, std::vector<Vec3>* keep) {
        MeshRelaxRequest r;
        r.w = w;
        r.h = h;
        r.smoothNormals = smooth;
        r.kinkConstraints = kinks;
        RelaxProduct p;
        relaxStitchMesh(prism, loops, faces, start, border, r, p, &err);
        MeshRelaxRequest plain = r;
        const WalkResidual wr = walkStationarity(prism, faces, p.pos, border, loops, plain, nullptr, 1e-3);
        double moved = 0;
        if (!facetPos.empty())
            for (size_t v = 0; v < p.pos.size(); ++v) moved = std::max(moved, norm(p.pos[v] - facetPos[v]));
        putf(csv, "%s,%s,%s,%s,%s,%d,%.4e,%.4e,%.4e,%d,%.4f\n", name, smooth ? "smooth" : "facet", startName, toString(p.status), p.stopReason.c_str(), p.steps, p.rmsForceEnd, wr.rms,
             wr.max, wr.onEdge, moved);
        putf(rep, "%-22s filter %-6s from %-14s: %s (%s) after %d steps, RMS force %.3e; walk residual rms %.3e max %.3e (vertex %d; %d on an edge); %.4f mm from the facet run\n",
             name, smooth ? "smooth" : "facet", startName, toString(p.status), p.stopReason.c_str(), p.steps, p.rmsForceEnd, wr.rms, wr.max, wr.maxVertex, wr.onEdge, moved);
        std::ofstream hist(out + "/" + name + "_history.csv", std::ios::binary);
        hist << "step,rms_force,max_force,max_force_vertex,max_move\n";
        for (const RelaxStep& st : p.history) putf(hist, "%d,%.9g,%.9g,%d,%.9g\n", st.step, st.rmsForce, st.maxForce, st.maxForceVertex, st.maxMove);
        if (keep) *keep = p.pos;
    };
    run("facet", false, tube.mesh.pos, "projected tube", &facetPos);
    run("smooth_R1", true, tube.mesh.pos, "projected tube", &smoothPos);
    run("facet_from_R1", false, smoothPos, "the R1 result", nullptr);
    run("smooth_from_facet", true, facetPos, "the facet result", nullptr);
    kinks = true;
    run("kinks_R2", false, tube.mesh.pos, "projected tube", nullptr);
    std::printf("relax-crease: written to %s\n", out.c_str());
    return 0;
}

int commandMarked(const Args& a) {
    MarkedFile f;
    std::string err;
    if (!readMarkedFile(a.get("--mesh"), f, err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    MarkingSet set;
    if (!loadMarkingFile(a.get("--markings"), set, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    if (set.digest != f.digest) {
        std::fprintf(stderr, "the markings were made on another mesh (digest %s, this file %s): stale, not applied\n", hex64(set.digest).c_str(),
                     hex64(f.digest).c_str());
        return 1;
    }
    const std::string out = a.get("--out");
    fs::create_directories(out);
    std::set<int> marked;
    for (const SurfaceMarking& m : set.markings) {
        const int p = pieceOfMarking(*f.pieces, m);
        if (p >= 0) marked.insert(p);
    }
    const int only = (int) a.number("--piece", -1);
    for (int pi : marked) {
        if (only >= 0 && pi != only) continue;
        MarkedBuildRequest rq;
        rq.pieces = f.pieces.get();
        rq.piece = pi;
        rq.markings = &set;
        rq.mmPerUnit = a.number("--unit-mm", 1);
        rq.w = a.number("--w", 5.5);
        rq.h = a.number("--h", 5.5);
        rq.depth = a.get("--depth", "relaxed") == "built" ? ReviewDepth::Built : ReviewDepth::Relaxed;
        rq.refine = static_cast<int>(a.number("--refine", -1));
        rq.options.adjacentShaping = a.get("--adjacent", "1") == "1";
        rq.options.boundaryShaping = a.get("--boundary-shaping", "1") == "1";
        // --anchor ends|path|none (RelaxAnchors; default: the app's)
        const std::string anchor = a.get("--anchor", "none");
        rq.options.anchors = anchor == "path" ? RelaxAnchors::RowOriginsOnPath : anchor == "ends" ? RelaxAnchors::YarnEnds : RelaxAnchors::None;
        // --border-rings placed|slide (BorderRings; default: the app's)
        const std::string rings = a.get("--border-rings", "first");
        rq.options.borderRings = rings == "slide" ? BorderRings::SlideInRelaxation : rings == "placed" ? BorderRings::PlacedOnLoop : BorderRings::PlacedFirst;
        const std::string order = a.get("--border-order", "kept");
        rq.options.borderOrder = order == "free" ? BorderOrder::Free : order == "ring" ? BorderOrder::RingSlidesAsOne : BorderOrder::Kept;
        rq.options.borderSprings = a.get("--border-springs", "stored") == "virtual" ? BorderSprings::VirtualFacesOnly : BorderSprings::StoredBorderEdges;
        rq.options.rowPhase = a.get("--row-phase", "wales") == "path" ? RowPhase::OnPath : RowPhase::WalesUpright;
        if (!a.get("--opening-rows").empty())
            rq.options.openingRows = a.get("--opening-rows") == "distance" ? OpeningRows::DistanceFromStart : OpeningRows::OpeningIsALevel;
        if (!a.get("--bind-off-count").empty())
            rq.options.bindOffCount = a.get("--bind-off-count") == "below" ? BindOffCount::OfRoundBelow : BindOffCount::FromLoopLength;
        rq.options.faceCorners = a.get("--face-corners", "kept") == "free" ? FaceCorners::Free : FaceCorners::Kept;
        if (!a.get("--step-control").empty())
            rq.options.stepControl = a.get("--step-control") == "fixed" ? StepControl::FixedTimeStep : StepControl::ForcesMayNotGrow;
        rq.options.smoothFilterNormals = a.get("--filter-normals", "facet") == "smooth";
        rq.options.kinkConstraints = a.get("--kinks", "free") == "held";
        // --preparation facets|loop_limit (default: the markings file's own); --relax-on loop_limit: the controlled
        // comparison (the layout as built, carried to that preparation of the same split and relaxed there)
        parsePreparation(set.preparation, rq.options.preparation);   // the reader admits only the two values
        if (!a.get("--preparation").empty() && !parsePreparation(a.get("--preparation"), rq.options.preparation)) {
            std::fprintf(stderr, "--preparation: facets or loop_limit\n");
            return 1;
        }
        if (!a.get("--relax-on").empty()) {
            rq.options.relaxOnOther = true;
            if (!parsePreparation(a.get("--relax-on"), rq.options.relaxPreparation)) {
                std::fprintf(stderr, "--relax-on: facets or loop_limit\n");
                return 1;
            }
            // --other-guides matched|carried|finish_point (MarkedBuildOptions::OtherGuides)
            const std::string og = a.get("--other-guides", "matched");
            rq.options.otherGuides = og == "carried" ? MarkedBuildOptions::OtherGuides::Carried
                                   : og == "finish_point" ? MarkedBuildOptions::OtherGuides::FinishPoint : MarkedBuildOptions::OtherGuides::MatchedLevel;
        }
        // --preview 3,5,7: offered proposals evaluated as if accepted (a research preview; the file is not changed)
        for (const std::string& id : splitList(a.get("--preview"))) rq.previewOffered.push_back(std::stoi(id));
        MarkedBuildProduct p;
        if (!buildFromMarkings(rq, p, &err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        const SurfacePiece& piece = f.pieces->pieces[pi];
        char name[64];
        std::snprintf(name, sizeof name, "piece_%03d", pi);
        const std::string dir = out + "/" + name;
        fs::create_directories(dir);
        const MarkingReview& r = p.review;
        const LayoutMeasures& m = r.measures;
        std::ofstream j(dir + "/review.json", std::ios::binary);
        putf(j, "{\n  \"piece\": %d, \"component\": %s, \"ordinal\": %d, \"triangles\": %zu,\n", pi, jsonText(nodeName(f, piece.node)).c_str(), piece.ordinal,
             piece.mesh.tri.size());
        putf(j, "  \"gauge\": {\"w\": %.6g, \"h\": %.6g, \"mm_per_unit\": %.9g}, \"depth\": \"%s\", \"valid\": %s,\n", rq.w, rq.h, rq.mmPerUnit,
             r.depth == ReviewDepth::Relaxed ? "relaxed" : "built", r.valid ? "true" : "false");
        putf(j, "  \"validity\": {\"contour_events\": %d, \"unsupported_cells\": %d},\n", m.contourEvents, m.unsupportedCells);
        putf(j, "  \"layout\": {\"rows\": %d, \"max_stitch_count_jump\": %d, \"exceptional_cells\": %d, \"distortion_p95\": %.6g},\n", m.rows, m.maxStitchCountJump,
             m.exceptionalCells, m.distortionP95);
        putf(j, "  \"coverage\": {\"covered_fraction\": %.6g, \"omitted_area_mm2\": %.6g, \"omitted_hidden_mm2\": %.6g, \"omitted_visible_mm2\": %.6g, \"pass\": %s},\n",
             m.coveredFraction, m.omittedAreaMm2, m.omittedHiddenMm2, m.omittedVisibleMm2, m.coveragePass ? "true" : "false");
        putf(j, "  \"route\": \"%s\", \"refine\": %d, \"preview_offered\": %zu,\n", p.evaluation.route.c_str(), p.evaluation.refineLevels, p.previewIds.size());
        putf(j, "  \"cost\": {\"seconds\": %.4f},\n  \"findings\": [\n", m.seconds);
        for (size_t i = 0; i < r.findings.size(); ++i) {
            const MarkingFinding& x = r.findings[i];
            putf(j, "    {\"marking\": %d, \"level\": \"%s\", \"cause\": \"%s\", \"message\": %s}%s\n", x.markingId, levelWord(x.level), causeWord(x.cause),
                 jsonText(x.message).c_str(), i + 1 < r.findings.size() ? "," : "");
        }
        putf(j, "  ],\n  \"recommendation\": %s\n}\n", jsonText(r.recommendation).c_str());
        writeObjPolygons(dir + "/stitch_mesh.obj", p.stitchMesh);
        writePolylinesObj(dir + "/rings.obj", p.rings);
        std::vector<std::vector<Vec3>> traces;
        for (const MarkingTrace& t : p.traces) {
            std::vector<Vec3> l;
            for (const TracePoint& x : t.points) l.push_back(x.pos);
            traces.push_back(l);
        }
        writePolylinesObj(dir + "/traces.obj", traces);
        std::ofstream rg(dir + "/regions.csv", std::ios::binary);
        rg << "piece_triangle,file_triangle,region\n";
        for (size_t t = 0; t < p.region.size(); ++t) putf(rg, "%zu,%d,%s\n", t, piece.fileTriangle[t], toString(p.region[t]));
        std::ofstream lg(dir + "/log.txt", std::ios::binary);
        for (const std::string& l : p.log) lg << l << "\n";
        writeEvaluation(dir, p, rq);
        std::printf("piece %d (%s %d): %s; rows %d, cells %d, exceptional %d, unsupported %d, covered %.4f, omitted %.0f mm2, p95 %.4f, %.2f s\n", pi,
                    nodeName(f, piece.node).c_str(), piece.ordinal, r.valid ? "built" : "NOT built", m.rows, p.stitchMesh.mesh.faces(), m.exceptionalCells,
                    m.unsupportedCells, m.coveredFraction, m.omittedAreaMm2, m.distortionP95, m.seconds);
        for (const MarkingFinding& x : r.findings) std::printf("  [%s, %s] marking %d: %s\n", levelWord(x.level), causeWord(x.cause), x.markingId, x.message.c_str());
    }
    return 0;
}

} // namespace

// The writers and readers the arm commands share (stitchmesh_lab_shared.h).
namespace knit::lab {
void writeTypedMeshFiles(const std::string& dir, const TypedStitchMesh& s, const TypedMeshCheck& c, const VirtualShift* vs, const VirtualShiftCheck* vc) {
    writeTypedMesh(dir, s, c, vs, vc);
}
bool readTypedMeshFiles(const std::string& dir, TypedStitchMesh& s, std::string& savedDigest) {
    bool strips = false;
    std::string cs;
    return readTypedFiles(dir, s, strips, savedDigest, cs);
}
void writeCellsFile(const std::string& path, const std::vector<CellMeasure>& cells) { writeCellsTo(path, cells); }
void writePolygonObj(const std::string& path, const TypedStitchMesh& s) { writeObjPolygons(path, s); }
void writeRelaxHistory(const std::string& path, const RelaxProduct& p) { writeHistory(path, p); }
std::string relaxSummary(const RelaxProduct& p) { return relaxText(p); }
std::string relaxSettings(const MeshRelaxRequest& r) { return requestText(r); }
} // namespace knit::lab

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: stitchmesh_lab surface|analyze|field|fixtures|build|load|relax-fixtures|relax|component|pieces|snap|marked|import|propose|markings|arm|rebind|propose-marks|propose-setups|prep-compare|smooth-fixture|shaping-fixtures|yarn-fixtures|yarn-embed|yarn-settle-fixtures|yarn-settle|ncp-replay|production-patch [--key value ...]\n");
        return 2;
    }
    const std::string cmd = argv[1];
    const Args a = parseArgs(argc, argv);
    if (cmd == "surface") return commandSurface(a);
    if (cmd == "analyze") return commandAnalyze(a);
    if (cmd == "field") return commandField(a);
    if (cmd == "fixtures") return commandFixtures(a);
    if (cmd == "build") return commandBuild(a);
    if (cmd == "load") return commandLoad(a);
    if (cmd == "relax") return commandRelax(a);
    if (cmd == "relax-fixtures") return commandRelaxFixtures(a);
    if (cmd == "component") return commandComponent(a);
    if (cmd == "pieces") return commandPieces(a);
    if (cmd == "snap") return commandSnap(a);
    if (cmd == "marked") return commandMarked(a);
    if (cmd == "relax-crease") return commandRelaxCrease(a);
    if (const int rc = smoothCommand(cmd, argc, argv); rc != 2) return rc;
    if (const int rc = yarnCommand(cmd, argc, argv); rc != 2) return rc;
    if (const int rc = settleCommand(cmd, argc, argv); rc != 2) return rc;
    return armCommand(cmd, argc, argv);
}
