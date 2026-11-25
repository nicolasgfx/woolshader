// stitchmesh_lab, the imported-arm commands of the stockinette experiment:
// the original model through the general import, saved markings on one component's prepared reference, and Method A's
// typed stitch mesh and its relaxation from those markings.
//
//   stitchmesh_lab import --obj <converted source .obj> --units <mm per file unit> --out <dir> [--component <id> --name <label>]
//       readMeshFile, finishNodes, analyzeImport, then prepareAndStore with the file's own surface (KeepFile) into
//       <dir>/import_cache; the stored entry is read back (loadImport) and every prepared triangle is traced to its file
//       triangle and node. Components = (node, material, welded connected piece). Writes import_report.txt,
//       components.csv, prepared_components.obj (one group per component), triangle_map.csv; with --component also
//       <dir>/reference/: reference.obj, reference.txt (identities, digest, topology), reference_map.csv, and every other
//       component of the same node as neighbour_<id>.obj (hidden-region and attachment analysis only).
//   stitchmesh_lab propose --ref <reference dir> --out <markings file> [--side inner|outer] [--w 5.5 --h 3.4]
//       automatic proposals from the reference geometry: start point, finish boundary, growth direction, row-start
//       path, attachment regions, saved as a marking record (provenance automatic_proposal)
//   stitchmesh_lab markings --ref <reference dir> --markings <file> --out <dir> [--w 5.5 --h 3.4]
//       reads the record back, validates it against the reference and the component, and previews its effect: the time
//       field from the start, the rings, their crossings with the row-start path, counts, hidden fractions
//   stitchmesh_lab arm --ref <reference dir> --markings <file> --out <dir> [--w 5.5 --h 3.4]
//                      [--variant baseline|adjacent] [--relax anchored|pole|free] [--filter-normals facet|smooth]
//                      [--solver-faces virtual|stored] [--timestep 1]
//       Method A from the markings: typed mesh, checks, virtual shift, coverage, cells; then the relaxation on the
//       virtual faces and the same measures after it (before/ and after/ products kept apart). variant adjacent: the
//       adjacent-shaping construction; relax anchored: the pole pinned, the bind-off ring on its level, every row's first
//       vertex on the row-start path (pole: only the pole and the bind-off ring; free: the CP3 solver as in CP3);
//       filter-normals smooth: solver note R1; solver-faces stored: a diagnostic without the helix
//   stitchmesh_lab rebind --ref <old reference dir> --new-ref <new reference dir> --markings <file> --out <file> [--tolerance 1e-6]
//       a checked closest-point rebinding of the markings onto another triangulation of the same component
//   stitchmesh_lab propose-marks --mesh <file> --piece <i> --out <file.markings> [--unit-mm 68.1 --w 5.5 --h 3.4]
//       the lab's nine proposals of a closed piece hidden at one end (stitchmesh/mark_proposals.h, labProposalFile)
//   stitchmesh_lab propose-setups --mesh <file> --piece <i> [--out <file.markings>] [--unit-mm 1 --w 5.5 --h 3.4]
//                                 [--preparation facets|loop_limit] [--refine <setup>]
//       the app's proposals for one piece (proposeMarkings): every setup offered, and a report of the setups and findings
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "import/import_analysis.h"
#include "import/import_pipeline.h"
#include "import/import_store.h"
#include "import/prepare_import.h"
#include "io/cache_paths.h"
#include "io/mesh_read.h"
#include "model/hash.h"
#include "model/mesh.h"
#include "model/mesh_check.h"
#include "model/mesh_query.h"
#include "stitchmesh/coverage.h"
#include "stitchmesh/guided_rows.h"
#include "jobs/marking_session.h"
#include "io/marking_file.h"
#include "model/surface_path.h"
#include "model/surface_pieces.h"
#include "stitchmesh/mark_proposals.h"
#include "stitchmesh/marked_build.h"
#include "stitchmesh/prepared_reference.h"
#include "stitchmesh/reference_markings.h"
#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/row_rings.h"
#include "stitchmesh/time_field.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"
#include "stitchmesh_lab_shared.h"

using namespace knit;
namespace fs = std::filesystem;

namespace {

// The historical reference-bound records (knit::refmark) by their short names inside this file; the Mark step's
// types are always written knit::SurfacePoint, knit::MarkingSet, ... here.
using knit::refmark::commonTriangle;
using knit::refmark::ComponentIdentity;
using knit::refmark::densifyPath;
using knit::refmark::Marking;
using knit::refmark::MarkingCheck;
using knit::refmark::MarkingSet;
using knit::refmark::MarkingType;
using knit::refmark::pointPosition;
using knit::refmark::Provenance;
using knit::refmark::readMarkings;
using knit::refmark::rebindMarkings;
using knit::refmark::RebindReport;
using knit::refmark::referenceDigest;
using knit::refmark::ReferenceIdentity;
using knit::refmark::SurfacePoint;
using knit::refmark::validateMarkings;
using knit::refmark::writeMarkings;
using knit::refmark::toString;
using knit::toString;

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

void putf(std::ostream& o, const char* fmt, ...) {
    char b[4096];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(b, sizeof b, fmt, ap);
    va_end(ap);
    o << b;
}

std::string hex64(uint64_t v) {
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(v));
    return b;
}

std::string contentDigest(const std::string& path) {
    // FNV-1a 64 of the file's bytes: the lab's content identity of the converted source (the conversion report holds the
    // FBX's sha256). Named sha256 in the record only where a sha256 was computed by the converter.
    std::ifstream f(path, std::ios::binary);
    std::string bytes((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    return "fnv1a64:" + hex64(fnv1a64Bytes(bytes.data(), bytes.size()));
}

bool writeTriObj(const std::string& path, const TriMesh& m, const char* header) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    f << "# " << header << "\n";
    for (const Vec3& p : m.pos) putf(f, "v %.17g %.17g %.17g\n", p.x, p.y, p.z);
    for (const auto& t : m.tri) putf(f, "f %d %d %d\n", t[0] + 1, t[1] + 1, t[2] + 1);
    return static_cast<bool>(f);
}

bool readTriObj(const std::string& path, TriMesh& m) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return false;
    m = TriMesh{};
    std::string line;
    while (std::getline(f, line)) {
        std::istringstream in(line);
        std::string w;
        in >> w;
        if (w == "v") {
            Vec3 p;
            in >> p.x >> p.y >> p.z;
            m.pos.push_back(p);
        } else if (w == "f") {
            std::array<int, 3> t{};
            for (int& x : t) {
                std::string tok;
                in >> tok;
                x = std::stoi(tok.substr(0, tok.find('/'))) - 1;
            }
            m.tri.push_back(t);
        }
    }
    m.computeNormals();
    return !m.tri.empty();
}

// ---------------------------------------------------------------------------------------------------------------------
// import

struct ImportedComponent {
    int id = -1;
    int node = -1, material = -1, piece = -1;
    std::vector<int> triangles;          // prepared triangle ids, ascending
    int firstFileTriangle = -1;
    double area = 0;
    Vec3 centroid, lo, hi;
    int vertices = 0, edges = 0, boundaryEdges = 0, nonManifoldEdges = 0, euler = 0;
};

struct ImportResult {
    std::shared_ptr<ImportedMesh> file;
    ImportAnalysis analysis;
    ImportOptions options;
    uint64_t identity = 0;
    std::string storedDir;
    PreparedMesh prepared;
    std::vector<std::string> changes;
    std::vector<int> fileTriangleOf;     // per prepared triangle
    double maxCorrespondence = 0;        // [file units] centroid distance of the matched file triangle
    int unmatched = 0;
    std::vector<ImportedComponent> components;
};

bool runImport(const Args& a, ImportResult& r, std::string* err) {
    const std::string obj = a.get("--obj"), out = a.get("--out");
    if (obj.empty() || out.empty()) {
        if (err) *err = "import: --obj and --out are required";
        return false;
    }
    const std::string cache = fs::absolute(fs::path(out) / "import_cache").string();
    fs::create_directories(cache);
#ifdef _WIN32
    _putenv_s("KNIT_IMPORT_CACHE", cache.c_str());
#else
    setenv("KNIT_IMPORT_CACHE", cache.c_str(), 1);
#endif
    r.file = std::make_shared<ImportedMesh>();
    if (!readMeshFile(obj, ImportLimits{}, *r.file, err)) return false;
    finishNodes(*r.file);
    r.analysis = analyzeImport(*r.file);
    r.options.unitMm = a.number("--units", 68.1);
    r.options.up = a.get("--up", "y") == "z" ? UpAxis::Z : UpAxis::Y;
    r.options.surface = ImportOptions::Surface::KeepFile;
    ImportRequest req;
    req.file = r.file;
    req.options = r.options;
    req.identity = importIdentity(*r.file, r.options);
    r.identity = req.identity;
    if (!prepareAndStore(req, r.changes, err, {})) return false;
    r.storedDir = (fs::path(importCacheRoot()) / identityHex(r.identity)).string();
    std::string src;
    if (!loadImport(r.storedDir, r.prepared, &src, err)) return false;
    // every prepared triangle back to its file triangle: the nearest file-triangle centroid, in file coordinates
    const ImportedMesh& f = *r.file;
    const int nF = static_cast<int>(f.tri.size());
    std::vector<Vec3> fc(nF);
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (int t = 0; t < nF; ++t) {
        fc[t] = (f.pos[f.tri[t][0]] + f.pos[f.tri[t][1]] + f.pos[f.tri[t][2]]) / 3.0;
        lo = Vec3{std::min(lo.x, fc[t].x), std::min(lo.y, fc[t].y), std::min(lo.z, fc[t].z)};
        hi = Vec3{std::max(hi.x, fc[t].x), std::max(hi.y, fc[t].y), std::max(hi.z, fc[t].z)};
    }
    const double cell = std::max(1e-9, norm(hi - lo) / 400.0);
    std::map<std::array<long long, 3>, std::vector<int>> grid;
    auto key = [&](const Vec3& p) { return std::array<long long, 3>{(long long) std::floor((p.x - lo.x) / cell), (long long) std::floor((p.y - lo.y) / cell), (long long) std::floor((p.z - lo.z) / cell)}; };
    for (int t = 0; t < nF; ++t) grid[key(fc[t])].push_back(t);
    const TriMesh& pm = r.prepared.mesh;
    r.fileTriangleOf.assign(pm.tri.size(), -1);
    for (size_t t = 0; t < pm.tri.size(); ++t) {
        Vec3 c{0, 0, 0};
        for (int k = 0; k < 3; ++k) c += toFileCoordinates(r.prepared, pm.pos[pm.tri[t][k]]);
        c = c / 3.0;
        const auto k0 = key(c);
        double best = 1e300;
        int bestT = -1;
        for (long long dx = -1; dx <= 1; ++dx)
            for (long long dy = -1; dy <= 1; ++dy)
                for (long long dz = -1; dz <= 1; ++dz) {
                    const auto it = grid.find({k0[0] + dx, k0[1] + dy, k0[2] + dz});
                    if (it == grid.end()) continue;
                    for (int ft : it->second) {
                        const double d = norm(fc[ft] - c);
                        if (d < best) {
                            best = d;
                            bestT = ft;
                        }
                    }
                }
        if (bestT < 0 || best > 1e-4) {
            ++r.unmatched;
            continue;
        }
        r.fileTriangleOf[t] = bestT;
        r.maxCorrespondence = std::max(r.maxCorrespondence, best);
    }
    // components: (node, material, connected piece over shared prepared vertices)
    std::map<std::pair<int, int>, std::vector<int>> groups;
    for (size_t t = 0; t < pm.tri.size(); ++t) {
        const int ft = r.fileTriangleOf[t];
        const int node = ft >= 0 ? f.nodeOfTriangle[ft] : -1;
        const int mat = t < r.prepared.triMaterial.size() ? r.prepared.triMaterial[t] : -1;
        groups[{node, mat}].push_back(static_cast<int>(t));
    }
    for (const auto& [nm, tris] : groups) {
        // union-find over vertices
        std::map<int, int> parent;
        std::function<int(int)> findRoot = [&](int v) {
            while (parent[v] != v) v = parent[v] = parent[parent[v]];
            return v;
        };
        for (int t : tris)
            for (int k = 0; k < 3; ++k) parent.emplace(pm.tri[t][k], pm.tri[t][k]);
        for (int t : tris)
            for (int k = 1; k < 3; ++k) {
                const int a0 = findRoot(pm.tri[t][0]), b0 = findRoot(pm.tri[t][k]);
                if (a0 != b0) parent[std::max(a0, b0)] = std::min(a0, b0);
            }
        std::map<int, std::vector<int>> pieces;
        for (int t : tris) pieces[findRoot(pm.tri[t][0])].push_back(t);
        std::vector<ImportedComponent> made;
        for (auto& [root, list] : pieces) {
            ImportedComponent c;
            c.node = nm.first;
            c.material = nm.second;
            c.triangles = list;
            c.firstFileTriangle = 1 << 30;
            for (int t : list) c.firstFileTriangle = std::min(c.firstFileTriangle, r.fileTriangleOf[t] < 0 ? (1 << 30) : r.fileTriangleOf[t]);
            made.push_back(c);
        }
        std::sort(made.begin(), made.end(), [](const ImportedComponent& x, const ImportedComponent& y) { return x.firstFileTriangle < y.firstFileTriangle; });
        for (size_t i = 0; i < made.size(); ++i) {
            made[i].piece = static_cast<int>(i);
            r.components.push_back(made[i]);
        }
    }
    std::sort(r.components.begin(), r.components.end(), [](const ImportedComponent& x, const ImportedComponent& y) { return x.firstFileTriangle < y.firstFileTriangle; });
    for (size_t i = 0; i < r.components.size(); ++i) {
        ImportedComponent& c = r.components[i];
        c.id = static_cast<int>(i);
        std::set<int> verts;
        std::map<std::pair<int, int>, int> edgeUse;
        Vec3 acc{0, 0, 0};
        c.lo = Vec3{1e300, 1e300, 1e300};
        c.hi = Vec3{-1e300, -1e300, -1e300};
        for (int t : c.triangles) {
            const auto& tr = pm.tri[t];
            const Vec3 p0 = pm.pos[tr[0]], p1 = pm.pos[tr[1]], p2 = pm.pos[tr[2]];
            const double ar = 0.5 * norm(cross(p1 - p0, p2 - p0));
            c.area += ar;
            acc += ar * (p0 + p1 + p2) / 3.0;
            for (int k = 0; k < 3; ++k) {
                verts.insert(tr[k]);
                const Vec3& p = pm.pos[tr[k]];
                c.lo = Vec3{std::min(c.lo.x, p.x), std::min(c.lo.y, p.y), std::min(c.lo.z, p.z)};
                c.hi = Vec3{std::max(c.hi.x, p.x), std::max(c.hi.y, p.y), std::max(c.hi.z, p.z)};
                const int u = tr[k], v = tr[(k + 1) % 3];
                edgeUse[{std::min(u, v), std::max(u, v)}] += 1;
            }
        }
        c.centroid = c.area > 0 ? acc / c.area : acc;
        c.vertices = static_cast<int>(verts.size());
        c.edges = static_cast<int>(edgeUse.size());
        for (const auto& [e, n] : edgeUse) {
            if (n == 1) ++c.boundaryEdges;
            if (n > 2) ++c.nonManifoldEdges;
        }
        c.euler = c.vertices - c.edges + static_cast<int>(c.triangles.size());
    }
    return true;
}

std::string nodeName(const ImportResult& r, int node) { return node >= 0 && node < (int) r.file->nodes.size() ? r.file->nodes[node].name : "?"; }
std::string materialName(const ImportResult& r, int mat) {
    return mat >= 0 && mat < (int) r.prepared.materials.size() ? r.prepared.materials[mat].name : "none";
}

// The component as a mesh of its own (library mm), with the maps back to the prepared mesh.
TriMesh componentMesh(const ImportResult& r, const ImportedComponent& c, std::vector<int>* preparedVertexOf) {
    TriMesh m;
    std::map<int, int> local;
    for (int t : c.triangles) {
        std::array<int, 3> tr{};
        for (int k = 0; k < 3; ++k) {
            const int v = r.prepared.mesh.tri[t][k];
            auto it = local.find(v);
            if (it == local.end()) {
                it = local.emplace(v, static_cast<int>(m.pos.size())).first;
                m.pos.push_back(r.prepared.mesh.pos[v]);
                if (preparedVertexOf) preparedVertexOf->push_back(v);
            }
            tr[k] = it->second;
        }
        m.tri.push_back(tr);
    }
    m.computeNormals();
    return m;
}

// Midpoint refinement, `levels` times: every triangle into four, the new vertices at edge midpoints. Every new vertex lies
// on its parent triangle's plane, so the surface is unchanged; only the field's sampling gets finer. origin: per output
// triangle, its triangle in the input; cornerBary: its three corners as barycentric coordinates in that triangle.
int commandImport(const Args& a) {
    ImportResult r;
    std::string err;
    const auto t0 = std::chrono::steady_clock::now();
    if (!runImport(a, r, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const std::string out = a.get("--out");
    const std::string obj = a.get("--obj");
    {
        std::ofstream o(out + "/import_report.txt", std::ios::binary);
        putf(o, "source %s (%s, %llu bytes, content %s)\n", obj.c_str(), r.file->format.c_str(), (unsigned long long) r.file->fileBytes, contentDigest(obj).c_str());
        putf(o, "read: %zu vertices, %zu triangles (%d faces as written, %d with more than 3 corners), %zu materials, %zu nodes\n", r.file->pos.size(), r.file->tri.size(),
             r.file->faces, r.file->polygons, r.file->materials.size(), r.file->nodes.size());
        for (const ImportNode& n : r.file->nodes) putf(o, "  node %-12s %8lld triangles, area %.6f file units^2\n", n.name.c_str(), n.triangles, n.area);
        putf(o, "analysis: %d connected components, %d boundary loops (%d edges), %d non-manifold edges, %d inconsistent edges, %lld degenerate, %lld duplicate faces, "
                "%lld unreferenced, %lld duplicate vertices; suggested unit %.4g mm (%s)\n",
             r.analysis.components, r.analysis.boundaryLoops, r.analysis.boundaryEdges, r.analysis.nonManifoldEdges, r.analysis.inconsistentEdges, r.analysis.degenerate,
             r.analysis.duplicateFaces, r.analysis.unreferenced, r.analysis.duplicateVertices, r.analysis.suggestedUnitMm, r.analysis.suggestionReason.c_str());
        for (const ImportIssue& i : r.analysis.issues) putf(o, "  issue %s (%lld): %s / %s\n", i.id.c_str(), i.count, i.text.c_str(), i.action.c_str());
        o << "options (importOptionsText):\n" << importOptionsText(r.options);
        putf(o, "import identity %s; stored entry %s\n", identityHex(r.identity).c_str(), r.storedDir.c_str());
        o << "preparation changes:\n";
        for (const std::string& c : r.changes) o << "  " << c << "\n";
        putf(o, "prepared (read back from mesh.kmesh): %zu vertices, %zu triangles; library = %.6g * zUp(file) + (%.6f, %.6f, %.6f) mm\n", r.prepared.mesh.pos.size(),
             r.prepared.mesh.tri.size(), r.prepared.scale, r.prepared.offset.x, r.prepared.offset.y, r.prepared.offset.z);
        putf(o, "correspondence prepared -> file triangle: %d unmatched, largest centroid distance %.3g file units\n", r.unmatched, r.maxCorrespondence);
        putf(o, "components (node, material, welded connected piece): %zu\n", r.components.size());
        putf(o, "time %.2f s\n", seconds);
    }
    {
        std::ofstream o(out + "/components.csv", std::ios::binary);
        o << "id,node,material,piece,triangles,vertices,area_mm2,cx,cy,cz,min_x,min_y,min_z,max_x,max_y,max_z,boundary_edges,nonmanifold_edges,euler,first_file_triangle\n";
        for (const ImportedComponent& c : r.components)
            putf(o, "%d,%s,%s,%d,%zu,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%d,%d,%d\n", c.id, nodeName(r, c.node).c_str(), materialName(r, c.material).c_str(),
                 c.piece, c.triangles.size(), c.vertices, c.area, c.centroid.x, c.centroid.y, c.centroid.z, c.lo.x, c.lo.y, c.lo.z, c.hi.x, c.hi.y, c.hi.z, c.boundaryEdges,
                 c.nonManifoldEdges, c.euler, c.firstFileTriangle);
    }
    {
        std::ofstream o(out + "/prepared_components.obj", std::ios::binary);
        o << "# the prepared import (library mm), one group per component c<id>\n";
        for (const Vec3& p : r.prepared.mesh.pos) putf(o, "v %.6f %.6f %.6f\n", p.x, p.y, p.z);
        for (const ImportedComponent& c : r.components) {
            putf(o, "g c%d\n", c.id);
            for (int t : c.triangles) putf(o, "f %d %d %d\n", r.prepared.mesh.tri[t][0] + 1, r.prepared.mesh.tri[t][1] + 1, r.prepared.mesh.tri[t][2] + 1);
        }
    }
    {
        std::ofstream o(out + "/triangle_map.csv", std::ios::binary);
        o << "prepared_triangle,file_triangle,node,material,component\n";
        std::vector<int> compOf(r.prepared.mesh.tri.size(), -1);
        for (const ImportedComponent& c : r.components)
            for (int t : c.triangles) compOf[t] = c.id;
        for (size_t t = 0; t < r.prepared.mesh.tri.size(); ++t) {
            const int ft = r.fileTriangleOf[t];
            putf(o, "%zu,%d,%s,%s,%d\n", t, ft, ft >= 0 ? nodeName(r, r.file->nodeOfTriangle[ft]).c_str() : "?", materialName(r, r.prepared.triMaterial[t]).c_str(), compOf[t]);
        }
    }
    std::printf("import: %zu prepared triangles, %zu components, %d unmatched, identity %s, %.2f s\n", r.prepared.mesh.tri.size(), r.components.size(), r.unmatched,
                identityHex(r.identity).c_str(), seconds);
    if (a.get("--component").empty()) return 0;
    // the chosen component's prepared reference
    const int id = std::stoi(a.get("--component"));
    if (id < 0 || id >= (int) r.components.size()) {
        std::fprintf(stderr, "no component %d\n", id);
        return 1;
    }
    const ImportedComponent& c = r.components[id];
    const std::string ref = out + "/reference";
    fs::create_directories(ref);
    std::vector<int> prepVertex;
    const TriMesh imported = componentMesh(r, c, &prepVertex);
    const int levels = static_cast<int>(a.number("--refine", 2));
    const RefinedMesh refined = knit::refineMidpoint(imported, levels);
    const TriMesh& m = refined.mesh;
    writeTriObj(ref + "/component_imported.obj", imported, "one imported component as the import prepared it, library mm (no refinement)");
    writeTriObj(ref + "/reference.obj", m, "prepared reference of one imported component: the import's triangles, midpoint-refined (surface unchanged), library mm");
    const MeshCheck mc = checkMesh(m);
    double planeDeviation = 0;   // every reference vertex from its origin triangle's plane: 0 up to rounding
    for (size_t t = 0; t < m.tri.size(); ++t) {
        const auto& p = imported.tri[refined.parent[t]];
        const Vec3 n = imported.faceNormal(refined.parent[t]);
        for (int k = 0; k < 3; ++k) planeDeviation = std::max(planeDeviation, std::fabs(dot(m.pos[m.tri[t][k]] - imported.pos[p[0]], n)));
    }
    {
        std::ofstream o(ref + "/reference_map.csv", std::ios::binary);
        o << "reference_triangle,component_triangle,prepared_triangle,file_triangle,b00,b01,b02,b10,b11,b12,b20,b21,b22\n";
        for (size_t t = 0; t < m.tri.size(); ++t) {
            const int ct = refined.parent[t];
            const auto& cb = refined.cornerBary[t];
            putf(o, "%zu,%d,%d,%d,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", t, ct, c.triangles[ct], r.fileTriangleOf[c.triangles[ct]], cb[0][0], cb[0][1], cb[0][2],
                 cb[1][0], cb[1][1], cb[1][2], cb[2][0], cb[2][1], cb[2][2]);
        }
        o << "component_vertex,prepared_vertex\n";
        for (size_t v = 0; v < prepVertex.size(); ++v) putf(o, "%zu,%d\n", v, prepVertex[v]);
    }
    {
        std::ofstream o(ref + "/reference.txt", std::ios::binary);
        putf(o, "name %s\n", a.get("--name", "component").c_str());
        putf(o, "source_obj %s\nsource_content %s\n", obj.c_str(), contentDigest(obj).c_str());
        putf(o, "import_identity %s\nstored_entry %s\n", identityHex(r.identity).c_str(), r.storedDir.c_str());
        putf(o, "node %s\nmaterial %s\npiece %d\ncomponent_id %d\n", nodeName(r, c.node).c_str(), materialName(r, c.material).c_str(), c.piece, c.id);
        putf(o, "frame library mm, Z up; library = %.6g * zUp(file) + (%.6f, %.6f, %.6f); reference coordinates = library coordinates (no component-local transform)\n",
             r.prepared.scale, r.prepared.offset.x, r.prepared.offset.y, r.prepared.offset.z);
        putf(o, "preparation import KeepFile (weld, drop degenerate and duplicate faces, compact, orient); then midpoint refinement x%d (no smoothing; largest vertex "
                "distance from its origin triangle's plane %.3g mm)\n", levels, planeDeviation);
        putf(o, "component_digest %s\ncomponent_vertices %zu\ncomponent_triangles %zu\n", referenceDigest(imported).c_str(), imported.pos.size(), imported.tri.size());
        putf(o, "digest %s\nvertices %zu\ntriangles %zu\narea_mm2 %.3f\ncomponent_area_mm2 %.3f\nboundary_edges %d\nnonmanifold_edges %d\neuler %d\n", referenceDigest(m).c_str(),
             m.pos.size(), m.tri.size(), m.area(), imported.area(), c.boundaryEdges, c.nonManifoldEdges, c.euler);
        putf(o, "mesh_check %s\n", mc.summary().c_str());
    }
    {
        std::ofstream o(ref + "/neighbours.csv", std::ios::binary);
        o << "id,node,material,piece,triangles,file\n";
        for (const ImportedComponent& d : r.components) {
            if (d.id == c.id || d.node != c.node) continue;
            const std::string name = "neighbour_" + std::to_string(d.id) + ".obj";
            writeTriObj(ref + "/" + name, componentMesh(r, d, nullptr), "another component of the same node (visibility and attachment analysis only)");
            putf(o, "%d,%s,%s,%d,%zu,%s\n", d.id, nodeName(r, d.node).c_str(), materialName(r, d.material).c_str(), d.piece, d.triangles.size(), name.c_str());
        }
    }
    std::printf("reference: component %d (%s / %s / piece %d), %zu triangles, digest %s\n", c.id, nodeName(r, c.node).c_str(), materialName(r, c.material).c_str(), c.piece,
                m.tri.size(), referenceDigest(m).c_str());
    return 0;
}

// ---------------------------------------------------------------------------------------------------------------------
// the prepared reference, its neighbours and the field from the start marking

struct Neighbour {
    int id = -1;
    std::string name;
    TriMesh mesh;
    bool closed = false;
    std::unique_ptr<MeshQuery> query;
};

struct ArmReference {
    std::string dir;
    TriMesh mesh;
    std::map<std::string, std::string> info;   // reference.txt
    ComponentIdentity component;
    ReferenceIdentity identity;
    std::vector<Neighbour> neighbours;
};

bool loadReference(const std::string& dir, ArmReference& r, std::string* err) {
    r.dir = dir;
    if (!readTriObj(dir + "/reference.obj", r.mesh)) {
        if (err) *err = "cannot read " + dir + "/reference.obj";
        return false;
    }
    std::ifstream f(dir + "/reference.txt", std::ios::binary);
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const size_t sp = line.find(' ');
        if (sp != std::string::npos) r.info[line.substr(0, sp)] = line.substr(sp + 1);
    }
    r.component.name = r.info["name"];
    for (char& ch : r.component.name)
        if (ch == ' ') ch = '_';
    r.component.sourceDigest = r.info["source_content"];
    r.component.importIdentity = r.info["import_identity"];
    r.component.node = r.info["node"];
    r.component.material = r.info["material"];
    r.component.piece = r.info.count("piece") ? std::stoi(r.info["piece"]) : -1;
    r.identity.digest = referenceDigest(r.mesh);
    r.identity.vertices = static_cast<int>(r.mesh.pos.size());
    r.identity.triangles = static_cast<int>(r.mesh.tri.size());
    r.identity.frame = r.info["frame"];
    if (r.info["digest"] != r.identity.digest) {
        if (err) *err = "reference.obj does not match the digest reference.txt records";
        return false;
    }
    std::ifstream nb(dir + "/neighbours.csv", std::ios::binary);
    std::getline(nb, line);
    while (std::getline(nb, line)) {
        std::vector<std::string> c;
        std::stringstream ss(line);
        std::string x;
        while (std::getline(ss, x, ',')) c.push_back(x);
        if (c.size() < 6) continue;
        Neighbour n;
        n.id = std::stoi(c[0]);
        n.name = c[1] + "/" + c[2] + "/" + c[3];
        if (!readTriObj(dir + "/" + c[5], n.mesh)) continue;
        std::map<std::pair<int, int>, int> use;
        for (const auto& t : n.mesh.tri)
            for (int k = 0; k < 3; ++k) use[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}] += 1;
        n.closed = true;
        for (const auto& [e, u] : use) n.closed = n.closed && u == 2;
        r.neighbours.push_back(std::move(n));
    }
    // the queries keep a reference to their mesh: made once the list no longer moves
    for (Neighbour& n : r.neighbours)
        if (n.closed) n.query = std::make_unique<MeshQuery>(n.mesh);
    return true;
}

// The closed neighbour containing p (the first in id order), or -1: a point of the arm hidden inside another component.
int hiddenIn(const ArmReference& r, const Vec3& p) {
    for (size_t i = 0; i < r.neighbours.size(); ++i)
        if (r.neighbours[i].query && r.neighbours[i].query->inside(p)) return static_cast<int>(i);
    return -1;
}

double pointValue(const TriMesh& m, const std::vector<double>& t, const SurfacePoint& p) {
    const auto& tr = m.tri[p.tri];
    return p.bary[0] * t[tr[0]] + p.bary[1] * t[tr[1]] + p.bary[2] * t[tr[2]];
}

// The row-time field from the start point: its triangle's corners seeded with their distances from the point.
bool fieldFrom(const TriMesh& m, const SurfacePoint& start, std::vector<double>& t, TimeFieldReport& rep, std::string* err) {
    const Vec3 x = pointPosition(m, start);
    std::vector<int> seeds;
    std::vector<double> values;
    for (int k = 0; k < 3; ++k) {
        seeds.push_back(m.tri[start.tri][k]);
        values.push_back(norm(m.pos[m.tri[start.tri][k]] - x));
    }
    return geodesicDistance(m, seeds, t, rep, err, &values);
}

SurfacePoint vertexPoint(const TriMesh& m, int v) {
    SurfacePoint p;
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
        for (int k = 0; k < 3; ++k)
            if (m.tri[t][k] == v) {
                p.tri = t;
                p.bary[0] = p.bary[1] = p.bary[2] = 0;
                p.bary[k] = 1;
                return p;
            }
    return p;
}

std::vector<std::vector<int>> neighboursOfVertices(const TriMesh& m) {
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

struct LevelHidden {
    double level = 0;
    int contours = 0;
    bool closed = false;
    double length = 0, hiddenFraction = 0;
    Vec3 centroid;
};

LevelHidden levelHidden(const ArmReference& r, const std::vector<double>& t, double level) {
    LevelHidden h;
    h.level = level;
    const RingLevel lv = extractLevel(r.mesh, t, level);
    h.contours = static_cast<int>(lv.contours.size());
    if (lv.contours.size() != 1) return h;
    const Contour& c = lv.contours.front();
    h.closed = c.closed;
    h.length = c.length;
    int hidden = 0;
    Vec3 acc{0, 0, 0};
    for (const RingPoint& p : c.points) {
        hidden += hiddenIn(r, p.p) >= 0 ? 1 : 0;
        acc += p.p;
    }
    h.hiddenFraction = c.points.empty() ? 0 : double(hidden) / c.points.size();
    h.centroid = c.points.empty() ? acc : acc / double(c.points.size());
    return h;
}

// ---------------------------------------------------------------------------------------------------------------------
// propose: automatic proposals, saved as a marking record

int commandPropose(const Args& a) {
    ArmReference r;
    std::string err;
    if (!loadReference(a.get("--ref"), r, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const std::string outFile = a.get("--out");
    const double h = a.number("--h", 3.4);
    const std::string side = a.get("--side", "inner");
    const double margin = a.number("--finish-margin", 0.5);   // x h beyond the first fully hidden level
    std::ostringstream rep;
    const TriMesh& m = r.mesh;
    // the ends: the geodesic diameter (two sweeps of the field)
    std::vector<double> t;
    TimeFieldReport tr;
    if (!geodesicDistance(m, {0}, t, tr, &err)) return 1;
    const int endA = static_cast<int>(std::max_element(t.begin(), t.end()) - t.begin());
    if (!geodesicDistance(m, {endA}, t, tr, &err)) return 1;
    const int endB = static_cast<int>(std::max_element(t.begin(), t.end()) - t.begin());
    const double diameter = t[endB];
    // which end is visible: the hidden share of the vertices within 15 mm of each end
    std::vector<double> tB;
    if (!geodesicDistance(m, {endB}, tB, tr, &err)) return 1;
    std::vector<int> hiddenBy(m.pos.size(), -1);
    for (size_t v = 0; v < m.pos.size(); ++v) hiddenBy[v] = hiddenIn(r, m.pos[v]);
    int nearA = 0, hidA = 0, nearB = 0, hidB = 0;
    for (size_t v = 0; v < m.pos.size(); ++v) {
        if (t[v] < 15) {
            ++nearA;
            hidA += hiddenBy[v] >= 0;
        }
        if (tB[v] < 15) {
            ++nearB;
            hidB += hiddenBy[v] >= 0;
        }
    }
    const double fracA = nearA ? double(hidA) / nearA : 0, fracB = nearB ? double(hidB) / nearB : 0;
    const int startVertex = fracA <= fracB ? endA : endB;
    putf(rep, "ends (geodesic diameter, %.2f mm): vertex %d hidden share near it %.3f, vertex %d hidden share %.3f -> start (visible end) = vertex %d at (%.3f, %.3f, %.3f)\n",
         diameter, endA, fracA, endB, fracB, startVertex, m.pos[startVertex].x, m.pos[startVertex].y, m.pos[startVertex].z);
    SurfacePoint start = vertexPoint(m, startVertex);
    std::string startNote = "the visible end of the geodesic diameter (the paw); rule: of the diameter's two ends, the one with the smaller hidden share within 15 mm";
    if (!a.get("--start-at").empty()) {
        // a given start position (x,y,z in library mm), put on the reference by its closest point
        Vec3 x;
        if (std::sscanf(a.get("--start-at").c_str(), "%lf,%lf,%lf", &x.x, &x.y, &x.z) != 3) {
            std::fprintf(stderr, "--start-at x,y,z\n");
            return 2;
        }
        const MeshQuery q(m);
        const MeshQuery::Hit hit = q.closest(x);
        const TrianglePoint tp = closestOnTriangle(x, m.pos[m.tri[hit.tri][0]], m.pos[m.tri[hit.tri][1]], m.pos[m.tri[hit.tri][2]]);
        start.tri = hit.tri;
        for (int k = 0; k < 3; ++k) start.bary[k] = tp.bary[k];
        const double sum = start.bary[0] + start.bary[1] + start.bary[2];
        for (int k = 0; k < 3; ++k) start.bary[k] /= sum;
        startNote = a.get("--start-note", "a given start position") + " (" + a.get("--start-at") + ", " + std::to_string(hit.dist) + " mm from the surface)";
        putf(rep, "start given at (%s): on the reference %.4f mm away, triangle %d\n", a.get("--start-at").c_str(), hit.dist, hit.tri);
    }
    std::vector<double> f;
    if (!fieldFrom(m, start, f, tr, &err)) return 1;
    const double tMax = *std::max_element(f.begin(), f.end());
    // the first level whose contour is entirely inside another component
    double fullyHidden = -1;
    int torso = -1;
    for (double L = 0.5 * h; L < tMax; L += 0.25) {
        const LevelHidden lh = levelHidden(r, f, L);
        if (lh.contours == 1 && lh.closed && lh.hiddenFraction >= 1.0) {
            fullyHidden = L;
            break;
        }
    }
    if (fullyHidden < 0) {
        putf(rep, "no level is entirely hidden: the finish needs a pole closure, which is not supported; no finish proposed\n");
        std::ofstream(outFile + ".report.txt", std::ios::binary) << rep.str();
        std::fprintf(stderr, "no fully hidden level\n");
        return 1;
    }
    const double finishLevel = std::min(fullyHidden + margin * h, tMax - h);
    const LevelHidden fin = levelHidden(r, f, finishLevel);
    // the component that hides the finish ring: the neighbour containing most of its points
    {
        const RingLevel lv = extractLevel(m, f, finishLevel);
        std::map<int, int> votes;
        for (const RingPoint& p : lv.contours.front().points) votes[hiddenIn(r, p.p)] += 1;
        int best = 0;
        for (const auto& [k, n] : votes)
            if (k >= 0 && n > best) {
                best = n;
                torso = k;
            }
    }
    putf(rep, "field from the start: max %.2f mm; first entirely hidden level %.2f mm; finish level = that + %.2f h = %.2f mm (hidden share %.3f, length %.2f mm); hidden by %s\n",
         tMax, fullyHidden, margin, finishLevel, fin.hiddenFraction, fin.length, torso >= 0 ? r.neighbours[torso].name.c_str() : "?");
    // the row-start path's end: the finish contour's point on the inner side (towards the hiding component's centroid,
    // across the arm's axis), or the opposite side
    const RingLevel lv = extractLevel(m, f, finishLevel);
    const Contour& fc = lv.contours.front();
    const Vec3 axis = fin.centroid - pointPosition(m, start);
    const Vec3 ax = axis / norm(axis);
    const Vec3 tc = r.neighbours[torso].mesh.centroid();
    Vec3 inward = tc - fin.centroid;
    inward = inward - dot(inward, ax) * ax;
    inward = inward / norm(inward);
    int pick = 0;
    double best = -1e300;
    for (size_t i = 0; i < fc.points.size(); ++i) {
        const double s = dot(fc.points[i].p - fin.centroid, inward) * (side == "outer" ? -1 : 1);
        if (s > best) {
            best = s;
            pick = static_cast<int>(i);
        }
    }
    // the vertex just above the finish level next to that point, then steepest descent of the field to the start triangle
    const RingPoint& rp = fc.points[pick];
    int v = f[rp.a] >= finishLevel ? rp.a : rp.b;
    const std::vector<std::vector<int>> nbv = neighboursOfVertices(m);
    std::set<int> seeds(m.tri[start.tri].begin(), m.tri[start.tri].end());
    std::vector<int> chain{v};
    for (int guard = 0; guard < 100000 && !seeds.count(v); ++guard) {
        int next = -1;
        double slope = 0;
        for (int w : nbv[v]) {
            const double s = (f[v] - f[w]) / norm(m.pos[v] - m.pos[w]);
            if (s > slope) {
                slope = s;
                next = w;
            }
        }
        if (next < 0) break;
        v = next;
        chain.push_back(v);
    }
    if (!seeds.count(v)) {
        std::fprintf(stderr, "the descent did not reach the start\n");
        return 1;
    }
    std::reverse(chain.begin(), chain.end());
    MarkingSet ms;
    ms.component = r.component;
    ms.reference = r.identity;
    Marking st;
    st.id = "start";
    st.type = MarkingType::StartPoint;
    st.provenance = Provenance::AutomaticProposal;
    st.points.push_back(start);
    st.settings = {{"construction", "gathered_pole_start"}, {"ring1", "auto"}};
    st.note = startNote;
    Marking fi;
    fi.id = "finish";
    fi.type = MarkingType::FinishBoundary;
    fi.provenance = Provenance::AutomaticProposal;
    // the finish point: where the path crosses the finish level
    SurfacePoint fp;
    {
        double lastT = -1;
        int lastV = -1;
        for (int c : chain) {
            if (lastV >= 0 && lastT < finishLevel && f[c] >= finishLevel) {
                const double u = (finishLevel - lastT) / (f[c] - lastT);
                const SurfacePoint pa = vertexPoint(m, lastV), pb = vertexPoint(m, c);
                const int tri = commonTriangle(m, pa, pb);
                fp.tri = tri;
                for (int k = 0; k < 3; ++k) fp.bary[k] = 0;
                for (int k = 0; k < 3; ++k) {
                    if (m.tri[tri][k] == lastV) fp.bary[k] = 1 - u;
                    if (m.tri[tri][k] == c) fp.bary[k] = u;
                }
            }
            lastT = f[c];
            lastV = c;
        }
    }
    fi.points.push_back(fp);
    char b[200];
    std::snprintf(b, sizeof b, "%.4f", finishLevel);
    fi.settings = {{"construction", "bind_off"}, {"level_rule", "first_entirely_hidden_level_plus_0.5h"}, {"proposed_level_mm", b}};
    fi.note = "the bind-off is the level contour of the row-time field through this point, inside " + (torso >= 0 ? r.neighbours[torso].name : std::string("?")) +
              "; the part of the arm beyond it stays hidden and is not knitted";
    Marking gd;
    gd.id = "growth";
    gd.type = MarkingType::GrowthDirection;
    gd.provenance = Provenance::AutomaticProposal;
    gd.settings = {{"from", "start"}, {"to", "finish"}};
    gd.note = "paw to shoulder";
    Marking pa;
    pa.id = "row_start";
    pa.type = MarkingType::RowStartPath;
    pa.provenance = Provenance::AutomaticProposal;
    pa.points.push_back(start);
    for (int c : chain) pa.points.push_back(vertexPoint(m, c));
    pa.settings = {{"side", side}};
    pa.note = "steepest descent of the row-time field from the finish ring's " + side + " point (" + side +
              " = towards the hiding component's centroid across the arm's axis, or away from it), reversed: start to finish";
    ms.markings = {st, fi, gd, pa};
    // attachment regions: triangles whose corners all lie inside one neighbour
    for (size_t k = 0; k < r.neighbours.size(); ++k) {
        Marking at;
        at.id = "attach_" + std::to_string(r.neighbours[k].id);
        at.type = MarkingType::AttachmentRegion;
        at.provenance = Provenance::AutomaticProposal;
        for (int tt = 0; tt < static_cast<int>(m.tri.size()); ++tt) {
            bool all = true;
            for (int q = 0; q < 3; ++q) all = all && hiddenBy[m.tri[tt][q]] == static_cast<int>(k);
            if (all) at.triangles.push_back(tt);
        }
        if (at.triangles.empty()) continue;
        at.settings = {{"other", "component_" + std::to_string(r.neighbours[k].id)}, {"join", "unresolved"}};
        at.note = "overlap with " + r.neighbours[k].name + ": a join candidate only, no correspondence";
        double area = 0;
        for (int tt : at.triangles) area += 0.5 * norm(cross(m.pos[m.tri[tt][1]] - m.pos[m.tri[tt][0]], m.pos[m.tri[tt][2]] - m.pos[m.tri[tt][0]]));
        putf(rep, "attachment candidate: %zu triangles (%.1f mm^2) inside %s\n", at.triangles.size(), area, r.neighbours[k].name.c_str());
        ms.markings.push_back(at);
    }
    putf(rep, "row-start path: %zu points, side %s\n", pa.points.size(), side.c_str());
    if (!writeMarkings(outFile, ms, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    std::ofstream(outFile + ".report.txt", std::ios::binary) << rep.str();
    std::printf("%s", rep.str().c_str());
    return 0;
}

// ---------------------------------------------------------------------------------------------------------------------
// markings: read back, validate, preview

struct ArmSetup {
    MarkingSet set;
    MarkingCheck check;
    std::vector<double> t;
    TimeFieldReport field;
    const Marking* start = nullptr;
    const Marking* finish = nullptr;
    const Marking* path = nullptr;
    double finishLevel = 0;
    std::vector<MarkedPathPoint> pathPoints;
    int pathDescents = 0;          // path segments along which the field decreases (backtracking)
    double pathStartGap = 0;       // [mm] first path point to the start point
};

bool loadSetup(const ArmReference& r, const std::string& file, ArmSetup& s, std::ostream& log, std::string* err) {
    if (!readMarkings(file, s.set, err)) return false;
    s.check = validateMarkings(s.set, r.mesh, r.component);
    log << s.check.text();
    if (s.check.status == Status::Fail) {
        if (err) *err = "the markings are refused (see the check)";
        return false;
    }
    s.start = s.set.find(MarkingType::StartPoint);
    s.finish = s.set.find(MarkingType::FinishBoundary);
    s.path = s.set.find(MarkingType::RowStartPath);
    if (!s.path) {
        if (err) *err = "no row_start_path: Method A needs one to place the row origin";
        return false;
    }
    if (!fieldFrom(r.mesh, s.start->points[0], s.t, s.field, err)) return false;
    s.finishLevel = pointValue(r.mesh, s.t, s.finish->points[0]);
    for (size_t i = 0; i < s.path->points.size(); ++i) {
        MarkedPathPoint p;
        p.p = pointPosition(r.mesh, s.path->points[i]);
        p.t = pointValue(r.mesh, s.t, s.path->points[i]);
        p.segmentTri = i + 1 < s.path->points.size() ? commonTriangle(r.mesh, s.path->points[i], s.path->points[i + 1]) : -1;
        if (!s.pathPoints.empty() && p.t < s.pathPoints.back().t) ++s.pathDescents;
        s.pathPoints.push_back(p);
    }
    s.pathStartGap = norm(s.pathPoints.front().p - pointPosition(r.mesh, s.start->points[0]));
    return true;
}

int commandMarkings(const Args& a) {
    ArmReference r;
    std::string err;
    if (!loadReference(a.get("--ref"), r, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const std::string out = a.get("--out");
    fs::create_directories(out);
    const double w = a.number("--w", 5.5), h = a.number("--h", 3.4);
    std::ofstream log(out + "/markings_check.txt", std::ios::binary);
    putf(log, "markings %s against reference %s (%s)\n", a.get("--markings").c_str(), r.identity.digest.c_str(), r.component.name.c_str());
    ArmSetup s;
    if (!loadSetup(r, a.get("--markings"), s, log, &err)) {
        log << "REFUSED: " << err << "\n";
        std::fprintf(stderr, "refused: %s\n", err.c_str());
        return 1;
    }
    const TriMesh& m = r.mesh;
    const double tMax = *std::max_element(s.t.begin(), s.t.end());
    const int R = std::max(1, static_cast<int>(std::lround(s.finishLevel / h)));
    const double hp = s.finishLevel / R;
    putf(log, "field from the start point: max %.3f mm, %d vertices (%d by a triangle update, %d by an edge update), %.1f ms\n", tMax, s.field.accepted, s.field.triangleUpdates,
         s.field.edgeUpdates, s.field.ms);
    putf(log, "finish level (the field at the finish point) %.4f mm -> R = %d rows, h' = %.4f mm (h %.3f)\n", s.finishLevel, R, hp, h);
    putf(log, "row-start path: %zu points, gap from the start point %.4f mm, field %.3f -> %.3f mm along it, %d segments where it decreases (backtracking)\n",
         s.pathPoints.size(), s.pathStartGap, s.pathPoints.front().t, s.pathPoints.back().t, s.pathDescents);
    {
        std::ofstream f(out + "/time_field.csv", std::ios::binary);
        f << "v,x,y,z,t,hidden_in\n";
        for (size_t v = 0; v < m.pos.size(); ++v) {
            const int hb = hiddenIn(r, m.pos[v]);
            putf(f, "%zu,%.5f,%.5f,%.5f,%.5f,%d\n", v, m.pos[v].x, m.pos[v].y, m.pos[v].z, s.t[v], hb >= 0 ? r.neighbours[hb].id : -1);
        }
    }
    {
        std::ofstream f(out + "/path.csv", std::ios::binary);
        f << "i,x,y,z,t,segment_triangle\n";
        for (size_t i = 0; i < s.pathPoints.size(); ++i)
            putf(f, "%zu,%.5f,%.5f,%.5f,%.5f,%d\n", i, s.pathPoints[i].p.x, s.pathPoints[i].p.y, s.pathPoints[i].p.z, s.pathPoints[i].t, s.pathPoints[i].segmentTri);
    }
    // the rings the markings ask for, with their crossings and hidden shares (a preview: no mesh)
    std::ofstream rf(out + "/rings_preview.csv", std::ios::binary);
    rf << "k,level,contours,closed,length,L_over_w,n,crossings,hidden_fraction,dn,mean_course_over_w\n";
    std::ofstream pf(out + "/ring_points.csv", std::ios::binary);
    pf << "k,i,x,y,z,hidden\n";
    int problems = 0, events = 0, prevN = 1;
    for (int k = 1; k <= R; ++k) {
        const double L = k * hp;
        const RingLevel lv = extractLevel(m, s.t, L);
        int crossings = 0;
        for (size_t i = 0; i + 1 < s.pathPoints.size(); ++i)
            if ((s.pathPoints[i].t - L < 0) != (s.pathPoints[i + 1].t - L < 0)) ++crossings;
        LevelHidden lh = levelHidden(r, s.t, L);
        const bool ok = lv.contours.size() == 1 && lv.contours.front().closed;
        if (!ok) ++events;
        if (crossings != 1) ++problems;
        const int n = ok ? std::max(k == 1 ? 3 : 1, static_cast<int>(std::lround(lh.length / w))) : 0;
        putf(rf, "%d,%.4f,%zu,%d,%.4f,%.4f,%d,%d,%.4f,%d,%.4f\n", k, L, lv.contours.size(), ok ? 1 : 0, lh.length, lh.length / w, n, crossings, lh.hiddenFraction, n - prevN,
             n > 0 ? lh.length / n / w : 0.0);
        if (ok)
            for (size_t i = 0; i < lv.contours.front().points.size(); ++i) {
                const Vec3& p = lv.contours.front().points[i].p;
                putf(pf, "%d,%zu,%.4f,%.4f,%.4f,%d\n", k, i, p.x, p.y, p.z, hiddenIn(r, p) >= 0 ? 1 : 0);
            }
        prevN = n;
    }
    // the part beyond the finish: area and hidden share (it is not knitted)
    double areaKnit = 0, areaBeyond = 0, hiddenBeyond = 0, hiddenKnit = 0;
    for (size_t tt = 0; tt < m.tri.size(); ++tt) {
        const auto& tr = m.tri[tt];
        const Vec3 c = (m.pos[tr[0]] + m.pos[tr[1]] + m.pos[tr[2]]) / 3.0;
        const double ar = 0.5 * norm(cross(m.pos[tr[1]] - m.pos[tr[0]], m.pos[tr[2]] - m.pos[tr[0]]));
        const double tc = (s.t[tr[0]] + s.t[tr[1]] + s.t[tr[2]]) / 3.0;
        const bool hid = hiddenIn(r, c) >= 0;
        if (tc <= s.finishLevel) {
            areaKnit += ar;
            hiddenKnit += hid ? ar : 0;
        } else {
            areaBeyond += ar;
            hiddenBeyond += hid ? ar : 0;
        }
    }
    putf(log, "rings 1..%d: %d contour events, %d rings the path does not cross exactly once\n", R, events, problems);
    putf(log, "area up to the finish level %.1f mm^2 (hidden %.1f), beyond it %.1f mm^2 (hidden %.1f, %.1f %%): not knitted\n", areaKnit, hiddenKnit, areaBeyond, hiddenBeyond,
         areaBeyond > 0 ? 100 * hiddenBeyond / areaBeyond : 0.0);
    log.close();
    std::ifstream back(out + "/markings_check.txt");
    std::cout << back.rdbuf();
    return problems + events > 0 ? 3 : 0;
}

// ---------------------------------------------------------------------------------------------------------------------
// arm: Method A from the markings, checked, relaxed, checked again

// The reference below a level: triangles split along the level contour (the same edge points extractLevel takes), so its
// border is the bind-off ring's level. parent: per output triangle, its reference triangle.
TriMesh clipBelow(const TriMesh& m, const std::vector<double>& t, double level, std::vector<int>& parent) {
    TriMesh out;
    out.pos = m.pos;
    std::map<std::pair<int, int>, int> cut;
    auto edgePoint = [&](int a, int b) {
        const std::pair<int, int> k{std::min(a, b), std::max(a, b)};
        const auto it = cut.find(k);
        if (it != cut.end()) return it->second;
        const double f = (level - t[a]) / (t[b] - t[a]);
        const int v = static_cast<int>(out.pos.size());
        out.pos.push_back((1 - f) * m.pos[a] + f * m.pos[b]);
        cut[k] = v;
        return v;
    };
    for (int tt = 0; tt < static_cast<int>(m.tri.size()); ++tt) {
        const auto& tr = m.tri[tt];
        int below = 0;
        for (int k = 0; k < 3; ++k) below += t[tr[k]] < level ? 1 : 0;
        if (below == 0) continue;
        if (below == 3) {
            out.tri.push_back(tr);
            parent.push_back(tt);
            continue;
        }
        // rotate so that corner 0 is the odd one out
        int r0 = 0;
        for (int k = 0; k < 3; ++k) {
            const bool b = t[tr[k]] < level;
            if ((below == 1 && b) || (below == 2 && !b)) r0 = k;
        }
        const int a = tr[r0], b = tr[(r0 + 1) % 3], c = tr[(r0 + 2) % 3];
        const int ab = edgePoint(a, b), ac = edgePoint(a, c);
        if (below == 1) {
            out.tri.push_back({a, ab, ac});
            parent.push_back(tt);
        } else {
            out.tri.push_back({ab, b, c});
            out.tri.push_back({ab, c, ac});
            parent.push_back(tt);
            parent.push_back(tt);
        }
    }
    // drop unused vertices
    std::vector<int> map(out.pos.size(), -1);
    TriMesh compact;
    for (auto& tr : out.tri)
        for (int& v : tr) {
            if (map[v] < 0) {
                map[v] = static_cast<int>(compact.pos.size());
                compact.pos.push_back(out.pos[v]);
            }
            v = map[v];
        }
    compact.tri = out.tri;
    compact.computeNormals();
    return compact;
}

std::vector<std::vector<int>> triBorderLoops(const TriMesh& m) {
    std::map<std::pair<int, int>, int> use;
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) use[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}] += 1;
    std::map<int, int> next;   // directed as the triangles walk it
    for (const auto& t : m.tri)
        for (int k = 0; k < 3; ++k) {
            const int x = t[k], y = t[(k + 1) % 3];
            if (use[{std::min(x, y), std::max(x, y)}] == 1) next[x] = y;
        }
    std::vector<std::vector<int>> loops;
    std::set<int> seen;
    for (const auto& [s, n0] : next) {
        if (seen.count(s)) continue;
        std::vector<int> loop;
        int v = s;
        while (!seen.count(v) && next.count(v)) {
            seen.insert(v);
            loop.push_back(v);
            v = next[v];
        }
        loops.push_back(loop);
    }
    return loops;
}

std::vector<char> polygonBorderVertices(const PolygonMesh& m) {
    std::vector<char> border(m.pos.size(), 0);
    for (const std::vector<int>& loop : polygonBorderLoops(m))
        for (int v : loop) border[v] = 1;
    return border;
}

struct CellStats {
    int measured = 0, within = 0, kCells = 0, kWithin = 0, poleCells = 0, poleWithin = 0, pentagons = 0, pentagonsWithin = 0, unresolved = 0;
    double medianK = 0, p95K = 0;   // |ratio - 1| over K widths and heights
    double widthLo = 0, widthHi = 0, heightLo = 0, heightHi = 0, angleLo = 0, angleHi = 0;
};

CellStats cellStats(const std::vector<CellMeasure>& cells, const TypedStitchMesh& s) {
    CellStats st;
    std::vector<double> dev, wr, hr, an;
    for (const CellMeasure& c : cells) {
        if (c.face < (int) s.reason.size() && s.reason[c.face] == "helix-start-step") continue;
        if (c.kind == StitchFaceKind::Unsupported) {
            ++st.unresolved;
            continue;
        }
        ++st.measured;
        st.within += c.withinBounds;
        if (c.kind == StitchFaceKind::Regular) {
            ++st.kCells;
            st.kWithin += c.withinBounds;
            dev.push_back(std::fabs(c.widthRatio - 1));
            dev.push_back(std::fabs(c.heightRatio - 1));
            wr.push_back(c.widthRatio);
            hr.push_back(c.heightRatio);
            an.push_back(c.minAngle);
            an.push_back(c.maxAngle);
        } else if (c.kind == StitchFaceKind::PoleStart) {
            ++st.poleCells;
            st.poleWithin += c.withinBounds;
        } else {
            ++st.pentagons;
            st.pentagonsWithin += c.withinBounds;
        }
    }
    auto q = [](std::vector<double> v, double p) {
        if (v.empty()) return 0.0;
        std::sort(v.begin(), v.end());
        return v[std::min(v.size() - 1, static_cast<size_t>(p * (v.size() - 1) + 0.5))];
    };
    st.medianK = q(dev, 0.5);
    st.p95K = q(dev, 0.95);
    st.widthLo = q(wr, 0);
    st.widthHi = q(wr, 1);
    st.heightLo = q(hr, 0);
    st.heightHi = q(hr, 1);
    st.angleLo = q(an, 0);
    st.angleHi = q(an, 1);
    return st;
}

std::string cellText(const char* what, const CellStats& c) {
    char b[600];
    std::snprintf(b, sizeof b,
                  "%s: cells measured %d, within frozen bounds %d (K %d/%d, pentagons %d/%d, pole start %d/%d), unresolved faces %d; K |ratio - 1| median %.4f p95 %.4f; "
                  "K width %.3f..%.3f w, height %.3f..%.3f h, wale-course angles %.1f..%.1f deg\n",
                  what, c.measured, c.within, c.kWithin, c.kCells, c.pentagonsWithin, c.pentagons, c.poleWithin, c.poleCells, c.unresolved, c.medianK, c.p95K, c.widthLo,
                  c.widthHi, c.heightLo, c.heightHi, c.angleLo, c.angleHi);
    return b;
}

// Bindings on the full reference (the relaxation binds to the clipped one).
void rebindTo(const TriMesh& ref, TypedStitchMesh& s) {
    const MeshQuery q(ref);
    s.binding.assign(s.mesh.pos.size(), VertexBinding{});
    for (size_t v = 0; v < s.mesh.pos.size(); ++v) {
        const MeshQuery::Hit h = q.closest(s.mesh.pos[v]);
        if (h.tri < 0) continue;
        const TrianglePoint tp = closestOnTriangle(s.mesh.pos[v], ref.pos[ref.tri[h.tri][0]], ref.pos[ref.tri[h.tri][1]], ref.pos[ref.tri[h.tri][2]]);
        s.binding[v].tri = h.tri;
        for (int k = 0; k < 3; ++k) s.binding[v].bary[k] = tp.bary[k];
        s.binding[v].distance = h.dist;
    }
}

// Per face: the largest distance of its centroid and edge midpoints from the reference (the geometry check's measure).
void writeDeviation(const std::string& path, const TypedStitchMesh& s, const TriMesh& ref) {
    const MeshQuery q(ref);
    std::ofstream f(path, std::ios::binary);
    f << "face,strip,slot,kind,deviation\n";
    for (int fc = 0; fc < s.mesh.faces(); ++fc) {
        double dev = q.closest(s.mesh.faceCentroid(fc)).dist;
        const int n = s.mesh.sides(fc);
        for (int k = 0; k < n; ++k) dev = std::max(dev, q.closest(0.5 * (s.mesh.pos[s.mesh.vertexAt(fc, k)] + s.mesh.pos[s.mesh.vertexAt(fc, (k + 1) % n)])).dist);
        putf(f, "%d,%d,%d,%s,%.4f\n", fc, fc < (int) s.strip.size() ? s.strip[fc] : -1, fc < (int) s.slot.size() ? s.slot[fc] : -1, toString(s.kind[fc]), dev);
    }
}

void writeShaping(const std::string& path, const MarkedRowsProduct& p) {
    std::ofstream f(path, std::ios::binary);
    f << "strip,kind,events,index,slots,sigma,target_fraction,target_slot,slot,moved,placed,face,reason\n";
    for (const ShapingPlacement& s : p.rows.shaping) {
        std::string reason = s.reason;
        for (char& c : reason)
            if (c == ',') c = ';';
        putf(f, "%d,%s,%d,%d,%d,%.4f,%.4f,%d,%d,%d,%d,%d,%s\n", s.strip, s.kind == StitchFaceKind::Increase ? "INC" : "DEC", s.events, s.index, s.slotsInStrip, s.sigma,
             s.targetFraction, s.targetSlot, s.slot, s.moved, s.placed ? 1 : 0, s.face, reason.c_str());
    }
}

int commandArm(const Args& a) {
    ArmReference r;
    std::string err;
    if (!loadReference(a.get("--ref"), r, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const std::string out = a.get("--out");
    fs::create_directories(out + "/before");
    fs::create_directories(out + "/after");
    const double w = a.number("--w", 5.5), h = a.number("--h", 3.4);
    std::ofstream rep(out + "/arm_report.txt", std::ios::binary);
    putf(rep, "arm: reference %s (%s), markings %s, gauge %.3g x %.3g mm, variant %s\n", r.identity.digest.c_str(), r.component.name.c_str(), a.get("--markings").c_str(), w, h,
         a.get("--variant", "baseline").c_str());
    ArmSetup s;
    if (!loadSetup(r, a.get("--markings"), s, rep, &err)) {
        rep << "REFUSED: " << err << "\n";
        std::fprintf(stderr, "refused: %s\n", err.c_str());
        return 1;
    }
    const TriMesh& m = r.mesh;
    MarkedRowsRequest rq;
    rq.w = w;
    rq.h = h;
    rq.pole = pointPosition(m, s.start->points[0]);
    rq.poleBinding.tri = s.start->points[0].tri;
    for (int k = 0; k < 3; ++k) rq.poleBinding.bary[k] = s.start->points[0].bary[k];
    rq.finishLevel = s.finishLevel;
    const std::string ring1 = s.start->setting("ring1", "auto");
    rq.ring1 = ring1 == "auto" ? 0 : std::stoi(ring1);
    rq.path = s.pathPoints;
    const std::string variant = a.get("--variant", "baseline");
    if (variant != "baseline" && variant != "adjacent") {
        std::fprintf(stderr, "unknown variant %s (baseline | adjacent)\n", variant.c_str());
        return 2;
    }
    rq.adjacentShaping = variant == "adjacent";
    MarkedRowsProduct mp;
    const auto t0 = std::chrono::steady_clock::now();
    if (!buildMarkedRows(m, s.t, rq, mp, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const double buildMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    putf(rep, "build: %s; %s (%.1f ms)\n", toString(mp.rows.status), mp.rows.detail.c_str(), buildMs);
    for (const std::string& l : mp.rows.log) rep << "  " << l << "\n";
    {
        std::ofstream f(out + "/rings.csv", std::ios::binary);
        f << "k,level,contours,length,L_over_w,n,crossings,contour_parameter,locate_error\n";
        for (size_t k = 0; k < mp.rows.rings.size(); ++k) {
            const RingSample& rs = mp.rows.rings[k];
            const PathCrossing* pc = k >= 1 && k - 1 < mp.crossings.size() ? &mp.crossings[k - 1] : nullptr;
            putf(f, "%d,%.4f,%d,%.4f,%.4f,%d,%d,%.4f,%.3g\n", rs.k, rs.level, rs.contours, rs.length, rs.length / w, rs.n, pc ? pc->crossings : 0, pc ? pc->contourParameter : -1.0,
                 pc ? pc->locateError : 0.0);
        }
    }
    if (mp.rows.status == Status::Fail) {
        putf(rep, "no mesh: %s\n", mp.rows.detail.c_str());
        std::printf("arm: FAIL, no mesh (%s)\n", mp.rows.detail.c_str());
        return 3;
    }
    writeShaping(out + "/shaping.csv", mp);
    TypedStitchMesh& stored = mp.rows.mesh;
    // the knit region of the reference (up to the bind-off level) and its border: the bind-off ring's declared loop
    std::vector<int> parent;
    const TriMesh knit = clipBelow(m, s.t, s.finishLevel, parent);
    const std::vector<std::vector<int>> loops = triBorderLoops(knit);
    writeTriObj(out + "/knit_region.obj", knit, "the reference up to the finish level (split along it), library mm");
    std::vector<DeclaredBoundary> declared;
    for (const auto& loop : loops) {
        DeclaredBoundary d;
        d.id = "bind-off";
        d.role = EdgeBoundary::BindOff;
        for (int v : loop) d.loop.push_back(knit.pos[v]);
        declared.push_back(d);
    }
    putf(rep, "knit region: %zu triangles, area %.1f mm^2 of %.1f, %zu border loop(s)\n", knit.tri.size(), knit.area(), m.area(), loops.size());
    TypedCheckOptions o;
    o.w = w;
    o.h = h;
    o.requireStrips = true;
    o.faceDeviation = 0.1 * w;
    o.bindingTolerance = 0.001 * w;
    TypedCheckOptions ov = o;
    ov.requireStrips = false;
    CoverageOptions co;
    co.coverDistance = 0.1 * w;
    co.boundaryBand = 0.5 * h;
    co.hausdorffLimit = 0.25 * h;
    VirtualShift vs;
    buildVirtualShift(stored, ShiftReading::PathSideDown, vs, &err);
    const VirtualShiftCheck vc = validateVirtualShift(stored, vs);
    const TypedStitchMesh virt = virtualTypedMesh(stored, vs);
    const TypedMeshCheck sb = checkTypedMesh(stored, knit, declared, o), vb = checkTypedMesh(virt, knit, declared, ov);
    const CoverageReport cvb = measureCoverage(knit, loops, virt.mesh, co), csb = measureCoverage(knit, loops, stored.mesh, co);
    const std::vector<CellMeasure> cellsB = measureCells(virt, w, h);
    lab::writeTypedMeshFiles(out + "/before", stored, sb, &vs, &vc);
    lab::writeCellsFile(out + "/before/cells_virtual.csv", cellsB);
    lab::writeCellsFile(out + "/before/cells_stored.csv", measureCells(stored, w, h));
    lab::writePolygonObj(out + "/before/virtual.obj", virt);
    writeDeviation(out + "/before/deviation_virtual.csv", virt, m);
    // reload: the typed files read back give the same digest
    TypedStitchMesh back;
    std::string savedDigest;
    const bool reread = lab::readTypedMeshFiles(out + "/before", back, savedDigest);
    putf(rep, "reload: %s, digest %s, stored %s -> %s\n", reread ? "read" : "NOT read", hex64(typedDigest(back)).c_str(), hex64(typedDigest(stored)).c_str(),
         reread && typedDigest(back) == typedDigest(stored) ? "equal" : "DIFFERENT");
    rep << "--- before relaxation: stored faces\n" << sb.text() << "--- before: virtual faces (the solver's)\n" << vb.text() << vc.text();
    rep << cellText("cells before (virtual faces)", cellStats(cellsB, virt));
    rep << "--- coverage before, stored tiling\n" << csb.text() << "--- coverage before, virtual faces\n" << cvb.text();
    // relaxation on the virtual faces, over the knit region; the bind-off ring slides on its level
    // the solver's faces: the virtual polygons (the helix, as in CP3); "stored" is a diagnostic only (closed rows, no helix)
    const std::string solverFaces = a.get("--solver-faces", "virtual");
    std::vector<RelaxFace> faces;
    for (int f = 0; f < stored.mesh.faces(); ++f) {
        std::vector<int> corners;
        for (int k = 0; k < stored.mesh.sides(f); ++k) corners.push_back(stored.mesh.vertexAt(f, k));
        faces.push_back(RelaxFace{solverFaces == "stored" ? corners : vs.corners[f], stored.role[f]});
    }
    putf(rep, "solver faces: %s\n", solverFaces == "stored" ? "STORED closed rows (diagnostic: no helix)" : "virtual (helix)");
    MeshRelaxRequest rr;
    rr.w = w;
    rr.h = h;
    rr.timestep = a.number("--timestep", rr.timestep);   // the frozen CP3 value unless a diagnostic run sets it
    rr.smoothNormals = a.get("--filter-normals", "facet") == "smooth";
    putf(rep, "filter normals: %s\n", rr.smoothNormals ? "smooth (interpolated vertex normals; solver note R1)" : "facet (CP3)");
    rr.maxSteps = static_cast<int>(a.number("--max-steps", rr.maxSteps));
    RelaxProduct rp;
    const std::string relaxMode = a.get("--relax", "anchored");
    if (relaxMode == "free") {
        // the CP3 solver as in CP3: the knit region as the surface, the bind-off ring sliding on its border loop, nothing else held
        relaxStitchMesh(knit, loops, faces, stored.mesh.pos, polygonBorderVertices(stored.mesh), rr, rp, &err);
    } else {
        // the markings' constraints: the pole stays at the start point, the bind-off ring slides on the finish level's contour,
        // the first vertex of every other ring slides on the row-start path; the full reference is the surface
        RelaxGuides g;
        const RingLevel fl = extractLevel(m, s.t, s.finishLevel);
        GuideLine finishLine;
        finishLine.closed = true;
        for (const RingPoint& p : fl.contours.front().points) {
            finishLine.points.push_back(p.p);
            finishLine.segmentTri.push_back(p.tri);
        }
        GuideLine pathLine;
        for (const MarkedPathPoint& p : s.pathPoints) {
            pathLine.points.push_back(p.p);
            if (p.segmentTri >= 0) pathLine.segmentTri.push_back(p.segmentTri);
        }
        g.lines = {finishLine, pathLine};
        g.lineOf.assign(stored.mesh.pos.size(), -1);
        g.pinned.assign(stored.mesh.pos.size(), 0);
        g.pinned[stored.rings[0][0]] = 1;
        const int R = static_cast<int>(stored.rings.size()) - 1;
        for (int v : stored.rings[R]) g.lineOf[v] = 0;
        if (relaxMode == "anchored")
            for (int k = 1; k < R; ++k) g.lineOf[stored.rings[k][0]] = 1;   // "pole": only the pole and the bind-off ring are held
        relaxStitchMesh(m, {}, faces, stored.mesh.pos, std::vector<char>(stored.mesh.pos.size(), 0), rr, rp, &err, &g);
        putf(rep, "anchored: %d pinned, %d guided vertices (largest snap onto a guide %.3g mm)\n", rp.pinnedVertices, rp.guidedVertices, rp.maxGuideSnap);
    }
    putf(rep, "relaxation mode %s\n", relaxMode.c_str());
    TypedStitchMesh storedAfter = stored;
    storedAfter.mesh.pos = rp.pos;
    rebindTo(m, storedAfter);
    const TypedStitchMesh virtAfter = virtualTypedMesh(storedAfter, vs);
    const TypedMeshCheck sa = checkTypedMesh(storedAfter, knit, declared, o), va = checkTypedMesh(virtAfter, knit, declared, ov);
    const VirtualShiftCheck vca = validateVirtualShift(storedAfter, vs);
    const CoverageReport cva = measureCoverage(knit, loops, virtAfter.mesh, co), csa = measureCoverage(knit, loops, storedAfter.mesh, co);
    const std::vector<CellMeasure> cellsA = measureCells(virtAfter, w, h);
    lab::writeTypedMeshFiles(out + "/after", storedAfter, sa, &vs, &vca);
    lab::writeCellsFile(out + "/after/cells_virtual.csv", cellsA);
    lab::writeCellsFile(out + "/after/cells_stored.csv", measureCells(storedAfter, w, h));
    lab::writePolygonObj(out + "/after/virtual.obj", virtAfter);
    writeDeviation(out + "/after/deviation_virtual.csv", virtAfter, m);
    {
        std::ofstream hf(out + "/after/relax_history.csv", std::ios::binary);
        hf << "step,rms_force,max_force,max_force_vertex,max_move,cg_iterations\n";
        for (const RelaxStep& st : rp.history)
            putf(hf, "%d,%.9g,%.9g,%d,%.9g,%d\n", st.step, st.rmsForce, st.maxForce, st.maxForceVertex, st.maxMove, st.cgIterations);
    }
    const double poleDrift = norm(rp.pos[0] - stored.mesh.pos[0]);
    rep << "--- relaxation (Y12 6.1 on the virtual faces, CP3 solver and settings): " << lab::relaxSettings(rr) << "\n" << lab::relaxSummary(rp);
    putf(rep, "pole vertex moved %.4f mm\n", poleDrift);
    // where the relaxed rows start: each ring's vertex 0 against the row-start path (the marking), as measured
    double worstOrigin = 0, sumOrigin = 0;
    {
        std::ofstream f(out + "/after/row_origins.csv", std::ios::binary);
        f << "ring,distance_to_path_mm\n";
        const int R = static_cast<int>(stored.rings.size()) - 1;
        for (int k = 1; k <= R; ++k) {
            const Vec3 x = rp.pos[stored.rings[k][0]];
            double best = 1e300;
            for (size_t i = 0; i + 1 < s.pathPoints.size(); ++i) {
                const ClosestSegmentPoints cp = closestSegmentPoints(x, x, s.pathPoints[i].p, s.pathPoints[i + 1].p);
                best = std::min(best, cp.distance);
            }
            putf(f, "%d,%.4f\n", k, best);
            worstOrigin = std::max(worstOrigin, best);
            sumOrigin += best;
        }
        putf(rep, "row origins after relaxation: distance from the row-start path mean %.3f mm, largest %.3f mm\n", sumOrigin / std::max(1, R), worstOrigin);
    }
    rep << "--- after relaxation: stored faces\n" << sa.text() << "--- after: virtual faces\n" << va.text() << vca.text();
    rep << cellText("cells after (virtual faces)", cellStats(cellsA, virtAfter));
    rep << "--- coverage after, stored tiling\n" << csa.text() << "--- coverage after, virtual faces\n" << cva.text();
    {
        std::ofstream f(out + "/summary.csv", std::ios::binary);
        const CellStats b0 = cellStats(cellsB, virt), b1 = cellStats(cellsA, virtAfter);
        f << "key,before,after\n";
        putf(f, "faces,%d,%d\nunsupported,%d,%d\n", sb.faces, sa.faces, sb.unsupported, sa.unsupported);
        putf(f, "connectivity_stored,%s,%s\nconnectivity_virtual,%s,%s\ngeometry_virtual,%s,%s\nvirtual_shift,%s,%s\n", toString(sb.connectivity), toString(sa.connectivity),
             toString(vb.connectivity), toString(va.connectivity), toString(vb.geometry), toString(va.geometry), toString(vc.status), toString(vca.status));
        putf(f, "coverage_stored,%s,%s\ncoverage_virtual,%s,%s\n", toString(csb.status), toString(csa.status), toString(cvb.status), toString(cva.status));
        putf(f, "uncovered_fraction_stored,%.5f,%.5f\nuncovered_fraction_virtual,%.5f,%.5f\ncovered_twice_virtual,%d,%d\nprojected_ratio_virtual,%.5f,%.5f\n", csb.uncoveredFraction,
             csa.uncoveredFraction, cvb.uncoveredFraction, cva.uncoveredFraction, cvb.coveredTwice, cva.coveredTwice, cvb.projectedRatio, cva.projectedRatio);
        putf(f, "hausdorff_stored,%.4f,%.4f\nhausdorff_virtual,%.4f,%.4f\n", csb.worstHausdorff, csa.worstHausdorff, cvb.worstHausdorff, cva.worstHausdorff);
        putf(f, "cells_measured,%d,%d\ncells_within,%d,%d\nk_cells,%d,%d\nk_within,%d,%d\npentagons,%d,%d\npentagons_within,%d,%d\npole_cells,%d,%d\npole_within,%d,%d\n", b0.measured,
             b1.measured, b0.within, b1.within, b0.kCells, b1.kCells, b0.kWithin, b1.kWithin, b0.pentagons, b1.pentagons, b0.pentagonsWithin, b1.pentagonsWithin, b0.poleCells,
             b1.poleCells, b0.poleWithin, b1.poleWithin);
        putf(f, "median_k,%.5f,%.5f\np95_k,%.5f,%.5f\nwidth_lo,%.4f,%.4f\nwidth_hi,%.4f,%.4f\nheight_lo,%.4f,%.4f\nheight_hi,%.4f,%.4f\nangle_lo,%.2f,%.2f\nangle_hi,%.2f,%.2f\n", b0.medianK,
             b1.medianK, b0.p95K, b1.p95K, b0.widthLo, b1.widthLo, b0.widthHi, b1.widthHi, b0.heightLo, b1.heightLo, b0.heightHi, b1.heightHi, b0.angleLo, b1.angleLo, b0.angleHi,
             b1.angleHi);
        putf(f, "max_face_deviation_virtual,%.4f,%.4f\ninverted_virtual,%d,%d\nintersecting_pairs_virtual,%d,%d\n", vb.maxFaceDeviation, va.maxFaceDeviation, vb.invertedFaces,
             va.invertedFaces, vb.intersectingPairs, va.intersectingPairs);
        putf(f, "relax_status,,%s\nrelax_stop,,%s\nrelax_steps,,%d\nrelax_ms,,%.1f\nrms_force,%.3e,%.3e\npole_drift_mm,,%.4f\nbuild_ms,%.1f,\n", toString(rp.status), rp.stopReason.c_str(),
             rp.steps, rp.ms, rp.rmsForceStart, rp.rmsForceEnd, poleDrift, buildMs);
        putf(f, "row_origin_mean_mm,,%.4f\nrow_origin_max_mm,,%.4f\nrelax_mode,,%s\nvariant,%s,\n", sumOrigin / std::max(1, (int) stored.rings.size() - 1), worstOrigin,
             relaxMode.c_str(), variant.c_str());
        putf(f, "rows,%d,\nh_prime,%.5f,\nfinish_level,%.4f,\nring1,%d,\nring1_auto,%d,\nstitches_total,%zu,\n", mp.rows.rows, mp.rows.hPrime, s.finishLevel,
             (int) stored.rings[1].size(), mp.ring1Auto, stored.mesh.pos.size() - 1);
    }
    std::printf("arm: build %s, %d faces (%d unsupported); virtual shift %s; relax %s (%s, %d steps, %.0f ms); stored coverage after %s; cells within after %d/%d\n",
                toString(mp.rows.status), sb.faces, sb.unsupported, toString(vc.status), toString(rp.status), rp.stopReason.c_str(), rp.steps, rp.ms, toString(csa.status),
                cellStats(cellsA, virtAfter).within, cellStats(cellsA, virtAfter).measured);
    return 0;
}

// rebind: markings made on one reference moved onto another triangulation of the same component, checked point by point
// (rebindMarkings); region markings cannot move by closest points and are dropped with a note (re-propose them).
int commandRebind(const Args& a) {
    ArmReference from, to;
    std::string err;
    if (!loadReference(a.get("--ref"), from, &err) || !loadReference(a.get("--new-ref"), to, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    MarkingSet in, kept, out;
    if (!readMarkings(a.get("--markings"), in, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const MarkingCheck before = validateMarkings(in, to.mesh, to.component);
    std::printf("the markings against the new reference, before rebinding: %s", before.text().c_str());
    kept = in;
    kept.markings.clear();
    for (const Marking& m : in.markings) {
        if (m.triangles.empty()) kept.markings.push_back(m);
        else std::printf("dropped region marking '%s' (%zu triangles): regions are not rebound by closest points; re-propose it\n", m.id.c_str(), m.triangles.size());
    }
    RebindReport rep;
    const bool ok = rebindMarkings(kept, from.mesh, to.mesh, a.number("--tolerance", 1e-6), out, rep);
    std::printf("rebinding: %s (%s)\n", ok ? "accepted" : "REFUSED", rep.detail.c_str());
    if (!ok) return 1;
    out.component = to.component;
    if (!writeMarkings(a.get("--out"), out, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    MarkingSet back;
    readMarkings(a.get("--out"), back, &err);
    const MarkingCheck after = validateMarkings(back, to.mesh, to.component);
    std::printf("the rebound markings read back against the new reference: %s", after.text().c_str());
    return after.status == Status::Fail ? 1 : 0;
}

// ---------------------------------------------------------------------------------------------------------------------
// propose-marks: automatic proposals for one piece of an imported file, as a Mark-step markings file (offered): the lab's
// nine-marking file, made by the same backend as the app's proposals (stitchmesh/mark_proposals.h, labProposalFile)

int commandProposeMarks(const Args& a) {
    MarkedFile file;
    std::string err;
    if (!readMarkedFile(a.get("--mesh"), file, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const int pieceId = static_cast<int>(a.number("--piece", -1));
    if (pieceId < 0 || pieceId >= static_cast<int>(file.pieces->pieces.size())) {
        std::fprintf(stderr, "no piece %d\n", pieceId);
        return 1;
    }
    knit::ProposalRequest rq;
    rq.pieces = file.pieces.get();
    rq.piece = pieceId;
    rq.fileDigest = file.digest;
    rq.source = a.get("--mesh");
    rq.componentName = file.mesh->nodes[file.pieces->pieces[pieceId].node].name;
    rq.mmPerUnit = a.number("--unit-mm", 68.1);
    rq.w = a.number("--w", 5.5);
    rq.h = a.number("--h", 3.4);
    knit::MarkingSet set;
    std::string report;
    if (!labProposalFile(rq, set, report, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    if (!saveMarkingFile(a.get("--out"), set, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    std::ofstream(a.get("--out") + ".report.txt", std::ios::binary) << report;
    std::printf("%s", report.c_str());
    return 0;
}

// propose-setups: the app's proposals for one piece (stitchmesh/mark_proposals.h, proposeMarkings): every setup's
// markings offered in one markings file, and a report of the shape, the setups (recipe, texts, measures, control points
// in file units) and the findings. --refine <setup>: rule 2 on that setup as well (one trial build).
int commandProposeSetups(const Args& a) {
    MarkedFile file;
    std::string err;
    if (!readMarkedFile(a.get("--mesh"), file, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    knit::ProposalRequest rq;
    rq.pieces = file.pieces.get();
    rq.piece = static_cast<int>(a.number("--piece", 0));
    rq.fileDigest = file.digest;
    rq.source = a.get("--mesh");
    rq.mmPerUnit = a.number("--unit-mm", 1);
    rq.w = a.number("--w", 5.5);
    rq.h = a.number("--h", 3.4);
    if (!parsePreparation(a.get("--preparation", "facets"), rq.preparation)) {
        std::fprintf(stderr, "unknown preparation\n");
        return 1;
    }
    knit::ProposalSet ps;
    if (!proposeMarkings(rq, ps, &err)) {
        std::fprintf(stderr, "%s\n", err.c_str());
        return 1;
    }
    const int refine = static_cast<int>(a.number("--refine", 0));
    if (refine > 0 && refine <= static_cast<int>(ps.setups.size())) {
        knit::ProposalSetup refined;
        if (refineProposalStart(rq, ps.setups[refine - 1], refined, &err)) {
            refined.setup = static_cast<int>(ps.setups.size()) + 1;
            for (knit::SurfaceMarking& m : refined.markings) m.setup = refined.setup;
            ps.setups.push_back(refined);
        } else {
            ps.log.push_back("refine " + std::to_string(refine) + " failed: " + err);
        }
    }
    std::ostringstream rep;
    putf(rep, "shape: %s (%d)\nbasis: %016llx\nseconds: %.3f\n", ps.shapeText.c_str(), static_cast<int>(ps.shape), static_cast<unsigned long long>(ps.basis), ps.seconds);
    for (const std::string& l : ps.log) rep << "log: " << l << "\n";
    const knit::SurfacePieces& pieces = *file.pieces;
    knit::MarkingSet set;
    set.source = rq.source;
    set.digest = file.digest;
    set.preparation = toString(rq.preparation);
    for (const knit::ProposalSetup& s : ps.setups) {
        putf(rep, "setup %d [%s]: %s\n  purpose: %s\n  reason: %s\n  limitation: %s\n  rows about %d, first ring %d, last ring %d, left out visible %.1f mm2, hidden %.1f mm2, partial %d, refinable %d\n",
             s.setup, s.recipe.c_str(), s.title.c_str(), s.purpose.c_str(), s.reason.c_str(), s.limitation.c_str(), s.rowsAbout, s.firstRingStitches, s.lastRingStitches,
             s.omittedVisibleMm2, s.omittedHiddenMm2, s.partial ? 1 : 0, s.refinable ? 1 : 0);
        for (const knit::SurfaceMarking& m : s.markings) {
            putf(rep, "  %s %s, %zu points:", toString(m.meaning), toString(m.shape), m.points.size());
            for (const knit::SurfacePoint& p : m.points) {
                int pi = -1;
                knit::SurfacePoint lp;
                if (!knit::pieceLocalPoint(pieces, p, pi, lp)) continue;
                const Vec3 x = knit::surfacePosition(pieces.pieces[pi].mesh, lp);
                putf(rep, " (%.4f %.4f %.4f)", x.x, x.y, x.z);
            }
            rep << "\n";
        }
    }
    for (const knit::MarkingFinding& f : ps.findings) putf(rep, "finding [level %d, cause %d, row %d, %zu positions]: %s\n", static_cast<int>(f.level), static_cast<int>(f.cause), f.row, f.at.size(), f.message.c_str());
    std::vector<knit::SurfaceMarking> all;
    for (const knit::ProposalSetup& s : ps.setups)
        for (const knit::SurfaceMarking& m : s.markings) all.push_back(m);
    knit::mergeProposals(set, pieces, rq.piece, all);
    if (!a.get("--out").empty()) {
        if (!saveMarkingFile(a.get("--out"), set, &err)) {
            std::fprintf(stderr, "%s\n", err.c_str());
            return 1;
        }
        std::ofstream(a.get("--out") + ".report.txt", std::ios::binary) << rep.str();
    }
    std::printf("%s", rep.str().c_str());
    return 0;
}

} // namespace

int armCommand(const std::string& cmd, int argc, char** argv) {
    const Args a = parseArgs(argc, argv);
    if (cmd == "import") return commandImport(a);
    if (cmd == "propose") return commandPropose(a);
    if (cmd == "markings") return commandMarkings(a);
    if (cmd == "arm") return commandArm(a);
    if (cmd == "rebind") return commandRebind(a);
    if (cmd == "propose-marks") return commandProposeMarks(a);
    if (cmd == "propose-setups") return commandProposeSetups(a);
    std::fprintf(stderr, "unknown command %s\n", cmd.c_str());
    return 2;
}
