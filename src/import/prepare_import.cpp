// The preparation run by Import (prepare_import.h): the repairs in their order, orientation, library space and the identity.
#include "import/prepare_import.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <set>

#include "import/repair_rules.h"
#include "io/files.h"
#include "model/cancel.h"
#include "model/colour.h"
#include "model/coordinates.h"
#include "model/hash.h"

namespace knit {

namespace fs = std::filesystem;

namespace {

// Frozen: every identity hashes this line. It named the chart generator (xatlas), which is gone; changing it would make
// every stored import, rebuilt ones too, miss its entry.
constexpr const char* kIdentityToken = "xatlas f700c7790aaa030e794b52ba7791a05c085faf0c";
constexpr int kImporterVersion = 1;
// The rebuild's stages (wrap_stages.h); bump when a change alters wrapped output; keep-file identities
// never see it. 5: the flatten's solver runs to convergence, and the surface file holds the flattened triangles. 6: a chart
// mostly outside the stretch bounds is split.
constexpr int kWrapVersion = 6;

} // namespace

std::string plainNodesText(const std::set<std::string>& nodes) {
    std::string text;
    for (const std::string& name : nodes) text += (text.empty() ? "" : "|") + name;
    return text;
}

std::set<std::string> parsePlainNodes(const std::string& text) {
    std::set<std::string> nodes;
    for (size_t at = 0; at < text.size();) {
        const size_t end = std::min(text.find('|', at), text.size());
        nodes.insert(text.substr(at, end - at));
        at = end + 1;
    }
    return nodes;
}

std::string importOptionsText(const ImportOptions& o) {
    char b[512];
    std::snprintf(b, sizeof b,
                  "unit_mm = %.17g\nup = %s\nweld = %d\nremove_degenerate = %d\nremove_duplicate_faces = %d\nremove_unreferenced = %d\n"
                  "fix_orientation = %d\ncharts = keep\n",
                  o.unitMm, o.up == UpAxis::Y ? "y" : "z", o.weld, o.removeDegenerate, o.removeDuplicateFaces, o.removeUnreferenced,
                  o.fixOrientation);
    // "charts = keep" is frozen: charts used to be a choice (the file's uv or generated ones), and every identity holds it.
    // Only a selection writes this line, and only a wrapped import the others: every existing identity stays as it is.
    const std::string selection = o.plainNodes.empty() ? "" : "plain_nodes = " + plainNodesText(o.plainNodes) + "\n";
    if (o.surface != ImportOptions::Surface::Wrapped) return b + selection;
    // A wrapped import names no engine: the import has one engine (since 2026-11-12), so an identity written with a
    // "wrap_engine" line is never matched again and its import is rebuilt once.
    char w[128];
    std::snprintf(w, sizeof w, "surface = wrapped\nwrap_version = %d\nwrap_stitch_mm = %.17g\n", kWrapVersion, o.wrapStitchMm);
    return std::string(b) + w + selection;
}

Vec3 toFileCoordinates(const PreparedMesh& p, const Vec3& lib) { return toFileCoordinates(p.scale, p.offset, p.up, lib); }

namespace {

void flipTri(PreparedMesh& p, size_t t) {
    std::swap(p.mesh.tri[t][1], p.mesh.tri[t][2]);
    if (!p.mesh.uv.empty()) std::swap(p.mesh.uv[t][1], p.mesh.uv[t][2]);
    if (!p.colourUv.empty()) std::swap(p.colourUv[t][1], p.colourUv[t][2]);
}

// The uses of each undirected edge: sortedEdges puts them in runs; runOf[3 t + k] is the start
// of the run of triangle t's corner edge k (corners k, k + 1), runEnd[start] the run's end.
struct EdgeRuns {
    std::vector<std::pair<uint64_t, int>> edges;
    std::vector<int> runOf, runEnd;
    int uses(int t, int k) const { return runEnd[runOf[3 * t + k]] - runOf[3 * t + k]; }
};

EdgeRuns edgeRuns(const std::vector<std::array<int, 3>>& tri) {
    EdgeRuns r;
    r.edges = sortedEdges(tri);
    const std::vector<std::pair<uint64_t, int>>& edges = r.edges;
    r.runOf.resize(tri.size() * 3);
    for (size_t i = 0; i < edges.size();) {
        if ((i & 0xffff) == 0) cancelCheckpoint();
        size_t j = i;
        while (j < edges.size() && edges[j].first == edges[i].first) ++j;
        for (size_t q = i; q < j; ++q) r.runOf[edges[q].second] = static_cast<int>(i);
        i = j;
    }
    r.runEnd.assign(edges.size(), 0);
    for (size_t i = 0; i < edges.size();) {
        size_t j = i;
        while (j < edges.size() && edges[j].first == edges[i].first) ++j;
        r.runEnd[i] = static_cast<int>(j);
        i = j;
    }
    return r;
}

// +1: triangle t uses the edge u -> v; -1 otherwise.
int edgeDirection(const TriMesh& m, int t, int u, int v) {
    for (int k = 0; k < 3; ++k)
        if (m.tri[t][k] == u && m.tri[t][(k + 1) % 3] == v) return 1;
    return -1;
}

// The breadth-first walk of orient over all components: each triangle's component (its seed
// triangle; -1: not reached yet) and whether it flips. queue and members are reused per component.
struct OrientationWalk {
    std::vector<int> comp;
    std::vector<char> flip;
    std::vector<int> queue, members;
    long long visited = 0;  // triangles visited, for the cancellation checkpoint
};

// Walks the component of `seed` across its manifold edges, giving every reached triangle the
// flip that makes it use each shared edge opposite to its neighbour. false: the component
// cannot be oriented (an edge with more than two uses, or a contradiction: non-orientable).
bool orientComponent(const TriMesh& m, const EdgeRuns& runs, int seed, OrientationWalk& walk) {
    walk.queue.assign(1, seed);
    walk.members.clear();
    walk.comp[seed] = seed;
    bool ok = true;
    for (size_t q = 0; q < walk.queue.size(); ++q) {
        if ((++walk.visited & 0xffff) == 0) cancelCheckpoint();
        const int t = walk.queue[q];
        walk.members.push_back(t);
        for (int k = 0; k < 3; ++k) {
            const int n = runs.uses(t, k);
            if (n != 2) {
                if (n > 2) ok = false;
                continue;
            }
            const int r = runs.runOf[3 * t + k];
            const int c0 = runs.edges[r].second, c1 = runs.edges[r + 1].second;
            const int o = (c0 / 3 == t ? c1 : c0) / 3;
            const int u = m.tri[t][k], v = m.tri[t][(k + 1) % 3];
            // t (after its flip) uses the edge one way; o must use it the other way.
            const int dt = edgeDirection(m, t, u, v) * (walk.flip[t] ? -1 : 1);
            const char needFlip = edgeDirection(m, o, u, v) != -dt;
            if (walk.comp[o] < 0) {
                walk.comp[o] = walk.comp[seed];
                walk.flip[o] = needFlip;
                walk.queue.push_back(o);
            } else if (walk.flip[o] != needFlip) {
                ok = false;  // non-orientable
            }
        }
    }
    return ok;
}

// A closed component (every edge used twice) that is inside-out (negative signed volume after
// its flips) flips as a whole.
void turnOutwards(const TriMesh& m, const EdgeRuns& runs, OrientationWalk& walk) {
    bool closed = true;
    double vol = 0;
    for (int t : walk.members) {
        for (int k = 0; k < 3; ++k) closed = closed && runs.uses(t, k) == 2;
        std::array<int, 3> tr = m.tri[t];
        if (walk.flip[t]) std::swap(tr[1], tr[2]);
        vol += dot(m.pos[tr[0]], cross(m.pos[tr[1]], m.pos[tr[2]])) / 6.0;
    }
    if (closed && vol < 0)
        for (int t : walk.members) walk.flip[t] = !walk.flip[t];
}

// Consistent orientation per component (breadth-first over manifold edges), then outward for
// closed components (positive signed volume). Returns the number of flipped triangles; counts
// components that cannot be oriented (non-orientable or non-manifold) in `unorientable`.
long long orient(PreparedMesh& p, int* unorientable) {
    const TriMesh& m = p.mesh;
    const EdgeRuns runs = edgeRuns(m.tri);
    OrientationWalk walk;
    walk.comp.assign(m.tri.size(), -1);
    walk.flip.assign(m.tri.size(), 0);
    int bad = 0;
    for (size_t seed = 0; seed < m.tri.size(); ++seed) {
        if (walk.comp[seed] >= 0) continue;
        if (!orientComponent(m, runs, static_cast<int>(seed), walk)) {
            ++bad;
            for (int t : walk.members) walk.flip[t] = 0;
            continue;
        }
        turnOutwards(m, runs, walk);
    }
    long long flipped = 0;
    for (size_t t = 0; t < m.tri.size(); ++t)
        if (walk.flip[t]) {
            flipTri(p, t);
            ++flipped;
        }
    if (unorientable) *unorientable = bad;
    return flipped;
}

// Keeps the triangles with keep[t] set, with their material and colour coordinates.
void keepTriangles(PreparedMesh& out, const std::vector<char>& keep) {
    TriMesh& m = out.mesh;
    const bool hasUv = !out.colourUv.empty();
    size_t w = 0;
    for (size_t t = 0; t < m.tri.size(); ++t) {
        if (!keep[t]) continue;
        m.tri[w] = m.tri[t];
        out.triMaterial[w] = out.triMaterial[t];
        if (hasUv) out.colourUv[w] = out.colourUv[t];
        ++w;
    }
    m.tri.resize(w);
    out.triMaterial.resize(w);
    if (hasUv) out.colourUv.resize(w);
}

// Coincident vertices (within kWeldRelative of the diagonal) merge: the triangles use each cluster's
// first vertex; the others become unreferenced.
void mergeCoincidentVertices(PreparedMesh& out, double diag) {
    TriMesh& m = out.mesh;
    const std::vector<int> remap = coincidentVertices(m.pos, kWeldRelative * diag);
    long long merged = 0;
    for (size_t v = 0; v < remap.size(); ++v) merged += remap[v] != static_cast<int>(v);
    for (auto& t : m.tri)
        for (int& v : t) v = remap[v];
    if (merged) out.changes.push_back("merged " + countText(merged) + " coincident vertices");
}

void removeZeroAreaTriangles(PreparedMesh& out, double diag) {
    const TriMesh& m = out.mesh;
    std::vector<char> keep(m.tri.size(), 1);
    long long n = 0;
    for (size_t t = 0; t < m.tri.size(); ++t) {
        const std::array<int, 3>& tr = m.tri[t];
        if (isZeroArea(m.pos, tr, diag)) keep[t] = 0, ++n;
    }
    keepTriangles(out, keep);
    if (n) out.changes.push_back("removed " + countText(n) + " zero-area faces");
}

// A triangle with the same three corners as an earlier one (in any order) is removed.
void removeDuplicateTriangles(PreparedMesh& out) {
    const TriMesh& m = out.mesh;
    std::set<std::array<int, 3>> seen;
    std::vector<char> keep(m.tri.size(), 1);
    long long n = 0;
    for (size_t t = 0; t < m.tri.size(); ++t) {
        if (!seen.insert(sortedCorners(m.tri[t])).second) keep[t] = 0, ++n;
    }
    keepTriangles(out, keep);
    if (n) out.changes.push_back("removed " + countText(n) + " duplicate faces");
}

// Compacts the vertex array to the referenced vertices, in order of first use (welding leaves
// the merged vertices unreferenced).
void compactVertices(const ImportedMesh& in, const ImportOptions& o, PreparedMesh& out) {
    TriMesh& m = out.mesh;
    std::vector<int> map(m.pos.size(), -1);
    std::vector<Vec3> pos;
    std::vector<Vec3> col;
    long long n = 0;
    for (auto& t : m.tri)
        for (int& v : t) {
            if ((++n & 0x3ffff) == 0) cancelCheckpoint();
            if (map[v] < 0) {
                map[v] = static_cast<int>(pos.size());
                pos.push_back(m.pos[v]);
                if (!out.vcolor.empty()) col.push_back(out.vcolor[v]);
            }
            v = map[v];
        }
    const long long dropped = static_cast<long long>(m.pos.size() - pos.size());
    if (o.removeUnreferenced && dropped) {
        // Report only the vertices that were unreferenced in the file (not the merged ones).
        const long long unref = unreferencedVertices(in.tri, in.pos.size());
        if (unref) out.changes.push_back("removed " + countText(unref) + " vertices that belong to no face");
    }
    m.pos = std::move(pos);
    out.vcolor = std::move(col);
}

// Orients every part consistently (orient) and notes the flips and the parts left as they are.
void orientFaces(PreparedMesh& out) {
    int bad = 0;
    const long long f = orient(out, &bad);
    if (f) out.changes.push_back("flipped " + countText(f) + " faces to a consistent orientation");
    if (bad) out.changes.push_back(std::to_string(bad) + " parts cannot be oriented consistently (non-manifold or one-sided); left as they are");
}

} // namespace

// Library space: Z-up, mm, centred on the x-y bounds with the lowest point at z = 0 (the offset
// is kept for toFileCoordinates).
LibraryPlacement toLibrarySpace(const ImportedMesh& in, const ImportOptions& o) {
    LibraryPlacement out;
    out.pos.resize(in.pos.size());
    for (size_t v = 0; v < in.pos.size(); ++v) {
        const Vec3& q = in.pos[v];
        const Vec3 z = o.up == UpAxis::Y ? Vec3{q.x, -q.z, q.y} : q;
        out.pos[v] = z * o.unitMm;
    }
    const BoundingBox box = boundingBoxOf(out.pos);
    out.offset = {-(box.lo.x + box.hi.x) * 0.5, -(box.lo.y + box.hi.y) * 0.5, -box.lo.z};
    for (Vec3& q : out.pos) q = q + out.offset;
    out.diagonal = box.diagonal;
    return out;
}

bool prepareImport(const ImportedMesh& in, const ImportOptions& o, PreparedMesh& out, std::string* err) {
    out = PreparedMesh{};
    if (!(o.unitMm > 0) || !std::isfinite(o.unitMm)) return fail(err, "the unit size must be a positive number of millimetres");
    out.scale = o.unitMm;
    out.up = o.up;
    out.materials = in.materials;
    LibraryPlacement placed = toLibrarySpace(in, o);
    out.mesh.pos = std::move(placed.pos);
    out.offset = placed.offset;
    const double diag = placed.diagonal;
    out.mesh.tri = in.tri;
    const bool hasUv = !in.uv.empty();
    if (hasUv) out.colourUv = in.uv;
    out.triMaterial = in.triMaterial;
    out.vcolor = in.vcolor;
    cancelCheckpoint();
    // Weld first (remap), then drop what became degenerate.
    if (o.weld) mergeCoincidentVertices(out, diag);
    cancelCheckpoint();
    if (o.removeDegenerate) removeZeroAreaTriangles(out, diag);
    if (o.removeDuplicateFaces) removeDuplicateTriangles(out);
    if (out.mesh.tri.empty()) return fail(err, "no faces are left after the repairs");
    if (o.removeUnreferenced || o.weld) compactVertices(in, o, out);
    cancelCheckpoint();
    if (o.fixOrientation) orientFaces(out);
    cancelCheckpoint();
    if (hasUv) out.mesh.uv = out.colourUv;  // the file's uv are its knitting charts
    out.mesh.computeNormals();
    return true;
}

uint64_t importIdentity(const ImportedMesh& m, const ImportOptions& o) {
    std::string s = "importer " + std::to_string(kImporterVersion) + "\n" + kIdentityToken + "\n" + importOptionsText(o);
    bool ok = false;
    const uint64_t src = fnv1aFile(m.path, &ok);
    s += "source " + std::to_string(src) + (ok ? "" : " unreadable") + "\n";
    for (const std::string& d : m.dependencies) {
        const uint64_t h = fnv1aFile(d, &ok);
        s += "dependency " + fs::path(d).filename().string() + " " + std::to_string(h) + (ok ? "" : " unreadable") + "\n";
    }
    for (const std::string& d : m.missingDependencies) s += "missing " + fs::path(d).filename().string() + "\n";
    return fnv1a(s.data(), s.size());
}

std::vector<Vec3> triangleColours(const PreparedMesh& p, const std::vector<ColourImage>& textures, const Vec3& fallback) {
    std::map<std::string, const ColourImage*> byPath;
    for (const ColourImage& t : textures) byPath[t.path] = &t;
    std::vector<Vec3> out(p.mesh.tri.size(), fallback);
    for (size_t t = 0; t < p.mesh.tri.size(); ++t) {
        const int mi = t < p.triMaterial.size() ? p.triMaterial[t] : -1;
        const ImportMaterial* mat = mi >= 0 && mi < static_cast<int>(p.materials.size()) ? &p.materials[mi] : nullptr;
        if (mat && !mat->texture.empty() && !p.colourUv.empty()) {
            auto it = byPath.find(mat->texture);
            if (it != byPath.end()) {
                out[t] = triangleTextureColour(*it->second, p.colourUv[t]);
                continue;
            }
        }
        if (!p.vcolor.empty()) {
            Vec3 c{0, 0, 0};
            for (int k = 0; k < 3; ++k) {
                const Vec3& v = p.vcolor[p.mesh.tri[t][k]];
                c = c + Vec3{srgbDecode(v.x), srgbDecode(v.y), srgbDecode(v.z)};
            }
            out[t] = c * (1.0 / 3.0);
            continue;
        }
        if (mat && mat->hasKd) out[t] = mat->kd;
    }
    return out;
}

} // namespace knit
