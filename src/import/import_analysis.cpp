// The analysis shown after Open (import_analysis.h): the file's defects counted by the repairs' own rules (repair_rules.h)
// and its nodes.
#include "import/import_analysis.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <numeric>
#include <set>

#include "import/prepare_import.h"
#include "import/repair_rules.h"
#include "model/cancel.h"

namespace knit {

namespace fs = std::filesystem;

bool ImportAnalysis::fatal() const {
    return std::any_of(issues.begin(), issues.end(), [](const ImportIssue& i) { return i.kind == ImportIssue::Kind::Fatal; });
}

namespace {

void addIssue(ImportAnalysis& a, ImportIssue::Kind k, const std::string& id, const std::string& text, const std::string& action,
              long long n = 0) {
    a.issues.push_back({k, id, text, action, n});
}

// The faces as the file has them: total area, degenerate faces (a repeated corner or zero
// area), duplicate faces (the same three corners) and vertices that belong to no face.
void countFaceDefects(const ImportedMesh& m, double diag, ImportAnalysis& a) {
    std::set<std::array<int, 3>> faces;
    for (const auto& t : m.tri) {
        a.area += triangleArea(m.pos, t);
        if (isZeroArea(m.pos, t, diag)) {
            ++a.degenerate;
            continue;
        }
        if (!faces.insert(sortedCorners(t)).second) ++a.duplicateFaces;
    }
    a.unreferenced = unreferencedVertices(m.tri, m.pos.size());
}

// The vertices a weld would merge into an earlier one.
long long countDuplicateVertices(const std::vector<int>& remap) {
    long long n = 0;
    for (size_t v = 0; v < remap.size(); ++v)
        if (remap[v] != static_cast<int>(v) && remap[v] >= 0) ++n;
    return n;
}

// The triangles after welding, without those that became degenerate (the topology's mesh).
TriMesh weldedTriangles(const ImportedMesh& m, const std::vector<int>& remap) {
    TriMesh w;
    w.pos = m.pos;
    for (const auto& t : m.tri) {
        const std::array<int, 3> r{remap[t[0]], remap[t[1]], remap[t[2]]};
        if (r[0] == r[1] || r[1] == r[2] || r[0] == r[2]) continue;
        w.tri.push_back(r);
    }
    return w;
}

// Union-find root with path halving.
int findRoot(std::vector<int>& parent, int x) {
    while (parent[x] != x) x = parent[x] = parent[parent[x]];
    return x;
}

// The edges by their uses: boundary (one face), non-manifold (more than two faces) and
// inconsistently oriented (two faces using it in the same direction); the boundary edges are
// joined into loops.
void countEdges(const TriMesh& w, ImportAnalysis& a) {
    const std::vector<std::pair<uint64_t, int>> edges = sortedEdges(w.tri);
    std::vector<int> parent(w.pos.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::vector<char> onBoundary(w.pos.size(), 0);
    for (size_t i = 0; i < edges.size();) {
        if ((i & 0xffff) == 0) cancelCheckpoint();
        size_t j = i;
        int forward = 0;
        while (j < edges.size() && edges[j].first == edges[i].first) {
            const int c = edges[j].second;
            forward += w.tri[c / 3][c % 3] < w.tri[c / 3][(c % 3 + 1) % 3];
            ++j;
        }
        const size_t uses = j - i;
        if (uses == 1) {
            ++a.boundaryEdges;
            const int u = int(edges[i].first >> 32), v = int(edges[i].first & 0xffffffffu);
            const int ru = findRoot(parent, u), rv = findRoot(parent, v);
            if (ru != rv) parent[ru] = rv;
            onBoundary[u] = onBoundary[v] = 1;
        } else if (uses > 2) {
            ++a.nonManifoldEdges;
        } else if (forward != 1) {
            ++a.inconsistentEdges;
        }
        i = j;
    }
    std::set<int> loops;
    for (size_t v = 0; v < onBoundary.size(); ++v)
        if (onBoundary[v]) loops.insert(findRoot(parent, static_cast<int>(v)));
    a.boundaryLoops = static_cast<int>(loops.size());
}

// The separate parts: each triangle's piece, where triangles sharing a vertex are one piece. Pieces are numbered by their first
// triangle.
std::vector<int> connectedPieces(const std::vector<std::array<int, 3>>& tri, size_t vertexCount) {
    std::vector<int> parent(vertexCount);
    std::iota(parent.begin(), parent.end(), 0);
    for (size_t t = 0; t < tri.size(); ++t) {
        if ((t & 0xffff) == 0) cancelCheckpoint();
        const int r0 = findRoot(parent, tri[t][0]);
        for (int k = 1; k < 3; ++k) {
            const int rk = findRoot(parent, tri[t][k]);
            if (rk != r0) parent[rk] = r0;
        }
    }
    std::vector<int> pieceOfRoot(vertexCount, -1), piece(tri.size());
    int pieces = 0;
    for (size_t t = 0; t < tri.size(); ++t) {
        int& p = pieceOfRoot[findRoot(parent, tri[t][0])];
        if (p < 0) p = pieces++;
        piece[t] = p;
    }
    return piece;
}

int pieceCount(const std::vector<int>& piece) { return piece.empty() ? 0 : *std::max_element(piece.begin(), piece.end()) + 1; }

// Where the colour comes from: textures (by file name), vertex colours, material colours or none.
std::string colourSource(const ImportedMesh& m) {
    bool tex = false, kd = false;
    for (const ImportMaterial& mat : m.materials) {
        tex |= !mat.texture.empty();
        kd |= mat.hasKd;
    }
    if (tex) {
        std::string names;
        for (const ImportMaterial& mat : m.materials)
            if (!mat.texture.empty()) names += (names.empty() ? "" : ", ") + fs::path(mat.texture).filename().string();
        return "texture " + names;
    }
    if (!m.vcolor.empty()) return "vertex colours";
    if (kd) return "material colours (Kd)";
    return "none";
}

// Suggested unit: the one that puts the largest extent closest to 300 mm (a suggestion).
void suggestUnit(const Vec3& ext, ImportAnalysis& a) {
    const double E = std::max({ext.x, ext.y, ext.z});
    struct U { double mm; const char* name; };
    const U units[] = {{1, "mm"}, {10, "cm"}, {1000, "m"}, {25.4, "in"}};
    double best = 1e300;
    for (const U& u : units) {
        const double d = std::fabs(std::log(E * u.mm / 300.0));
        if (d < best) {
            best = d;
            a.suggestedUnitMm = u.mm;
            char b[160];
            std::snprintf(b, sizeof b, "the largest extent is %.4g file units: %.4g mm if the unit is %s", E, E * u.mm, u.name);
            a.suggestionReason = b;
        }
    }
}

// The defects the import can repair (each one an option of prepareImport).
void addRepairIssues(ImportAnalysis& a) {
    if (a.degenerate)
        addIssue(a, ImportIssue::Kind::Repairable, "degenerate", countText(a.degenerate) + " faces have zero area.", "Remove zero-area faces", a.degenerate);
    if (a.duplicateFaces)
        addIssue(a, ImportIssue::Kind::Repairable, "duplicate_faces", countText(a.duplicateFaces) + " faces are duplicates of others.",
                 "Remove duplicate faces", a.duplicateFaces);
    if (a.duplicateVertices)
        addIssue(a, ImportIssue::Kind::Repairable, "weld", countText(a.duplicateVertices) + " vertices coincide with others (split seams).",
                 "Merge coincident vertices (connects the surface across split seams)", a.duplicateVertices);
    if (a.unreferenced)
        addIssue(a, ImportIssue::Kind::Repairable, "unreferenced", countText(a.unreferenced) + " vertices belong to no face.", "Remove unused vertices",
                 a.unreferenced);
    if (a.inconsistentEdges)
        addIssue(a, ImportIssue::Kind::Repairable, "orientation", countText(a.inconsistentEdges) + " edges join faces that face opposite ways.",
                 "Make each part face one way (outwards for closed parts)", a.inconsistentEdges);
}

// The shape's supported features and its limits for knitting.
void addShapeIssues(const ImportedMesh& m, ImportAnalysis& a) {
    if (a.boundaryLoops)
        addIssue(a, ImportIssue::Kind::Feature, "open", std::to_string(a.boundaryLoops) + " open boundaries (holes or rims).",
                 "A rebuild closes them; with the knitting skipped, the file's surface stays as it is.");
    if (a.components > 1)
        addIssue(a, ImportIssue::Kind::Feature, "components", std::to_string(a.components) + " separate parts.",
                 "A rebuild joins the ticked parts into one closed body; with the knitting skipped, each part keeps its own surface.");
    if (m.polygons)
        addIssue(a, ImportIssue::Kind::Feature, "polygons", countText(m.polygons) + " faces with more than 3 corners.", "Supported: split into triangles.");
    if (a.nonManifoldEdges)
        addIssue(a, ImportIssue::Kind::Warning, "non_manifold", countText(a.nonManifoldEdges) + " edges are shared by more than two faces.",
                 "Such parts do not hide fibres inside them.");
}

// What the file references or contains that the import does not use.
void addFileIssues(const ImportedMesh& m, ImportAnalysis& a) {
    for (const std::string& d : m.missingDependencies)
        addIssue(a, ImportIssue::Kind::Warning, "missing_file", "Referenced file not found: " + d, "The colour falls back to the material colour.");
    if (!m.ignored.empty()) {
        std::string s;
        for (const std::string& i : m.ignored) s += (s.empty() ? "" : "; ") + i;
        addIssue(a, ImportIssue::Kind::Feature, "ignored", "Not used: " + s, "Only geometry, texture coordinates and the base colour are imported.");
    }
}

} // namespace

std::string defaultYarnColourSource(const ImportAnalysis& a) {
    if (a.colourSource.rfind("texture", 0) == 0) return "texture";
    return a.colourSource == "none" ? "single" : "object";
}

void finishNodes(ImportedMesh& m) {
    if (m.nodes.size() < 2) {
        const std::vector<int> remap = coincidentVertices(m.pos, kWeldRelative * boundingBoxOf(m.pos).diagonal);
        std::vector<std::array<int, 3>> welded = m.tri;
        for (std::array<int, 3>& t : welded)
            for (int& v : t) v = remap[v];
        std::vector<int> piece = connectedPieces(welded, m.pos.size());
        const int pieces = pieceCount(piece);
        if (pieces != 1 || m.nodes.empty()) {
            m.nodes.clear();
            for (int p = 0; p < pieces; ++p) m.nodes.push_back({"Part " + std::to_string(p + 1)});
            m.nodeOfTriangle = std::move(piece);
        }
    }
    for (size_t t = 0; t < m.tri.size(); ++t) {
        ImportNode& node = m.nodes[m.nodeOfTriangle[t]];
        ++node.triangles;
        node.area += triangleArea(m.pos, m.tri[t]);
    }
}

ImportAnalysis analyzeImport(const ImportedMesh& m) {
    ImportAnalysis a;
    const BoundingBox box = boundingBoxOf(m.pos);
    a.lo = box.lo;
    a.hi = box.hi;
    const double diag = box.diagonal;
    const Vec3 ext = a.hi - a.lo;
    if (!(diag > 0)) {
        addIssue(a, ImportIssue::Kind::Fatal, "zero_extent", "All vertices are at one point: the mesh has no size.", "Check the file in a 3D tool.");
        return a;
    }
    countFaceDefects(m, diag, a);
    if (a.degenerate == static_cast<long long>(m.tri.size())) {
        addIssue(a, ImportIssue::Kind::Fatal, "all_degenerate", "Every face has zero area.", "Check the file in a 3D tool.");
        return a;
    }
    const std::vector<int> remap = coincidentVertices(m.pos, kWeldRelative * diag);
    a.duplicateVertices = countDuplicateVertices(remap);
    // Topology on the welded, non-degenerate triangles.
    const TriMesh w = weldedTriangles(m, remap);
    cancelCheckpoint();
    countEdges(w, a);
    cancelCheckpoint();
    a.components = pieceCount(connectedPieces(w.tri, w.pos.size()));
    a.colourSource = colourSource(m);
    suggestUnit(ext, a);
    addRepairIssues(a);
    addShapeIssues(m, a);
    addFileIssues(m, a);
    return a;
}

bool readSurfaceFile(const std::string& path, SurfaceFile& out, std::string* err) {
    std::shared_ptr<ImportedMesh> mesh = std::make_shared<ImportedMesh>();
    if (!readMeshFile(path, ImportLimits{}, *mesh, err)) return false;
    finishNodes(*mesh);
    out.pieces = std::make_shared<SurfacePieces>(findSurfacePieces(mesh->pos, mesh->tri, mesh->nodeOfTriangle, pieceWeldTolerance(mesh->pos)));
    out.digest = surfaceDigest(mesh->pos, mesh->tri);
    out.mesh = mesh;
    return true;
}

} // namespace knit
