#include "stitchmesh/polygon_mesh.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <utility>

#include "model/hash.h"

namespace knit {

// --- PolygonMesh

void PolygonMesh::addFace(const std::vector<int>& vertices) {
    corner.insert(corner.end(), vertices.begin(), vertices.end());
    faceStart.push_back(static_cast<int>(corner.size()));
}

Vec3 PolygonMesh::faceNormal(int f) const {
    Vec3 n;
    const int s = sides(f);
    for (int k = 0; k < s; ++k) n += cross(pos[vertexAt(f, k)], pos[vertexAt(f, (k + 1) % s)]);
    return normalized(n);
}

Vec3 PolygonMesh::faceCentroid(int f) const {
    Vec3 c;
    const int s = sides(f);
    for (int k = 0; k < s; ++k) c += pos[vertexAt(f, k)];
    return c / s;
}

double PolygonMesh::faceArea(int f) const {
    Vec3 n;
    const int s = sides(f);
    for (int k = 0; k < s; ++k) n += cross(pos[vertexAt(f, k)], pos[vertexAt(f, (k + 1) % s)]);
    return 0.5 * norm(n);
}

uint64_t PolygonMesh::digest() const {
    const uint64_t sizes[3] = {pos.size(), faceStart.size(), corner.size()};
    uint64_t h = fnv1a64Bytes(sizes, sizeof sizes);
    h = fnv1a64Bytes(pos.data(), pos.size() * sizeof(Vec3), h);
    h = fnv1a64Bytes(faceStart.data(), faceStart.size() * sizeof(int), h);
    return fnv1a64Bytes(corner.data(), corner.size() * sizeof(int), h);
}

// --- HalfEdgeMesh

bool buildHalfEdges(const PolygonMesh& m, HalfEdgeMesh& out, std::string* why) {
    out = HalfEdgeMesh{};
    const int corners = static_cast<int>(m.corner.size());
    out.next.assign(corners, -1);
    out.twin.assign(corners, -1);
    out.origin.assign(corners, -1);
    out.face.assign(corners, -1);
    out.halfEdgeOfVertex.assign(m.pos.size(), -1);
    out.boundaryVertex.assign(m.pos.size(), 0);
    std::map<std::pair<int, int>, int> directed;  // (origin, destination) -> half-edge
    int inconsistent = 0, nonManifold = 0;
    for (int f = 0; f < m.faces(); ++f) {
        const int s = m.sides(f);
        for (int k = 0; k < s; ++k) {
            const int h = m.faceStart[f] + k;
            out.next[h] = m.faceStart[f] + (k + 1) % s;
            out.origin[h] = m.corner[h];
            out.face[h] = f;
            if (out.halfEdgeOfVertex[m.corner[h]] < 0) out.halfEdgeOfVertex[m.corner[h]] = h;
            const std::pair<int, int> key(m.corner[h], m.corner[m.faceStart[f] + (k + 1) % s]);
            if (directed.count(key)) ++inconsistent;  // the same direction twice: the two faces disagree
            else directed[key] = h;
        }
    }
    for (const auto& [key, h] : directed) {
        const auto t = directed.find({key.second, key.first});
        if (t == directed.end()) continue;
        out.twin[h] = t->second;
    }
    // three or more half-edges on one undirected edge
    std::map<std::pair<int, int>, int> uses;
    for (int h = 0; h < corners; ++h) {
        const int a = out.origin[h], b = out.destination(h);
        ++uses[{std::min(a, b), std::max(a, b)}];
    }
    for (const auto& [edge, n] : uses)
        if (n > 2) ++nonManifold;
    if (inconsistent || nonManifold) {
        if (why) *why = std::to_string(nonManifold) + " edges with more than two faces, " + std::to_string(inconsistent) + " used twice in one direction";
        return false;
    }
    // boundary loops: each boundary half-edge's successor along the boundary leaves its destination
    std::map<int, int> boundaryOut;  // vertex -> the boundary half-edge leaving it
    for (int h = 0; h < corners; ++h) {
        if (out.twin[h] >= 0) continue;
        boundaryOut[out.origin[h]] = h;
        out.boundaryVertex[out.origin[h]] = 1;
        out.boundaryVertex[out.destination(h)] = 1;
    }
    std::vector<char> seen(corners, 0);
    for (int h = 0; h < corners; ++h) {
        if (out.twin[h] >= 0 || seen[h]) continue;
        ++out.boundaryLoops;
        int g = h;
        for (int guard = 0; guard < corners && !seen[g]; ++guard) {
            seen[g] = 1;
            const auto it = boundaryOut.find(out.destination(g));
            if (it == boundaryOut.end()) break;
            g = it->second;
        }
    }
    return true;
}

// --- checkPolygonMesh

std::string PolygonMeshCheck::text() const {
    char b[512];
    std::snprintf(b, sizeof b,
                  "polygon mesh %s: %d vertices, %d edges, %d faces (%d triangles, %d quads, %d pentagons, %d larger); "
                  "boundary edges %d in %d loops, non-manifold %d, inconsistent %d, degenerate %d; %d components, Euler %d",
                  toString(status), vertices, edges, faces, triangles, quads, pentagons, larger, boundaryEdges, boundaryLoops,
                  nonManifoldEdges, inconsistentEdges, degenerateFaces, components, euler);
    std::string t = b;
    if (!detail.empty()) t += "; " + detail;
    return t;
}

PolygonMeshCheck checkPolygonMesh(const PolygonMesh& m) {
    PolygonMeshCheck c;
    c.faces = m.faces();
    if (c.faces == 0) {
        c.status = Status::Fail;
        c.detail = "no faces";
        return c;
    }
    std::vector<char> used(m.pos.size(), 0);
    std::map<std::pair<int, int>, int> undirected;
    std::map<std::pair<int, int>, int> directed;
    double areaSum = 0;
    std::vector<double> areas(c.faces);
    for (int f = 0; f < c.faces; ++f) {
        const int s = m.sides(f);
        if (s == 3) ++c.triangles;
        else if (s == 4) ++c.quads;
        else if (s == 5) ++c.pentagons;
        else ++c.larger;
        std::set<int> distinct;
        for (int k = 0; k < s; ++k) {
            const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % s);
            used[a] = 1;
            distinct.insert(a);
            ++undirected[{std::min(a, b), std::max(a, b)}];
            if (directed.count({a, b})) ++c.inconsistentEdges;
            else directed[{a, b}] = f;
        }
        if (static_cast<int>(distinct.size()) < s) ++c.degenerateFaces;
        areas[f] = m.faceArea(f);
        areaSum += areas[f];
    }
    const double meanArea = areaSum / c.faces;
    for (int f = 0; f < c.faces; ++f)
        if (areas[f] < 1e-6 * meanArea) ++c.degenerateFaces;
    for (char u : used) c.vertices += u;
    c.edges = static_cast<int>(undirected.size());
    for (const auto& [edge, n] : undirected) {
        if (n == 1) ++c.boundaryEdges;
        if (n > 2) ++c.nonManifoldEdges;
    }
    // components by union-find over the edges
    std::vector<int> parent(m.pos.size());
    for (size_t v = 0; v < parent.size(); ++v) parent[v] = static_cast<int>(v);
    auto root = [&](int v) {
        while (parent[v] != v) v = parent[v] = parent[parent[v]];
        return v;
    };
    for (const auto& [edge, n] : undirected) parent[root(edge.first)] = root(edge.second);
    std::set<int> roots;
    for (size_t v = 0; v < parent.size(); ++v)
        if (used[v]) roots.insert(root(static_cast<int>(v)));
    c.components = static_cast<int>(roots.size());
    c.euler = c.vertices - c.edges + c.faces;
    HalfEdgeMesh he;
    std::string why;
    if (buildHalfEdges(m, he, &why)) c.boundaryLoops = he.boundaryLoops;
    c.status = (c.nonManifoldEdges == 0 && c.inconsistentEdges == 0 && c.degenerateFaces == 0) ? Status::Pass : Status::Fail;
    if (c.status == Status::Fail) c.detail = why.empty() ? std::to_string(c.degenerateFaces) + " degenerate faces" : why;
    return c;
}

// --- normalizePolygons

namespace {

double cornerAngleDeg(const std::vector<int>& f, size_t k, const std::vector<Vec3>& pos) {
    const size_t n = f.size();
    const Vec3 a = pos[f[(k + n - 1) % n]] - pos[f[k]], b = pos[f[(k + 1) % n]] - pos[f[k]];
    const double la = norm(a), lb = norm(b);
    if (la <= 0 || lb <= 0) return 0;
    return std::acos(std::clamp(dot(a, b) / (la * lb), -1.0, 1.0)) * 180 / kPi;
}

struct SplitChoice {
    bool found = false;
    size_t j = 0, k = 0;
    double alignment = -1;
};

// The diagonal of f most aligned with the field at its first vertex; existing edges and adjacent corners are no diagonals.
SplitChoice bestDiagonal(const std::vector<int>& f, const std::vector<Vec3>& pos, const std::vector<Vec3>& axis, const std::vector<Vec3>& normal,
                         const std::set<std::pair<int, int>>& edges) {
    SplitChoice best;
    const size_t n = f.size();
    for (size_t j = 0; j < n; ++j)
        for (size_t k = j + 2; k < n; ++k) {
            if (j == 0 && k == n - 1) continue;  // adjacent corners
            const int a = f[j], b = f[k];
            if (edges.count({std::min(a, b), std::max(a, b)})) continue;
            const Vec3 nj = normal[a];
            Vec3 d = pos[b] - pos[a];
            d -= nj * dot(nj, d);
            if (!(norm(d) > 0)) continue;
            d = normalized(d);
            Vec3 q = axis[a];
            q -= nj * dot(nj, q);
            if (!(norm(q) > 0)) continue;
            q = normalized(q);
            const double s = std::max(std::fabs(dot(d, q)), std::fabs(dot(d, cross(nj, q))));
            if (s > best.alignment) {
                best.found = true;
                best.j = j;
                best.k = k;
                best.alignment = s;
            }
        }
    return best;
}

} // namespace

PolygonMesh normalizePolygons(const PolygonMesh& m, const std::vector<Vec3>& axis, const std::vector<Vec3>& normal, NormalizeLog& log) {
    log = NormalizeLog{};
    PolygonMesh out;
    out.pos = m.pos;
    std::set<std::pair<int, int>> edges;
    for (int f = 0; f < m.faces(); ++f)
        for (int k = 0; k < m.sides(f); ++k) {
            const int a = m.vertexAt(f, k), b = m.vertexAt(f, (k + 1) % m.sides(f));
            edges.insert({std::min(a, b), std::max(a, b)});
        }
    char b[200];
    for (int f = 0; f < m.faces(); ++f) {
        std::vector<int> face(m.corner.begin() + m.faceStart[f], m.corner.begin() + m.faceStart[f + 1]);
        // (2) polygons of more than five sides: split recursively along the best diagonal
        std::vector<std::vector<int>> work{face}, pieces;
        while (!work.empty()) {
            std::vector<int> p = std::move(work.back());
            work.pop_back();
            if (p.size() <= 5) {
                pieces.push_back(std::move(p));
                continue;
            }
            const SplitChoice d = bestDiagonal(p, m.pos, axis, normal, edges);
            if (!d.found) {
                std::snprintf(b, sizeof b, "face %d (%zu sides): no usable diagonal, kept", f, p.size());
                log.reasons.push_back(b);
                pieces.push_back(std::move(p));
                continue;
            }
            std::snprintf(b, sizeof b, "face %d (%zu sides): split along %d-%d, field alignment %.2f", f, p.size(), p[d.j], p[d.k], d.alignment);
            log.reasons.push_back(b);
            ++log.largeSplits;
            edges.insert({std::min(p[d.j], p[d.k]), std::max(p[d.j], p[d.k])});
            std::vector<int> first(p.begin() + d.j, p.begin() + d.k + 1), second(p.begin() + d.k, p.end());
            second.insert(second.end(), p.begin(), p.begin() + d.j + 1);
            work.push_back(std::move(second));
            work.push_back(std::move(first));
        }
        // pentagons: three triangles from the widest corner; quads with a reflex corner: two triangles from it
        for (const std::vector<int>& p : pieces) {
            if (p.size() == 4) {
                // the turn at every corner about the face's reference normal (the corners' vertex normals, else Newell)
                Vec3 n{0, 0, 0};
                for (int v : p) if (v < static_cast<int>(normal.size())) n += normal[v];
                if (norm(n) < 1e-9) n = cross(m.pos[p[1]] - m.pos[p[0]], m.pos[p[2]] - m.pos[p[0]]) + cross(m.pos[p[2]] - m.pos[p[0]], m.pos[p[3]] - m.pos[p[0]]);
                int reflex = -1;
                double worst = 0;
                for (int k = 0; k < 4; ++k) {
                    const Vec3& a = m.pos[p[(k + 3) % 4]];
                    const Vec3& c = m.pos[p[k]];
                    const Vec3& d = m.pos[p[(k + 1) % 4]];
                    const double turn = dot(cross(c - a, d - c), n);
                    if (turn <= 0 && turn <= worst) { worst = turn; reflex = k; }
                }
                if (reflex >= 0) {
                    const int v[4] = {p[reflex], p[(reflex + 1) % 4], p[(reflex + 2) % 4], p[(reflex + 3) % 4]};
                    out.addFace({v[0], v[1], v[2]});
                    out.addFace({v[0], v[2], v[3]});
                    std::snprintf(b, sizeof b, "face %d (quad): a reflex corner at vertex %d, two triangles from it", f, v[0]);
                    log.reasons.push_back(b);
                    ++log.reflexQuads;
                    continue;
                }
            }
            if (p.size() != 5) {
                out.addFace(p);
                continue;
            }
            size_t widest = 0;
            double widestDeg = -1;
            for (size_t k = 0; k < 5; ++k) {
                const double a = cornerAngleDeg(p, k, m.pos);
                if (a > widestDeg) widestDeg = a, widest = k;
            }
            const int v[5] = {p[widest], p[(widest + 1) % 5], p[(widest + 2) % 5], p[(widest + 3) % 5], p[(widest + 4) % 5]};
            out.addFace({v[0], v[1], v[2]});
            out.addFace({v[0], v[2], v[3]});
            out.addFace({v[0], v[3], v[4]});
            std::snprintf(b, sizeof b, "face %d (pentagon): three triangles from its widest corner, vertex %d (%.0f degrees)", f, v[0], widestDeg);
            log.reasons.push_back(b);
            ++log.pentagons;
        }
    }
    return out;
}

} // namespace knit
