#pragma once
// The coarse mesh of the stitch-mesh layout: polygons of three or more sides on the reference surface, its half-edge
// connectivity, its validity, and the normalization of what the remesher extracts into triangles and quads
// (Wu et al. 2018, section 4, after the extraction). No file IO, no solver: records and the operations on them.
#include <cstdint>
#include <string>
#include <vector>

#include "model/status.h"
#include "model/vec.h"

namespace knit {

// Faces as one corner list: face f has the corners [faceStart[f], faceStart[f + 1]), counter-clockwise seen from
// outside the surface. faceStart has faces() + 1 entries (one entry, 0, for no face).
struct PolygonMesh {
    std::vector<Vec3> pos;
    std::vector<int> faceStart{0};
    std::vector<int> corner;

    int faces() const { return static_cast<int>(faceStart.size()) - 1; }
    int sides(int f) const { return faceStart[f + 1] - faceStart[f]; }
    int vertexAt(int f, int k) const { return corner[faceStart[f] + k]; }
    void addFace(const std::vector<int>& vertices);
    Vec3 faceNormal(int f) const;    // Newell's method, unit (zero for a degenerate face)
    Vec3 faceCentroid(int f) const;  // the mean of the corners
    double faceArea(int f) const;
    // FNV-1a 64 over the three array sizes, then the positions bit for bit, then faceStart and corner: two meshes with
    // the same corners partitioned into different faces differ.
    uint64_t digest() const;
};

// One half-edge per corner: h runs from corner h to the next corner of its face.
struct HalfEdgeMesh {
    std::vector<int> next, twin, origin, face;  // twin -1: a boundary half-edge
    std::vector<int> halfEdgeOfVertex;          // one outgoing half-edge per vertex (-1: unused vertex)
    std::vector<char> boundaryVertex;           // per vertex: some half-edge out of or into it has no twin
    int boundaryLoops = 0;
    int destination(int h) const { return origin[next[h]]; }
};
// False, with the reason, when an edge is used by more than two half-edges or twice in the same direction.
bool buildHalfEdges(const PolygonMesh& m, HalfEdgeMesh& out, std::string* why);

// Is the polygon mesh something a later stage can read at all?
struct PolygonMeshCheck {
    Status status = Status::Unresolved;  // PASS: manifold edges, one consistent orientation, no degenerate face
    int vertices = 0, edges = 0, faces = 0, triangles = 0, quads = 0, pentagons = 0, larger = 0;
    int boundaryEdges = 0, nonManifoldEdges = 0, inconsistentEdges = 0, boundaryLoops = 0, components = 0, euler = 0;
    int degenerateFaces = 0;  // a repeated vertex, or an area below 1e-6 of the mean face area
    std::string detail;
    std::string text() const;
};
PolygonMeshCheck checkPolygonMesh(const PolygonMesh& m);

// Wu 2018 section 4, after the extraction: every polygon with more than five sides is split along the diagonal
// most aligned with the field at the diagonal's first vertex (the alignment max(|d . q|, |d . n x q|) of the
// diagonal's direction d in that vertex's tangent plane; ties go to the lowest corner pair), recursively; then
// every pentagon is split from its widest corner (the T-junction of the extraction) into three triangles.
// axis and normal are per vertex of m. Each operation is logged with its reason.
struct NormalizeLog {
    int largeSplits = 0, pentagons = 0;
    int reflexQuads = 0;   // quads with a corner turning the other way (a dart), split from that corner into two triangles:
                           // the subdivision of stage 4 would make an inverted cell there (the car's coarse face 923)
    std::vector<std::string> reasons;
};
PolygonMesh normalizePolygons(const PolygonMesh& m, const std::vector<Vec3>& axis, const std::vector<Vec3>& normal, NormalizeLog& log);

} // namespace knit
