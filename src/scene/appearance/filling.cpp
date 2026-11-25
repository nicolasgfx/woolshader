// The filling behind a knitted part; see filling.h.
#include "scene/appearance/filling.h"

#include <cstdio>

namespace knit {

Filling makeFilling(const TriMesh& mesh, const ComponentCheck& component, double offset) {
    Filling f;
    f.offset = offset;
    if (!component.consistent || component.nonManifoldEdges > 0) {
        f.reason = component.nonManifoldEdges > 0 ? "non-manifold edges" : "inconsistent orientation";
        return f;
    }
    if (component.closed && component.signedVolume <= 0) {
        f.reason = "closed with non-positive volume (inward-facing)";
        return f;
    }
    f.kind = component.closed ? Filling::Kind::Stuffing : Filling::Kind::Backing;
    TriMesh moved = mesh;
    moved.uv.clear();
    moved.computeNormals();
    for (size_t v = 0; v < moved.pos.size(); ++v) moved.pos[v] = moved.pos[v] - moved.vnormal[v] * offset;
    moved.computeNormals();
    for (size_t q = 0; q < moved.tri.size(); ++q)
        f.inverted += dot(moved.faceNormal(static_cast<int>(q)), mesh.faceNormal(static_cast<int>(q))) <= 0;
    const MeshCheck check = checkMesh(moved);
    f.selfIntersection = Status::Pass;
    for (const ComponentCheck& c : check.components) {
        f.selfIntersection = combine(f.selfIntersection, c.selfIntersection);
        f.intersectingPairs += c.intersectingPairs;
        f.unresolvedPairs += c.unresolvedPairs;
        f.signedVolume += c.signedVolume;
    }
    f.mesh = std::move(moved);
    return f;
}

std::string Filling::summary() const {
    char b[320];
    if (kind == Kind::None) {
        std::snprintf(b, sizeof b, "no filling (%s)", reason.c_str());
        return b;
    }
    std::snprintf(b, sizeof b, "%s at %.3f mm: %d of %zu moved triangles inverted, self-intersection %s (%d pairs, %d unresolved)%s",
                  kind == Kind::Stuffing ? "stuffing" : "backing", offset, inverted, mesh.tri.size(), toString(selfIntersection),
                  intersectingPairs, unresolvedPairs, kind == Kind::Stuffing && signedVolume <= 0 ? ", turned inside out (volume <= 0)" : "");
    return std::string(b) + " (rendering only, not a check of the knit)";
}

} // namespace knit
