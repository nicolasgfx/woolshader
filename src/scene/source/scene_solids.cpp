// The scene's closed components and their inside tests, made when the source loads (scene_solids.h).
#include "scene/source/scene_solids.h"

#include <map>
#include "model/mesh_query.h"

namespace knit {

TriMesh extractTriangles(const TriMesh& mesh, const std::vector<int>& triangles, Vec3& centre) {
    TriMesh out;
    std::map<int, int> newIndex;
    for (int t : triangles) {
        std::array<int, 3> corners;
        for (int k = 0; k < 3; ++k) {
            auto it = newIndex.find(mesh.tri[t][k]);
            if (it == newIndex.end()) {
                it = newIndex.insert({mesh.tri[t][k], static_cast<int>(out.pos.size())}).first;
                out.pos.push_back(mesh.pos[mesh.tri[t][k]]);
            }
            corners[k] = it->second;
        }
        out.tri.push_back(corners);
        if (!mesh.uv.empty()) out.uv.push_back(mesh.uv[t]);
    }
    centre = out.centroid();
    for (Vec3& p : out.pos) p = p - centre;
    out.computeNormals();
    return out;
}

SceneSolids::SceneSolids(const SourceScene& source) {
    for (size_t pi = 0; pi < source.prims.size(); ++pi) {
        const SourcePrimitive& prim = source.prims[pi];
        if (prim.cls == SourcePrimitive::Class::Skip) continue;
        for (size_t k = 0; k < prim.check.components.size(); ++k) {
            if (!prim.check.components[k].insideEligible) continue;
            Vec3 centre;
            auto mesh = std::make_unique<TriMesh>(extractTriangles(prim.welded, prim.check.components[k].triangles, centre));
            for (Vec3& p : mesh->pos) p = p + centre;  // back in library coordinates; not the source's bits
            mesh->computeNormals();
            auto query = std::make_unique<MeshQuery>(*mesh);
            queries_.push_back(query.get());
            solids_.push_back({static_cast<int>(pi), static_cast<int>(k), std::move(mesh), std::move(query)});
        }
    }
}

LatticeSolids SceneSolids::otherThan(const ScenePart& part) const {
    LatticeSolids others;
    for (const Solid& solid : solids_)
        if (!(solid.prim == part.prim && solid.component == part.component)) others.meshes.push_back(solid.query.get());
    others.offset = part.centre;
    return others;
}

bool SceneSolids::insideAny(const Vec3& pointMm) const {
    for (const MeshQuery* query : queries_)
        if (query->inside(pointMm)) return true;
    return false;
}

} // namespace knit
