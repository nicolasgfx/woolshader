#pragma once
// The scene's closed components, in library mm, with their inside tests. Knit and Fibres read one set
// (Scene::solids, made when the source loads): the lattice nodes and the fibres must stay out of them.
#include <memory>
#include <vector>

#include "model/mesh_query.h"
#include "knitting/layout/chart_lattice.h"
#include "scene/scene.h"

namespace knit {

class SceneSolids {
public:
    SceneSolids() = default;  // none: the sewn pieces touch by design
    explicit SceneSolids(const SourceScene& source);

    const std::vector<const MeshQuery*>& queries() const { return queries_; }
    // Every solid but that component, for its lattice, in the frame of the part's centre.
    LatticeSolids otherThan(const ScenePart& part) const;
    bool insideAny(const Vec3& pointMm) const;

private:
    struct Solid {
        int prim, component;
        std::unique_ptr<TriMesh> mesh;  // MeshQuery refers to it
        std::unique_ptr<MeshQuery> query;
    };
    std::vector<Solid> solids_;
    std::vector<const MeshQuery*> queries_;
};

// The triangles of a mesh as a compact mesh (uv kept), re-centred at their centroid.
TriMesh extractTriangles(const TriMesh& mesh, const std::vector<int>& triangles, Vec3& centre);

} // namespace knit
