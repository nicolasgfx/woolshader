#pragma once
// Which parts of a GLB file a scene takes, and at which size (io/glb_load.h reads the file). Apart from
// the reader, so the scene registry describes a scene without the mesh and image types.
#include <string>
#include <vector>

namespace knit {

// What to take from a file, and how large it is.
struct GlbSelection {
    // One scene-wide similarity: sizeMm is the bbox measure `sizeRef` of ALL meshes, after Y-up -> Z-up.
    enum class SizeRef { Height, HorizontalExtent } sizeRef = SizeRef::Height;
    double sizeMm = 0;
    // The mesh selection, as an import's (import.plain_nodes): the glTF nodes, by name, kept as the file has
    // them (plain: shown and exported as their mesh, never knitted). Decided first; a preset knits every other
    // node (knitAllMaterials).
    std::vector<std::string> plainNodes;
    // Per primitive, by its material name: knitted, plain or skipped. A material in none of the lists is an
    // error (no silent default).
    std::vector<std::string> knitMaterials, plainMaterials, skipMaterials;
    bool knitAllMaterials = false;  // every material is knitted
    // Per glTF mesh name, before the material lists: plain (as above) or skipped. For scenes whose
    // meshes share one material.
    std::vector<std::string> plainMeshes, skipMeshes;
};

// A preset's components as the Import step lists them: each glTF node with a mesh, its triangles and
// whether the selection knits it (glb_load.cpp; the file's JSON only, no buffer is read). False with *err
// for a file that is no glTF, or a material the selection does not name.
struct GlbNode {
    std::string name;
    long long triangles = 0;
    bool knitted = true;
};
bool listGlbNodes(const std::string& path, const GlbSelection& selection, std::vector<GlbNode>& out, std::string* err);

} // namespace knit
