#pragma once
// The prepared import (prepare_import.h) as a self-contained cache entry:
//   out/cache/imports/<identity>/mesh.kmesh     the repaired mesh in library mm, Z-up, knit uv,
//                                                   colour uv, vertex colours, materials
//                                   /plain.kmesh    the same for the nodes that stay plain (none: no file)
//                                   /textures/...   copies of the base-colour images
//                                   /import.txt     source path, options, what the repairs changed
// published through publish.h (cache policy). The identity (importIdentity) covers the importer's
// version, a frozen token, the options, the source file's bytes and every dependency's bytes, so a
// changed MTL or image is another entry; the source file itself is never written.
#include <cstdint>
#include <string>
#include <vector>

#include "import/prepare_import.h"
#include "io/glb_load.h"

namespace knit {

// Stores the prepared meshes and copies of their textures as the entry r.identity names (published as one directory). r.file
// is the source: its path, format and dependencies go into import.txt with r.options.
bool storeImport(const std::string& root, const ImportRequest& r, const PreparedImport& p, std::string* err);
// Reads a stored entry; texture paths in the materials point into the entry's directory.
bool loadImport(const std::string& dir, PreparedMesh& p, std::string* sourcePath, std::string* err);
// The stored entry exists and is complete.
bool importStored(const std::string& root, uint64_t identity);
// What preparing a stored entry changed, as PreparedMesh::changes said then.
std::vector<std::string> storedImportChanges(const std::string& root, uint64_t identity);


// A stored import as a GLB-like scene for the tiled pipeline: one knitted primitive per material
// (a mesh without materials: one), then one plain primitive per material of plain.kmesh, welded,
// checked (checkMesh); baseColor = the area-weighted
// mean base colour of the primitive's triangles; triColor per triangle (the per-part yarn colour
// averages it over the part); the colour uv and the texture path for the viewer.
bool loadImportedScene(const std::string& dir, SourceScene& out, std::string* err, const SourceOptions& o = {});

} // namespace knit
