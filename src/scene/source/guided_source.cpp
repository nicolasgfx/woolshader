// The guided route's scene: the source file's pieces, one knitted, the others plain (guided_source.h).
#include "scene/source/guided_source.h"

#include <cmath>
#include <cstdio>
#include <filesystem>

#include "import/import_analysis.h"
#include "io/marking_file.h"
#include "model/hash.h"
#include "model/mesh_check.h"
#include "scene/source/scene_solids.h"

namespace knit {

namespace {

// One piece as a primitive of the scene, in library mm; Knit for the selected piece, Plain for the context.
SourcePrimitive piecePrimitive(const SurfaceFile& file, const GuidedInputs& in, int index, SourcePrimitive::Class cls) {
    const SurfacePiece& piece = file.pieces->pieces[index];
    const ImportedMesh& mesh = *file.mesh;
    SourcePrimitive g;
    const std::string node = piece.node >= 0 && piece.node < static_cast<int>(mesh.nodes.size()) ? mesh.nodes[piece.node].name : std::string("object");
    g.node = node;
    g.mesh = node + "_piece" + std::to_string(piece.ordinal);
    g.nodeIndex = 0;  // the file is one object (as the import's primitives)
    g.cls = cls;
    const int firstTriangle = piece.fileTriangle.empty() ? -1 : piece.fileTriangle.front();
    const int material = firstTriangle >= 0 && firstTriangle < static_cast<int>(mesh.triMaterial.size()) ? mesh.triMaterial[firstTriangle] : -1;
    g.material = material >= 0 && !mesh.materials[material].name.empty() ? mesh.materials[material].name : std::string("default");
    if (material >= 0 && mesh.materials[material].hasKd) g.baseColor = mesh.materials[material].kd;
    g.welded.tri = piece.mesh.tri;
    for (const Vec3& p : piece.mesh.pos) g.welded.pos.push_back(libraryPointOfFileMm(in, p * in.mmPerUnit));
    g.welded.computeNormals();
    for (const std::array<int, 3>& t : g.welded.tri) g.cornerNormals.push_back({g.welded.vnormal[t[0]], g.welded.vnormal[t[1]], g.welded.vnormal[t[2]]});
    g.check = checkMesh(g.welded, SelfIntersectionTest::Skip);
    return g;
}

std::string hexOf(uint64_t v) {
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(v));
    return b;
}

} // namespace

Vec3 libraryPointOfFileMm(const GuidedInputs& in, const Vec3& q) {
    const Vec3 z = in.yUp ? Vec3{q.x, -q.z, q.y} : q;  // the import's axes (toLibrarySpace, import/prepare_import.cpp)
    return z + in.offset;
}

bool loadGuidedParts(const SceneSettings& s, Scene& out, GuidedInputs& in, std::string* err) {
    auto refuse = [&](const std::string& why) {
        if (err) *err = "guided stitch mesh: " + why;
        return false;
    };
    const SourceScene& placed = out.source;
    if (!placed.imported) return refuse("the scene is not an import");
    if (std::fabs(placed.scale - s.importUnitMm) > 1e-9 * std::max(1.0, s.importUnitMm)) return refuse("the import's unit differs from the settings' (import again)");
    SurfaceFile file;
    std::string why;
    if (!readSurfaceFile(s.importSource, file, &why)) return refuse("cannot read " + s.importSource + ": " + why);
    if (s.guidedMarkings.empty()) return refuse("no markings file is chosen (Mark step: Knit this piece)");
    if (!loadMarkingFile(s.guidedMarkings, in.markings, &why)) return refuse("cannot read the markings " + s.guidedMarkings + ": " + why);
    if (in.markings.digest != file.digest)
        return refuse("the markings belong to another version of " + std::filesystem::path(s.importSource).filename().string() + " (markings " + hexOf(in.markings.digest) +
                      ", file " + hexOf(file.digest) + "): they are not rebound, mark the piece again");
    const int nPieces = static_cast<int>(file.pieces->pieces.size());
    if (s.guidedPiece < 0 || s.guidedPiece >= nPieces) return refuse("piece " + std::to_string(s.guidedPiece) + " does not exist (the file has " + std::to_string(nPieces) + ")");
    if (s.guidedKnitComponent < 0) return refuse("knitting component " + std::to_string(s.guidedKnitComponent) + " does not exist");
    const uint64_t effective = effectiveMarkingsDigest(in.markings, *file.pieces, s.guidedPiece, s.guidedKnitComponent);
    if (hexOf(effective) != s.guidedMarkingsDigest)
        return refuse("the markings file changed after it was chosen for knitting (" + hexOf(effective) + " against " + s.guidedMarkingsDigest + "): choose the piece again");
    in.pieces = file.pieces;
    in.piece = s.guidedPiece;
    in.knitComponent = s.guidedKnitComponent;
    in.mmPerUnit = placed.scale;
    in.yUp = placed.importYUp;
    in.offset = placed.offset;
    const SurfacePiece& piece = file.pieces->pieces[in.piece];
    in.pieceTriangles = fnv1a64Bytes(piece.fileTriangle.data(), piece.fileTriangle.size() * sizeof(int));
    const std::string node = piece.node >= 0 && piece.node < static_cast<int>(file.mesh->nodes.size()) ? file.mesh->nodes[piece.node].name : std::string("object");
    in.pieceName = node + " piece " + std::to_string(piece.ordinal);

    // The pieces replace the prepared import's primitives: the same placement, one primitive per piece.
    out.source.prims.clear();
    out.parts.clear();
    int knitPrim = -1;
    for (int i = 0; i < nPieces; ++i) {
        if (i == in.piece) knitPrim = static_cast<int>(out.source.prims.size());
        out.source.prims.push_back(piecePrimitive(file, in, i, i == in.piece ? SourcePrimitive::Class::Knit : SourcePrimitive::Class::Plain));
    }
    out.source.log.push_back("guided stitch mesh: " + in.pieceName + " knitted from " + s.guidedMarkings + "; " + std::to_string(nPieces - 1) + " other pieces plain");
    out.solids = std::make_shared<const SceneSolids>(out.source);
    const SourcePrimitive& g = out.source.prims[knitPrim];
    ScenePart part;
    part.prim = knitPrim;
    part.component = 0;
    part.name = g.mesh;
    part.componentName = part.name;
    std::vector<int> all(g.welded.tri.size());
    for (size_t t = 0; t < all.size(); ++t) all[t] = static_cast<int>(t);
    part.mesh = std::make_shared<const TriMesh>(extractTriangles(g.welded, all, part.centre));
    part.layout = KnitLayout::Guided;
    out.parts.push_back(std::move(part));
    return true;
}

} // namespace knit
