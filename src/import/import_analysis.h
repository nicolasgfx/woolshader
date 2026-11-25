#pragma once
// What a file is, read just after Open (the Import step's stage 2): the sanity check in three groups and its nodes. Nothing
// here changes the mesh except finishNodes, which fills its node table. Defects are counted by the same rules the repairs
// apply (repair_rules.h), on the file as read.
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "io/mesh_read.h"
#include "model/surface_pieces.h"
#include "model/vec.h"

namespace knit {

// The sanity check, in three groups (the Import dialog's sections).
struct ImportIssue {
    enum class Kind { Fatal, Repairable, Feature, Warning } kind = Kind::Warning;
    std::string id;       // stable identifier (tests, repair toggles)
    std::string text;     // what was found
    std::string action;   // what the user can do / what the repair changes
    long long count = 0;
};

struct ImportAnalysis {
    std::vector<ImportIssue> issues;
    Vec3 lo, hi;                                 // file units, file axes
    double area = 0;                             // surface area, file units^2 (the stitch estimate)
    int components = 0, boundaryLoops = 0, boundaryEdges = 0, nonManifoldEdges = 0, inconsistentEdges = 0;
    long long degenerate = 0, duplicateFaces = 0, unreferenced = 0, duplicateVertices = 0;
    std::string colourSource;                    // "texture bunny.png", "vertex colours", "material colours (Kd)", "none"
    double suggestedUnitMm = 1;                  // from the bounding box; a suggestion only
    std::string suggestionReason;
    bool fatal() const;
};

// The pass after reading (Assimp's post-processing): a file naming fewer than two nodes gets its connected pieces on the welded mesh as nodes, unless there
// is only one piece; then every node's triangles and area.
void finishNodes(ImportedMesh& m);

ImportAnalysis analyzeImport(const ImportedMesh& m);

// The yarn colour source (material.yarn_colour_source) a new object opens with: a texture is knitted
// ("texture"), other colours in the file are one colour per part ("object"), else one colour ("single").
std::string defaultYarnColourSource(const ImportAnalysis& a);

// The file as the Mark step and the guided Knit route read it: the mesh with its components named (finishNodes), its pieces
// (model/surface_pieces.h, at pieceWeldTolerance) and its surfaceDigest, the identity a markings file records. Both read it
// here, so a piece index and a digest mean the same file and the same piece in both. False with a message when the file
// cannot be read.
struct SurfaceFile {
    std::shared_ptr<const ImportedMesh> mesh;
    std::shared_ptr<const SurfacePieces> pieces;
    uint64_t digest = 0;
};
[[nodiscard]] bool readSurfaceFile(const std::string& path, SurfaceFile& out, std::string* err);

} // namespace knit
