#pragma once
// The guided route's scene (knit.route = guided, settings/scene_settings.h): an import whose one selected piece is knitted
// as a stockinette stitch mesh from the Mark step's markings, while the file's other pieces stay plain meshes as context.
//
// The pieces come from the source file as the Mark step reads it (readSurfaceFile, import/import_analysis.h), so a piece
// index and the file digest mean the same thing here and there. They are placed as the import placed the file
// (library = scale * zUp(file) + offset, the prepared import's own numbers), so the knitted piece sits exactly where the
// import showed it. Nothing is rebound by name or position: the markings file must name this very file (its digest), the
// piece must exist, and the markings must still be the ones the settings recorded (knit.guided_markings_digest).
#include <memory>
#include <string>

#include "model/surface_marking.h"
#include "model/surface_pieces.h"
#include "scene/scene.h"

namespace knit {

// What the guided Knit operation (knit/guided_part.cpp) reads beside the settings.
struct GuidedInputs {
    std::shared_ptr<const SurfacePieces> pieces;  // of import.source, file units and axes
    MarkingSet markings;                          // the markings file as read
    int piece = -1;                               // the knitted piece
    int knitComponent = 0;                        // its knitting component (SurfaceMarking::knitComponent)
    std::string pieceName;                        // "<component> piece <n>" (the Mark step's label)
    double mmPerUnit = 1;                         // the import's unit
    bool yUp = true;                              // the import's up axis
    Vec3 offset;                                  // library = mmPerUnit * zUp(file) + offset [mm]
    uint64_t pieceTriangles = 0;                  // FNV-1a 64 of the piece's file triangles (its identity in the file)
};

// A point in file axes, already in mm, as a library point [mm] (Z up).
Vec3 libraryPointOfFileMm(const GuidedInputs& in, const Vec3& fileMm);

// Replaces the prepared import's primitives in `out` (loaded by loadMeshScene) by the source file's pieces: the selected
// one a Knit primitive with one Guided part, every other one a Plain primitive. False with the reason, for people, when
// the file, the markings or the piece do not fit the settings.
bool loadGuidedParts(const SceneSettings& s, Scene& out, GuidedInputs& in, std::string* err);

} // namespace knit
