#pragma once
// The spiral layout's source (SceneEntry::knitsAsPieces: the teddy): a scene of sewn pieces, delivered
// as the common knit result (Scene, KnitLayout::Spiral) that the display, the fibres, the knit
// cache and the export consume for every scene.
//
// The source scene has one primitive per piece (knitted) and per node the scene's selection keeps (plain,
// from its GLB: eyes, noses, mouths, the bow), in library mm (Z-up) about the package origin (SpiralSource),
// each under its bear's object node. Each piece is one part (part k = piece k of the manifest): its mesh in library mm about its
// area centroid, exactly as the spiral knitting knits it, so the part's knit and its mesh agree.
#include <map>
#include <string>
#include <vector>

#include "knitting/knit/pieces.h"
#include "scene/scene.h"
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"

namespace knit {

// The pieces as the knitting places them: the manifest, the gauge sigma [mm per scene unit], the package
// origin (scene units, Z-up: the centre of all bounds over the floor) and the parts the scene's selection
// keeps (an object without pieces, the bow, goes with the bear it overlaps most). Library mm of a piece
// point p: (p - origin) sigma about the origin.
struct SpiralSource {
    std::string manifestPath;
    Manifest man;
    double sigma = 0;
    Vec3 origin;
    // The parts the scene's selection keeps (glb.plainNodes), in library mm about the origin, and the bear each
    // goes with.
    std::vector<SourcePrimitive> kept;
    std::vector<std::string> keptOwner;
};
// Scene e's pieces manifest, its gauge and every piece's mesh. The gauge is the stitch count around
// the longest ring of the entry's gauge piece; the bears keep their calibrated size at every stitch
// width (strand.wale).
bool loadSpiralSource(const SceneEntry& e, const SceneSettings& s, SpiralSource& out, std::string* err);

// The settings the spiral pieces cannot honour: REFUSED with a reason, never ignored.
bool checkSpiralSettings(const SceneEntry& e, const SceneSettings& s, std::string* err);

// The source scene and its parts, none knitted: one spiral part per piece, carrying its mesh with the
// chart as uv, its rings and its bear (spiral_part.h).
bool loadPiecesScene(const SceneEntry& e, const SceneSettings& s, SpiralSource& src, Scene& out, std::string* err);

// What the pieces take from the scene's GLB, read once per process: the base colour (linear) of each
// material, the primitives its selection keeps (library mm of the GLB) and the box of its knitted fabric.
struct PresetMesh {
    std::map<std::string, Vec3> colours;
    std::vector<SourcePrimitive> kept;
    Vec3 fabricLo{1e300, 1e300, 1e300}, fabricHi{-1e300, -1e300, -1e300};
};
bool presetMeshOf(const SceneEntry& e, const PresetMesh*& out, std::string* err);

// The base colours (linear) of the scene's GLB materials that colour its pieces
// (SceneEntry::pieceColourMaterials), by GLB material name.
bool piecesObjectColours(const SceneEntry& e, std::map<std::string, Vec3>& out, std::string* err);

} // namespace knit
