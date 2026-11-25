// The spiral layout's source: the teddy's pieces, placed and gauged, as a scene of parts.
#include "scene/source/sewn_pieces.h"

#include <algorithm>
#include <cstdio>
#include <map>
#include <mutex>

#include "model/cancel.h"
#include "model/choices.h"
#include "model/mesh_check.h"
#include "knitting/knit/pieces.h"
#include "knitting/layout/uv_atlas.h"
#include "io/glb_load.h"
#include "scene/source/scene_parts.h"
#include "scene/source/scene_solids.h"
#include "fibres/fibre.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"
#include "io/files.h"
#include "scene/scene.h"

namespace knit {

namespace {

// An axis-aligned box in scene units.
struct SceneBox {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};  // empty
};

SceneBox boxOf(const TriMesh& m) {
    SceneBox b;
    m.bounds(b.lo, b.hi);
    return b;
}

void growBox(SceneBox& b, const SceneBox& c) {
    b.lo = {std::min(b.lo.x, c.lo.x), std::min(b.lo.y, c.lo.y), std::min(b.lo.z, c.lo.z)};
    b.hi = {std::max(b.hi.x, c.hi.x), std::max(b.hi.y, c.hi.y), std::max(b.hi.z, c.hi.z)};
}

// Length of the overlap of the intervals [l0, h0] and [l1, h1] (0 if disjoint).
double overlapLength(double l0, double h0, double l1, double h1) { return std::max(0.0, std::min(h0, h1) - std::max(l0, l1)); }

// The point a box is placed by: its centre over its floor.
Vec3 floorCentre(const Vec3& lo, const Vec3& hi) { return {0.5 * (lo.x + hi.x), 0.5 * (lo.y + hi.y), lo.z}; }

// The object (bear) whose box overlaps b most, else the one whose box centre is nearest.
std::string ownerObject(const SceneBox& b, const std::map<std::string, SceneBox>& objectBox) {
    double bestOverlap = 0, bestDist = 1e300;
    std::string byOverlap, byDist;
    for (const auto& [obj, ob] : objectBox) {
        const double v = overlapLength(b.lo.x, b.hi.x, ob.lo.x, ob.hi.x) * overlapLength(b.lo.y, b.hi.y, ob.lo.y, ob.hi.y) *
                         overlapLength(b.lo.z, b.hi.z, ob.lo.z, ob.hi.z);
        const double d = norm((b.lo + b.hi) * 0.5 - (ob.lo + ob.hi) * 0.5);
        if (v > bestOverlap) bestOverlap = v, byOverlap = obj;
        if (d < bestDist) bestDist = d, byDist = obj;
    }
    return byOverlap.empty() ? byDist : byOverlap;
}

// The package's origin and the kept parts. The kept parts are the scene's GLB nodes its selection keeps
// (glb.plainNodes: the eyes, noses, mouths, the bow). The GLB is placed on the pieces by its knitted
// fabric's box (centre over the floor), which must have the pieces' size: both are the same posed bears.
// The origin is the bounds of every piece and kept part centred over the floor, as the tiled scenes
// (glb_load.h), so an imported asset pivots at its feet. Each kept part goes with the bear whose box
// overlaps its own most (else the nearest box centre).
bool placePackageOrigin(const SceneEntry& e, SpiralSource& src, std::string* err) {
    SceneBox all;
    std::map<std::string, SceneBox> objectBox;  // piece objects
    for (size_t k = 0; k < src.man.pieces.size(); ++k) {
        TriMesh piece;
        if (!readObj(piecePath(src.man, src.man.pieces[k].file), piece, err)) return false;
        const SceneBox b = boxOf(piece);
        growBox(all, b);
        if (auto [it, added] = objectBox.emplace(src.man.pieces[k].object, b); !added) growBox(it->second, b);
    }
    const SceneBox pieces = all;
    const PresetMesh* glb = nullptr;
    if (!e.path.empty() && !e.glb.plainNodes.empty() && !presetMeshOf(e, glb, err)) return false;
    Vec3 glbFromScene;  // GLB library mm = scene units * sigma + glbFromScene
    std::vector<SceneBox> keptBox;
    if (glb && !glb->kept.empty()) {
        const Vec3 piecesSize = (pieces.hi - pieces.lo) * src.sigma, fabricSize = glb->fabricHi - glb->fabricLo;
        const double worst = std::max({std::abs(fabricSize.x / piecesSize.x - 1), std::abs(fabricSize.y / piecesSize.y - 1),
                                       std::abs(fabricSize.z / piecesSize.z - 1)});
        std::printf("kept parts: %zu nodes of %s placed on the pieces (fabric %.1f x %.1f x %.1f mm, pieces %.1f x %.1f x %.1f mm)\n",
                    glb->kept.size(), e.path.c_str(), fabricSize.x, fabricSize.y, fabricSize.z, piecesSize.x, piecesSize.y, piecesSize.z);
        if (worst > 0.02) {
            if (err) *err = "scene " + e.id + ": the knitted fabric of " + e.path + " is not the pieces' size (" + std::to_string(worst * 100) +
                            " % apart): the kept parts cannot be placed";
            return false;
        }
        glbFromScene = floorCentre(glb->fabricLo, glb->fabricHi) - floorCentre(pieces.lo, pieces.hi) * src.sigma;
        for (const SourcePrimitive& p : glb->kept) {
            SceneBox b;
            p.welded.bounds(b.lo, b.hi);
            b.lo = (b.lo - glbFromScene) * (1 / src.sigma);
            b.hi = (b.hi - glbFromScene) * (1 / src.sigma);
            growBox(all, b);
            keptBox.push_back(b);
        }
    }
    src.origin = floorCentre(all.lo, all.hi);
    std::printf("export: true metres (%.4f m per scene unit), origin (%.4f, %.4f, %.4f) scene units\n", src.sigma / 1000,
                src.origin.x, src.origin.y, src.origin.z);
    if (glb) {
        src.kept = glb->kept;
        for (size_t k = 0; k < src.kept.size(); ++k) {
            src.keptOwner.push_back(ownerObject(keptBox[k], objectBox));
            for (Vec3& q : src.kept[k].welded.pos) q = q - glbFromScene - src.origin * src.sigma;
        }
    }
    return true;
}

std::string fileStem(const std::string& file) { return file.substr(0, file.find_last_of('.')); }

// The node of a manifest object (a bear): numbered in the order its first primitive appears, so the
// export groups every piece and extra of one bear under one object.
int objectNode(std::map<std::string, int>& nodes, const std::string& object) {
    return nodes.emplace(object, static_cast<int>(nodes.size())).first->second;
}

// The source scene of the loaded pieces: one part per piece (in manifest order) and the kept parts as
// plain primitives.
bool piecesSceneOf(const SceneEntry& e, const SceneSettings& s, const SpiralSource& src, Scene& out, std::string* err) {
    out = Scene{};
    // A piece's colour is its GLB fabric's (manifest material -> fabric -> base colour).
    std::map<std::string, Vec3> fabricColour;
    if (!piecesObjectColours(e, fabricColour, err)) return false;
    const std::map<std::string, std::string> fabricOf(e.pieceColourMaterials.begin(), e.pieceColourMaterials.end());
    std::map<std::string, int> nodes;
    SourceScene& source = out.source;
    source.path = src.manifestPath;
    source.fileHash = sceneAssetHash(e);
    for (size_t k = 0; k < src.man.pieces.size(); ++k) {
        cancelCheckpoint();
        const PieceEntry& pe = src.man.pieces[k];
        LoadedPiece lp;
        if (!loadPiece(src.man, static_cast<int>(k), src.sigma, lp, err)) return false;
        ScenePart r;
        r.prim = static_cast<int>(source.prims.size());
        r.component = 0;
        r.name = fileStem(pe.file);
        r.centre = (lp.map.centre - src.origin) * src.sigma;  // library mm about the package origin
        r.mesh = std::make_shared<const TriMesh>(lp.mesh);
        r.layout = KnitLayout::Spiral;
        r.rings = pe.rings;
        SourcePrimitive g;
        g.node = pe.object;
        g.nodeIndex = objectNode(nodes, pe.object);
        r.object = g.nodeIndex;
        g.mesh = r.name;
        const auto fabric = fabricOf.find(pe.material);
        g.material = fabric != fabricOf.end() ? fabric->second : pe.material;
        const auto colour = fabricColour.find(g.material);
        g.baseColor = colour != fabricColour.end() ? colour->second : s.yarnWhite.color;
        g.welded = lp.mesh;
        for (Vec3& q : g.welded.pos) q = q + r.centre;
        g.welded.computeNormals();
        g.check = checkMesh(g.welded);
        if (g.check.components.size() != 1) {
            if (err) *err = r.name + ": " + std::to_string(g.check.components.size()) + " components (a piece is one surface)";
            return false;
        }
        source.prims.push_back(std::move(g));
        out.parts.push_back(std::move(r));
    }
    // The kept parts (the GLB's unselected nodes: eyes, noses, mouths, the bow) as the file has them, each under
    // the bear it belongs to (SpiralSource::keptOwner).
    for (size_t k = 0; k < src.kept.size(); ++k) {
        SourcePrimitive g = src.kept[k];
        g.nodeIndex = objectNode(nodes, src.keptOwner[k]);
        source.prims.push_back(std::move(g));
    }
    source.lo = {1e300, 1e300, 1e300};
    source.hi = {-1e300, -1e300, -1e300};
    for (const SourcePrimitive& g : source.prims)
        for (const Vec3& q : g.welded.pos) {
            source.lo = {std::min(source.lo.x, q.x), std::min(source.lo.y, q.y), std::min(source.lo.z, q.z)};
            source.hi = {std::max(source.hi.x, q.x), std::max(source.hi.y, q.y), std::max(source.hi.z, q.z)};
        }
    // The pieces of a bear are sewn together: they touch along their seams by design, so the
    // clearance between parts is not measured on this layout.
    out.cross.lines = {"not measured: the pieces of a bear are sewn together and touch along their seams"};
    return true;
}

} // namespace

bool loadSpiralSource(const SceneEntry& e, const SceneSettings& s, SpiralSource& out, std::string* err) {
    out = SpiralSource{};
    out.manifestPath = findProjectFile(e.piecesManifest);
    if (out.manifestPath.empty()) {
        if (err) *err = "scene " + e.id + ": pieces manifest " + e.piecesManifest + " not found (the scene preparation tools make it)";
        return false;
    }
    if (!(s.wale > 0) || !(e.gaugeStitches > 0) || !(e.piecesW > 0)) {
        if (err) *err = "scene " + e.id + ": stitch width and the pieces' gauge must be positive";
        return false;
    }
    if (!readManifest(out.manifestPath, out.man, err)) return false;
    std::printf("manifest %s: %zu pieces\n", out.manifestPath.c_str(), out.man.pieces.size());
    // The bears keep their calibrated size at every stitch width: sigma = stitches w / ring is
    // gaugeStitches piecesW / ring when stitches = gaugeStitches piecesW / w (not rounded: a whole
    // count would change the size).
    Gauge g;
    if (!computeGauge(out.man, {e.gaugePiece, e.gaugeStitches * e.piecesW / s.wale, s.wale}, g, err)) return false;
    out.sigma = g.sigma;
    std::printf("gauge: %s longest ring %.4f scene units at a = %.4f -> %.0f stitches of w = %.2f mm: "
                "sigma = %.3f mm per scene unit\n", g.file.c_str(), g.ringLength, g.aAt, g.stitches, g.w, g.sigma);
    return placePackageOrigin(e, out, err);
}

bool checkSpiralSettings(const SceneEntry& e, const SceneSettings& s, std::string* err) {
    auto refuse = [&](const std::string& why) {
        if (err) *err = "scene " + e.id + ": " + why;
        return false;
    };
    if (s.pattern != "stockinette") return refuse("pattern " + s.pattern + ": spiral pieces knit stockinette only (purl INC/DEC not derived)");
    if (s.angleDeg != 0) return refuse("row direction: fixed by each piece's axis on spiral pieces");
    if (s.coverage != "component") return refuse("knit.coverage " + s.coverage + ": spiral pieces have no charts to continue");
    // The look of the knitted pieces on the common export: the pieces have no texture coordinates
    // (no colourwork), and their plies are round.
    if (s.yarnColourSource != "object" && s.yarnColourSource != "single")
        return refuse("yarn colour " + s.yarnColourSource + ": spiral pieces take one colour per bear (object) or one colour (single)");
    if (s.flatten && s.plies > 1) return refuse("contact flattening with explicit plies is not supported (the plies are round)");
    if (s.yarnWhite.model == MaterialSpec::Model::Fibre || s.stuffing.model == MaterialSpec::Model::Fibre)
        return refuse("the fibre (hair) BSDF applies to strands only, not to yarn tubes or the filling");
    int texPlies = 1;
    PlyProfileModel texModel = PlyProfileModel::Swept;
    if (!parseTextureProfile(s.textureProfile, texPlies, texModel, err)) return false;
    if (s.diameterAttempts < 0) return refuse("knit.diameter_attempts: 0 (all) or a positive count");
    if (std::string why; !tubeSettingsValid(s, &why)) return refuse(why);
    // The fibre settings are checked with the knitting, as before the common fibres (a teddy whose
    // fibre settings the generator refuses is refused at its knit).
    FibreParams fp;
    if (!fibreParamsFromSettings(s, fp, err)) return refuse(err ? *err : "fibre settings");
    return true;
}

bool presetMeshOf(const SceneEntry& e, const PresetMesh*& out, std::string* err) {
    // bears.glb is large (the fabric subdivided): read once per process and scene, and only what the pieces take.
    static std::mutex m;
    static std::map<std::string, PresetMesh> loaded;
    std::lock_guard<std::mutex> lk(m);
    auto it = loaded.find(e.id);
    if (it == loaded.end()) {
        const std::string path = findProjectFile(e.path);
        if (path.empty()) {
            if (err) *err = "scene " + e.id + ": " + e.path + " not found (the pieces' colours and kept parts)";
            return false;
        }
        SourceScene g;
        SourceOptions go;
        go.selfIntersection = false;
        if (!loadGlb(path, e.glb, g, err, go)) return false;
        PresetMesh p;
        for (SourcePrimitive& q : g.prims) {
            p.colours.emplace(q.material, q.baseColor);
            if (q.cls == SourcePrimitive::Class::Plain) p.kept.push_back(std::move(q));
            else if (q.cls == SourcePrimitive::Class::Knit)
                for (const Vec3& v : q.welded.pos) {
                    p.fabricLo = {std::min(p.fabricLo.x, v.x), std::min(p.fabricLo.y, v.y), std::min(p.fabricLo.z, v.z)};
                    p.fabricHi = {std::max(p.fabricHi.x, v.x), std::max(p.fabricHi.y, v.y), std::max(p.fabricHi.z, v.z)};
                }
        }
        it = loaded.emplace(e.id, std::move(p)).first;
    }
    out = &it->second;
    return true;
}

bool piecesObjectColours(const SceneEntry& e, std::map<std::string, Vec3>& out, std::string* err) {
    out.clear();
    if (e.pieceColourMaterials.empty()) return true;  // nothing coloured from the GLB
    const PresetMesh* glb = nullptr;
    if (!presetMeshOf(e, glb, err)) return false;
    for (const auto& [id, material] : e.pieceColourMaterials) {
        auto c = glb->colours.find(material);
        if (c == glb->colours.end()) {
            if (err) *err = "scene " + e.id + ": material " + material + " (the colour of " + id + ") is not in " + e.path;
            return false;
        }
        out[material] = c->second;
    }
    return true;
}

bool loadPiecesScene(const SceneEntry& e, const SceneSettings& s, SpiralSource& src, Scene& out, std::string* err) {
    if (!loadSpiralSource(e, s, src, err) || !piecesSceneOf(e, s, src, out, err)) return false;
    // The pieces' closed surfaces: what the fibres hide in (the knit uses none, the pieces touch by design).
    out.solids = std::make_shared<const SceneSolids>(out.source);
    return true;
}

} // namespace knit
