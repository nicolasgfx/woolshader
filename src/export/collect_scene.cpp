// collectScene: a knitted scene and its fibres as an ExportScene.
//
// Stage by stage, in the order the output lists them: the materials, the surfaces that are not
// knitted, each part's tubes and filling, the close-up cameras, the fibres. It computes no knit, no
// fibres and no colour: those are results it is given, and PartMaterials says which material each
// one takes. What it does compute is the form exports want: tube meshes, fillings and metres.
//
// What has no yarn (plain primitives, parts too small for the gauge, failed or skipped parts, the
// body under hair) is written as its source mesh by addPlainMesh, so no part disappears.
#include "export/collect_scene.h"

#include <algorithm>
#include <array>
#include <functional>
#include <cmath>
#include <cstdio>
#include <filesystem>

#include "io/files.h"
#include "model/cancel.h"
#include "model/choices.h"
#include "model/coordinates.h"
#include "model/hash.h"
#include "model/log.h"
#include "scene/appearance/baked_maps.h"
#include "scene/appearance/filling.h"
#include "scene/appearance/ply_maps.h"
#include "scene/appearance/rendered_part.h"
#include "scene/appearance/yarn_tube.h"
#include "scene/fibres/scene_fibres.h"
#include "scene/run_log.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"
#include "settings/yarn_look.h"
#include "fibres/fur_loops.h"

namespace knit {

namespace {

// The export being collected: what is read, what is filled, and what the stages share.
struct Collector {
    const SceneEntry& entry;
    const Scene& input;
    const SceneSettings& settings;
    ExportScene& output;
    PartMaterials materials;
    std::string& log;
    bool stageProfile = false;
    bool fur = false;
    TubeStyle tube;        // ring and sides from the settings, compression, then the map set's uvs
    double mapPartD = 0;   // the part diameter the texture map set is made for (0: no map set)
    std::set<std::pair<int, int>> plainComponents;           // (primitive, component) already written as a plain mesh
    std::map<const ColourImage*, std::string> texturePath;   // a decoded texture -> its PNG
    StrandTaper taper;     // the late taper (fibres.tail_taper*, ExportOptions::taperPoints)
    long long bodyFileSegments = ExportScene::kBodyFileSegments;  // ExportOptions::bodyFileSegments
    std::function<void(const std::string& stage, int done, int total)> onStage;  // ExportOptions::onStage
    BakeInterrupt bake = BakeInterrupt::Cancellable;  // ExportOptions::finishBakes
};

// How a part's tubes are shaded: the ply map set is made for one part diameter.
enum class TubeShading { Plain, PlyMaps };

// An unknitted surface (a plain primitive, a part not knitted, the fur's base mesh).
struct PlainSurface {
    const std::string& name;
    const SourcePrimitive& primitive;
    const std::vector<int>* tris;  // the primitive's triangles in mesh's order (null: all of them)
    const TriMesh& mesh;
};

// One population of a part's fibres, as the strand parts it becomes.
struct FibrePart {
    const PartStrands& strands;
    FibrePopulation population;
    const std::string& name;
};

// The fibres of one colour, with the budget key of each.
struct StrandGroup {
    StrandSetF strands;
    std::vector<float> keys;
};

void addStrandToGroup(StrandGroup& group, const PartStrands& from, size_t strand) {
    const StrandSetF& set = from.strands;
    const uint32_t a = set.first[strand], b = set.first[strand + 1];
    group.strands.add(std::vector<float>(set.xyz.begin() + 3 * a, set.xyz.begin() + 3 * b),
                      std::vector<float>(set.radius.begin() + a, set.radius.begin() + b));
    group.keys.push_back(from.keys[strand]);
}

TriMesh toScene(const TriMesh& lib, const SceneMapping& map) {
    TriMesh t = lib;
    for (Vec3& p : t.pos) p = map.toScene(p);
    t.uv.clear();
    t.computeNormals();
    return t;
}

int partIndex(const Collector& c, const ScenePart& r) { return static_cast<int>(&r - &c.input.parts[0]); }

// The strand's tube (strand.tube_*); ply tubes bring their own. Compressed (flattened) rendered tubes
// are exported at the resolution gate (d) checks, as the pieces' --compress. Knitted scenes compute
// on the GPU, the tubes' surface model and mesh too: a device failure fails the export, never a CPU
// stand-in, because a compressed tube's model at that resolution takes minutes per part on the CPU
// (the teddy's body: 120 s).
TubeStyle tubeStyleOf(const SceneSettings& s) {
    TubeStyle style;
    style.ring = s.tubeRing;
    style.sides = s.tubeSides;
    const SurfaceParams sp = yarnSurfaceParams(s);  // the one owner of the surface the body and the maps follow
    style.compression = sp.compression;
    style.maxCompression = sp.maxCompression;
    style.modelRing = sp.ringSpacing;
    return style;
}

// Texture profile (plies as maps on the single tube, S5): one set for the scene, from the
// content-addressed cache (ply_map_sets.h), made at the rendered diameter of the scene's most common
// part diameter. A set is never reused across D / pitch ratios: parts knitted at another D (fallback
// chain) stay plain. Explicit plies replace the texture profile (it applies to a single tube only).
bool attachTextureProfile(Collector& c, std::string* err) {
    const SceneSettings& s = c.settings;
    if (s.yarnBody.enabled) {  // the yarn body's fibres are the relief the maps would draw
        appendFormat(c.log, "  texture profile %s: not attached (the yarn body draws the fibres)\n", s.textureProfile.c_str());
        return true;
    }
    int texPlies = 1;
    PlyProfileModel texModel = PlyProfileModel::Swept;
    if (!parseTextureProfile(s.textureProfile, texPlies, texModel, err)) return false;
    if (!(s.plyPhaseJitter >= 0 && s.plyPhaseJitter <= 1) || !(s.plyPitchJitter >= 0 && s.plyPitchJitter <= 0.5) ||
        !(s.plyCrownTiltDeg > 0 && s.plyCrownTiltDeg < 90)) {
        if (err) *err = "strand.ply_phase_jitter 0..1, strand.ply_pitch_jitter 0..0.5, strand.ply_crown_tilt_deg between 0 and 90";
        return false;
    }
    PlyModeOptions plyOpt;
    PlyGeometry plyGeom;
    plyGeom.count = texPlies;
    plyGeom.pitch = s.plyPitch;
    plyGeom.fill = s.plyFill;
    plyGeom.handedness = s.plyHandedness;
    plyOpt.mode = texPlies > 1 && s.plies <= 1 ? PlyModeOptions::Mode::Texture : PlyModeOptions::Mode::None;
    plyOpt.spec.pp.model = texModel;
    plyOpt.spec.pp.blendTiltDeg = s.plyCrownTiltDeg;
    plyOpt.spec.g = plyGeom;
    if (plyOpt.textured()) {
        std::map<double, int> nD;
        for (const ScenePart& r : c.input.parts)
            if (r.hasStitches()) ++nD[r.yarnD];
        for (const auto& [D, n] : nD)
            if (c.mapPartD == 0 || n > nD[c.mapPartD]) c.mapPartD = D;
    }
    if (c.mapPartD <= 0) return true;
    plyOpt.spec.yarnD = finalSurface(s, c.mapPartD).renderedDiameter();
    PlyMapSetResult ms;
    if (!ensurePlyMapSet(plyMapCacheRoot(), plyOpt.spec, ms, err)) return false;
    const std::optional<PlyShading> shading = attachPlySet(ms, plyOpt, c.output.folder(), 1000.0, err);
    if (!shading) return false;
    c.output.plyShading = *shading;
    c.tube.uvRepeat = shading->uvRepeat;
    c.tube.uvWrap = shading->uvWrap;
    c.tube.uvPhaseJitter = s.plyPhaseJitter;
    c.tube.uvPitchJitter = s.plyPitchJitter;
    c.output.plyShading.normalScale = s.normalStrength;
    appendFormat(c.log, "  texture profile %s: map set %s (%s, %.0f ms) at D_r %.3f mm (parts at D %.2f)%s\n", s.textureProfile.c_str(),
            ms.dir.c_str(), ms.hit ? "cache hit" : "generated", ms.ms, plyOpt.spec.yarnD, c.mapPartD,
            renderFlattened(s) ? "; relief approximated on the elliptical section" : "");
    return true;
}

bool addMaterials(Collector& c, std::string* err) {
    c.materials.addSceneMaterials();
    return attachTextureProfile(c, err);
}

// The uvs a primitive's colour texture lies on, one set per welded triangle.
const std::vector<std::array<Vec2, 3>>& colourUvsOf(const SourcePrimitive& g) {
    return !g.colourUv.empty() ? g.colourUv : g.welded.uv;
}

// An unknitted surface in the source's base colour texture on its colour uvs: one PNG per decoded
// texture, one material per material and texture (the base colour factor times the texture).
bool addTexturedPlainMesh(Collector& c, const PlainSurface& surface, std::string* err) {
    const SourcePrimitive& g = surface.primitive;
    const std::vector<std::array<Vec2, 3>>& uvs = colourUvsOf(g);
    auto known = c.texturePath.find(g.albedo.get());
    if (known == c.texturePath.end()) {
        const std::string name = std::to_string(c.texturePath.size());
        const std::optional<std::string> path = c.output.addColourTexture(name, *g.albedo, err);
        if (!path) return false;
        known = c.texturePath.emplace(g.albedo.get(), *path).first;
    }
    const std::string material = c.materials.texturedPlainMaterial(g, known->second);
    std::vector<std::array<Vec2, 3>> cornerUv(surface.mesh.tri.size());
    for (size_t q = 0; q < cornerUv.size(); ++q) cornerUv[q] = uvs[surface.tris ? (*surface.tris)[q] : static_cast<int>(q)];
    return c.output.addTexturedMesh(surface.name, material, surface.mesh, cornerUv, err);
}

// An unknitted surface, coloured by material.plain_colour_source: the texture, the part's
// area-weighted mean base colour, or the material's base colour (glTF: the factor; imports: the mean colour).
bool addPlainMesh(Collector& c, const PlainSurface& surface, std::string* err) {
    const SourcePrimitive& g = surface.primitive;
    const size_t nt = surface.tris ? surface.tris->size() : g.welded.tri.size();
    const std::string& colourSource = c.settings.plainColourSource;
    if (colourSource == "texture" && g.albedo && colourUvsOf(g).size() == g.welded.tri.size() && surface.mesh.tri.size() == nt)
        return addTexturedPlainMesh(c, surface, err);
    Vec3 mean;
    if (colourSource != "material" && meanTriangleColour(g, surface.tris, mean))
        return c.output.addMesh(surface.name, c.materials.meanColourMaterial(surface.name, mean), surface.mesh, err);
    return c.output.addMesh(surface.name, c.materials.plain[g.material], surface.mesh, err);
}

// Plain primitives as they are (their own colour).
bool addPlainPrimitives(Collector& c, std::string* err) {
    const SourceScene& source = c.input.source;
    const std::vector<std::string> names = plainPrimitiveNames(source);
    for (size_t pi = 0; pi < source.prims.size(); ++pi) {
        const SourcePrimitive& g = source.prims[pi];
        const std::string& name = names[pi];
        if (name.empty()) continue;
        cancelCheckpoint();
        if (!addPlainMesh(c, {name, g, nullptr, toScene(g.welded, sceneMetres({}))}, err)) return false;
        c.output.place(name, exportObjectOf(source, static_cast<int>(pi)));
    }
    return true;
}

// A knitted part's tubes: one mesh, or with texture colourwork one per stitch colour.
bool addPartTubes(Collector& c, const ScenePart& r, const RenderedPart& rp, TubeShading shading, std::string* err) {
    const int part = partIndex(c, r);
    const PartMaterials& m = c.materials;
    const auto materialOf = [&](const std::string& yarn) {
        return shading == TubeShading::PlyMaps ? ExportScene::yarnMaterial(yarn) : yarn;
    };
    TubeStyle style = c.tube;
    if (rp.plies) {
        style.ring = rp.ring;
        style.sides = rp.sides;
    }
    const std::optional<ExportTube> tube = buildExportTube(*rp.strands, rp.strandPp, sceneMetres(r.centre), style, err);
    if (!tube) return false;
    const auto pc = m.stitchColours.find(part);
    if (pc == m.stitchColours.end() || rp.yarn->topo.empty())
        return c.output.addMeshPart(r.name, materialOf(m.yarnOf(part)), tube->mesh, err);
    // Colourwork: one tube mesh per stitch colour, each triangle with its stitch's colour.
    const std::vector<uint32_t>& key = pc->second.key;
    const StitchAlongYarn along(*rp.yarn, rp.plies.get());
    TubeGroups groups;
    for (uint32_t colour : m.colourKeys) {
        const std::string id = stitchColourId(colour);
        groups.groups.push_back({"__" + id.substr(5), materialOf(id)});
    }
    const std::vector<uint32_t>& tri = tube->mesh.tri;
    for (size_t t = 0; t + 2 < tri.size(); t += 3) {
        const uint32_t a = tri[t], b = tri[t + 1], d = tri[t + 2];
        const int st = along.at(tube->yarn[a], (tube->s[a] + tube->s[b] + tube->s[d]) / 3);
        const uint32_t colour = st >= 0 && st < static_cast<int>(key.size()) ? key[st] : m.dominant.at(part);
        groups.groupOfTriangle.push_back(m.colourIndex.at(colour));
    }
    const std::optional<std::vector<std::string>> parts = c.output.addTubeGroups(r.name, *tube, groups, err);
    if (!parts) return false;
    for (const std::string& n : *parts) c.output.place(n, exportObjectOf(c.input.source, r.prim));
    appendFormat(c.log, "  %s: %zu stitch colours (one tube mesh each)\n", r.name.c_str(), parts->size());
    return true;
}

// The filling behind a knitted part (filling.h): its source surface moved behind the yarn, in the
// material "stuffing".
bool addPartFilling(Collector& c, const ScenePart& r, std::string* err) {
    const SourcePrimitive& g = c.input.source.prims[r.prim];
    const Filling f = makeFilling(*r.mesh, g.check.components[r.component], fillingOffset(c.settings, r.knit->t.pp.yarnDiameter));
    appendFormat(c.log, "  %s: %s\n", r.name.c_str(), f.summary().c_str());
    if (f.kind == Filling::Kind::None) return true;
    const std::string name = r.name + (f.kind == Filling::Kind::Stuffing ? "_stuffing" : "_backing");
    if (!c.output.addMesh(name, "stuffing", toScene(f.mesh, sceneMetres(r.centre)), err)) return false;
    c.output.place(name, exportObjectOf(c.input.source, r.prim));
    return true;
}

// The yarn of fibres rendered as maps (yarn_body.render = maps; fibres/yarn_bake.h, scene/appearance/
// baked_maps.h): one 16-sided tube per yarn under the maps baked from its fibres, in the yarn body's group, no
// body strands and no cores. The maps come from the cache, baked once per knit, body and texel; the material's
// values (the dye, cavity, roughness floor, highlight) are the settings', applied as the maps are attached, so
// a material edit never bakes. Progress: a yarn at a time while baking.
bool addPartMappedTubes(Collector& c, const ScenePart& r, const RenderedPart& rp, const YarnBodyParams& p, std::string* err) {
    const SceneSettings::YarnBody& b = c.settings.yarnBody;
    const std::string yarn = c.materials.yarnOf(partIndex(c, r));
    const auto base = std::find_if(c.output.baseMaterials.begin(), c.output.baseMaterials.end(), [&](const BaseMaterial& m) { return m.id == yarn; });
    if (base == c.output.baseMaterials.end()) return fail(err, r.name + ": no base material " + yarn + " for the yarn maps");
    const MaterialSpec spec = base->spec;  // a copy: addBakedYarn grows baseMaterials
    SurfaceModel surface;
    if (!buildYarnSurface(*rp.strands, rp.strandPp, c.tube, surface, err)) return false;
    const BakeParams bp = bakeParams(c.settings);
    BakedMapSpec ms;
    if (!bakedMapSpecOf(c.settings, r.name, *rp.strands, rp.strandPp, ms, err)) return false;
    BakeProgress progress;
    if (c.onStage)
        progress.report = [&](int done, int yarns) {
            c.onStage("baking " + r.name + ": yarn " + std::to_string(done) + " of " + std::to_string(yarns), done, yarns);
        };
    BakedMapEntry entry;
    if (!ensureBakedMapSet(bakeCacheRoot(), ms, *rp.strands, surface, c.bake, progress, entry, err)) return false;
    if (entry.yarns.size() != surface.yarns.size())
        return fail(err, r.name + ": the cached maps hold " + std::to_string(entry.yarns.size()) + " yarns, the knit " +
                             std::to_string(surface.yarns.size()) + " (" + entry.dir + ")");
    const std::optional<ExportTube> tube = buildBakedExportTube(*rp.strands, rp.strandPp, sceneMetres(r.centre), c.tube, bp, err);
    if (!tube) return false;
    std::map<int, int> indexOfYarnId;
    TubeGroups groups;
    for (size_t y = 0; y < surface.yarns.size(); ++y) {
        cancelCheckpoint();
        const YarnSurface& ys = surface.yarns[y];
        indexOfYarnId[ys.yarnId] = static_cast<int>(y);
        if (!bakedPlanMatches(entry.yarns[y], ys.length, ys.radius, ms))
            return fail(err, r.name + ": the cached maps of yarn " + std::to_string(y) + " do not fit its tube (" + entry.dir + ")");
        char suffix[16];
        std::snprintf(suffix, sizeof suffix, "_y%03zu", y);
        const std::string id = "baked_" + yarn + "_" + materialIdOf(r.name) + suffix;
        BakedYarnFiles files;
        if (!attachBakedYarn(entry, static_cast<int>(y), c.output.folder(), "textures/" + id, b.mapCavity, b.mapRoughness, files, err)) return false;
        BakedMaps shading;
        shading.fabric = b.mapMaterial != "shiny";
        shading.cavity = files.cavity;
        shading.normal = files.normal;
        shading.specular = files.specular;
        shading.roughness = files.roughness;
        shading.height = files.height;
        shading.specularScale = b.mapSpecular;
        shading.roughnessAlong = b.mapAnisotropy;
        shading.anisoRotation = -entry.yarns[y].fibreAngleDeg / 180;  // half-turns; the renderer's v axis is -theta (BakedMaps)
        c.output.addBakedYarn(id, spec, shading);
        groups.groups.push_back({suffix, id});
    }
    const std::vector<uint32_t>& tri = tube->mesh.tri;
    for (size_t k = 0; k + 2 < tri.size(); k += 3) groups.groupOfTriangle.push_back(indexOfYarnId.at(tube->yarn[tri[k]]));
    const std::optional<std::vector<std::string>> parts = c.output.addTubeGroups(r.name + "_baked", *tube, groups, err);
    if (!parts) return false;
    for (const std::string& n : *parts) c.output.place(n, exportObjectOf(c.input.source, r.prim, ExportScene::kBodyGroup));
    appendFormat(c.log, "  %s yarn maps (%d plies): %zu yarns, texel %.0f um, %s in %.0f s (%s)\n", r.name.c_str(), p.plies, parts->size(),
                 1e3 * bp.texel, entry.hit ? "from the cache" : "baked", entry.ms / 1000, entry.dir.c_str());
    return true;
}

// The yarn body (fibres/yarn_body.h, a prototype): the part's yarn as structural fibres on the surface
// model its tube is meshed from, in place of the tube, and the core proxy when asked. Colourwork gives
// the body its part's yarn colour only (logged).
bool addPartBody(Collector& c, const ScenePart& r, const RenderedPart& rp, std::string* err) {
    YarnBodyParams p;
    if (!yarnBodyParams(c.settings, p, err)) return false;
    SurfaceModel surface;
    if (!buildYarnSurface(*rp.strands, rp.strandPp, c.tube, surface, err)) return false;
    YarnBody body;
    if (!buildYarnBody(*rp.strands, surface, p, fnv1a64(r.name), body, err, ExecPolicy{true})) return false;
    const int part = partIndex(c, r);
    const std::string yarn = c.materials.yarnOf(part);
    if (c.materials.stitchColours.count(part))
        appendFormat(c.log, "  %s: the yarn body takes the part's yarn colour (%s), not its stitch colours\n", r.name.c_str(), yarn.c_str());
    // The strand files: whole strands, at most c.bodyFileSegments segments each, all in the body's group;
    // one file keeps the name <part>_body, several are <part>_body_00, _01, ...
    const SceneMapping map = sceneMetres(r.centre);
    const std::string material = c.output.bodyFibreMaterial(yarn);
    const Placement where = exportObjectOf(c.input.source, r.prim, ExportScene::kBodyGroup);
    const std::vector<size_t> cuts = ExportScene::strandFileCuts(body.strands.first, c.bodyFileSegments);
    std::vector<float> pts, radii;
    for (size_t k = 0; k + 1 < cuts.size(); ++k) {
        cancelCheckpoint();
        StrandSetF strands;
        for (size_t i = cuts[k]; i < cuts[k + 1]; ++i) {
            pts.clear();
            radii.clear();
            for (uint64_t v = body.strands.first[i]; v < body.strands.first[i + 1]; ++v) {
                const Vec3 q = map.toExport(body.strands.points[v]);
                pts.insert(pts.end(), {static_cast<float>(q.x), static_cast<float>(q.y), static_cast<float>(q.z)});
                radii.push_back(static_cast<float>(body.strands.radius[v] / map.sigma));
            }
            strands.add(pts, radii);
        }
        char suffix[16];
        std::snprintf(suffix, sizeof suffix, "_body_%02zu", k);
        const std::string name = cuts.size() == 2 ? r.name + "_body" : r.name + suffix;
        if (!c.output.addStrands(name, material, strands, err)) return false;
        c.output.place(name, where);
    }
    appendFormat(c.log, "  %s yarn body (%d plies, %s core): %s, %zu strand file(s)\n", r.name.c_str(), p.plies,
                 c.settings.yarnBody.core.c_str(), body.stats.line().c_str(), cuts.size() - 1);
    if (p.core == BodyCore::None) return true;
    const std::string core = r.name + "_core";
    if (!c.output.addMesh(core, c.output.bodyCoreMaterial(yarn, p.core == BodyCore::Translucent), toScene(body.core, map), err)) return false;
    c.output.place(core, exportObjectOf(c.input.source, r.prim, ExportScene::kBodyGroup));
    return true;
}

bool addKnittedPart(Collector& c, const ScenePart& r, std::string* err) {
    const PatternParams& pp = r.knit->t.pp;
    StageClock clock;
    const std::optional<RenderedPart> rendered = renderedPart(r, c.settings, err);
    if (!rendered) return false;
    const RenderedPart& rp = *rendered;
    if (c.stageProfile) appendFormat(c.log, "  %s profile: %s\n", r.name.c_str(), clock.line("renderedPart").c_str());
    clock = StageClock{};
    const bool mapped = c.output.plyShading.on && r.yarnD == c.mapPartD;
    if (c.output.plyShading.on && !mapped)
        appendFormat(c.log, "  %s: plain tube (knitted at D %.2f, the map set is for D %.2f)\n", r.name.c_str(), r.yarnD, c.mapPartD);
    if (!c.settings.yarnBody.enabled) {
        if (!addPartTubes(c, r, rp, mapped ? TubeShading::PlyMaps : TubeShading::Plain, err)) return false;
    } else {
        // A yarn made of fibres: its strands, or its tubes under the maps baked from them.
        YarnBodyParams body;
        if (!yarnBodyParams(c.settings, body, err)) return false;
        switch (yarnBodyRenderOf(c.settings)) {
        case YarnBodyRender::Fibres:
            if (!addPartBody(c, r, rp, err)) return false;
            break;
        case YarnBodyRender::Maps:
            if (!addPartMappedTubes(c, r, rp, body, err)) return false;
            break;
        }
    }
    if (c.stageProfile) appendFormat(c.log, "  %s profile: %s\n", r.name.c_str(), clock.line(c.settings.yarnBody.enabled ? "addBody" : "addTubes").c_str());
    if (rp.plies)
        appendFormat(c.log, "  %s: %d plies, r_p %.3f mm, containment e %.2e mm (%s c/2), %d folds FAIL, %d UNRESOLVED\n", r.name.c_str(),
                c.settings.plies, rp.plies->layout.rp, rp.plies->containment, rp.plies->containment < 0.5 * pp.clearance ? "<" : ">=",
                rp.plyFolds, rp.plyFoldsUnresolved);
    if (!c.settings.filling) return true;
    clock = StageClock{};
    if (!addPartFilling(c, r, err)) return false;
    if (c.stageProfile) appendFormat(c.log, "  %s profile: %s\n", r.name.c_str(), clock.line("filling").c_str());
    return true;
}

// One part. Knitted: its tubes and filling. Not knitted (too small for the gauge, or no stitch fits):
// its surface stays visible as a plain mesh in its own colour.
bool addPart(Collector& c, const ScenePart& r, std::string* err) {
    // The chart parts of a component share its mesh: an unknitted component is written once.
    if (!r.hasStitches() && !c.plainComponents.insert({r.prim, r.component}).second) return true;
    const SourcePrimitive& g = c.input.source.prims[r.prim];
    c.output.place(r.name, exportObjectOf(c.input.source, r.prim));
    if (r.hasStitches()) return addKnittedPart(c, r, err);
    const std::vector<int>& tris = g.check.components[r.component].triangles;
    if (!addPlainMesh(c, {r.name, g, &tris, toScene(*r.mesh, sceneMetres(r.centre))}, err)) return false;
    appendFormat(c.log, "  %s: plain mesh (not knitted: %s)\n", r.name.c_str(), r.note.c_str());
    return true;
}

// The close-up cameras of a part (closeupPartOf): the mesh vertex within 15 % of the part's height
// around its middle whose normal faces the front camera most (scene front = -Y, Blender's front view;
// export +Z), seen along that normal from 40 mm ("closeup_<part>"), 115 mm ("mid"), 180 mm
// ("medium_<part>"), all with a 20 degree field of view, and 400 mm with 30 degrees ("far"). The
// renderer selects one by name (--camera mid), the Blender export makes each a camera object.
void addCloseupCameras(Collector& c, int part) {
    const ScenePart& r = c.input.parts[static_cast<size_t>(part)];
    TriMesh mesh = *r.mesh;
    mesh.computeNormals();
    Vec3 lo, hi;
    mesh.bounds(lo, hi);
    double best = -2;
    Vec3 at, n;
    // band: the middle 15 % of the height; a coarse mesh with no vertex there takes every vertex.
    for (double band : {0.15, 1.0}) {
        for (size_t v = 0; v < mesh.pos.size(); ++v) {
            const Vec3& p = mesh.pos[v];
            if (std::fabs(p.z - 0.5 * (lo.z + hi.z)) > band * (hi.z - lo.z)) continue;
            const Vec3 nv = normalized(mesh.vnormal[v]);
            if (-nv.y > best) {
                best = -nv.y;
                at = p;
                n = nv;
            }
        }
        if (best > -2) break;
    }
    const SceneMapping map = sceneMetres(r.centre);
    c.output.extraCameras.push_back({"closeup_" + r.name, map.toExport(at + n * 40.0), map.toExport(at), 20.0});
    c.output.extraCameras.push_back({"medium_" + r.name, map.toExport(at + n * 180.0), map.toExport(at), 20.0});
    c.output.extraCameras.push_back({"mid", map.toExport(at + n * 115.0), map.toExport(at), 20.0});
    c.output.extraCameras.push_back({"far", map.toExport(at + n * 400.0), map.toExport(at), 30.0});
    appendFormat(c.log, "  cameras closeup_%s, mid, medium_%s and far (part %d)\n", r.name.c_str(), r.name.c_str(), part);
}

// Colourwork: the fibres of each stitch colour as their own part (the fibre's anchor on its strand:
// the stitch it grows from). Fibres on separate ply tubes find it through the ply set.
bool addFibresByStitchColour(Collector& c, const PartFibres& pf, const FibrePart& fibres, std::string* err) {
    const ScenePart& part = c.input.parts[static_cast<size_t>(pf.part)];
    const StitchColours& colours = c.materials.stitchColours.at(pf.part);
    const bool random = fibres.population == FibrePopulation::Random;
    std::optional<RenderedPart> plied;
    if (c.settings.plies > 1) {
        plied = renderedPart(part, c.settings, err);
        if (!plied) return false;
    }
    const StitchAlongYarn along(part.knit->t.curves, plied ? plied->plies.get() : nullptr);
    std::map<int, StrandGroup> byColour;
    size_t j = 0;
    for (const Fibre& f : pf.fs.fibres) {
        if (f.population != fibres.population) continue;
        const int st = along.at(f.yarnId, f.s);
        const uint32_t key = st >= 0 && st < static_cast<int>(colours.key.size()) ? colours.key[st] : c.materials.dominant.at(pf.part);
        addStrandToGroup(byColour[c.materials.colourIndex.at(key)], fibres.strands, j);
        ++j;
    }
    for (const auto& [g, group] : byColour) {
        const std::string id = stitchColourId(c.materials.colourKeys[g]), name = fibres.name + "__" + id.substr(5);
        const std::string material = random ? c.output.randomFibreMaterial(id) : c.output.orientedFibreMaterial(id);
        if (!c.output.addStrands(name, material, keepStrands(group.strands, group.keys, pf.keyBelow), err)) return false;
        c.output.place(name, exportObjectOf(c.input.source, part.prim, random ? ExportScene::kRandomGroup : ExportScene::kOrientedGroup));
    }
    return true;
}

// Fur in texture colours: the hairs of each palette colour as their own strand part, in a material of
// that colour shared by the parts.
bool addFurByColour(Collector& c, const PartFibres& pf, const PartStrands& strands, std::string* err) {
    const ScenePart& part = c.input.parts[static_cast<size_t>(pf.part)];
    std::map<uint32_t, StrandGroup> byColour;
    for (size_t j = 0; j < pf.fibreColour.size(); ++j) addStrandToGroup(byColour[pf.fibreColour[j]], strands, j);
    for (const auto& [key, group] : byColour) {
        const std::string id = c.materials.furColourMaterial(key), name = pf.name + "_fibres__" + id.substr(4);
        if (!c.output.addStrands(name, id, keepStrands(group.strands, group.keys, pf.keyBelow), err)) return false;
        c.output.place(name, exportObjectOf(c.input.source, part.prim, "fur"));
    }
    return true;
}

// One part's fibres below the budget's key: fur as one strand part in its own material, or one per
// texture colour; wool as one part per population, each with its own material per yarn colour.
bool addPartFibres(Collector& c, const PartFibres& pf, long long& exported, std::string* err) {
    const ScenePart& part = c.input.parts[static_cast<size_t>(pf.part)];
    const SourceScene& source = c.input.source;
    const StrandTaper& taper = c.taper;
    const PartStrands orientedStrands = partStrands(pf, part.centre, FibrePopulation::Directed, taper);
    const PartStrands randomStrands = partStrands(pf, part.centre, FibrePopulation::Random, taper);
    if (taper.tail > 0 && taper.points == TaperPoints::Dense)
        appendFormat(c.log, "  %s late taper (tail %.2f, at least %.1f diameters, tip %.2f): oriented %lld tapered, %lld kept their taper; "
                     "random %lld tapered, %lld kept their taper\n", pf.name.c_str(), taper.tail, taper.minDiameters, taper.tip,
                     orientedStrands.tapered, orientedStrands.untapered, randomStrands.tapered, randomStrands.untapered);
    else if (taper.tail > 0)
        appendFormat(c.log, "  %s late taper (tail %.2f, at least %.1f diameters, tip %.2f, points %s, tolerance %.3f): oriented %lld "
                     "tapered, %lld kept their taper; random %lld tapered, %lld kept their taper\n", pf.name.c_str(), taper.tail,
                     taper.minDiameters, taper.tip, taperPointsName(taper.points), taper.tolerance, orientedStrands.tapered,
                     orientedStrands.untapered, randomStrands.tapered, randomStrands.untapered);
    const StrandSetF d = keepStrands(orientedStrands.strands, orientedStrands.keys, pf.keyBelow);
    const StrandSetF rr = keepStrands(randomStrands.strands, randomStrands.keys, pf.keyBelow);
    exported += static_cast<long long>(d.vertices() - d.strands() + rr.vertices() - rr.strands());
    if (c.fur) {
        if (c.settings.pileColourSource == "texture" && !pf.fibreColour.empty()) return addFurByColour(c, pf, orientedStrands, err);
        const auto own = c.materials.fur.find(pf.part);
        const std::string name = pf.name + "_fibres";
        if (!c.output.addStrands(name, own != c.materials.fur.end() ? own->second : std::string("fur"), d, err)) return false;
        c.output.place(name, exportObjectOf(source, part.prim, "fur"));
        return true;
    }
    const std::string oriented = pf.name + "_fibres_oriented", random = pf.name + "_fibres_random";
    const std::string yarnId = c.materials.yarnOf(pf.part);
    const bool byStitch = c.materials.stitchColours.count(pf.part) && part.knit;
    if (byStitch) {
        if (!addFibresByStitchColour(c, pf, {orientedStrands, FibrePopulation::Directed, oriented}, err)) return false;
    } else {
        if (!c.output.addStrands(oriented, c.output.orientedFibreMaterial(yarnId), d, err)) return false;
        c.output.place(oriented, exportObjectOf(source, part.prim, ExportScene::kOrientedGroup));
    }
    if (byStitch && c.settings.fibreRandomInheritsYarn)
        return addFibresByStitchColour(c, pf, {randomStrands, FibrePopulation::Random, random}, err);
    if (!c.output.addStrands(random, c.output.randomFibreMaterial(yarnId), rr, err)) return false;
    c.output.place(random, exportObjectOf(source, part.prim, ExportScene::kRandomGroup));
    return true;
}

// The given fibre result as it is, each part below its object's segment budget threshold.
bool addFibres(Collector& c, const SceneFibres& fibres, std::string* err) {
    appendFormat(c.log, "fibres: the given fibre result (%zu parts), not generated again\n", fibres.parts.size());
    c.log += fibres.log;
    long long before = 0, after = 0;
    for (const PartFibres& pf : fibres.parts) before += segmentCount(pf.fs);
    for (const PartFibres& pf : fibres.parts) {
        cancelCheckpoint();
        if (!addPartFibres(c, pf, after, err)) return false;
    }
    appendFormat(c.log, "strand budget %lld per object: %lld segments generated -> %lld exported\n", c.settings.segmentBudget, before, after);
    return true;
}

// Validate what was written before anything is published: every strand file reads back.
bool readStrandsBack(const ExportScene& scene, std::string* err) {
    for (const PartFile& part : scene.strandParts()) {
        cancelCheckpoint();
        StrandSetF back;
        if (!readStrandsFile(scene.pathOf(part), back, err)) return false;
    }
    return true;
}

} // namespace

namespace {

// A tube mesh in Y-up metres, wound for the renderer.
ExportTube exportTubeOf(TubeMesh& tm, const SceneMapping& map) {
    ExportTube out;
    SceneMesh& mesh = out.mesh;
    mesh.pos.resize(tm.pos.size());
    mesh.nrm.resize(tm.nrm.size());
    for (size_t i = 0; i < tm.pos.size(); ++i) {
        mesh.pos[i] = map.toExport(tm.pos[i]);
        mesh.nrm[i] = yUp(tm.nrm[i]);
    }
    // Winding: counter-clockwise seen along the (outward) vertex normals, as the renderer expects
    // (buildTubeMesh's strips wind the other way: its import warned "winding opposes normals").
    mesh.tri = tm.tri;
    for (size_t t = 0; t + 2 < mesh.tri.size(); t += 3) {
        const Vec3 &a = mesh.pos[mesh.tri[t]], &b = mesh.pos[mesh.tri[t + 1]], &c = mesh.pos[mesh.tri[t + 2]];
        if (dot(cross(b - a, c - a), mesh.nrm[mesh.tri[t]] + mesh.nrm[mesh.tri[t + 1]] + mesh.nrm[mesh.tri[t + 2]]) < 0)
            std::swap(mesh.tri[t + 1], mesh.tri[t + 2]);
    }
    mesh.uv = std::move(tm.uv);
    out.yarn = std::move(tm.yarn);
    out.s = std::move(tm.s);
    return out;
}

} // namespace

std::optional<ExportTube> buildExportTube(const CurveSet& cs, const PatternParams& pp, const SceneMapping& map,
                                          const TubeStyle& style, std::string* err) {
    TubeMesh tm;
    if (!buildYarnTube(cs, pp, style, tm, err)) return std::nullopt;
    return exportTubeOf(tm, map);
}

std::optional<ExportTube> buildBakedExportTube(const CurveSet& cs, const PatternParams& pp, const SceneMapping& map,
                                               const TubeStyle& style, const BakeParams& bake, std::string* err) {
    TubeStyle baked = style;
    baked.ring = kBakeRing;
    baked.sides = kBakeSides;
    baked.uvRepeat = bake.texel * bake.atlasWidth;  // one strip of yarn per unit of u
    baked.uvWrap = 1;                               // the ring doubled where a strip ends
    baked.uvPhaseJitter = baked.uvPitchJitter = 0;
    TubeMesh tm;
    if (!buildYarnTube(cs, pp, baked, tm, err)) return std::nullopt;
    relabelAtlasUv(tm, bake);
    return exportTubeOf(tm, map);
}

int closeupPartOf(const std::vector<ScenePart>& parts, int requested) {
    if (requested >= 0) return requested;
    int best = -1;
    double bestArea = 0;
    for (size_t i = 0; i < parts.size(); ++i) {
        if (!parts[i].hasStitches() || !parts[i].mesh) continue;
        const double area = parts[i].mesh->area();
        if (best < 0 || area > bestArea) best = static_cast<int>(i), bestArea = area;
    }
    return best;
}

Placement exportObjectOf(const SourceScene& source, int prim, const std::string& group) {
    const SourcePrimitive& g = source.prims[prim];
    const int object = sceneObjectOf(source, prim);
    char key[32];
    std::snprintf(key, sizeof key, object >= 0 ? "node%06d" : "prim%06d", object >= 0 ? object : prim);
    const std::string name = !g.node.empty()   ? g.node
                             : !g.mesh.empty() ? g.mesh
                             : g.nodeIndex >= 0 ? "node_" + std::to_string(g.nodeIndex)
                                                : "primitive_" + std::to_string(prim);
    return {key, name, group};
}

bool collectScene(const ExportSource& source, const ExportOptions& options, ExportScene& scene, std::string& log,
                  std::string* err) {
    const SceneEntry& e = source.entry;
    const Scene& t = source.scene;
    const SceneSettings& s = source.settings;
    const bool fur = furWorkflow(s);
    Collector c{e, t, s, scene, PartMaterials(t, s, scene, log), log, options.stageProfile, fur, tubeStyleOf(s)};
    scene.projectName = e.projectName;
    scene.fibreTipScale = s.tailTaper > 0 ? 1.0 : scene.fibreTipScale;
    c.taper = StrandTaper{s.tailTaper, s.tailTaperMinDiameters, s.tailTaperTip, TaperPoints::Dense, options.taperTolerance};
    if (!options.taperPoints.empty() && !parseTaperPoints(options.taperPoints, c.taper.points)) {
        if (err) *err = "unknown taper points '" + options.taperPoints + "' (dense, existing, onset, refined)";
        return false;
    }
    if (c.taper.points == TaperPoints::Refined && !(options.taperTolerance > 0 && options.taperTolerance < 1)) {
        if (err) *err = "the taper tolerance must lie above 0 and below 1 (a fraction of the root radius)";
        return false;
    }
    StrandRules rules = taperRules(c.taper);
    if (furWorkflow(s)) {
        // The pile's thick loops keep the renderer's own joint rule (fur_loops.h kDiscMinSegOverR), never wool.
        rules.minSegOverR = std::min(rules.minSegOverR, kDiscMinSegOverR);
        if (rules.minSegOverRMax > 0) rules.minSegOverRMax = std::min(rules.minSegOverRMax, kDiscMinSegOverR);
    }
    scene.setStrandRules(rules);
    if (e.kind == SceneEntry::Kind::Mesh && t.source.imported)
        scene.importFrame = ImportFrame{t.source.scale, t.source.offset, t.source.importYUp ? UpAxis::Y : UpAxis::Z};
    StageClock clock;
    const auto profileLine = [&](const char* label) {
        if (options.stageProfile) appendFormat(log, "profile: %s\n", clock.line(label).c_str());
        clock = StageClock{};
    };
    c.bodyFileSegments = options.bodyFileSegments;
    c.onStage = options.onStage;
    c.bake = options.finishBakes ? BakeInterrupt::FinishOnce : BakeInterrupt::Cancellable;
    if (YarnBodyParams body; s.yarnBody.enabled && !yarnBodyParams(s, body, err)) return false;  // refused before anything is written
    if (!addMaterials(c, err)) return false;
    profileLine("export: materials and stitch colours");
    if (!addPlainPrimitives(c, err)) return false;
    profileLine("export: plain surfaces");
    for (const ScenePart& r : t.parts) {
        cancelCheckpoint();
        if (!addPart(c, r, err)) return false;
    }
    profileLine("export: parts");
    if (const int part = closeupPartOf(t.parts, options.closeupPart); part >= 0) addCloseupCameras(c, part);
    if (source.fibres && !addFibres(c, *source.fibres, err)) return false;
    profileLine("export: fibres");
    scene.sortParts();
    cancelCheckpoint();
    if (!readStrandsBack(scene, err)) return false;
    profileLine("export: read-back");
    return true;
}

namespace {

// The shading of a map set whose files are at <package>/<prefix><name>.
PlyShading plyShadingOf(const PlyModeOptions& o, const PlyMapSet& set, const std::string& prefix) {
    PlyShading ps;
    ps.on = true;
    ps.uvRepeat = set.uvRepeat;  // the maps' u period is one ply pitch
    ps.uvWrap = 64;
    ps.frame = o.shading == PlyModeOptions::Shading::Bump ? PlyShading::Frame::Bump : PlyShading::Frame::Normal;
    ps.normal = prefix + (o.shading == PlyModeOptions::Shading::NormalFlipGreen ? set.normalFlipGreen : set.normal);
    ps.height = prefix + set.height;
    ps.cavity = o.cavity ? prefix + set.cavity : "";
    ps.bumpScale = set.bumpScale;
    return ps;
}

} // namespace

std::optional<PlyShading> attachPlySet(const PlyMapSetResult& r, const PlyModeOptions& o, const std::string& packageDir, double sigma,
                                       std::string* err) {
    namespace fs = std::filesystem;
    const std::string dst = packageDir + "/textures";
    std::error_code ec;
    fs::create_directories(dst, ec);
    PlyMapSet set = r.set;
    // The maps the material can use. The AO map stays in the cache (the viewer's preview reads it; the
    // renderer's material has no AO slot), and the cavity map is copied only when the material uses it.
    std::vector<std::string> files = {set.height, set.depth, set.normal, set.normalFlipGreen};
    if (o.cavity) files.push_back(set.cavity);
    for (const std::string& f : files) {
        fs::copy_file(r.dir + "/" + f, dst + "/" + f, fs::copy_options::overwrite_existing, ec);
        if (ec) {
            if (err) *err = "map set: cannot copy " + f + " to " + dst + ": " + ec.message();
            return std::nullopt;
        }
    }
    PlyMapSpec scene = o.spec;
    scene.sigma = sigma;
    set.bumpScale = (set.prof.rhoMax - set.prof.rhoMin) / sigma;
    writePlyMapsJson(dst, scene, set);
    return plyShadingOf(o, set, "textures/");
}

std::optional<PlyShading> setupPlyTextures(const PlyModeOptions& o, const std::string& packageDir, double sigma, std::string* err) {
    if (!o.textured()) return PlyShading{};
    PlyMapSetResult cachedSet;
    if (!ensurePlyMapSet(plyMapCacheRoot(), o.spec, cachedSet, err)) return std::nullopt;
    const std::optional<PlyShading> shading = attachPlySet(cachedSet, o, packageDir, sigma, err);
    if (!shading) return std::nullopt;
    const PlyShadingProfile& prof = cachedSet.set.prof;
    std::printf("ply textures (%s%s): %d plies, pitch %.2f mm, map D %.2f mm, %s, rho [%.4f, %.4f] mm, bump_scale %.6g; "
                "tube uvRepeat %.3f mm, uvWrap %d\n",
                o.shading == PlyModeOptions::Shading::Bump ? "bump" : o.shading == PlyModeOptions::Shading::Normal ? "normal" : "normal, green flipped",
                o.cavity ? " + cavity" : "", o.spec.g.count, o.spec.g.pitch, o.spec.yarnD, prof.gap ? "gap" : "crease", prof.rhoMin,
                prof.rhoMax, shading->bumpScale, shading->uvRepeat, shading->uvWrap);
    return shading;
}

} // namespace knit
