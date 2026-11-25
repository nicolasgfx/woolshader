// The photon_path_tracer package's scene.json: the one place the renderer's scene format is written
// (schema 12, one concrete type per material).
//
// Tubes are triangle meshes, the renderer has no curve primitives: one binary little-endian PLY per
// part with x y z nx ny nz u v (u = s / uvRepeat along the yarn, less an integer multiple of uvWrap
// if set, v = theta / 2 pi). Fibres are strand parts (.pptstrands, round linear;
// fibres/strands_file.h) after the meshes. scene.json is a v2 PhotonScene, or v5 when it has strand
// parts (the renderer requires v5 for root strands[]). One instance per mesh or strand part, all
// with identity transforms, placed in a groups[] tree per source object and fibre type when the
// caller places them. The default camera frames every part; a grey ground quad lies under them, lit
// by the uniform analytical sky.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <ostream>
#include <set>
#include <sstream>

#include "export/export_scene.h"
#include "io/files.h"

namespace knit {

namespace {

// The two precisions of a colour in the file: the OBJ follows them (colourAsWritten).
constexpr const char* kHairColourFormat = "[%.6f, %.6f, %.6f]";
constexpr const char* kAlbedoFormat = "[%.3f, %.3f, %.3f]";

std::string tripleText(const char* format, const Vec3& c) {
    char b[128];
    std::snprintf(b, sizeof b, format, c.x, c.y, c.z);
    return b;
}

std::string vec(const Vec3& v) { return tripleText(kHairColourFormat, v); }

// JSON numbers as the first exporter's literals printed them: shortest form, ".0" on integers.
std::string number(double v) {
    char b[40];
    std::snprintf(b, sizeof b, "%.6g", v);
    std::string s(b);
    if (s.find_first_of(".eEn") == std::string::npos) s += ".0";
    return s;
}

// Part entries of one list ("meshes" or "strands").
void writePartList(std::ostream& o, const char* key, const std::vector<PartFile>& parts) {
    o << "  \"" << key << "\": [\n";
    for (size_t i = 0; i < parts.size(); ++i)
        o << "    {\"name\": \"" << parts[i].name << "\", \"path\": \"" << parts[i].path << "\", \"material\": \"" << parts[i].material
          << "\"}" << (i + 1 < parts.size() ? "," : "") << "\n";
    o << "  ],\n";
}

void writeCamera(std::ostream& o, const char* id, const ExportCamera& camera) {
    o << "    {\"id\": \"" << id << "\", \"name\": \"" << camera.name << "\", \"position\": " << vec(camera.eye)
      << ", \"look_at\": " << vec(camera.lookAt) << ", \"up\": [0, 1, 0], \"fov_deg\": " << camera.fov
      << ", \"width\": 1600, \"height\": 1000, \"spp\": 64, \"max_bounces\": 8, \"view\": {\"kind\": \"look_at\", \"eye\": "
      << vec(camera.eye) << ", \"look_at\": " << vec(camera.lookAt) << ", \"up\": [0, 1, 0]}}";
}

// The renderer's groups[]: one root, a group per placed source object (by key), then a child group
// per (object, group) that holds a part; ids g0001.. in that order. Unplaced parts stay at the root.
struct GroupTree {
    std::map<std::string, Placement> placed;                 // part -> placement, the ground's included
    std::map<std::string, std::string> objectName;           // object key -> display name
    std::map<std::string, Box> objectBounds;
    std::string rootGroup;
    std::map<std::string, std::string> objectGroup;          // object key -> group id
    std::map<std::pair<std::string, std::string>, std::string> childGroup;  // (object key, group) -> group id
    int lastId = 0;

    const Placement* placementOf(const PartFile& part) const {
        const auto p = placed.find(part.name);
        return p != placed.end() && !p->second.object.empty() ? &p->second : nullptr;
    }
    std::string nextId() {
        char b[16];
        std::snprintf(b, sizeof b, "g%04d", ++lastId);
        return b;
    }
};

// A +Y quad (counter-clockwise from above) just below the parts: the offset keeps a flat bottom off
// the plane (no coplanar faces). A floor the renderer's editor frames, picks and deletes as an
// ordinary object (50 diagonals made a ~775 m quad for the teddies: frame-all left the bears specks).
bool writeGroundMesh(const ExportScene& scene, const SceneFrame& frame, std::vector<PartFile>& meshList,
                     std::string* err) {
    const Vec3& c = frame.centre;
    const double extent = 2 * frame.radius, y = scene.bounds().lo.y - 1e-4 * extent;
    const double h = ExportScene::kGroundHalfSize * extent;
    SceneMesh quad;
    quad.pos = {{c.x - h, y, c.z - h}, {c.x + h, y, c.z - h}, {c.x + h, y, c.z + h}, {c.x - h, y, c.z + h}};
    quad.nrm.assign(4, Vec3{0, 1, 0});
    quad.tri = {0, 2, 1, 0, 3, 2};
    const std::string file = std::string("meshes/") + ExportScene::kGround + ".ply";
    if (!writePartMesh(scene.folder() + "/" + file, quad)) return fail(err, "cannot write " + file);
    meshList.push_back({ExportScene::kGround, file, ExportScene::kGround});
    return true;
}

// A hair material: hair is its only layer.
void writeHairMaterial(std::ostringstream& t, const MaterialSpec& sp) {
    t << "\"type\": \"hair\", \"fibre\": {\"color\": " << tripleText(kHairColourFormat, sp.color)
      << ", \"longitudinal_roughness\": " << number(sp.fibreLongRoughness) << ", \"azimuthal_roughness\": " << number(sp.fibreAzimRoughness)
      << ", \"ior\": " << number(sp.fibreIor) << ", \"cuticle_tilt_deg\": " << number(sp.fibreCuticleDeg) << "}}";
}

// The maps of a ply-mapped yarn. The cavity multiplier is the diffuse block's albedo_tex (sRGB,
// multiplied by the albedo); the normal or height map is material-level (linear). Schema 12 has no AO map.
void writePlyMapKeys(std::ostringstream& t, const PlyShading& ps) {
    if (ps.frame == PlyShading::Frame::Normal && !ps.normal.empty())
        t << ", \"normal_tex\": \"" << ps.normal << "\", \"normalmap_scale\": " << ps.normalScale;
    if (ps.frame == PlyShading::Frame::Bump && !ps.height.empty()) {
        char b[32];
        std::snprintf(b, sizeof b, "%.8g", ps.bumpScale);
        t << ", \"bump_tex\": \"" << ps.height << "\", \"bump_scale\": " << b;
    }
}

// A fabric's sheen block and, when the cloth lets light through, its translucency block (the renderer's
// TranslucencyParams: fabric only), each with its leading ", ". The sheen layer's weight is its coverage.
void writeSheenKeys(std::ostringstream& t, const MaterialSpec& sp) {
    t << ", \"sheen\": {\"color\": [1.0, 1.0, 1.0], \"tint\": " << number(sp.sheenTint) << ", \"roughness\": " << number(sp.sheenRoughness);
    if (sp.sheen != 1) t << ", \"coverage\": " << number(sp.sheen);
    t << "}";
    if (sp.translucency > 0) t << ", \"translucency\": {\"transmit_share\": " << number(sp.translucency) << "}";
}

// The highlight of a baked yarn tube, with its leading ", ": the renderer's specular block (the plastic's, and since
// photon_path_tracer dd8870fda also a fabric's). Its x axis follows the measured fibre direction; the roughness map
// is the value across the yarn (the authored roughness across is 1) and scales the one along it; the specular map,
// the share of rays that met a fibre, scales the reflectance (the renderer's MAT-32: a coverage read as an F0 factor).
void writeBakedHighlight(std::ostringstream& t, const BakedMaps& b) {
    t << ", \"specular\": {\"specular\": " << number(b.specularScale) << ", \"roughness\": " << number(b.roughnessAlong)
      << ", \"roughness_y\": 1.0, \"aniso_rotation\": " << number(b.anisoRotation) << ", \"roughness_tex\": \"" << b.roughness
      << "\", \"reflectance_tex\": \"" << b.specular << "\"}";
}

// A baked yarn tube (fibres/yarn_bake.h), the dye under its grey cavity map (albedo_tex) and the normal map, in
// one of two materials (BakedMaps::fabric), both with the highlight of the roughness and specular maps:
//   fabric  the yarn's own fabric: its sheen, the highlight under it, its translucency. Always a fabric, the only
//           type with the highlight and translucency: a yarn without a sheen gets the least one (0.01), as the
//           core proxy does. The body and its glow receive only what the highlight does not reflect.
//   shiny   plastic: the highlight alone over the body, brighter (no sheen in front of it). Opaque.
void writeBakedYarnMaterial(std::ostringstream& t, const ExportScene& scene, const ExportMaterial& m) {
    const BakedMaps& b = scene.bakedMaps.at(m.id);
    if (b.fabric) {
        MaterialSpec sp = m.spec;
        sp.sheen = std::max(sp.sheen, 0.01);
        t << "\"type\": \"fabric\", \"diffuse\": {\"albedo\": " << tripleText(kAlbedoFormat, sp.color) << ", \"albedo_tex\": \"" << b.cavity
          << "\"}";
        writeSheenKeys(t, sp);
        writeBakedHighlight(t, b);
        t << ", \"normal_tex\": \"" << b.normal << "\", \"normalmap_scale\": 1.0}";
        return;
    }
    t << "\"type\": \"plastic\", \"diffuse\": {\"albedo\": " << tripleText(kAlbedoFormat, m.spec.color) << ", \"albedo_tex\": \"" << b.cavity
      << "\"}";
    writeBakedHighlight(t, b);
    t << ", \"normal_tex\": \"" << b.normal << "\", \"normalmap_scale\": 1.0}";
}

// One scene.json material entry, without the separator: the one place the material schema is
// written. Schema 12, one concrete type per material (the renderer's RDR-320): hair for a fibre
// material, else fabric (diffuse + sheen), coated_diffuse (diffuse + coat) or diffuse, by the
// rule of the renderer's material migration script (migrate_materials.py).
void writeMaterial(std::ostream& o, const ExportScene& scene, const ExportMaterial& m) {
    const MaterialSpec& sp = m.spec;
    std::ostringstream t;
    t << "    {\"id\": \"" << m.id << "\", \"name\": \"" << m.id << "\", \"scope\": \"scene\", \"schema_version\": 12, "
      << "\"provenance\": {\"kind\": \"created_in_scene\"}, ";
    if (m.kind == MaterialKind::Hair) {
        writeHairMaterial(t, sp);
        o << t.str();
        return;
    }
    if (m.kind == MaterialKind::BakedYarn) {
        writeBakedYarnMaterial(t, scene, m);
        o << t.str();
        return;
    }
    const bool plyMapped = m.kind == MaterialKind::PlyMappedYarn;
    const bool sheen = sp.sheen > 0, coat = sp.clearcoat && !sheen;
    t << "\"type\": \"" << (sheen ? "fabric" : coat ? "coated_diffuse" : "diffuse") << "\", ";
    t << "\"diffuse\": {\"albedo\": " << tripleText(kAlbedoFormat, sp.color);
    if (plyMapped && !scene.plyShading.cavity.empty()) t << ", \"albedo_tex\": \"" << scene.plyShading.cavity << "\"";
    else if (!m.colourTexture.empty()) t << ", \"albedo_tex\": \"" << m.colourTexture << "\"";
    t << "}";
    if (sheen) writeSheenKeys(t, sp);
    if (coat) {
        t << ", \"coat\": {\"roughness\": " << number(sp.clearcoatRoughness) << ", \"ior\": " << number(sp.clearcoatIor);
        if (sp.clearcoatWeight != 1) t << ", \"coverage\": " << number(sp.clearcoatWeight);
        t << "}";
    }
    if (plyMapped) writePlyMapKeys(t, scene.plyShading);
    t << "}";
    o << t.str();
}

// With a scene graph the ground is an object of its own (top-level group "Ground"): the renderer
// selects, deletes and (Import Scene as Asset, per top-level group) imports it apart. One scene node
// (projectName) over every object and the ground: the package imports as one object.
GroupTree buildGroups(const ExportScene& scene, const std::vector<const PartFile*>& parts) {
    GroupTree tree;
    tree.placed = scene.placements();
    bool anyPlaced = false;
    for (const auto& [part, at] : scene.placements()) anyPlaced = anyPlaced || !at.object.empty();
    if (scene.groundPlane && anyPlaced) tree.placed[ExportScene::kGround] = {ExportScene::kGroundObject, "Ground", ""};
    std::set<std::pair<std::string, std::string>> children;  // (object key, group)
    for (const PartFile* e : parts) {
        const Placement* p = tree.placementOf(*e);
        if (!p) continue;
        tree.objectName.emplace(p->object, p->objectName.empty() ? p->object : p->objectName);
        if (!p->group.empty()) children.insert({p->object, p->group});
        if (const auto b = scene.partBounds().find(e->name); b != scene.partBounds().end()) tree.objectBounds[p->object].grow(b->second);
    }
    if (tree.objectName.empty()) return tree;
    tree.rootGroup = tree.nextId();
    for (const auto& [object, name] : tree.objectName)
        if (object != ExportScene::kGroundObject) tree.objectGroup[object] = tree.nextId();
    for (const auto& child : children) tree.childGroup[child] = tree.nextId();
    if (tree.objectName.count(ExportScene::kGroundObject)) tree.objectGroup[ExportScene::kGroundObject] = tree.nextId();  // last
    return tree;
}

void writeInstances(std::ostream& o, const ExportScene& scene, const std::vector<const PartFile*>& parts,
                    size_t meshCount, const GroupTree& tree) {
    o << "  \"instances\": [\n";
    for (size_t i = 0; i < parts.size(); ++i) {
        char id[32];
        std::snprintf(id, sizeof id, "inst_%04zu", i);
        o << "    {\"id\": \"" << id << "\", \"name\": \"" << parts[i]->name << "\"";
        if (const Placement* p = tree.placementOf(*parts[i]))
            o << ", \"group_id\": \"" << (p->group.empty() ? tree.objectGroup.at(p->object) : tree.childGroup.at({p->object, p->group})) << "\"";
        o << ", \"asset_mesh_index\": " << i
          << ", \"material\": {\"kind\": \"scene_material\", \"id\": \"" << parts[i]->material
          << "\"}, \"transform\": {\"translate\": [0, 0, 0], \"rotate_euler_deg\": [0, 0, 0], \"scale\": [1, 1, 1]}";
        // The wool look's thin fibre ends. The renderer applies this shape to the exported radii at
        // compile, so the radii stay as they are, or the taper would apply twice. Its key is the
        // part's own asset_mesh_index; root and tip diameters are unused with "imported" thickness.
        if (i >= meshCount && scene.isWoolFibreMaterial(parts[i]->material))
            o << ", \"part_fibre_shapes\": {\"" << i << "\": {\"thickness\": \"imported\", \"diameter_scale\": 1.0, \"tip_scale\": "
              << scene.fibreTipScale << ", \"root_diameter_m\": 7e-05, \"tip_diameter_m\": 7e-05, \"diameter_variation\": 0.0, \"seed\": "
              << 1000 + i << "}}";
        o << "}" << (i + 1 < parts.size() ? "," : "") << "\n";
    }
}

void writeCameras(std::ostream& o, const ExportScene& scene, const SceneFrame& frame) {
    o << "  ],\n  \"cameras\": [\n";
    writeCamera(o, "cam_0000", {"front", frame.eye, frame.centre, frame.fov});
    for (size_t i = 0; i < scene.extraCameras.size(); ++i) {
        char id[32];
        std::snprintf(id, sizeof id, "cam_%04zu", i + 1);
        o << ",\n";
        writeCamera(o, id, scene.extraCameras[i]);
    }
    o << "\n  ],\n  \"default_camera_id\": \"cam_0000\",\n";
}

void writeGroups(std::ostream& o, const ExportScene& scene, const GroupTree& tree) {
    if (tree.objectName.empty()) return;
    // The pivot of a group is its object's bottom centre (the editor rotates and scales a group
    // about it: a bear about its feet); child groups share it.
    const auto bottomCentre = [](const Box& box) { return vec({0.5 * (box.lo.x + box.hi.x), box.lo.y, 0.5 * (box.lo.z + box.hi.z)}); };
    const auto pivotOf = [&](const std::string& object) {
        const auto b = tree.objectBounds.find(object);
        return b == tree.objectBounds.end() ? std::string("[0, 0, 0]") : bottomCentre(b->second);
    };
    std::vector<std::string> rows;
    const auto addRow = [&](const std::string& id, const std::string& name, const std::string& parent, const std::string& pivot) {
        rows.push_back("    {\"id\": \"" + id + "\", \"name\": \"" + jsonEscaped(name) + "\"" +
                       (parent.empty() ? "" : ", \"parent_id\": \"" + parent + "\"") +
                       ", \"origin\": \"User\", \"is_light\": false, \"primary_camera_film\": false, \"pivot\": " + pivot + "}");
    };
    addRow(tree.rootGroup, scene.projectName, "", bottomCentre(scene.bounds()));
    for (const auto& [object, id] : tree.objectGroup) addRow(id, tree.objectName.at(object), tree.rootGroup, pivotOf(object));
    for (const auto& [child, id] : tree.childGroup) addRow(id, child.second, tree.objectGroup.at(child.first), pivotOf(child.first));
    std::sort(rows.begin(), rows.end());  // by id (g%04d)
    o << "  \"groups\": [\n";
    for (size_t i = 0; i < rows.size(); ++i) o << rows[i] << (i + 1 < rows.size() ? "," : "") << "\n";
    // The renderer mints new group ids from next_group_seq (g%04d), so it follows the last.
    o << "  ],\n  \"next_group_seq\": " << tree.lastId + 1 << ",\n";
}

// The light: the renderer's analytical env is one constant colour (no sun; there is no directional
// light type). Its sky_* keys are read only with sky_color present, `intensity` scales HDRIs only,
// and lower_hemisphere defaults to black for it, so all are written.
void writeEnvironment(std::ostream& o, const ExportScene& scene) {
    char sky[32];
    std::snprintf(sky, sizeof sky, scene.skyIntensity == 1.0 ? "1.0" : "%.6g", scene.skyIntensity);
    o << "  \"env\": {\"mode\": \"analytical\", \"sky_enabled\": true, \"sky_color\": [1.0, 1.0, 1.0], \"sky_intensity\": " << sky
      << ", \"intensity\": 1.0, \"lower_hemisphere\": \"full\", \"visible_to_camera\": true},\n";
    o << "  \"render\": {\"width\": 1600, \"height\": 1000, \"spp\": 64, \"max_bounces\": 8}\n}\n";
}

} // namespace

Vec3 colourAsWritten(const ExportMaterial& material) {
    const char* format = material.kind == MaterialKind::Hair ? kHairColourFormat : kAlbedoFormat;
    std::istringstream text(tripleText(format, material.spec.color));  // "[x, y, z]"
    Vec3 written;
    char bracket = 0, comma = 0;
    text >> bracket >> written.x >> comma >> written.y >> comma >> written.z;
    return written;
}

bool writeSceneJson(const ExportScene& scene, std::string* err) {
    const std::optional<SceneFrame> frame = frameScene(scene, err);
    if (!frame) return false;
    std::vector<PartFile> meshList = scene.meshParts();
    if (scene.groundPlane && !writeGroundMesh(scene, *frame, meshList, err)) return false;
    const std::optional<std::vector<ExportMaterial>> materials = scene.materials(err);
    if (!materials) return false;

    // One stream for the whole file: its default number format is part of the file.
    std::ostringstream o;
    const std::vector<PartFile>& strands = scene.strandParts();
    const int version = strands.empty() ? 2 : 5;  // root strands[] needs format version 5
    o << "{\n  \"format\": {\"name\": \"PhotonScene\", \"version\": " << version << ", \"minimum_reader\": " << version
      << ", \"generator\": {\"name\": \"wool_shader knit_check --export-ppt\"}},\n";
    o << "  \"units\": \"meters\",\n  \"up_axis\": \"Y\",\n  \"project_name\": \"" << scene.projectName << "\",\n";
    o << "  \"materials\": [\n";
    for (size_t i = 0; i < materials->size(); ++i) {
        writeMaterial(o, scene, (*materials)[i]);
        o << (i + 1 < materials->size() ? "," : "") << "\n";
    }
    o << "  ],\n";
    writePartList(o, "meshes", meshList);
    if (!strands.empty()) writePartList(o, "strands", strands);
    // One instance per part; meshes and strands share one part index space, meshes first.
    std::vector<const PartFile*> parts;
    for (const PartFile& e : meshList) parts.push_back(&e);
    for (const PartFile& e : strands) parts.push_back(&e);
    const GroupTree tree = buildGroups(scene, parts);
    writeInstances(o, scene, parts, meshList.size(), tree);
    writeCameras(o, scene, *frame);
    writeGroups(o, scene, tree);
    writeEnvironment(o, scene);

    std::ofstream f(scene.folder() + "/scene.json");
    if (!f) return fail(err, "cannot write scene.json");
    f << o.str();
    return static_cast<bool>(f);
}

} // namespace knit
