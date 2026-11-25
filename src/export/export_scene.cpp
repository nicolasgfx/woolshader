// ExportScene: the parts on disk, the records beside them, and the materials resolved. See export_scene.h.
#include "export/export_scene.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <set>
#include <tuple>

#include "io/files.h"
#include "model/cancel.h"
#include "model/colour.h"

namespace knit {

namespace {

constexpr double kMinHairColour = 1e-4;  // the renderer reads no hair colour below this

// A refusal from a function that returns an optional.
std::nullopt_t refuse(std::string* err, const std::string& why) {
    (void)fail(err, why);
    return std::nullopt;
}

// --- hair colours: what a hair material must say so that fibres show a yarn's colour

// A hair colour the renderer reads: every channel at least its floor.
Vec3 atHairColourFloor(const Vec3& colour) {
    return {std::max(colour.x, kMinHairColour), std::max(colour.y, kMinHairColour), std::max(colour.z, kMinHairColour)};
}

// The hair colour channel that makes one fibre show `channel`. The renderer's hair colour is that of
// a dense mass, so one thin fibre of that colour is almost clear. A ray through the fibre's centre
// crosses two radii, so keeping the fraction c takes sigma_a = -ln(c) / 2 per radius; the colour of
// that absorption is exp(-sqrt(sigma_a) den(beta_n)), Chiang et al. 2016 Eq. 9, as the renderer's
// fibre_absorption.h.
double hairChannelOfFibre(double channel, double azimuthalRoughness) {
    const double b = azimuthalRoughness;
    const double den = 5.969 - 0.215 * b + 2.532 * b * b - 10.73 * b * b * b + 5.574 * b * b * b * b + 0.245 * b * b * b * b * b;
    const double sigmaA = -std::log(std::clamp(channel, kMinHairColour, 1.0)) / 2;
    return std::max(std::exp(-std::sqrt(sigmaA) * den), kMinHairColour);
}

Vec3 hairColourOfFibre(const Vec3& colour, double azimuthalRoughness) {
    return {hairChannelOfFibre(colour.x, azimuthalRoughness), hairChannelOfFibre(colour.y, azimuthalRoughness),
            hairChannelOfFibre(colour.z, azimuthalRoughness)};
}

// Between the renderer's dense mass (the colour as it stands) and one fibre of it: the geometric mean.
Vec3 hairColourOfLayer(const Vec3& colour, double azimuthalRoughness) {
    const Vec3 one = hairColourOfFibre(colour, azimuthalRoughness);
    return {std::sqrt(colour.x * one.x), std::sqrt(colour.y * one.y), std::sqrt(colour.z * one.z)};
}

// A dense coat in its yarn's colour. Hair adds an untinted shine to every channel, which washes out
// channels near zero: red turned salmon and black grey on the car. So the shine comes off first, never
// more than the lowest channel: the same amount off all three keeps the hue. Pale wool hardly changes.
Vec3 hairColourOfCoat(const Vec3& colour) {
    constexpr double kShine = 0.06;  // matched on the car (red, tyres) and on cream and tan wool
    const double shine = std::min(kShine, std::min({colour.x, colour.y, colour.z}));
    return {std::max(colour.x - shine, kMinHairColour), std::max(colour.y - shine, kMinHairColour),
            std::max(colour.z - shine, kMinHairColour)};
}

// --- materials

// A material as the formats see it. Hair has no diffuse layer to hold a texture, and on a ply-mapped
// yarn the cavity map takes the base colour's place.
ExportMaterial resolvedMaterial(const ExportScene& scene, const std::string& id, const MaterialSpec& spec,
                                MaterialKind kind) {
    ExportMaterial m{id, spec, kind, ""};
    const bool hasCavityMap = kind == MaterialKind::PlyMappedYarn && !scene.plyShading.cavity.empty();
    const auto texture = scene.albedoTextures.find(id);
    if (kind != MaterialKind::Hair && !hasCavityMap && texture != scene.albedoTextures.end())
        m.colourTexture = texture->second;
    // No type has a coat and a sheen: a coated sheen material is fabric, its coat dropped (as the
    // renderer's migration drops it).
    if (kind != MaterialKind::Hair && spec.clearcoat && spec.sheen > 0)
        std::fprintf(stderr, "warning: material %s: the renderer has no coat over a sheen; written as fabric without the coat\n",
                     id.c_str());
    return m;
}

MaterialKind kindOfSpec(const MaterialSpec& spec) {
    return spec.model == MaterialSpec::Model::Fibre ? MaterialKind::Hair : MaterialKind::Surface;
}

bool finite(const Vec3& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

// A mesh's vertices in Y-up with their normals, triangles as they stand.
SceneMesh sceneMeshOf(const TriMesh& sceneMesh) {
    TriMesh m = sceneMesh;
    m.computeNormals();
    SceneMesh out;
    out.pos.resize(m.pos.size());
    out.nrm.resize(m.pos.size());
    for (size_t i = 0; i < m.pos.size(); ++i) {
        out.pos[i] = yUp(m.pos[i]);
        out.nrm[i] = yUp(m.vnormal[i]);
    }
    for (const auto& t : m.tri) out.tri.insert(out.tri.end(), {uint32_t(t[0]), uint32_t(t[1]), uint32_t(t[2])});
    return out;
}

Box boxOf(const std::vector<Vec3>& points) {
    Box box;
    for (const Vec3& p : points) box.grow(p);
    return box;
}

} // namespace

Vec3 yUp(const Vec3& p) { return {p.x, p.z, -p.y}; }

bool writePartMesh(const std::string& path, const SceneMesh& mesh) {
    std::ofstream f(path, std::ios::binary);
    if (!f) return false;
    const std::vector<Vec3>& pos = mesh.pos;
    const std::vector<Vec3>& nrm = mesh.nrm;
    const bool hasUv = mesh.uv.size() == 2 * pos.size();
    f << "ply\nformat binary_little_endian 1.0\nelement vertex " << pos.size()
      << "\nproperty float x\nproperty float y\nproperty float z\nproperty float nx\nproperty float ny\nproperty float nz\n";
    if (hasUv) f << "property float u\nproperty float v\n";
    f << "element face " << mesh.tri.size() / 3 << "\nproperty list uchar uint vertex_indices\nend_header\n";
    for (size_t i = 0; i < pos.size(); ++i) {
        float v[8] = {float(pos[i].x), float(pos[i].y), float(pos[i].z), float(nrm[i].x), float(nrm[i].y), float(nrm[i].z),
                      hasUv ? mesh.uv[2 * i] : 0.f, hasUv ? mesh.uv[2 * i + 1] : 0.f};
        f.write(reinterpret_cast<const char*>(v), (hasUv ? 8 : 6) * sizeof(float));
    }
    for (size_t t = 0; t + 2 < mesh.tri.size(); t += 3) {
        unsigned char n = 3;
        f.write(reinterpret_cast<const char*>(&n), 1);
        f.write(reinterpret_cast<const char*>(&mesh.tri[t]), 3 * sizeof(uint32_t));
    }
    return static_cast<bool>(f);
}

std::optional<StoredMesh> readPartMesh(const std::string& path, std::string* err) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return refuse(err, "cannot open " + path);
    std::string line;
    size_t vertices = 0, faces = 0;
    int properties = 0;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind("element vertex ", 0) == 0) vertices = std::stoull(line.substr(15));
        else if (line.rfind("element face ", 0) == 0) faces = std::stoull(line.substr(13));
        else if (line.rfind("property float ", 0) == 0) ++properties;
        else if (line == "end_header") break;
    }
    if (properties != 6 && properties != 8) return refuse(err, path + ": unexpected vertex layout");
    StoredMesh mesh;
    mesh.stride = properties;
    mesh.vertices.resize(vertices * properties);
    f.read(reinterpret_cast<char*>(mesh.vertices.data()), static_cast<std::streamsize>(mesh.vertices.size() * sizeof(float)));
    mesh.tri.resize(faces * 3);
    for (size_t t = 0; t < faces; ++t) {
        unsigned char n = 0;
        f.read(reinterpret_cast<char*>(&n), 1);
        if (n != 3) return refuse(err, path + ": a face is not a triangle");
        f.read(reinterpret_cast<char*>(&mesh.tri[3 * t]), 3 * sizeof(uint32_t));
        if ((t & 0xffff) == 0) cancelCheckpoint();
    }
    if (!f) return refuse(err, path + ": truncated");
    for (const uint32_t i : mesh.tri)
        if (i >= vertices) return refuse(err, path + ": a face index is out of range");
    return mesh;
}

std::vector<BaseMaterial> defaultMaterials() {
    const auto material = [](const char* id, const Vec3& colour, bool sheen) {
        MaterialSpec m;
        m.color = colour;
        m.sheen = sheen ? 1.0 : 0.0;
        m.sheenTint = kYarnSheenTint;
        return BaseMaterial{id, m};
    };
    return {material("wool_white", kDefaultYarnColour, true), material("wool_brown", {0.42, 0.20, 0.08}, true),
            material("black", {0.02, 0.02, 0.02}, false), material("bow", {0.57, 0.43, 0.29}, true)};
}

MaterialSpec defaultRandomFibre() {
    MaterialSpec m;
    m.color = {ExportScene::kRandomFibreGrey, ExportScene::kRandomFibreGrey, ExportScene::kRandomFibreGrey};
    return m;
}

void Box::grow(const Vec3& p) {
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
}

void Box::grow(const Box& b) {
    lo = {std::min(lo.x, b.lo.x), std::min(lo.y, b.lo.y), std::min(lo.z, b.lo.z)};
    hi = {std::max(hi.x, b.hi.x), std::max(hi.y, b.hi.y), std::max(hi.z, b.hi.z)};
}

ExportScene::ExportScene(std::string folder) : folder_(std::move(folder)) {}

bool ExportScene::addMeshFile(const std::string& name, const std::string& material, const Box& box, const SceneMesh& mesh,
                              std::string* err) {
    const std::string file = "meshes/" + name + ".ply";
    if (!writePartMesh(folder_ + "/" + file, mesh)) return fail(err, "cannot write " + file);
    bounds_.grow(box);
    partBounds_[name] = box;
    meshes_.push_back({name, file, material});
    triangles_ += mesh.tri.size() / 3;
    return true;
}

bool ExportScene::addMeshPart(const std::string& name, const std::string& material, const SceneMesh& mesh,
                              std::string* err) {
    return addMeshFile(name, material, boxOf(mesh.pos), mesh, err);
}

bool ExportScene::addMesh(const std::string& name, const std::string& material, const TriMesh& sceneMesh,
                          std::string* err) {
    return addMeshPart(name, material, sceneMeshOf(sceneMesh), err);
}

bool ExportScene::addTexturedMesh(const std::string& name, const std::string& material, const TriMesh& sceneMesh,
                                  const std::vector<std::array<Vec2, 3>>& cornerUv, std::string* err) {
    if (cornerUv.size() != sceneMesh.tri.size()) return fail(err, name + ": one set of corner uvs per triangle");
    const SceneMesh whole = sceneMeshOf(sceneMesh);
    // One output vertex per distinct (vertex, uv): shared where the uvs agree, split on seams.
    SceneMesh split;
    std::map<std::tuple<int, float, float>, uint32_t> index;
    for (size_t t = 0; t < sceneMesh.tri.size(); ++t)
        for (int k = 0; k < 3; ++k) {
            const int v = sceneMesh.tri[t][k];
            const float u = float(cornerUv[t][k].x), w = float(cornerUv[t][k].y);
            const auto [it, added] = index.emplace(std::make_tuple(v, u, w), uint32_t(split.pos.size()));
            if (added) {
                split.pos.push_back(whole.pos[v]);
                split.nrm.push_back(whole.nrm[v]);
                split.uv.insert(split.uv.end(), {u, w});
            }
            split.tri.push_back(it->second);
        }
    return addMeshFile(name, material, boxOf(whole.pos), split, err);
}

std::optional<std::vector<std::string>> ExportScene::addTubeGroups(const std::string& name, const ExportTube& tube,
                                                                   const TubeGroups& groups, std::string* err) {
    const SceneMesh& whole = tube.mesh;
    std::map<int, std::vector<uint32_t>> byGroup;  // group -> triangles (corner indices into the tube)
    for (size_t t = 0; t + 2 < whole.tri.size(); t += 3) {
        const uint32_t a = whole.tri[t], b = whole.tri[t + 1], c = whole.tri[t + 2];
        const int g = groups.groupOfTriangle[t / 3];
        if (g < 0) return refuse(err, name + ": a tube triangle without a group (yarn " + std::to_string(tube.yarn[a]) + ")");
        std::vector<uint32_t>& v = byGroup[g];
        v.insert(v.end(), {a, b, c});
    }
    std::vector<std::string> written;
    for (const auto& [g, triangles] : byGroup) {
        // The group's vertices, compacted.
        std::map<uint32_t, uint32_t> index;
        SceneMesh part;
        part.tri.reserve(triangles.size());
        for (uint32_t v : triangles) {
            auto it = index.find(v);
            if (it == index.end()) {
                it = index.emplace(v, static_cast<uint32_t>(part.pos.size())).first;
                part.pos.push_back(whole.pos[v]);
                part.nrm.push_back(whole.nrm[v]);
                if (whole.uv.size() == 2 * whole.pos.size()) part.uv.insert(part.uv.end(), {whole.uv[2 * v], whole.uv[2 * v + 1]});
            }
            part.tri.push_back(it->second);
        }
        const std::string partName = name + groups.groups[g].suffix;
        if (!addMeshPart(partName, groups.groups[g].material, part, err)) return std::nullopt;
        written.push_back(partName);
    }
    return written;
}

bool ExportScene::addStrands(const std::string& name, const std::string& material, const StrandSetF& s,
                             std::string* err) {
    if (s.strands() == 0) return true;  // omitted: the renderer rejects empty strand parts
    const std::string file = "strands/" + name + ".pptstrands";
    if (!writeStrandsFile(folder_ + "/" + file, s, strandRules_, err)) return false;
    Box box;
    for (size_t i = 0; i + 2 < s.xyz.size(); i += 3) box.grow(Vec3{s.xyz[i], s.xyz[i + 1], s.xyz[i + 2]});
    bounds_.grow(box);
    partBounds_[name] = box;
    strands_.push_back({name, file, material});
    strandSegments_ += s.vertices() - s.strands();
    return true;
}

std::optional<std::string> ExportScene::addColourTexture(const std::string& name, const ColourImage& image,
                                                         std::string* err) {
    const size_t pixels = size_t(std::max(image.width, 0)) * size_t(std::max(image.height, 0));
    if (image.width <= 0 || image.height <= 0 || image.rgba.size() != pixels * 4)
        return refuse(err, "colour texture " + name + ": unexpected size");
    std::vector<uint8_t> rgb(pixels * 3);
    for (size_t i = 0; i < pixels; ++i)
        for (int c = 0; c < 3; ++c) rgb[3 * i + c] = image.rgba[4 * i + c];
    const std::string path = kColourTexturePrefix + name + ".png";
    std::error_code ec;
    std::filesystem::create_directories(folder_ + "/textures", ec);
    if (!writePngRgb8(folder_ + "/" + path, image.width, image.height, rgb)) return refuse(err, "cannot write " + path);
    return path;
}

void ExportScene::place(const std::string& name, const Placement& at) { placements_[name] = at; }

std::vector<size_t> ExportScene::strandFileCuts(const std::vector<uint64_t>& first, long long maxSegments) {
    const size_t strands = first.empty() ? 0 : first.size() - 1;
    std::vector<size_t> cuts{0};
    long long segments = 0;  // in the file being filled
    for (size_t i = 0; i < strands; ++i) {
        const long long own = static_cast<long long>(first[i + 1] - first[i]) - 1;
        if (i > cuts.back() && segments + own > maxSegments) {
            cuts.push_back(i);
            segments = 0;
        }
        segments += own;
    }
    cuts.push_back(strands);
    return cuts;
}

std::string ExportScene::orientedFibreMaterial(const std::string& yarnId) {
    const std::string id = "fibre_oriented_" + yarnId;
    fibreMaterials_[id] = {yarnId, false};
    return id;
}

std::string ExportScene::randomFibreMaterial(const std::string& yarnId) {
    const std::string id = "fibre_random_" + yarnId;
    fibreMaterials_[id] = {yarnId, true};
    return id;
}

std::string ExportScene::bodyFibreMaterial(const std::string& yarnId) {
    const std::string id = "fibre_body_" + yarnId;
    const auto has = [&](const std::string& m) {
        return std::any_of(baseMaterials.begin(), baseMaterials.end(), [&](const BaseMaterial& b) { return b.id == m; });
    };
    const auto yarn = std::find_if(baseMaterials.begin(), baseMaterials.end(), [&](const BaseMaterial& b) { return b.id == yarnId; });
    if (!has(id) && yarn != baseMaterials.end()) {
        MaterialSpec m = yarn->spec;
        m.model = MaterialSpec::Model::Fibre;
        m.color = atHairColourFloor(yarn->spec.color);
        baseMaterials.push_back({id, m});
    }
    return id;
}

std::string ExportScene::bodyCoreMaterial(const std::string& yarnId, bool translucent) {
    const std::string id = (translucent ? "core_translucent_" : "core_") + yarnId;
    const auto yarn = std::find_if(baseMaterials.begin(), baseMaterials.end(), [&](const BaseMaterial& b) { return b.id == yarnId; });
    const bool known = std::any_of(baseMaterials.begin(), baseMaterials.end(), [&](const BaseMaterial& b) { return b.id == id; });
    if (!known && yarn != baseMaterials.end()) {
        MaterialSpec m = yarn->spec;
        m.model = MaterialSpec::Model::DiffuseSheen;
        m.clearcoat = false;
        if (translucent) {
            m.sheen = std::max(m.sheen, 0.01);  // fabric: the only type with translucency
            m.translucency = kBodyCoreTranslucency;
        } else {
            m.sheen = 0;  // plain diffuse: the fibres carry the gloss
            m.translucency = 0;
        }
        baseMaterials.push_back({id, m});
    }
    return id;
}

void ExportScene::addBakedYarn(const std::string& id, const MaterialSpec& yarnSpec, const BakedMaps& maps) {
    const bool known = std::any_of(baseMaterials.begin(), baseMaterials.end(), [&](const BaseMaterial& m) { return m.id == id; });
    if (!known) baseMaterials.push_back({id, yarnSpec});
    bakedMaps[id] = maps;
}

void ExportScene::sortParts() {
    const auto byName = [](const PartFile& a, const PartFile& b) { return a.name < b.name; };
    std::sort(meshes_.begin(), meshes_.end(), byName);
    std::sort(strands_.begin(), strands_.end(), byName);
    const auto cameraByName = [](const ExportCamera& a, const ExportCamera& b) { return a.name < b.name; };
    std::sort(extraCameras.begin(), extraCameras.end(), cameraByName);
}

std::optional<MaterialSpec> ExportScene::fibreMaterialSpec(const std::string& id, std::string* err) const {
    using Model = MaterialSpec::Model;
    const auto it = fibreMaterials_.find(id);
    if (it == fibreMaterials_.end()) return refuse(err, "fibre material " + id + ": not registered");
    const FibreMaterial& fm = it->second;
    const auto yarn = std::find_if(baseMaterials.begin(), baseMaterials.end(), [&](const BaseMaterial& m) { return m.id == fm.yarn; });
    const bool yarnFound = yarn != baseMaterials.end();
    const std::string noYarn = "fibre material " + id + ": no base material " + fm.yarn;
    // A hair BSDF from fibreLook.hair in the given colour (its model may still say diffuse when only
    // hairMaterials asks for hair, e.g. the strand demo).
    const auto hair = [&](const Vec3& colour) {
        MaterialSpec m = fibreLook.hair;
        m.model = Model::Fibre;
        m.color = colour;
        return m;
    };
    if (fm.random) {
        // Its own hair model (material.fibre_random.model = fibre: translucent stray fibres beside
        // diffuse oriented ones), or the oriented fibres' hair parameters with its own colour.
        const MaterialSpec& random = fibreLook.random;
        MaterialSpec out = random.model == Model::Fibre ? random : fibreLook.hairMaterials ? hair(random.color) : random;
        if (fibreLook.randomInheritsYarn) {
            if (!yarnFound) return refuse(err, noYarn);
            const Vec3 c = yarn->spec.color;
            // A stray is a fibre of that yarn and carries its dye, and a mass of hair in the yarn's
            // colour reads as that yarn (the renderer's hair colour is a dense mass's).
            out.color = out.model == Model::Fibre ? atHairColourFloor(c) : c;
        }
        return out;
    }
    if (fibreLook.orientedOverride) return *fibreLook.orientedOverride;
    if (!yarnFound) return refuse(err, noYarn);
    if (!fibreLook.hairMaterials) return yarn->spec;
    // Hair: the fuzz lies between the renderer's dense mass and a single fibre. At the yarn's colour as
    // it stands (a dense mass) it washes out over the yarn, at one fibre's colour a layer of them is
    // too dark, so the geometric mean of the two. Matched on the car against the yarn's red.
    // A dense coat (hairKeepsYarnColour) stays close to the yarn: that mean would darken and saturate it.
    const Vec3 c = yarn->spec.color;
    return hair(fibreLook.hairKeepsYarnColour ? hairColourOfCoat(c) : hairColourOfLayer(c, fibreLook.hair.fibreAzimRoughness));
}

std::optional<std::vector<ExportMaterial>> ExportScene::materials(std::string* err) const {
    std::set<std::string> referenced;
    for (const PartFile& e : meshes_) referenced.insert(e.material);
    for (const PartFile& e : strands_) referenced.insert(e.material);
    std::vector<ExportMaterial> out;
    for (const BaseMaterial& m : baseMaterials)
        if (referenced.count(m.id))
            out.push_back(resolvedMaterial(*this, m.id, m.spec, bakedMaps.count(m.id) ? MaterialKind::BakedYarn : kindOfSpec(m.spec)));
    if (plyShading.on)
        for (int k = 0; k < yarnEntries && k < static_cast<int>(baseMaterials.size()); ++k) {
            const std::string id = yarnMaterial(baseMaterials[k].id);
            if (referenced.count(id))
                out.push_back(resolvedMaterial(*this, id, baseMaterials[k].spec, MaterialKind::PlyMappedYarn));
        }
    for (const auto& [id, fibre] : fibreMaterials_) {
        if (!referenced.count(id)) continue;
        const std::optional<MaterialSpec> spec = fibreMaterialSpec(id, err);
        if (!spec) return std::nullopt;
        out.push_back(resolvedMaterial(*this, id, *spec, kindOfSpec(*spec)));
    }
    if (groundPlane) {
        MaterialSpec ground;
        ground.color = {groundAlbedo, groundAlbedo, groundAlbedo};
        ground.sheen = 0;
        out.push_back(resolvedMaterial(*this, kGround, ground, MaterialKind::Surface));
    }
    return out;
}

std::optional<SceneFrame> frameScene(const ExportScene& scene, std::string* err) {
    const Vec3& lo = scene.bounds().lo;
    const Vec3& hi = scene.bounds().hi;
    if (!(lo.x <= hi.x && lo.y <= hi.y && lo.z <= hi.z) || !finite(lo) || !finite(hi))
        return refuse(err, "scene.json: no finite part bounds to frame the camera");
    SceneFrame frame;
    frame.centre = (lo + hi) * 0.5;
    frame.radius = 0.5 * norm(hi - lo);
    const double dist = 1.15 * frame.radius / std::tan(frame.fov * kPi / 360.0);
    frame.eye = frame.centre + Vec3{0, 0.15 * frame.radius, dist};
    bool camerasFinite = finite(frame.eye) && finite(frame.centre);
    for (const ExportCamera& cam : scene.extraCameras)
        camerasFinite = camerasFinite && finite(cam.eye) && finite(cam.lookAt) && std::isfinite(cam.fov);
    if (!camerasFinite) return refuse(err, "scene.json: a camera is not finite");
    return frame;
}

std::string jsonEscaped(const std::string& text) {
    std::string out;
    for (const char c : text) {
        if (c == '"' || c == '\\') out += {'\\', c};
        else if (c == '\n') out += "\\n";
        else if (static_cast<unsigned char>(c) < 0x20) {
            char b[8];
            std::snprintf(b, sizeof b, "\\u%04x", static_cast<unsigned>(c));
            out += b;
        } else out += c;
    }
    return out;
}

} // namespace knit
