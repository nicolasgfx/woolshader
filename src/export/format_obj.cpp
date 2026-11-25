// Wavefront OBJ: one .obj with every part as a named object, fibres as polylines, and an .mtl with one
// base colour per material. It reads the parts back one at a time: a scene does not fit in memory.
#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>

#include "export/export_scene.h"
#include "io/files.h"
#include "model/cancel.h"
#include "model/coordinates.h"
#include "model/log.h"

namespace knit {

namespace fs = std::filesystem;

namespace {

// Text output in bounded chunks with a cancellation checkpoint per chunk.
class ChunkedWriter {
public:
    explicit ChunkedWriter(const std::string& path) : f_(path, std::ios::binary) {}
    bool ok() const { return static_cast<bool>(f_); }
    void put(const char* s, size_t n) {
        buf_.append(s, n);
        if (buf_.size() >= (8u << 20)) flush();
    }
    void put(const std::string& s) { put(s.data(), s.size()); }
    bool close() {
        flush();
        f_.close();
        return static_cast<bool>(f_);
    }

private:
    void flush() {
        f_.write(buf_.data(), static_cast<std::streamsize>(buf_.size()));
        buf_.clear();
        cancelCheckpoint();
    }
    std::ofstream f_;
    std::string buf_;
};

// Running OBJ numbers. Every vertex writes v and vn, only meshes with uv write vt, so a mix of
// meshes with and without uv needs its own vt numbering.
struct ObjNumbers {
    size_t vertex = 1, texture = 1, triangles = 0, lines = 0;
};

std::string objName(const std::string& s) {
    std::string o;
    for (char c : s) o += (c == ' ' || c == '\t') ? '_' : c;
    return o;
}

// Package (Y-up metres) -> library (Z-up mm) -> the imported file's coordinates. Without an import
// the numbers pass untouched, so a package's file and its OBJ hold the same values.
Vec3 outputPoint(const ExportScene& scene, const Vec3& p) {
    if (!scene.importFrame) return p;
    const ImportFrame& frame = *scene.importFrame;
    return toFileCoordinates(frame.scale, frame.offset, frame.up, Vec3{p.x, -p.z, p.y} * 1000.0);
}

Vec3 outputNormal(const ExportScene& scene, const Vec3& n) {
    if (!scene.importFrame) return n;
    const Vec3 z{n.x, -n.z, n.y};
    return normalized(scene.importFrame->up == UpAxis::Y ? Vec3{z.x, z.z, -z.y} : z);
}

// One output unit per metre, for the fibres' radii.
double outputUnit(const ExportScene& scene) {
    return scene.importFrame ? norm(outputPoint(scene, {1, 0, 0}) - outputPoint(scene, {0, 0, 0})) : 1.0;
}

std::string coordinateNote(const ExportScene& scene) {
    if (!scene.importFrame) return "metres, Y up (as the photon_path_tracer package)";
    const std::string up = scene.importFrame->up == UpAxis::Y ? "Y" : "Z";
    return "the imported file's units and up axis (" + up + " up)";
}

// The model's meshes. The name decides, because a source mesh with no name and a material called
// "ground" takes the ground quad's name, and neither is written.
std::vector<PartFile> modelMeshes(const ExportScene& scene) {
    std::vector<PartFile> meshes;
    for (const PartFile& part : scene.meshParts())
        if (part.name != ExportScene::kGround) meshes.push_back(part);
    return meshes;
}

// One base colour per material the parts use, sorted by id; a duplicate id takes the last entry, an
// unknown id a light grey. Sheen, coat and fibre shading are reduced to the base colour. A base
// colour texture goes beside the OBJ (the package's uvs are its vt, v up).
bool writeMtlFile(const ExportScene& scene, const std::vector<PartFile>& meshes, const FormatTarget& target,
                  std::string* err) {
    const std::optional<std::vector<ExportMaterial>> materials = scene.materials(err);
    if (!materials) return false;
    std::map<std::string, ExportMaterial> byId;
    for (const ExportMaterial& m : *materials) byId[m.id] = m;
    std::set<std::string> used;
    for (const PartFile& part : meshes) used.insert(part.material);
    for (const PartFile& part : scene.strandParts()) used.insert(part.material);
    const std::string path = target.folder + "/" + target.baseName + ".mtl";
    std::ofstream m(path);
    m << "# wool_shader knit export: base colours of the photon_path_tracer package's materials\n";
    for (const std::string& id : used) {
        const auto it = byId.find(id);
        const Vec3 k = it != byId.end() ? colourAsWritten(it->second) : Vec3{0.8, 0.8, 0.8};
        char b[256];
        std::snprintf(b, sizeof b, "\nnewmtl %s\nKd %.4f %.4f %.4f\nKa 0 0 0\nKs 0 0 0\nillum 1\n", objName(id).c_str(), k.x, k.y, k.z);
        m << b;
        if (it == byId.end()) continue;
        // A baked yarn: its cavity map as the colour map, its normal map as the bump map; the others are not OBJ's.
        const auto baked = scene.bakedMaps.find(id);
        const std::string colourMap = baked != scene.bakedMaps.end() ? baked->second.cavity : it->second.colourTexture;
        const std::string bumpMap = baked != scene.bakedMaps.end() ? baked->second.normal : "";
        const std::pair<const char*, const std::string*> maps[2] = {{"map_Kd ", &colourMap}, {"map_bump ", &bumpMap}};
        for (const auto& [statement, texture] : maps) {
            if (texture->empty()) continue;
            const std::string copy = target.folder + "/" + *texture;
            std::error_code ec;
            fs::create_directories(fs::path(copy).parent_path(), ec);
            if (!fs::copy_file(scene.folder() + "/" + *texture, copy, fs::copy_options::overwrite_existing, ec) && ec)
                return fail(err, "cannot copy " + *texture);
            m << statement << *texture << "\n";
        }
    }
    return m ? true : fail(err, "cannot write " + path);
}

bool writeMeshObject(ChunkedWriter& w, const ExportScene& scene, const PartFile& part, ObjNumbers& at,
                     std::string* err) {
    const std::optional<StoredMesh> mesh = readPartMesh(scene.pathOf(part), err);
    if (!mesh) return false;
    w.put("o " + objName(part.name) + "\nusemtl " + objName(part.material) + "\n");
    const bool hasUv = mesh->stride == 8;
    const size_t n = mesh->vertices.size() / mesh->stride;
    char b[192];
    for (size_t i = 0; i < n; ++i) {
        const float* q = &mesh->vertices[i * mesh->stride];
        const Vec3 x = outputPoint(scene, {q[0], q[1], q[2]}), nn = outputNormal(scene, {q[3], q[4], q[5]});
        int k = std::snprintf(b, sizeof b, "v %.7g %.7g %.7g\nvn %.5f %.5f %.5f\n", x.x, x.y, x.z, nn.x, nn.y, nn.z);
        w.put(b, static_cast<size_t>(k));
        if (hasUv) {
            k = std::snprintf(b, sizeof b, "vt %.6g %.6g\n", q[6], q[7]);
            w.put(b, static_cast<size_t>(k));
        }
    }
    const std::vector<uint32_t>& tri = mesh->tri;
    for (size_t t = 0; t + 2 < tri.size(); t += 3) {
        const size_t v0 = at.vertex + tri[t], v1 = at.vertex + tri[t + 1], v2 = at.vertex + tri[t + 2];
        const size_t t0 = at.texture + tri[t], t1 = at.texture + tri[t + 1], t2 = at.texture + tri[t + 2];
        const int k = hasUv ? std::snprintf(b, sizeof b, "f %zu/%zu/%zu %zu/%zu/%zu %zu/%zu/%zu\n", v0, t0, v0, v1, t1, v1, v2, t2, v2)
                            : std::snprintf(b, sizeof b, "f %zu//%zu %zu//%zu %zu//%zu\n", v0, v0, v1, v1, v2, v2);
        w.put(b, static_cast<size_t>(k));
    }
    at.vertex += n;
    if (hasUv) at.texture += n;
    at.triangles += tri.size() / 3;
    return true;
}

// Each strand one line element through its vertices. OBJ lines carry no radius; the comment gives the
// group's radii in the output's units.
bool writeFibreObject(ChunkedWriter& w, const ExportScene& scene, const PartFile& part, ObjNumbers& at,
                      std::string* err) {
    StrandSetF s;
    if (!readStrandsFile(scene.pathOf(part), s, err)) return false;
    const double unit = outputUnit(scene);
    double rMin = 1e300, rMax = 0, rSum = 0;
    for (const float r : s.radius) {
        rMin = std::min(rMin, double(r));
        rMax = std::max(rMax, double(r));
        rSum += r;
    }
    const size_t n = s.vertices();
    w.put("o " + objName(part.name) + "\nusemtl " + objName(part.material) + "\n");
    char b[192];
    std::snprintf(b, sizeof b, "# %zu fibres, radius %.4g (min %.4g, max %.4g)\n", s.strands(), n ? rSum / double(n) * unit : 0.0,
                  n ? rMin * unit : 0.0, rMax * unit);
    w.put(b);
    for (size_t i = 0; i < n; ++i) {
        const Vec3 x = outputPoint(scene, {s.xyz[3 * i], s.xyz[3 * i + 1], s.xyz[3 * i + 2]});
        const int k = std::snprintf(b, sizeof b, "v %.7g %.7g %.7g\n", x.x, x.y, x.z);
        w.put(b, static_cast<size_t>(k));
    }
    std::string line;
    for (size_t f = 0; f < s.strands(); ++f) {
        line = "l";
        for (uint32_t v = s.first[f]; v < s.first[f + 1]; ++v) line += " " + std::to_string(at.vertex + v);
        line += "\n";
        w.put(line);
    }
    at.vertex += n;
    at.lines += s.strands();
    return true;
}

} // namespace

bool writeObjFile(const ExportScene& scene, const FormatTarget& target, std::string& log, std::string* err) {
    const std::vector<PartFile> meshes = modelMeshes(scene);
    if (meshes.empty()) return fail(err, "the package lists no meshes");
    if (!writeMtlFile(scene, meshes, target, err)) return false;
    const std::string base = target.baseName, obj = target.folder + "/" + base + ".obj";
    ChunkedWriter w(obj);
    if (!w.ok()) return fail(err, "cannot write " + obj);
    w.put("# wool_shader knit export (Wavefront OBJ)\n# coordinates: " + coordinateNote(scene) +
          "\n# fibres and hair: polylines (l), one object per fibre group; their radius is in each object's comment\n"
          "# not included: the fibres' per-vertex radius, ply shading maps (see the photon_path_tracer package)\nmtllib " + base + ".mtl\n");
    ObjNumbers at;
    for (const PartFile& part : meshes)
        if (!writeMeshObject(w, scene, part, at, err)) return false;
    for (const PartFile& part : scene.strandParts())
        if (!writeFibreObject(w, scene, part, at, err)) return false;
    if (!w.close()) return fail(err, "cannot write " + obj);
    appendFormat(log, "OBJ: %zu objects, %zu triangles, %zu fibre lines, %zu vertices -> %s.obj\n",
                 meshes.size() + scene.strandParts().size(), at.triangles, at.lines, at.vertex - 1, base.c_str());
    return true;
}

} // namespace knit
