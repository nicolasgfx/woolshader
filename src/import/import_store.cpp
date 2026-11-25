// The prepared import as a cache entry: mesh.kmesh, plain.kmesh, textures, and the way back as a scene (import_store.h).
#include "import/import_store.h"
#include "io/cache_paths.h"
#include "model/colour.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>

#include "io/image_io.h"
#include "model/cancel.h"
#include "model/coordinates.h"
#include "model/mesh_check.h"
#include "io/publish.h"
#include "io/files.h"

namespace knit {

namespace fs = std::filesystem;

namespace {

constexpr char kMagic[8] = {'K', 'M', 'E', 'S', 'H', '1', 0, 0};
constexpr uint32_t kKmeshVersion = 1;
constexpr const char* kFile = "mesh.kmesh";
constexpr const char* kPlainFile = "plain.kmesh";

std::string safeFileName(const std::string& name, int index) {
    std::string s;
    for (char c : name) s += (std::isalnum(static_cast<unsigned char>(c)) || c == '.' || c == '-' || c == '_') ? c : '_';
    if (s.empty()) s = "texture";
    return std::to_string(index) + "_" + s;
}

bool writeKmesh(const std::string& path, const PreparedMesh& p, std::string* err) {
    ByteWriter w;
    w.bytes(kMagic, 8);
    w.u32(kKmeshVersion);
    const TriMesh& m = p.mesh;
    w.u64(m.pos.size());
    for (const Vec3& v : m.pos) w.vec3(v);
    w.u64(m.tri.size());
    for (const auto& t : m.tri) w.i32(t[0]), w.i32(t[1]), w.i32(t[2]);
    w.u64(m.uv.size());
    for (const auto& c : m.uv)
        for (const Vec2& v : c) w.vec2(v);
    w.u64(p.colourUv.size());
    for (const auto& c : p.colourUv)
        for (const Vec2& v : c) w.vec2(v);
    w.u64(p.vcolor.size());
    for (const Vec3& v : p.vcolor) w.vec3(v);
    w.u64(p.triMaterial.size());
    for (int i : p.triMaterial) w.i32(i);
    w.u64(p.materials.size());
    for (size_t i = 0; i < p.materials.size(); ++i) {
        const ImportMaterial& mat = p.materials[i];
        w.string32(mat.name);
        w.vec3(mat.kd);
        w.u8(uint8_t(mat.hasKd));
        w.string32(mat.textureRef);
        w.string32(i < p.textures.size() ? p.textures[i] : std::string());  // file name inside the entry
        w.u8(uint8_t(mat.textureMissing));
    }
    w.f64(p.scale);
    w.vec3(p.offset);
    w.u8(uint8_t(p.up == UpAxis::Y ? 0 : 1));
    w.u8(uint8_t(p.chartsGenerated));
    w.u64(p.changes.size());
    for (const std::string& c : p.changes) w.string32(c);
    return replaceFile(path, w.written(), err);
}

bool readKmesh(const std::string& path, PreparedMesh& p, std::string* err) {
    std::string b;
    if (!readWholeFile(path, b, err)) return false;
    if (b.size() < 12 || std::memcmp(b.data(), kMagic, 8) != 0) return fail(err, path + ": not a prepared mesh");
    ByteReader r(b, 8);
    if (r.u32() != kKmeshVersion) return fail(err, path + ": another version");
    p = PreparedMesh{};
    TriMesh& m = p.mesh;
    m.pos.resize(r.count(24));
    for (Vec3& v : m.pos) v = r.vec3();
    m.tri.resize(r.count(12));
    for (auto& t : m.tri) t = {r.i32(), r.i32(), r.i32()};
    m.uv.resize(r.count(48));
    for (auto& c : m.uv)
        for (Vec2& v : c) v = r.vec2();
    p.colourUv.resize(r.count(48));
    for (auto& c : p.colourUv)
        for (Vec2& v : c) v = r.vec2();
    p.vcolor.resize(r.count(24));
    for (Vec3& v : p.vcolor) v = r.vec3();
    p.triMaterial.resize(r.count(4));
    for (int& i : p.triMaterial) i = r.i32();
    const uint64_t nm = r.count(40);
    p.materials.resize(nm);
    p.textures.resize(nm);
    for (size_t i = 0; i < nm; ++i) {
        ImportMaterial& mat = p.materials[i];
        mat.name = r.string32();
        mat.kd = r.vec3();
        mat.hasKd = r.u8() != 0;
        mat.textureRef = r.string32();
        p.textures[i] = r.string32();
        mat.textureMissing = r.u8() != 0;
    }
    p.scale = r.f64();
    p.offset = r.vec3();
    p.up = r.u8() == 0 ? UpAxis::Y : UpAxis::Z;
    p.chartsGenerated = r.u8() != 0;
    p.changes.resize(r.count(4));
    for (std::string& c : p.changes) c = r.string32();
    if (!r.ok()) return fail(err, path + ": truncated or corrupt");
    // Consistency: indices in range, per-triangle arrays of the triangle count.
    for (const auto& t : m.tri)
        for (int v : t)
            if (v < 0 || v >= static_cast<int>(m.pos.size())) return fail(err, path + ": corrupt (index out of range)");
    if ((!m.uv.empty() && m.uv.size() != m.tri.size()) || (!p.colourUv.empty() && p.colourUv.size() != m.tri.size()) ||
        (!p.vcolor.empty() && p.vcolor.size() != m.pos.size()) || p.triMaterial.size() != m.tri.size())
        return fail(err, path + ": corrupt (array sizes)");
    m.computeNormals();
    return true;
}

// Copies of the mesh's textures into the entry's staging folder, each once, and writes the mesh.
bool writeMeshWithTextures(const PreparedMesh& p, const std::string& staging, const std::string& file, std::map<std::string, std::string>& copied, std::string* err) {
    PreparedMesh q = p;
    q.textures.assign(q.materials.size(), "");
    std::error_code ec;
    for (size_t i = 0; i < q.materials.size(); ++i) {
        const std::string& src = q.materials[i].texture;
        if (src.empty()) continue;
        cancelCheckpoint();
        auto it = copied.find(src);
        if (it == copied.end()) {
            const std::string name = "textures/" + safeFileName(fs::path(src).filename().string(), static_cast<int>(copied.size()));
            fs::create_directories(fs::path(staging) / "textures", ec);
            fs::copy_file(src, fs::path(staging) / name, fs::copy_options::overwrite_existing, ec);
            if (ec) return fail(err, "cannot copy " + src + " into the import cache: " + ec.message());
            it = copied.emplace(src, name).first;
        }
        q.textures[i] = it->second;
    }
    return writeKmesh(staging + "/" + file, q, err);
}

} // namespace

bool storeImport(const std::string& root, const ImportRequest& r, const PreparedImport& p, std::string* err) {
    const ImportedMesh& source = *r.file;
    const uint64_t identity = r.identity;
    const std::string target = importCacheDir(root, identity);
    std::error_code ec;
    fs::create_directories(root, ec);
    recoverPublished(target);
    Publisher pub;
    if (!pub.begin(target, err, PublishPolicy::CacheEntry)) return false;
    std::map<std::string, std::string> copied;
    if (!writeMeshWithTextures(p.knitted, pub.staging(), kFile, copied, err)) return false;
    const bool hasPlain = !p.plain.mesh.tri.empty();
    if (hasPlain && !writeMeshWithTextures(p.plain, pub.staging(), kPlainFile, copied, err)) return false;
    {
        std::ofstream f(pub.staging() + "/import.txt", std::ios::binary);
        f << "source = " << source.path << "\nformat = " << source.format << "\nidentity = " << identityHex(identity) << "\n" << importOptionsText(r.options);
        for (const std::string& d : source.dependencies) f << "dependency = " << d << "\n";
        for (const std::string& c : p.knitted.changes) f << "change = " << c << "\n";
        if (!f) return fail(err, "cannot write import.txt");
    }
    // Validate before publishing: the prepared meshes read back.
    PreparedMesh back;
    if (!readKmesh(pub.staging() + "/" + kFile, back, err)) return false;
    if (hasPlain && !readKmesh(pub.staging() + "/" + kPlainFile, back, err)) return false;
    return pub.commit("import " + identityHex(identity), err);
}

bool importStored(const std::string& root, uint64_t identity) {
    const std::string dir = importCacheDir(root, identity);
    recoverPublished(dir);
    return isPublished(dir, "import " + identityHex(identity));
}

std::vector<std::string> storedImportChanges(const std::string& root, uint64_t identity) {
    std::vector<std::string> changes;
    std::ifstream f(importCacheDir(root, identity) + "/import.txt");
    std::string line;
    while (std::getline(f, line))
        if (line.rfind("change = ", 0) == 0) changes.push_back(line.substr(9));
    return changes;
}

namespace {

// One stored mesh; texture paths in its materials point into the entry's directory.
bool loadStoredMesh(const std::string& dir, const std::string& file, PreparedMesh& p, std::string* err) {
    if (!readKmesh(dir + "/" + file, p, err)) return false;
    for (size_t i = 0; i < p.materials.size(); ++i) p.materials[i].texture = i < p.textures.size() && !p.textures[i].empty() ? dir + "/" + p.textures[i] : "";
    return true;
}

} // namespace

bool loadImport(const std::string& dir, PreparedMesh& p, std::string* sourcePath, std::string* err) {
    if (!loadStoredMesh(dir, kFile, p, err)) return false;
    if (sourcePath) {
        std::ifstream f(dir + "/import.txt");
        std::string line;
        while (std::getline(f, line))
            if (line.rfind("source = ", 0) == 0) *sourcePath = line.substr(9);
    }
    return true;
}

namespace {

// One material region of a prepared mesh as a primitive: welded, checked, with its base colour.
SourcePrimitive primitiveOfMaterial(const PreparedMesh& p, int material, const std::vector<int>& tris, const std::vector<Vec3>& triColour,
                                 const std::string& name, const std::map<std::string, std::shared_ptr<const ColourImage>>& textures, const SourceOptions& o) {
    const bool named = material >= 0 && material < static_cast<int>(p.materials.size());
    SourcePrimitive g;
    g.node = g.mesh = name;
    g.nodeIndex = 0;  // the file is one object
    g.material = named && !p.materials[material].name.empty() ? p.materials[material].name : std::string("default");
    g.primitive = material;
    g.cls = SourcePrimitive::Class::Knit;
    std::vector<int> map(p.mesh.pos.size(), -1);
    double area = 0;
    Vec3 col{0, 0, 0};
    for (int t : tris) {
        std::array<int, 3> tr;
        for (int k = 0; k < 3; ++k) {
            int& v = map[p.mesh.tri[t][k]];
            if (v < 0) {
                v = static_cast<int>(g.welded.pos.size());
                g.welded.pos.push_back(p.mesh.pos[p.mesh.tri[t][k]]);
            }
            tr[k] = v;
        }
        g.welded.tri.push_back(tr);
        if (!p.mesh.uv.empty()) g.welded.uv.push_back(p.mesh.uv[t]);
        if (!p.colourUv.empty()) g.colourUv.push_back(p.colourUv[t]);
        g.triColor.push_back(triColour[t]);
        const double a = 0.5 * norm(cross(p.mesh.pos[p.mesh.tri[t][1]] - p.mesh.pos[p.mesh.tri[t][0]], p.mesh.pos[p.mesh.tri[t][2]] - p.mesh.pos[p.mesh.tri[t][0]]));
        area += a;
        col = col + triColour[t] * a;
    }
    if (!p.vcolor.empty()) {
        g.vertexColor.resize(g.welded.pos.size());
        for (size_t v = 0; v < map.size(); ++v)
            if (map[v] >= 0) g.vertexColor[map[v]] = p.vcolor[v];
    }
    g.baseColor = area > 0 ? col * (1.0 / area) : kDefaultYarnColour;
    if (named && !p.materials[material].texture.empty()) {
        g.baseColorImage = p.materials[material].texture;
        const auto it = textures.find(p.materials[material].texture);
        if (it != textures.end()) g.albedo = it->second;
    }
    g.welded.computeNormals();
    g.check = checkMesh(g.welded, o.selfIntersection ? SelfIntersectionTest::Run : SelfIntersectionTest::Skip);
    return g;
}

// One primitive per material used of a stored mesh (in material order; triangles without a material last).
void addPrimitives(const PreparedMesh& p, SourcePrimitive::Class cls, const std::string& name, const SourceOptions& o, SourceScene& out) {
    ImportedMesh materials;
    materials.materials = p.materials;
    std::vector<std::string> warnings;
    const std::vector<ColourImage> textures = loadImportTextures(materials, ImportLimits{}, &warnings);
    for (const std::string& w : warnings) out.log.push_back(w);
    const std::vector<Vec3> triColour = triangleColours(p, textures, kDefaultYarnColour);
    std::map<std::string, std::shared_ptr<const ColourImage>> texByPath;
    for (const ColourImage& t : textures) texByPath[t.path] = std::make_shared<const ColourImage>(t);
    std::map<int, std::vector<int>> byMaterial;
    for (size_t t = 0; t < p.mesh.tri.size(); ++t) byMaterial[p.triMaterial[t]].push_back(static_cast<int>(t));
    for (const auto& [material, tris] : byMaterial) {
        cancelCheckpoint();
        out.prims.push_back(primitiveOfMaterial(p, material, tris, triColour, name, texByPath, o));
        out.prims.back().cls = cls;
    }
}

} // namespace

bool loadImportedScene(const std::string& dir, SourceScene& out, std::string* err, const SourceOptions& o) {
    PreparedMesh p;
    std::string source;
    if (!loadImport(dir, p, &source, err)) return false;
    out = SourceScene{};
    out.path = dir + "/" + kFile;
    bool ok = false;
    out.fileHash = fnv1aFile(out.path, &ok);
    out.scale = p.scale;
    out.offset = p.offset;
    out.imported = true;
    out.importYUp = p.up == UpAxis::Y;
    if (p.mesh.uv.empty() && o.requireCharts)
        return fail(err, "the import has no knitting charts: change the import, and it rebuilds the surface");
    const std::string stem = fs::path(source).stem().string().empty() ? std::string("mesh") : fs::path(source).stem().string();
    addPrimitives(p, SourcePrimitive::Class::Knit, stem, o, out);
    std::error_code ec;
    if (fs::exists(dir + "/" + kPlainFile, ec)) {
        PreparedMesh plain;
        if (!loadStoredMesh(dir, kPlainFile, plain, err)) return false;
        SourceOptions asPlain = o;
        asPlain.selfIntersection = false;  // only a knitted surface needs its inside test
        addPrimitives(plain, SourcePrimitive::Class::Plain, stem + "_plain", asPlain, out);
    }
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const SourcePrimitive& g : out.prims)
        for (const Vec3& q : g.welded.pos) {
            lo = {std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z)};
            hi = {std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z)};
        }
    out.lo = lo;
    out.hi = hi;
    out.log.push_back("import " + source + ": " + std::to_string(out.prims.size()) + " material regions, " + std::to_string(p.mesh.tri.size()) +
                      " triangles" + (p.chartsGenerated ? ", generated knitting charts" : ""));
    return true;
}

} // namespace knit
