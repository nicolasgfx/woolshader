// The persistent cache of what the Knit job showed, so a known scene opens at once.
#include "jobs/scene_cache.h"

#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>

#include "io/cache_paths.h"
#include "io/files.h"
#include "io/publish.h"
#include "knitting/tube/yarn_surface.h"
#include "model/hash.h"
#include "scene/result_keys.h"

namespace knit {

namespace fs = std::filesystem;

namespace {

constexpr char kMagic[8] = {'K', 'N', 'I', 'T', 'S', 'C', 'N', '1'};
// 2: part colours, base-colour textures and per-vertex colours of the meshes, fibre populations;
// 3: BuiltScene::displayNote; 4: DisplayLines::ownColor;
// 5: without what nothing read: the yarn lines, the lines' default colours, the meshes' colour slots and the
//    flat pattern references.
constexpr uint32_t kFormat = 5;
// Bumped whenever the Knit job's display output changes for the same settings and assets:
// 2: the source surfaces with their colour (the viewer's Mesh view), the knit stored with it;
// 3: the teddy's source surfaces at their true scale; 4: knit.bin v2 (with its checksum);
// 5: the display tube at the strand tessellation (strand.tube_*), not 0.8 mm / 6 sides;
// 6: the teddy's pieces drawn at the drawn diameter D, not the knit diameter k D;
// 7: the teddy's separate ply tubes drawn (strand.plies > 1), not the yarn's single tube;
// 8: the yarn display carries its stitch colours whatever the colour source (stitch_colours.h);
// 9: overlap coverage measures the yarn around separate plies (no patch over a whole component);
// 10: base colours through the texture's KHR_texture_transform (the man);
// 11: gap coverage at a fixed tube, not strand.tube_* (an entry knitted at another tessellation held
//     other gap patches under the same key).
constexpr int kGenerator = 11;
constexpr const char* kFile = "scene.bin";

// The display arrays as they are in memory (little-endian): an i64 count, then the elements. A count the remaining
// bytes cannot hold refuses the read (ByteReader::count), so a corrupt entry never allocates beyond its size.
void writeFloats(ByteWriter& w, const std::vector<float>& v) {
    w.i64(static_cast<int64_t>(v.size()));
    if (!v.empty()) w.bytes(v.data(), v.size() * sizeof(float));
}
std::vector<float> readFloats(ByteReader& r) {
    std::vector<float> v(r.count(sizeof(float)));
    if (!v.empty()) r.bytes(v.data(), v.size() * sizeof(float));
    return v;
}
void writeIndices(ByteWriter& w, const std::vector<uint32_t>& v) {
    w.i64(static_cast<int64_t>(v.size()));
    if (!v.empty()) w.bytes(v.data(), v.size() * sizeof(uint32_t));
}
std::vector<uint32_t> readIndices(ByteReader& r) {
    std::vector<uint32_t> v(r.count(sizeof(uint32_t)));
    if (!v.empty()) r.bytes(v.data(), v.size() * sizeof(uint32_t));
    return v;
}
void writeBytes(ByteWriter& w, const std::vector<uint8_t>& v) {
    w.i64(static_cast<int64_t>(v.size()));
    if (!v.empty()) w.bytes(v.data(), v.size());
}
std::vector<uint8_t> readBytes(ByteReader& r) {
    std::vector<uint8_t> v(r.count(1));
    if (!v.empty()) r.bytes(v.data(), v.size());
    return v;
}
// A colour channel: 4 bytes, the float as it is (never widened to double).
void writeFloat(ByteWriter& w, float f) { w.bytes(&f, sizeof f); }
float readFloat(ByteReader& r) {
    float f = 0;
    r.bytes(&f, sizeof f);
    return f;
}

} // namespace

std::string sceneCacheRoot() { return outputCacheDir("scenes"); }

std::string sceneCacheKey(const SceneSettings& s, uint64_t assetHash) {
    char b[360];
    // The display's tube is the strand's tessellation (a Fibres-stage setting: not in the knit or
    // surface key).
    std::snprintf(b, sizeof b, "scene-cache v%u gen %d scene=%s knit=%016llx surface=%016llx asset=%016llx tube=%d/%a/%d/%a", kFormat,
                  kGenerator, s.scene.c_str(), static_cast<unsigned long long>(stageKey(s, Stage::Knit, assetHash)),
                  static_cast<unsigned long long>(stageKey(s, Stage::Surface, assetHash)), static_cast<unsigned long long>(assetHash),
                  s.tubeSides, s.tubeRing, s.plyTubeSides, s.plyTubeRing);
    return b;
}

std::string serializeBuiltScene(const BuiltScene& b) {
    ByteWriter w;
    w.bytes(kMagic, sizeof kMagic);
    w.u32(kFormat);
    w.u32(static_cast<uint32_t>(b.meshes.size()));
    for (const DisplayMesh& d : b.meshes) {
        w.string64(d.name);
        writeFloats(w, d.pos);
        writeFloats(w, d.nrm);
        writeIndices(w, d.idx);
        writeFloats(w, d.uv);
        for (float c : d.color) writeFloat(w, c);
        w.i32(static_cast<int32_t>(d.kind));
        w.u8(d.hasPartColor);
        for (float c : d.partColor) writeFloat(w, c);
        writeFloats(w, d.vcol);
        writeFloats(w, d.albedoUv);
        w.i32(d.albedo);
    }
    w.u32(static_cast<uint32_t>(b.albedo.size()));
    for (const auto& t : b.albedo) {
        w.string64(t->name);
        w.i32(t->width);
        w.i32(t->height);
        writeBytes(w, t->rgba);
    }
    writeFloats(w, b.fibres.pos);
    writeBytes(w, b.fibres.population);
    writeFloats(w, b.fibres.partColor);
    writeFloats(w, b.fibres.ownColor);
    w.u32(static_cast<uint32_t>(b.parts.size()));
    for (const PartRow& r : b.parts) {
        w.string64(r.name);
        w.u8(r.knitted);
        w.i32(static_cast<int32_t>(r.verdict));
        w.i32(r.stitches);
        w.f64(r.yarnD);
        w.string64(r.note);
        w.u32(static_cast<uint32_t>(r.lines.size()));
        for (const std::string& l : r.lines) w.string64(l);
        w.vec3(r.lo);
        w.vec3(r.hi);
    }
    w.string64(b.report);
    w.string64(b.displayNote);
    w.vec3(b.lo);
    w.vec3(b.hi);
    w.f64(b.yarnLength);
    w.f64(b.ms);
    return w.written();
}

bool deserializeBuiltScene(const std::string& bytes, BuiltScene& b) {
    if (bytes.size() < sizeof kMagic || std::memcmp(bytes.data(), kMagic, sizeof kMagic) != 0) return false;
    ByteReader r(bytes, sizeof kMagic);
    if (r.u32() != kFormat) return false;
    BuiltScene o;
    const uint32_t nm = r.u32();
    for (uint32_t i = 0; r.ok() && i < nm; ++i) {
        DisplayMesh d;
        d.name = r.string64();
        d.pos = readFloats(r);
        d.nrm = readFloats(r);
        d.idx = readIndices(r);
        d.uv = readFloats(r);
        for (float& c : d.color) c = readFloat(r);
        d.kind = static_cast<DisplayMesh::Kind>(r.i32());
        d.hasPartColor = r.u8() != 0;
        for (float& c : d.partColor) c = readFloat(r);
        d.vcol = readFloats(r);
        d.albedoUv = readFloats(r);
        d.albedo = r.i32();
        o.meshes.push_back(std::move(d));
    }
    const uint32_t nt = r.u32();
    for (uint32_t i = 0; r.ok() && i < nt; ++i) {
        auto t = std::make_shared<DisplayTexture>();
        t->name = r.string64();
        t->width = r.i32();
        t->height = r.i32();
        t->rgba = readBytes(r);
        if (t->rgba.size() != size_t(std::max(0, t->width)) * size_t(std::max(0, t->height)) * 4) return false;
        o.albedo.push_back(std::move(t));
    }
    o.fibres.pos = readFloats(r);
    o.fibres.population = readBytes(r);
    o.fibres.partColor = readFloats(r);
    o.fibres.ownColor = readFloats(r);
    const uint32_t np = r.u32();
    for (uint32_t i = 0; r.ok() && i < np; ++i) {
        PartRow p;
        p.name = r.string64();
        p.knitted = r.u8() != 0;
        p.verdict = static_cast<Status>(r.i32());
        p.stitches = r.i32();
        p.yarnD = r.f64();
        p.note = r.string64();
        const uint32_t nl = r.u32();
        for (uint32_t k = 0; r.ok() && k < nl; ++k) p.lines.push_back(r.string64());
        p.lo = r.vec3();
        p.hi = r.vec3();
        o.parts.push_back(std::move(p));
    }
    o.report = r.string64();
    o.displayNote = r.string64();
    o.lo = r.vec3();
    o.hi = r.vec3();
    o.yarnLength = r.f64();
    o.ms = r.f64();
    if (!r.ok() || r.at() != bytes.size()) return false;  // truncated, or bytes left over
    b = std::move(o);
    return true;
}

std::string sceneCacheEntryDir(const std::string& root, const SceneSettings& s, uint64_t assetHash) {
    const std::string key = sceneCacheKey(s, assetHash);
    char id[24];
    std::snprintf(id, sizeof id, "%016llx", static_cast<unsigned long long>(fnv1a(key.data(), key.size())));
    return root + "/" + s.scene + "-" + id;
}

bool loadSceneCacheKnit(const std::string& root, const SceneSettings& s, uint64_t assetHash, std::string& bytes) {
    const std::string target = sceneCacheEntryDir(root, s, assetHash);
    if (!isPublished(target, sceneCacheKey(s, assetHash))) return false;
    std::ifstream f(target + "/knit.bin", std::ios::binary);
    if (!f) return false;
    bytes.assign((std::istreambuf_iterator<char>(f)), {});
    return true;
}

bool storeSceneCache(const std::string& root, const SceneSettings& s, uint64_t assetHash, const BuiltScene& b, std::string* err,
                     const std::string* knitBytes) {
    const std::string key = sceneCacheKey(s, assetHash);
    const std::string target = sceneCacheEntryDir(root, s, assetHash);
    std::error_code ec;
    fs::create_directories(root, ec);
    recoverPublished(target);
    Publisher pub;
    if (!pub.begin(target, err, PublishPolicy::CacheEntry)) return false;
    const std::string bytes = serializeBuiltScene(b);
    {
        std::ofstream f(pub.staging() + "/" + kFile, std::ios::binary);
        f.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
        if (!f) {
            if (err) *err = "cannot write the scene cache in " + pub.staging();
            return false;
        }
    }
    if (knitBytes) {
        std::ofstream f(pub.staging() + "/knit.bin", std::ios::binary);
        f.write(knitBytes->data(), static_cast<std::streamsize>(knitBytes->size()));
        if (!f) {
            if (err) *err = "cannot write the knit cache in " + pub.staging();
            return false;
        }
    }
    {
        // Validate before publishing: the file reads back.
        std::ifstream f(pub.staging() + "/" + kFile, std::ios::binary);
        const std::string back((std::istreambuf_iterator<char>(f)), {});
        BuiltScene chk;
        if (back != bytes || !deserializeBuiltScene(back, chk)) {
            if (err) *err = "the scene cache does not read back";
            return false;
        }
    }
    return pub.commit(key, err);
}

std::shared_ptr<BuiltScene> loadSceneCache(const std::string& root, const SceneSettings& s, uint64_t assetHash) {
    const std::string key = sceneCacheKey(s, assetHash);
    const std::string target = sceneCacheEntryDir(root, s, assetHash);
    recoverPublished(target);
    if (!isPublished(target, key)) return nullptr;
    std::ifstream f(target + "/" + kFile, std::ios::binary);
    const std::string bytes((std::istreambuf_iterator<char>(f)), {});
    auto b = std::make_shared<BuiltScene>();
    if (!deserializeBuiltScene(bytes, *b)) return nullptr;
    b->cached = true;
    return b;
}

} // namespace knit
