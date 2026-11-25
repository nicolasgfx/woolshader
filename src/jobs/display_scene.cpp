// Builds a BuiltScene from the results: tubes, source meshes and fibre lines as CPU arrays.
#include "jobs/display_scene.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <map>
#include <set>

#include "fibres/fibre_budget.h"
#include "fibres/yarn_bake.h"
#include "io/files.h"
#include "io/image_io.h"
#include "knitting/knit/mesh_knit.h"
#include "model/cancel.h"
#include "model/choices.h"
#include "model/mesh.h"
#include "scene/appearance/ply_maps.h"
#include "scene/appearance/yarn_tube.h"
#include "scene/knit/knit_report.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"
#include "settings/yarn_look.h"

namespace knit {

// A display copy of a base-colour texture, at most 2048 texels on a side (box filter).
std::shared_ptr<const DisplayTexture> makeDisplayTexture(const ColourImage& t) {
    auto d = std::make_shared<DisplayTexture>();
    d->name = t.path;
    int f = 1;
    while (t.width / f > 2048 || t.height / f > 2048) f *= 2;
    d->width = std::max(1, t.width / f);
    d->height = std::max(1, t.height / f);
    d->rgba.resize(size_t(d->width) * d->height * 4);
    for (int y = 0; y < d->height; ++y)
        for (int x = 0; x < d->width; ++x)
            for (int c = 0; c < 4; ++c) {
                unsigned sum = 0;
                for (int j = 0; j < f; ++j)
                    for (int i = 0; i < f; ++i) sum += t.rgba[(size_t(y * f + j) * t.width + size_t(x * f + i)) * 4 + c];
                d->rgba[(size_t(y) * d->width + x) * 4 + c] = uint8_t(sum / unsigned(f * f));
            }
    return d;
}

void addSourceMesh(BuiltScene& out, std::map<const ColourImage*, int>& textureSlots, const std::string& name, const SourcePrimitive& g,
                   const std::vector<int>* triangles, const SourceLook& look) {
    const Vec3& tint = look.tint;
    DisplayMesh d;
    d.name = name;
    d.kind = look.kind;
    TriMesh m = g.welded;
    if (m.vnormal.size() != m.pos.size()) m.computeNormals();
    const bool tex = g.albedo != nullptr && (!g.colourUv.empty() || !m.uv.empty());
    if (tex) {
        auto it = textureSlots.find(g.albedo.get());
        if (it == textureSlots.end()) {
            it = textureSlots.emplace(g.albedo.get(), static_cast<int>(out.albedo.size())).first;
            out.albedo.push_back(makeDisplayTexture(*g.albedo));
        }
        d.albedo = it->second;
    }
    auto corner = [&](int t, int k) {
        const int v = m.tri[t][k];
        const Vec3& n = m.vnormal[v];
        const Vec3 p = m.pos[v] - n * look.inset;
        d.pos.insert(d.pos.end(), {float(p.x), float(p.y), float(p.z)});
        d.nrm.insert(d.nrm.end(), {float(n.x), float(n.y), float(n.z)});
        if (tex) {
            // Colour uv in the OBJ convention (v up); the display texture's rows are top first.
            const Vec2 uv = !g.colourUv.empty() ? g.colourUv[t][k] : m.uv[t][k];
            d.albedoUv.insert(d.albedoUv.end(), {float(uv.x), float(1.0 - uv.y)});
        }
        if (!g.vertexColor.empty()) {
            const Vec3& c = g.vertexColor[v];
            d.vcol.insert(d.vcol.end(), {float(c.x * tint.x), float(c.y * tint.y), float(c.z * tint.z)});
        }
        d.idx.push_back(static_cast<uint32_t>(d.idx.size()));
    };
    if (triangles)
        for (int t : *triangles)
            for (int k = 0; k < 3; ++k) corner(t, k);
    else
        for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
            for (int k = 0; k < 3; ++k) corner(t, k);
    const Vec3 c{g.baseColor.x * tint.x, g.baseColor.y * tint.y, g.baseColor.z * tint.z};
    d.color[0] = float(tex ? tint.x : c.x), d.color[1] = float(tex ? tint.y : c.y), d.color[2] = float(tex ? tint.z : c.z);
    out.meshes.push_back(std::move(d));
}

KnitStandIns knitStandIns(const Scene& scene, double inset) {
    BuiltScene b;
    std::map<const ColourImage*, int> textureSlots;  // this call's own: indices into b.albedo
    for (const SourcePrimitive& g : scene.source.prims)
        if (g.cls == SourcePrimitive::Class::Plain) addSourceMesh(b, textureSlots, g.mesh, g, nullptr, kPlainLook);
    const long fixed = static_cast<long>(b.meshes.size());
    for (const ScenePart& r : scene.parts) {
        const SourcePrimitive& g = scene.source.prims[r.prim];
        addSourceMesh(b, textureSlots, r.name, g, &g.check.components[r.component].triangles, standInLook(inset));
    }
    KnitStandIns out;
    out.fixed.assign(b.meshes.begin(), b.meshes.begin() + fixed);
    out.base.assign(b.meshes.begin() + fixed, b.meshes.end());
    out.albedo = std::move(b.albedo);
    return out;
}

// --- Yarn tubes

// The display tube of a part's strands with `sides` sides and `ring` mm between rings (the strand tessellation, or
// the display's coarser one); with uvRepeat > 0 the export's uv convention (format_ppt.cpp: u = s / pitch, wrapped).
static TubeMesh displayTube(const std::string& name, const RenderedPart& rp, double ring, int sides, double uvRepeat) {
    TubeStyle style;
    style.modelRing = ring;
    style.ring = ring;
    style.sides = sides;
    if (uvRepeat > 0) {
        style.uvRepeat = uvRepeat;
        style.uvWrap = 64;
    }
    TubeMesh tm;
    std::string err;
    if (!buildYarnTube(*rp.strands, rp.strandPp, style, tm, &err)) throw std::runtime_error(name + ": display tube on the GPU: " + err);
    return tm;
}

// A tube as a yarn display mesh at `offset`, in the yarn colour.
static DisplayMesh tubeDisplayMesh(const std::string& name, const TubeMesh& tm, const Vec3& offset, bool withUv, const Vec3& color) {
    DisplayMesh d;
    d.name = name;
    d.kind = DisplayMesh::Kind::Yarn;
    for (size_t v = 0; v < tm.pos.size(); ++v) {
        const Vec3 p = tm.pos[v] + offset;
        d.pos.insert(d.pos.end(), {float(p.x), float(p.y), float(p.z)});
        d.nrm.insert(d.nrm.end(), {float(tm.nrm[v].x), float(tm.nrm[v].y), float(tm.nrm[v].z)});
    }
    d.idx = tm.tri;
    if (withUv) d.uv = tm.uv;
    d.color[0] = float(color.x), d.color[1] = float(color.y), d.color[2] = float(color.z);
    return d;
}

// The display colour (sRGB) of each tube vertex: its stitch's colour, found from its strand and arc length; a yarn
// without a stitch takes its part's. Empty when the part has no stitch colours. The viewer shows these colours while
// the texture source is chosen, so they are there whatever the current source.
static std::vector<float> stitchVertexColours(const TubeMesh& tm, const StitchColours& colours, const RenderedPart& rp, const Vec3& partColour) {
    std::vector<float> out;
    if (colours.key.empty() || rp.yarn->topo.empty() || tm.yarn.size() != tm.pos.size() || tm.s.size() != tm.pos.size()) return out;
    const StitchAlongYarn along(*rp.yarn, rp.plies.get());
    const uint32_t partKey = stitchColourKey(partColour);
    for (size_t v = 0; v < tm.pos.size(); ++v) {
        const int st = along.at(tm.yarn[v], tm.s[v]);
        const uint32_t k = st >= 0 && st < static_cast<int>(colours.key.size()) ? colours.key[st] : partKey;
        const Vec3 c{(k >> 16 & 0xff) / 255.0, (k >> 8 & 0xff) / 255.0, (k & 0xff) / 255.0};
        out.insert(out.end(), {float(c.x), float(c.y), float(c.z)});
    }
    return out;
}

Vec3 objectColourOf(const Scene& t, const ScenePart& part, const SceneSettings& s) {
    SceneSettings obj = s;
    obj.yarnColourSource = "object";
    return partYarnColour(t, part, obj);
}

bool addKnittedPartTube(BuiltScene& out, const Scene& t, const ScenePart& part, const SceneSettings& s, const KnitTubeLook& look,
                        std::string* err) {
    const std::optional<RenderedPart> rendered = renderedPart(part, s, err);
    if (!rendered) return false;
    const RenderedPart& rp = *rendered;
    const Vec3 pc = objectColourOf(t, part, s);
    const StitchColours& colours = sceneStitchColours(t).parts[static_cast<size_t>(&part - t.parts.data())];
    const double pitch = look.texture ? look.texture->pitch : 0.0;
    const double uvRepeat = part.yarnD == look.textureD ? pitch : 0.0;
    for (const YarnPath& y : rp.strands->yarns) out.yarnLength += y.length();
    const TubeMesh tm = displayTube(part.name, rp, look.ring, rp.sides, uvRepeat);
    DisplayMesh d = tubeDisplayMesh(part.name, tm, part.centre, uvRepeat > 0, s.yarnWhite.color);
    d.vcol = stitchVertexColours(tm, colours, rp, pc);
    d.hasPartColor = true;
    d.partColor[0] = float(pc.x), d.partColor[1] = float(pc.y), d.partColor[2] = float(pc.z);
    out.meshes.push_back(std::move(d));
    return true;
}

bool addKnittedPartMappedTubes(BuiltScene& out, const Scene& t, const ScenePart& part, const SceneSettings& s, const KnitTubeLook& look,
                               const RenderedPart& rp, const BakedMapEntry& entry, std::string* err) {
    const Vec3 pc = objectColourOf(t, part, s);
    for (const YarnPath& y : rp.strands->yarns) out.yarnLength += y.length();
    const BakeParams bp = bakeParams(s);
    TubeStyle style;
    style.modelRing = look.ring;
    style.ring = look.ring;
    style.sides = rp.sides;
    style.uvRepeat = bp.texel * bp.atlasWidth;  // one atlas strip of yarn per unit of u
    style.uvWrap = 1;
    TubeMesh tm;
    if (!buildYarnTube(*rp.strands, rp.strandPp, style, tm, err)) return false;
    relabelAtlasUv(tm, bp);
    if (tm.yarnFirstVertex.size() != entry.yarns.size() + 1) {
        if (err) *err = part.name + ": the baked maps hold " + std::to_string(entry.yarns.size()) + " yarns, the tube " +
                        std::to_string(tm.yarnFirstVertex.size() - 1);
        return false;
    }
    for (size_t y = 0; y + 1 < tm.yarnFirstVertex.size(); ++y) {
        const uint32_t first = tm.yarnFirstVertex[y], end = tm.yarnFirstVertex[y + 1];
        char suffix[16];
        std::snprintf(suffix, sizeof suffix, "#y%03zu", y);
        DisplayMesh d;
        d.name = part.name + suffix;
        d.kind = DisplayMesh::Kind::Yarn;
        for (uint32_t v = first; v < end; ++v) {
            const Vec3 p = tm.pos[v] + part.centre;
            d.pos.insert(d.pos.end(), {float(p.x), float(p.y), float(p.z)});
            d.nrm.insert(d.nrm.end(), {float(tm.nrm[v].x), float(tm.nrm[v].y), float(tm.nrm[v].z)});
            d.albedoUv.insert(d.albedoUv.end(), {tm.uv[2 * v], 1.0f - tm.uv[2 * v + 1]});  // image space: the first row on top
        }
        for (size_t k = 0; k + 2 < tm.tri.size(); k += 3) {
            if (tm.tri[k] < first || tm.tri[k] >= end) continue;  // a strip never crosses two yarns
            d.idx.insert(d.idx.end(), {tm.tri[k] - first, tm.tri[k + 1] - first, tm.tri[k + 2] - first});
        }
        d.color[0] = float(s.yarnWhite.color.x), d.color[1] = float(s.yarnWhite.color.y), d.color[2] = float(s.yarnWhite.color.z);
        d.hasPartColor = true;
        d.partColor[0] = float(pc.x), d.partColor[1] = float(pc.y), d.partColor[2] = float(pc.z);
        ColourImage cavity;
        std::vector<uint8_t> rgb;
        if (!bakedCavityMap(entry, static_cast<int>(y), s.yarnBody.mapCavity, cavity.width, cavity.height, rgb, err)) return false;
        cavity.rgba.resize(rgb.size() / 3 * 4);
        for (size_t k = 0; k < rgb.size() / 3; ++k) {
            cavity.rgba[4 * k] = rgb[3 * k], cavity.rgba[4 * k + 1] = rgb[3 * k + 1], cavity.rgba[4 * k + 2] = rgb[3 * k + 2];
            cavity.rgba[4 * k + 3] = 255;
        }
        out.albedo.push_back(makeDisplayTexture(cavity));
        d.albedo = static_cast<int>(out.albedo.size()) - 1;
        out.meshes.push_back(std::move(d));
    }
    return true;
}

bool addKnittedPartTubes(BuiltScene& out, const Scene& t, const ScenePart& part, const SceneSettings& s, const KnitTubeLook& look,
                         std::string* err) {
    if (yarnBodyRenderOf(s) != YarnBodyRender::Maps) return addKnittedPartTube(out, t, part, s, look, err);
    const std::optional<RenderedPart> rendered = renderedPart(part, s, err);
    if (!rendered) return false;
    BakedMapSpec spec;
    if (!bakedMapSpecOf(s, part.name, *rendered->strands, rendered->strandPp, spec, err)) return false;
    BakedMapEntry entry;
    if (loadBakedMapSet(bakeCacheRoot(), spec, entry)) return addKnittedPartMappedTubes(out, t, part, s, look, *rendered, entry, err);
    if (!addKnittedPartTube(out, t, part, s, look, err)) return false;
    out.displayNote = "The yarn maps are not baked yet: the viewport shows the plain tubes (Knit step, Bake the yarn maps).";
    return true;
}

bool carriesYarnMaps(const BuiltScene& b) {
    for (const DisplayMesh& d : b.meshes)
        if (d.kind == DisplayMesh::Kind::Yarn && d.albedo >= 0) return true;
    return false;
}

// The part a yarn mesh belongs to: its name, or the part before "#y" (one mesh per yarn under the maps).
static std::string partOfYarnMesh(const std::string& name) {
    const size_t hash = name.rfind("#y");
    return hash == std::string::npos ? name : name.substr(0, hash);
}

// Puts the rebuilt tubes where the shown parts' tubes were, keeping each part's colour choices: every yarn mesh
// of the shown scene goes (one per part, or one per yarn under the maps), the rebuilt ones take their place,
// their cavity textures appended to the scene's. Every shown part must get tubes again.
bool replaceYarnMeshes(BuiltScene& out, BuiltScene& tubes, std::string* err) {
    struct PartColour {
        bool has = false;
        float rgb[3] = {0, 0, 0};
    };
    std::map<std::string, PartColour> colourOf;
    std::vector<DisplayMesh> kept;
    size_t firstYarn = out.meshes.size();
    for (DisplayMesh& d : out.meshes) {
        if (d.kind != DisplayMesh::Kind::Yarn) {
            kept.push_back(std::move(d));
            continue;
        }
        firstYarn = std::min(firstYarn, kept.size());
        PartColour& c = colourOf[partOfYarnMesh(d.name)];
        c.has = d.hasPartColor;
        std::copy(std::begin(d.partColor), std::end(d.partColor), std::begin(c.rgb));
    }
    const int albedoOffset = static_cast<int>(out.albedo.size());
    std::vector<DisplayMesh> placed;
    std::set<std::string> partsPlaced;
    for (DisplayMesh& n : tubes.meshes) {
        if (n.pos.empty()) continue;
        const auto c = colourOf.find(partOfYarnMesh(n.name));
        if (c == colourOf.end()) continue;
        n.hasPartColor = c->second.has;
        std::copy(std::begin(c->second.rgb), std::end(c->second.rgb), std::begin(n.partColor));
        if (n.albedo >= 0) n.albedo += albedoOffset;
        partsPlaced.insert(partOfYarnMesh(n.name));
        placed.push_back(std::move(n));
    }
    if (partsPlaced.size() != colourOf.size()) {
        if (err) *err = "the rebuilt strand mesh does not match the shown parts (" + std::to_string(partsPlaced.size()) + " of " +
                        std::to_string(colourOf.size()) + "): knit again";
        return false;
    }
    kept.insert(kept.begin() + static_cast<long long>(std::min(firstYarn, kept.size())), std::make_move_iterator(placed.begin()),
                std::make_move_iterator(placed.end()));
    out.meshes = std::move(kept);
    out.albedo.insert(out.albedo.end(), tubes.albedo.begin(), tubes.albedo.end());
    return true;
}

bool showCachedYarnMaps(BuiltScene& c, const Scene& t, const SceneSettings& s, std::string* err) {
    if (yarnBodyRenderOf(s) != YarnBodyRender::Maps || carriesYarnMaps(c)) return true;
    const KnitTubeLook look = knitTubeLook(t, s, mostCommonD(c.parts));
    BuiltScene tubes;
    for (const ScenePart& p : t.parts) {
        cancelCheckpoint();
        if (p.knit && !addKnittedPartTubes(tubes, t, p, s, look, err)) return false;
    }
    if (!replaceYarnMeshes(c, tubes, err)) return false;
    c.displayNote = tubes.displayNote.empty() ? look.displayNote : tubes.displayNote;  // the stored note may say "not baked yet"
    return true;
}

BuiltScene withoutYarnMaps(const BuiltScene& b) {
    BuiltScene plain = b;
    std::vector<bool> used(plain.albedo.size(), false);
    for (DisplayMesh& d : plain.meshes) {
        if (d.kind == DisplayMesh::Kind::Yarn && d.albedo >= 0) {
            d.albedo = -1;
            d.albedoUv.clear();
        }
        if (d.albedo >= 0 && d.albedo < static_cast<int>(used.size())) used[d.albedo] = true;
    }
    for (size_t k = 0; k < plain.albedo.size(); ++k)  // an empty texture keeps the other meshes' indices
        if (!used[k]) plain.albedo[k] = std::make_shared<DisplayTexture>();
    return plain;
}

std::optional<DisplayMesh> coarsePartTube(const ScenePart& part, const SceneSettings& s) {
    std::string why;
    const std::optional<RenderedPart> rendered = renderedPart(part, s, &why);
    if (!rendered) return std::nullopt;
    const RenderedPart& rp = *rendered;
    const TubeMesh tm = displayTube(part.name, rp, std::max(rp.ring, 0.8), std::min(rp.sides, 6), 0.0);
    return tubeDisplayMesh(part.name, tm, part.centre, false, s.yarnWhite.color);
}

// The length [mm] of the tubes a tiled knit shows (plies: one tube each).
static double knitStrandLength(const Scene& t, const SceneSettings& s) {
    double L = 0;
    for (const ScenePart& r : t.parts)
        if (r.knit)
            for (const YarnPath& y : r.knit->t.curves.yarns) L += y.length();
    return L * std::max(1, s.plies);
}

KnitTubeLook knitTubeLook(const Scene& t, const SceneSettings& s, double textureD) {
    KnitTubeLook look;
    look.textureD = textureD;
    if (textureD > 0) look.texture = texturePreview(s, textureD);
    // The strand tessellation (with plies: the ply tube's), or coarser where the scene's tubes would exceed
    // kDisplayTriangles; the export always uses the setting.
    const int sides = s.plies > 1 ? s.plyTubeSides : s.tubeSides;
    const double ring = s.plies > 1 ? s.plyTubeRing : s.tubeRing;
    const double lengthMm = knitStrandLength(t, s);
    const double triangles = 2.0 * sides * lengthMm / ring;
    look.ring = ring;
    if (triangles <= double(kDisplayTriangles)) return look;
    look.ring = 2.0 * sides * lengthMm / double(kDisplayTriangles);
    char b[256];
    std::snprintf(b, sizeof b,
                  "Display: a ring every %.2f mm (%.0f million triangles at %.2f mm exceed the display budget of %lld million); "
                  "the export uses %.2f mm.",
                  look.ring, triangles / 1e6, ring, kDisplayTriangles / 1000000, ring);
    look.displayNote = b;
    return look;
}

// --- Fibre lines

void addFibreLine(DisplayLines& out, std::span<const Vec3f> points, const Vec3& offset, uint8_t population, const Vec3& partColour,
                  const Vec3& ownColour) {
    // Polyline fibres start with the buried root: draw from the emergence point.
    if (points.size() > 2) out.fibreStarts.push_back(static_cast<uint32_t>(out.pos.size() / 3));
    for (size_t k = 2; k < points.size(); ++k) {
        const Vec3 a = toVec3(points[k - 1]) + offset, q = toVec3(points[k]) + offset;
        out.pos.insert(out.pos.end(), {float(a.x), float(a.y), float(a.z), float(q.x), float(q.y), float(q.z)});
        for (int e = 0; e < 2; ++e) {
            out.population.push_back(population);
            out.partColor.insert(out.partColor.end(), {float(partColour.x), float(partColour.y), float(partColour.z)});
            out.ownColor.insert(out.ownColor.end(), {float(ownColour.x), float(ownColour.y), float(ownColour.z)});
        }
    }
}

// --- Fibre display chunks (arrangeFibreLines)

// About 2000 segments per chunk: some 1000 chunks for 2 million segments, few enough to test
// every frame, small enough to cull a useful share in a close-up.
constexpr uint32_t kChunkVertices = 4096;

// The 10-bit cell of x in [lo, hi].
static uint32_t mortonCell(double x, double lo, double hi) {
    const double t = hi > lo ? (x - lo) / (hi - lo) : 0;
    return static_cast<uint32_t>(std::clamp(t, 0.0, 1.0) * 1023.0);
}

// The 10 bits of v spread to every third bit (bit i -> bit 3 i).
static uint32_t spreadBits(uint32_t v) {
    v = (v | (v << 16)) & 0x030000FF;
    v = (v | (v << 8)) & 0x0300F00F;
    v = (v | (v << 4)) & 0x030C30C3;
    v = (v | (v << 2)) & 0x09249249;
    return v;
}

// The 30-bit Morton (Z-order) code of p in the box [lo, hi]: points sorted by it stay near their
// neighbours, so consecutive runs of the sorted fibres are compact chunks.
static uint32_t mortonCode(const Vec3& p, const Vec3& lo, const Vec3& hi) {
    return spreadBits(mortonCell(p.x, lo.x, hi.x)) | spreadBits(mortonCell(p.y, lo.y, hi.y)) << 1 |
           spreadBits(mortonCell(p.z, lo.z, hi.z)) << 2;
}

// A fibre's place in its chunk: a fixed integer hash of its index (lowbias32, C. Wellons), so a
// chunk's fibres are in random order and any prefix of the chunk is an even thinning.
static uint32_t fibreRank(uint32_t i) {
    i ^= i >> 16;
    i *= 0x7feb352dU;
    i ^= i >> 15;
    i *= 0x846ca68bU;
    i ^= i >> 16;
    return i;
}

static Vec3 linePoint(const DisplayLines& l, uint32_t v) { return {l.pos[3 * v], l.pos[3 * v + 1], l.pos[3 * v + 2]}; }

// Appends the vertices [first, end) of `from` to `to` (each per-vertex array that `from` has).
static void appendLineVertices(const DisplayLines& from, uint32_t first, uint32_t end, DisplayLines& to) {
    to.pos.insert(to.pos.end(), from.pos.begin() + 3 * first, from.pos.begin() + 3 * end);
    if (from.partColor.size() == from.pos.size())
        to.partColor.insert(to.partColor.end(), from.partColor.begin() + 3 * first, from.partColor.begin() + 3 * end);
    if (from.ownColor.size() == from.pos.size())
        to.ownColor.insert(to.ownColor.end(), from.ownColor.begin() + 3 * first, from.ownColor.begin() + 3 * end);
    if (from.population.size() * 3 == from.pos.size())
        to.population.insert(to.population.end(), from.population.begin() + first, from.population.begin() + end);
}

// The chunk's bounds.
static void boundChunk(const DisplayLines& lines, FibreChunk& chunk) {
    chunk.lo = {1e300, 1e300, 1e300};
    chunk.hi = {-1e300, -1e300, -1e300};
    for (uint32_t v = chunk.firstVertex; v < chunk.firstVertex + chunk.vertexCount; ++v) {
        const Vec3 p = linePoint(lines, v);
        chunk.lo = {std::min(chunk.lo.x, p.x), std::min(chunk.lo.y, p.y), std::min(chunk.lo.z, p.z)};
        chunk.hi = {std::max(chunk.hi.x, p.x), std::max(chunk.hi.y, p.y), std::max(chunk.hi.z, p.z)};
    }
}

void addToBounds(DisplayBounds& box, const DisplayMesh& d) {
    for (size_t i = 0; i + 2 < d.pos.size(); i += 3) {
        box.lo = {std::min(box.lo.x, double(d.pos[i])), std::min(box.lo.y, double(d.pos[i + 1])), std::min(box.lo.z, double(d.pos[i + 2]))};
        box.hi = {std::max(box.hi.x, double(d.pos[i])), std::max(box.hi.y, double(d.pos[i + 1])), std::max(box.hi.z, double(d.pos[i + 2]))};
    }
}

void updateSceneBounds(BuiltScene& b) {
    DisplayBounds box;
    for (const DisplayMesh& d : b.meshes) addToBounds(box, d);
    if (box.empty()) return;
    b.lo = box.lo;
    b.hi = box.hi;
}

std::shared_ptr<const TexturePreview> texturePreview(const SceneSettings& s, double partD) {
    int texPlies = 1;
    PlyProfileModel model = PlyProfileModel::Swept;
    if (!parseTextureProfile(s.textureProfile, texPlies, model, nullptr) || texPlies <= 1 || s.plies > 1) return nullptr;
    PlyMapSpec spec;
    spec.g.count = texPlies;
    spec.g.pitch = s.plyPitch;
    spec.g.fill = s.plyFill;
    spec.g.handedness = s.plyHandedness;
    spec.pp.model = model;
    spec.yarnD = finalSurface(s, partD).renderedDiameter();
    PlyMapSetResult ms;
    std::string e2;
    auto tp = std::make_shared<TexturePreview>();
    if (!ensurePlyMapSet(plyMapCacheRoot(), spec, ms, &e2) ||
        !readPngRgb8(ms.dir + "/" + ms.set.ao, tp->width, tp->height, tp->rgb)) {
        return tp;  // empty: the map set could not be made
    }
    tp->pitch = spec.g.pitch;
    return tp;
}

float displayThreshold(const SceneFibres& tf) {
    std::vector<std::pair<float, int>> all;
    for (const PartFibres& pf : tf.parts) addKeySegments(pf.fs, all);
    return budgetThreshold(all, kDisplaySegments);
}

double mostCommonD(const std::vector<PartRow>& parts) {
    std::map<double, int> nD;
    for (const PartRow& r : parts)
        if (r.knitted && r.yarnD > 0) ++nD[r.yarnD];
    double mapD = 0;
    for (const auto& [D, n] : nD)
        if (mapD == 0 || n > nD[mapD]) mapD = D;
    return mapD;
}

// The yarn diameter most of the knit's parts have, counted over every part with a knit (before the rows exist).
// Not mostCommonD: that counts rows with stitches and a diameter above 0, and the two may differ, which would move
// which parts carry the texture uvs.
static double mostCommonKnitD(const Scene& t) {
    std::map<double, int> nD;
    for (const ScenePart& r : t.parts)
        if (r.knit) ++nD[r.yarnD];
    double mapD = 0;
    for (const auto& [D, n] : nD)
        if (mapD == 0 || n > nD[mapD]) mapD = D;
    return mapD;
}

bool knitDisplay(BuiltScene& out, const SceneEntry& e, const Scene& t, const SceneSettings& s, std::string* err) {
    out.report = formatKnitReport(t, s, false);
    const KnitTubeLook look = knitTubeLook(t, s, mostCommonKnitD(t));
    if (look.texture) out.texture = look.texture;
    out.displayNote = look.displayNote;
    std::map<const ColourImage*, int> textureSlots;
    for (const SourcePrimitive& g : t.source.prims)
        if (g.cls == SourcePrimitive::Class::Plain) addSourceMesh(out, textureSlots, g.mesh, g, nullptr, kPlainLook);
    // A guided part's generated yarn (GuidedPart::previewTubes) as a plain mesh at the part's place, named for what it is.
    const auto addGeneratedYarn = [&](const ScenePart& r, const SourcePrimitive& g, const std::string& label) {
        SourcePrimitive preview = g;
        preview.welded = *r.guided->previewTubes;
        for (Vec3& p : preview.welded.pos) p = p + r.centre;
        preview.welded.vnormal.clear();
        preview.albedo = nullptr;
        preview.colourUv.clear();
        preview.triColor.clear();
        preview.vertexColor.clear();
        addSourceMesh(out, textureSlots, r.name + " " + label, preview, nullptr, kPlainLook);
    };
    for (const ScenePart& r : t.parts) {
        cancelCheckpoint();
        PartRow row;
        row.name = r.name;
        row.knitted = r.hasStitches();
        row.note = r.note;
        const SourcePrimitive& g = t.source.prims[r.prim];
        const size_t firstMesh = out.meshes.size();
        if (r.knit) {
            row.verdict = r.knit->verdict;
            row.stitches = static_cast<int>(r.knit->t.prog.stitches.size());
            row.yarnD = r.yarnD;
            row.lines.push_back(r.knit->detail);
            for (const CheckLine& l : r.knit->after.lines)
                if (l.status != Status::Pass) row.lines.push_back(std::string(toString(l.status)) + " " + l.name + ": " + l.detail);
            if (r.guided)
                for (const std::string& l : guidedPartLines(*r.guided)) row.lines.push_back(l);
            // A guided knit whose conversion to curves added defects shows the generated yarn (the faithful one); the
            // knit's curves are still what the export writes (GuidedPart::yarnShown, yarnExported say so).
            if (r.guided && r.guided->previewTubes) addGeneratedYarn(r, g, "(generated yarn, before its conversion)");
            else if (!addKnittedPartTubes(out, t, r, s, look, err)) return false;
        } else if (r.guided && r.guided->previewTubes) {
            // A guided part that is not knitted: its generated yarn as a preview, so what was refused or stopped is visible.
            for (const std::string& l : guidedPartLines(*r.guided)) row.lines.push_back(l);
            addGeneratedYarn(r, g, "(generated yarn: preview, not knitted)");
        } else {
            if (r.guided)
                for (const std::string& l : guidedPartLines(*r.guided)) row.lines.push_back(l);
            addSourceMesh(out, textureSlots, r.name, g, &g.check.components[r.component].triangles, kPlainLook);
        }
        DisplayBounds box;
        for (size_t k = firstMesh; k < out.meshes.size(); ++k) addToBounds(box, out.meshes[k]);  // several meshes per part under maps
        row.lo = box.lo;
        row.hi = box.hi;
        out.parts.push_back(std::move(row));
    }
    // The source surfaces of the knitted parts (the viewer's "Mesh" view), with their colour.
    for (const ScenePart& r : t.parts) {
        if (!r.knit) continue;
        const SourcePrimitive& g = t.source.prims[r.prim];
        addSourceMesh(out, textureSlots, r.name + " (source)", g, &g.check.components[r.component].triangles, kSourceLook);
    }
    return true;
}

void arrangeFibreLines(DisplayLines& lines) {
    const size_t fibres = lines.fibreStarts.size();
    const uint32_t vertices = static_cast<uint32_t>(lines.pos.size() / 3);
    auto fibreEnd = [&](size_t i) { return i + 1 < fibres ? lines.fibreStarts[i + 1] : vertices; };
    // 1. The fibres in Z order of their first point.
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (uint32_t v = 0; v < vertices; ++v) {
        const Vec3 p = linePoint(lines, v);
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    std::vector<std::pair<uint32_t, uint32_t>> order(fibres);  // (sort key, fibre)
    for (size_t i = 0; i < fibres; ++i) order[i] = {mortonCode(linePoint(lines, lines.fibreStarts[i]), lo, hi), static_cast<uint32_t>(i)};
    std::sort(order.begin(), order.end());
    // 2. The grid cells a chunk stays in: Z order jumps between cells, so a run of codes across a
    // cell boundary may join distant fibres. The finest level (2^level cells per axis) with at most
    // twice as many occupied cells as the chunks the vertices make.
    const size_t targetChunks = std::max<size_t>(1, vertices / kChunkVertices);
    int cellShift = 0;  // a fibre's cell: its code >> cellShift
    for (int level = 10; level >= 1; --level) {
        cellShift = 3 * (10 - level);
        size_t cells = 0;
        for (size_t i = 0; i < order.size(); ++i) cells += i == 0 || order[i].first >> cellShift != order[i - 1].first >> cellShift;
        if (cells <= 2 * targetChunks) break;
    }
    // 3. Chunks of consecutive fibres in one cell, each chunk's fibres in random order.
    DisplayLines out;
    for (size_t a = 0; a < order.size();) {
        size_t b = a;
        uint32_t chunkVertices = 0;
        const uint32_t cell = order[a].first >> cellShift;
        for (; b < order.size() && chunkVertices < kChunkVertices && order[b].first >> cellShift == cell; ++b)
            chunkVertices += fibreEnd(order[b].second) - lines.fibreStarts[order[b].second];
        for (size_t k = a; k < b; ++k) order[k].first = fibreRank(order[k].second);
        std::sort(order.begin() + static_cast<long long>(a), order.begin() + static_cast<long long>(b));
        FibreChunk chunk;
        chunk.firstVertex = static_cast<uint32_t>(out.pos.size() / 3);
        for (size_t k = a; k < b; ++k) appendLineVertices(lines, lines.fibreStarts[order[k].second], fibreEnd(order[k].second), out);
        chunk.vertexCount = static_cast<uint32_t>(out.pos.size() / 3) - chunk.firstVertex;
        if (chunk.vertexCount > 0) {
            boundChunk(out, chunk);
            out.chunks.push_back(chunk);
        }
        a = b;
    }
    lines = std::move(out);
}

} // namespace knit
