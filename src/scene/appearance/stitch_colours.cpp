// One colour per stitch from the base colour texture, reduced to the scene's palette; see stitch_colours.h.
#include "scene/appearance/stitch_colours.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <limits>
#include <memory>
#include <mutex>
#include <set>

#include "io/image_io.h"
#include "knitting/tube/plies.h"
#include "model/cancel.h"
#include "model/mesh_query.h"
#include "model/parallel.h"
#include "scene/scene.h"

namespace knit {

namespace {

// A part's surface, as its colour is read from it.
struct TexturedSurface {
    const TriMesh& mesh;
    const SourcePrimitive& prim;
    const std::vector<int>& sourceTriangles;  // part mesh triangle -> the primitive's triangle
    const ColourImage& texture;
};

// The texture uv (OBJ convention, v up) of a point of part mesh triangle `tri`: the primitive's colour
// uv of the source triangle, else its knitting uv (glTF's, v down). Never the part mesh's own uv: a
// spiral part's holds its chart (a, b).
Vec2 textureUvAt(const TexturedSurface& surface, int tri, const double bary[3]) {
    const int source = surface.sourceTriangles[tri];
    const bool hasColourUv = source < static_cast<int>(surface.prim.colourUv.size());
    const std::array<Vec2, 3>& uv = hasColourUv ? surface.prim.colourUv[source] : surface.prim.welded.uv[source];
    const Vec2 q{bary[0] * uv[0].x + bary[1] * uv[1].x + bary[2] * uv[2].x, bary[0] * uv[0].y + bary[1] * uv[1].y + bary[2] * uv[2].y};
    return hasColourUv ? q : Vec2{q.x, 1.0 - q.y};
}

// The linear mean of the texture over the stitch's cell; false: no sample is near the surface.
bool meanColourOverCell(const TexturedSurface& surface, const MeshQuery& query, const StitchRecord& stitch, double w, double h, Vec3& mean) {
    const double reach = 2 * std::max(w, h);
    Vec3 sum{0, 0, 0};
    int samples = 0;
    for (int i = 0; i < kCellSamples; ++i)
        for (int j = 0; j < kCellSamples; ++j) {
            const double a = (i + 0.5) / kCellSamples - 0.5, b = (j + 0.5) / kCellSamples;
            const MeshQuery::Hit hit = query.closest(stitch.origin + stitch.U * (a * w) + stitch.V * (b * h), reach);
            if (hit.tri < 0) continue;
            const std::array<int, 3>& v = surface.mesh.tri[hit.tri];
            const TrianglePoint onTri =
                closestOnTriangle(hit.p, surface.mesh.pos[v[0]], surface.mesh.pos[v[1]], surface.mesh.pos[v[2]]);
            sum = sum + bilinearTexel(surface.texture, textureUvAt(surface, hit.tri, onTri.bary));
            ++samples;
        }
    if (samples == 0) return false;
    mean = sum * (1.0 / samples);
    return true;
}

bool partHasStitches(const ScenePart& r) { return r.knit && !r.knit->t.curves.topo.stitches.empty(); }

// Whether a texture can be read on the part's stitches.
bool partHasTexture(const Scene& t, const ScenePart& r) { return partHasStitches(r) && partSurfaceHasTexture(t, r); }

Vec3 partObjectColour(const Scene& t, const ScenePart& r) {
    SceneSettings objectColour;
    objectColour.yarnColourSource = "object";
    return partYarnColour(t, r, objectColour);
}

// Part r's stitch colours; returns how many stitches took the part's mean colour.
int colourPartStitches(const Scene& t, const ScenePart& r, const MeshQuery* query, StitchColours& out) {
    const std::vector<StitchRecord>& stitches = r.knit->t.curves.topo.stitches;
    const uint32_t partKey = stitchColourKey(partObjectColour(t, r));
    out.key.assign(stitches.size(), partKey);
    if (!query) return 0;
    const SourcePrimitive& g = t.source.prims[r.prim];
    const TexturedSurface surface{*r.mesh, g, g.check.components[r.component].triangles, *g.albedo};
    const double w = r.knit->t.pp.waleSpacing, h = r.knit->t.pp.courseSpacing;
    std::atomic<int> offSurface{0};
    ExecPolicy parallel;
    parallel.parallel = true;
    parallelFor(static_cast<int64_t>(stitches.size()), 256, [&](int64_t begin, int64_t end, int64_t) {
        cancelCheckpoint();
        for (int64_t k = begin; k < end; ++k) {
            const StitchRecord& stitch = stitches[k];
            Vec3 mean;
            if (stitch.id >= 0 && stitch.id < static_cast<int>(out.key.size()) && meanColourOverCell(surface, *query, stitch, w, h, mean))
                out.key[stitch.id] = stitchColourKey(mean);
            else
                offSurface.fetch_add(1, std::memory_order_relaxed);
        }
    }, parallel);
    return offSurface.load();
}

// A colour of the median cut and the stitches that have it.
using WeightedColour = std::pair<uint32_t, long long>;

// One box of the median cut.
struct ColourBox {
    std::vector<WeightedColour> colours;
};

int channelOf(uint32_t key, int channel) { return key >> (16 - 8 * channel) & 0xff; }

// The channel with the widest range in the box, and that range.
int widestChannel(const ColourBox& box, int& range) {
    int best = 0;
    range = -1;
    for (int c = 0; c < 3; ++c) {
        int lo = 255, hi = 0;
        for (const WeightedColour& colour : box.colours) {
            lo = std::min(lo, channelOf(colour.first, c));
            hi = std::max(hi, channelOf(colour.first, c));
        }
        if (hi - lo > range) range = hi - lo, best = c;
    }
    return best;
}

// Moves the upper half of the box's stitches along `channel` into a new box.
ColourBox splitAtMedian(ColourBox& box, int channel) {
    std::sort(box.colours.begin(), box.colours.end(), [channel](const WeightedColour& a, const WeightedColour& b) {
        const int ca = channelOf(a.first, channel), cb = channelOf(b.first, channel);
        return ca != cb ? ca < cb : a.first < b.first;
    });
    long long total = 0;
    for (const WeightedColour& colour : box.colours) total += colour.second;
    long long below = 0;
    size_t cut = 1;
    for (; cut + 1 < box.colours.size(); ++cut) {
        below += box.colours[cut - 1].second;
        if (2 * below >= total) break;
    }
    ColourBox upper;
    upper.colours.assign(box.colours.begin() + cut, box.colours.end());
    box.colours.resize(cut);
    return upper;
}

uint32_t meanColourOfBox(const ColourBox& box) {
    Vec3 sum{0, 0, 0};
    double stitches = 0;
    for (const WeightedColour& colour : box.colours) {
        sum = sum + stitchColourLinear(colour.first) * double(colour.second);
        stitches += double(colour.second);
    }
    return stitchColourKey(sum * (1.0 / stitches));
}

} // namespace

bool meanTriangleColour(const SourcePrimitive& g, const std::vector<int>* tris, Vec3& out) {
    if (g.triColor.size() != g.welded.tri.size()) return false;
    const size_t count = tris ? tris->size() : g.welded.tri.size();
    double area = 0;
    Vec3 sum{0, 0, 0};
    for (size_t q = 0; q < count; ++q) {
        const int t = tris ? (*tris)[q] : static_cast<int>(q);
        const std::array<int, 3>& v = g.welded.tri[t];
        const double a = 0.5 * norm(cross(g.welded.pos[v[1]] - g.welded.pos[v[0]], g.welded.pos[v[2]] - g.welded.pos[v[0]]));
        area += a;
        sum = sum + g.triColor[t] * a;
    }
    if (area <= 0) return false;
    out = sum * (1.0 / area);
    return true;
}

Vec3 partYarnColour(const Scene& t, const ScenePart& r, const SceneSettings& s) {
    if (s.yarnColourSource != "object" || r.prim < 0 || r.prim >= static_cast<int>(t.source.prims.size())) return s.yarnWhite.color;
    const SourcePrimitive& g = t.source.prims[r.prim];
    if (r.component < 0 || r.component >= static_cast<int>(g.check.components.size())) return g.baseColor;
    Vec3 mean = g.baseColor;
    meanTriangleColour(g, &g.check.components[r.component].triangles, mean);
    return mean;
}

Vec3 stitchColourLinear(uint32_t key) {
    return {srgbByteToLinear(key >> 16 & 0xff), srgbByteToLinear(key >> 8 & 0xff), srgbByteToLinear(key & 0xff)};
}

uint32_t stitchColourKey(const Vec3& c) { return linearToSrgbByte(c.x) << 16 | linearToSrgbByte(c.y) << 8 | linearToSrgbByte(c.z); }

std::map<uint32_t, uint32_t> reducePalette(const std::map<uint32_t, long long>& stitchesOfColour, int maxColours) {
    std::map<uint32_t, uint32_t> entryOf;
    if (static_cast<int>(stitchesOfColour.size()) <= maxColours) {
        for (const std::pair<const uint32_t, long long>& colour : stitchesOfColour) entryOf[colour.first] = colour.first;
        return entryOf;
    }
    std::vector<ColourBox> boxes(1);
    boxes[0].colours.assign(stitchesOfColour.begin(), stitchesOfColour.end());
    while (static_cast<int>(boxes.size()) < maxColours) {
        size_t widest = 0;
        int widestRange = -1, widestAxis = 0;
        for (size_t b = 0; b < boxes.size(); ++b) {
            int range = 0;
            const int axis = widestChannel(boxes[b], range);
            if (range > widestRange) widestRange = range, widestAxis = axis, widest = b;
        }
        if (widestRange <= 0) break;
        boxes.push_back(splitAtMedian(boxes[widest], widestAxis));
    }
    for (const ColourBox& box : boxes) {
        const uint32_t entry = meanColourOfBox(box);
        for (const WeightedColour& colour : box.colours) entryOf[colour.first] = entry;
    }
    return entryOf;
}

bool partSurfaceHasTexture(const Scene& t, const ScenePart& r) {
    if (r.prim < 0 || r.prim >= static_cast<int>(t.source.prims.size())) return false;
    const SourcePrimitive& g = t.source.prims[r.prim];
    if (!g.albedo || r.component < 0 || r.component >= static_cast<int>(g.check.components.size())) return false;
    return g.colourUv.size() == g.welded.tri.size() || g.welded.uv.size() == g.welded.tri.size();
}

Vec3 partSurfaceColour(const Scene& t, const ScenePart& r, int tri, const double bary[3]) {
    const SourcePrimitive& g = t.source.prims[r.prim];
    const TexturedSurface surface{*r.mesh, g, g.check.components[r.component].triangles, *g.albedo};
    return bilinearTexel(*g.albedo, textureUvAt(surface, tri, bary));
}

namespace {

SceneStitchColours computeStitchColours(const Scene& t) {
    SceneStitchColours out;
    out.parts.resize(t.parts.size());
    bool anyTexture = false;
    for (const ScenePart& r : t.parts) anyTexture |= partHasTexture(t, r);
    if (!anyTexture) return out;
    // Chart parts share their component's mesh.
    std::map<const TriMesh*, std::unique_ptr<MeshQuery>> queries;
    std::map<uint32_t, long long> stitchesOfColour;
    for (size_t i = 0; i < t.parts.size(); ++i) {
        const ScenePart& r = t.parts[i];
        if (!partHasStitches(r)) continue;
        cancelCheckpoint();
        const MeshQuery* query = nullptr;
        if (partHasTexture(t, r)) {
            std::unique_ptr<MeshQuery>& q = queries[r.mesh.get()];
            if (!q) q = std::make_unique<MeshQuery>(*r.mesh);
            query = q.get();
        }
        out.stitchesOffSurface += colourPartStitches(t, r, query, out.parts[i]);
        for (uint32_t k : out.parts[i].key) ++stitchesOfColour[k];
    }
    const std::map<uint32_t, uint32_t> entryOf = reducePalette(stitchesOfColour, kYarnPaletteColours);
    std::set<uint32_t> palette;
    for (StitchColours& part : out.parts)
        for (uint32_t& k : part.key) {
            k = entryOf.at(k);
            palette.insert(k);
        }
    out.palette.assign(palette.begin(), palette.end());
    return out;
}

} // namespace

StitchAlongYarn::StitchAlongYarn(const CurveSet& cs, const PlySet* plies) : cs_(cs), plies_(plies) {}


int StitchAlongYarn::at(int strand, double s) const {
    int yarn = strand;
    if (plies_) {
        if (strand < 0 || strand >= static_cast<int>(plies_->curves.yarns.size())) return -1;
        const YarnPath& ply = plies_->curves.yarns[strand];
        const std::pair<int, double> parent = plyParent(*plies_, strand, ply.paramAtArcLength(std::clamp(s, 0.0, ply.length())));
        yarn = parent.first;
        s = parent.second;
    }
    if (yarn < 0 || yarn >= static_cast<int>(cs_.yarns.size())) return -1;
    const YarnPath& y = cs_.yarns[yarn];
    return cs_.topo.stitchNear(yarn, y.paramAtArcLength(std::clamp(s, 0.0, y.length())));
}

const SceneStitchColours& sceneStitchColours(const Scene& t) {
    static std::mutex one;  // two jobs must not compute the same knit's colours twice
    std::lock_guard<std::mutex> lock(one);
    if (!t.stitchColours) t.stitchColours = std::make_shared<const SceneStitchColours>(computeStitchColours(t));
    return *t.stitchColours;
}

} // namespace knit
