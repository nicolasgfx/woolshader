// The Fur stage: a loop pile and fine fibres grown on a scene's meshes (fur_stage.h has the steps).
#include "scene/fibres/fur_stage.h"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>

#include "fibres/fibre_budget.h"
#include "fibres/fur_field.h"
#include "fibres/fur_loops.h"
#include "io/files.h"
#include "io/image_io.h"
#include "model/cancel.h"
#include "model/hash.h"
#include "model/log.h"
#include "model/mesh_query.h"
#include "scene/appearance/stitch_colours.h"
#include "scene/scene.h"
#include "scene/source/scene_parts.h"
#include "scene/source/scene_solids.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"

namespace knit {

bool loadFurMask(const std::string& path, FurMask& out, std::string* err) {
    out = FurMask{};
    std::string why;
    if (!decodeGreyImageFile(path, out.width, out.height, out.value, &why)) {
        if (err) *err = "fur mask " + path + ": cannot read the image (" + why + ")";
        return false;
    }
    out.path = path;
    bool ok = false;
    out.hash = fnv1aFile(path, &ok);
    return true;
}

double sampleFurMask(const FurMask& m, const Vec2& uv) {
    if (m.value.empty()) return 1.0;
    const double u = uv.x - std::floor(uv.x), v = uv.y - std::floor(uv.y);
    const int x = std::clamp(int(u * m.width), 0, m.width - 1), y = std::clamp(int((1.0 - v) * m.height), 0, m.height - 1);
    return m.value[size_t(y) * size_t(m.width) + size_t(x)] / 255.0;
}

bool furMaskPath(const SceneSettings& s, std::string& path, std::string* err) {
    path.clear();
    if (s.pile.coverage == "entire") return true;
    // A mask image the user loaded: refused if the file changed since (its hash is the result key).
    if (s.pile.maskPath.empty()) {
        if (err) *err = "no mask image is loaded: Load one, or choose Entire surface";
        return false;
    }
    const std::string h = fileHashHex(s.pile.maskPath);
    if (h.empty()) {
        if (err) *err = "the mask image " + s.pile.maskPath + " cannot be read";
        return false;
    }
    if (!s.pile.maskHash.empty() && h != s.pile.maskHash) {
        if (err) *err = "the mask image " + s.pile.maskPath + " changed since it was loaded: Reload it";
        return false;
    }
    path = s.pile.maskPath;
    return true;
}

bool furPartGrows(const SceneSettings& s, const ScenePart& part) {
    const std::vector<std::string> excluded = furExcludedParts(s);
    return std::find(excluded.begin(), excluded.end(), part.name) == excluded.end();
}

bool furPartCoverage(const Scene& t, const ScenePart& r, const FurMask& mask, FurCoverage& out, std::string* err) {
    out = nullptr;
    if (mask.value.empty()) return true;
    const SourcePrimitive& g = t.source.prims[r.prim];
    const std::vector<int>& primTri = g.check.components[r.component].triangles;
    const std::vector<std::array<Vec2, 3>>& uvs = !g.colourUv.empty() ? g.colourUv : g.welded.uv;
    if (uvs.empty()) {
        if (err) *err = "part " + r.name + ": the mask image needs colour texture coordinates, which this mesh lacks";
        return false;
    }
    out = [&uvs, &primTri, &mask](int tri, const Vec3& b) {
        const std::array<Vec2, 3>& uv = uvs[static_cast<size_t>(primTri[static_cast<size_t>(tri)])];
        return sampleFurMask(mask, {b.x * uv[0].x + b.y * uv[1].x + b.z * uv[2].x, b.x * uv[0].y + b.y * uv[1].y + b.z * uv[2].y});
    };
    return true;
}

FurPartModifiers furPartModifiers(const Scene& t, size_t part, const FurLook& look) {
    const ScenePart& r = t.parts[part];
    TriMesh m = *r.mesh;
    for (Vec3& p : m.pos) p = p + r.centre;
    return buildFurModifiers(m, look, fnv1a(r.name.data(), r.name.size()));
}

FurPartRoots sampleFurPartRoots(const Scene& t, size_t part, const FurLook& look, double loopsCapacity, double fibresCapacity, const ExecPolicy& exec) {
    const ScenePart& r = t.parts[part];
    const SceneSolids& solids = *t.solids;
    const FurHiddenTest hidden = [&solids](const Vec3& q) { return solids.insideAny(q); };
    TriMesh m = *r.mesh;
    for (Vec3& p : m.pos) p = p + r.centre;
    FurPartRoots out;
    out.part = static_cast<int>(part);
    const uint64_t nameHash = fnv1a(r.name.data(), r.name.size());
    out.loops = sampleFurRoots(m, FurLayer::Loops, loopsCapacity, furPartSeed(look.seed, FurLayer::Loops, nameHash), hidden, exec);
    if (look.fibres.enabled && fibresCapacity > 0)
        out.fibres = sampleFurRoots(m, FurLayer::Fibres, fibresCapacity, furPartSeed(look.seed, FurLayer::Fibres, nameHash), hidden, exec);
    return out;
}

void addHairFibre(FibreSet& fs, const Hair& h, const Vec3& centre) {
    Fibre f;
    f.yarnId = f.parentYarn = h.tri;
    f.candidate = h.index;
    f.population = FibrePopulation::Directed;
    std::vector<Vec3> pts;
    double length = 0;
    for (const Vec3& p : h.points) {
        pts.push_back(p - centre);
        if (pts.size() > 1) length += norm(pts.back() - pts[pts.size() - 2]);
    }
    const double emergeL = norm(pts[1] - pts[0]);
    f.anchor = toVec3f(pts[1]);
    f.normal = toVec3f(h.normal);
    f.emergeL = float(emergeL);
    f.radius0 = float(h.rootRadius);
    f.radius1 = float(h.tipRadius);
    f.key = h.key;
    f.targetLength = float(length - emergeL);
    fs.addFibre(f, pts);
}

namespace {

// The base colour at each fibre's root (its points[1], on triangle h.tri of the part's mesh).
std::vector<uint32_t> hairRootColours(const Scene& t, const ScenePart& r, const std::vector<Hair>& hairs) {
    std::vector<uint32_t> keys;
    keys.reserve(hairs.size());
    for (const Hair& h : hairs) {
        const std::array<int, 3>& v = r.mesh->tri[static_cast<size_t>(h.tri)];
        const TrianglePoint onTri = closestOnTriangle(h.points[1] - r.centre, r.mesh->pos[v[0]], r.mesh->pos[v[1]], r.mesh->pos[v[2]]);
        keys.push_back(stitchColourKey(partSurfaceColour(t, r, h.tri, onTri.bary)));
    }
    return keys;
}

// The fibre colours of all parts reduced to one palette (the renderer has one colour per material); returns
// the palette's size.
size_t reduceHairColours(std::vector<PartFibres>& parts) {
    std::map<uint32_t, long long> hairsOfColour;
    for (const PartFibres& pf : parts)
        for (uint32_t k : pf.fibreColour) ++hairsOfColour[k];
    const std::map<uint32_t, uint32_t> entryOf = reducePalette(hairsOfColour, kYarnPaletteColours);
    std::set<uint32_t> palette;
    for (PartFibres& pf : parts)
        for (uint32_t& k : pf.fibreColour) {
            k = entryOf.at(k);
            palette.insert(k);
        }
    return palette.size();
}

// The body fibres collide with: every non-skipped primitive (library mm).
TriMesh furBody(const Scene& t) {
    TriMesh body;
    for (const SourcePrimitive& g : t.source.prims) {
        if (g.cls == SourcePrimitive::Class::Skip) continue;
        const int base = static_cast<int>(body.pos.size());
        body.pos.insert(body.pos.end(), g.welded.pos.begin(), g.welded.pos.end());
        for (const auto& tr : g.welded.tri) body.tri.push_back({tr[0] + base, tr[1] + base, tr[2] + base});
    }
    body.computeNormals();
    return body;
}

// One part's fur field (the root caches are built for it and dropped: the viewer's live preview keeps them).
struct PartField {
    int part = -1;
    FurField field;
};

// The fields of every part but the excluded ones: the budget is global, so they are all built first.
bool buildPartFields(const Scene& scene, const SceneSettings& s, const FurLook& look, const FurMask& mask, std::vector<PartField>& out,
                     int& excludedParts, FurReport& rep, std::string* err) {
    const ExecPolicy exec{true};
    for (size_t i = 0; i < scene.parts.size(); ++i) {
        cancelCheckpoint();
        if (!furPartGrows(s, scene.parts[i])) {
            ++excludedParts;
            continue;
        }
        FurCoverage coverage;
        if (!furPartCoverage(scene, scene.parts[i], mask, coverage, err)) return false;
        const auto t0 = std::chrono::steady_clock::now();
        const FurPartRoots roots =
            sampleFurPartRoots(scene, i, look, furDensityBound(look, FurLayer::Loops), furDensityBound(look, FurLayer::Fibres), exec);
        rep.msSample += msSince(t0);
        const auto t1 = std::chrono::steady_clock::now();
        const FurPartModifiers mods = furPartModifiers(scene, i, look);
        PartField pf;
        pf.part = static_cast<int>(i);
        pf.field = buildFurField({roots.loops, roots.fibres, look, coverage, &mods}, exec);
        rep.msField += msSince(t1);
        for (FurLayer layer : {FurLayer::Loops, FurLayer::Fibres}) {
            FurLayerCounts &to = layer == FurLayer::Loops ? rep.loops : rep.fibres;
            const FurLayerCounts& from = layer == FurLayer::Loops ? pf.field.loops : pf.field.fibres;
            to.sampled += from.sampled;
            to.hidden += from.hidden;
            to.masked += from.masked;
            to.thinned += from.thinned;
            to.kept += from.kept;
        }
        out.push_back(std::move(pf));
    }
    return true;
}

// The segment budget before the settle: the export's threshold over (key, polyline segments) of every disc of
// both layers, in part order, the loops before the fine fibres. Returns the threshold.
float applySegmentBudget(std::vector<PartField>& parts, long long budget, FurReport& rep) {
    std::vector<std::pair<float, int>> keySegments;
    for (const PartField& pf : parts)
        for (const FurDisc& d : pf.field.discs) keySegments.push_back({d.budgetKey, furSegmentCount(d) + 1});
    const float q = budgetThreshold(keySegments, budget);
    for (PartField& pf : parts) {
        std::vector<FurDisc>& discs = pf.field.discs;
        for (const FurDisc& d : discs)
            if (!(d.budgetKey < q)) ++(d.layer == FurLayer::Loops ? rep.loops : rep.fibres).overBudget;
        discs.erase(std::remove_if(discs.begin(), discs.end(), [q](const FurDisc& d) { return !(d.budgetKey < q); }), discs.end());
    }
    return q;
}

std::string layerLine(const char* name, const FurLayerCounts& c) {
    char b[400];
    const std::string thinned = c.thinned ? std::to_string(c.thinned) + " thinned by the local density, " : std::string();
    std::snprintf(b, sizeof b,
                  "%s: %lld sampled, %lld inside another solid, %lld outside the coverage, %s%lld kept, %lld over the budget, %lld "
                  "stuck in the body (rejected), %lld grown (%lld lengths raised to the minimum), %lld segments",
                  name, static_cast<long long>(c.sampled), static_cast<long long>(c.hidden), static_cast<long long>(c.masked),
                  thinned.c_str(), static_cast<long long>(c.kept), static_cast<long long>(c.overBudget), static_cast<long long>(c.stuck),
                  static_cast<long long>(c.grown), static_cast<long long>(c.clamped), static_cast<long long>(c.segments));
    return b;
}

} // namespace

std::optional<SceneFibres> growSceneFur(const SceneEntry&, const Scene& t, const SceneSettings& s, std::string* err) {
    const auto started = std::chrono::steady_clock::now();
    if (!pileSettingsValid(s, err)) return std::nullopt;
    const FurLook look = furLookFromSettings(s);
    std::string maskPath;
    if (!furMaskPath(s, maskPath, err)) return std::nullopt;
    FurMask mask;
    if (!maskPath.empty() && !loadFurMask(maskPath, mask, err)) return std::nullopt;
    const TriMesh body = furBody(t);
    const MeshQuery bodyQuery(body);
    FurReport rep;
    std::vector<PartField> parts;
    int excluded = 0;
    if (!buildPartFields(t, s, look, mask, parts, excluded, rep, err)) return std::nullopt;
    const float q = applySegmentBudget(parts, s.segmentBudget, rep);
    SceneFibres out;
    for (const PartField& pf : parts) {
        cancelCheckpoint();
        const ScenePart& r = t.parts[static_cast<size_t>(pf.part)];
        const FurReport before = rep;
        const std::vector<Hair> hairs = growLoops(pf.field, bodyQuery, look.settle, ExecPolicy{true}, &rep);
        FibreSet fs;
        fs.rep.candidates = static_cast<int>(pf.field.discs.size());
        for (const Hair& h : hairs) addHairFibre(fs, h, r.centre);
        fs.rep.accepted = static_cast<int>(fs.fibres.size());
        PartFibres part{pf.part, r.name, std::move(fs), {}};
        if (partSurfaceHasTexture(t, r)) part.fibreColour = hairRootColours(t, r, hairs);
        appendFormat(out.log, "  %s fur: %lld loops, %lld fine fibres grown\n", r.name.c_str(), static_cast<long long>(rep.loops.grown - before.loops.grown),
                     static_cast<long long>(rep.fibres.grown - before.fibres.grown));
        out.parts.push_back(std::move(part));
    }
    const size_t colours = reduceHairColours(out.parts);
    // Fur keeps one segment budget over all its parts (wool's holds per object): every kept fibre is below q,
    // so the export writes them all.
    for (PartFibres& pf : out.parts) pf.keyBelow = q;
    appendFormat(out.log, "fur: %s\n", layerLine("loops", rep.loops).c_str());
    appendFormat(out.log, "fur: %s\n", layerLine("fine fibres", rep.fibres).c_str());
    appendFormat(out.log, "fur: coverage %s%s%s, %d excluded part%s; segment budget %lld (key < %.6f); largest last-step move %.4f mm\n",
                 s.pile.coverage.c_str(), maskPath.empty() ? "" : " ", maskPath.c_str(), excluded, excluded == 1 ? "" : "s", s.segmentBudget,
                 double(q), rep.maxResidual);
    appendFormat(out.log, "fur: sample %.0f ms, field %.0f ms, settle %.0f ms; %zu texture colours at the roots; wall time %.1f s\n", rep.msSample,
                 rep.msField, rep.msSettle, colours, secondsSince(started));
    return out;
}

} // namespace knit
