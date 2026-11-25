// From a scene entry to its parts: load the source, make one part per component or chart.
#include "scene/source/scene_parts.h"

#include <algorithm>

#include "model/hash.h"
#include "import/import_store.h"
#include "knitting/knit/pieces.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"
#include "knitting/knit/mesh_knit.h"
#include "io/files.h"
#include "scene/scene.h"
#include "scene/source/scene_solids.h"

namespace knit {

namespace {

// The settings a knitting run refuses before it starts (the reason in err). `pattern` is parsed here.
bool checkKnitSettings(const SceneEntry& e, const SceneSettings& s, StitchPattern& pattern, std::string* err) {
    const std::optional<StitchPattern> parsed = parsePattern(s.pattern);
    if (!parsed) {
        if (err) *err = "unknown pattern " + s.pattern;
        return false;
    }
    pattern = *parsed;
    if (std::string why; !finalSurface(s, s.diameter).valid(&why)) {
        if (err) *err = "scene " + e.id + ": " + why;
        return false;
    }
    return true;
}

bool hasPartNamed(const std::vector<ScenePart>& parts, const std::string& name) {
    // A component's name too: with knit.coverage = overlap its parts are named after it (<name>_ch<k>).
    for (const ScenePart& q : parts)
        if (q.name == name || q.componentName == name) return true;
    return false;
}

// The parts to knit: one per component of every Knit primitive, or with knit.coverage = overlap one
// per usable chart, knitted on the whole component. Every part of an import is knitted on its charts.
bool collectParts(const SceneSettings& s, bool requireCharts, Scene& out, std::string* err) {
    for (size_t pi = 0; pi < out.source.prims.size(); ++pi) {
        const SourcePrimitive& g = out.source.prims[pi];
        if (g.cls != SourcePrimitive::Class::Knit) continue;
        for (size_t k = 0; k < g.check.components.size(); ++k) {
            ScenePart r;
            r.prim = static_cast<int>(pi);
            r.component = static_cast<int>(k);
            r.name = g.mesh + "_" + g.material + "_c" + std::to_string(k);
            // Instances of one mesh (several nodes) give equal names: the export files are
            // named after the part, so the name must be unique.
            for (int n = 2; hasPartNamed(out.parts, r.name); ++n)
                r.name = g.mesh + "_" + g.material + "_c" + std::to_string(k) + "_" + std::to_string(n);
            r.mesh = std::make_shared<const TriMesh>(extractTriangles(g.welded, g.check.components[k].triangles, r.centre));
            r.componentName = r.name;
            if (g.welded.uv.empty() && requireCharts) {
                if (err) *err = "primitive " + g.mesh + " has no TEXCOORD_0: cannot be tiled";
                return false;
            }
            if (overlapCoverage(s) && !r.mesh->uv.empty()) {
                const UvAtlas atlas(*r.mesh, TileParams{}.up);
                size_t chartParts = 0;
                for (size_t ci = 0; ci < atlas.charts().size(); ++ci) {
                    if (atlas.charts()[ci].status != Status::Pass) continue;
                    ScenePart rc;
                    rc.prim = r.prim;
                    rc.component = r.component;
                    rc.name = r.name + "_ch" + std::to_string(ci);
                    rc.componentName = r.name;
                    rc.centre = r.centre;
                    rc.mesh = r.mesh;
                    rc.chart = static_cast<int>(ci);
                    out.parts.push_back(std::move(rc));
                    ++chartParts;
                }
                // No chart passes the atlas checks: the component stays one part (knitted inside its
                // charts, or left out of the export), so its surface never disappears from the knit.
                if (chartParts > 0) continue;
            }
            out.parts.push_back(std::move(r));
        }
    }
    return true;
}

// The source scene of a mesh layout: the GLB or the import, and its solids (the parts come next,
// collectParts).
bool loadMeshScene(const SceneEntry& e, const SceneSettings& s, const SourceOptions& options, Scene& out, std::string* err) {
    const std::string path = findProjectFile(e.path);
    if (path.empty()) {
        if (err) *err = "scene " + e.id + ": " + e.path + " not found";
        return false;
    }
    if (e.kind == SceneEntry::Kind::Mesh) {
        if (!loadImportedScene(path, out.source, err, options)) return false;
    } else if (!loadGlb(path, e.glb, out.source, err, options)) {
        return false;
    }
    out.solids = std::make_shared<const SceneSolids>(out.source);
    return true;
}

void hashFile(const std::string& path, uint64_t& h, bool& ok) {
    bool fileOk = false;
    const uint64_t f = fnv1aFile(path, &fileOk);
    ok = ok && fileOk;
    h = fnv1a(&f, sizeof f, h);
}

} // namespace

uint64_t sceneAssetHash(const SceneEntry& e, std::string* err) {
    const bool manifestOnly = e.path.empty() && e.knitsAsPieces();  // knit_check --pieces: no GLB beside the pieces
    const std::string path = manifestOnly ? std::string() : findProjectFile(e.path);
    if (path.empty() && !manifestOnly) {
        if (err) *err = e.path + " not found";
        return 0;
    }
    uint64_t h = fnv1a(e.id.data(), e.id.size());
    bool ok = true;
    // An import: the prepared mesh (its directory name is the import identity, which covers
    // the source, its MTL / images and the options; the textures are copies under it).
    if (!manifestOnly) hashFile(e.kind == SceneEntry::Kind::Mesh ? path + "/mesh.kmesh" : path, h, ok);
    // Knitted as spiral pieces (the teddy): the pieces are the knit's input and the GLB the colour
    // source, the kept parts' source and the skipped-knit mesh, so every file of both; the routing token
    // keeps these entries apart from the ones of the same GLB knitted on its UV charts. v2: the entry holds
    // the common knit result (knit.bin; v1 entries referred to knit cache folders and now miss). v5: the kept
    // parts are the GLB's unselected nodes, no longer the manifest's extras files (a v4 entry holds those).
    if (e.knitsAsPieces()) {
        constexpr char kRouting[] = "knit: spiral pieces v5 (kept parts from the GLB selection)";
        h = fnv1a(kRouting, sizeof kRouting - 1, h);
        const std::string manifest = findProjectFile(e.piecesManifest);
        if (manifest.empty()) {
            if (err) *err = e.piecesManifest + " not found";
            return 0;
        }
        hashFile(manifest, h, ok);
        Manifest m;
        if (!readManifest(manifest, m, err)) return 0;
        for (const PieceEntry& p : m.pieces) hashFile(piecePath(m, p.file), h, ok);
    }
    if (!ok) {
        if (err) *err = "cannot read the inputs of scene " + e.id;
        return 0;
    }
    return h;
}

std::optional<KnitSource> loadPartsToKnit(const SceneEntry& e, const SceneSettings& s, const SourceOptions& options, Scene& out,
                                          std::string* err) {
    KnitSource source;
    if (knitsGuidedPiece(s)) {
        SourceOptions bare = options;
        bare.requireCharts = false;  // the stitch mesh needs no charts: the import may keep the file's surface
        std::shared_ptr<GuidedInputs> guided = std::make_shared<GuidedInputs>();
        if (!loadMeshScene(e, s, bare, out, err) || !checkKnitSettings(e, s, source.pattern, err) || !loadGuidedParts(s, out, *guided, err)) return std::nullopt;
        source.guided = guided;
        return source;
    }
    if (knitsSewnPieces(e, s)) {
        if (!checkKnitSettings(e, s, source.pattern, err) || !checkSpiralSettings(e, s, err) || !loadPiecesScene(e, s, source.pieces, out, err))
            return std::nullopt;
        return source;
    }
    SourceOptions charted = options;
    charted.requireCharts = true;
    if (!loadMeshScene(e, s, charted, out, err) || !checkKnitSettings(e, s, source.pattern, err) || !collectParts(s, true, out, err))
        return std::nullopt;
    return source;
}

bool loadBaseParts(const SceneEntry& e, const SceneSettings& s, const SourceOptions& options, Scene& out, std::string* err) {
    SourceOptions bare = options;
    bare.requireCharts = false;
    if (!loadMeshScene(e, s, bare, out, err) || !collectParts(s, false, out, err)) return false;
    for (ScenePart& r : out.parts) r.note = "Knit skipped: the part stays its mesh";
    return true;
}

} // namespace knit
