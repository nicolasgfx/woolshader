// GLB scenes into library space, one scene-wide similarity for every mesh (glb_load.h has the units contract).
#include "io/glb_load.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>

#include "cgltf.h"
#include "io/files.h"
#include "io/image_io.h"
#include "io/mesh_read.h"

namespace knit {

namespace {

// glTF world (Y-up) -> library Z-up, before scaling.
Vec3 zUp(const Vec3& p) { return {p.x, -p.z, p.y}; }

bool has(const std::vector<std::string>& v, const std::string& s) { return std::find(v.begin(), v.end(), s) != v.end(); }

// A base-colour image of the file (embedded buffer view or a file next to the GLB), decoded within
// the import texture limits (mesh_read.h) and kept at most 2048 px on a side: enough for the
// display and the per-triangle base colour, a fraction of an 8k image's memory. Null if it cannot be.
std::shared_ptr<const ColourImage> decodeImage(const cgltf_image* im, const std::string& glbPath, std::vector<std::string>& log) {
    constexpr int kKeptSide = 2048;
    const int maxSide = ImportLimits{}.maxTextureSide;
    ColourImage image;
    bool decoded = false;
    std::string name = im->name ? im->name : "";
    if (im->buffer_view && im->buffer_view->buffer && im->buffer_view->buffer->data) {
        const auto* bytes = static_cast<const unsigned char*>(im->buffer_view->buffer->data) + im->buffer_view->offset;
        decoded = decodeImageBytes(bytes, im->buffer_view->size, maxSide, image, nullptr);
        if (name.empty()) name = "(embedded)";
    } else if (im->uri && std::string(im->uri).rfind("data:", 0) != 0) {
        const std::string p = (std::filesystem::path(glbPath).parent_path() / im->uri).string();
        int w = 0, h = 0;
        decoded = readImageSize(p, w, h, nullptr) && w <= maxSide && h <= maxSide && decodeImageFile(p, image, nullptr);
        name = p;
    }
    if (!decoded) {
        log.push_back("base colour image " + name + ": not decoded (missing, unreadable or larger than 8192 px); its material colour is used");
        return nullptr;
    }
    image.path = name;
    return std::make_shared<const ColourImage>(shrunkToSide(image, kKeptSide));
}

// The parsed file, freed on every way out.
class CgltfFile {
public:
    CgltfFile() = default;
    CgltfFile(const CgltfFile&) = delete;
    CgltfFile& operator=(const CgltfFile&) = delete;
    ~CgltfFile() {
        if (data) cgltf_free(data);
    }
    cgltf_data* data = nullptr;
};

struct RawPrim {
    SourcePrimitive meta;
    TriMesh mesh;                              // unscaled Z-up glTF units, unwelded
    std::vector<std::array<Vec3, 3>> normals;  // per triangle corner
    TextureTransform colourTransform;          // of the base colour texture
};

// A node's world transform (column-major), and its normal matrix: the inverse transpose of the 3x3
// (cofactors / det).
struct NodeTransform {
    double M[16];
    double C[9];
    double det;

    explicit NodeTransform(const cgltf_node* n) {
        float world[16];
        cgltf_node_transform_world(n, world);
        for (int i = 0; i < 16; ++i) M[i] = world[i];
        det = M[0] * (M[5] * M[10] - M[9] * M[6]) - M[4] * (M[1] * M[10] - M[9] * M[2]) + M[8] * (M[1] * M[6] - M[5] * M[2]);
        const double cofactors[9] = {M[5] * M[10] - M[6] * M[9], M[6] * M[8] - M[4] * M[10], M[4] * M[9] - M[5] * M[8],
                                     M[2] * M[9] - M[1] * M[10], M[0] * M[10] - M[2] * M[8], M[1] * M[8] - M[0] * M[9],
                                     M[1] * M[6] - M[2] * M[5], M[2] * M[4] - M[0] * M[6], M[0] * M[5] - M[1] * M[4]};
        std::copy(cofactors, cofactors + 9, C);
    }
    Vec3 point(const Vec3& p) const {
        return {M[0] * p.x + M[4] * p.y + M[8] * p.z + M[12], M[1] * p.x + M[5] * p.y + M[9] * p.z + M[13], M[2] * p.x + M[6] * p.y + M[10] * p.z + M[14]};
    }
    Vec3 normal(const Vec3& v) const {
        const Vec3 r{C[0] * v.x + C[3] * v.y + C[6] * v.z, C[1] * v.x + C[4] * v.y + C[7] * v.z, C[2] * v.x + C[5] * v.y + C[8] * v.z};
        const double l = norm(r);
        return l > 0 ? r * (det < 0 ? -1.0 / l : 1.0 / l) : r;
    }
};

// What walking the file's nodes gathers. A problem stops the walk.
struct GlbWalk {
    const cgltf_data* data;
    const std::string& path;
    std::vector<std::string>& log;
    std::map<const cgltf_image*, std::shared_ptr<const ColourImage>> images;  // decoded once per image
    std::vector<RawPrim> prims;
    std::string problem;
};

bool openGlb(const std::string& path, CgltfFile& file, std::vector<std::string>& log, std::string* err) {
    cgltf_options opt{};
    if (cgltf_parse_file(&opt, path.c_str(), &file.data) != cgltf_result_success ||
        cgltf_load_buffers(&opt, file.data, path.c_str()) != cgltf_result_success || cgltf_validate(file.data) != cgltf_result_success)
        return fail(err, path + ": not a valid glTF / GLB file");
    const cgltf_data* data = file.data;
    for (size_t k = 0; k < data->extensions_used_count; ++k) {
        const std::string x = data->extensions_used[k];
        log.push_back("extension " + x + (x == "KHR_texture_transform" ? ": applied to the base colour's texture coordinates"
                                          : x == "KHR_materials_sheen" ? ": ignored (materials come from the settings)"
                                                                       : ": ignored"));
    }
    if (data->extensions_required_count > 0) return fail(err, path + ": requires extension " + data->extensions_required[0] + " (not supported)");
    return true;
}

// The base colour only (metallic, roughness, normal and other maps are not read).
void readBaseColour(const cgltf_primitive& pr, GlbWalk& walk, RawPrim& rp) {
    if (!pr.material || !pr.material->has_pbr_metallic_roughness) return;
    const cgltf_pbr_metallic_roughness& pbr = pr.material->pbr_metallic_roughness;
    rp.meta.baseColor = {pbr.base_color_factor[0], pbr.base_color_factor[1], pbr.base_color_factor[2]};
    const cgltf_texture* tx = pbr.base_color_texture.texture;
    if (!tx || !tx->image) return;
    rp.meta.baseColorImage = tx->image->uri ? tx->image->uri : tx->image->name ? tx->image->name : "(embedded)";
    auto it = walk.images.find(tx->image);
    if (it == walk.images.end()) it = walk.images.emplace(tx->image, decodeImage(tx->image, walk.path, walk.log)).first;
    rp.meta.albedo = it->second;
    if (pbr.base_color_texture.has_transform) {
        const cgltf_texture_transform& x = pbr.base_color_texture.transform;
        rp.colourTransform = {{x.offset[0], x.offset[1]}, x.rotation, {x.scale[0], x.scale[1]}};
    }
}

void readPrimitive(const cgltf_node* n, size_t index, const NodeTransform& world, GlbWalk& walk) {
    const cgltf_primitive& pr = n->mesh->primitives[index];
    RawPrim rp;
    rp.meta.node = n->name ? n->name : "";
    rp.meta.mesh = n->mesh->name ? n->mesh->name : "";
    rp.meta.primitive = static_cast<int>(index);
    rp.meta.nodeIndex = static_cast<int>(n - walk.data->nodes);
    rp.meta.material = pr.material && pr.material->name ? pr.material->name : "";
    rp.meta.mirrored = world.det < 0;
    const std::string where = "node " + rp.meta.node + " mesh " + rp.meta.mesh + " primitive " + std::to_string(index);
    if (pr.type != cgltf_primitive_type_triangles) {
        walk.problem = where + ": not a triangle list";
        return;
    }
    readBaseColour(pr, walk, rp);
    const cgltf_accessor *P = nullptr, *N = nullptr, *UV = nullptr;
    for (size_t a = 0; a < pr.attributes_count; ++a) {
        const cgltf_attribute& at = pr.attributes[a];
        if (at.type == cgltf_attribute_type_position) P = at.data;
        else if (at.type == cgltf_attribute_type_normal) N = at.data;
        else if (at.type == cgltf_attribute_type_texcoord && at.index == 0) UV = at.data;
    }
    if (!P) {
        walk.problem = where + ": no POSITION";
        return;
    }
    const size_t nv = P->count;
    rp.mesh.pos.resize(nv);
    std::vector<Vec3> vn(nv);
    std::vector<Vec2> vuv(UV ? nv : 0);
    for (size_t v = 0; v < nv; ++v) {
        float f[3] = {0, 0, 0};
        cgltf_accessor_read_float(P, v, f, 3);
        rp.mesh.pos[v] = zUp(world.point({f[0], f[1], f[2]}));
        if (N) {
            cgltf_accessor_read_float(N, v, f, 3);
            vn[v] = zUp(world.normal({f[0], f[1], f[2]}));
        }
        if (UV) {
            cgltf_accessor_read_float(UV, v, f, 2);
            vuv[v] = {f[0], f[1]};
        }
    }
    const size_t ni = pr.indices ? pr.indices->count : nv;
    if (ni % 3) {
        walk.problem = where + ": index count not a multiple of 3";
        return;
    }
    for (size_t t = 0; t < ni / 3; ++t) {
        std::array<int, 3> tri;
        for (int k = 0; k < 3; ++k) tri[k] = static_cast<int>(pr.indices ? cgltf_accessor_read_index(pr.indices, 3 * t + k) : 3 * t + k);
        if (world.det < 0) std::swap(tri[1], tri[2]);
        rp.mesh.tri.push_back(tri);
        if (UV) rp.mesh.uv.push_back({vuv[tri[0]], vuv[tri[1]], vuv[tri[2]]});
        rp.normals.push_back({vn[tri[0]], vn[tri[1]], vn[tri[2]]});
    }
    walk.prims.push_back(std::move(rp));
}

void readNode(const cgltf_node* n, GlbWalk& walk) {
    if (n->mesh) {
        const NodeTransform world(n);
        for (size_t p = 0; p < n->mesh->primitives_count; ++p) {
            readPrimitive(n, p, world, walk);
            if (!walk.problem.empty()) return;
        }
    }
    for (size_t c = 0; c < n->children_count && walk.problem.empty(); ++c) readNode(n->children[c], walk);
}

// A primitive's class: the node selection first, then the mesh and material lists, with no silent default.
bool classOf(const GlbSelection& selection, const SourcePrimitive& p, SourcePrimitive::Class& cls, const std::string& path, std::string* err) {
    const std::string& m = p.material;
    if (has(selection.plainNodes, p.node)) cls = SourcePrimitive::Class::Plain;
    else if (has(selection.skipMeshes, p.mesh)) cls = SourcePrimitive::Class::Skip;
    else if (has(selection.plainMeshes, p.mesh)) cls = SourcePrimitive::Class::Plain;
    else if (has(selection.skipMaterials, m)) cls = SourcePrimitive::Class::Skip;
    else if (has(selection.plainMaterials, m)) cls = SourcePrimitive::Class::Plain;
    else if (selection.knitAllMaterials || has(selection.knitMaterials, m)) cls = SourcePrimitive::Class::Knit;
    else return fail(err, path + ": material \"" + m + "\" is in none of the scene's lists");
    return true;
}

bool classifyPrimitives(std::vector<RawPrim>& raw, const GlbSelection& selection, const std::string& path, std::string* err) {
    for (RawPrim& rp : raw)
        if (!classOf(selection, rp.meta, rp.meta.cls, path, err)) return false;
    return true;
}

// The nodes below n (n included) that have a mesh, for listGlbNodes.
bool listNode(const cgltf_node* n, const GlbSelection& selection, const std::string& path, std::vector<GlbNode>& out, std::string* err) {
    if (n->mesh) {
        GlbNode row;
        row.name = n->name ? n->name : "";
        bool anyKnitted = false;
        for (size_t k = 0; k < n->mesh->primitives_count; ++k) {
            const cgltf_primitive& pr = n->mesh->primitives[k];
            SourcePrimitive p;
            p.node = row.name;
            p.mesh = n->mesh->name ? n->mesh->name : "";
            p.material = pr.material && pr.material->name ? pr.material->name : "";
            SourcePrimitive::Class cls;
            if (!classOf(selection, p, cls, path, err)) return false;
            if (cls == SourcePrimitive::Class::Skip) continue;
            anyKnitted |= cls == SourcePrimitive::Class::Knit;
            const cgltf_size corners = pr.indices ? pr.indices->count : pr.attributes_count ? pr.attributes[0].data->count : 0;
            row.triangles += static_cast<long long>(corners / 3);
        }
        row.knitted = anyKnitted;
        if (row.triangles > 0) out.push_back(row);
    }
    for (size_t c = 0; c < n->children_count; ++c)
        if (!listNode(n->children[c], selection, path, out, err)) return false;
    return true;
}

// The scene-wide scale from the bounds of all included meshes (Z-up glTF units).
bool scaleToLibrary(const std::vector<RawPrim>& raw, const GlbSelection& selection, const std::string& path, SourceScene& out, std::string* err) {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const RawPrim& rp : raw) {
        if (rp.meta.cls == SourcePrimitive::Class::Skip) continue;
        for (const auto& t : rp.mesh.tri)
            for (int v : t) {
                const Vec3& p = rp.mesh.pos[v];
                lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
                hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
            }
    }
    const bool byHeight = selection.sizeRef == GlbSelection::SizeRef::Height;
    const double measure = byHeight ? hi.z - lo.z : std::max(hi.x - lo.x, hi.y - lo.y);
    if (!(measure > 0) || !(selection.sizeMm > 0)) return fail(err, path + ": empty scene or no size");
    out.scale = selection.sizeMm / measure;
    out.offset = Vec3{-0.5 * (lo.x + hi.x), -0.5 * (lo.y + hi.y), -lo.z} * out.scale;
    out.lo = lo * out.scale + out.offset;
    out.hi = hi * out.scale + out.offset;
    char b[256];
    std::snprintf(b, sizeof b, "scale %.6g mm per glTF unit (%s %.6g units -> %.1f mm)", out.scale, byHeight ? "height" : "horizontal extent", measure, selection.sizeMm);
    out.log.push_back(b);
    return true;
}

// The base colour per triangle (texture x factor; the per-part yarn colour averages it). The colour
// uv follows the OBJ convention (v up): glTF's v runs down the image.
void colourFromTexture(const RawPrim& rp, SourcePrimitive& g) {
    g.colourUv.resize(g.welded.uv.size());
    g.triColor.resize(g.welded.uv.size());
    for (size_t t = 0; t < g.welded.uv.size(); ++t) {
        for (int k = 0; k < 3; ++k) {
            const Vec2 q = transformTextureUv(rp.colourTransform, g.welded.uv[t][k]);
            g.colourUv[t][k] = {q.x, 1.0 - q.y};
        }
        const Vec3 c = triangleTextureColour(*g.albedo, g.colourUv[t]);
        g.triColor[t] = {c.x * g.baseColor.x, c.y * g.baseColor.y, c.z * g.baseColor.z};
    }
}

// Into library space, welded per source (coincident vertices of different sources are never merged).
void addPrimitive(RawPrim& rp, const SourceOptions& o, SourceScene& out) {
    SourcePrimitive g = std::move(rp.meta);
    for (Vec3& p : rp.mesh.pos) p = p * out.scale + out.offset;
    Vec3 plo, phi;
    rp.mesh.bounds(plo, phi);
    WeldedMesh weld = weldVertices(rp.mesh, o.weldRel * norm(phi - plo));
    g.welded = std::move(weld.mesh);
    g.droppedTriangles = weld.droppedTriangles;
    const std::vector<int>& remap = weld.remap;
    // The corner normals stay aligned with the welded triangles (degenerates dropped).
    for (size_t t = 0; t < rp.mesh.tri.size(); ++t) {
        const std::array<int, 3>& tr = rp.mesh.tri[t];
        const int a = remap[tr[0]], b = remap[tr[1]], c = remap[tr[2]];
        if (a != b && b != c && a != c) g.cornerNormals.push_back(rp.normals[t]);
    }
    // A skipped primitive is never pair-tested.
    const bool pairTest = o.selfIntersection && g.cls != SourcePrimitive::Class::Skip;
    g.check = checkMesh(g.welded, pairTest ? SelfIntersectionTest::Run : SelfIntersectionTest::Skip);
    if (g.albedo && !g.welded.uv.empty()) colourFromTexture(rp, g);
    char b[256];
    std::snprintf(b, sizeof b, "%s / %s / %s [%s]: %zu verts, %zu tris (%d degenerate dropped)%s; %zu components", g.node.c_str(), g.mesh.c_str(),
                  g.material.c_str(), g.cls == SourcePrimitive::Class::Knit ? "knit" : g.cls == SourcePrimitive::Class::Plain ? "plain" : "skip",
                  g.welded.pos.size(), g.welded.tri.size(), g.droppedTriangles, g.mirrored ? ", mirrored node" : "", g.check.components.size());
    out.log.push_back(b);
    out.prims.push_back(std::move(g));
}

} // namespace

Vec2 transformTextureUv(const TextureTransform& t, const Vec2& uv) {
    const double c = std::cos(t.rotation), s = std::sin(t.rotation), u = uv.x * t.scale.x, v = uv.y * t.scale.y;
    return {t.offset.x + c * u + s * v, t.offset.y - s * u + c * v};
}

bool loadGlb(const std::string& path, const GlbSelection& selection, SourceScene& out, std::string* err, const SourceOptions& o) {
    out = SourceScene{};
    out.path = path;
    bool ok = false;
    out.fileHash = fnv1aFile(path, &ok);
    if (!ok) return fail(err, "cannot read " + path);
    std::vector<RawPrim> raw;
    {   // The meshes are copied: the file's buffers are not held through welding.
        CgltfFile file;
        if (!openGlb(path, file, out.log, err)) return false;
        GlbWalk walk{file.data, path, out.log, {}, {}, {}};
        const cgltf_data* data = file.data;
        const cgltf_scene* scene = data->scene ? data->scene : data->scenes_count ? &data->scenes[0] : nullptr;
        if (scene)
            for (size_t k = 0; k < scene->nodes_count && walk.problem.empty(); ++k) readNode(scene->nodes[k], walk);
        if (!walk.problem.empty()) return fail(err, path + ": " + walk.problem);
        raw = std::move(walk.prims);
    }
    if (!classifyPrimitives(raw, selection, path, err) || !scaleToLibrary(raw, selection, path, out, err)) return false;
    for (RawPrim& rp : raw) addPrimitive(rp, o, out);
    return true;
}

bool listGlbNodes(const std::string& path, const GlbSelection& selection, std::vector<GlbNode>& out, std::string* err) {
    out.clear();
    CgltfFile file;
    cgltf_options opt{};
    if (cgltf_parse_file(&opt, path.c_str(), &file.data) != cgltf_result_success) return fail(err, path + ": not a valid glTF / GLB file");
    const cgltf_data* data = file.data;
    const cgltf_scene* scene = data->scene ? data->scene : data->scenes_count ? &data->scenes[0] : nullptr;
    if (scene)
        for (size_t k = 0; k < scene->nodes_count; ++k)
            if (!listNode(scene->nodes[k], selection, path, out, err)) return false;
    return true;
}

} // namespace knit
