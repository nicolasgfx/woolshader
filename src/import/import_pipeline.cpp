// The preparation of an import, read top to bottom at the end of this file: prepareAndStore.
#include "import/import_pipeline.h"

#include <filesystem>

#include "model/coordinates.h"
#include "import/import_store.h"
#include "import/wrap_surface.h"
#include "io/cache_paths.h"
#include "io/files.h"

namespace knit {

namespace fs = std::filesystem;

namespace {

// The folder a rebuild works in. Its destructor removes it, however the job ends: the stored import keeps its own copies.
class RebuildFolder {
public:
    explicit RebuildFolder(const fs::path& path) : path_(path) {
        std::error_code ec;
        fs::remove_all(path_, ec);
        fs::create_directories(path_, ec);
    }
    RebuildFolder(const RebuildFolder&) = delete;
    RebuildFolder& operator=(const RebuildFolder&) = delete;
    ~RebuildFolder() {
        std::error_code ec;
        fs::remove_all(path_, ec);
    }
    const fs::path& path() const { return path_; }

private:
    fs::path path_;
};

// The rebuilt surface's three files (wrapped.obj, its MTL and texture) from this application's own stages (wrap_surface.h).
bool rebuildSurface(const ImportedMesh& in, const ImportRequest& r, const fs::path& folder, std::vector<std::string>& summary, std::string* err,
                    const ImportProgress& progress) {
    ImportedMesh source = in;  // loadImportTextures updates the materials
    std::vector<std::string> warnings;
    const std::vector<ColourImage> textures = loadImportTextures(source, r.limits, &warnings);
    return wrapFileToFolder(source, textures, r.options, folder.string(), summary, err, progress);
}

} // namespace

bool readRebuiltSurface(const fs::path& folder, const ImportOptions& o, const std::vector<std::string>& summary, PreparedMesh& out, std::string* err,
                        const ImportProgress& progress) {
    if (progress) progress("reading the rebuilt surface", 0);
    ImportedMesh surface;
    if (!readMeshFile((folder / (std::string(kRebuiltSurfaceName) + ".obj")).string(), ImportLimits{}, surface, err, {})) return false;
    ImportOptions asWritten;
    asWritten.unitMm = 1;
    asWritten.up = UpAxis::Y;
    if (!prepareImport(surface, asWritten, out, err)) return false;
    out.scale = o.unitMm;  // the way back to the source file's own coordinates (toFileCoordinates)
    out.up = o.up;
    out.chartsGenerated = true;
    out.changes = summary;
    return true;
}

namespace {

// The surface to knit, in library space.
bool prepareSurface(const ImportedMesh& in, const ImportRequest& r, const fs::path& folder, PreparedMesh& out, std::string* err, const ImportProgress& progress) {
    if (r.options.surface == ImportOptions::Surface::KeepFile) return prepareImport(in, r.options, out, err);
    std::vector<std::string> summary;
    return rebuildSurface(in, r, folder, summary, err, progress) && readRebuiltSurface(folder, r.options, summary, out, err, progress);
}

void moveToOffset(PreparedMesh& p, const Vec3& offset) {
    for (Vec3& q : p.mesh.pos) q = q + (offset - p.offset);
    p.offset = offset;
}

// The plain nodes as the file has them, then both meshes where the whole object stands, as its preview shows it. A rebuilt surface
// keeps the file's frame, so one offset lines the two up.
bool preparePlainNodes(const ImportedMesh& plainNodes, const ImportRequest& r, PreparedImport& out, std::string* err) {
    if (!prepareImport(plainNodes, r.options, out.plain, err)) return false;
    const Vec3 offset = toLibrarySpace(*r.file, r.options).offset;
    moveToOffset(out.knitted, offset);
    moveToOffset(out.plain, offset);
    return true;
}

// The triangles of the wanted nodes, with only the vertices they use.
ImportedMesh meshOfNodes(const ImportedMesh& file, const std::vector<char>& wanted) {
    ImportedMesh m;
    m.path = file.path;
    m.format = file.format;
    m.fileBytes = file.fileBytes;
    m.materials = file.materials;
    m.dependencies = file.dependencies;
    m.missingDependencies = file.missingDependencies;
    std::vector<int> newVertex(file.pos.size(), -1), newTriangle(file.tri.size() + 1, -1);
    for (size_t t = 0; t < file.tri.size(); ++t) {
        if (!wanted[file.nodeOfTriangle[t]]) continue;
        std::array<int, 3> corners;
        for (int k = 0; k < 3; ++k) {
            int& v = newVertex[file.tri[t][k]];
            if (v < 0) {
                v = static_cast<int>(m.pos.size());
                m.pos.push_back(file.pos[file.tri[t][k]]);
                if (!file.vcolor.empty()) m.vcolor.push_back(file.vcolor[file.tri[t][k]]);
            }
            corners[k] = v;
        }
        newTriangle[t] = static_cast<int>(m.tri.size());
        m.tri.push_back(corners);
        if (!file.uv.empty()) m.uv.push_back(file.uv[t]);
        if (t < file.triMaterial.size()) m.triMaterial.push_back(file.triMaterial[t]);
    }
    for (int first : file.quadFirstTriangle)
        if (newTriangle[first] >= 0 && newTriangle[first + 1] >= 0) m.quadFirstTriangle.push_back(newTriangle[first]);
    return m;
}

} // namespace

NodeSplit splitByNodes(const ImportedMesh& file, const std::set<std::string>& plainNodes) {
    std::vector<char> plain(file.nodes.size(), 0), knitted(file.nodes.size(), 1);
    for (size_t n = 0; n < file.nodes.size(); ++n) {
        plain[n] = plainNodes.count(file.nodes[n].name) ? 1 : 0;
        knitted[n] = !plain[n];
    }
    return {meshOfNodes(file, knitted), meshOfNodes(file, plain)};
}

bool prepareAndStore(const ImportRequest& r, std::vector<std::string>& changes, std::string* err, const ImportProgress& progress) {
    const std::string root = importCacheRoot();
    const bool rebuilt = r.options.surface == ImportOptions::Surface::Wrapped;
    if (rebuilt && importStored(root, r.identity)) {  // a rebuild costs minutes, and its identity fixes its result
        changes = storedImportChanges(root, r.identity);
        return true;
    }
    const NodeSplit split = r.options.plainNodes.empty() ? NodeSplit{*r.file, {}} : splitByNodes(*r.file, r.options.plainNodes);
    if (split.knitted.tri.empty()) return fail(err, "every component is left out: tick at least one to knit");
    // Outlives the store: the rebuilt surface's texture is copied from here.
    const RebuildFolder folder(fs::absolute(root) / (".rebuild_" + identityHex(r.identity)));
    PreparedImport prepared;
    if (!prepareSurface(split.knitted, r, folder.path(), prepared.knitted, err, progress)) return false;
    if (!split.plain.tri.empty() && !preparePlainNodes(split.plain, r, prepared, err)) return false;
    if (progress) progress("storing the prepared mesh", 0);
    if (!importStored(root, r.identity) && !storeImport(root, r, prepared, err)) return false;
    changes = prepared.knitted.changes;
    return true;
}

} // namespace knit
