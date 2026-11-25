// The Import step's work on its thread: read the file, then prepare and store it.
#include "jobs/import_session.h"

#include <algorithm>
#include <filesystem>
#include <map>
#include <sstream>

#include "model/coordinates.h"
#include "import/import_pipeline.h"
#include "import/import_store.h"
#include "io/cache_paths.h"
#include "io/files.h"
#include "jobs/display_scene.h"

namespace knit {

struct ImportSession::Task {
    enum class Kind { Open, Prepare } kind = Kind::Open;
    std::string path;
    ImportLimits lim;
    ImportOptions options;
    std::shared_ptr<const ImportedMesh> mesh;  // Prepare: the analysed file
    uint64_t request = 0;
};

std::string repairsText(const ImportOptions& o) {
    std::string s;
    auto add = [&](bool on, const char* name) {
        if (on) s += (s.empty() ? "" : ",") + std::string(name);
    };
    add(o.weld, "weld");
    add(o.removeDegenerate, "degenerate");
    add(o.removeDuplicateFaces, "duplicates");
    add(o.removeUnreferenced, "unreferenced");
    add(o.fixOrientation, "orientation");
    return s;
}

void parseRepairs(const std::string& text, ImportOptions& o) {
    auto has = [&](const char* w) { return ("," + text + ",").find("," + std::string(w) + ",") != std::string::npos; };
    o.weld = has("weld");
    o.removeDegenerate = has("degenerate");
    o.removeDuplicateFaces = has("duplicates");
    o.removeUnreferenced = has("unreferenced");
    o.fixOrientation = has("orientation");
}

ImportOptions importOptionsFrom(const SceneSettings& s) {
    ImportOptions o;
    o.unitMm = s.importUnitMm;
    o.up = s.importUp == "z" ? UpAxis::Z : UpAxis::Y;
    parseRepairs(s.importRepairs, o);
    o.surface = s.importSurface == "wrapped" ? ImportOptions::Surface::Wrapped : ImportOptions::Surface::KeepFile;
    o.wrapStitchMm = s.importWrapStitchMm;
    o.plainNodes = parsePlainNodes(s.importPlainNodes);
    return o;
}

namespace {

// The file as it looks: one display mesh per node and material, named after the node, per corner (flat normals, colour uv,
// vertex colours), in file coordinates.
std::shared_ptr<const BuiltScene> makePreview(const ImportedMesh& m, const std::vector<ColourImage>& tex) {
    auto b = std::make_shared<BuiltScene>();
    std::map<std::string, int> texIndex;
    for (const ColourImage& t : tex) {
        texIndex[t.path] = static_cast<int>(b->albedo.size());
        b->albedo.push_back(makeDisplayTexture(t));
    }
    std::map<std::pair<int, int>, std::vector<int>> trianglesOf;  // by (node, material)
    for (size_t t = 0; t < m.tri.size(); ++t) trianglesOf[{m.nodeOfTriangle[t], t < m.triMaterial.size() ? m.triMaterial[t] : -1}].push_back(static_cast<int>(t));
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& [nodeAndMaterial, tris] : trianglesOf) {
        cancelCheckpoint();
        const int mi = nodeAndMaterial.second;
        DisplayMesh d;
        d.kind = DisplayMesh::Kind::Base;
        const ImportMaterial* mat = mi >= 0 && mi < static_cast<int>(m.materials.size()) ? &m.materials[mi] : nullptr;
        d.name = m.nodes[nodeAndMaterial.first].name;  // the viewer greys an unticked node by this name
        if (mat && !mat->texture.empty() && !m.uv.empty()) {
            auto it = texIndex.find(mat->texture);
            if (it != texIndex.end()) d.albedo = it->second;
        }
        const Vec3 kd = mat && mat->hasKd ? mat->kd : kDefaultYarnColour;
        d.color[0] = d.albedo >= 0 ? 1.0f : float(kd.x);
        d.color[1] = d.albedo >= 0 ? 1.0f : float(kd.y);
        d.color[2] = d.albedo >= 0 ? 1.0f : float(kd.z);
        for (int t : tris) {
            const std::array<int, 3>& tr = m.tri[t];
            Vec3 n = cross(m.pos[tr[1]] - m.pos[tr[0]], m.pos[tr[2]] - m.pos[tr[0]]);
            const double l = norm(n);
            n = l > 0 ? n * (1.0 / l) : Vec3{0, 0, 1};
            for (int k = 0; k < 3; ++k) {
                const Vec3& p = m.pos[tr[k]];
                d.pos.insert(d.pos.end(), {float(p.x), float(p.y), float(p.z)});
                d.nrm.insert(d.nrm.end(), {float(n.x), float(n.y), float(n.z)});
                if (d.albedo >= 0) d.albedoUv.insert(d.albedoUv.end(), {float(m.uv[t][k].x), float(1.0 - m.uv[t][k].y)});  // rows top first
                if (!m.vcolor.empty()) d.vcol.insert(d.vcol.end(), {float(m.vcolor[tr[k]].x), float(m.vcolor[tr[k]].y), float(m.vcolor[tr[k]].z)});
                d.idx.push_back(static_cast<uint32_t>(d.idx.size()));
                lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
                hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
            }
        }
        b->meshes.push_back(std::move(d));
    }
    b->lo = lo;
    b->hi = hi;
    return b;
}

// The file, read and looked at: mesh, analysis and preview in `out`. Failed with a message in *err.
ImportSession::Phase readAndAnalyse(const std::string& path, const ImportLimits& lim, ImportSession::Snapshot& out, std::string* err, const ImportProgress& progress) {
    auto mesh = std::make_shared<ImportedMesh>();
    if (!readMeshFile(path, lim, *mesh, err, progress)) return ImportSession::Phase::Failed;
    progress("checking the mesh", 0);
    finishNodes(*mesh);
    if (static_cast<int>(mesh->nodes.size()) > lim.maxNodes) {
        (void)fail(err, std::filesystem::path(path).filename().string() + " holds " + countText(mesh->nodes.size()) + " components; the import limit is " +
                      countText(lim.maxNodes) + ". Join small parts in a modelling program and open it again.");
        return ImportSession::Phase::Failed;
    }
    out.analysis = std::make_shared<ImportAnalysis>(analyzeImport(*mesh));
    progress("decoding textures", 0);
    const std::vector<ColourImage> tex = loadImportTextures(*mesh, lim, &out.textureWarnings);
    out.preview = makePreview(*mesh, tex);
    out.mesh = mesh;
    return ImportSession::Phase::Analysed;
}

} // namespace

ImportSession::ImportSession() { worker_ = std::thread([this] { loop(); }); }

ImportSession::~ImportSession() {
    {
        std::lock_guard<std::mutex> lk(m_);
        stop_ = true;
        pending_.reset();
        if (running_) running_->cancel();
    }
    cv_.notify_all();
    worker_.join();
}

void ImportSession::dropFile() {
    state_.mesh.reset();
    state_.analysis.reset();
    state_.preview.reset();
    state_.textureWarnings.clear();
    state_.error.clear();
    state_.stage.clear();
    state_.fraction = 0;
}

void ImportSession::open(const std::string& path, const ImportLimits& lim) {
    {
        std::lock_guard<std::mutex> lk(m_);
        pending_ = std::make_unique<Task>();
        pending_->kind = Task::Kind::Open;
        pending_->path = path;
        pending_->lim = lim;
        pending_->request = ++request_;
        limits_ = lim;
        if (running_) running_->cancel();
        dropFile();
        state_.phase = Phase::Reading;
        state_.path = path;
        state_.stage = "reading";
    }
    cv_.notify_all();
}

void ImportSession::prepare(const ImportOptions& o) {
    {
        std::lock_guard<std::mutex> lk(m_);
        if (!state_.mesh) return;
        pending_ = std::make_unique<Task>();
        pending_->kind = Task::Kind::Prepare;
        pending_->options = o;
        pending_->mesh = state_.mesh;
        pending_->lim = limits_;
        pending_->request = ++request_;
        if (running_) running_->cancel();
        state_.phase = Phase::Preparing;
        state_.stage = "preparing";
        state_.fraction = 0;
        state_.error.clear();
    }
    cv_.notify_all();
}

void ImportSession::cancel() {
    std::lock_guard<std::mutex> lk(m_);
    ++request_;  // the queued or running request publishes nothing any more
    pending_.reset();
    if (running_) running_->cancel();
    if (state_.phase == Phase::Preparing) {
        state_.phase = state_.mesh ? Phase::Analysed : Phase::Idle;  // back to the analysed file (the Import step stays)
    } else if (state_.phase == Phase::Reading) {
        state_.phase = Phase::Cancelled;
    } else {
        return;  // nothing was in progress
    }
    state_.stage.clear();
    state_.fraction = 0;
    ++state_.generation;
}

void ImportSession::reset() {
    std::lock_guard<std::mutex> lk(m_);
    ++request_;
    pending_.reset();
    if (running_) running_->cancel();
    state_.phase = Phase::Idle;
    dropFile();
    ++state_.generation;
}

bool ImportSession::busy() const {
    std::lock_guard<std::mutex> lk(m_);
    return running_ != nullptr || pending_ != nullptr;
}

void ImportSession::setTestHook(std::function<void(const char*)> hook) {
    std::lock_guard<std::mutex> lk(m_);
    hook_ = std::move(hook);
}

void ImportSession::loop() {
    for (;;) {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return stop_ || pending_; });
        if (stop_) return;
        const Task task = *pending_;
        pending_.reset();
        CancelToken token;
        running_ = &token;
        const std::function<void(const char*)> hook = hook_;
        lk.unlock();
        if (hook) hook("start");
        auto progress = [this, &task](const std::string& st, double f) {
            std::lock_guard<std::mutex> l(m_);
            if (task.request != request_) return;  // retired: a newer state stands
            state_.stage = st;
            state_.fraction = f;
        };
        Phase result = Phase::Failed;
        std::string err;
        Snapshot done;  // what this task publishes
        try {
            ThreadCancelScope scope(&token);
            if (task.kind == Task::Kind::Open) {
                result = readAndAnalyse(task.path, task.lim, done, &err, progress);
            } else {
                done.identity = importIdentity(*task.mesh, task.options);
                if (prepareAndStore({task.mesh, task.options, task.lim, done.identity}, done.changes, &err, progress)) result = Phase::Prepared;
            }
        } catch (const Cancelled&) {
            result = Phase::Cancelled;
        } catch (const std::exception& e) {
            err = e.what();
        }
        if (hook) hook("finish");
        lk.lock();
        running_ = nullptr;
        if (task.request != request_) continue;  // superseded, cancelled or reset: the newer state stands
        state_.phase = result;
        state_.error = result == Phase::Failed ? (err.empty() ? std::string("failed") : err) : std::string();
        if (task.kind == Task::Kind::Open && result == Phase::Analysed) {
            state_.mesh = done.mesh;
            state_.analysis = done.analysis;
            state_.preview = done.preview;
            state_.textureWarnings = done.textureWarnings;
        }
        if (task.kind == Task::Kind::Prepare && result == Phase::Prepared) {
            state_.identity = done.identity;
            state_.options = task.options;
            state_.changes = done.changes;
        }
        if (task.kind == Task::Kind::Prepare && result == Phase::Cancelled && state_.mesh) state_.phase = Phase::Analysed;
        ++state_.generation;
    }
}

ImportSession::Snapshot ImportSession::snapshot() const {
    std::lock_guard<std::mutex> lk(m_);
    return state_;
}

bool ImportSession::storedImportMatches(const ImportOptions& o, const std::string& storedIdentity) const {
    std::shared_ptr<const ImportedMesh> mesh;
    {
        std::lock_guard<std::mutex> lk(m_);
        mesh = state_.mesh;
    }
    const uint64_t id = mesh ? importIdentity(*mesh, o) : 0;
    uint64_t want = 0;
    parseIdentityHex(storedIdentity, want);  // a bad text leaves what it read: then no stored import matches
    return id == want && importStored(importCacheRoot(), id);
}

SceneSettings ImportSession::applyTo(const SceneSettings& base) const {
    std::lock_guard<std::mutex> lk(m_);
    SceneSettings s = base;
    s.scene = "import";
    s.importSource = state_.path;
    s.importIdentity = identityHex(state_.identity);
    s.importUnitMm = state_.options.unitMm;
    s.importUp = state_.options.up == UpAxis::Y ? "y" : "z";
    s.importRepairs = repairsText(state_.options);
    s.importSurface = state_.options.surface == ImportOptions::Surface::Wrapped ? "wrapped" : "file";
    s.importWrapStitchMm = state_.options.wrapStitchMm;
    s.importPlainNodes = plainNodesText(state_.options.plainNodes);
    return s;
}

} // namespace knit
