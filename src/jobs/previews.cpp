// The Fibres step's reference preview: a small sample the controls act on in about a second.
#include "jobs/previews.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "jobs/display_scene.h"
#include "model/choices.h"
#include "model/hash.h"
#include "knitting/tube/ply_profile.h"
#include "knitting/tube/yarn_surface.h"
#include "scene/fibres/fur_stage.h"
#include "scene/fibres/reference_samples.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"

namespace knit {

namespace {

void pushVertex(DisplayMesh& d, const Vec3& p, const Vec3& n) {
    d.pos.insert(d.pos.end(), {float(p.x), float(p.y), float(p.z)});
    d.nrm.insert(d.nrm.end(), {float(n.x), float(n.y), float(n.z)});
}

// Quads between consecutive rings of `sides` vertices starting at vertex `first`.
void connectRings(DisplayMesh& d, uint32_t first, int rings, int sides) {
    for (int k = 0; k + 1 < rings; ++k)
        for (int j = 0; j < sides; ++j) {
            const uint32_t a = first + uint32_t(k * sides + j), b = first + uint32_t(k * sides + (j + 1) % sides);
            const uint32_t c = a + uint32_t(sides), e = b + uint32_t(sides);
            d.idx.insert(d.idx.end(), {a, b, e, a, e, c});
        }
}

// Vertex normals from the triangles (area-weighted), for displaced geometry.
void smoothNormals(DisplayMesh& d, size_t firstVertex) {
    const size_t nv = d.pos.size() / 3;
    std::vector<Vec3> acc(nv - firstVertex);
    auto P = [&](uint32_t i) { return Vec3{d.pos[3 * i], d.pos[3 * i + 1], d.pos[3 * i + 2]}; };
    for (size_t t = 0; t + 2 < d.idx.size(); t += 3) {
        const uint32_t a = d.idx[t], b = d.idx[t + 1], c = d.idx[t + 2];
        if (a < firstVertex || b < firstVertex || c < firstVertex) continue;
        const Vec3 n = cross(P(b) - P(a), P(c) - P(a));
        for (uint32_t v : {a, b, c}) acc[v - firstVertex] += n;
    }
    for (size_t v = firstVertex; v < nv; ++v) {
        const Vec3 n = normalized(acc[v - firstVertex]);
        d.nrm[3 * v] = float(n.x), d.nrm[3 * v + 1] = float(n.y), d.nrm[3 * v + 2] = float(n.z);
    }
}

// The yarn tubes of a surface model, sampled every ~0.1 mm with 48 sides; with a profile, every
// point is moved radially by rho(phi) / r (the twisted-ply relief the maps encode).
void appendYarnSurface(const SurfaceModel& sm, const PlyShadingProfile* prof, DisplayMesh& d) {
    const int sides = 48;
    for (const YarnSurface& ys : sm.yarns) {
        if (ys.rings.size() < 2) continue;
        const size_t step = std::max<size_t>(1, static_cast<size_t>(std::lround(0.1 / std::max(1e-6, ys.ds))));
        const size_t first = d.pos.size() / 3;
        int rings = 0;
        for (size_t i = 0; i < ys.rings.size(); i += step) {
            const Ring& R = ys.rings[i];
            const Vec3 ea = R.n1 * std::cos(R.phi) + R.n2 * std::sin(R.phi), eb = R.n2 * std::cos(R.phi) - R.n1 * std::sin(R.phi);
            for (int j = 0; j < sides; ++j) {
                const double t = 2 * kPi * j / sides;
                Vec3 off = ea * (R.a * std::cos(t)) + eb * (R.b * std::sin(t));
                if (prof) off = off * (prof->eval(prof->phiAt(R.s, t + R.phi)) / prof->r);
                pushVertex(d, R.c + off, normalized(off));
            }
            ++rings;
        }
        connectRings(d, uint32_t(first), rings, sides);
        if (prof) smoothNormals(d, first);
    }
}

// Appends a mesh (or its triangles `tris`) with its vertex normals.
void appendMesh(DisplayMesh& d, const TriMesh& m, const std::vector<int>* tris = nullptr) {
    const uint32_t first = uint32_t(d.pos.size() / 3);
    for (size_t v = 0; v < m.pos.size(); ++v) pushVertex(d, m.pos[v], m.vnormal[v]);
    auto add = [&](int t) {
        for (int k = 0; k < 3; ++k) d.idx.push_back(first + uint32_t(m.tri[t][k]));
    };
    if (tris)
        for (int t : *tris) add(t);
    else
        for (size_t t = 0; t < m.tri.size(); ++t) add(int(t));
}

void setColour(DisplayMesh& d, const Vec3& c) { d.color[0] = float(c.x), d.color[1] = float(c.y), d.color[2] = float(c.z); }

// Lit tubes along the visible part of every fibre (a polyline's buried root segment is
// skipped), with Fibre::vertexRadius per vertex; `sides` around.
void appendStrandTubes(const FibreSet& fs, const Vec3& offset, int sides, DisplayMesh& out) {
    for (const Fibre& f : fs.fibres) {
        const std::vector<Vec3> pts = fibrePoints(fs, f);
        const size_t k0 = 1;  // the tube starts at the curve start: the root stays inside the yarn
        if (pts.size() < k0 + 2) continue;
        const uint32_t first = uint32_t(out.pos.size() / 3);
        Vec3 u;  // transported normal
        int rings = 0;
        for (size_t k = k0; k < pts.size(); ++k) {
            const Vec3 t = normalized(pts[std::min(k + 1, pts.size() - 1)] - pts[k > k0 ? k - 1 : k]);
            if (k == k0) u = normalized(cross(t, std::fabs(t.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}));
            u = normalized(u - t * dot(u, t));
            const Vec3 v = cross(t, u);
            const double r = vertexRadius(f, static_cast<int>(k));
            for (int j = 0; j < sides; ++j) {
                const double a = 2 * kPi * j / sides;
                const Vec3 n = u * std::cos(a) + v * std::sin(a);
                pushVertex(out, pts[k] + offset + n * r, n);
            }
            ++rings;
        }
        connectRings(out, first, rings, sides);
    }
}

bool woolPreview(const SceneSettings& s, FibrePreview& out, std::string* err) {
    std::optional<WoolSample> grown = growWoolSample(s, err);
    if (!grown) return false;
    WoolSample& sample = *grown;
    // The yarn, with the twisted-ply appearance as relief (single tube only).
    DisplayMesh yarn;
    yarn.name = "reference patch";
    yarn.kind = DisplayMesh::Kind::Yarn;
    appendYarnSurface(sample.surface, sample.relief.get(), yarn);
    out.meshes.push_back(std::move(yarn));
    out.fibres = static_cast<long long>(sample.along.fibres.size() + sample.stray.fibres.size());
    DisplayMesh d, r;
    d.name = "fibres along the yarn";
    d.kind = r.kind = DisplayMesh::Kind::Plain;
    setColour(d, kOrientedFibreColour);  // the two populations told apart (the export keeps the materials)
    appendStrandTubes(sample.along, {}, 5, d);
    r.name = "stray fibres";
    setColour(r, kStrayFibreColour);
    appendStrandTubes(sample.stray, {}, 5, r);
    out.meshes.push_back(std::move(d));
    out.meshes.push_back(std::move(r));
    char b[400];
    std::snprintf(b, sizeof b,
                  "Reference patch: fibres on a %.0f mm square of %s, D %.2f mm at %.1f x %.1f mm%s; %lld fibres (yellow: along "
                  "the yarn, blue: stray). Approximate: a flat patch (no hidden test, no export rounding).",
                  sample.square, s.pattern.c_str(), s.diameter, s.wale, s.course,
                  sample.plyTubes ? (", " + std::to_string(s.plies) + " ply tubes").c_str()
                                  : sample.relief ? (", " + std::to_string(sample.reliefPlies) + "-ply relief shown as geometry").c_str() : "",
                  out.fibres);
    out.note = b;
    return true;
}

bool furPreview(const SceneSettings& s, FibrePreview& out, std::string* err) {
    (void)err;
    const FurSample sample = growFurSample(s);
    DisplayMesh skin;
    skin.name = "reference sphere";
    skin.kind = DisplayMesh::Kind::Plain;
    setColour(skin, {0.45, 0.42, 0.40});
    appendMesh(skin, sample.sphere);
    DisplayMesh fur;
    fur.name = "fur";
    fur.kind = DisplayMesh::Kind::Plain;
    setColour(fur, s.pileMaterial.color);
    appendStrandTubes(sample.hairs, {}, 5, fur);
    out.meshes.push_back(std::move(skin));
    out.meshes.push_back(std::move(fur));
    out.fibres = static_cast<long long>(sample.hairs.fibres.size());
    char b[400];
    std::snprintf(b, sizeof b,
                  "Reference sphere, radius %.0f mm: %lld fibres%s, %lld stuck in the body (rejected). Approximate: no coverage mask, a "
                  "static groom (no animation).",
                  sample.radius, out.fibres,
                  sample.thinning < 1 ? (" (density thinned to " + std::to_string(int(std::lround(100 * sample.thinning))) + " % for the preview)").c_str() : "",
                  sample.stuck);
    out.note = b;
    return true;
}

} // namespace

std::vector<Vec2> plyCrossSection(int plies, PlyProfileModel model, double fill, int points) {
    PlyGeometry g;
    g.count = plies;
    g.fill = fill;
    PlyProfileParams p;
    p.model = model;
    const PlyShadingProfile prof = buildPlyShadingProfile(g, 1.0, p);
    std::vector<Vec2> out;
    for (int i = 0; i < points; ++i) {
        const double phi = 2 * kPi * i / points, rho = prof.eval(phi);
        out.push_back({rho * std::cos(phi), rho * std::sin(phi)});
    }
    return out;
}

uint64_t fibrePreviewKey(const SceneSettings& s) {
    std::stringstream in(writeSettings(s));
    std::string line, kept;
    while (std::getline(in, line))
        if (line.rfind("scene ", 0) != 0 && line.rfind("import.", 0) != 0 && line.rfind("display.", 0) != 0 &&
            line.rfind("fibres.preset ", 0) != 0)  // a label of the values (fibre_presets.h)
            kept += line + "\n";
    return fnv1a(kept.data(), kept.size());
}

bool buildFibrePreview(const SceneSettings& s, FibrePreview& out, std::string* err) {
    const auto t0 = std::chrono::steady_clock::now();
    out = FibrePreview{};
    out.key = fibrePreviewKey(s);
    out.workflow = workflowOf(s);
    const bool ok = out.workflow == Workflow::Fur ? furPreview(s, out, err) : woolPreview(s, out, err);
    if (!ok) return false;
    DisplayBounds box;  // assigned even when empty (no mesh): the preview has no other bounds
    for (const DisplayMesh& d : out.meshes) addToBounds(box, d);
    out.lo = box.lo;
    out.hi = box.hi;
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return true;
}

// --- the hair coverage overlay

std::shared_ptr<BuiltScene> furCoverageOverlay(const BuiltScene& b, const SceneSettings& s, double* covered, std::string* err) {
    auto out = std::make_shared<BuiltScene>(b);
    std::string maskPath;
    if (!furMaskPath(s, maskPath, err)) return nullptr;
    FurMask mask;
    if (!maskPath.empty() && !loadFurMask(maskPath, mask, err)) return nullptr;
    const Vec3 red{0.95, 0.12, 0.10};
    long long kept = 0, texels = 0;
    if (!mask.value.empty())
        for (auto& tex : out->albedo) {
            auto t = std::make_shared<DisplayTexture>(*tex);
            for (int y = 0; y < t->height; ++y)
                for (int x = 0; x < t->width; ++x) {
                    // Display textures are image space (row 0 at the top); the mask's uv is v up.
                    const double m = sampleFurMask(mask, {(x + 0.5) / t->width, 1.0 - (y + 0.5) / t->height});
                    ++texels;
                    if (m >= 0.5) {
                        ++kept;
                        continue;
                    }
                    uint8_t* p = &t->rgba[(size_t(y) * size_t(t->width) + size_t(x)) * 4];
                    p[0] = uint8_t(0.4 * p[0] + 0.6 * 255 * red.x);
                    p[1] = uint8_t(0.4 * p[1] + 0.6 * 255 * red.y);
                    p[2] = uint8_t(0.4 * p[2] + 0.6 * 255 * red.z);
                }
            tex = t;
        }
    const std::vector<std::string> excluded = furExcludedParts(s);
    for (DisplayMesh& d : out->meshes)
        if (std::find(excluded.begin(), excluded.end(), d.name) != excluded.end()) {
            d.albedo = -1;
            d.vcol.clear();
            d.color[0] = float(red.x), d.color[1] = float(red.y), d.color[2] = float(red.z);
        }
    if (covered) *covered = texels ? double(kept) / double(texels) : 1.0;
    return out;
}

// --- the workers

FibrePreviewWorker::FibrePreviewWorker() { worker_ = std::thread([this] { loop(); }); }

FibrePreviewWorker::~FibrePreviewWorker() {
    {
        std::lock_guard<std::mutex> lk(m_);
        stop_ = true;
        if (token_) token_->cancel();
    }
    cv_.notify_all();
    worker_.join();
}

void FibrePreviewWorker::request(const SceneSettings& s) {
    const uint64_t key = fibrePreviewKey(s);
    std::lock_guard<std::mutex> lk(m_);
    requestedKey_ = key;
    if (key == doneKey_ || key == pendingKey_ || (running_ && key == runningKey_)) return;
    pending_ = std::make_unique<SceneSettings>(s);
    pendingKey_ = key;
    if (running_ && token_) token_->cancel();  // latest wins
    cv_.notify_all();
}

std::shared_ptr<const FibrePreview> FibrePreviewWorker::latest() const {
    std::lock_guard<std::mutex> lk(m_);
    return latest_;
}

std::string FibrePreviewWorker::error() const {
    std::lock_guard<std::mutex> lk(m_);
    return error_;
}

uint64_t FibrePreviewWorker::requestedKey() const {
    std::lock_guard<std::mutex> lk(m_);
    return requestedKey_;
}

bool FibrePreviewWorker::busy() const {
    std::lock_guard<std::mutex> lk(m_);
    return running_ || pending_ != nullptr;
}

void FibrePreviewWorker::loop() {
    for (;;) {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return stop_ || pending_ != nullptr; });
        if (stop_) return;
        // Debounce: wait until the edits pause for 150 ms (a slider drag sends many).
        for (uint64_t k = pendingKey_;;) {
            cv_.wait_for(lk, std::chrono::milliseconds(150));
            if (stop_) return;
            if (pendingKey_ == k) break;
            k = pendingKey_;
        }
        std::unique_ptr<SceneSettings> s = std::move(pending_);
        runningKey_ = pendingKey_;
        pendingKey_ = 0;
        running_ = true;
        token_ = std::make_shared<CancelToken>();
        std::shared_ptr<CancelToken> token = token_;
        lk.unlock();
        auto p = std::make_shared<FibrePreview>();
        std::string err;
        bool ok = false, cancelled = false;
        try {
            ThreadCancelScope scope(token.get());
            ok = buildFibrePreview(*s, *p, &err);
        } catch (const Cancelled&) {
            cancelled = true;
        } catch (const std::exception& e) {
            err = e.what();
        }
        lk.lock();
        running_ = false;
        if (!cancelled) {
            if (ok) {
                latest_ = p;
                doneKey_ = p->key;
                error_.clear();
            } else {
                error_ = err.empty() ? "the reference could not be built" : err;
                doneKey_ = runningKey_;
            }
        }
    }
}

FlatReferences::FlatReferences() { worker_ = std::thread([this] { loop(); }); }

FlatReferences::~FlatReferences() {
    {
        std::lock_guard<std::mutex> lk(m_);
        stop_ = true;
    }
    cv_.notify_all();
    worker_.join();
}

std::string FlatReferences::key(const SceneSettings& s) {
    char b[160];
    std::snprintf(b, sizeof b, "%.17g %.17g %.17g %.17g %.17g", knitDiameter(s, s.diameter), s.wale, s.course, s.depthAmplitude,
                  s.clearance);
    return b;
}

std::shared_ptr<const std::vector<PatternReference>> FlatReferences::get(const SceneSettings& s) {
    std::lock_guard<std::mutex> lk(m_);
    auto it = done_.find(key(s));
    return it == done_.end() ? nullptr : it->second;
}

void FlatReferences::request(const SceneSettings& s) {
    const std::string k = key(s);
    {
        std::lock_guard<std::mutex> lk(m_);
        if (done_.count(k) || running_ == k || pending_ == k) return;
        pending_ = k;
        pendingSettings_ = s;  // the swatch's parameters are the key's values (flatPatternParams)
    }
    cv_.notify_all();
}

void FlatReferences::loop() {
    for (;;) {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return stop_ || !pending_.empty(); });
        if (stop_) return;
        running_ = pending_;
        const SceneSettings s = pendingSettings_;
        pending_.clear();
        lk.unlock();
        auto refs = std::make_shared<std::vector<PatternReference>>();
        for (StitchPattern p : {StitchPattern::Stockinette, StitchPattern::ReverseStockinette, StitchPattern::Garter,
                                      StitchPattern::Rib1x1, StitchPattern::Rib2x2}) {
            refs->push_back(flatPatternReference(s, p));
            std::lock_guard<std::mutex> l(m_);
            if (stop_) return;
        }
        lk.lock();
        done_[running_] = refs;
        running_.clear();
    }
}

} // namespace knit
