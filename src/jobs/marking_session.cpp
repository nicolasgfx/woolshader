// The Mark step's background work (marking_session.h).
#include "jobs/marking_session.h"

#include <algorithm>
#include <chrono>
#include <exception>

#include "import/import_analysis.h"
#include "model/hash.h"
#include "stitchmesh/mark_proposals.h"
#include "stitchmesh/marked_build.h"
#include "stitchmesh/quick_review.h"

namespace knit {

namespace {

// The display colours of the Mark step's own meshes (diagnostic colours, not materials).
constexpr float kPieceGrey[3] = {0.72f, 0.72f, 0.70f};
constexpr float kRegionColour[4][3] = {
    {0.86f, 0.82f, 0.70f},  // kept: the band Method A knits
    {0.42f, 0.56f, 0.80f},  // omitted beyond the start curve
    {0.80f, 0.56f, 0.40f},  // omitted beyond the finish curve
    {0.46f, 0.46f, 0.46f},  // outside the boundary
};
constexpr float kHiddenOmitted[3] = {0.36f, 0.38f, 0.42f};  // omitted, but inside another piece of the component: hidden
// Stitch faces by kind, in StitchKind's order. Regular faces are light blue so they stand apart from the beige piece and
// the grey neighbours.
constexpr float kKindColour[6][3] = {
    {0.62f, 0.78f, 0.95f},  // regular K
    {0.30f, 0.72f, 0.34f},  // increase
    {0.86f, 0.46f, 0.20f},  // decrease
    {0.62f, 0.42f, 0.84f},  // special quad
    {0.92f, 0.12f, 0.12f},  // unsupported
    {0.95f, 0.80f, 0.20f},  // pole start
};
// Stitch faces by what the checks found on them (the findings view).
constexpr float kNoFinding[3] = {0.84f, 0.87f, 0.90f};
constexpr float kOverFit[3] = {0.98f, 0.82f, 0.25f};       // over the fit target, up to twice it
constexpr float kFarOverFit[3] = {0.93f, 0.42f, 0.16f};    // more than twice the fit target
constexpr float kStructural[3] = {0.85f, 0.20f, 0.75f};    // inverted, degenerate or not triangulable
constexpr float kUnsupported[3] = {0.92f, 0.12f, 0.12f};

void addTriangle(DisplayMesh& d, const Vec3& a, const Vec3& b, const Vec3& c, const float colour[3]) {
    const Vec3 n = normalized(cross(b - a, c - a));
    for (const Vec3& p : {a, b, c}) {
        d.idx.push_back(static_cast<uint32_t>(d.pos.size() / 3));
        d.pos.insert(d.pos.end(), {float(p.x), float(p.y), float(p.z)});
        d.nrm.insert(d.nrm.end(), {float(n.x), float(n.y), float(n.z)});
        d.vcol.insert(d.vcol.end(), {colour[0], colour[1], colour[2]});
    }
}

void growBounds(BuiltScene& s, const TriMesh& m, bool first) {
    Vec3 lo, hi;
    m.bounds(lo, hi);
    s.lo = first ? lo : componentMin(s.lo, lo);
    s.hi = first ? hi : componentMax(s.hi, hi);
}

// One display mesh per piece, indexed, smooth normals, in file coordinates.
std::shared_ptr<const BuiltScene> piecesPreview(const SurfacePieces& pieces) {
    std::shared_ptr<BuiltScene> s = std::make_shared<BuiltScene>();
    for (size_t i = 0; i < pieces.pieces.size(); ++i) {
        const TriMesh& m = pieces.pieces[i].mesh;
        DisplayMesh d;
        d.name = "piece " + std::to_string(i);
        for (size_t v = 0; v < m.pos.size(); ++v) {
            d.pos.insert(d.pos.end(), {float(m.pos[v].x), float(m.pos[v].y), float(m.pos[v].z)});
            d.nrm.insert(d.nrm.end(), {float(m.vnormal[v].x), float(m.vnormal[v].y), float(m.vnormal[v].z)});
        }
        for (const std::array<int, 3>& t : m.tri) d.idx.insert(d.idx.end(), {uint32_t(t[0]), uint32_t(t[1]), uint32_t(t[2])});
        std::copy(kPieceGrey, kPieceGrey + 3, d.color);
        s->meshes.push_back(std::move(d));
        growBounds(*s, m, i == 0);
    }
    return s;
}

int regionIndex(MarkedRegion r) {
    switch (r) {
    case MarkedRegion::Kept: return 0;
    case MarkedRegion::OmittedStart: return 1;
    case MarkedRegion::OmittedFinish: return 2;
    case MarkedRegion::OutsideBoundary: return 3;
    }
    return 0;
}

StitchKind stitchKind(StitchFaceKind k) {
    switch (k) {
    case StitchFaceKind::Regular: return StitchKind::Regular;
    case StitchFaceKind::Increase: return StitchKind::Increase;
    case StitchFaceKind::Decrease: return StitchKind::Decrease;
    case StitchFaceKind::SpecialQuad: return StitchKind::SpecialQuad;
    case StitchFaceKind::Unsupported: return StitchKind::Unsupported;
    case StitchFaceKind::PoleStart: return StitchKind::PoleStart;
    }
    return StitchKind::Unsupported;
}

// The built polygons as model data, with the per-face findings of the final check (`check` measured `mesh` in mm, before
// it was scaled to file units, so its faces are the same).
StitchMeshData stitchMeshData(const TypedStitchMesh& mesh, const TypedMeshCheck& check) {
    StitchMeshData d;
    const PolygonMesh& m = mesh.mesh;
    d.pos = m.pos;
    d.faceStart.push_back(0);
    for (int f = 0; f < m.faces(); ++f) {
        for (int k = 0; k < m.sides(f); ++k) d.corner.push_back(m.vertexAt(f, k));
        d.faceStart.push_back(static_cast<int>(d.corner.size()));
        d.kind.push_back(f < static_cast<int>(mesh.kind.size()) ? stitchKind(mesh.kind[f]) : StitchKind::Unsupported);
        d.row.push_back(f < static_cast<int>(mesh.strip.size()) ? mesh.strip[f] : -1);
        uint8_t finding = 0;
        float deviation = 0;
        if (f < static_cast<int>(check.perFace.size())) {
            const FaceCheck& fc = check.perFace[f];
            finding |= fc.overDeviation ? kFaceOverFit : 0;
            finding |= fc.inverted ? kFaceInverted : 0;
            finding |= fc.degenerate ? kFaceDegenerate : 0;
            finding |= fc.untriangulable ? kFaceUntriangulable : 0;
            deviation = fc.deviation;
        }
        if (d.kind.back() == StitchKind::Unsupported) finding |= kFaceUnsupported;
        d.finding.push_back(finding);
        d.deviationMm.push_back(deviation);
    }
    if (d.faceStart.size() == 1) d.faceStart.clear();   // no faces
    return d;
}

// What the build's checks measured, copied from its evaluation (nothing measured again).
StitchMeshQuality stitchMeshQuality(const MarkedBuildProduct& p, const StitchMeshData& d, ReviewDepth depth) {
    const MarkedEvaluation& ev = p.evaluation;
    const TypedMeshCheck& c = ev.storedAfter;
    const CoverageReport& cov = ev.coverStoredAfter;
    StitchMeshQuality q;
    q.relaxed = depth == ReviewDepth::Relaxed;
    q.faces = d.faces();
    for (int f = 0; f < d.faces(); ++f) {
        switch (d.kind[f]) {
        case StitchKind::Regular: ++q.regular; break;
        case StitchKind::Increase: ++q.increase; break;
        case StitchKind::Decrease: ++q.decrease; break;
        case StitchKind::SpecialQuad: ++q.special; break;
        case StitchKind::Unsupported: ++q.unsupported; break;
        case StitchKind::PoleStart: ++q.poleStart; break;
        }
        if (d.row[f] < 0) continue;
        if (d.row[f] >= static_cast<int>(q.facesPerRow.size())) q.facesPerRow.resize(d.row[f] + 1, 0);
        ++q.facesPerRow[d.row[f]];
    }
    q.rows = p.review.measures.rowsKnown == Known::Current ? p.review.measures.rows : static_cast<int>(q.facesPerRow.size());
    q.fitTargetMm = ev.checkOptions.faceDeviation;
    q.maxDeviationMm = c.maxFaceDeviation;
    q.facesOverFit = c.facesOverDeviation;
    q.inverted = c.invertedFaces;
    q.degenerate = c.degenerateFaces;
    q.untriangulable = c.untriangulableFaces;
    q.intersectingPairs = c.intersectingPairs;
    q.connectivityPass = c.connectivity == Status::Pass;
    if (!q.connectivityPass) {
        const size_t open = c.detail.find_first_of("[{");
        q.connectivityDetail = open == std::string::npos ? c.detail : c.detail.substr(open + 1, c.detail.find_last_of("]}") - open - 1);
    }
    q.uncoveredFraction = cov.uncoveredFraction;
    q.uncoveredMeasured = cov.uncoveredMeasured || cov.samples == 0;   // no samples: no build reached the coverage
    q.uncoveredLimit = ev.coverOptions.uncoveredLimit;
    q.coveredTwice = cov.coveredTwice;
    q.areaRatio = cov.projectedRatio;
    q.areaLow = ev.coverOptions.areaLow;
    q.areaHigh = ev.coverOptions.areaHigh;
    q.boundaryMm = cov.worstHausdorff;
    q.boundaryLimitMm = ev.coverOptions.hausdorffLimit;
    q.coveragePass = cov.status == Status::Pass;
    q.relaxRan = ev.relaxed;
    q.relaxStop = ev.relaxed ? ev.relax.stopReason : "";
    for (const MarkingFinding& f : p.review.findings)
        if (f.cause == FindingCause::SolverLimitation && f.message.rfind("The relaxation did not run: ", 0) == 0) q.relaxNote = f.message.substr(28);
    q.clean = p.review.measures.checksKnown == Known::Current && p.review.measures.checksPass;
    return q;
}

// Why a build that returned made no stitch faces: the backend's own words (the first problem, else the first incomplete
// marking, else the first limitation, else its recommendation).
std::string noMeshReason(const MarkingReview& review) {
    for (const FindingLevel level : {FindingLevel::Problem, FindingLevel::Incomplete, FindingLevel::Warning})
        for (const MarkingFinding& f : review.findings)
            if (f.level == level) return f.message;
    if (!review.recommendation.empty()) return review.recommendation;
    return "the build made no stitch faces";
}

// A face's colour in the findings view: the most serious finding wins (structure, then a missing construction, then fit).
const float* findingColour(const StitchMeshData& d, int f, double fitTargetMm) {
    const uint8_t finding = d.finding[f];
    if (finding & kFaceStructural) return kStructural;
    if (finding & kFaceUnsupported) return kUnsupported;
    if (finding & kFaceOverFit) return d.deviationMm[f] > 2 * fitTargetMm ? kFarOverFit : kOverFit;
    return kNoFinding;
}

// Which build a view answers: its piece and knitting component, its request, its inputs, its unit and depth.
struct BuildIdentity {
    int piece = -1, knitComponent = 0;
    uint64_t request = 0, inputs = 0;
    double mmPerUnit = 1;
    ReviewDepth depth = ReviewDepth::Relaxed;
};

// The product as display data: the piece coloured by region (left-out surface another piece hides drawn darker); when the
// mesh can be drawn, the stitch faces by kind and by finding (lifted a little off the piece so they do not fight its
// depth) and the face outlines; the rings and the traces. `failure`: the build failed with this message (no product).
std::shared_ptr<const MarkedPieceView> pieceView(const SurfacePiece& piece, const BuildIdentity& id, const MarkedBuildProduct& p, const std::string& failure) {
    const ReviewDepth depth = id.depth;
    std::shared_ptr<MarkedPieceView> v = std::make_shared<MarkedPieceView>();
    v->piece = id.piece;
    v->knitComponent = id.knitComponent;
    v->request = id.request;
    v->inputs = id.inputs;
    v->review = p.review;
    v->review.depth = depth;
    v->markingIds = p.markingIds;
    v->previewIds = p.previewIds;
    v->traces = p.traces;
    v->rings = p.rings;
    v->log = p.log;
    v->stitch = stitchMeshData(p.stitchMesh, p.evaluation.storedAfter);
    v->quality = stitchMeshQuality(p, v->stitch, depth);
    v->use = stitchMeshUse(v->stitch);
    if (!failure.empty()) v->use = StitchMeshUse{false, "the build failed: " + failure};
    else if (!p.review.valid) v->use = StitchMeshUse{false, noMeshReason(p.review)};
    std::shared_ptr<BuiltScene> s = std::make_shared<BuiltScene>();
    DisplayMesh regions;
    regions.name = "marked piece";
    const TriMesh& m = piece.mesh;
    for (size_t t = 0; t < m.tri.size(); ++t) {
        const bool hidden = t < p.regionHidden.size() && p.regionHidden[t];
        const float* colour = p.region.empty() ? kPieceGrey : hidden ? kHiddenOmitted : kRegionColour[regionIndex(p.region[t])];
        addTriangle(regions, m.pos[m.tri[t][0]], m.pos[m.tri[t][1]], m.pos[m.tri[t][2]], colour);
    }
    s->meshes.push_back(std::move(regions));
    growBounds(*s, m, true);
    if (v->use.drawable) {
        const StitchMeshData& d = v->stitch;
        const double lift = 0.15 / id.mmPerUnit;  // 0.15 mm off the surface
        DisplayMesh byKind, byFinding;
        byKind.name = "stitch mesh";
        byFinding.name = "stitch mesh findings";
        for (int f = 0; f < d.faces(); ++f) {
            std::vector<Vec3> corners;
            for (int k = d.faceStart[f]; k < d.faceStart[f + 1]; ++k) corners.push_back(d.pos[d.corner[k]]);
            Vec3 normal;   // Newell's normal of the polygon
            for (size_t k = 0; k < corners.size(); ++k) normal = normal + cross(corners[k], corners[(k + 1) % corners.size()]);
            const Vec3 n = normalized(normal) * lift;
            std::vector<Vec3> outline;
            for (const Vec3& c : corners) outline.push_back(c + n);
            const float* kind = kKindColour[static_cast<int>(d.kind[f])];
            const float* found = findingColour(d, f, v->quality.fitTargetMm);
            for (size_t k = 1; k + 1 < outline.size(); ++k) {
                addTriangle(byKind, outline[0], outline[k], outline[k + 1], kind);
                addTriangle(byFinding, outline[0], outline[k], outline[k + 1], found);
            }
            outline.push_back(outline.front());
            v->stitchEdges.push_back(outline);
        }
        s->meshes.push_back(std::move(byKind));
        s->meshes.push_back(std::move(byFinding));
    }
    v->scene = s;
    return v;
}

// The entry of `piece`'s knitting component in a snapshot's builds (-1: none).
int buildIndex(const std::vector<std::shared_ptr<const MarkedPieceView>>& builds, int piece, int knitComponent) {
    for (size_t i = 0; i < builds.size(); ++i)
        if (builds[i]->piece == piece && builds[i]->knitComponent == knitComponent) return static_cast<int>(i);
    return -1;
}

uint64_t hashInts(const std::vector<int>& ids, uint64_t hash) { return ids.empty() ? hash : fnv1a64Bytes(ids.data(), ids.size() * sizeof(int), hash); }

} // namespace

uint64_t markBuildInputs(const MarkingSet& set, const SurfacePieces& pieces, int piece, double mmPerUnit, double w, double h, ReviewDepth depth,
                         const std::vector<int>& previewOffered, int knitComponent) {
    uint64_t hash = effectiveMarkingsDigest(set, pieces, piece, knitComponent);
    if (knitComponent != 0) hash = fnv1a64Bytes(&knitComponent, sizeof knitComponent, hash);   // component 0 hashes as before
    const double numbers[3] = {mmPerUnit, w, h};
    hash = fnv1a64Bytes(numbers, sizeof numbers, hash);
    const int ints[2] = {piece, static_cast<int>(depth)};
    hash = fnv1a64Bytes(ints, sizeof ints, hash);
    if (previewOffered.empty()) return hash;
    // a preview builds offered markings: their content counts too
    for (const SurfaceMarking& m : set.markings) {
        if (std::find(previewOffered.begin(), previewOffered.end(), m.id) == previewOffered.end()) continue;
        const int head[3] = {m.id, static_cast<int>(m.meaning), static_cast<int>(m.shape)};
        hash = fnv1a64Bytes(head, sizeof head, hash);
        for (const SurfacePoint& p : m.points) {
            hash = fnv1a64Bytes(&p.tri, sizeof p.tri, hash);
            hash = fnv1a64Bytes(p.bary, sizeof p.bary, hash);
        }
    }
    return hashInts(previewOffered, hash);
}

uint64_t quickReviewInputs(const MarkingSet& set, const SurfacePieces& pieces, int piece, double mmPerUnit, double w, double h, const std::vector<int>& previewOffered,
                           int knitComponent) {
    uint64_t hash = fnv1a64Bytes(&set.digest, sizeof set.digest);
    hash = fnv1a64(set.preparation, hash);
    const double numbers[3] = {mmPerUnit, w, h};
    hash = fnv1a64Bytes(numbers, sizeof numbers, hash);
    hash = fnv1a64Bytes(&piece, sizeof piece, hash);
    if (knitComponent != 0) hash = fnv1a64Bytes(&knitComponent, sizeof knitComponent, hash);   // component 0 hashes as before
    for (const SurfaceMarking& m : set.markings) {
        if (!markingInKnitComponent(pieces, m, piece, knitComponent)) continue;
        const int head[5] = {m.id, static_cast<int>(m.meaning), static_cast<int>(m.shape), static_cast<int>(m.origin), static_cast<int>(m.proposal)};
        hash = fnv1a64Bytes(head, sizeof head, hash);
        for (const std::pair<std::string, std::string>& kv : m.settings) hash = fnv1a64(kv.first + "=" + kv.second + ";", hash);
        for (const SurfacePoint& p : m.points) {
            hash = fnv1a64Bytes(&p.tri, sizeof p.tri, hash);
            hash = fnv1a64Bytes(p.bary, sizeof p.bary, hash);
        }
    }
    return hashInts(previewOffered, hash);
}

bool readMarkedFile(const std::string& path, MarkedFile& out, std::string* err) {
    SurfaceFile file;
    if (!readSurfaceFile(path, file, err)) return false;
    out.digest = file.digest;
    out.preview = piecesPreview(*file.pieces);
    out.mesh = file.mesh;
    out.pieces = file.pieces;
    return true;
}

// What the session keeps between requests, in stitchmesh's own types (never in the header).
struct MarkingSession::Kept {
    std::shared_ptr<const PreparedPiece> prepared;   // the quick lane's prepared piece, shared with builds of the same key
    std::shared_ptr<const MarkedField> field;        // the row field of the last build that reached the rows
    ProposalRequest proposalRequest;                 // the last proposals' request (refine repeats it)
};

struct MarkingSession::Task {
    enum class Kind { Open, Build, Propose, Refine } kind = Kind::Open;
    uint64_t request = 0;
    std::string path;
    MarkingSet set;
    int piece = -1;
    int knitComponent = 0;
    double mmPerUnit = 1, w = 5.5, h = 3.4;
    ReviewDepth depth = ReviewDepth::Relaxed;
    std::vector<int> previewOffered;
    std::string preparation = "facets";
    int setup = 0;
    uint64_t digest = 0;   // the open file's (proposals record it in their basis)
    std::shared_ptr<const SurfacePieces> pieces;
    std::shared_ptr<const PreparedPiece> prepared;
    std::shared_ptr<const ProposalSet> proposals;
    ProposalRequest proposalRequest;
};

struct MarkingSession::QuickTask {
    uint64_t serial = 0, inputs = 0;
    int fileGeneration = 0;
    MarkingSet set;
    int piece = -1;
    int knitComponent = 0;
    double mmPerUnit = 1, w = 5.5, h = 3.4;
    std::vector<int> previewOffered;
    std::shared_ptr<const SurfacePieces> pieces;
};

MarkingSession::MarkingSession() : kept_(std::make_unique<Kept>()) {
    worker_ = std::thread([this] { loop(); });
    quickWorker_ = std::thread([this] { quickLoop(); });
}

MarkingSession::~MarkingSession() {
    {
        std::lock_guard<std::mutex> lk(m_);
        stop_ = true;
        ++request_;
        ++quickSerial_;
        if (running_) running_->cancel();
    }
    cv_.notify_all();
    quickCv_.notify_all();
    worker_.join();
    quickWorker_.join();
}

void MarkingSession::open(const std::string& path) {
    std::lock_guard<std::mutex> lk(m_);
    if (running_) running_->cancel();
    pending_ = std::make_unique<Task>();
    pending_->kind = Task::Kind::Open;
    pending_->request = ++request_;
    pending_->path = path;
    quickPending_.reset();
    ++quickSerial_;   // a review of the previous file can no longer publish
    *kept_ = Kept{};
    const int generation = state_.generation + 1, fileGeneration = state_.fileGeneration + 1;
    state_ = Snapshot{};
    state_.phase = Phase::Reading;
    state_.path = path;
    state_.generation = generation;
    state_.fileGeneration = fileGeneration;
    cv_.notify_all();
}

uint64_t MarkingSession::build(const MarkingSet& set, int piece, double mmPerUnit, double w, double h, ReviewDepth depth, const std::vector<int>& previewOffered,
                               int knitComponent) {
    std::lock_guard<std::mutex> lk(m_);
    if (!state_.pieces) return 0;  // nothing read yet
    if (running_) running_->cancel();
    pending_ = std::make_unique<Task>();
    pending_->kind = Task::Kind::Build;
    pending_->request = ++request_;
    pending_->set = set;
    pending_->piece = piece;
    pending_->knitComponent = knitComponent;
    pending_->mmPerUnit = mmPerUnit;
    pending_->w = w;
    pending_->h = h;
    pending_->depth = depth;
    pending_->previewOffered = previewOffered;
    pending_->pieces = state_.pieces;
    pending_->prepared = kept_->prepared;
    state_.phase = Phase::Building;
    state_.error.clear();
    ++state_.generation;
    cv_.notify_all();
    return request_;
}

void MarkingSession::propose(int piece, double mmPerUnit, double w, double h, const std::string& preparation) {
    std::lock_guard<std::mutex> lk(m_);
    if (!state_.pieces || piece < 0 || piece >= static_cast<int>(state_.pieces->pieces.size())) return;
    if (running_) running_->cancel();
    pending_ = std::make_unique<Task>();
    pending_->kind = Task::Kind::Propose;
    pending_->request = ++request_;
    pending_->piece = piece;
    pending_->mmPerUnit = mmPerUnit;
    pending_->w = w;
    pending_->h = h;
    pending_->preparation = preparation;
    pending_->digest = state_.digest;
    pending_->pieces = state_.pieces;
    state_.phase = Phase::Proposing;
    state_.error.clear();
    state_.progress = 0;
    state_.progressText = "Suggesting markings ...";
    step_ = 0;
    stepsTotal_ = kProposeSteps;
    ++state_.generation;
    cv_.notify_all();
}

void MarkingSession::refine(int setup) {
    std::lock_guard<std::mutex> lk(m_);
    if (!state_.pieces || !state_.proposals || setup <= 0 || setup > static_cast<int>(state_.proposals->setups.size())) return;
    if (running_) running_->cancel();
    pending_ = std::make_unique<Task>();
    pending_->kind = Task::Kind::Refine;
    pending_->request = ++request_;
    pending_->setup = setup;
    pending_->pieces = state_.pieces;
    pending_->proposals = state_.proposals;
    pending_->proposalRequest = kept_->proposalRequest;
    state_.phase = Phase::Refining;
    state_.error.clear();
    state_.progress = 0;
    state_.progressText = "Refining the start (one trial layout with relaxation) ...";
    step_ = 0;
    stepsTotal_ = kRefineSteps;
    ++state_.generation;
    cv_.notify_all();
}

void MarkingSession::cancelDeep() {
    std::lock_guard<std::mutex> lk(m_);
    pending_.reset();
    ++request_;
    if (running_) running_->cancel();
    if (state_.phase == Phase::Building || state_.phase == Phase::Proposing || state_.phase == Phase::Refining) state_.phase = state_.pieces ? Phase::Ready : Phase::Idle;
    state_.progressText.clear();
    ++state_.generation;
}

uint64_t MarkingSession::review(const MarkingSet& set, int piece, double mmPerUnit, double w, double h, const std::vector<int>& previewOffered, int knitComponent) {
    std::lock_guard<std::mutex> lk(m_);
    if (!state_.pieces || piece < 0 || piece >= static_cast<int>(state_.pieces->pieces.size())) return 0;
    quickPending_ = std::make_unique<QuickTask>();
    quickPending_->serial = ++quickSerial_;
    quickPending_->fileGeneration = state_.fileGeneration;
    quickPending_->inputs = quickReviewInputs(set, *state_.pieces, piece, mmPerUnit, w, h, previewOffered, knitComponent);
    quickPending_->set = set;
    quickPending_->piece = piece;
    quickPending_->knitComponent = knitComponent;
    quickPending_->mmPerUnit = mmPerUnit;
    quickPending_->w = w;
    quickPending_->h = h;
    quickPending_->previewOffered = previewOffered;
    quickPending_->pieces = state_.pieces;
    state_.quickBusy = true;
    quickCv_.notify_all();
    return quickSerial_;
}

void MarkingSession::reset() {
    std::lock_guard<std::mutex> lk(m_);
    pending_.reset();
    quickPending_.reset();
    ++request_;
    ++quickSerial_;
    if (running_) running_->cancel();
    *kept_ = Kept{};
    const int generation = state_.generation + 1, fileGeneration = state_.fileGeneration + 1;
    state_ = Snapshot{};
    state_.generation = generation;
    state_.fileGeneration = fileGeneration;
}

namespace {

StitchLegendEntry legendEntry(const char* label, const float rgb[3]) {
    StitchLegendEntry e;
    e.label = label;
    for (int k = 0; k < 3; ++k) e.rgb[k] = rgb[k];
    return e;
}

} // namespace

std::vector<StitchLegendEntry> stitchKindLegend() {
    const char* labels[6] = {"Regular stitch", "Increase", "Decrease", "Special quad", "No construction", "Pole start"};
    std::vector<StitchLegendEntry> legend;
    for (int k = 0; k < 6; ++k) legend.push_back(legendEntry(labels[k], kKindColour[k]));
    return legend;
}

std::vector<StitchLegendEntry> stitchFindingLegend() {
    return {legendEntry("No finding", kNoFinding), legendEntry("Over the fit target", kOverFit), legendEntry("Over twice the fit target", kFarOverFit),
            legendEntry("Inverted, degenerate or not triangulable", kStructural), legendEntry("No construction", kUnsupported)};
}

std::shared_ptr<const MarkedPieceView> findBuild(const MarkingSession::Snapshot& snapshot, int piece, int knitComponent) {
    const int at = buildIndex(snapshot.builds, piece, knitComponent);
    return at < 0 ? nullptr : snapshot.builds[at];
}

MarkingSession::Snapshot MarkingSession::snapshot() const {
    std::lock_guard<std::mutex> lk(m_);
    Snapshot s = state_;
    if (s.phase == Phase::Proposing || s.phase == Phase::Refining) s.progress = std::min(1.0f, float(step_.load()) / float(std::max(1, stepsTotal_)));
    return s;
}

void MarkingSession::delayQuickForTest(int ms) {
    std::lock_guard<std::mutex> lk(m_);
    quickDelayMs_ = ms;
}

void MarkingSession::loop() {
    for (;;) {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return stop_ || pending_; });
        if (stop_) return;
        const Task task = *pending_;
        pending_.reset();
        std::shared_ptr<CancelToken> token = std::make_shared<CancelToken>();
        running_ = token;
        lk.unlock();
        std::string err;
        BuildIdentity id;   // a build's: known before it runs, so a failure can be published for its piece
        if (task.kind == Task::Kind::Build)
            id = BuildIdentity{task.piece, task.knitComponent, task.request,
                               markBuildInputs(task.set, *task.pieces, task.piece, task.mmPerUnit, task.w, task.h, task.depth, task.previewOffered, task.knitComponent),
                               task.mmPerUnit, task.depth};
        MarkedFile file;
        std::shared_ptr<const MarkedPieceView> view;
        std::shared_ptr<const MarkedField> field;
        std::shared_ptr<ProposalSet> proposals;
        ProposalRequest proposalRequest;
        bool ok = false, cancelled = false;
        try {
            ThreadCancelScope scope(token.get());
            if (task.kind == Task::Kind::Open) {
                ok = readMarkedFile(task.path, file, &err);
            } else if (task.kind == Task::Kind::Build) {
                MarkedBuildRequest rq;
                rq.pieces = task.pieces.get();
                rq.piece = task.piece;
                rq.knitComponent = task.knitComponent;
                rq.markings = &task.set;
                rq.mmPerUnit = task.mmPerUnit;
                rq.w = task.w;
                rq.h = task.h;
                rq.depth = task.depth;
                rq.previewOffered = task.previewOffered;
                rq.prepared = task.prepared;
                parsePreparation(task.set.preparation, rq.options.preparation);   // the file's choice (v3); facets otherwise
                MarkedBuildProduct product;
                ok = buildFromMarkings(rq, product, &err);
                if (ok) field = product.field;
                // a failed build publishes its piece's view too, without a mesh: the earlier mesh is never shown for it
                view = pieceView(task.pieces->pieces[task.piece], id, ok ? product : MarkedBuildProduct{}, ok ? std::string() : (err.empty() ? std::string("unknown error") : err));
            } else if (task.kind == Task::Kind::Propose) {
                proposalRequest.pieces = task.pieces.get();
                proposalRequest.piece = task.piece;
                proposalRequest.fileDigest = task.digest;
                proposalRequest.mmPerUnit = task.mmPerUnit;
                proposalRequest.w = task.w;
                proposalRequest.h = task.h;
                parsePreparation(task.preparation, proposalRequest.preparation);
                proposalRequest.step = &step_;
                proposals = std::make_shared<ProposalSet>();
                ok = proposeMarkings(proposalRequest, *proposals, &err);
            } else {
                ProposalRequest rq = task.proposalRequest;
                rq.pieces = task.pieces.get();
                rq.step = &step_;
                ProposalSetup refined;
                ok = refineProposalStart(rq, task.proposals->setups[task.setup - 1], refined, &err);
                if (ok) {
                    proposals = std::make_shared<ProposalSet>(*task.proposals);
                    refined.setup = static_cast<int>(proposals->setups.size()) + 1;
                    for (SurfaceMarking& m : refined.markings) {
                        m.setup = refined.setup;
                        m.basis = proposals->basis;
                    }
                    proposals->setups.push_back(refined);
                    proposalRequest = task.proposalRequest;
                }
            }
        } catch (const Cancelled&) {
            cancelled = true;
        } catch (const std::exception& e) {
            err = e.what();
            ok = false;
        }
        if (task.kind == Task::Kind::Build && !cancelled && !view) view = pieceView(task.pieces->pieces[task.piece], id, MarkedBuildProduct{}, err.empty() ? std::string("unknown error") : err);
        lk.lock();
        running_.reset();
        if (cancelled) ++state_.cancelledTasks;
        if (cancelled || task.request != request_) continue;  // superseded, cancelled or reset: the newer state stands
        ++state_.generation;
        state_.progressText.clear();
        if (task.kind == Task::Kind::Build) {   // also a failed build: it replaces its piece's entry
            state_.built = view;
            const int at = buildIndex(state_.builds, view->piece, view->knitComponent);
            if (at >= 0) state_.builds[at] = view;
            else state_.builds.push_back(view);
        }
        if (!ok) {
            state_.phase = Phase::Failed;
            state_.error = err.empty() ? std::string("failed") : err;
            if (task.kind != Task::Kind::Open && state_.pieces) state_.phase = Phase::Ready;  // the file stays usable
            continue;
        }
        state_.phase = Phase::Ready;
        state_.error.clear();
        if (task.kind == Task::Kind::Open) {
            state_.mesh = file.mesh;
            state_.pieces = file.pieces;
            state_.preview = file.preview;
            state_.digest = file.digest;
            state_.built.reset();
            state_.builds.clear();
        } else if (task.kind == Task::Kind::Build) {
            if (field) {
                kept_->field = field;
                state_.fieldKey = field->key;
                state_.fieldPiece = task.piece;
            }
        } else {
            state_.proposals = proposals;
            kept_->proposalRequest = proposalRequest;
            kept_->proposalRequest.step = nullptr;
        }
    }
}

void MarkingSession::quickLoop() {
    for (;;) {
        std::unique_lock<std::mutex> lk(m_);
        quickCv_.wait(lk, [&] { return stop_ || quickPending_; });
        if (stop_) return;
        const QuickTask task = *quickPending_;
        quickPending_.reset();
        const uint64_t digest = state_.digest;
        std::shared_ptr<const PreparedPiece> prepared = kept_->prepared;
        const std::shared_ptr<const MarkedField> field = kept_->field;
        const int delay = quickDelayMs_;
        ReferencePreparation preparation = ReferencePreparation::Facets;
        parsePreparation(task.set.preparation, preparation);
        const uint64_t key = preparedPieceKey(digest, task.piece, task.mmPerUnit, task.h, preparation);
        const bool cold = !prepared || prepared->key != key || prepared->piece != task.piece;
        if (cold) state_.quickPreparing = true;
        lk.unlock();
        std::shared_ptr<QuickView> view = std::make_shared<QuickView>();
        try {
            if (cold) {
                std::shared_ptr<PreparedPiece> made = std::make_shared<PreparedPiece>();
                preparePiece(*task.pieces, task.piece, digest, task.mmPerUnit, task.h, preparation, *made);
                view->preparedMs = made->seconds * 1000;
                prepared = made;
            }
            const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
            QuickRequest q;
            q.reference = prepared;
            q.field = field;
            q.markings = &task.set;
            q.pieces = task.pieces.get();
            q.piece = task.piece;
            q.knitComponent = task.knitComponent;
            q.mmPerUnit = task.mmPerUnit;
            q.w = task.w;
            q.h = task.h;
            q.preparation = preparation;
            q.previewOffered = task.previewOffered;
            view->review = quickReview(q, &view->markingIds, &view->traces);
            view->ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        } catch (const std::exception&) {
            view.reset();
        }
        if (delay > 0) std::this_thread::sleep_for(std::chrono::milliseconds(delay));
        lk.lock();
        if (cold && prepared && state_.fileGeneration == task.fileGeneration) kept_->prepared = prepared;
        state_.quickPreparing = false;
        state_.quickBusy = quickPending_ != nullptr;
        // only the newest request of this file publishes
        if (!view || task.serial != quickSerial_ || task.fileGeneration != state_.fileGeneration) continue;
        view->serial = task.serial;
        view->fileGeneration = task.fileGeneration;
        view->inputs = task.inputs;
        view->piece = task.piece;
        state_.quick = view;
        ++state_.generation;
    }
}

} // namespace knit
