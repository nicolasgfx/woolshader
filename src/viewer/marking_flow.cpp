// The Mark step's state and actions (marking_flow.h).
#include "viewer/marking_flow.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <filesystem>

#include "io/cache_paths.h"
#include "io/marking_file.h"
#include "io/stitch_mesh_obj.h"
#include "model/hash.h"
#include "viewer/view_content.h"
#include "viewer/workflow.h"

namespace fs = std::filesystem;

namespace viewer {

namespace {

// The overlay's colours per meaning (RGB), chosen to stay apart on a cream or grey surface.
constexpr uint8_t kMeaningColour[6][3] = {
    {40, 190, 80},    // start
    {240, 130, 30},   // finish
    {220, 50, 200},   // row-start path
    {235, 205, 40},   // boundary
    {225, 35, 35},    // cut
    {40, 190, 230},   // attachment
};
constexpr uint8_t kRingColour[3] = {250, 250, 250};
constexpr uint8_t kStitchEdgeColour[3] = {40, 40, 40};
constexpr uint8_t kProblemColour[3] = {255, 0, 0};
constexpr uint8_t kHighlightColour[3] = {255, 220, 0};   // a clicked finding: its row and places
constexpr uint8_t kShadowColour[3] = {120, 120, 120};    // in-use markings a previewed setup stands in for
constexpr uint8_t kHiddenAlpha = 70;  // what lies behind the surface: faint, still findable

int meaningIndex(knit::MarkingMeaning m) {
    switch (m) {
    case knit::MarkingMeaning::Start: return 0;
    case knit::MarkingMeaning::Finish: return 1;
    case knit::MarkingMeaning::RowStartPath: return 2;
    case knit::MarkingMeaning::Boundary: return 3;
    case knit::MarkingMeaning::Cut: return 4;
    case knit::MarkingMeaning::Attachment: return 5;
    }
    return 0;
}

void addVertex(std::vector<float>& pos, std::vector<uint8_t>& rgba, std::vector<uint8_t>& hidden, const knit::Vec3& p, const uint8_t c[3], uint8_t alpha) {
    pos.insert(pos.end(), {float(p.x), float(p.y), float(p.z)});
    rgba.insert(rgba.end(), {c[0], c[1], c[2], alpha});
    hidden.insert(hidden.end(), {c[0], c[1], c[2], kHiddenAlpha});
}

void addPolyline(OverlayGeometry& g, const std::vector<knit::Vec3>& line, const uint8_t c[3], uint8_t alpha) {
    for (size_t i = 1; i < line.size(); ++i) {
        addVertex(g.linePos, g.lineRgba, g.lineRgbaHidden, line[i - 1], c, alpha);
        addVertex(g.linePos, g.lineRgba, g.lineRgbaHidden, line[i], c, alpha);
    }
}

// A marking's colour: its meaning's; an offered proposal paler, a rejected one grey.
void markingColour(const knit::SurfaceMarking& m, uint8_t out[3]) {
    const uint8_t* c = kMeaningColour[meaningIndex(m.meaning)];
    for (int k = 0; k < 3; ++k) {
        out[k] = c[k];
        if (m.proposal == knit::ProposalState::Offered) out[k] = static_cast<uint8_t>((c[k] + 255) / 2);
        if (m.proposal == knit::ProposalState::Rejected) out[k] = 140;
    }
}

double nowSeconds() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// A finding's highlighted part, drawn on the trace of the review that made it (its trace ranges index that review's own
// traces, traced on the refined reference: the viewer's piece traces can have other point counts).
void addFindingRange(OverlayGeometry& g, const knit::MarkingFinding& f, const std::vector<int>& markingIds, const std::vector<knit::MarkingTrace>& traces,
                     const uint8_t colour[3]) {
    if (f.firstTracePoint < 0) return;
    for (size_t i = 0; i < markingIds.size() && i < traces.size(); ++i) {
        if (markingIds[i] != f.markingId) continue;
        const knit::MarkingTrace& tr = traces[i];
        const int last = std::min(f.lastTracePoint, static_cast<int>(tr.points.size()) - 1);
        for (int k = std::max(1, f.firstTracePoint + 1); k <= last; ++k) {
            addVertex(g.linePos, g.lineRgba, g.lineRgbaHidden, tr.points[k - 1].pos, colour, 255);
            addVertex(g.linePos, g.lineRgba, g.lineRgbaHidden, tr.points[k].pos, colour, 255);
        }
    }
}

void rebuildOverlay(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    OverlayGeometry g;
    g.revision = mk.overlay.revision + 1;
    const knit::MarkedPieceView* built = shownBuild(ui);
    if (built && mk.view == MarkView::Rows)
        for (const std::vector<knit::Vec3>& ring : built->rings) addPolyline(g, ring, kRingColour, 255);
    if (built && mk.view == MarkView::StitchMesh)
        for (const std::vector<knit::Vec3>& edge : built->stitchEdges) addPolyline(g, edge, kStitchEdgeColour, 255);
    // a preview: the setup's markings drawn in full colour, the in-use markings it stands in for greyed
    std::vector<int> shadowed;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    if (mk.previewSetup > 0 && snap.pieces) {
        const knit::SetupConflicts c = knit::setupConflicts(mk.set, *snap.pieces, mk.piece, mk.previewSetup);
        shadowed = c.acceptedProposals;
        shadowed.insert(shadowed.end(), c.userMarkings.begin(), c.userMarkings.end());
    }
    // two passes: offered and rejected proposals first, then the markings in use (and a previewed setup) on top of them
    for (int pass = 0; pass < 2; ++pass)
    for (size_t i = 0; i < mk.set.markings.size() && i < mk.traces.size(); ++i) {
        const knit::SurfaceMarking& m = mk.set.markings[i];
        if (mk.tracePiece[i] != mk.piece && mk.piece >= 0) continue;  // only the isolated piece's markings
        const bool inUse = m.proposal != knit::ProposalState::Offered && m.proposal != knit::ProposalState::Rejected;
        const bool front = inUse || (mk.previewSetup > 0 && m.setup == mk.previewSetup);
        if (front != (pass == 1)) continue;
        uint8_t c[3];
        markingColour(m, c);
        if (mk.previewSetup > 0 && m.setup == mk.previewSetup) std::copy(kMeaningColour[meaningIndex(m.meaning)], kMeaningColour[meaningIndex(m.meaning)] + 3, c);
        if (std::find(shadowed.begin(), shadowed.end(), m.id) != shadowed.end()) std::copy(kShadowColour, kShadowColour + 3, c);
        const knit::MarkingTrace& tr = mk.traces[i];
        for (size_t k = 1; k < tr.points.size(); ++k) {
            addVertex(g.linePos, g.lineRgba, g.lineRgbaHidden, tr.points[k - 1].pos, c, 255);
            addVertex(g.linePos, g.lineRgba, g.lineRgbaHidden, tr.points[k].pos, c, 255);
        }
        const bool active = m.id == mk.activeMarking;
        for (size_t k = 0; k < tr.controlAt.size(); ++k) {
            const uint8_t white[3] = {255, 255, 255}, accent[3] = {255, 230, 0};
            const uint8_t* colour = active && int(k) == mk.selectedPoint ? accent : active ? white : c;
            addVertex(g.pointPos, g.pointRgba, g.pointRgbaHidden, tr.points[tr.controlAt[k]].pos, colour, 255);
        }
    }
    // the places a suggestion's refusal names (where the rows split), for the isolated piece
    if (mk.proposals && mk.proposals->piece == mk.piece)
        for (const knit::MarkingFinding& f : mk.proposals->findings)
            for (const knit::Vec3& x : f.at) addVertex(g.pointPos, g.pointRgba, g.pointRgbaHidden, x, kProblemColour, 255);
    // problems of the current review, on its own traces; the clicked finding's range, row and places
    const std::shared_ptr<const knit::MarkingSession::QuickView> q = currentQuick(ui);
    const std::vector<knit::MarkingFinding>* findings = q ? &q->review.findings : built && !shownBuildOutdated(ui) ? &built->review.findings : nullptr;
    const std::vector<int>* ids = q ? &q->markingIds : built ? &built->markingIds : nullptr;
    const std::vector<knit::MarkingTrace>* traces = q ? &q->traces : built ? &built->traces : nullptr;
    if (findings)
        for (const knit::MarkingFinding& f : *findings)
            if (f.level == knit::FindingLevel::Problem) addFindingRange(g, f, *ids, *traces, kProblemColour);
    const knit::MarkingFinding* h = nullptr;
    const std::vector<int>* hIds = nullptr;
    const std::vector<knit::MarkingTrace>* hTraces = nullptr;
    if (mk.highlightFinding >= 0 && mk.highlightQuick && q && mk.highlightFinding < int(q->review.findings.size())) {
        h = &q->review.findings[mk.highlightFinding];
        hIds = &q->markingIds;
        hTraces = &q->traces;
    } else if (mk.highlightFinding >= 0 && !mk.highlightQuick && built && mk.highlightFinding < int(built->review.findings.size())) {
        h = &built->review.findings[mk.highlightFinding];
        hIds = &built->markingIds;
        hTraces = &built->traces;
    }
    if (h) {
        addFindingRange(g, *h, *hIds, *hTraces, kHighlightColour);
        if (built && h->row >= 0 && built->rings.size() == static_cast<size_t>(built->review.measures.rows) + 1 && h->row < int(built->rings.size()))
            addPolyline(g, built->rings[h->row], kHighlightColour, 255);
        for (const knit::Vec3& x : h->at) addVertex(g.pointPos, g.pointRgba, g.pointRgbaHidden, x, kHighlightColour, 255);
    }
    mk.overlay = std::move(g);
}

// Every marking's path and piece again (after loading, undo, a new file).
void retraceAll(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    mk.traces.assign(mk.set.markings.size(), knit::MarkingTrace{});
    mk.tracePiece.assign(mk.set.markings.size(), -1);
    if (snap.pieces) mk.markingsOnPiece.assign(snap.pieces->pieces.size(), 0);
    mk.stale = snap.pieces && mk.set.digest != snap.digest;  // made on another mesh: listed, never applied
    if (!snap.pieces || mk.stale) return;
    for (size_t i = 0; i < mk.set.markings.size(); ++i) {
        const knit::SurfaceMarking& m = mk.set.markings[i];
        const int piece = knit::pieceOfMarking(*snap.pieces, m);
        mk.tracePiece[i] = piece;
        if (piece < 0) continue;
        ++mk.markingsOnPiece[piece];
        const knit::TriMesh& mesh = snap.pieces->pieces[piece].mesh;
        std::vector<knit::SurfacePoint> local;
        for (const knit::SurfacePoint& p : m.points) {
            int pi = -1;
            knit::SurfacePoint lp;
            if (knit::pieceLocalPoint(*snap.pieces, p, pi, lp) && pi == piece) local.push_back(lp);
        }
        if (mk.neighbours.size() != snap.pieces->pieces.size()) mk.neighbours.assign(snap.pieces->pieces.size(), knit::TriangleNeighbours{});
        if (mk.neighbours[piece].across.empty()) mk.neighbours[piece] = knit::triangleNeighbours(mesh);  // once per piece
        // a closed start or termination on one border loop is that opening: drawn along the rim, as the mesher takes it
        std::vector<int> loop;
        const bool opening = m.shape == knit::MarkingShape::ClosedCurve && local.size() == m.points.size() &&
                             (knit::markingRole(m.meaning) == knit::MarkingRole::Start || knit::markingRole(m.meaning) == knit::MarkingRole::Termination) &&
                             knit::pointsOnOneBorderLoop(mesh, local, loop);
        mk.traces[i] = opening ? knit::borderLoopTrace(mesh, loop, local) : knit::traceMarking(mesh, mk.neighbours[piece], local, m.shape == knit::MarkingShape::ClosedCurve);
    }
}

// The guided Knit route's markings as a file of their own, named by the file's digest and the effective markings' digest
// (out/cache/markings/<stem>-<file digest>-guided-<markings digest>.markings): written once and never changed, so saved
// settings always find the markings they were knitted from, while the automatic save keeps following every edit.
// "" when it cannot be written.
std::string writeGuidedSnapshot(SceneUi& ui, const std::string& digest) {
    MarkingFlow& mk = ui.mark;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    const std::string stem = fs::path(snap.path).stem().string();
    const fs::path path = fs::absolute(fs::path(knit::markingsCacheRoot()) / (stem + "-" + knit::identityHex(snap.digest) + "-guided-" + digest + ".markings"));
    std::string err;
    fs::create_directories(path.parent_path());
    if (!fs::exists(path) && !knit::saveMarkingFile(path.string(), mk.set, &err)) {
        mk.note = "Not saved: " + err;
        return "";
    }
    return path.generic_string();
}

// The guided Knit route follows the markings: when an edit changes what its piece's build uses, the route points at a new
// snapshot with the new digest, so the knit made from the old markings shows as outdated. An edit that changes nothing
// the build uses (another piece, an offered or rejected proposal) leaves the route as it is.
void syncGuidedRoute(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    if (ui.s.knitRoute != "guided" || ui.s.guidedPiece < 0 || mk.stale) return;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    if (!snap.pieces || ui.s.guidedPiece >= int(snap.pieces->pieces.size()) || ui.s.importSource != snap.path) return;
    const std::string digest = knit::identityHex(knit::effectiveMarkingsDigest(mk.set, *snap.pieces, ui.s.guidedPiece, ui.s.guidedKnitComponent));
    if (digest == ui.s.guidedMarkingsDigest) return;
    const std::string path = writeGuidedSnapshot(ui, digest);
    if (path.empty()) return;
    ui.s.guidedMarkings = path;
    ui.s.guidedMarkingsDigest = digest;
}

void saveAutomatically(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    if (mk.stale || mk.savedTo.empty()) return;
    std::string err;
    fs::create_directories(fs::path(mk.savedTo).parent_path());
    if (!knit::saveMarkingFile(mk.savedTo, mk.set, &err)) {
        mk.note = "Not saved: " + err;
        return;
    }
    syncGuidedRoute(ui);
}

void requestQuickReview(SceneUi& ui);

// After an edit: undo kept by the caller, paths, overlay, the file, the quick review.
void markingsChanged(SceneUi& ui) {
    retraceAll(ui);
    ui.mark.highlightFinding = -1;
    rebuildOverlay(ui);
    saveAutomatically(ui);
    ui.mark.layoutDueAt = 0;
    requestQuickReview(ui);
}

// Every edit keeps the set for undo; an edit always applies to the saved intent, so it ends a preview first.
void keepForUndo(MarkingFlow& mk) {
    mk.previewSetup = 0;
    mk.undo.push_back(mk.set);
    if (mk.undo.size() > kUndoLevels) mk.undo.erase(mk.undo.begin());
    mk.redo.clear();
}

// A user's edit of a proposal makes it the user's marking (the proposal id stays, so the two can be compared).
void editedByUser(knit::SurfaceMarking& m) {
    if (m.origin != knit::MarkingOrigin::Proposal) return;
    if (m.proposalId == 0) m.proposalId = m.id;
    m.origin = knit::MarkingOrigin::User;
    m.proposal = knit::ProposalState::None;
}

std::string pieceLabel(const knit::MarkingSession::Snapshot& snap, int piece) {
    const knit::SurfacePiece& p = snap.pieces->pieces[piece];
    const std::string node = p.node >= 0 && p.node < int(snap.mesh->nodes.size()) ? snap.mesh->nodes[p.node].name : "component";
    return node + ", piece " + std::to_string(p.ordinal);
}

// The file was read: the pieces' labels, and the set saved for this file (or a new one).
void fileRead(SceneUi& ui, const knit::MarkingSession::Snapshot& snap) {
    MarkingFlow& mk = ui.mark;
    mk.pieceLabels.clear();
    for (size_t i = 0; i < snap.pieces->pieces.size(); ++i) mk.pieceLabels.push_back(pieceLabel(snap, int(i)));
    mk.savedTo = autosavePath(ui);
    mk.set = knit::MarkingSet{};
    std::string err;
    // Settings that knit this file on the guided route name their markings file: that set opens (copied into the automatic
    // save from then on, so the named file is never overwritten). Otherwise the automatic save of this file, if any.
    const bool guidedFile = ui.s.knitRoute == "guided" && !ui.s.guidedMarkings.empty() && ui.s.importSource == snap.path && fs::exists(ui.s.guidedMarkings);
    knit::MarkingSet named;
    if (guidedFile && knit::loadMarkingFile(ui.s.guidedMarkings, named, &err) && named.digest == snap.digest) {
        mk.set = named;
        mk.note = "Markings of the guided knit opened from " + fs::path(ui.s.guidedMarkings).filename().string();
    } else if (fs::exists(mk.savedTo) && knit::loadMarkingFile(mk.savedTo, mk.set, &err)) {
        mk.note = "Markings restored from " + fs::path(mk.savedTo).filename().string();
    }
    mk.set.source = snap.path;
    mk.set.digest = snap.digest;
    mk.undo.clear();
    mk.redo.clear();
    mk.neighbours.clear();
    mk.piece = -1;
    mk.activeMarking = 0;
    mk.selectedPoint = -1;
    retraceAll(ui);
    rebuildOverlay(ui);
}

// What a baseline is comparable on: the file, the reference preparation, the unit, the gauge, the piece and its knitting
// component.
uint64_t baselineKeyOf(const SceneUi& ui, uint64_t digest, int piece, int knitComponent) {
    uint64_t h = knit::fnv1a64Bytes(&digest, sizeof digest);
    h = knit::fnv1a64(ui.mark.set.preparation, h);
    const double numbers[3] = {ui.s.importUnitMm, ui.s.wale, ui.s.course};
    h = knit::fnv1a64Bytes(numbers, sizeof numbers, h);
    h = knit::fnv1a64Bytes(&piece, sizeof piece, h);
    return knit::fnv1a64Bytes(&knitComponent, sizeof knitComponent, h);
}

// New suggestions: merged into the set as one undo step (a refinement adds only its new setup).
void integrateProposals(SceneUi& ui, const knit::MarkingSession::Snapshot& snap) {
    MarkingFlow& mk = ui.mark;
    if (!snap.proposals || snap.proposals == mk.proposals || !snap.pieces || mk.stale) return;
    const knit::ProposalSet& ps = *snap.proposals;
    const bool refinement = mk.proposals && mk.proposals->piece == ps.piece && mk.proposals->basis == ps.basis && ps.setups.size() > mk.proposals->setups.size() &&
                            mk.cardSetup.size() == mk.proposals->setups.size();
    keepForUndo(mk);
    std::vector<knit::SurfaceMarking> offered;
    for (size_t i = refinement ? mk.proposals->setups.size() : 0; i < ps.setups.size(); ++i)
        for (const knit::SurfaceMarking& m : ps.setups[i].markings) {
            offered.push_back(m);
            offered.back().knitComponent = mk.knitComponent;   // suggestions go to the component being marked
        }
    std::vector<int> numbers;
    if (refinement) knit::appendProposals(mk.set, offered, &numbers);
    else knit::mergeProposals(mk.set, *snap.pieces, ps.piece, offered, &numbers);
    if (!refinement) mk.cardSetup.clear();
    for (int n : numbers) mk.cardSetup.push_back(n);
    mk.proposals = snap.proposals;
    mk.note = refinement ? "The refined start was added as suggestion " + std::to_string(ps.setups.size()) + "." : ps.shapeText + ": " + std::to_string(ps.setups.size()) + " suggestions.";
    retraceAll(ui);
    saveAutomatically(ui);
    requestQuickReview(ui);
}

// The arrived build of `view`'s piece and component in a list (null: none).
const MarkingFlow::ArrivedBuild* arrivedBuild(const std::vector<MarkingFlow::ArrivedBuild>& builds, const knit::MarkedPieceView& view) {
    for (const MarkingFlow::ArrivedBuild& a : builds)
        if (a.view.get() == &view) return &a;
    return nullptr;
}

// The session's builds taken over: for each one that arrived since the last frame, when it arrived; the explicit build
// selects the Stitch Mesh view when it answers the piece and component still shown with a mesh to draw; a build with a
// mesh may become the baseline of its own piece and component (a preview or a build without a mesh never does).
void integrateBuild(SceneUi& ui, const knit::MarkingSession::Snapshot& snap) {
    MarkingFlow& mk = ui.mark;
    std::vector<MarkingFlow::ArrivedBuild> builds;
    for (const std::shared_ptr<const knit::MarkedPieceView>& view : snap.builds) {
        if (const MarkingFlow::ArrivedBuild* known = arrivedBuild(mk.builds, *view)) {
            builds.push_back(*known);
            continue;
        }
        builds.push_back(MarkingFlow::ArrivedBuild{view, nowSeconds()});
        const bool shownPiece = view->piece == mk.piece && view->knitComponent == mk.knitComponent;
        if (shownPiece) mk.highlightFinding = -1;
        if (mk.explicitBuild != 0 && view->request == mk.explicitBuild) {
            mk.explicitBuild = 0;
            if (shownPiece && view->use.drawable) mk.view = MarkView::StitchMesh;
        }
        if (!view->previewIds.empty() || !view->use.drawable) continue;
        const uint64_t key = baselineKeyOf(ui, snap.digest, view->piece, view->knitComponent);
        if ((mk.baselineFromNextBuild && shownPiece) || !mk.baseline || mk.baselineKey != key) {
            mk.baseline = view;
            mk.baselineKey = key;
            mk.baselineAt = builds.back().at;
            if (shownPiece) mk.baselineFromNextBuild = false;
        }
    }
    mk.builds = builds;
}

void requestQuickReview(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    if (!snap.pieces || mk.stale || mk.piece < 0) return;
    std::vector<int> previewIds;
    const knit::MarkingSet set = effectiveMarkingSet(ui, previewIds);
    mk.quickInputs = knit::quickReviewInputs(set, *snap.pieces, mk.piece, ui.s.importUnitMm, ui.s.wale, ui.s.course, previewIds, mk.knitComponent);
    mk.quickFieldKey = snap.fieldKey;
    mk.quickFileGeneration = snap.fileGeneration;
    mk.quickSerial = mk.session.review(set, mk.piece, ui.s.importUnitMm, ui.s.wale, ui.s.course, previewIds, mk.knitComponent);
    using Phase = knit::MarkingSession::Phase;
    const bool background = snap.phase == Phase::Reading || snap.phase == Phase::Building || snap.phase == Phase::Proposing ||
                            snap.phase == Phase::Refining || ui.job.active;
    if (mk.quickSerial != 0) mk.pendingStamps.push_back({mk.quickSerial, nowSeconds(), background});
}

} // namespace

knit::MarkingSet effectiveMarkingSet(const SceneUi& ui, std::vector<int>& previewIds) {
    previewIds.clear();
    const MarkingFlow& mk = ui.mark;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    if (mk.previewSetup <= 0 || !snap.pieces) return mk.set;
    return knit::previewMarkings(mk.set, *snap.pieces, mk.piece, mk.previewSetup, previewIds);
}

std::shared_ptr<const knit::MarkingSession::QuickView> currentQuick(const SceneUi& ui) {
    const MarkingFlow& mk = ui.mark;
    const std::shared_ptr<const knit::MarkingSession::QuickView> q = mk.session.snapshot().quick;
    if (!q || mk.stale || q->serial != mk.quickSerial || q->inputs != mk.quickInputs || q->fileGeneration != mk.quickFileGeneration || q->piece != mk.piece) return nullptr;
    return q;
}

void suggestMarkings(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    if (mk.piece < 0 || mk.stale) return;
    mk.session.propose(mk.piece, ui.s.importUnitMm, ui.s.wale, ui.s.course, mk.set.preparation);
}

void refineSuggestion(SceneUi& ui, int card) {
    if (card <= 0 || !ui.mark.proposals || card > static_cast<int>(ui.mark.proposals->setups.size())) return;
    ui.mark.session.refine(card);
}

void cancelSuggestionWork(SceneUi& ui) { ui.mark.session.cancelDeep(); }

namespace {

int setupOfCard(const SceneUi& ui, int card) {
    const MarkingFlow& mk = ui.mark;
    if (card <= 0 || card > static_cast<int>(mk.cardSetup.size())) return 0;
    return mk.cardSetup[card - 1];
}

uint64_t currentBasis(const SceneUi& ui) {
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    const int piece = ui.mark.proposals ? ui.mark.proposals->piece : ui.mark.piece;
    return knit::proposalBasis(snap.digest, piece, ui.mark.set.preparation, ui.s.importUnitMm, ui.s.wale, ui.s.course);
}

} // namespace

void previewSuggestion(SceneUi& ui, int card) {
    MarkingFlow& mk = ui.mark;
    const int setup = setupOfCard(ui, card);
    mk.previewSetup = card <= 0 || setup == mk.previewSetup || suggestionState(ui, card) != "offered" ? 0 : setup;
    mk.highlightFinding = -1;
    rebuildOverlay(ui);
    requestQuickReview(ui);
}

bool acceptSuggestion(SceneUi& ui, int card, bool replaceUser) {
    MarkingFlow& mk = ui.mark;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    const int setup = setupOfCard(ui, card);
    if (!snap.pieces || setup <= 0 || mk.stale || !mk.proposals) return false;
    knit::MarkingSet trial = mk.set;
    if (!knit::acceptSetup(trial, *snap.pieces, mk.proposals->piece, setup, currentBasis(ui), replaceUser)) return false;
    keepForUndo(mk);
    mk.set = trial;
    mk.baselineFromNextBuild = true;
    markingsChanged(ui);
    return true;
}

void rejectSuggestion(SceneUi& ui, int card) {
    MarkingFlow& mk = ui.mark;
    const int setup = setupOfCard(ui, card);
    if (setup <= 0 || mk.stale) return;
    keepForUndo(mk);
    knit::rejectSetup(mk.set, setup);
    markingsChanged(ui);
}

std::string suggestionState(const SceneUi& ui, int card) {
    const MarkingFlow& mk = ui.mark;
    const int setup = setupOfCard(ui, card);
    if (setup <= 0 || !mk.proposals) return "";
    const size_t expected = mk.proposals->setups[card - 1].markings.size();
    size_t found = 0, accepted = 0, rejected = 0, outdated = 0;
    for (const knit::SurfaceMarking& m : mk.set.markings) {
        if (m.setup != setup) continue;
        ++found;
        if (m.origin == knit::MarkingOrigin::User) return "edited";
        accepted += m.proposal == knit::ProposalState::Accepted ? 1 : 0;
        rejected += m.proposal == knit::ProposalState::Rejected ? 1 : 0;
        outdated += m.proposal == knit::ProposalState::Offered && m.basis != currentBasis(ui) ? 1 : 0;
    }
    if (found != expected) return "edited";
    if (accepted == found) return "accepted";
    if (rejected == found) return "rejected";
    if (outdated > 0) return "outdated";
    return "offered";
}

std::string suggestionReplaces(const SceneUi& ui, int card) {
    const MarkingFlow& mk = ui.mark;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    const int setup = setupOfCard(ui, card);
    if (!snap.pieces || setup <= 0 || !mk.proposals) return "";
    std::string out;
    for (int id : knit::setupConflicts(mk.set, *snap.pieces, mk.proposals->piece, setup).userMarkings)
        if (const knit::SurfaceMarking* m = knit::findMarking(mk.set, id)) out += (out.empty() ? "" : ", ") + std::string(knit::meaningLabel(m->meaning));
    return out;
}

void setLayoutRefresh(SceneUi& ui, bool on) {
    ui.mark.layoutRefresh = on;
    ui.mark.layoutDueAt = 0;
}

std::string baselineComparison(const SceneUi& ui) {
    const MarkingFlow& mk = ui.mark;
    const knit::MarkedPieceView* b = shownBuild(ui);
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    if (!b || !mk.baseline || b == mk.baseline.get() || mk.baselineKey != baselineKeyOf(ui, snap.digest, b->piece, b->knitComponent)) return "";
    const knit::LayoutMeasures& now = b->review.measures;
    const knit::LayoutMeasures& was = mk.baseline->review.measures;
    auto problems = [](const knit::MarkingReview& r) {
        int n = 0;
        for (const knit::MarkingFinding& f : r.findings) n += f.level == knit::FindingLevel::Problem || f.level == knit::FindingLevel::Warning ? 1 : 0;
        return n;
    };
    char line[400];
    std::string t;
    std::snprintf(line, sizeof line, "Compared with the layout of the accepted markings (%s, built %.0f s ago):\n", mk.baseline->review.depth == knit::ReviewDepth::Relaxed ? "relaxed" : "unrelaxed",
                  nowSeconds() - mk.baselineAt);
    t += line;
    std::snprintf(line, sizeof line, "rows %d -> %d; rows not crossed once or split %d -> %d; cells without a construction %d -> %d\n", was.rows, now.rows, was.contourEvents, now.contourEvents,
                  was.unsupportedCells, now.unsupportedCells);
    t += line;
    std::snprintf(line, sizeof line, "visible surface left out %.0f -> %.0f mm2; hidden %.0f -> %.0f mm2; coverage %s -> %s\n", was.omittedVisibleMm2, now.omittedVisibleMm2, was.omittedHiddenMm2,
                  now.omittedHiddenMm2, was.coveragePass ? "PASS" : "FAIL", now.coveragePass ? "PASS" : "FAIL");
    t += line;
    const int fewer = problems(mk.baseline->review) - problems(b->review);
    if (fewer > 0 && now.omittedVisibleMm2 > was.omittedVisibleMm2 + 1) {
        std::snprintf(line, sizeof line, "%d warning%s fewer, but %.0f mm2 more visible surface left out: a trade, not an improvement.\n", fewer, fewer == 1 ? "" : "s",
                      now.omittedVisibleMm2 - was.omittedVisibleMm2);
        t += line;
    }
    return t;
}

namespace {

std::string knownText(knit::Known k, const std::string& value) {
    switch (k) {
    case knit::Known::No: return "not computed";
    case knit::Known::Outdated: return value + " (outdated)";
    case knit::Known::Current: return value;
    }
    return value;
}

} // namespace

std::string quickReviewText(const SceneUi& ui) {
    const std::shared_ptr<const knit::MarkingSession::QuickView> q = currentQuick(ui);
    if (!q) return "no current quick review (checking)\n";
    const knit::LayoutMeasures& m = q->review.measures;
    std::string t = "quick review of piece " + std::to_string(q->piece) + " (" + std::to_string(q->ms) + " ms, preparation " + std::to_string(q->preparedMs) + " ms)\n";
    t += "rows: " + knownText(m.rowsKnown, std::to_string(m.rows)) + "\n";
    t += "rows not crossed once or split: " + knownText(m.contourEventsKnown, std::to_string(m.contourEvents)) + "\n";
    t += "cells without a construction: " + knownText(m.unsupportedKnown, std::to_string(m.unsupportedCells)) + "\n";
    t += "coverage: " + knownText(m.coverageKnown, m.coveragePass ? "PASS" : "FAIL") + "\n";
    t += "layout checks: " + knownText(m.checksKnown, m.checksPass ? "PASS" : "FAIL") + "\n";
    for (const knit::MarkingFinding& f : q->review.findings) {
        const char* level = f.level == knit::FindingLevel::Ok ? "ok" : f.level == knit::FindingLevel::Warning ? "warning" : f.level == knit::FindingLevel::Problem ? "problem" : "incomplete";
        char line[64];
        std::snprintf(line, sizeof line, "  [%s] marking %d, row %d, %zu places: ", level, f.markingId, f.row, f.at.size());
        t += line + f.message + "\n";
    }
    return t;
}

namespace {

// One line of latency statistics: answered and coalesced counts, median, p95 and worst of the answered latencies.
std::string latencyLine(const char* label, std::vector<double> v, int coalesced) {
    if (v.empty()) return std::string(label) + "answered 0, coalesced " + std::to_string(coalesced) + ": no latency measured\n";
    std::sort(v.begin(), v.end());
    auto at = [&](double share) { return v.empty() ? 0.0 : v[static_cast<size_t>(share * double(v.size() - 1))]; };
    char line[300];
    std::snprintf(line, sizeof line, "%sanswered %zu, coalesced %d, median %.1f ms, p95 %.1f ms, worst %.1f ms\n", label, v.size(), coalesced, at(0.5),
                  at(0.95), v.empty() ? 0.0 : v.back());
    return line;
}

} // namespace

std::string latencyText(const SceneUi& ui) {
    const MarkingFlow& mk = ui.mark;
    std::vector<double> during;
    for (size_t i = 0; i < mk.latencyMs.size() && i < mk.latencyBackground.size(); ++i)
        if (mk.latencyBackground[i]) during.push_back(mk.latencyMs[i]);
    // all edits; then only those made while background work ran (the per-edit truth, whatever the burst's start)
    std::string t = latencyLine("", mk.latencyMs, mk.coalesced);
    t += latencyLine("while background work ran: ", during, mk.coalescedBackground);
    char line[64];
    for (size_t i = 0; i < mk.latencyMs.size(); ++i) {
        std::snprintf(line, sizeof line, "%.2f%s\n", mk.latencyMs[i], i < mk.latencyBackground.size() && mk.latencyBackground[i] ? " background" : "");
        t += line;
    }
    return t;
}

void queueBurst(SceneUi& ui, int id, int point, const std::vector<knit::Vec3>& positions) {
    for (const knit::Vec3& p : positions) ui.mark.burst.push_back({id, point, p});
}

void markFrameShown(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    const std::shared_ptr<const knit::MarkingSession::QuickView> q = currentQuick(ui);
    if (!q || q->serial == mk.drawnQuick) return;
    mk.drawnQuick = q->serial;
    const double now = nowSeconds();
    std::vector<MarkingFlow::EditStamp> waiting;
    for (const MarkingFlow::EditStamp& s : mk.pendingStamps) {
        if (s.serial == q->serial) {
            mk.latencyMs.push_back(1000 * (now - s.at));
            mk.latencyBackground.push_back(s.background ? 1 : 0);
        } else if (s.serial < q->serial) {   // superseded before an answer: not a response
            ++mk.coalesced;
            if (s.background) ++mk.coalescedBackground;
        }
        else waiting.push_back(s);
    }
    mk.pendingStamps = waiting;
}

void highlightFinding(SceneUi& ui, int index, bool quick) {
    ui.mark.highlightFinding = index;
    ui.mark.highlightQuick = quick;
    rebuildOverlay(ui);
}

std::string autosavePath(const SceneUi& ui) {
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    if (!snap.pieces) return "";
    const std::string stem = fs::path(snap.path).stem().string();
    return (fs::path(knit::markingsCacheRoot()) / (stem + "-" + knit::identityHex(snap.digest) + ".markings")).string();
}

const knit::MarkedPieceView* shownBuild(const SceneUi& ui) {
    for (const MarkingFlow::ArrivedBuild& a : ui.mark.builds)
        if (a.view->piece == ui.mark.piece && a.view->knitComponent == ui.mark.knitComponent) return a.view.get();
    return nullptr;
}

double shownBuildArrivedAt(const SceneUi& ui) {
    for (const MarkingFlow::ArrivedBuild& a : ui.mark.builds)
        if (a.view->piece == ui.mark.piece && a.view->knitComponent == ui.mark.knitComponent) return a.at;
    return 0;
}

bool shownBuildOutdated(const SceneUi& ui) {
    const knit::MarkedPieceView* built = shownBuild(ui);
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    if (!built || !snap.pieces) return false;
    std::vector<int> previewIds;
    const knit::MarkingSet set = effectiveMarkingSet(ui, previewIds);
    // either depth's build of today's markings is current (an automatic layout is unrelaxed and says so)
    for (const knit::ReviewDepth depth : {knit::ReviewDepth::Built, knit::ReviewDepth::Relaxed})
        if (built->inputs == knit::markBuildInputs(set, *snap.pieces, ui.mark.piece, ui.s.importUnitMm, ui.s.wale, ui.s.course, depth, previewIds, ui.mark.knitComponent))
            return false;
    return true;
}

void integrateMarking(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    const bool importScene = ui.s.scene == "import" && !ui.s.importSource.empty();
    if (ui.stage == Stage::Mark && importScene && mk.openedPath != ui.s.importSource) {
        mk.openedPath = ui.s.importSource;
        mk.session.open(ui.s.importSource);
    }
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    if (snap.generation != mk.generation) {
        const bool newFile = snap.pieces && (mk.pieceLabels.size() != snap.pieces->pieces.size() || mk.set.digest != snap.digest);
        mk.generation = snap.generation;
        if (newFile) fileRead(ui, snap);
        integrateProposals(ui, snap);
        integrateBuild(ui, snap);
        rebuildOverlay(ui);  // a build, proposals or a quick review arrived (or failed)
    }
    if (!snap.pieces || mk.stale || mk.piece < 0) return;
    // a scripted burst: one point move per frame, as a drag would
    if (!mk.burst.empty()) {
        const MarkingFlow::BurstMove b = mk.burst.front();
        mk.burst.erase(mk.burst.begin());
        moveMarkingPoint(ui, b.id, b.point, b.position);
    }
    // anything else the quick review depends on changed (the gauge, the reference, a new field): ask again
    std::vector<int> previewIds;
    const knit::MarkingSet set = effectiveMarkingSet(ui, previewIds);
    if (knit::quickReviewInputs(set, *snap.pieces, mk.piece, ui.s.importUnitMm, ui.s.wale, ui.s.course, previewIds, mk.knitComponent) != mk.quickInputs ||
        snap.fieldKey != mk.quickFieldKey || snap.fileGeneration != mk.quickFileGeneration)
        requestQuickReview(ui);
    // "Update rows after edits": the rows unknown for complete markings, 0.7 s after the newest quick answer, once per input
    const std::shared_ptr<const knit::MarkingSession::QuickView> q = currentQuick(ui);
    if (!mk.layoutRefresh || !q || snap.phase != knit::MarkingSession::Phase::Ready) return;
    if (shownBuild(ui) && !shownBuildOutdated(ui)) return;   // a build of today's markings is shown: its rows are known
    bool complete = q->review.measures.rowsKnown != knit::Known::Current;
    for (const knit::MarkingFinding& f : q->review.findings) complete = complete && f.level != knit::FindingLevel::Problem && f.level != knit::FindingLevel::Incomplete;
    const uint64_t inputs = knit::markBuildInputs(set, *snap.pieces, mk.piece, ui.s.importUnitMm, ui.s.wale, ui.s.course, knit::ReviewDepth::Built, previewIds, mk.knitComponent);
    if (!complete || inputs == mk.autoLayoutInputs) return;
    if (mk.layoutDueAt == 0) mk.layoutDueAt = nowSeconds() + 0.7;
    if (nowSeconds() < mk.layoutDueAt) return;
    mk.layoutDueAt = 0;
    mk.autoLayoutInputs = inputs;
    mk.session.build(set, mk.piece, ui.s.importUnitMm, ui.s.wale, ui.s.course, knit::ReviewDepth::Built, previewIds, mk.knitComponent);
}

void isolatePiece(SceneUi& ui, int piece) {
    MarkingFlow& mk = ui.mark;
    if (piece == mk.piece) return;
    mk.previewSetup = 0;
    mk.highlightFinding = -1;
    mk.piece = piece;
    mk.knitComponent = 0;   // another piece starts at its first knitting component
    mk.activeMarking = 0;
    mk.selectedPoint = -1;
    ui.camera.frame = piece >= 0;  // the view frames the isolated piece (viewContent's bounds)
    rebuildOverlay(ui);
}

void setKnitComponent(SceneUi& ui, int knitComponent) {
    MarkingFlow& mk = ui.mark;
    if (knitComponent < 0 || knitComponent == mk.knitComponent) return;
    mk.knitComponent = knitComponent;
    mk.previewSetup = 0;     // a previewed setup belongs to the component it was offered for
    mk.highlightFinding = -1;
    mk.activeMarking = 0;
    mk.selectedPoint = -1;
    rebuildOverlay(ui);
    requestQuickReview(ui);
}

int knitComponentsOnPiece(const SceneUi& ui) {
    const MarkingFlow& mk = ui.mark;
    int highest = mk.knitComponent;
    for (size_t i = 0; i < mk.set.markings.size() && i < mk.tracePiece.size(); ++i)
        if (mk.tracePiece[i] == mk.piece) highest = std::max(highest, mk.set.markings[i].knitComponent);
    return highest + 1;
}

void addMarking(SceneUi& ui, knit::MarkingMeaning meaning) {
    MarkingFlow& mk = ui.mark;
    if (mk.piece < 0 || mk.stale) return;
    keepForUndo(mk);
    knit::SurfaceMarking m;
    m.id = mk.set.nextId++;
    m.meaning = meaning;
    m.shape = knit::defaultShape(meaning);
    m.knitComponent = mk.knitComponent;
    mk.set.markings.push_back(m);
    mk.activeMarking = m.id;
    mk.selectedPoint = -1;
    markingsChanged(ui);
}

void selectMarking(SceneUi& ui, int id) {
    ui.mark.activeMarking = id;
    ui.mark.selectedPoint = -1;
    rebuildOverlay(ui);
}

void setMarkingMeaning(SceneUi& ui, int id, knit::MarkingMeaning meaning) {
    knit::SurfaceMarking* m = knit::findMarking(ui.mark.set, id);
    if (!m || m->meaning == meaning) return;
    keepForUndo(ui.mark);
    m = knit::findMarking(ui.mark.set, id);
    m->meaning = meaning;
    if (!knit::shapeAllowed(meaning, m->shape)) m->shape = knit::defaultShape(meaning);
    editedByUser(*m);
    markingsChanged(ui);
}

void setMarkingShape(SceneUi& ui, int id, knit::MarkingShape shape) {
    knit::SurfaceMarking* m = knit::findMarking(ui.mark.set, id);
    if (!m || m->shape == shape || !knit::shapeAllowed(m->meaning, shape)) return;
    keepForUndo(ui.mark);
    m = knit::findMarking(ui.mark.set, id);
    m->shape = shape;
    if (shape == knit::MarkingShape::Point && m->points.size() > 1) m->points.resize(1);
    editedByUser(*m);
    markingsChanged(ui);
}

void deleteMarking(SceneUi& ui, int id) {
    MarkingFlow& mk = ui.mark;
    if (!knit::findMarking(mk.set, id)) return;
    keepForUndo(mk);
    mk.set.markings.erase(std::remove_if(mk.set.markings.begin(), mk.set.markings.end(), [id](const knit::SurfaceMarking& m) { return m.id == id; }),
                          mk.set.markings.end());
    if (mk.activeMarking == id) mk.activeMarking = 0;
    mk.selectedPoint = -1;
    markingsChanged(ui);
}

void deleteSelectedPoint(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    knit::SurfaceMarking* m = knit::findMarking(mk.set, mk.activeMarking);
    if (!m || mk.selectedPoint < 0 || mk.selectedPoint >= int(m->points.size())) return;
    keepForUndo(mk);
    m = knit::findMarking(mk.set, mk.activeMarking);
    m->points.erase(m->points.begin() + mk.selectedPoint);
    editedByUser(*m);
    mk.selectedPoint = -1;
    markingsChanged(ui);
}

void acceptProposal(SceneUi& ui, int id) {
    knit::SurfaceMarking* m = knit::findMarking(ui.mark.set, id);
    if (!m || m->origin != knit::MarkingOrigin::Proposal) return;
    keepForUndo(ui.mark);
    knit::findMarking(ui.mark.set, id)->proposal = knit::ProposalState::Accepted;
    markingsChanged(ui);
}

bool moveMarkingPoint(SceneUi& ui, int id, int point, const knit::Vec3& filePosition) {
    MarkingFlow& mk = ui.mark;
    knit::SurfaceMarking* m = knit::findMarking(mk.set, id);
    if (!m || mk.stale || point < 0 || point >= int(m->points.size())) return false;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    if (!snap.pieces) return false;
    const int piece = knit::pieceOfMarking(*snap.pieces, *m);
    if (piece < 0) return false;
    const knit::SurfacePiece& pc = snap.pieces->pieces[piece];
    const knit::SurfacePoint local = knit::closestSurfacePoint(pc.mesh, filePosition);
    if (local.tri < 0) return false;
    keepForUndo(mk);
    m = knit::findMarking(mk.set, id);
    knit::SurfacePoint filePoint = local;
    filePoint.tri = pc.fileTriangle[local.tri];   // a piece keeps every triangle's corner order: the weights carry over
    m->points[point] = filePoint;
    editedByUser(*m);
    markingsChanged(ui);
    return true;
}

namespace {

// One measure for the text report: its value, or "not computed" when the build did not reach it.
std::string measureText(knit::Known known, const char* format, double value) {
    if (known == knit::Known::No) return "not computed";
    char b[64];
    std::snprintf(b, sizeof b, format, value);
    return known == knit::Known::Outdated ? std::string(b) + " (outdated)" : std::string(b);
}

} // namespace

std::string markReviewText(const SceneUi& ui) {
    const knit::MarkedPieceView* b = shownBuild(ui);
    if (!b) return "no build of the isolated piece\n";
    const knit::MarkingReview& r = b->review;
    const knit::LayoutMeasures& m = r.measures;
    char line[512];
    std::string t;
    // every measure as computed, or "not computed" (never a 0 for a measure the build did not reach)
    std::snprintf(line, sizeof line, "piece %d, stitch mesh built %s, checks pass %s\n", b->piece, r.valid ? "yes" : "no",
                  m.checksKnown == knit::Known::No ? "not computed" : m.checksPass ? "yes" : "no");
    t += line;
    // the Stitch Mesh view's headline and the mesh as drawn and saved (file units)
    t += "stitch mesh: " + knit::stitchMeshHeadline(b->use, b->quality) + (shownBuildOutdated(ui) ? " (outdated)" : "") + "\n";
    if (b->use.drawable) {
        const knit::StitchMeshQuality& q = b->quality;
        knit::Vec3 lo, hi;
        knit::stitchMeshBounds(b->stitch, lo, hi);
        std::snprintf(line, sizeof line, "stitch mesh faces %d (regular %d, increase %d, decrease %d, special %d, unsupported %d, pole start %d), vertices %zu, rows %d\n",
                      q.faces, q.regular, q.increase, q.decrease, q.special, q.unsupported, q.poleStart, b->stitch.pos.size(), q.rows);
        t += line;
        std::snprintf(line, sizeof line, "stitch mesh bounds (file units) %.9g %.9g %.9g .. %.9g %.9g %.9g\n", lo.x, lo.y, lo.z, hi.x, hi.y, hi.z);
        t += line;
        t += "faces per row";
        for (int n : q.facesPerRow) t += " " + std::to_string(n);
        t += "\n";
        std::snprintf(line, sizeof line, "fit: %d faces over %.3f mm, max %.3f mm; inverted %d, degenerate %d, untriangulable %d, intersecting pairs %d; connectivity %s\n",
                      q.facesOverFit, q.fitTargetMm, q.maxDeviationMm, q.inverted, q.degenerate, q.untriangulable, q.intersectingPairs, q.connectivityPass ? "pass" : "fail");
        t += line;
        if (!q.uncoveredMeasured) t += "coverage: uncovered share not measurable (every sample within the boundary band)\n";
        std::snprintf(line, sizeof line, "coverage: uncovered %.4f (limit %.4f), covered twice %d, area ratio %.4f, boundary %.3f mm (limit %.3f mm); relaxation %s\n",
                      q.uncoveredFraction, q.uncoveredLimit, q.coveredTwice, q.areaRatio, q.boundaryMm, q.boundaryLimitMm,
                      !q.relaxed ? "not requested" : q.relaxRan ? q.relaxStop.c_str() : "did not run");
        t += line;
        for (const std::string& f : knit::stitchMeshFindings(q)) t += "  finding: " + f + "\n";
        if (!knit::stitchMeshRelaxationNote(q).empty()) t += "  note: " + knit::stitchMeshRelaxationNote(q) + "\n";
    }
    t += "validity: contour events " + measureText(m.contourEventsKnown, "%.0f", m.contourEvents) + ", unsupported cells " +
         measureText(m.unsupportedKnown, "%.0f", m.unsupportedCells) + "\n";
    t += "layout: rows " + measureText(m.rowsKnown, "%.0f", m.rows) + ", max stitch-count jump " + measureText(m.jumpKnown, "%.0f", m.maxStitchCountJump) +
         ", exceptional cells " + measureText(m.exceptionalKnown, "%.0f", m.exceptionalCells) + ", distortion p95 " +
         measureText(m.distortionKnown, "%.4f", m.distortionP95) + "\n";
    t += "coverage: covered " + measureText(m.coverageKnown, "%.4f", m.coveredFraction) + ", pass " +
         (m.coverageKnown == knit::Known::No ? std::string("not computed") : m.coveragePass ? std::string("yes") : std::string("no")) + ", omitted " +
         measureText(m.omittedKnown, "%.1f", m.omittedAreaMm2) + " mm2 (hidden " + measureText(m.omittedKnown, "%.1f", m.omittedHiddenMm2) + ", visible " +
         measureText(m.omittedKnown, "%.1f", m.omittedVisibleMm2) + ")\n";
    t += "cost: " + measureText(m.secondsKnown, "%.2f", m.seconds) + " s\n";
    t += "reference: " + ui.mark.set.preparation + "\n";
    for (const knit::MarkingFinding& f : r.findings) {
        const char* cause = f.cause == knit::FindingCause::Marking               ? "marking"
                            : f.cause == knit::FindingCause::AlgorithmLimitation ? "algorithm limitation"
                            : f.cause == knit::FindingCause::SolverLimitation    ? "solver limitation"
                                                                                 : "input";
        const char* level = f.level == knit::FindingLevel::Ok        ? "ok"
                            : f.level == knit::FindingLevel::Warning ? "warning"
                            : f.level == knit::FindingLevel::Problem ? "problem"
                                                                     : "incomplete";
        std::snprintf(line, sizeof line, "  [%s, %s] marking %d: %s\n", level, cause, f.markingId, f.message.c_str());
        t += line;
    }
    return t;
}

std::string markInventoryText(const SceneUi& ui) {
    const knit::MarkingSession::Snapshot snap = ui.mark.session.snapshot();
    if (!snap.pieces || !snap.mesh) return "file not read\n";
    char line[512];
    std::string t;
    const double mm = ui.s.importUnitMm;
    std::snprintf(line, sizeof line, "file %s\nunit %.9g mm per file unit, up %s\npieces %zu\n", snap.path.c_str(), mm, ui.s.importUp.c_str(), snap.pieces->pieces.size());
    t += line;
    knit::Vec3 lo, hi;
    for (size_t i = 0; i < snap.pieces->pieces.size(); ++i) {
        knit::Vec3 a, b;
        snap.pieces->pieces[i].mesh.bounds(a, b);
        lo = i == 0 ? a : knit::componentMin(lo, a);
        hi = i == 0 ? b : knit::componentMax(hi, b);
        const std::string label = i < ui.mark.pieceLabels.size() ? ui.mark.pieceLabels[i] : "";
        std::snprintf(line, sizeof line, "piece %zu \"%s\" triangles %zu size %.1f x %.1f x %.1f mm\n", i, label.c_str(), snap.pieces->pieces[i].mesh.tri.size(), (b.x - a.x) * mm,
                      (b.y - a.y) * mm, (b.z - a.z) * mm);
        t += line;
    }
    std::snprintf(line, sizeof line, "bounds (file units) %.9g %.9g %.9g .. %.9g %.9g %.9g\nsize %.1f x %.1f x %.1f mm\n", lo.x, lo.y, lo.z, hi.x, hi.y, hi.z, (hi.x - lo.x) * mm,
                  (hi.y - lo.y) * mm, (hi.z - lo.z) * mm);
    t += line;
    return t;
}

void rejectProposal(SceneUi& ui, int id) {
    knit::SurfaceMarking* m = knit::findMarking(ui.mark.set, id);
    if (!m || m->origin != knit::MarkingOrigin::Proposal) return;
    keepForUndo(ui.mark);
    knit::findMarking(ui.mark.set, id)->proposal = knit::ProposalState::Rejected;
    markingsChanged(ui);
}

void setMarkingPreparation(SceneUi& ui, const std::string& preparation) {
    if ((preparation != "facets" && preparation != "loop_limit") || ui.mark.set.preparation == preparation) return;
    keepForUndo(ui.mark);
    ui.mark.set.preparation = preparation;
    markingsChanged(ui);
}

void undoMarking(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    if (mk.undo.empty()) return;
    mk.redo.push_back(mk.set);
    mk.set = mk.undo.back();
    mk.undo.pop_back();
    if (!knit::findMarking(mk.set, mk.activeMarking)) mk.activeMarking = 0;
    mk.selectedPoint = -1;
    markingsChanged(ui);
}

void redoMarking(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    if (mk.redo.empty()) return;
    mk.undo.push_back(mk.set);
    mk.set = mk.redo.back();
    mk.redo.pop_back();
    if (!knit::findMarking(mk.set, mk.activeMarking)) mk.activeMarking = 0;
    mk.selectedPoint = -1;
    markingsChanged(ui);
}

void saveMarkingsAs(SceneUi& ui, const std::string& path) {
    std::string err;
    ui.mark.note = knit::saveMarkingFile(path, ui.mark.set, &err) ? "Saved to " + path : "Not saved: " + err;
}

void loadMarkingsFrom(SceneUi& ui, const std::string& path) {
    MarkingFlow& mk = ui.mark;
    knit::MarkingSet loaded;
    std::string err;
    if (!knit::loadMarkingFile(path, loaded, &err)) {
        mk.note = "Not loaded: " + err;
        return;
    }
    keepForUndo(mk);
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    mk.set = loaded;
    mk.note = snap.pieces && loaded.digest != snap.digest ? "These markings were made on another mesh: shown in the list, not applied" : "Loaded " + fs::path(path).filename().string();
    mk.activeMarking = 0;
    mk.selectedPoint = -1;
    retraceAll(ui);
    rebuildOverlay(ui);
    if (!mk.stale) saveAutomatically(ui);
}

void setMarkView(SceneUi& ui, MarkView view) {
    ui.mark.view = view;
    rebuildOverlay(ui);
}

void setNeighbourLook(SceneUi& ui, NeighbourLook look) { ui.mark.around = look; }

void buildMarkedPiece(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    if (mk.piece < 0 || mk.stale) return;
    std::vector<int> previewIds;
    const knit::MarkingSet set = effectiveMarkingSet(ui, previewIds);
    mk.layoutDueAt = 0;
    mk.explicitBuild = mk.session.build(set, mk.piece, ui.s.importUnitMm, ui.s.wale, ui.s.course, mk.relax ? knit::ReviewDepth::Relaxed : knit::ReviewDepth::Built,
                                        previewIds, mk.knitComponent);
}

std::string stitchMeshFileName(const SceneUi& ui) {
    const MarkingFlow& mk = ui.mark;
    std::string piece = mk.piece >= 0 && mk.piece < static_cast<int>(mk.pieceLabels.size()) ? mk.pieceLabels[mk.piece] : "piece";
    if (mk.knitComponent > 0) piece += " component " + std::to_string(mk.knitComponent + 1);
    // the piece label as one word: "left, piece 0" -> "left_piece_0"
    std::string word;
    for (const char c : piece) {
        const bool keep = std::isalnum(static_cast<unsigned char>(c)) || c == '-';
        if (keep) word += c;
        else if (!word.empty() && word.back() != '_') word += '_';
    }
    return fs::path(mk.openedPath).stem().string() + "-" + word + "-stitch-mesh." + knit::kStitchMeshExtension;
}

bool exportStitchMesh(SceneUi& ui, const std::string& path) {
    MarkingFlow& mk = ui.mark;
    const knit::MarkedPieceView* b = shownBuild(ui);
    if (!b || !b->use.drawable) {
        const std::string reason = b ? b->use.reason : std::string("build the isolated piece first");
        mk.note = "No stitch mesh to save: " + reason + (!reason.empty() && reason.back() == '.' ? "" : ".");
        return false;
    }
    const std::string label = b->piece < static_cast<int>(mk.pieceLabels.size()) ? mk.pieceLabels[b->piece] : "piece " + std::to_string(b->piece);
    std::vector<std::string> comments;
    comments.push_back("stitch mesh of " + mk.openedPath);
    comments.push_back("piece " + std::to_string(b->piece) + " (" + label + "), knitting component " + std::to_string(b->knitComponent + 1));
    comments.push_back("coordinates: the imported file's units and axes (" + std::to_string(ui.s.importUnitMm) + " mm per unit)");
    comments.push_back(knit::stitchMeshHeadline(b->use, b->quality));
    for (const std::string& f : knit::stitchMeshFindings(b->quality)) comments.push_back("finding: " + f);
    if (!knit::stitchMeshRelaxationNote(b->quality).empty()) comments.push_back("note: " + knit::stitchMeshRelaxationNote(b->quality));
    if (shownBuildOutdated(ui)) comments.push_back("outdated: the markings, the gauge or a preview changed after this build");
    if (!b->previewIds.empty()) comments.push_back("a preview of a suggestion, not the saved markings");
    std::string err;
    if (!knit::writeStitchMeshObj(path, b->stitch, comments, label, &err)) {
        mk.note = "Stitch mesh not saved: " + err;
        return false;
    }
    mk.note = "Saved " + std::to_string(b->stitch.faces()) + " stitch faces to " + path;
    return true;
}

bool knitThisPiece(SceneUi& ui) {
    MarkingFlow& mk = ui.mark;
    mk.previewSetup = 0;   // the Knit input is the saved intent, never a preview
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    if (mk.piece < 0 || mk.stale || !snap.pieces || mk.savedTo.empty()) return false;
    saveAutomatically(ui);
    const std::string digest = knit::identityHex(knit::effectiveMarkingsDigest(mk.set, *snap.pieces, mk.piece, mk.knitComponent));
    const std::string path = writeGuidedSnapshot(ui, digest);
    if (path.empty()) return false;
    if (!knit::knitWorkflow(ui.s)) chooseWorkflow(ui, knit::Workflow::Knit);  // a guided piece is knitted: the Knit workflow
    ui.s.knitRoute = "guided";
    ui.s.guidedPiece = mk.piece;
    ui.s.guidedKnitComponent = mk.knitComponent;
    ui.s.guidedMarkings = path;
    ui.s.guidedMarkingsDigest = digest;
    ui.stage = Stage::Knit;
    return true;
}

namespace {

// A point of the active marking placed: it moves the selected point, replaces a point marking's point, or extends a path or curve.
void placePoint(SceneUi& ui, const knit::SurfacePoint& point) {
    MarkingFlow& mk = ui.mark;
    keepForUndo(mk);
    knit::SurfaceMarking* m = knit::findMarking(mk.set, mk.activeMarking);
    if (mk.selectedPoint >= 0 && mk.selectedPoint < int(m->points.size())) {
        m->points[mk.selectedPoint] = point;
        mk.selectedPoint = -1;
    } else if (m->shape == knit::MarkingShape::Point) {
        m->points.assign(1, point);
    } else {
        m->points.push_back(point);
    }
    editedByUser(*m);
    markingsChanged(ui);
}

} // namespace

void markClick(SceneUi& ui, const Camera& camera, const ViewPixel& pixel, ClickIntent intent) {
    MarkingFlow& mk = ui.mark;
    const std::optional<SurfaceHit> hit = surfaceUnderPixel(ui, camera, pixel);
    if (!hit) return;
    if (hit->piece != mk.piece) {  // a click on another piece isolates it
        isolatePiece(ui, hit->piece);
        return;
    }
    if (mk.view == MarkView::StitchMesh) return;   // the Stitch Mesh view inspects: a click on the piece never edits its markings
    knit::SurfaceMarking* m = knit::findMarking(mk.set, mk.activeMarking);
    if (!m || mk.stale) return;
    if (intent == ClickIntent::Select) {
        // the active marking's control point nearest to the click, within 4 % of the piece's size
        const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
        knit::Vec3 lo, hi;
        snap.pieces->pieces[mk.piece].mesh.bounds(lo, hi);
        const double reach = 0.04 * knit::norm(hi - lo);
        int best = -1;
        double bestDistance = reach;
        for (size_t i = 0; i < mk.set.markings.size(); ++i) {
            if (mk.set.markings[i].id != mk.activeMarking) continue;
            const knit::MarkingTrace& tr = mk.traces[i];
            for (size_t k = 0; k < tr.controlAt.size(); ++k) {
                const double d = knit::norm(tr.points[tr.controlAt[k]].pos - hit->position);
                if (d >= bestDistance) continue;
                best = int(k);
                bestDistance = d;
            }
        }
        mk.selectedPoint = best;
        rebuildOverlay(ui);
        return;
    }
    placePoint(ui, hit->point);
}

bool placeMarkingPoint(SceneUi& ui, const knit::Vec3& filePosition) {
    MarkingFlow& mk = ui.mark;
    const knit::MarkingSession::Snapshot snap = mk.session.snapshot();
    if (!knit::findMarking(mk.set, mk.activeMarking) || mk.stale || mk.piece < 0 || !snap.pieces) return false;
    const knit::SurfacePiece& pc = snap.pieces->pieces[mk.piece];
    knit::SurfacePoint point = knit::closestSurfacePoint(pc.mesh, filePosition);
    if (point.tri < 0) return false;
    point.tri = pc.fileTriangle[point.tri];  // a piece keeps every triangle's corner order: the weights carry over
    placePoint(ui, point);
    return true;
}

} // namespace viewer
