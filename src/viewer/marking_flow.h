#pragma once
// The Mark step's state and what changes it (no ImGui, GL or GLFW).
// The step reads the imported file again through its own session (jobs/marking_session.h),
// lists its pieces, and lets the user isolate one and draw markings on it by clicking points. Every edit keeps the
// previous set for undo, re-traces the edited marking, rebuilds the overlay and saves the set next to the caches;
// "Build stitch mesh" hands the set to Method A on the session's worker (at the Knit settings' gauge) and shows its review;
// "Knit this piece" makes the isolated piece the guided Knit route's (knit.route = guided) with these markings.
//
// Assisted stockinette: "Suggest markings" asks the session for setups of the
// isolated piece (offered markings, merged into the set as one undo step); a setup is previewed on an ephemeral set in
// which it shadows the in-use markings of its roles (never saved, never the Knit input) and accepted as one undo step.
// Every change requests a quick review; one is shown as current only when it answers the newest request. After an edit
// that leaves the rows unknown, an unrelaxed layout is built 0.7 s later ("Update rows after edits").
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "jobs/marking_session.h"
#include "model/surface_marking.h"
#include "model/surface_path.h"
#include "viewer/camera.h"

namespace viewer {

struct SceneUi;

enum class NeighbourLook { Transparent, Hidden, Shown };  // the pieces around the isolated one
enum class MarkView { Markings, Rows, StitchMesh };        // what the viewport shows: markings, rows, or the Stitch Mesh view
enum class ClickIntent { Place, Select };                  // a click: place or move a point; Ctrl+click: select a point

// Lines and points drawn over the viewport's meshes, in file coordinates (the import preview's). Built when the
// markings or the build change (revision), never per frame; hidden parts are drawn faintly (the *Hidden colours).
struct OverlayGeometry {
    std::vector<float> linePos;           // GL_LINES pairs, xyz
    std::vector<uint8_t> lineRgba, lineRgbaHidden;
    std::vector<float> pointPos;
    std::vector<uint8_t> pointRgba, pointRgbaHidden;
    int revision = 0;
};

struct MarkingFlow {
    knit::MarkingSession session;
    int generation = -1;              // the session's last integrated generation
    std::string openedPath;           // the file the session was opened for
    knit::MarkingSet set;
    bool stale = false;               // the set was made on another mesh: shown in the list, never applied
    std::string savedTo;              // where every edit is saved (out/markings/<file>-<digest>.markings)
    std::string note;                 // the last save / load message
    int piece = -1;                   // the isolated piece (-1: none)
    // The isolated piece's knitting component (SurfaceMarking::knitComponent, 0: the first) that new markings and
    // suggestions go to and that reviews, builds and Knit this piece use. A piece knitted as several overlapping parts
    // (a collar over a body) has one component per part.
    int knitComponent = 0;
    int activeMarking = 0;            // the marking new points go to (0: none)
    int selectedPoint = -1;           // a point of the active marking the next click moves (-1: none)
    NeighbourLook around = NeighbourLook::Transparent;
    MarkView view = MarkView::Markings;
    bool relax = true;                // the build relaxes the stitch mesh (ReviewDepth::Relaxed)
    std::vector<knit::MarkingSet> undo, redo;
    std::vector<knit::MarkingTrace> traces;  // per marking of the set: its path on its piece
    std::vector<int> tracePiece;             // per marking of the set: its piece (-1: none)
    std::vector<std::string> pieceLabels;    // per piece: "white_bear_fabric, piece 4"
    std::vector<int> markingsOnPiece;        // per piece: how many markings it has
    std::vector<knit::TriangleNeighbours> neighbours;  // per piece, made when a marking first lies on it
    OverlayGeometry overlay;
    // suggestions
    std::shared_ptr<const knit::ProposalSet> proposals;  // the session's suggestions merged into the set (cards below)
    std::vector<int> cardSetup;      // per card (proposals->setups index): its setup number in the set
    int previewSetup = 0;            // a setup number previewed (0: none): reviewed and built on its ephemeral set
    // quick reviews: the newest request, and what the shown view must match to be drawn as current
    uint64_t quickSerial = 0, quickInputs = 0, quickFieldKey = 0;
    int quickFileGeneration = -1;
    struct EditStamp {
        uint64_t serial = 0;
        double at = 0;               // steady clock [s] when the request was made
        bool background = false;     // deep work (reading, a build, suggestions, a refinement) or a job ran at that moment
    };
    std::vector<EditStamp> pendingStamps;   // requests not answered yet
    std::vector<double> latencyMs;          // answered requests: request -> the frame that drew the answer
    std::vector<char> latencyBackground;    // per answered request: made while background work ran (EditStamp::background)
    int coalesced = 0;                      // requests superseded before an answer (never counted as answered)
    int coalescedBackground = 0;            // of those, made while background work ran
    uint64_t drawnQuick = 0;                // the serial of the quick view last drawn as current
    // the layout: automatic refresh, the builds as they arrived, and the baseline a build is compared with
    bool layoutRefresh = true;
    double layoutDueAt = 0;          // steady clock [s] when an automatic layout build is due (0: none)
    uint64_t autoLayoutInputs = 0;   // the last automatic build's inputs (never twice for the same, even without rows)
    struct ArrivedBuild {
        std::shared_ptr<const knit::MarkedPieceView> view;
        double at = 0;               // steady clock [s] when it arrived
    };
    // the session's last build of every piece and knitting component, taken over once per frame (integrateMarking), so
    // what a frame draws stays alive for the whole frame
    std::vector<ArrivedBuild> builds;
    // the explicit "Build stitch mesh" waiting for its result: when it arrives for the piece and component still shown,
    // with a mesh to draw, the Stitch Mesh view is selected (automatic layouts and other pieces' results never select it)
    uint64_t explicitBuild = 0;      // its request number (0: none)
    // the Stitch Mesh view's switches
    bool stitchSource = false;       // the piece's own surface under the stitch faces (off: the curved piece would poke
                                     // through the flat faces; on, it shows where the mesh lies inside or outside it)
    bool stitchHighlight = false;    // the faces coloured by what the checks found on them, instead of by kind
    std::shared_ptr<const knit::MarkedPieceView> baseline;   // the accepted (or first saved) setup's layout
    uint64_t baselineKey = 0;        // file digest, piece, unit, gauge, preparation of the baseline
    double baselineAt = 0;
    bool baselineFromNextBuild = false;     // a setup was accepted: its next layout becomes the baseline
    int highlightFinding = -1;              // a clicked finding (index into highlightReview's findings)
    bool highlightQuick = false;            // the finding is the quick view's (else the shown build's)
    // a scripted burst of point moves (mark-burst), one per frame
    struct BurstMove {
        int id = 0, point = 0;
        knit::Vec3 position;
    };
    std::vector<BurstMove> burst;
};

constexpr size_t kUndoLevels = 100;

// --- what changes it (marking_flow.cpp)

// Once per frame: opens the imported file when the Mark step is shown, integrates what the session published.
void integrateMarking(SceneUi& ui);
void isolatePiece(SceneUi& ui, int piece);              // -1: none
// The isolated piece's knitting component that new markings, suggestions, reviews, builds and Knit this piece use (0: the
// first); the number of components the piece has so far (its markings' highest component + 1, at least the active one).
void setKnitComponent(SceneUi& ui, int knitComponent);
int knitComponentsOnPiece(const SceneUi& ui);
void addMarking(SceneUi& ui, knit::MarkingMeaning meaning);
void selectMarking(SceneUi& ui, int id);
void setMarkingMeaning(SceneUi& ui, int id, knit::MarkingMeaning meaning);
void setMarkingShape(SceneUi& ui, int id, knit::MarkingShape shape);  // only shapes the meaning allows
void deleteMarking(SceneUi& ui, int id);
void deleteSelectedPoint(SceneUi& ui);
void acceptProposal(SceneUi& ui, int id);
// Moves control point `point` of marking `id` to the surface point of its piece nearest `filePosition` (file coordinates):
// the same edit as a click that moves a selected point (undo, retrace, autosave; a proposal becomes the user's marking,
// keeping its proposal id). False when there is no such marking, point or piece.
bool moveMarkingPoint(SceneUi& ui, int id, int point, const knit::Vec3& filePosition);
// The shown build's review as text (measures by group, findings with their causes): the script's mark-report.
std::string markReviewText(const SceneUi& ui);
// The file as the Mark step read it: its pieces (label, triangles, bounds in file units and mm) and the whole file's size
// in mm at the import's unit: the script's mark-inventory (a prepared file checked after its reimport).
std::string markInventoryText(const SceneUi& ui);
void rejectProposal(SceneUi& ui, int id);
// The reference the mesher prepares for the set ("facets" or "loop_limit", saved with the markings, v3): an edit like any other.
void setMarkingPreparation(SceneUi& ui, const std::string& preparation);
void undoMarking(SceneUi& ui);
void redoMarking(SceneUi& ui);
void saveMarkingsAs(SceneUi& ui, const std::string& path);
void loadMarkingsFrom(SceneUi& ui, const std::string& path);
void setMarkView(SceneUi& ui, MarkView view);
void setNeighbourLook(SceneUi& ui, NeighbourLook look);
void buildMarkedPiece(SceneUi& ui);                      // the footer: Method A on the isolated piece
// Places a point of the active marking on the isolated piece at the surface point nearest `filePosition` (file
// coordinates): what a click there does (undo, retrace, autosave). False without an active marking or piece.
bool placeMarkingPoint(SceneUi& ui, const knit::Vec3& filePosition);
// The footer's "Knit this piece (guided)": knit.route = guided on the isolated piece with the current markings file and
// their digest (Knit is turned on), then the Knit step. False when no piece is isolated or the markings are stale.
bool knitThisPiece(SceneUi& ui);
// A click in the viewport: isolates the piece under it, places or moves a point of the active marking, or (Select)
// selects the active marking's point nearest to it.
void markClick(SceneUi& ui, const Camera& camera, const ViewPixel& pixel, ClickIntent intent);

// --- what reads it

// The last build of the isolated piece's knitting component (null: none since the file was opened). A failed build is a
// build too: it has no mesh and says why (MarkedPieceView::use).
const knit::MarkedPieceView* shownBuild(const SceneUi& ui);
// When the shown build arrived (steady clock [s]; 0 without one).
double shownBuildArrivedAt(const SceneUi& ui);
// Writes the shown build's stitch polygons to an OBJ file at `path` (file units and axes, polygons kept); the result goes
// to the step's note. False without a drawable mesh or when the file cannot be written.
bool exportStitchMesh(SceneUi& ui, const std::string& path);
// The file name Export Stitch Mesh offers: "<file stem>-<piece label>-stitch-mesh.obj".
std::string stitchMeshFileName(const SceneUi& ui);
// The shown build was made from other markings (or another preview), gauge or unit than today's (it stays shown, labelled
// outdated). Its depth is not part of it: an unrelaxed layout of today's markings is current, and says it is unrelaxed.
bool shownBuildOutdated(const SceneUi& ui);
// The file a set of this file is saved to automatically ("" before the file is read).
std::string autosavePath(const SceneUi& ui);

// --- assisted stockinette (marking_flow.cpp)

// "Suggest markings for this piece": proposals for the isolated piece on the session's deep lane.
void suggestMarkings(SceneUi& ui);
// Rule 2 on a card (1-based): a trial relaxation of that setup's layout, in the background.
void refineSuggestion(SceneUi& ui, int card);
void cancelSuggestionWork(SceneUi& ui);   // stops a running suggestion, refinement or build
// Preview a card's setup (0: stop): its markings shadow the in-use ones of their roles in the quick review and the
// preview build; nothing is saved and the Knit input does not change. Any edit ends a preview first.
void previewSuggestion(SceneUi& ui, int card);
// Accept a card's setup as one undo step (replaceUser: the user's own markings of its roles are replaced). False when the
// setup is outdated or conflicts with user markings and replaceUser is false.
bool acceptSuggestion(SceneUi& ui, int card, bool replaceUser);
void rejectSuggestion(SceneUi& ui, int card);
// What a card's setup is now: "offered", "accepted", "rejected", "edited" (a marking edited or deleted), "outdated" (made
// for another gauge, unit or reference).
std::string suggestionState(const SceneUi& ui, int card);
// The user's markings a card's acceptance would replace (by role), as "Start, Boundary" ("" when none).
std::string suggestionReplaces(const SceneUi& ui, int card);
// The set the quick review and a build see now: the preview's ephemeral set, or the saved set; previewIds: the offered
// markings named for a research preview.
knit::MarkingSet effectiveMarkingSet(const SceneUi& ui, std::vector<int>& previewIds);
// The quick view when it answers the newest request (serial, file generation and inputs): null otherwise ("Checking ...").
std::shared_ptr<const knit::MarkingSession::QuickView> currentQuick(const SceneUi& ui);
void setLayoutRefresh(SceneUi& ui, bool on);
// The comparison of the shown build with the baseline, as lines for people ("" without a matching baseline).
std::string baselineComparison(const SceneUi& ui);
// The quick view's review as text (the script's mark-quick-report), and the latency record (mark-latency).
std::string quickReviewText(const SceneUi& ui);
std::string latencyText(const SceneUi& ui);
// Queue point moves (one per frame, like a drag): marking id, control point, positions in file coordinates.
void queueBurst(SceneUi& ui, int id, int point, const std::vector<knit::Vec3>& positions);
// Called after the frame is shown: the quick view drawn in it answers its request (latency).
void markFrameShown(SceneUi& ui);
// Highlight a finding (index into the shown review's findings; quick: the quick view's): its range, row and places.
void highlightFinding(SceneUi& ui, int index, bool quick);

} // namespace viewer
