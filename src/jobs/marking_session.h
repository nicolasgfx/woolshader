#pragma once
// The Mark step's background work (GL-free, tested without a window), on two lanes.
//
//   deep lane   read the imported file as the Import step reads it, find its pieces (model/surface_pieces.h) and make their
//               display; on request build one marked piece's stitch mesh with Method A (stitchmesh/marked_build.h,
//               private), propose setups for a piece (stitchmesh/mark_proposals.h) or refine one setup's start (rule 2: one
//               trial build). One worker; a newer deep request cancels the running one (CancelToken), whose result is
//               never published (ImportSession's rule: every request carries a number, only the current one publishes).
//   quick lane  the quick review while drawing (stitchmesh/quick_review.h) on its own worker and pending slot, so it never
//               waits behind a long build and never cancels one. The piece is prepared once (PreparedPiece, kept and shared
//               with builds); the rows are checked against the field of the last build when that field matches. A quick
//               result is published only when it answers the newest request of the open file.
//
// Everything the viewer reads is model data: proposals are ProposalSet (model/marking_review.h); the row field and the
// prepared piece stay inside (stitchmesh is a private dependency of jobs).
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "io/mesh_read.h"
#include "jobs/built_scene.h"
#include "model/cancel.h"
#include "model/marking_review.h"
#include "model/stitch_mesh_data.h"
#include "model/surface_marking.h"
#include "model/surface_path.h"
#include "model/surface_pieces.h"

namespace knit {

// A built piece as the viewer shows it, in file units and axes (the import preview's coordinates). Every finished build
// publishes one, also a build that made no stitch faces or failed: then `use` says why there is no mesh, and the piece's
// earlier mesh is gone (never shown in its place).
struct MarkedPieceView {
    int piece = -1;
    int knitComponent = 0;                     // the piece's knitting component built
    uint64_t request = 0;                      // the session request that made it (MarkingSession::build's return value)
    uint64_t inputs = 0;                       // markBuildInputs of the request: the viewer compares it with today's to say "outdated"
    std::vector<int> previewIds;               // offered proposals the build previewed (a setup's preview), empty for saved intent
    MarkingReview review;
    std::vector<int> markingIds;               // the markings used, as the review's findings name them
    std::vector<MarkingTrace> traces;          // per used marking, on the piece's triangles (the findings' trace ranges index these)
    std::vector<std::vector<Vec3>> rings;      // row contours
    std::vector<std::vector<Vec3>> stitchEdges;  // the stitch mesh's faces as closed outlines (lifted like the display faces)
    StitchMeshData stitch;                     // the polygons as built, not lifted (what Export Stitch Mesh writes)
    StitchMeshQuality quality;                 // what the build's checks measured
    StitchMeshUse use;                         // whether the mesh can be drawn and saved, and why not
    // meshes[0]: the piece, per triangle coloured by region (kept / omitted visible / omitted hidden / outside). When
    // use.drawable: meshes[1] the stitch faces coloured by kind (regular, increase, decrease, special, unsupported, pole
    // start), meshes[2] the same faces coloured by what the checks found on them (faces without a finding paled).
    std::shared_ptr<const BuiltScene> scene;
    std::vector<std::string> log;
};

class MarkingSession {
public:
    enum class Phase { Idle, Reading, Ready, Building, Proposing, Refining, Failed };

    // The last quick review, as the viewer draws it.
    struct QuickView {
        uint64_t serial = 0;           // the request it answers (review's return value)
        int fileGeneration = 0;        // the open() it belongs to
        uint64_t inputs = 0;           // quickReviewInputs of that request
        int piece = -1;
        MarkingReview review;
        std::vector<int> markingIds;
        std::vector<MarkingTrace> traces;   // the findings' trace ranges index these
        double ms = 0;                 // the review itself
        double preparedMs = 0;         // the piece's preparation when this request had to make it (0: warm)
    };

    struct Snapshot {
        Phase phase = Phase::Idle;
        std::string path, error;
        int generation = 0;  // increases with every published change (the viewer polls it)
        int fileGeneration = 0;  // increases with every open() and reset()
        // Ready and later: the file as read, its pieces and digest, and one display mesh per piece (named "piece <i>").
        std::shared_ptr<const ImportedMesh> mesh;
        std::shared_ptr<const SurfacePieces> pieces;
        std::shared_ptr<const BuiltScene> preview;
        uint64_t digest = 0;
        std::shared_ptr<const MarkedPieceView> built;  // the last build (null: none since the file was opened)
        // the last build of every piece and knitting component built since the file was opened (findBuild); a newer
        // build of the same piece and component replaces its entry, also when it failed
        std::vector<std::shared_ptr<const MarkedPieceView>> builds;
        std::shared_ptr<const ProposalSet> proposals;  // the last proposals (null: none since the file was opened)
        uint64_t fieldKey = 0;                         // the key of the row field kept from the last build (0: none)
        int fieldPiece = -1;
        std::shared_ptr<const QuickView> quick;        // the last quick review published
        bool quickBusy = false, quickPreparing = false;
        float progress = 0;                            // of a running proposal or refinement, 0 .. 1
        std::string progressText;
        int cancelledTasks = 0;   // deep tasks stopped before their end (superseded or cancelDeep); nothing of theirs is published
    };

    MarkingSession();
    ~MarkingSession();  // cancels the running work and joins both lanes
    MarkingSession(const MarkingSession&) = delete;
    MarkingSession& operator=(const MarkingSession&) = delete;

    void open(const std::string& path);
    // Ready: Method A on one piece with the set's markings on it; mmPerUnit is the import's unit, w and h the gauge.
    // previewOffered: offered proposals built as if accepted (a setup's preview; the set is then the preview's own).
    // knitComponent: the piece's knitting component built (SurfaceMarking::knitComponent; 0: the piece's first).
    // Returns the request number its MarkedPieceView carries (0: nothing read yet, no build started).
    uint64_t build(const MarkingSet& set, int piece, double mmPerUnit, double w, double h, ReviewDepth depth, const std::vector<int>& previewOffered = {},
                   int knitComponent = 0);
    // Ready: proposals for one piece on the reference preparation ("facets" or "loop_limit").
    void propose(int piece, double mmPerUnit, double w, double h, const std::string& preparation);
    // Ready, after proposals: rule 2 on setup number `setup` of the last proposals; the refined setup is added to them.
    void refine(int setup);
    void cancelDeep();   // stops a running build, proposal or refinement (nothing is published)
    // The quick review of `piece` (quick lane). Returns the request's serial; a result for an older serial is never
    // published.
    uint64_t review(const MarkingSet& set, int piece, double mmPerUnit, double w, double h, const std::vector<int>& previewOffered, int knitComponent = 0);
    void reset();  // back to Idle (no file)
    Snapshot snapshot() const;
    // Tests only: hold the next quick reviews this long before publishing, to show that an older one never overtakes.
    void delayQuickForTest(int ms);

private:
    struct Task;
    struct QuickTask;
    struct Kept;
    void loop();
    void quickLoop();
    mutable std::mutex m_;
    std::condition_variable cv_, quickCv_;
    std::unique_ptr<Task> pending_;
    std::unique_ptr<QuickTask> quickPending_;
    uint64_t request_ = 0, quickSerial_ = 0;
    Snapshot state_;
    std::unique_ptr<Kept> kept_;   // the prepared piece, the row field, the last proposal request (stitchmesh types, .cpp only)
    std::atomic<int> step_{0};
    int stepsTotal_ = 1;
    int quickDelayMs_ = 0;
    bool stop_ = false;
    std::shared_ptr<CancelToken> running_;   // the running deep task's token, cancelled by a newer deep request
    std::thread worker_, quickWorker_;
};

// The colours the Stitch Mesh view draws with (MarkedPieceView meshes 1 and 2), for the legend beside it.
struct StitchLegendEntry {
    std::string label;
    float rgb[3] = {0, 0, 0};
};
std::vector<StitchLegendEntry> stitchKindLegend();      // by kind: regular, increase, decrease, special, unsupported, pole start
std::vector<StitchLegendEntry> stitchFindingLegend();   // by finding: none, over the fit target, over twice it, structure, unsupported

// The last build of `piece`'s knitting component in a snapshot (null: none since the file was opened).
std::shared_ptr<const MarkedPieceView> findBuild(const MarkingSession::Snapshot& snapshot, int piece, int knitComponent);

// What a build of `piece`'s knitting component depends on: the markings it uses (effectiveMarkingsDigest: offered and rejected
// proposals, other components' and other pieces' markings change nothing), the unit, the gauge, the depth, the previewed
// proposals and the component. Equal values: the same build.
uint64_t markBuildInputs(const MarkingSet& set, const SurfacePieces& pieces, int piece, double mmPerUnit, double w, double h, ReviewDepth depth,
                         const std::vector<int>& previewOffered = {}, int knitComponent = 0);

// What a quick review of `piece`'s knitting component depends on: every marking of that component (in use or offered, with
// its state, so a preview and an edit both count), the set's digest and preparation, the unit, the gauge, the previewed
// proposals and the component.
uint64_t quickReviewInputs(const MarkingSet& set, const SurfacePieces& pieces, int piece, double mmPerUnit, double w, double h, const std::vector<int>& previewOffered,
                           int knitComponent = 0);

// The file read as the Mark step and stitchmesh_lab read it: its pieces, digest and per-piece display. False with a
// message when the file cannot be read.
struct MarkedFile {
    std::shared_ptr<const ImportedMesh> mesh;
    std::shared_ptr<const SurfacePieces> pieces;
    std::shared_ptr<const BuiltScene> preview;
    uint64_t digest = 0;
};
[[nodiscard]] bool readMarkedFile(const std::string& path, MarkedFile& out, std::string* err);

} // namespace knit
