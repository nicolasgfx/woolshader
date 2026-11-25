#pragma once
// Method A on a marked piece of an imported file: the markings a
// user drew (or accepted) on one piece become the mesher's input, and the product says, per marking, whether it served
// its purpose. One backend for the app (jobs/marking_session.h) and the lab (stitchmesh_lab marked).
//
//   1. the piece's markings are resolved to the piece's triangles, then to the prepared reference: the piece in mm,
//      midpoint-refined (prepared_reference.h; the surface does not move, every point maps back exactly); paths traced
//      (model/surface_path.h); shapes checked per meaning. Input: user markings and accepted proposals; offered proposals
//      only when named in previewOffered (a labelled research preview, never saved as acceptance); rejected ones never.
//   2. the route follows the start's shape:
//      band  a start curve: the start and finish (or boundary) curves trim the piece (surface_trim.h); the band between
//            them is knitted, the caps beyond are omitted on purpose; the field runs from the start curve.
//      pole  a start point (experimental): the field runs from the point; ring 0 is the
//            pole, strip 0 the gathered pole start; the finish is a curve (the region it bounds around the start is
//            knitted) or a point (the bind-off is the field's level through it; beyond it is omitted).
//   3. the row-start path places every row's vertex 0 (guided_rows.h) — it never cuts;
//   4. the typed stitch mesh, its virtual shift and (ReviewDepth::Relaxed) the CP3 relaxation on the virtual faces, with
//      the markings' constraints by RelaxAnchors (pole pinned; row origins on the path or the yarn's ends pinned; bind-off on its curve/level);
//   5. coverage of the knit region by the stored tiling and by the solver's (virtual) faces; omitted area split into
//      what another piece of the same component hides and what stays visible.
//
// A construction the mesher does not have (cuts, attachments, a finish point after a start curve) is reported as
// FindingCause::AlgorithmLimitation, never as a marking to move. ReviewDepth::Quick has no cheaper path yet: it runs
// as Built.
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include "model/marking_review.h"
#include "model/mesh_query.h"
#include "model/surface_marking.h"
#include "model/surface_path.h"
#include "model/surface_pieces.h"
#include "stitchmesh/coverage.h"
#include "stitchmesh/guided_rows.h"
#include "stitchmesh/prepared_reference.h"
#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/row_rings.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"
namespace knit {
enum class MarkedRegion : char { Kept, OmittedStart, OmittedFinish, OutsideBoundary };
const char* toString(MarkedRegion r);

// One piece's prepared reference (prepared_reference.h), made once and shared by quick reviews and builds of the piece: the
// piece in mm, refined (levels: the automatic rule, mean edge <= h / 2), and its triangle neighbours.
struct PreparedPiece {
    uint64_t key = 0;                  // preparedPieceKey of what it was made from
    int piece = -1;
    RefinedMesh refined;
    TriangleNeighbours nb;
    double seconds = 0;
};
uint64_t preparedPieceKey(uint64_t fileDigest, int piece, double mmPerUnit, double h, ReferencePreparation p);
void preparePiece(const SurfacePieces& pieces, int piece, uint64_t fileDigest, double mmPerUnit, double h, ReferencePreparation p, PreparedPiece& out);

// The row field a build laid out, kept so that a quick review can check a row-start path against it without building
// (quick_review.h). It depends on the start, the finish, the unit, the row height and the reference only; the path, ring1
// and the stitch width do not change the field, its levels or its contours. Made in place and shared read-only (the
// closest-point query refers to `mesh`): never copied.
struct MarkedField {
    uint64_t key = 0;                          // markedFieldKey of the build's start and finish
    std::shared_ptr<const PreparedPiece> reference;
    std::string route;                         // "band" | "pole"
    TriMesh mesh;                              // the mesh the rows builder ran on (the band, the disc, or the reference)
    std::vector<int> parentTri;                // per mesh triangle: its reference triangle (empty: the reference itself)
    std::vector<double> t;                     // per mesh vertex
    std::unique_ptr<MeshQuery> query;          // pole route: pathOnField's lookup on `mesh`
    int rows = 0;                              // R
    double rowHeight = 0;                      // h' = T / R
    bool loopFinish = false;                   // pole route: the bind-off is a curve (ring R is not crossed)
    std::vector<Contour> contours;             // band route: the contour per ring the builder tested (empty: a contour event)
    std::vector<int> contoursPerRing;          // per ring k = 0 .. R: contours at the level (1: fine), path-independent
};
// file digest, preparation, piece, unit, h, and the resolved start and finish (id, meaning, shape, points bit for bit)
struct ResolvedMarkings;
uint64_t markedFieldKey(const ResolvedMarkings& m, uint64_t fileDigest, const std::string& preparation, int piece, double mmPerUnit, double h);
// What the mesh relaxation holds of the layout the markings gave (relax_mesh.h: pinned vertices and guide lines).
enum class RelaxAnchors {
    // The pole pinned and each row's first vertex on the row-start path (the first integration's rule). The path is not a
    // wale of the relaxed mesh, so the rows slide against their held first vertex: on the teddy body whole rings moved
    // 12-17 mm, with inverted cells.
    RowOriginsOnPath,
    // The yarn's ends on an open ring pinned: the first vertex of the bind-off ring, and of the cast-on ring on a band. The
    // solver's faces are the virtual (helical) ones, which leave the closing edge of those rings without a face, so
    // nothing else keeps these vertices from sliding into their neighbours. The rows between are free, and so is a
    // pole: pinned, it turned a cell over beside it on the teddy arms (MarkedEvaluation::poleMoved says how far it went).
    YarnEnds,
    // Nothing held by the layout (the default since the teddy mesh repair):
    // the border rings are held by BorderRings::PlacedFirst instead, and the border edge the virtual faces leave out
    // has its spring (BorderSprings::StoredBorderEdges), so no single vertex needs holding. Also the proposals' probe of
    // where the pole would go.
    None,
};
// Where the vertices of a border ring (the bind-off ring; on a band also the cast-on ring) end up on their boundary loop.
enum class BorderRings {
    // They slide along the loop under the relaxation's spring forces (CP3). The forces know nothing of the loop's bends:
    // sides grew from 5.6 to 13.7 mm on the teddy body and the border polygon left the loop by up to 3 mm.
    SlideInRelaxation,
    // After the relaxation, placed on the loop for the smallest distance between the border polygon and the loop
    // (border_placement.h), then held while the other vertices are relaxed once more around them. That second relaxation
    // is kept only when it ends with smaller forces than it started with: on the teddy arms it runs away along the
    // row-start column (a wale edge there shrank from 2.3 to 0.15 mm) and is dropped, leaving the placed rings on the
    // first relaxation's mesh. With RelaxAnchors::YarnEnds only.
    PlacedOnLoop,
    // Placed on the loop as laid out (border_placement.h, the first vertex within half a stitch of its place) and held
    // through the one relaxation: the rows below arrange themselves under a border that already fits. With
    // RowPhase::WalesUpright the rows hardly turn against the ring (under 1 mm on the teddy body), which is what makes
    // holding it from the start possible; with RowPhase::OnPath it left 15 inverted faces there.
    PlacedFirst,
};
// The rows of a piece that is knitted from a start point to its own opening.
enum class OpeningRows {
    // Levels of the distance from the start point, as far as the nearest point of the opening; the opening itself is the
    // last ring. Where the opening is not at one distance from the start, everything beyond its nearest point falls into
    // the last strip: on the inner half of a teddy ear cut along its rim the opening is 28.9 to 36.0 mm from the start, the
    // last strip is up to three rows tall and 22 % of the surface is left uncovered.
    DistanceFromStart,
    // Levels of t' = T t / (t + d), with d the distance to the opening and T the mean of t along the opening: the start is
    // still the lowest point and the opening is a level, so the rows lie between the two and each strip is as many
    // stitch heights tall as the others along any line from the start to the opening.
    OpeningIsALevel,
};
// The springs along the stitch mesh's border in the relaxation.
enum class BorderSprings {
    // Those of the solver's (virtual, helical) faces only. The edge between the first and second vertex of the bind-off
    // ring has none: on the teddy body the second vertex slid through the held first one in four steps.
    VirtualFacesOnly,
    // Every border edge of the stored mesh has its course spring (RelaxGuides::courseEdges): the border ring is a closed
    // chain, as it is in the stored mesh and in Method B.
    StoredBorderEdges,
};
// The experimental constructions and solver options; the defaults are this integration's (marked_build.cpp states each).
struct MarkedBuildOptions {
    bool adjacentShaping = true;       // shaping next to shaping where the frozen spacing cannot be met
    bool boundaryShaping = true;       // kp / d12k in a cast-on or bind-off strip
    RelaxAnchors anchors = RelaxAnchors::None;       // what the relaxation holds of the layout
    BorderRings borderRings = BorderRings::PlacedFirst;
    RowPhase rowPhase = RowPhase::WalesUpright;      // guided_rows.h: rings turned so that their wale edges stand upright
    OpeningRows openingRows = OpeningRows::OpeningIsALevel;
    BindOffCount bindOffCount = BindOffCount::OfRoundBelow;   // guided_rows.h: no shaping in the bind-off strip
    BorderOrder borderOrder = BorderOrder::Kept;     // relax_mesh.h: border vertices keep their order around the loop
    BorderSprings borderSprings = BorderSprings::StoredBorderEdges;
    FaceCorners faceCorners = FaceCorners::Kept;     // relax_mesh.h: no step turns a face corner over
    // relax_mesh.h. The fixed time step: with the rings turned for upright wales the forces no longer grow, and taking
    // steps back ended the two teddy arms in different states (the right one with two faces at 0.56 and 0.59 mm, the left
    // at 0.50; with the fixed step both at 0.48 to 0.52).
    StepControl stepControl = StepControl::FixedTimeStep;
    bool smoothFilterNormals = false;  // solver note R1 (relax_mesh.h); off: the facet normals of CP3
    bool kinkConstraints = false;      // solver note R2 (relax_mesh.h); off: CP3's one-sided projection at kinks
    // The prepared reference (prepared_reference.h): the layout is built and
    // relaxed on it. Facets keeps the imported surface; LoopLimit is the smooth-reference variant.
    ReferencePreparation preparation = ReferencePreparation::Facets;
    // Controlled comparison: the layout as built on `preparation`, carried to this preparation of the same piece at the same
    // levels (every stitch vertex, guide and border keeps its refined triangle and weights) and relaxed there.
    bool relaxOnOther = false;
    ReferencePreparation relaxPreparation = ReferencePreparation::LoopLimit;
    // With relaxOnOther, the markings' constraint lines on the other preparation (only the pole route with a finish point walks
    // on the whole reference; the other routes always carry):
    //   Carried      the lines as laid out, carried point by point (piecewise linear per triangle: a straight contour gains
    //                kinks on the other surface)
    //   MatchedLevel the row-start path traced again between the same control points; the finish contour the level of the
    //                other surface's field that the carried bind-off ring sits on on average, the knit region clipped there
    //                (the controlled comparison's default: the constraints by the same rules, fitting the carried layout)
    //   FinishPoint  as MatchedLevel, but the level through the finish point (the rule of a full build; geodesic lengths
    //                change with the surface, so the carried bind-off ring is off it)
    enum class OtherGuides : char { Carried, MatchedLevel, FinishPoint };
    OtherGuides otherGuides = OtherGuides::MatchedLevel;
};
struct MarkedBuildRequest {
    const SurfacePieces* pieces = nullptr;
    int piece = -1;
    int knitComponent = 0;                 // the piece's knitting component (SurfaceMarking::knitComponent) built
    double layerOffsetMm = 0;              // the reference moved this far out along its normals (a component over another)
    const MarkingSet* markings = nullptr;  // the file's markings: those on this piece and component are used
    double mmPerUnit = 1;                  // the import's unit
    double w = 5.5, h = 5.5;               // gauge [mm]
    ReviewDepth depth = ReviewDepth::Relaxed;
    std::vector<int> previewOffered;       // offered proposals evaluated as if accepted (research preview)
    int refine = -1;                       // midpoint refinement levels of the piece; -1: automatic (mean edge <= h / 2)
    MarkedBuildOptions options;
    // An already prepared reference of this piece (the quick review's); used when its key matches (automatic refinement
    // only), else the build prepares its own. It never changes the result.
    std::shared_ptr<const PreparedPiece> prepared;
};
// What a lab report needs beyond the app's review, in mm on the prepared reference.
struct MarkedEvaluation {
    std::string route;                         // "band" or "pole"
    int refineLevels = 0;
    TriMesh reference;                         // the prepared reference (refined piece, mm)
    TriMesh knitRegion;                        // the part of it the layout must cover
    std::vector<std::vector<int>> regionLoops; // its border loops
    TypedStitchMesh before, after;             // stored tiling before and after relaxation (mm); after == before without relaxation
    VirtualShift shift;
    VirtualShiftCheck shiftBefore, shiftAfter;
    TypedMeshCheck storedBefore, storedAfter, virtualBefore, virtualAfter;
    CoverageReport coverStoredBefore, coverStoredAfter, coverVirtualBefore, coverVirtualAfter;
    RelaxProduct relax;
    bool relaxed = false;
    std::vector<RelaxFace> solverFaces;
    RelaxGuides guides;                        // the constraints the relaxation used
    std::vector<std::vector<int>> relaxLoops;  // reference loops the border slid on
    std::vector<char> border;
    double finishLevel = 0;                    // pole route with a finish point: the bind-off level [mm]
    double poleMoved = 0;                      // pole route: how far the relaxation moved the pole [mm]
    double rowOriginMean = 0, rowOriginMax = 0;  // each ring's vertex 0 against the row-start path after relaxation [mm]
    // the preparation (smooth-reference report): which surface, and how far it is from the imported one
    std::string preparation = "facets";        // of the layout's reference
    std::string relaxPreparation = "facets";   // of the reference the relaxation ran on (and the checks measure against)
    uint64_t referenceDigest = 0;              // meshDigest of `reference`
    TriMesh original;                          // the imported piece at the same refinement (facets, mm)
    double changeMax = 0, changeP95 = 0;       // reference vertices -> the original surface [mm] (0 for facets)
    double changeBackMax = 0, changeBackP95 = 0;  // original vertices -> the reference [mm]
    double areaOriginal = 0, areaReference = 0;   // [mm2]
    TypedMeshCheck storedAfterVsOriginal;      // the relaxed stored faces measured against the original surface
    TypedCheckOptions checkOptions;            // the limits the mesh checks used (gauge-scaled)
    CoverageOptions coverOptions;              // the limits the coverage checks used (gauge-scaled)
    struct MarkingShift { int id = 0; double controlMax = 0; double traceHausdorff = 0; };
    std::vector<MarkingShift> markingShift;    // per used marking: control points moved by the preparation, trace difference [mm]
};
struct MarkedBuildProduct {
    MarkingReview review;
    std::vector<int> markingIds;               // the markings on this piece used as input, in set order
    std::vector<int> previewIds;               // of those, the offered proposals used by request (research preview)
    std::vector<MarkingTrace> traces;          // per used marking: its path on the piece, file units (piece triangles)
    std::vector<MarkedRegion> region;          // per piece triangle (empty: no trim)
    std::vector<char> regionHidden;            // per piece triangle: an omitted triangle inside another piece of the component
    std::vector<std::vector<Vec3>> rings;      // row contours for the preview, file units
    TypedStitchMesh stitchMesh;                // file units: the stored tiling, relaxed at ReviewDepth::Relaxed
    MarkedEvaluation evaluation;
    std::shared_ptr<const MarkedField> field;  // the row field, once the rows were laid out (null before that)
    std::vector<std::string> log;
};
// False only for an unusable request (no pieces, no such piece, no markings set, a gauge or unit <= 0); every problem
// of the markings is a finding in out.review.
[[nodiscard]] bool buildFromMarkings(const MarkedBuildRequest& r, MarkedBuildProduct& out, std::string* err);

// Step 1 of the build on its own, shared with the quick review (quick_review.h): the piece's markings in use resolved to
// the prepared reference, traced, and checked per marking (piece, shape, completeness, the trace), and the one start,
// termination and row-start path chosen. Offered proposals only when named in previewOffered, rejected ones never.
struct MeaningChoice {
    int start = -1, finish = -1, rowStart = -1;   // index into ResolvedMarkings::used (-1: none)
    bool finishIsBoundary = false;
    bool incompleteStart = false, incompleteFinish = false, incompleteRowStart = false;   // that role's marking is still being drawn
};
struct ResolvedMarking {
    const SurfaceMarking* marking = nullptr;
    std::vector<SurfacePoint> local;   // its points on the prepared reference
    MarkingTrace trace;                // on the prepared reference (mm); a border loop: the loop itself
    std::vector<int> borderLoop;       // a closed Start or Finish lying on one border loop of the piece: that loop on the
                                       // prepared reference (vertex ids, walked as the control points go), else empty
    bool resolved = false;             // every point lies on this piece
    bool preview = false;              // an offered proposal named for a research preview
};
struct ResolvedMarkings {
    std::vector<ResolvedMarking> used;
    MeaningChoice choice;
    int offeredSkipped = 0;            // offered proposals on this piece that are not input
    std::vector<int> markingIds, previewIds;
    std::vector<MarkingTrace> traces;  // per used marking: its path on the piece, file units (piece triangles)
};
// The per-marking findings, in set order.
void resolvePieceMarkings(const MarkedBuildRequest& r, const RefinedMesh& refined, const TriangleNeighbours& nb, MarkingReview& review, ResolvedMarkings& out);
// The findings about the choice as a whole: a research preview, offered proposals not used, a missing start, finish or
// path, a finish point after a start curve.
void reviewMeaningChoice(const ResolvedMarkings& m, MarkingReview& review);

// The row-start path's crossings with each row of a field, counted by the layout's own rules: band route pathCrossing on
// the contours the builder tested; pole route pathOnField and levelCrossings. Rows 1 .. R - 1 (and R for a finish level).
struct RowCrossing {
    int row = 0;
    int count = 0;                     // crossings
    bool split = false;                // the row is not one closed contour (path-independent)
    std::vector<double> pathParameters;// along the trace: point index + fraction of each crossing
    std::vector<Vec3> at;              // their positions (mm)
};
std::vector<RowCrossing> rowCrossings(const MarkedField& f, const MarkingTrace& path);
// The findings those crossings give, the same in a build and a quick review: split rows (on the start), rows the path
// crosses twice or more (per row, with positions and the trace range), rows it misses, where it turns back, or "crosses
// every row once". toFile: mm to file units for the positions.
void crossingFindings(const std::vector<RowCrossing>& rows, const std::string& route, int startId, int pathId, double toFile, MarkingReview& review);
} // namespace knit
