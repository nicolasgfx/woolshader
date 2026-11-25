#pragma once
// The quick review while drawing (ReviewDepth::Quick):
// what can be said about a piece's markings without building a layout, in milliseconds on a prepared piece.
//
//   always      the markings resolved and traced on the prepared reference (marked_build.h's resolvePieceMarkings: the
//               piece, the shape, an unfinished marking as Incomplete, a trace that cannot be drawn), the one start,
//               termination and path, and a path or curve that crosses itself
//   with a current field
//               the row-start path's crossings with every row, counted by the layout's own code (rowCrossings: the same
//               count the build makes on the same field and trace), so rows crossed twice, missed rows and turn-backs
//               are named with their places
//   otherwise   "row checks pending": the rows are not known for this start, finish, gauge and reference, and no row
//               finding is given (never a pass)
//
// Not quick: the trim, the typed mesh, coverage and relaxation (a layout build). The review marks which measures it knows.
#include <memory>
#include <vector>

#include "model/marking_review.h"
#include "model/surface_marking.h"
#include "model/surface_path.h"
#include "model/surface_pieces.h"
#include "stitchmesh/marked_build.h"
#include "stitchmesh/prepared_reference.h"

namespace knit {

struct QuickRequest {
    std::shared_ptr<const PreparedPiece> reference;   // the piece prepared once (required)
    std::shared_ptr<const MarkedField> field;         // the last layout's field (null, or of other inputs: checked by its key)
    const MarkingSet* markings = nullptr;
    const SurfacePieces* pieces = nullptr;
    int piece = -1;
    int knitComponent = 0;                            // the piece's knitting component reviewed (SurfaceMarking::knitComponent)
    double mmPerUnit = 1, w = 5.5, h = 3.4;
    ReferencePreparation preparation = ReferencePreparation::Facets;
    std::vector<int> previewOffered;                  // offered proposals reviewed as if accepted (a preview's setup)
};

// markingIds, traces (may be null): the markings used and their paths on the piece's triangles in file units, as a build
// gives them, so a finding's trace range lands on the drawn path.
MarkingReview quickReview(const QuickRequest& q, std::vector<int>* markingIds, std::vector<MarkingTrace>* traces);

} // namespace knit
