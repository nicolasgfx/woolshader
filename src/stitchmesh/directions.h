#pragma once
// Stage 3 of the stitch-mesh layout: the knitting (wale) direction on every face (Wu et al. 2018, section 6), from the
// labelled coarse mesh of stage 2.
//
// Faces joined across their wale edges form ROWS: a path (it begins and ends with a triangle of one wale edge, or at
// the mesh border) or a cycle. Along a row the wale direction is one choice: which of a face's course edges are its
// bottom edges (the previous row's loops come in there) and which its top edges. A face's course edges fall on two
// SIDES of the row, A and B (travelling the row from its entry wale edge to its exit, side B's edges come after the
// entry counter-clockwise). One binary per row says whether side A is the top. Two rows that share course edges
// MATCH when one calls the edge top and the other bottom; the program maximizes the number of matched course edges
// (weighted by count per row pair, a meta-graph with one node per row), solved with HiGHS. Mismatches that remain
// are kept and listed: they are connectivity the later stages must support or name unsupported, never a seam to hide.
//
// The objective is the same when every row of a connected meta-graph component flips, so each component is oriented
// by a recorded anchor: its area-weighted mean wale direction points along `up` (a repository tie-break), or the
// caller's explicit sign per component (the user's growth control).
#include <string>
#include <vector>

#include "model/status.h"
#include "model/vec.h"
#include "stitchmesh/labels.h"
#include "stitchmesh/stage_report.h"

namespace knit {

struct DirectionOptions {
    Vec3 up{0, 0, 1};                 // the components' anchor: rows climb this way unless reversed
    std::vector<char> reverseComponent;  // per meta-graph component (index), 1: flip the anchor's choice; empty: none
    double timeLimitS = 60;
    int maxNodes = 200000;            // branch-and-bound node limit: a deterministic stop
};

struct DirectedMesh {
    const PolygonMesh* mesh = nullptr;  // the labelled mesh (LabelProduct::mesh); must outlive the product
    std::vector<int> rowOfFace;         // per face
    std::vector<std::vector<int>> rows; // faces in row order
    std::vector<char> rowClosed;        // per row: a cycle
    std::vector<int> componentOfRow;    // meta-graph component
    // per face, per corner edge k (face-local): the role of a course edge, 1 top, 0 bottom; -1 for a wale edge
    std::vector<std::vector<char>> top;
    // per face: the wale direction, from the bottom course edges' midpoints to the top ones' (zero where one side is empty)
    std::vector<Vec3> waleDirection;
    std::vector<Vec3> faceCentre;
    // interior course edges whose two faces disagree: both bottom or both top. Two kinds, kept apart: between two rows
    // (the row program's residue) and within one row (a row that meets itself, which no orientation can match)
    struct Mismatch {
        int a = -1, b = -1, faceA = -1, faceB = -1;
        bool bothTop = false;
        bool sameRow = false;
    };
    std::vector<Mismatch> mismatches;
};

struct DirectionProduct {
    DirectedMesh directed;
    int rows = 0, closedRows = 0, openRows = 0, metaEdges = 0, components = 0;
    double matchedWeight = 0, totalWeight = 0;   // course edges between different rows: matched, all
    int sameRowMatched = 0, sameRowMismatched = 0;  // course edges between the two sides of one row / one side
    int interRowMismatched = 0;                     // mismatches between different rows (the list holds both kinds)
    int boundaryCourseEdges = 0;
    std::string solverStatus;
    bool optimal = false;
    double msSolve = 0;
    StageReport report;  // PASS iff every row is coherent and the program optimal; UNRESOLVED when feasible but not optimal,
                         // or when mismatches remain (they are listed; a later stage decides); FAIL when a face has no valid row
};

// False only when the labelled mesh is unusable (a face with no or more than two wale edges, an edge without a role).
bool assignDirections(const PolygonMesh& m, const EdgeLabels& labels, const DirectionOptions& o, DirectionProduct& out, std::string* err);

} // namespace knit
