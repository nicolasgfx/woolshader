#pragma once
// What the stitch mesher says about a piece's markings: findings tied to a marking (and a highlighted part of its
// path), and the layout's measures, kept in four separate groups (validity, layout quality, coverage, cost) instead
// of one score. A finding names its cause, so a problem of the algorithm (an unresolved pole construction) or of the
// solver (a relaxation that does not converge) is never presented as something moving a line would cure.
// Every measure says whether it was computed: a quick review knows row crossings but not the
// cells, and a measure that was not computed is never shown as 0 or PASS.
//
// The proposals the mesher makes for a piece are plain data here too (ProposalSetup, ProposalSet), so the app sees them
// without seeing the mesher: a setup is a compatible start, termination and row-start path made from one field.
#include <cstdint>
#include <string>
#include <vector>

#include "model/surface_marking.h"
#include "model/vec.h"

namespace knit {

// How deep a review goes: cheap checks while drawing, a rebuilt layout, or a rebuilt and relaxed one.
enum class ReviewDepth { Quick, Built, Relaxed };
// Incomplete (appended, so the older values keep their numbers): a marking still being drawn, not an error.
enum class FindingLevel { Ok, Warning, Problem, Incomplete };
enum class FindingCause {
    Marking,              // the marking itself: moving or changing it can help
    AlgorithmLimitation,  // a construction the mesher does not have yet (splits, joins, a pole end): no marking cures it
    Input,                // the surface or settings (no markings, a stale set, a degenerate piece)
    SolverLimitation      // a solver that did not converge or certify (the relaxation): never a reason to move a marking
};

struct MarkingFinding {
    int markingId = 0;  // 0: the piece as a whole
    FindingLevel level = FindingLevel::Ok;
    FindingCause cause = FindingCause::Marking;
    std::string message;
    int firstTracePoint = -1, lastTracePoint = -1;  // the highlighted part of the marking's path in the review's own trace (-1: none)
    int row = -1;                                   // a row to highlight (-1: none)
    std::vector<Vec3> at;                           // positions to highlight, file units (the crossings, a split contour)
};

// Whether a measure was computed, and for the current inputs.
enum class Known : char { No, Current, Outdated };

struct LayoutMeasures {
    // validity
    int contourEvents = 0;       // rows that are not one closed contour, or not crossed once by the row-start path
    int unsupportedCells = 0;    // typed-mesh faces of kind Unsupported
    bool checksPass = false;     // the built layout's own checks pass: connectivity and geometry, coverage, no contour event,
                                 // no unsupported cell, no Problem finding (a built mesh alone, MarkingReview::valid, is not that)
    // layout quality
    int rows = 0;
    int maxStitchCountJump = 0;  // largest |n(k+1) - n(k)|
    int exceptionalCells = 0;    // shaping, helix steps and other non-regular faces
    double distortionP95 = 0;    // relaxed cell distortion, 95th percentile (0: not relaxed)
    // coverage
    double coveredFraction = 0;  // of the kept surface
    double omittedAreaMm2 = 0;   // removed on purpose by start, finish and boundary curves (or beyond a finish level)
    double omittedHiddenMm2 = 0; // of that, inside another closed piece of the same component (hidden)
    double omittedVisibleMm2 = 0;// of that, not inside another piece: visible surface left unknitted
    bool coveragePass = false;
    // cost
    double seconds = 0;
    // per measure: computed for these inputs, for older ones, or not at all
    Known contourEventsKnown = Known::No, unsupportedKnown = Known::No, checksKnown = Known::No, rowsKnown = Known::No, jumpKnown = Known::No,
          exceptionalKnown = Known::No, distortionKnown = Known::No, coverageKnown = Known::No, omittedKnown = Known::No, secondsKnown = Known::No;
};

struct MarkingReview {
    ReviewDepth depth = ReviewDepth::Built;
    bool valid = false;                    // a stitch mesh was built (not: the layout passes its checks; see checksPass)
    std::vector<MarkingFinding> findings;
    LayoutMeasures measures;
    // The comparison that gives "better" a meaning: the accepted setup's layout on the same surface and gauge.
    bool hasBaseline = false;
    LayoutMeasures baseline;
    std::string compared;                  // what the baseline is, its depth and age, for people ("" without one)
    std::string recommendation;
};

// What kind of piece the proposals were made for, from its own geometry.
enum class PieceShape { OpenTube, ClosedHidden, ClosedIsolated, Branching, Unsupported, Cup };

// One compatible setup: a start, a termination and a row-start path made from the same field, offered together.
struct ProposalSetup {
    int setup = 0;                         // its markings' setup number (SurfaceMarking::setup)
    std::string recipe;                    // how it was made ("tube/inner", "pole/outer", "band", ...): for refinement and tests
    std::string title, purpose, reason, limitation;   // plain words; every setup has all four
    std::vector<SurfaceMarking> markings;  // offered, on the file's triangles, setup and basis set
    // measured on the prepared reference before any build; "about" where rounded (rows = round(T / h), n = round(L / w))
    int rowsAbout = 0, firstRingStitches = 0, lastRingStitches = 0;
    double omittedVisibleMm2 = 0, omittedHiddenMm2 = 0;
    bool partial = false;                  // leaves visible surface unknitted on purpose
    bool refinable = false;                // a pole start: the deeper trial-relaxation refinement applies
};

struct ProposalSet {
    int piece = -1;
    uint64_t basis = 0;                    // proposalBasis (surface_pieces.h) of the request
    PieceShape shape = PieceShape::Unsupported;
    std::string shapeText;                 // "An open tube with two openings", ...
    std::vector<ProposalSetup> setups;
    std::vector<MarkingFinding> findings;  // refusals and partial guidance, with positions
    std::vector<std::string> log;
    double seconds = 0;
};

} // namespace knit
