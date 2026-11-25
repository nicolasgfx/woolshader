#pragma once
// Proposals for one piece of an imported file:
// a few compatible setups, each a start, a termination and a row-start path made from one field, from the piece's own
// geometry only (its border loops, its topology, the closed pieces of the same component that hide it). No coordinates,
// axes, names or manifest data enter. One backend for the app (jobs/marking_session.h) and the lab (stitchmesh_lab
// propose-marks).
//
//   open tube        two border loops: the loops themselves are the cast-on and the bind-off, in either direction; the
//                    row-start path by steepest descent, on the inner and the outer side of a bend (one path when straight)
//   closed, hidden   a start point at the free end of the geodesic diameter (rule 1), the finish level where the piece
//                    first lies entirely inside another piece plus half a row; inner and outer paths; a band; and, when
//                    that finish ring is narrow, a finish where the ring still has 8 stitches (the visible area it costs)
//   closed, alone    nothing hides an end: every setup leaves a visible cap, disclosed with its area (no closed finish)
//   branching        rows split: refused at the split, with a partial setup up to it
//
// Rule 2 (where an unconstrained relaxation of a setup's layout settles the pole) is a deeper calculation, run only on
// request (refineProposalStart): one layout build with relaxation, never final yarn settling, never a verdict on rule 1.
#include <atomic>
#include <cstdint>
#include <string>

#include "model/marking_review.h"
#include "model/surface_marking.h"
#include "model/surface_pieces.h"
#include "stitchmesh/prepared_reference.h"

namespace knit {

struct ProposalRequest {
    const SurfacePieces* pieces = nullptr;
    int piece = -1;
    uint64_t fileDigest = 0;              // the markings set's digest (meshDigest of the file as read)
    std::string source;                   // the file's path (the lab's markings file names it)
    std::string componentName;            // the piece's component, for the lab's report
    double mmPerUnit = 1, w = 5.5, h = 3.4;
    ReferencePreparation preparation = ReferencePreparation::Facets;
    std::atomic<int>* step = nullptr;     // advanced after every step, for progress (kProposeSteps, kRefineSteps); may be null
};
constexpr int kProposeSteps = 4;   // the prepared piece, its shape, the fields, the setups
constexpr int kRefineSteps = 2;    // the trial build, the setup re-derived

// The setups and findings for one piece (out.basis = proposalBasis of the request). False only for an unusable request;
// a piece that cannot be proposed for gets a refusal finding and no setup (or a partial one).
[[nodiscard]] bool proposeMarkings(const ProposalRequest& r, ProposalSet& out, std::string* err);

// Rule 2 on a setup with a start point: the setup's layout built and relaxed once with nothing holding the pole, the
// start moved to where the pole settled, the finish and path derived again in that field by the setup's own recipe.
[[nodiscard]] bool refineProposalStart(const ProposalRequest& r, const ProposalSetup& from, ProposalSetup& out, std::string* err);

// The lab's nine-marking proposals file of a closed piece hidden at one end, exactly as stitchmesh_lab propose-marks wrote it
// before assisted stockinette (ids 1-9, alternative groups 1-3, rule 2 offered first), and its text report. False with a
// message when no level of the piece lies inside another piece, or the trial build does not relax.
[[nodiscard]] bool labProposalFile(const ProposalRequest& r, MarkingSet& out, std::string& report, std::string* err);

} // namespace knit
