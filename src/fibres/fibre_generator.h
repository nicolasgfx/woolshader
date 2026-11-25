#pragma once
// Fibres on knitted yarn, on the CPU: the Fibres step's reference preview and the research tools.
// The scenes grow theirs with gpu::generateFibres; growYarnFibres runs that kernel's light rules
// (fibres.cu growFibre, finishLight) stage by stage, under the same names and in the same order.
//
// Per candidate: place the anchor on its strand, keep it if it faces out of the fabric, is
// unoccluded and survives the thinning; shape a smooth curve (fibre_shape.h) into a polyline that
// starts buried inside the own tube; cut the fibre before anything goes deeper than 20 % of the
// yarn radius into any strand.
#include <cstdint>
#include <vector>

#include "fibres/fibre.h"

namespace knit {

// The budget's view of one part: every candidate's identity and key, with no geometry. A candidate
// is (supporting strand = its position in SurfaceModel::yarns, index on that strand); strand y holds
// FibreFrame::candidates[y] or floor(density L_y). Exactly the candidates the GPU enumerates.
struct CandidateSource {
    uint64_t seed = 0;            // the hash seed of anchors, draws and keys
    std::vector<int> yarnId;      // per strand (hash identity)
    std::vector<int> count;       // per strand
    int strands() const { return static_cast<int>(count.size()); }
    int64_t total() const;
    float key(int strand, int i) const;  // float(Fibre::key) of the candidate if it is accepted
};

uint64_t fibreSeed(const FibreParams& fp, uint64_t pieceSeed);
CandidateSource candidateSource(const SurfaceModel& sm, const FibreParams& fp, const FibreFrame* frame);

FibreSet growYarnFibres(const CurveSet& cs, const SurfaceModel& sm, const PatternParams& pp, const FibreParams& fp,
                        const FibreFrame* frame = nullptr);

} // namespace knit
