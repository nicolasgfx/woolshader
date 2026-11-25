#pragma once
// Production's stitch records (model/topology.h) for a yarn embedded in a typed stitch mesh (yarn_embed.h) and converted
// to production curves (yarn_curves.h): what the yarn surface, the fibres, the stitch colours and the export read from a
// knit (counts, the stitch owning a parameter, frames, parents). Derived from the stitch mesh, never assumed.
//
// A loop is a face's own row yarn through one of its top edges (Yuksel et al. 2012 section 4: a face holds the top of one
// loop and the bottom of the next). Along the strand a face's row pieces and the heads or caps between them are
// contiguous: row 0 (from the W- port) -> the head above, through the top edge's two ports -> row 1 ... -> the last row
// (to the W+ port). A face with n top edges has n + 1 row pieces and makes n loops; loop j runs from row j's middle (row 0's
// start for j = 0) to row j + 1's middle (the last row's end for j = n - 1). So K and d12k make one loop, kp two, and the
// ownership on the strand is half-open and gap-free inside a face.
//
// Parents: the loops a loop is drawn through are the heads in its face (one per bottom edge, in head order): each head
// belongs to the loop of the face below whose top edge it passes. K: 1 parent; each loop of a kp: the same 1 (op Inc, slots
// 0 and 1); d12k: 2 (op Dec); a cast-on face and a pole-start face: none (the first turn, as a flat patch's first course).
// The total of parent links must equal the construction check's expected interlocks (YarnCheck::interlocksExpected): the
// records and the yarn check count the same heads.
//
// Frame: origin the face centroid, N the reference normal there, U the loop's chord (W- to W+, knitting order) projected
// off N, V = N x U. The local signature check is NOT applicable to these records: its expected crossings are production's
// 16-control template, not these templates; the guided route checks the construction with checkYarn instead.
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/yarn_path.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/yarn_curves.h"
#include "stitchmesh/yarn_embed.h"

namespace knit {

struct GuidedRecordCounts {
    int loops = 0, increases = 0, decreases = 0, firstTurn = 0, parentLinks = 0;
};

// Fills curves.topo (one record per loop, ids = indices) and every yarn's stitch spans, then finalizes the curves again.
// s and faces as embedYarn took them; e: the embedding whose pieces cv converted; curves: cv.curves or the same curves moved
// (same parameters); expectedInterlocks: YarnCheck::interlocksExpected of e. False with the reason when the pieces do not
// form loops as described (a face without row pieces, overlapping ownership, a head through no known top edge, a count
// that disagrees).
[[nodiscard]] bool guidedStitchRecords(const TypedStitchMesh& s, const std::vector<std::vector<int>>& faces, const YarnEmbedding& e, const YarnConversion& cv,
                                       const TriMesh& reference, int expectedInterlocks, CurveSet& curves, GuidedRecordCounts& counts, std::string* err);

} // namespace knit
