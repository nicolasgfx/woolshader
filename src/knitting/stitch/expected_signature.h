#pragma once
// Expected crossing signature of plain stockinette (jersey), front projection.
//
// Source: the knit-stitch construction (TopoKnit, arXiv:2101.04560): a new loop is
// formed by drawing yarn through the previous course's loop of the SAME wale,
// from back to front. Viewed from +N (the technical face, the side the loops are
// drawn towards), the newer yarn of an intermeshed pair:
//   * enters the old loop from behind: at the two OUTER crossings (newer foot over
//     the old loop's side) the newer yarn is BEHIND;
//   * exits in front: at the two INNER crossings (newer legs over the old head)
//     the newer yarn is IN FRONT.
// Each leg/head "twist" (one TopoKnit contact neighbourhood) contributes one outer
// and one inner crossing, so a pair has 4 crossings.
//
// This file is deliberately independent of the stitch template. Viewed from -N
// every sign flips (the purl face). Only RELATIVE depths are constrained here.
#include <cstdlib>

namespace knit::signature {

enum class Depth { Behind, Front };

// Crossings of one intermeshed pair, sorted along the NEWER yarn's parameter,
// expressed as the newer yarn's depth. Palindromic, hence independent of the
// traversal direction (flat knitting alternates direction per course).
inline constexpr Depth kNewerAlongNewer[4] = {Depth::Behind, Depth::Front, Depth::Front, Depth::Behind};

// The same four crossings sorted along the OLDER yarn's parameter (old shoulder,
// old crown, old crown, old shoulder), still expressed as the newer yarn's depth.
inline constexpr Depth kNewerAlongOlder[4] = {Depth::Behind, Depth::Front, Depth::Front, Depth::Behind};

// PURL, derived from the construction, not from a template: a purl loop is drawn
// through its parent from FRONT to back. Seen from +N that is the knit construction
// seen from -N, so the newer yarn enters the old loop from in front (outer crossings:
// newer in FRONT) and leaves it behind (inner crossings: newer BEHIND). Both orders stay
// palindromic. Only the CHILD's face enters: how a pair crosses is fixed by the way the
// child is drawn through; the parent's own face governs the parent's pair with ITS
// parent, not this one. (Mixed patterns - rib, garter - therefore need no further
// table: every pair is a knit pair or a purl pair by its child.)
inline constexpr Depth kPurlNewerAlongNewer[4] = {Depth::Front, Depth::Behind, Depth::Behind, Depth::Front};
inline constexpr Depth kPurlNewerAlongOlder[4] = {Depth::Front, Depth::Behind, Depth::Behind, Depth::Front};

// Shaping, derived from the constructions, not from the templates:
//
// INC = "k2 into the same loop": two adjacent children are drawn back to front
// through ONE untwisted parent loop. Along the children's yarn the strand enters the
// loop from behind once (outer crossing of slot 0's outer leg with the parent's
// head side), comes out over the crown with each of the four legs (inner crossings),
// and leaves behind the other head side (slot 1's outer leg). The two inner legs and
// the sinker between the children stay INSIDE the parent loop and cross nothing.
// Per child, sorted along the child (and, restricted to that child, along the
// parent, which runs the same way):
inline constexpr Depth kIncSlot0[3] = {Depth::Behind, Depth::Front, Depth::Front};
inline constexpr Depth kIncSlot1[3] = {Depth::Front, Depth::Front, Depth::Behind};
//
// DEC = k2tog: one child is drawn back to front through TWO consecutive loops lying
// on top of each other, the later one (p1) in front. Each leg passes through both
// loops, so the child makes with EACH parent the four crossings of a plain pair
// (B, F, F, B). The two parent loops overlap in projection: every crossing between
// p0 and p1 must have p1 in front (their count depends on the outlines; it is even
// for two closed-up loops but is not constrained further).

// Every other pair of yarn portions has ZERO projected crossings: different wales,
// courses r and r+2, and the same yarn with itself (within a span, between
// adjacent spans, and between distant portions).

} // namespace knit::signature
