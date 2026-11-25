// Whether yarn crossings are the crossings of a correct knit: the crossing
// signature of a flat patch seen along one projection, the local signature of every stitch of a
// surface in its own frame, and the attribution of crossings near the ends of a stitch's span.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <map>
#include <sstream>

#include "model/cancel.h"
#include "knitting/checks/checks.h"
#include "knitting/stitch/expected_signature.h"
#include "knitting/stitch/stockinette.h"
#include "model/log.h"

namespace knit {

namespace {
using Clock = std::chrono::steady_clock;
// --- crossing signature (shaping and purl expectations: expected_signature.h)

// The expected crossings of one intermeshed pair: the newer (child's) yarn's depth at each,
// sorted along the newer and along the older yarn, and the sequence's name for messages.
struct ExpectedSequence {
    std::vector<signature::Depth> alongNewer, alongOlder;
    const char* name = "";
};

// Purl shaping stitches have no derived expectation (patterns are K only).
bool hasDerivedExpectation(const StitchRecord& child) {
    return !(child.face == StitchFace::Purl && child.op != StitchOp::K);
}

// The child's op decides the sequence (K: knit or purl face; INC: its slot; DEC: the K sequence
// with each parent). Only called when hasDerivedExpectation(child).
ExpectedSequence expectedSequence(const StitchRecord& child) {
    using namespace signature;
    ExpectedSequence e;
    e.alongNewer.assign(kNewerAlongNewer, kNewerAlongNewer + 4);
    e.alongOlder.assign(kNewerAlongOlder, kNewerAlongOlder + 4);
    e.name = "B,F,F,B";
    if (child.face == StitchFace::Purl) {
        e.alongNewer.assign(kPurlNewerAlongNewer, kPurlNewerAlongNewer + 4);
        e.alongOlder.assign(kPurlNewerAlongOlder, kPurlNewerAlongOlder + 4);
        e.name = "F,B,B,F (purl)";
    }
    if (child.op == StitchOp::Inc) {
        const Depth* inc = child.incSlot == 0 ? kIncSlot0 : kIncSlot1;
        e.alongNewer.assign(inc, inc + 3);
        e.alongOlder.assign(inc, inc + 3);
        e.name = child.incSlot == 0 ? "B,F,F (INC slot 0)" : "F,F,B (INC slot 1)";
    }
    return e;
}

using PairIndex = std::map<std::pair<int, int>, int>;  // (child, parent) -> index into the pairs

// One pair per (child, parent) of the knit topology. For a flat patch these are the
// (wale, course) / (wale, course - 1) pairs in (course, wale) order.
PairIndex listIntermeshedPairs(const KnitTopology& topo, std::vector<PairSignature>& pairs) {
    PairIndex pairIndex;
    for (const StitchRecord& c : topo.stitches)
        for (int q = 0; q < c.nParents; ++q) {
            pairIndex[{c.id, c.parents[q]}] = static_cast<int>(pairs.size());
            PairSignature ps;
            ps.child = c.id;
            ps.parent = c.parents[q];
            ps.wale = c.index;
            ps.olderCourse = topo.stitches[c.parents[q]].turn;
            pairs.push_back(ps);
        }
    return pairIndex;
}

// A same-yarn crossing closer than pi (D + c) along the yarn is a local portion (e.g. a
// buckled yarn overlapping itself in projection) and is not subject to D + c (the
// curvature check covers it).
bool isLocalPortion(const CurveSet& cs, const Crossing& x, double Dc) {
    return x.yarnA == x.yarnB &&
           std::fabs(cs.yarns[x.yarnA].arcLengthAt(x.uB) - cs.yarns[x.yarnA].arcLengthAt(x.uA)) < kPi * Dc;
}

// The depth gap at a crossing is the 3D distance of the two points, so |depth| >= D + c is a
// necessary condition of the clearance requirement - for NONLOCAL pairs only.
void checkDepthGaps(const CurveSet& cs, const CrossingResult& cr, const PatternParams& pp, SignatureResult& res) {
    res.minDepthGapLo = 1e300;
    const double Dc = pp.yarnDiameter + pp.clearance;
    Status gapStatus = Status::Pass;
    for (const Crossing& x : cr.crossings) {
        if (isLocalPortion(cs, x, Dc)) continue;
        const double gapLo = x.depth.containsZero() ? 0.0 : std::min(std::fabs(x.depth.lo), std::fabs(x.depth.hi));
        const double gapHi = x.depth.mag();
        res.minDepthGapLo = std::min(res.minDepthGapLo, gapLo);
        if (gapHi < Dc) gapStatus = combine(gapStatus, Status::Fail);
        else if (gapLo < Dc) gapStatus = combine(gapStatus, Status::Unresolved);
    }
    if (cr.crossings.empty()) res.minDepthGapLo = 0;
    res.depthGapVsClearance = gapStatus;
}

// The pair of topology stitches ta and tb (in either role), -1 if they are not intermeshed.
int intermeshedPairOf(const PairIndex& pairIndex, int ta, int tb) {
    if (ta < 0 || tb < 0) return -1;
    PairIndex::const_iterator it = pairIndex.find({ta, tb});
    if (it == pairIndex.end()) it = pairIndex.find({tb, ta});
    return it == pairIndex.end() ? -1 : it->second;
}

// The DEC child drawn through both stitches ta and tb, -1 if there is none.
int decChildOfParents(const KnitTopology& topo, int ta, int tb) {
    int decChild = -1;
    if (ta >= 0 && tb >= 0 && ta != tb)
        for (int q = 0; q < topo.stitches[ta].nChildren && decChild < 0; ++q) {
            const StitchRecord& c = topo.stitches[topo.stitches[ta].children[q]];
            if (c.op == StitchOp::Dec && (topo.isParent(tb, c.id))) decChild = c.id;
        }
    return decChild;
}

// Two stitches next to each other on one yarn (the same row) with different faces: a rib's knit and
// purl wale. The relaxed rib is two-layered and its loops are wider than the wale spacing, so these
// neighbours overlap in the front projection, the knit loop in front.
bool ribLayerNeighbours(const KnitTopology& topo, int ta, int tb) {
    if (ta < 0 || tb < 0 || ta == tb) return false;
    const StitchRecord &a = topo.stitches[ta], &b = topo.stitches[tb];
    return a.yarn == b.yarn && a.turn == b.turn && std::abs(a.index - b.index) == 1 && a.face != b.face;
}

// Crossings are expected only between a stitch and a loop it is drawn through; that includes
// portions of the SAME yarn (a spiral piece interlocks with itself). The two parents of one DEC
// child may also overlap, with the later parent in front, and so may a rib's
// neighbouring knit and purl stitches, with the knit one in front. Every other crossing is
// unexpected.
void assignCrossingsToPairs(const CurveSet& cs, const KnitTopology& topo, const CrossingResult& cr,
                            const PairIndex& pairIndex, SignatureResult& res) {
    for (size_t k = 0; k < cr.crossings.size(); ++k) {
        const Crossing& x = cr.crossings[k];
        const int ta = cs.yarns[x.yarnA].stitchId(x.stitchA), tb = cs.yarns[x.yarnB].stitchId(x.stitchB);
        const int pair = intermeshedPairOf(pairIndex, ta, tb);
        if (pair >= 0) {
            res.pairs[pair].crossingIdx.push_back(static_cast<int>(k));
            continue;
        }
        const int decChild = decChildOfParents(topo, ta, tb);
        if (decChild >= 0 && x.status == Status::Pass) {
            const StitchRecord& c = topo.stitches[decChild];
            const bool aIsLater = c.parents[1] == ta;
            const bool aFront = x.depth.lo > 0;
            if (aFront != aIsLater) ++res.overlapWrongOrder;
            res.overlap.push_back(static_cast<int>(k));
            continue;
        }
        if (ribLayerNeighbours(topo, ta, tb) && x.status == Status::Pass) {
            const bool aFront = x.depth.lo > 0;
            const bool aKnit = topo.stitches[ta].face == StitchFace::Knit;
            if (aFront != aKnit) ++res.layerWrongOrder;
            res.layerOverlap.push_back(static_cast<int>(k));
            continue;
        }
        res.unexpected.push_back(static_cast<int>(k));
    }
}

// The newer yarn's depth at crossing x, from the sign of z_A - z_B (well defined also when A and B
// are the same yarn, where the front yarn id alone would not tell the two apart). An unsigned or
// unverified crossing makes s UNRESOLVED.
signature::Depth childDepthAt(const CurveSet& cs, const Crossing& x, int child, Status& s) {
    using namespace signature;
    if (x.front < 0 || x.status != Status::Pass) { s = combine(s, Status::Unresolved); return Depth::Behind; }
    const bool aFront = x.depth.lo > 0;
    const bool childIsA = cs.yarns[x.yarnA].stitchId(x.stitchA) == child;
    return aFront == childIsA ? Depth::Front : Depth::Behind;
}

// The pair's crossings against the expected sequence of the child's op: the count, then the
// newer yarn's depth in order along the newer and along the older yarn.
// Sets ps.status and ps.detail; ps.crossingIdx ends sorted along the newer yarn.
void checkPairOrder(const CurveSet& cs, const CrossingResult& cr, PairSignature& ps) {
    const StitchRecord& chRec = cs.topo.stitches[ps.child];
    if (!hasDerivedExpectation(chRec)) {
        ps.status = Status::Unresolved;
        ps.detail = std::string("purl ") + toString(chRec.op) + ": no derived expectation";
        return;
    }
    const ExpectedSequence expected = expectedSequence(chRec);
    if (ps.crossingIdx.size() != expected.alongNewer.size()) {
        ps.status = Status::Fail;
        ps.detail = "expected " + std::to_string(expected.alongNewer.size()) + " crossings (" + toString(chRec.op) +
                    "), found " + std::to_string(ps.crossingIdx.size());
        return;
    }
    // A crossing's parameter on the newer (child's) and on the older (parent's) yarn.
    auto newerU = [&](int k) {
        const Crossing& x = cr.crossings[k];
        return cs.yarns[x.yarnA].stitchId(x.stitchA) == ps.child ? x.uA : x.uB;
    };
    auto olderU = [&](int k) {
        const Crossing& x = cr.crossings[k];
        return cs.yarns[x.yarnA].stitchId(x.stitchA) == ps.child ? x.uB : x.uA;
    };
    Status ps_st = Status::Pass;
    std::vector<int> byNewer = ps.crossingIdx, byOlder = ps.crossingIdx;
    std::sort(byNewer.begin(), byNewer.end(), [&](int a, int b) { return newerU(a) < newerU(b); });
    std::sort(byOlder.begin(), byOlder.end(), [&](int a, int b) { return olderU(a) < olderU(b); });
    for (size_t q = 0; q < expected.alongNewer.size(); ++q) {
        Status s1 = Status::Pass, s2 = Status::Pass;
        const signature::Depth d1 = childDepthAt(cs, cr.crossings[byNewer[q]], ps.child, s1);
        const signature::Depth d2 = childDepthAt(cs, cr.crossings[byOlder[q]], ps.child, s2);
        ps_st = combine(ps_st, combine(s1, s2));
        if (s1 == Status::Pass && d1 != expected.alongNewer[q]) {
            ps_st = Status::Fail;
            ps.detail = std::string("order along newer yarn does not match ") + expected.name;
        }
        if (s2 == Status::Pass && d2 != expected.alongOlder[q]) {
            ps_st = Status::Fail;
            ps.detail = std::string("order along older yarn does not match ") + expected.name;
        }
    }
    ps.crossingIdx = byNewer;
    ps.status = ps_st;
    if (ps_st == Status::Unresolved && ps.detail.empty()) ps.detail = "a crossing sign or ownership is unresolved";
}

// The one-line summary of a signature result (knit_check and the report print it).
std::string describeSignature(const SignatureResult& res, const CrossingResult& cr) {
    std::ostringstream os;
    int okPairs = 0;
    for (const PairSignature& ps : res.pairs) okPairs += ps.status == Status::Pass;
    os << okPairs << "/" << res.pairs.size() << " intermeshed pairs match; " << res.unexpected.size()
       << " unexpected crossings; " << cr.unresolved.size() << " unresolved search regions";
    if (!res.overlap.empty())
        os << "; " << res.overlap.size() << " DEC parent overlap crossings (" << res.overlapWrongOrder
           << " with the wrong layer in front)";
    if (!res.layerOverlap.empty())
        os << "; " << res.layerOverlap.size() << " rib layer overlap crossings (" << res.layerWrongOrder << " with the purl loop in front)";
    return os.str();
}
} // namespace

// The crossing signature of a curve set seen along one front projection:
// list the intermeshed pairs, check the depth gaps, sort every crossing into its pair, then
// compare each pair with its expected sequence.
SignatureResult checkSignature(const CurveSet& cs, const CrossingResult& cr, const PatternParams& pp) {
    SignatureResult res;
    const KnitTopology& topo = cs.topo;
    const PairIndex pairIndex = listIntermeshedPairs(topo, res.pairs);
    checkDepthGaps(cs, cr, pp, res);
    assignCrossingsToPairs(cs, topo, cr, pairIndex, res);

    Status st = Status::Pass;
    for (PairSignature& ps : res.pairs) {
        checkPairOrder(cs, cr, ps);
        st = combine(st, ps.status);
    }
    // An unexpected crossing is a FAIL once its existence and depth sign are verified.
    for (int k : res.unexpected) {
        Status s = cr.crossings[k].status == Status::Pass ? Status::Fail : Status::Unresolved;
        st = combine(st, s);
    }
    if (res.overlapWrongOrder > 0 || res.layerWrongOrder > 0) st = combine(st, Status::Fail);
    if (!cr.unresolved.empty()) st = combine(st, Status::Unresolved);
    res.status = st;
    res.detail = describeSignature(res, cr);
    return res;
}

namespace {

// --- local signature: each stitch checked in its own frame

// The exact sub-B-spline of yarn y over the span range [u0, u1] (integers: whole spans), its
// controls expressed in frame F (origin, U, V, N), so that the front projection along N of the
// mini set is the projection seen from the stitch's own surface normal.
YarnPath extractInFrame(const CurveSet& cs, int y, double u0, double u1, const StitchRecord& F) {
    const BSpline& c = cs.yarns[y].curve;
    int k0 = static_cast<int>(std::floor(u0)), k1 = static_cast<int>(std::ceil(u1)) - 1;
    YarnPath out;
    for (int k = k0; k <= k1 + 3; ++k) {
        Vec3 d = c.ctrl[k] - F.origin;
        out.curve.ctrl.push_back({dot(d, F.U), dot(d, F.V), dot(d, F.N)});
    }
    return out;
}

// Stitches lo..hi along one yarn, from the stitch before `first` to the stitch after `last`
// (where the yarn has them), and [u0, u1], the union of their spans.
struct YarnNeighbourhood {
    int lo = 0, hi = 0;
    double u0 = 1e300, u1 = -1e300;
};

YarnNeighbourhood yarnNeighbourhood(const KnitTopology& topo, const StitchRecord& first, const StitchRecord& last) {
    YarnNeighbourhood nb;
    nb.lo = first.id;
    nb.hi = last.id;
    if (nb.lo > 0 && topo.stitches[nb.lo - 1].yarn == first.yarn) --nb.lo;
    if (nb.hi + 1 < static_cast<int>(topo.stitches.size()) && topo.stitches[nb.hi + 1].yarn == last.yarn) ++nb.hi;
    for (int q = nb.lo; q <= nb.hi; ++q) {
        nb.u0 = std::min(nb.u0, topo.stitches[q].uBegin);
        nb.u1 = std::max(nb.u1, topo.stitches[q].uEnd);
    }
    return nb;
}

// The last parent whose neighbourhood is extracted. Parents on both ends of a cut ring (the
// last before the first): the first parent's neighbourhood only.
const StitchRecord& lastParentInNeighbourhood(const KnitTopology& topo, const StitchRecord& ch) {
    const StitchRecord& p0 = topo.stitches[ch.parents[0]];
    const StitchRecord& last = topo.stitches[ch.parents[ch.nParents - 1]];
    return last.id < p0.id || last.yarn != p0.yarn ? p0 : last;
}

bool startsEarlier(const StitchSpan& x, const StitchSpan& y) { return x.uBegin < y.uBegin; }

// One child's two-yarn mini set in the child's frame: yarn 0 = the child's stitch, yarn 1 = the
// parent neighbourhood nb (on p0's yarn), with a topology in which the child keeps its op, face
// and parents, so checkSignature applies unchanged.
CurveSet childMiniSet(const CurveSet& cs, const KnitTopology& topo, const StitchRecord& ch, const StitchRecord& p0,
                      const YarnNeighbourhood& nb) {
    CurveSet mini;
    mini.kind = CurveKind::FlatPatch;
    YarnPath a = extractInFrame(cs, ch.yarn, ch.uBegin, ch.uEnd, ch);
    a.id = 0;
    YarnPath b = extractInFrame(cs, p0.yarn, nb.u0, nb.u1, ch);
    b.id = 1;
    const double bOff = std::floor(nb.u0);
    StitchRecord mc;
    mc.id = 0;
    mc.yarn = 0;
    mc.uBegin = 0;
    mc.uEnd = a.curve.uMax();
    mc.op = ch.op;
    mc.incSlot = ch.incSlot;
    mc.face = ch.face;
    a.stitches.push_back({ch.index, ch.turn, mc.uBegin, mc.uEnd, 0});
    mini.topo.stitches.push_back(mc);
    for (int q = nb.lo; q <= nb.hi; ++q) {
        const StitchRecord& src = topo.stitches[q];
        StitchRecord m;
        m.id = static_cast<int>(mini.topo.stitches.size());
        m.yarn = 1;
        m.uBegin = src.uBegin - bOff;
        m.uEnd = src.uEnd - bOff;
        m.turn = src.turn;
        m.index = src.index;
        for (int k = 0; k < ch.nParents; ++k)
            if (ch.parents[k] == q) mini.topo.stitches[0].parents[mini.topo.stitches[0].nParents++] = m.id;
        b.stitches.push_back({src.index, src.turn, m.uBegin, m.uEnd, m.id});
        mini.topo.stitches.push_back(m);
    }
    std::sort(b.stitches.begin(), b.stitches.end(), startsEarlier);
    a.finalize();
    b.finalize();
    mini.yarns.push_back(std::move(a));
    mini.yarns.push_back(std::move(b));
    mini.topo.buildIndex({mini.yarns[0].uEnd(), mini.yarns[1].uEnd()});
    return mini;
}

// Appends one reason to a check's detail line, "; "-separated.
void appendDetail(std::string& detail, const std::string& reason) {
    detail += (detail.empty() ? "" : "; ") + reason;
}

// The child's verdict from its mini set's signature. Only crossings involving the child (mini
// yarn 0) speak about its interlock. Crossings of the parent neighbourhood with itself
// (neighbouring heads of the parent row overlapping in this projection, which happens where the
// surface bends within the neighbourhood) are reported but not counted.
LocalPairResult childVerdict(const StitchRecord& ch, const CrossingResult& cr, const SignatureResult& sg) {
    LocalPairResult lp;
    lp.child = ch.id;
    lp.parent = ch.parents[0];
    lp.crossings = 0;
    for (const PairSignature& ps : sg.pairs) lp.crossings += static_cast<int>(ps.crossingIdx.size());
    Status st = Status::Pass;
    for (const PairSignature& ps : sg.pairs) st = combine(st, ps.status);
    if (sg.overlapWrongOrder || sg.layerWrongOrder) st = combine(st, Status::Fail);
    if (!cr.unresolved.empty()) st = combine(st, Status::Unresolved);
    lp.unexpected = 0;
    for (int k : sg.unexpected) {
        const Crossing& x = cr.crossings[k];
        if (x.yarnA != 0 && x.yarnB != 0) { ++lp.neighbourhoodSelf; continue; }
        ++lp.unexpected;
        st = combine(st, x.status == Status::Pass ? Status::Fail : Status::Unresolved);
    }
    lp.overlap = static_cast<int>(sg.overlap.size());
    lp.status = st;
    for (const PairSignature& ps : sg.pairs)
        if (!ps.detail.empty()) appendDetail(lp.detail, ps.detail);
    if (sg.overlapWrongOrder) appendDetail(lp.detail, "DEC layers in the wrong order");
    if (sg.layerWrongOrder) appendDetail(lp.detail, "rib layers in the wrong order");
    if (lp.unexpected) appendDetail(lp.detail, std::to_string(lp.unexpected) + " unexpected crossings");
    return lp;
}

// One child (with at least one parent) against its parent(s) in its own frame: the child's span
// and the parent neighbourhood as a mini set, its crossing search and its signature.
LocalPairResult checkChildInOwnFrame(const CurveSet& cs, const KnitTopology& topo, const PatternParams& pp,
                                     const StitchRecord& ch) {
    const StitchRecord& p0 = topo.stitches[ch.parents[0]];
    const YarnNeighbourhood nb = yarnNeighbourhood(topo, p0, lastParentInNeighbourhood(topo, ch));
    const CurveSet mini = childMiniSet(cs, topo, ch, p0, nb);
    const CrossingResult cr = findCrossings(mini, Allowances::forExtent(mini.extent()));
    const SignatureResult sg = checkSignature(mini, cr, pp);
    return childVerdict(ch, cr, sg);
}

void countVerdict(LocalSignatureResult& res, Status st) {
    res.pass += st == Status::Pass;
    res.fail += st == Status::Fail;
    res.unresolved += st == Status::Unresolved;
    res.status = combine(res.status, st);
}

// Every child is independent of the others, so they run in parallel and are accumulated in
// stitch order: the result is identical for every ExecPolicy.
void checkAllChildren(const CurveSet& cs, const KnitTopology& topo, const PatternParams& pp, const ExecPolicy& exec,
                      LocalSignatureResult& res) {
    const int64_t nSt = static_cast<int64_t>(topo.stitches.size());
    std::vector<LocalPairResult> per(static_cast<size_t>(nSt));
    std::vector<char> has(static_cast<size_t>(nSt), 0);
    parallelFor(nSt, exec.chunk > 0 ? exec.chunk : 16, [&](int64_t begin, int64_t end, int64_t) {
        for (int64_t i = begin; i < end; ++i) {
            if ((topo.stitches[i].id & 63) == 0) cancelCheckpoint();
            const StitchRecord& ch = topo.stitches[i];
            has[i] = ch.nParents > 0;  // the cast-on loops have no parent to check
            if (has[i]) per[i] = checkChildInOwnFrame(cs, topo, pp, ch);
        }
    }, exec);
    for (int64_t i = 0; i < nSt; ++i) {
        if (!has[i]) continue;
        LocalPairResult& lp = per[i];
        res.neighbourhoodSelf += lp.neighbourhoodSelf;
        countVerdict(res, lp.status);
        res.pairs.push_back(std::move(lp));
    }
}

// One CLOSE pass's mini set in the gathered head's frame: yarn 0 = the drawstring over the
// strand range, yarn 1 = the head and its yarn neighbours nb. headMini is the head's stitch id in
// the mini topology.
struct StrandRange {
    int yarn;
    double s0, s1;
};

struct CloseMiniSet {
    CurveSet mini;
    int headMini = -1;
};

CloseMiniSet closeMiniSet(const CurveSet& cs, const StrandRange& strand, const StitchRecord& head, const YarnNeighbourhood& nb) {
    const KnitTopology& topo = cs.topo;
    CloseMiniSet out;
    CurveSet& mini = out.mini;
    YarnPath a = extractInFrame(cs, strand.yarn, strand.s0, strand.s1, head);
    a.id = 0;
    YarnPath b = extractInFrame(cs, head.yarn, nb.u0, nb.u1, head);
    b.id = 1;
    const double bOff = std::floor(nb.u0);
    a.stitches.push_back({0, 0, 0.0, a.curve.uMax(), 0});
    StitchRecord ms0;
    ms0.id = 0;
    ms0.uEnd = a.curve.uMax();
    mini.topo.stitches.push_back(ms0);
    for (int q = nb.lo; q <= nb.hi; ++q) {
        StitchRecord m;
        m.id = static_cast<int>(mini.topo.stitches.size());
        m.yarn = 1;
        m.uBegin = topo.stitches[q].uBegin - bOff;
        m.uEnd = topo.stitches[q].uEnd - bOff;
        if (q == head.id) out.headMini = m.id;
        b.stitches.push_back({0, 0, m.uBegin, m.uEnd, m.id});
        mini.topo.stitches.push_back(m);
    }
    std::sort(b.stitches.begin(), b.stitches.end(), startsEarlier);
    a.finalize();
    b.finalize();
    mini.yarns.push_back(std::move(a));
    mini.yarns.push_back(std::move(b));
    mini.topo.buildIndex({mini.yarns[0].uEnd(), mini.yarns[1].uEnd()});
    return out;
}

// The CLOSE signature: exactly two crossings with the gathered loop, B then F
// along the strand on even loops, F then B on odd ones. withHead: (u on the strand, strand in
// front), sorted along the strand. Returns the mismatch, empty when the pass matches.
std::string closeOrderMismatch(const std::vector<std::pair<double, bool>>& withHead, bool even) {
    if (withHead.size() != 2)
        return "CLOSE: expected 2 crossings with the gathered loop, found " + std::to_string(withHead.size());
    if (withHead[0].second == even || withHead[1].second != even)
        return std::string("CLOSE: order along the strand does not match ") + (even ? "B,F" : "F,B");
    return "";
}

// The drawstring's pass through gathered loop k in the head's frame: it
// crosses that loop twice in the expected order and nothing else of the neighbourhood.
LocalPairResult checkCloseLoopInOwnFrame(const CurveSet& cs, const KnitTopology& topo, const TerminalRecord& tr, size_t k) {
    const StitchRecord& head = topo.stitches[tr.stitches[k]];
    const YarnNeighbourhood nb = yarnNeighbourhood(topo, head, head);
    const double s0 = std::max(0.0, std::floor(tr.segments[k].first));
    const double s1 = std::min(cs.yarns[tr.yarn].uEnd(), std::ceil(tr.segments[k].second));
    const CloseMiniSet m = closeMiniSet(cs, {tr.yarn, s0, s1}, head, nb);
    const CrossingResult cr = findCrossings(m.mini, Allowances::forExtent(m.mini.extent()));
    LocalPairResult lp;
    lp.child = head.id;
    lp.parent = -1;
    lp.terminal = 1;
    std::vector<std::pair<double, bool>> withHead;  // (u on the strand, strand in front)
    Status st = cr.unresolved.empty() ? Status::Pass : Status::Unresolved;
    for (const Crossing& x : cr.crossings) {
        if (x.status != Status::Pass) { st = combine(st, Status::Unresolved); continue; }
        if (x.yarnA != x.yarnB) {
            const bool strandA = x.yarnA == 0;
            // Only crossings inside this loop's nominal strand range count;
            // the extracted span range reaches into the neighbours' passes.
            const double us = (strandA ? x.uA : x.uB) + s0;
            if (us < tr.segments[k].first || us > tr.segments[k].second) continue;
            const int other = m.mini.topo.stitchAt(1, strandA ? x.uB : x.uA);
            const bool strandFront = strandA ? x.depth.lo > 0 : x.depth.hi < 0;
            if (other == m.headMini) withHead.push_back({strandA ? x.uA : x.uB, strandFront});
            else ++lp.unexpected;
        } else {
            ++lp.unexpected;  // the strand or the heads crossing themselves
        }
    }
    std::sort(withHead.begin(), withHead.end());
    lp.crossings = static_cast<int>(withHead.size());
    if (lp.unexpected > 0) st = Status::Fail;
    const std::string mismatch = closeOrderMismatch(withHead, k % 2 == 0);
    if (!mismatch.empty()) {
        st = Status::Fail;
        lp.detail = mismatch;
    }
    if (lp.unexpected) appendDetail(lp.detail, std::to_string(lp.unexpected) + " unexpected crossings");
    lp.status = st;
    return lp;
}

// Every CLOSE drawstring pass, in loop order; returns the number of loops checked.
int checkCloseLoops(const CurveSet& cs, const KnitTopology& topo, LocalSignatureResult& res) {
    int closeChecks = 0;
    for (const TerminalRecord& tr : topo.terminals) {
        if (tr.op != TerminalOp::Close || tr.segments.size() != tr.stitches.size()) continue;
        for (size_t k = 0; k < tr.stitches.size(); ++k) {
            const LocalPairResult lp = checkCloseLoopInOwnFrame(cs, topo, tr, k);
            countVerdict(res, lp.status);
            res.pairs.push_back(lp);
            ++closeChecks;
        }
    }
    return closeChecks;
}

} // namespace

// The local signature: every child against its parents in its own frame, every
// CLOSE drawstring pass against its gathered loop, then the crossing attribution correction.
LocalSignatureResult checkLocalSignature(const CurveSet& cs, const PatternParams& pp, const ExecPolicy& exec, double window) {
    auto t0 = Clock::now();
    LocalSignatureResult res;
    const KnitTopology& topo = cs.topo;
    checkAllChildren(cs, topo, pp, exec, res);
    const int closeChecks = checkCloseLoops(cs, topo, res);
    char buf[240];
    std::snprintf(buf, sizeof buf, "%d/%zu stitches match in their own frame (%d FAIL, %d UNRESOLVED)%s", res.pass,
                  res.pairs.size(), res.fail, res.unresolved,
                  closeChecks ? (", incl. " + std::to_string(closeChecks) + " CLOSE loops").c_str() : "");
    res.detail = buf;
    applyCrossingAttribution(cs, res, window);
    res.ms = msSince(t0);
    return res;
}

// ------------------------------------------------------- crossing attribution ----

CrossingOwner attributeCrossing(const KnitTopology& topo, int yarn, double u, int parentStitch, double window) {
    CrossingOwner r;
    const int own = topo.stitchAt(yarn, u);
    if (own < 0) return r;  // outside every stitch span: not attributed (the original rule)
    if (topo.isParent(parentStitch, own)) {  // the own span first (the original rule)
        r.owner = own;
        return r;
    }
    // window < the shortest stitch span: only the stitches at u -+ window can be within reach
    int cand[2], n = 0;
    for (const double v : {u - window, u + window}) {
        const int c = topo.stitchAt(yarn, v);
        if (c < 0 || c == own || (n > 0 && cand[0] == c)) continue;
        const StitchRecord& s = topo.stitches[c];
        if (topo.isParent(parentStitch, c) && u >= s.uBegin - window && u < s.uEnd + window) cand[n++] = c;
    }
    if (n == 1) r.owner = cand[0];
    else if (n > 1) r.ambiguous = true;
    return r;
}

CrossingOwner attributeCloseCrossing(const TerminalRecord& tr, int headStitch, double u, double window) {
    CrossingOwner r;
    bool owned = false;  // inside some loop's nominal range (outside every range: not attributed)
    for (size_t q = 0; q < tr.stitches.size() && q < tr.segments.size(); ++q)
        owned = owned || (u >= tr.segments[q].first && u <= tr.segments[q].second);
    if (!owned) return r;
    int n = 0, cand = -1;
    for (size_t q = 0; q < tr.stitches.size() && q < tr.segments.size(); ++q) {
        if (tr.stitches[q] != headStitch) continue;
        if (u >= tr.segments[q].first && u <= tr.segments[q].second) {
            r.owner = static_cast<int>(q);
            return r;
        }
        if (u >= tr.segments[q].first - window && u <= tr.segments[q].second + window) {
            ++n;
            cand = static_cast<int>(q);
        }
    }
    if (n == 1) r.owner = cand;
    else if (n > 1) r.ambiguous = true;
    return r;
}

namespace {

// A range of B-spline parameters on one yarn.
struct YarnRange {
    int yarn;
    double u0, u1;
};

// One crossing of two ranges in a stitch's frame, in global parameters.
struct FramedCrossing {
    double ua, ub;
    bool resolved, aFront;
};

struct FramedCrossings {
    std::vector<FramedCrossing> crossings;  // between the two ranges
    size_t unresolvedRegions = 0;
    int selfCrossings = 0;                  // of a range with itself: not listed
};

// The range widened to whole spans, its controls in frame F; off: the first span.
YarnPath framedRange(const CurveSet& cs, const StitchRecord& F, const YarnRange& range, double& off) {
    const BSpline& c = cs.yarns[range.yarn].curve;
    const int k0 = std::max(0, static_cast<int>(std::floor(range.u0)));
    const int k1 = std::min(static_cast<int>(c.ctrl.size()) - 4, static_cast<int>(std::ceil(range.u1)) - 1);
    off = k0;
    YarnPath out;
    for (int k = k0; k <= k1 + 3; ++k) {
        const Vec3 d = c.ctrl[k] - F.origin;
        out.curve.ctrl.push_back({dot(d, F.U), dot(d, F.V), dot(d, F.N)});
    }
    out.stitches.push_back({0, 0, 0.0, out.curve.uMax(), 0});
    out.finalize();
    return out;
}

// Crossings of two ranges projected along F.N at F.origin: global parameters, depth order, status.
FramedCrossings framedCrossings(const CurveSet& cs, const StitchRecord& F, const YarnRange& rangeA, const YarnRange& rangeB) {
    double offA, offB;
    CurveSet mini;
    mini.kind = CurveKind::FlatPatch;
    YarnPath a = framedRange(cs, F, rangeA, offA), b = framedRange(cs, F, rangeB, offB);
    a.id = 0;
    b.id = 1;
    mini.yarns.push_back(std::move(a));
    mini.yarns.push_back(std::move(b));
    const CrossingResult cr = findCrossings(mini, Allowances::forExtent(mini.extent()));
    FramedCrossings out;
    out.unresolvedRegions = cr.unresolved.size();
    for (const Crossing& x : cr.crossings) {
        if (x.yarnA == x.yarnB) {
            ++out.selfCrossings;
            continue;
        }
        const bool aIsA = x.yarnA == 0;
        out.crossings.push_back({(aIsA ? x.uA : x.uB) + offA, (aIsA ? x.uB : x.uA) + offB, x.status == Status::Pass, aIsA ? x.depth.lo > 0 : x.depth.hi < 0});
    }
    return out;
}

struct Claim {
    int yarnA, yarnB;
    double ua, ub;
    int check;  // index into the pairs
};

// One check's re-evaluation: the window and where its claims go.
struct AttributionRun {
    double window;
    int check;
    std::vector<Claim>& claims;
};

// The child's framed crossings sorted by the attribution rule: the ones it owns,
// per parent; the ones in its own span that no stitch can take (unexpected); and the ambiguous ones
// that concern it. status: UNRESOLVED for an owned crossing that is not verified, FAIL for a verified
// unexpected one. Every owned crossing is also recorded as a claim of the check.
struct ChildCrossings {
    std::vector<std::vector<const FramedCrossing*>> perParent;  // indexed like ch.parents
    int unexpected = 0, ambiguous = 0;
    Status status = Status::Pass;
};

ChildCrossings attributeChildCrossings(const KnitTopology& topo, const StitchRecord& ch, int parentYarn,
                                       const std::vector<FramedCrossing>& xs, const AttributionRun& run) {
    const double w = run.window;
    ChildCrossings out;
    out.perParent.resize(ch.nParents);
    for (const FramedCrossing& x : xs) {
        const int sp = topo.stitchAt(parentYarn, x.ub);
        if (sp < 0) continue;
        const CrossingOwner o = attributeCrossing(topo, ch.yarn, x.ua, sp, w);
        const bool inOwn = x.ua >= ch.uBegin && x.ua < ch.uEnd;
        if (o.ambiguous) {
            // this child is one of the candidates, or the crossing lies in its own span
            if (inOwn || (topo.isParent(sp, ch.id) && x.ua >= ch.uBegin - w && x.ua < ch.uEnd + w)) ++out.ambiguous;
            continue;
        }
        if (o.owner == ch.id) {
            if (!x.resolved) out.status = combine(out.status, Status::Unresolved);
            for (int q = 0; q < ch.nParents; ++q)
                if (ch.parents[q] == sp) out.perParent[q].push_back(&x);
            run.claims.push_back({ch.yarn, parentYarn, x.ua, x.ub, run.check});
        } else if (o.owner < 0 && inOwn) {
            ++out.unexpected;  // in the child's own span, no stitch can take it
            out.status = combine(out.status, x.resolved ? Status::Fail : Status::Unresolved);
        }
    }
    return out;
}

// One parent's attributed crossings against the expected sequence of the child's op, as
// checkPairOrder does on the mini set: the count, then the child's depth in order along the
// child's and along the parent's yarn. Unverified crossings are skipped (already UNRESOLVED).
// Appends the first mismatch to detail and returns false on one.
bool parentCrossingsMatch(const StitchRecord& ch, std::vector<const FramedCrossing*> v, std::string& detail) {
    using namespace signature;
    const ExpectedSequence expected = expectedSequence(ch);
    if (v.size() != expected.alongNewer.size()) {
        appendDetail(detail, std::string("expected ") + std::to_string(expected.alongNewer.size()) + " crossings (" +
                                 toString(ch.op) + "), found " + std::to_string(v.size()));
        return false;
    }
    std::vector<const FramedCrossing*> byOlder = v;
    std::sort(v.begin(), v.end(), [](const FramedCrossing* a, const FramedCrossing* b) { return a->ua < b->ua; });
    std::sort(byOlder.begin(), byOlder.end(), [](const FramedCrossing* a, const FramedCrossing* b) { return a->ub < b->ub; });
    for (size_t k = 0; k < v.size(); ++k) {
        if (!v[k]->resolved || !byOlder[k]->resolved) continue;
        if ((v[k]->aFront ? Depth::Front : Depth::Behind) != expected.alongNewer[k]) {
            appendDetail(detail, std::string("order along newer yarn does not match ") + expected.name);
            return false;
        }
        if ((byOlder[k]->aFront ? Depth::Front : Depth::Behind) != expected.alongOlder[k]) {
            appendDetail(detail, std::string("order along older yarn does not match ") + expected.name);
            return false;
        }
    }
    return true;
}

// One child's pairs (K, INC, DEC; knit or purl face) with the attribution rule. Frame, parent
// neighbourhood (parents and one stitch on either side) and expected sequences as
// checkLocalSignature; the child's yarn widened by the window on either side.
LocalPairResult reevaluatePair(const CurveSet& cs, const LocalPairResult& prod, const AttributionRun& run) {
    const KnitTopology& topo = cs.topo;
    const double w = run.window;
    const StitchRecord& ch = topo.stitches[prod.child];
    LocalPairResult lp = prod;
    lp.detail.clear();
    lp.unexpected = 0;
    lp.crossings = 0;
    if (!hasDerivedExpectation(ch)) return prod;  // unchanged
    const StitchRecord& p0 = topo.stitches[ch.parents[0]];
    const StitchRecord& p1 = topo.stitches[ch.parents[ch.nParents - 1]];
    const YarnNeighbourhood nb = yarnNeighbourhood(topo, p0, p1);
    const FramedCrossings framed = framedCrossings(cs, ch, {ch.yarn, ch.uBegin - w, ch.uEnd + w}, {p0.yarn, nb.u0, nb.u1});
    const ChildCrossings owned = attributeChildCrossings(topo, ch, p0.yarn, framed.crossings, run);
    lp.unexpected = owned.unexpected;
    Status st = framed.unresolvedRegions ? Status::Unresolved : Status::Pass;
    st = combine(st, owned.status);
    if (owned.ambiguous) st = combine(st, Status::Unresolved);
    for (int q = 0; q < ch.nParents; ++q) {
        lp.crossings += static_cast<int>(owned.perParent[q].size());
        if (!parentCrossingsMatch(ch, owned.perParent[q], lp.detail)) st = Status::Fail;
    }
    // The DEC parents' layer order does not depend on the attribution: kept.
    if (prod.detail.find("DEC layers in the wrong order") != std::string::npos) {
        st = Status::Fail;
        appendDetail(lp.detail, "DEC layers in the wrong order");
    }
    if (lp.unexpected) appendDetail(lp.detail, std::to_string(lp.unexpected) + " unexpected crossings");
    if (owned.ambiguous) appendDetail(lp.detail, std::to_string(owned.ambiguous) + " crossings of ambiguous ownership");
    lp.status = st;
    return lp;
}

// One CLOSE loop with the attribution rule on the drawstring's nominal segments. Frame and
// neighbourhood as checkLocalSignature; the strand widened by the window; crossings of the strand
// or the heads with themselves counted on the original strand range (as before).
LocalPairResult reevaluateClose(const CurveSet& cs, const TerminalRecord& tr, size_t k, const LocalPairResult& prod, const AttributionRun& run) {
    const KnitTopology& topo = cs.topo;
    const double w = run.window;
    const StitchRecord& h = topo.stitches[tr.stitches[k]];
    const YarnNeighbourhood nb = yarnNeighbourhood(topo, h, h);
    const double hu0 = topo.stitches[nb.lo].uBegin, hu1 = topo.stitches[nb.hi].uEnd;
    const double s0 = std::max(0.0, std::floor(tr.segments[k].first)), s1 = std::min(cs.yarns[tr.yarn].uEnd(), std::ceil(tr.segments[k].second));
    const YarnRange heads{h.yarn, hu0, hu1};
    const FramedCrossings nominal = framedCrossings(cs, h, {tr.yarn, s0, s1}, heads);
    const FramedCrossings widened =
        framedCrossings(cs, h, {tr.yarn, std::max(0.0, tr.segments[k].first - w), std::min(cs.yarns[tr.yarn].uEnd(), tr.segments[k].second + w)}, heads);
    LocalPairResult lp = prod;
    Status st = widened.unresolvedRegions || nominal.unresolvedRegions ? Status::Unresolved : Status::Pass;
    std::vector<std::pair<double, bool>> withHead;
    int unexpected = nominal.selfCrossings, ambiguous = 0;
    for (const FramedCrossing& x : widened.crossings) {
        const int j = topo.stitchAt(h.yarn, x.ub);
        const CrossingOwner o = attributeCloseCrossing(tr, j, x.ua, w);
        const bool inK = x.ua >= tr.segments[k].first && x.ua <= tr.segments[k].second;
        if (o.ambiguous) {
            if (inK || j == h.id) ++ambiguous;
            continue;
        }
        if (o.owner == static_cast<int>(k)) {
            if (!x.resolved) st = combine(st, Status::Unresolved);
            withHead.push_back({x.ua, x.aFront});
            run.claims.push_back({tr.yarn, h.yarn, x.ua, x.ub, run.check});
        } else if (o.owner < 0 && inK) {
            ++unexpected;
        }
    }
    std::sort(withHead.begin(), withHead.end());
    lp.detail.clear();
    if (unexpected) st = Status::Fail;
    if (ambiguous) st = combine(st, Status::Unresolved);
    const std::string mismatch = closeOrderMismatch(withHead, k % 2 == 0);
    if (!mismatch.empty()) {
        st = Status::Fail;
        lp.detail = mismatch;
    }
    if (unexpected) appendDetail(lp.detail, std::to_string(unexpected) + " unexpected crossings");
    if (ambiguous) appendDetail(lp.detail, std::to_string(ambiguous) + " crossings of ambiguous ownership");
    lp.crossings = static_cast<int>(withHead.size());
    lp.unexpected = unexpected;
    lp.status = st;
    return lp;
}

// The CLOSE terminal whose loops are checked (the last one with a strand range per gathered
// loop), nullptr if the piece has none.
const TerminalRecord* checkedCloseTerminal(const KnitTopology& topo) {
    const TerminalRecord* close = nullptr;
    for (const TerminalRecord& tr : topo.terminals)
        if (tr.op == TerminalOp::Close && tr.segments.size() == tr.stitches.size() && !tr.segments.empty()) close = &tr;
    return close;
}

// The checks the attribution can change, as flags per stitch and per CLOSE loop.
// A crossing moves to a stitch only from a yarn neighbour's span (window < a stitch span, never
// from outside every span), where it was unexpected, so that neighbour failed. CLOSE likewise:
// failing loops and their neighbouring loops.
struct ChecksToRedo {
    std::vector<char> stitch, closeLoop;
};

ChecksToRedo checksThatCanChange(const KnitTopology& topo, const LocalSignatureResult& sig, const TerminalRecord* close) {
    const int N = static_cast<int>(topo.stitches.size());
    ChecksToRedo redo;
    redo.stitch.assign(N, 0);
    if (close) redo.closeLoop.assign(close->stitches.size(), 0);
    for (const LocalPairResult& lp : sig.pairs) {
        if (lp.status == Status::Pass || lp.child < 0) continue;
        if (!lp.terminal) {
            for (int d = -1; d <= 1; ++d) {
                const int c = lp.child + d;
                if (c >= 0 && c < N && topo.stitches[c].yarn == topo.stitches[lp.child].yarn) redo.stitch[c] = 1;
            }
        } else if (close) {
            for (size_t k = 0; k < close->stitches.size(); ++k)
                if (close->stitches[k] == lp.child)
                    for (int d = -1; d <= 1; ++d)
                        if (static_cast<int>(k) + d >= 0 && static_cast<int>(k) + d < static_cast<int>(redo.closeLoop.size()))
                            redo.closeLoop[k + d] = 1;
        }
    }
    return redo;
}

// The CLOSE loop that gathers head stitch `head`; close.stitches.size() if none does.
size_t closeLoopOf(const TerminalRecord& close, int head) {
    size_t k = 0;
    while (k < close.stitches.size() && close.stitches[k] != head) ++k;
    return k;
}

// The marked checks re-evaluated with the attribution rule, in check order; every crossing a
// re-evaluated check takes is recorded as a claim.
struct Reevaluation {
    std::vector<Claim> claims;
    int reevaluated = 0, changed = 0;
};

Reevaluation reevaluateMarkedChecks(const CurveSet& cs, const TerminalRecord* close, const ChecksToRedo& redo, double window, LocalSignatureResult& sig) {
    const KnitTopology& topo = cs.topo;
    Reevaluation out;
    for (size_t i = 0; i < sig.pairs.size(); ++i) {
        LocalPairResult& lp = sig.pairs[i];
        if (lp.child < 0) continue;
        if ((i & 63) == 0) cancelCheckpoint();  // a widely failing knit re-evaluates many checks
        const AttributionRun run{window, static_cast<int>(i), out.claims};
        LocalPairResult nl;
        if (!lp.terminal) {
            if (!redo.stitch[lp.child] || topo.stitches[lp.child].nParents == 0) continue;
            nl = reevaluatePair(cs, lp, run);
        } else {
            if (!close) continue;
            const size_t k = closeLoopOf(*close, lp.child);
            if (k == close->stitches.size() || !redo.closeLoop[k]) continue;
            nl = reevaluateClose(cs, *close, k, lp, run);
        }
        ++out.reevaluated;
        out.changed += nl.status != lp.status;
        lp = nl;
    }
    return out;
}

// Each crossing counted once: the same crossing (same yarns, both parameters within 0.05
// spans: one crossing seen in two frames) claimed by two checks makes both UNRESOLVED.
void markDoubleClaims(std::vector<Claim>& claims, LocalSignatureResult& sig) {
    std::sort(claims.begin(), claims.end(), [](const Claim& a, const Claim& b) { return a.ua < b.ua; });
    for (size_t i = 0; i < claims.size(); ++i) {
        if ((i & 1023) == 0) cancelCheckpoint();
        for (size_t j = i + 1; j < claims.size() && claims[j].ua - claims[i].ua < 0.05; ++j)
            if (claims[j].check != claims[i].check && claims[j].yarnA == claims[i].yarnA && claims[j].yarnB == claims[i].yarnB &&
                std::fabs(claims[j].ub - claims[i].ub) < 0.05)
                for (int c : {claims[i].check, claims[j].check}) {
                    LocalPairResult& lp = sig.pairs[c];
                    if (lp.detail.find("claimed by two checks") == std::string::npos)
                        appendDetail(lp.detail, "a crossing claimed by two checks");
                    lp.status = combine(lp.status, Status::Unresolved);
                }
    }
}

// The totals recounted from the corrected checks; returns the number of CLOSE checks.
int recountVerdicts(LocalSignatureResult& sig) {
    sig.pass = sig.fail = sig.unresolved = 0;
    sig.status = Status::Pass;
    int closeChecks = 0;
    for (const LocalPairResult& lp : sig.pairs) {
        countVerdict(sig, lp.status);
        closeChecks += lp.terminal;
    }
    return closeChecks;
}

} // namespace

// The crossing attribution correction: mark the checks it can change,
// re-evaluate them with the rule, make crossings claimed twice UNRESOLVED, recount.
void applyCrossingAttribution(const CurveSet& cs, LocalSignatureResult& sig, double window) {
    if (window <= 0 || sig.pairs.empty()) return;
    const KnitTopology& topo = cs.topo;
    const TerminalRecord* close = checkedCloseTerminal(topo);
    const ChecksToRedo redo = checksThatCanChange(topo, sig, close);
    Reevaluation ev = reevaluateMarkedChecks(cs, close, redo, window, sig);
    markDoubleClaims(ev.claims, sig);
    const int closeChecks = recountVerdicts(sig);
    sig.attributionReevaluated = ev.reevaluated;
    char buf[320];
    std::snprintf(buf, sizeof buf, "%d/%zu stitches match in their own frame (%d FAIL, %d UNRESOLVED)%s; crossing attribution window %g spans (%d re-evaluated, %d changed)",
                  sig.pass, sig.pairs.size(), sig.fail, sig.unresolved, closeChecks ? (", incl. " + std::to_string(closeChecks) + " CLOSE loops").c_str() : "",
                  window, ev.reevaluated, ev.changed);
    sig.detail = buf;
}

} // namespace knit
