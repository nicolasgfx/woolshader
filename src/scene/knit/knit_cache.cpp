// The per-piece .knitcache file: one knitted piece with its attempts, and the curve encodings (knit_cache.h).
#include "scene/knit/knit_cache.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <sstream>


namespace knit {

namespace {

const char kMagic[] = "knit-cache v2\n";  // v2: PatternParams::ribLayerOffset, ribPairSpread, garterWave

void putStitch(ByteWriter& w, const StitchRecord& r) {
    w.i64(r.id); w.i64(r.yarn); w.f64(r.uBegin); w.f64(r.uEnd); w.i64(r.turn); w.i64(r.index); w.i64(static_cast<int>(r.op));
    w.i64(r.parents[0]); w.i64(r.parents[1]); w.i64(r.nParents); w.i64(r.children[0]); w.i64(r.children[1]); w.i64(r.nChildren);
    w.i64(r.incSlot); w.i64(static_cast<int>(r.face));
    w.i64(r.lattice.chart); w.i64(r.lattice.i); w.i64(r.lattice.r); w.i64(r.lattice.run);
    w.vec3(r.origin); w.vec3(r.U); w.vec3(r.V); w.vec3(r.N);
}
StitchRecord getStitch(ByteReader& q) {
    StitchRecord r;
    r.id = int(q.i64()); r.yarn = int(q.i64()); r.uBegin = q.f64(); r.uEnd = q.f64(); r.turn = int(q.i64()); r.index = int(q.i64());
    r.op = static_cast<StitchOp>(q.i64());
    r.parents[0] = int(q.i64()); r.parents[1] = int(q.i64()); r.nParents = int(q.i64());
    r.children[0] = int(q.i64()); r.children[1] = int(q.i64()); r.nChildren = int(q.i64());
    r.incSlot = int(q.i64()); r.face = static_cast<StitchFace>(q.i64());
    r.lattice.chart = int(q.i64()); r.lattice.i = int(q.i64()); r.lattice.r = int(q.i64()); r.lattice.run = int(q.i64());
    r.origin = q.vec3(); r.U = q.vec3(); r.V = q.vec3(); r.N = q.vec3();
    return r;
}
void putNode(ByteWriter& w, const LayoutNode& n) {
    w.vec3(n.p); w.vec3(n.n); w.vec3(n.c); w.vec3(n.v); w.f64(n.w); w.f64(n.h); w.f64(n.wL); w.f64(n.wR); w.f64(n.a); w.f64(n.b); w.f64(n.f);
}
LayoutNode getNode(ByteReader& q) {
    LayoutNode n;
    n.p = q.vec3(); n.n = q.vec3(); n.c = q.vec3(); n.v = q.vec3();
    n.w = q.f64(); n.h = q.f64(); n.wL = q.f64(); n.wR = q.f64(); n.a = q.f64(); n.b = q.f64(); n.f = q.f64();
    return n;
}
void putNodes(ByteWriter& w, const std::vector<LayoutNode>& v) { w.i64(int64_t(v.size())); for (const LayoutNode& n : v) putNode(w, n); }
std::vector<LayoutNode> getNodes(ByteReader& q) {
    std::vector<LayoutNode> v(q.count(8));
    for (LayoutNode& n : v) n = getNode(q);
    return v;
}

void putYarn(ByteWriter& w, const YarnPath& y) {
    w.i64(y.id); w.string64(y.name); w.i64(y.course); w.i64(y.reversed);
    w.i64(int64_t(y.curve.ctrl.size()));
    for (const Vec3& c : y.curve.ctrl) w.vec3(c);
    w.i64(int64_t(y.stitches.size()));
    for (const StitchSpan& sp : y.stitches) { w.i64(sp.wale); w.i64(sp.course); w.f64(sp.uBegin); w.f64(sp.uEnd); w.i64(sp.id); }
    const YarnPath::Tables& t = y.tables();
    w.doubles(t.tu); w.doubles(t.ts); w.doubles(t.preLo); w.doubles(t.preHi); w.f64(t.length); w.f64(t.lengthErr); w.i64(t.stitchesSorted);
    w.i64(static_cast<int64_t>(t.sortedCount));
}
YarnPath getYarn(ByteReader& q) {
    YarnPath y;
    y.id = int(q.i64()); y.name = q.string64(); y.course = int(q.i64()); y.reversed = q.i64() != 0;
    y.curve.ctrl.resize(q.count(24));
    for (Vec3& c : y.curve.ctrl) c = q.vec3();
    y.stitches.resize(q.count(8));
    for (StitchSpan& sp : y.stitches) { sp.wale = int(q.i64()); sp.course = int(q.i64()); sp.uBegin = q.f64(); sp.uEnd = q.f64(); sp.id = int(q.i64()); }
    YarnPath::Tables t;
    t.tu = q.doubles(); t.ts = q.doubles(); t.preLo = q.doubles(); t.preHi = q.doubles(); t.length = q.f64(); t.lengthErr = q.f64();
    t.stitchesSorted = q.i64() != 0;
    t.sortedCount = static_cast<size_t>(q.i64());
    y.setTables(std::move(t));
    return y;
}

void putProjState(ByteWriter& w, const ProjState& s) { w.i64(s.tri); w.f64(s.b[0]); w.f64(s.b[1]); w.f64(s.b[2]); w.vec3(s.foot); w.vec3(s.n); w.vec3(s.query); }
ProjState getProjState(ByteReader& q) {
    ProjState s;
    s.tri = int(q.i64()); s.b[0] = q.f64(); s.b[1] = q.f64(); s.b[2] = q.f64(); s.foot = q.vec3(); s.n = q.vec3(); s.query = q.vec3();
    return s;
}

} // namespace

void writeStitches(ByteWriter& w, const std::vector<StitchRecord>& v) { w.i64(int64_t(v.size())); for (const StitchRecord& r : v) putStitch(w, r); }
std::vector<StitchRecord> readStitches(ByteReader& q) {
    std::vector<StitchRecord> v(q.count(8));
    for (StitchRecord& r : v) r = getStitch(q);
    return v;
}

void writeCurves(ByteWriter& w, const CurveSet& cs) {
    w.i64(int64_t(cs.yarns.size()));
    for (const YarnPath& y : cs.yarns) putYarn(w, y);
    w.i64(0);  // a retired field: the file keeps its slot
    w.i64(static_cast<int>(cs.kind));
    writeStitches(w, cs.topo.stitches);
    w.i64(int64_t(cs.topo.terminals.size()));
    for (const TerminalRecord& t : cs.topo.terminals) {
        w.i64(static_cast<int>(t.op)); w.i64(t.yarn); w.f64(t.uBegin); w.f64(t.uEnd); w.ints(t.stitches);
        w.i64(int64_t(t.segments.size()));
        for (const auto& [a, b] : t.segments) { w.f64(a); w.f64(b); }
        w.vec3(t.centre); w.vec3(t.axis); w.f64(t.radius);
    }
    w.i64(cs.wales); w.i64(cs.courses); w.f64(cs.waleSpacing); w.f64(cs.courseSpacing);
}
void readCurves(ByteReader& q, CurveSet& cs) {
    cs.yarns.resize(q.count(8));
    for (YarnPath& y : cs.yarns) y = getYarn(q);
    q.i64();  // a retired field: the file keeps its slot
    cs.kind = static_cast<CurveKind>(q.i64());
    cs.topo.stitches = readStitches(q);
    cs.topo.terminals.resize(q.count(8));
    for (TerminalRecord& t : cs.topo.terminals) {
        t.op = static_cast<TerminalOp>(q.i64()); t.yarn = int(q.i64()); t.uBegin = q.f64(); t.uEnd = q.f64(); t.stitches = q.ints();
        t.segments.resize(q.count(16));
        for (auto& [a, b] : t.segments) { a = q.f64(); b = q.f64(); }
        t.centre = q.vec3(); t.axis = q.vec3(); t.radius = q.f64();
    }
    cs.wales = int(q.i64()); cs.courses = int(q.i64()); cs.waleSpacing = q.f64(); cs.courseSpacing = q.f64();
    // The topology's per-yarn index (derived): rebuilt as the knitting built it.
    std::vector<double> ends;
    for (const YarnPath& y : cs.yarns) ends.push_back(y.uEnd());
    if (!cs.topo.stitches.empty()) cs.topo.buildIndex(ends);
}

void writePatternParams(ByteWriter& w, const PatternParams& p) {
    w.i64(p.wales); w.i64(p.courses); w.f64(p.waleSpacing); w.f64(p.courseSpacing); w.f64(p.yarnDiameter); w.f64(p.clearance);
    w.f64(p.depthAmplitude); w.f64(p.headHalfWidth); w.f64(p.legHalfWidth); w.f64(p.sampleTolerance); w.i64(static_cast<int>(p.pattern));
    w.f64(p.transNeckY); w.f64(p.transNeckZ); w.f64(p.transLegY); w.f64(p.ribLayerOffset); w.f64(p.ribPairSpread); w.f64(p.garterWave);
}
PatternParams readPatternParams(ByteReader& q) {
    PatternParams p;
    p.wales = int(q.i64()); p.courses = int(q.i64()); p.waleSpacing = q.f64(); p.courseSpacing = q.f64(); p.yarnDiameter = q.f64();
    p.clearance = q.f64(); p.depthAmplitude = q.f64(); p.headHalfWidth = q.f64(); p.legHalfWidth = q.f64(); p.sampleTolerance = q.f64();
    p.pattern = static_cast<StitchPattern>(q.i64());
    p.transNeckY = q.f64(); p.transNeckZ = q.f64(); p.transLegY = q.f64(); p.ribLayerOffset = q.f64(); p.ribPairSpread = q.f64(); p.garterWave = q.f64();
    return p;
}

bool writeKnitCache(const std::string& path, const KnitCacheEntry& e, std::string* err) {
    ByteWriter w;
    w.bytes(kMagic, sizeof kMagic - 1);
    w.string64(e.key);
    const PieceKnit& k = e.knit;
    // Program.
    w.i64(int64_t(k.prog.turns.size()));
    for (const TurnPlan& t : k.prog.turns) { w.i64(t.input); w.i64(t.k); w.i64(t.inc); w.i64(t.dec); w.i64(t.output); }
    writeStitches(w, k.prog.stitches);
    w.ints(k.prog.turnStart);
    w.i64(k.prog.castRing); w.i64(k.prog.closeDrawstring);
    // Layout.
    putNodes(w, k.lay.nodes);
    putNode(w, k.lay.before);
    putNode(w, k.lay.after);
    w.i64(k.lay.helical); w.f64(k.lay.aSign);
    putNodes(w, k.lay.virtualNodes);
    w.ints(k.lay.virtualChildOf);
    w.ints(k.lay.rowPrev);
    w.ints(k.lay.rowNext);
    w.i64(int64_t(k.lay.waleHint.size()));
    for (const Vec3& v : k.lay.waleHint) w.vec3(v);
    w.i64(int64_t(k.lay.nodeState.size()));
    for (const ProjState& s : k.lay.nodeState) putProjState(w, s);
    // Curves and the rest of the piece knit.
    writeCurves(w, k.curves);
    w.i64(k.startAtA); w.f64(k.dir); w.ints(k.counts); w.f64(k.startArc); w.f64(k.endArc); w.string64(k.detail);
    writePatternParams(w, e.pp);
    w.i64(e.poleAInside); w.i64(e.poleBInside);
    w.i64(int64_t(e.attempts.size()));
    for (const KnitCacheAttempt& a : e.attempts) { w.f64(a.D); w.i64(static_cast<int>(a.verdict)); }
    w.i64(e.residue); w.i64(e.generatorChecks);
    w.i64(static_cast<int64_t>(e.meshDigest)); w.i64(static_cast<int64_t>(e.validationDigest));
    w.string64(e.log); w.f64(e.msKnit);
    w.bytes(kMagic, sizeof kMagic - 1);  // end marker (a truncated file is refused)
    return replaceFile(path, w.written(), err);
}

bool readKnitCache(const std::string& path, KnitCacheEntry& e, std::string* err) {
    std::string b;
    if (!readWholeFile(path, b, err)) return false;
    const size_t m = sizeof kMagic - 1;
    if (b.size() < 2 * m || b.compare(0, m, kMagic) != 0) { if (err) *err = path + ": not a knit cache (knit-cache v2)"; return false; }
    ByteReader q(b, m);
    e = KnitCacheEntry{};
    e.key = q.string64();
    PieceKnit& k = e.knit;
    k.prog.turns.resize(q.count(40));
    for (TurnPlan& t : k.prog.turns) { t.input = int(q.i64()); t.k = int(q.i64()); t.inc = int(q.i64()); t.dec = int(q.i64()); t.output = int(q.i64()); }
    k.prog.stitches = readStitches(q);
    k.prog.turnStart = q.ints();
    k.prog.castRing = q.i64() != 0; k.prog.closeDrawstring = q.i64() != 0;
    k.lay.nodes = getNodes(q);
    k.lay.before = getNode(q);
    k.lay.after = getNode(q);
    k.lay.helical = q.i64() != 0; k.lay.aSign = q.f64();
    k.lay.virtualNodes = getNodes(q);
    k.lay.virtualChildOf = q.ints();
    k.lay.rowPrev = q.ints();
    k.lay.rowNext = q.ints();
    k.lay.waleHint.resize(q.count(24));
    for (Vec3& v : k.lay.waleHint) v = q.vec3();
    k.lay.nodeState.resize(q.count(8));
    for (ProjState& s : k.lay.nodeState) s = getProjState(q);
    readCurves(q, k.curves);
    k.startAtA = q.i64() != 0; k.dir = q.f64(); k.counts = q.ints(); k.startArc = q.f64(); k.endArc = q.f64(); k.detail = q.string64();
    e.pp = readPatternParams(q);
    e.poleAInside = q.i64() != 0; e.poleBInside = q.i64() != 0;
    e.attempts.resize(q.count(16));
    for (KnitCacheAttempt& a : e.attempts) { a.D = q.f64(); a.verdict = static_cast<Status>(q.i64()); }
    e.residue = q.i64() != 0; e.generatorChecks = q.i64() != 0;
    e.meshDigest = static_cast<uint64_t>(q.i64()); e.validationDigest = static_cast<uint64_t>(q.i64());
    e.log = q.string64(); e.msKnit = q.f64();
    if (!q.ok() || q.at() + m != b.size() || b.compare(q.at(), m, kMagic) != 0) {
        if (err) *err = path + ": truncated or corrupt knit cache";
        return false;
    }
    return true;
}

std::string staleKeyField(const std::string& have, const std::string& want) {
    auto lines = [](const std::string& s) {
        std::vector<std::string> v;
        std::istringstream in(s);
        for (std::string l; std::getline(in, l);) v.push_back(l);
        return v;
    };
    const std::vector<std::string> h = lines(have), w = lines(want);
    for (size_t i = 0; i < std::max(h.size(), w.size()); ++i) {
        const std::string a = i < h.size() ? h[i] : "", b = i < w.size() ? w[i] : "";
        if (a != b) {
            const std::string& l = b.empty() ? a : b;
            const size_t eq = l.find(" = ");
            return eq == std::string::npos ? l : l.substr(0, eq);
        }
    }
    return "";
}

bool sameKnitDiameter(double a, double b) { return std::abs(a - b) <= kKnitDiameterToleranceMm; }

bool selectCachedAttempt(const KnitCacheEntry& e, const std::vector<double>& dListIn, double requestedD, double* missing) {
    if (e.attempts.empty()) { if (missing) *missing = requestedD; return false; }
    // The same walk as the pipeline's loop: the first D whose gate (a) passes is kept; if every
    // D fails, the requested D is appended once and kept (residue).
    std::vector<double> dList = dListIn;
    const KnitCacheAttempt& kept = e.attempts.back();
    for (size_t di = 0; di < dList.size(); ++di) {
        const double D = dList[di];
        const KnitCacheAttempt* a = nullptr;
        for (const KnitCacheAttempt& x : e.attempts)
            if (sameKnitDiameter(x.D, D)) { a = &x; break; }
        if (!a) { if (missing) *missing = D; return false; }
        if (a->verdict == Status::Pass) {
            // Only the kept attempt's result is stored (the original loop stopped at it).
            if (!sameKnitDiameter(kept.D, D) || kept.verdict != Status::Pass) { if (missing) *missing = D; return false; }
            return true;
        }
        if (di + 1 == dList.size() && dList.size() > 1 && dList.back() != requestedD) dList.push_back(requestedD);
    }
    // Every D failed: the loop keeps its last attempt (the residue at the requested D when
    // fallbacks were tried).
    const double last = dList.back();
    if (!sameKnitDiameter(kept.D, last) || kept.verdict == Status::Pass) { if (missing) *missing = last; return false; }
    return true;
}

} // namespace knit
