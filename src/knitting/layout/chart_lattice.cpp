// The stitch lattice on the UV charts of a welded mesh: nodes, yarns, program, layout and curves (chart_lattice.h).
#include "knitting/layout/chart_lattice.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <tuple>

#include "model/cancel.h"
#include "model/mesh_query.h"
#include "model/spatial_hash.h"
#include "knitting/stitch/stockinette.h"

namespace knit {

namespace {

double dot2(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }

// The wale direction in the chart's uv': its up direction turned by the orientation on the mesh.
Vec2 waleDirection(const UvChart& ch, double angleDeg) {
    const double th = angleDeg * kPi / 180;
    return {ch.upDir.x * std::cos(th) - ch.upDir.y * std::sin(th), ch.upDir.x * std::sin(th) + ch.upDir.y * std::cos(th)};
}

struct UvRect {
    double x0, x1, y0, y1;
};

// Does segment (a, b) meet the closed axis-aligned rectangle? (Liang-Barsky clipping.)
bool segmentMeetsRect(Vec2 a, Vec2 b, const UvRect& rect) {
    double t0 = 0, t1 = 1;
    const double dx = b.x - a.x, dy = b.y - a.y;
    auto clip = [&](double p, double q) {
        if (p == 0) return q >= 0;
        const double r = q / p;
        if (p < 0) {
            if (r > t1) return false;
            t0 = std::max(t0, r);
        } else {
            if (r < t0) return false;
            t1 = std::min(t1, r);
        }
        return true;
    };
    return clip(-dx, a.x - rect.x0) && clip(dx, rect.x1 - a.x) && clip(-dy, a.y - rect.y0) && clip(dy, rect.y1 - a.y) && t0 <= t1;
}

struct Node {
    int i = 0, r = 0, tri = -1;
    double bary[3] = {0, 0, 0};
    Vec3 p;
    double x = 0, y = 0;  // lattice-frame uv' coordinates
    Vec3 c3, v3;          // overhang: the node's course / wale directions (unit, tangent)
    bool walked = false;  // reached by a surface walk beyond the chart (its wale hint is v3)
    // The run-end footprint tests, if already known (classified on the device, overhang): 1 clear, 0 not, -1 test here.
    signed char clearStart = -1, clearEnd = -1;
};

// A step on the lattice: the 3D displacement to project and the change of (i, r).
struct LatticeStep {
    Vec3 move;
    int di, dr;
};

// One walk step of the overlap mode from node a on the surface (continueFrom); the frame is carried
// into the new tangent plane.
bool walkStep(const MeshSurface& surf, const ContinuationParams& cp, const Node& a, const LatticeStep& step, Node& nx) {
    const ProjState s0 = surf.stateAt(a.tri, a.bary, a.p);
    ProjState s1;
    if (surf.continueFrom(s0, a.p + step.move, s1, cp) != ProjResult::Owned) return false;
    nx = a;
    nx.i = a.i + step.di;
    nx.r = a.r + step.dr;
    nx.tri = s1.tri;
    for (int k = 0; k < 3; ++k) nx.bary[k] = s1.b[k];
    nx.p = s1.foot;
    auto tangent = [&](const Vec3& v) { return normalized(v - s1.n * dot(v, s1.n)); };
    nx.c3 = tangent(a.c3);
    nx.v3 = tangent(a.v3);
    nx.walked = true;
    return true;
}

double pointSegmentDistance(const Vec3& p, const Vec3& a, const Vec3& b) {
    const Vec3 ab = b - a;
    const double len2 = dot(ab, ab);
    const double t = len2 > 0 ? std::clamp(dot(p - a, ab) / len2, 0.0, 1.0) : 0.0;
    return norm(p - (a + ab * t));
}

// Without the footprint test (the overlap mode: a chart's stitches cross its seams onto the
// neighbouring charts) a stitch must still not reach over an OPEN edge of the mesh (an edge of
// one triangle): there is no surface beyond it, and yarn placed there collapses onto the edge
// (near-coincident curves, whose local signature check then subdivides for minutes).
class OpenEdgeTest {
public:
    OpenEdgeTest(const TriMesh& m, double reach) : hash_(std::max(reach, 1e-6)), reach_(reach) {
        std::map<std::pair<int, int>, int> uses;
        for (const auto& t : m.tri)
            for (int k = 0; k < 3; ++k) ++uses[{std::min(t[k], t[(k + 1) % 3]), std::max(t[k], t[(k + 1) % 3])}];
        for (const auto& [e, n] : uses) {
            if (n != 1) continue;
            const Vec3 a = m.pos[e.first], b = m.pos[e.second];
            edges_.push_back({a, b});
            hash_.insert({std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}, {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)});
        }
    }
    // The point is closer than the reach to an open edge.
    bool near(const Vec3& p) const {
        bool found = false;
        const Vec3 r{reach_, reach_, reach_};
        hash_.query(p - r, p + r, [&](int e) { found = found || pointSegmentDistance(p, edges_[e][0], edges_[e][1]) < reach_; });
        return found;
    }

private:
    std::vector<std::array<Vec3, 2>> edges_;
    SpatialHash3 hash_;
    double reach_;
};

} // namespace

bool LatticeSolids::contains(const Vec3& p) const {
    const Vec3 x = p + offset;
    for (const MeshQuery* q : meshes)
        if (q->inside(x)) return true;
    return false;
}

LatticeFrame latticeFrame(const UvAtlas& A, int chart, const LatticeSpacing& spacing) {
    LatticeFrame f;
    const UvChart& ch = A.charts()[chart];
    const double w = spacing.w, h = spacing.h;
    f.v2 = waleDirection(ch, spacing.angleDeg);
    f.c2 = {f.v2.y, -f.v2.x};  // course: wale rotated by -90 deg (c x v . n > 0 in uv')
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
    for (int t : ch.tris) {
        if (!A.usable(t)) continue;
        for (const Vec2& q : A.uvp(t)) {
            xmin = std::min(xmin, dot2(q, f.c2)), xmax = std::max(xmax, dot2(q, f.c2));
            ymin = std::min(ymin, dot2(q, f.v2)), ymax = std::max(ymax, dot2(q, f.v2));
        }
    }
    if (xmin > xmax) return f;
    f.ok = true;
    f.x0 = xmin;
    f.y0 = ymin;
    f.ni = std::max(0, static_cast<int>(std::floor((xmax - f.x0) / w)) + 1);
    f.nr = std::max(0, static_cast<int>(std::floor((ymax - f.y0) / h)) + 1);
    return f;
}

LatticeCheck checkLatticeProgram(const StitchProgram& prog) {
    LatticeCheck c;
    c.stitches = static_cast<int>(prog.stitches.size());
    auto fail = [&](const std::string& why) {
        if (c.status != Status::Fail) c.detail = why;
        c.status = Status::Fail;
    };
    for (size_t k = 0; k < prog.stitches.size(); ++k) {
        const StitchRecord& s = prog.stitches[k];
        if (s.id != static_cast<int>(k)) fail("ids not in yarn order at " + std::to_string(k));
        if (s.lattice.chart < 0) fail("stitch " + std::to_string(k) + " without lattice coordinates");
        if (k > 0) {
            const StitchRecord& q = prog.stitches[k - 1];
            if (s.turn < q.turn) fail("runs out of order at stitch " + std::to_string(k));
            if (s.turn == q.turn && (s.lattice.chart != q.lattice.chart || s.lattice.r != q.lattice.r || s.lattice.i != q.lattice.i + 1))
                fail("run " + std::to_string(s.turn) + " not contiguous in i at stitch " + std::to_string(k));
        }
        if (s.nParents > 1) fail("stitch " + std::to_string(k) + " with " + std::to_string(s.nParents) + " parents");
        if (s.nParents == 1) {
            const StitchRecord& p = prog.stitches[s.parents[0]];
            if (p.lattice.chart != s.lattice.chart || p.lattice.i != s.lattice.i || p.lattice.r != s.lattice.r - 1)
                fail("stitch " + std::to_string(k) + ": parent not at (i, r - 1) of the same chart");
        } else {
            ++c.parentless;
        }
        if (s.nChildren > 1) fail("loop " + std::to_string(k) + " with " + std::to_string(s.nChildren) + " children");
        if (s.nChildren == 0) ++c.uncaught;
    }
    if (c.status == Status::Pass) {
        char b[160];
        std::snprintf(b, sizeof b, "%d stitches: lattice graph consistent; %d parentless (cast-on omitted), %d uncaught loops (bind-off omitted)",
                      c.stitches, c.parentless, c.uncaught);
        c.detail = b;
    }
    return c;
}

TemplateExtent templateExtent(const PatternParams& pp) {
    TemplateExtent e{1e300, -1e300, 0};
    for (int f = 0; f < 8; ++f) {
        StitchFaces sf;
        sf.own = f & 1 ? StitchFace::Purl : StitchFace::Knit;
        sf.child = f & 2 ? StitchFace::Purl : StitchFace::Knit;
        sf.left = f & 4 ? StitchFace::Purl : StitchFace::Knit;
        const std::vector<TemplateControl> t = stitchTemplate(pp, sf);
        for (const TemplateControl& c : t) e.Ymin = std::min(e.Ymin, c.Y), e.Ymax = std::max(e.Ymax, c.Y);
        e.Xpad = std::max({e.Xpad, 1 - t[15].X, 1 + t[1].X});  // padding controls beyond the run's end nodes
    }
    return e;
}

namespace {

using LatticeSite = std::pair<int, int>;  // (r, i)

// A run: consecutive stitches of one row, knitted as one yarn.
struct LatticeRun {
    int chart, r;
    std::vector<int> nodes;  // indices into the collected nodes
};

// Collects the surface nodes of the lattice and groups them into runs, chart by chart or on one
// walked patch. A node is kept if its stitch fits on the surface (chart_lattice.h).
class LatticeRunBuilder {
public:
    LatticeRunBuilder(const ChartedMesh& on, const TileParams& p, KnitNumerics& numerics);

    const PatternParams& pattern() const { return pp_; }
    const std::vector<Node>& nodes() const { return nodes_; }
    const std::vector<LatticeRun>& runs() const { return runs_; }

    // TileParams::patch: the walked local lattice (no chart; its stitches are filed under chart 0).
    bool walkPatch(std::string* err);
    // One chart's nodes, overhang and runs, counted in ct.
    void tileChart(int chart, ChartTiling& ct);

private:
    using NodeMap = std::map<LatticeSite, int>;  // (r, i) -> nodes_ index

    bool nearOpenEdge(const Vec3& q) const { return openEdges_ && openEdges_->near(q); }
    LatticeChartQuery chartQuery(int chart, const LatticeFrame& f) const;
    bool nodeClear(const LatticeChartQuery& q, const NodeExtent& node) const;
    bool takeClassifiedNode(const LatticeNodeClass& c, const LatticeChartQuery& q, Node& nd, ChartTiling& ct) const;
    void placeChartNodes(const LatticeChartQuery& q, NodeMap& kept, ChartTiling& ct);
    void walkOverhang(const LatticeChartQuery& q, NodeMap& kept, ChartTiling& ct);
    void collectChartRuns(const LatticeChartQuery& q, const NodeMap& kept, ChartTiling& ct);
    void collectPatchRuns(const NodeMap& at);

    const UvAtlas& A_;
    const MeshSurface& surf_;
    const TileParams& p_;
    KnitNumerics& numerics_;
    PatternParams pp_;                     // p.pp at the gauge
    double w_, h_;                         // gauge [mm]
    double margin_;                        // support radius + half the clearance
    TemplateExtent extent_{};              // over every face combination (transitions move heights)
    ContinuationParams walkParams_;        // of the surface walks
    std::unique_ptr<OpenEdgeTest> openEdges_;  // without the footprint test: the open edges still bound the stitches
    std::vector<Node> nodes_;
    std::vector<LatticeRun> runs_;
};

ContinuationParams surfaceWalkParams(const TileParams& p) {
    ContinuationParams cp = p.cp;
    if (cp.lambdaW <= 0) cp.lambdaW = 0.5 * p.g.w;
    return cp;
}

LatticeRunBuilder::LatticeRunBuilder(const ChartedMesh& on, const TileParams& p, KnitNumerics& numerics)
    : A_(on.atlas), surf_(on.surface), p_(p), numerics_(numerics), pp_(p.pp), w_(p.g.w), h_(p.g.h), walkParams_(surfaceWalkParams(p)) {
    pp_.waleSpacing = w_;
    pp_.courseSpacing = h_;
    const double Rs = p.supportRadius > 0 ? p.supportRadius : pp_.radius();
    margin_ = Rs + 0.5 * pp_.clearance;
    extent_ = templateExtent(pp_);
    // In a chart the footprint test against its open edges (footprintClear); beyond it (overhang
    // walks, gap patches) this 3D test: the footprint's largest half extent plus the margin, isotropic.
    const double openReach = std::max({0.5 * w_, -extent_.Ymin * h_, extent_.Ymax * h_}) + margin_;
    if (!p.footprintCheck && (p.overhang > 0 || p.patch.tri >= 0)) openEdges_ = std::make_unique<OpenEdgeTest>(A_.mesh(), openReach);
}

// The patch's seed node: the walk starts here with wale `wale` projected into the tangent plane.
Node patchSeed(const TriMesh& m, const MeshSurface& surf, const TileParams::WalkPatch& patch) {
    Node seed;
    seed.tri = patch.tri;
    for (int k = 0; k < 3; ++k) seed.bary[k] = patch.bary[k];
    const std::array<int, 3>& tv = m.tri[seed.tri];
    seed.p = m.pos[tv[0]] * seed.bary[0] + m.pos[tv[1]] * seed.bary[1] + m.pos[tv[2]] * seed.bary[2];
    const Vec3 n0 = surf.normalAt(seed.tri, seed.bary);
    seed.v3 = patch.wale - n0 * dot(patch.wale, n0);
    // A seed facing exactly along the wale has no projected wale, and any tangent direction tiles it.
    if (norm(seed.v3) < 1e-6) {
        const Vec3 axis = std::fabs(n0.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
        seed.v3 = axis - n0 * dot(axis, n0);
    }
    seed.v3 = normalized(seed.v3);
    seed.c3 = normalized(cross(seed.v3, n0));
    seed.walked = true;
    return seed;
}

void LatticeRunBuilder::collectPatchRuns(const NodeMap& at) {
    std::map<int, std::vector<int>> rowsOf;
    for (const auto& kv : at) rowsOf[kv.first.first].push_back(kv.first.second);
    for (auto& [r, is] : rowsOf) {
        std::sort(is.begin(), is.end());
        std::vector<int> cur;
        auto flush = [&] {
            if (static_cast<int>(cur.size()) >= p_.minRun) runs_.push_back({0, r, cur});
            cur.clear();
        };
        for (size_t k = 0; k < is.size(); ++k) {
            if (k > 0 && is[k] != is[k - 1] + 1) flush();
            cur.push_back(at.at({r, is[k]}));
        }
        flush();
    }
}

bool LatticeRunBuilder::walkPatch(std::string* err) {
    if (A_.charts().empty()) {
        if (err) *err = "walk patch: the atlas has no chart";
        return false;
    }
    const TileParams::WalkPatch& patch = p_.patch;
    const Node seed = patchSeed(A_.mesh(), surf_, patch);
    NodeMap at;
    auto inside = [&](const Node& nd) {
        return norm(nd.p - seed.p) <= patch.radius && !p_.solids.contains(nd.p) && !nearOpenEdge(nd.p);
    };
    auto add = [&](const Node& nd) {
        at[{nd.r, nd.i}] = static_cast<int>(nodes_.size());
        nodes_.push_back(nd);
    };
    // The seed's wale, then a course through each of its nodes.
    const int nr = static_cast<int>(std::ceil(patch.radius / h_)) + 1, ni = static_cast<int>(std::ceil(patch.radius / w_)) + 1;
    std::vector<Node> column{seed};
    for (const int dir : {-1, 1}) {
        Node cur = seed;
        for (int k = 0; k < nr; ++k) {
            Node nx;
            if (!walkStep(surf_, walkParams_, cur, {cur.v3 * dir * h_, 0, dir}, nx) || norm(nx.p - seed.p) > patch.radius + h_) break;
            column.push_back(nx);
            cur = nx;
        }
    }
    for (const Node& base : column) {
        if (inside(base)) add(base);
        for (const int dir : {-1, 1}) {
            Node cur = base;
            for (int k = 0; k < ni; ++k) {
                Node nx;
                if (!walkStep(surf_, walkParams_, cur, {cur.c3 * dir * w_, dir, 0}, nx) || norm(nx.p - seed.p) > patch.radius + w_) break;
                if (inside(nx)) add(nx);
                cur = nx;
            }
        }
    }
    collectPatchRuns(at);
    return true;
}

LatticeChartQuery LatticeRunBuilder::chartQuery(int chart, const LatticeFrame& f) const {
    LatticeChartQuery q;
    q.chart = chart;
    q.c2 = f.c2;
    q.v2 = f.v2;
    q.xmin = f.x0;
    q.ymin = f.y0;
    q.w = w_;
    q.h = h_;
    q.ni = f.ni;
    q.nr = f.nr;
    q.margin = margin_;
    q.Ymin = extent_.Ymin;
    q.Ymax = extent_.Ymax;
    q.Xpad = extent_.Xpad;
    q.footprintCheck = p_.footprintCheck;
    return q;
}

// The footprint test, where one applies: not on a part that runs past its chart with no open edge to respect.
bool LatticeRunBuilder::nodeClear(const LatticeChartQuery& q, const NodeExtent& node) const {
    return (!p_.footprintCheck && !openEdges_) || footprintClear(A_, q, node);
}

// The node as the batch classified it; false (counted in ct) if it is not kept.
bool LatticeRunBuilder::takeClassifiedNode(const LatticeNodeClass& c, const LatticeChartQuery& q, Node& nd, ChartTiling& ct) const {
    if (c.code == 1) { ++ct.outside; return false; }
    if (c.code == 2) { ++ct.ambiguous; return false; }
    if (c.code == 3) { ++ct.insideOther; return false; }
    if (c.code == 4) { ++ct.boundary; return false; }
    nd.tri = c.tri;
    for (int k = 0; k < 3; ++k) nd.bary[k] = c.bary[k];
    nd.p = c.p;
    // The overlap mode: the open edges here (the classification ran without a footprint test).
    if (openEdges_ && !nodeClear(q, {nd.x, nd.y, 0.5 * w_, 0.5 * w_})) { ++ct.boundary; return false; }
    nd.clearStart = openEdges_ ? -1 : c.clearStart;
    nd.clearEnd = openEdges_ ? -1 : c.clearEnd;
    return true;
}

void LatticeRunBuilder::placeChartNodes(const LatticeChartQuery& q, NodeMap& kept, ChartTiling& ct) {
    std::vector<LatticeNodeClass> cls;
    numerics_.classifyLattice(A_, q, p_.solids, cls);
    if (cls.size() != static_cast<size_t>(q.ni) * q.nr) throw PartFailure("lattice node classification: size mismatch");
    for (int r = 0; r < q.nr && (cancelCheckpoint(), true); ++r)  // per lattice row
        for (int i = 0; i < q.ni; ++i) {
            Node nd;
            nd.i = i;
            nd.r = r;
            nd.x = q.xmin + i * w_;
            nd.y = q.ymin + r * h_;
            ++ct.candidates;
            if (!takeClassifiedNode(cls[static_cast<size_t>(r) * q.ni + i], q, nd, ct)) continue;
            kept[{r, i}] = static_cast<int>(nodes_.size());
            nodes_.push_back(nd);
        }
}

// TileParams::overhang: walks beyond the chart, first along the columns, then along the runs.
void LatticeRunBuilder::walkOverhang(const LatticeChartQuery& q, NodeMap& kept, ChartTiling& ct) {
    for (auto& kv : kept) {
        Node& nd = nodes_[kv.second];
        nd.c3 = normalized(A_.pushForward(nd.tri, q.c2));
        nd.v3 = normalized(A_.pushForward(nd.tri, q.v2));
    }
    auto walk = [&](const Node& a, const LatticeStep& step, Node& nx) {
        if (!walkStep(surf_, walkParams_, a, step, nx)) return false;
        nx.x = q.xmin + nx.i * w_;
        nx.y = q.ymin + nx.r * h_;
        return !p_.solids.contains(nx.p) && !nearOpenEdge(nx.p);
    };
    auto add = [&](Node nd) {
        nd.clearStart = nd.clearEnd = 1;  // beyond the chart: its border does not apply
        kept[{nd.r, nd.i}] = static_cast<int>(nodes_.size());
        nodes_.push_back(nd);
        ++ct.overhang;
    };
    std::vector<LatticeSite> base;
    for (const auto& kv : kept) base.push_back(kv.first);
    for (const auto& [r0, i0] : base)
        for (const int dir : {-1, 1}) {
            Node cur = nodes_[kept.at({r0, i0})];
            for (int k = 0; k < p_.overhang && !kept.count({cur.r + dir, i0}); ++k) {
                Node nx;
                if (!walk(cur, {cur.v3 * dir * h_, 0, dir}, nx)) break;
                add(nx);
                cur = nx;
            }
        }
    base.clear();
    for (const auto& kv : kept) base.push_back(kv.first);
    for (const auto& [r0, i0] : base)
        for (const int dir : {-1, 1}) {
            if (kept.count({r0, i0 + dir})) continue;  // not a run end on this side
            Node cur = nodes_[kept.at({r0, i0})];
            for (int k = 0; k < p_.overhang && !kept.count({r0, cur.i + dir}); ++k) {
                Node nx;
                if (!walk(cur, {cur.c3 * dir * w_, dir, 0}, nx)) break;
                add(nx);
                cur = nx;
            }
        }
}

// Runs: consecutive i per row; the end nodes also need room for the padding controls.
void LatticeRunBuilder::collectChartRuns(const LatticeChartQuery& q, const NodeMap& kept, ChartTiling& ct) {
    auto clearAtStart = [&](int n) {
        const Node& nd = nodes_[n];
        return nd.clearStart >= 0 ? nd.clearStart != 0 : nodeClear(q, {nd.x, nd.y, extent_.Xpad * w_, 0.5 * w_});
    };
    auto clearAtEnd = [&](int n) {
        const Node& nd = nodes_[n];
        return nd.clearEnd >= 0 ? nd.clearEnd != 0 : nodeClear(q, {nd.x, nd.y, 0.5 * w_, extent_.Xpad * w_});
    };
    int rLo = 0, rHi = q.nr - 1, iLo = 0, iHi = q.ni - 1;
    for (const auto& kv : kept) {
        rLo = std::min(rLo, kv.first.first);
        rHi = std::max(rHi, kv.first.first);
        iLo = std::min(iLo, kv.first.second);
        iHi = std::max(iHi, kv.first.second);
    }
    for (int r = rLo; r <= rHi; ++r) {
        cancelCheckpoint();  // per row of runs
        std::vector<int> cur;
        auto flush = [&] {
            while (!cur.empty() && !clearAtStart(cur.front())) {
                cur.erase(cur.begin());
                ++ct.runEnd;
            }
            while (!cur.empty() && !clearAtEnd(cur.back())) {
                cur.pop_back();
                ++ct.runEnd;
            }
            if (static_cast<int>(cur.size()) >= p_.minRun) {
                runs_.push_back({q.chart, r, cur});
                ++ct.runs;
            } else if (!cur.empty()) {
                ct.shortRun += static_cast<int>(cur.size());
                ++ct.droppedRuns;
            }
            cur.clear();
        };
        for (int i = iLo; i <= iHi; ++i) {
            auto it = kept.find({r, i});
            if (it == kept.end()) flush();
            else cur.push_back(it->second);
        }
        flush();
    }
}

void LatticeRunBuilder::tileChart(int chart, ChartTiling& ct) {
    cancelCheckpoint();  // per chart (lattice, node classification, runs)
    const UvChart& ch = A_.charts()[chart];
    ct.chart = chart;
    ct.status = ch.status;
    ct.upFromWorld = ch.upFromWorld;
    ct.baseAngleDeg = std::atan2(ch.upDir.y, ch.upDir.x) * 180 / kPi;
    if (ch.status != Status::Pass) return;
    if (p_.onlyChart >= 0 && chart != p_.onlyChart) {
        ct.skipped = true;
        return;
    }
    const LatticeFrame f = latticeFrame(A_, chart, {p_.angleDeg, w_, h_});
    if (!f.ok) return;
    const LatticeChartQuery q = chartQuery(chart, f);
    NodeMap kept;
    placeChartNodes(q, kept, ct);
    if (p_.overhang > 0) walkOverhang(q, kept, ct);
    collectChartRuns(q, kept, ct);
}

std::string noStitchFitsMessage(const UvAtlas& A, const std::vector<ChartTiling>& charts) {
    std::string message = "no stitch fits the charts at this gauge:";
    for (const ChartTiling& ct : charts) {
        const UvChart& ch = A.charts()[ct.chart];
        char b[320];
        std::snprintf(b, sizeof b,
                      " [chart %d %s, %.0f mm2 (folded %.0f, overlap %.0f, degenerate %.0f): nodes %d, outside %d, ambiguous %d, "
                      "boundary %d, inside other %d, run ends %d, short runs %d]",
                      ct.chart, toString(ct.status), ch.area, ch.foldedArea, ch.overlapArea, ch.degenerateArea, ct.candidates,
                      ct.outside, ct.ambiguous, ct.boundary, ct.insideOther, ct.runEnd, ct.shortRun);
        message += b;
    }
    return message;
}

// The stitch program: stitches in (chart, r, i) order, one turn per run; a stitch's parent is the
// one at (i, r - 1) of its chart. Returns the node of each stitch id.
std::vector<int> buildStitchProgram(const std::vector<LatticeRun>& runs, const std::vector<Node>& nodes, StitchPattern pattern, StitchProgram& prog) {
    std::map<std::tuple<int, int, int>, int> idOf;  // (chart, i, r) -> stitch id
    std::vector<int> nodeOf;
    for (size_t ru = 0; ru < runs.size(); ++ru) {
        if ((ru & 0x3ff) == 0) cancelCheckpoint();
        const LatticeRun& run = runs[ru];
        prog.turnStart.push_back(static_cast<int>(prog.stitches.size()));
        TurnPlan tp;
        tp.k = tp.output = static_cast<int>(run.nodes.size());
        prog.turns.push_back(tp);
        for (int n : run.nodes) {
            const Node& nd = nodes[n];
            StitchRecord s;
            s.id = static_cast<int>(prog.stitches.size());
            s.turn = static_cast<int>(ru);
            s.index = nd.i;
            s.lattice = {run.chart, nd.i, nd.r, static_cast<int>(ru)};
            s.face = stitchFaceIsPurl(pattern, nd.i, nd.r) ? StitchFace::Purl : StitchFace::Knit;
            idOf[{run.chart, nd.i, nd.r}] = s.id;
            nodeOf.push_back(n);
            prog.stitches.push_back(s);
        }
    }
    for (StitchRecord& s : prog.stitches) {
        auto it = idOf.find({s.lattice.chart, s.lattice.i, s.lattice.r - 1});
        if (it != idOf.end()) {
            s.parents[0] = it->second;
            s.nParents = 1;
        }
    }
    linkChildren(prog.stitches);
    return nodeOf;
}

// The surface nodes grouped into runs: on the walked patch, or chart by chart. false: nothing to knit.
bool buildRuns(LatticeRunBuilder& lattice, const UvAtlas& A, const TileParams& p, std::vector<ChartTiling>& charts, std::string* err) {
    if (p.patch.tri >= 0) {
        if (!lattice.walkPatch(err)) return false;
    } else {
        try {
            for (size_t ci = 0; ci < A.charts().size(); ++ci) {
                ChartTiling ct;
                lattice.tileChart(static_cast<int>(ci), ct);
                charts.push_back(ct);
            }
        } catch (const PartFailure& failure) {
            if (err) *err = failure.what();
            return false;
        }
    }
    if (lattice.runs().empty()) {
        if (err) *err = noStitchFitsMessage(A, charts);
        return false;
    }
    return true;
}

// The layout's nodes on the surface and each stitch's three faces, then the frames.
void layOutNodes(const ChartedMesh& on, const TileParams& p, const std::vector<Node>& nodes, const std::vector<int>& nodeOf, LatticeKnit& out) {
    const UvAtlas& A = on.atlas;
    const MeshSurface& surf = on.surface;
    const StitchProgram& prog = out.prog;
    PieceLayout& lay = out.lay;
    const int N = static_cast<int>(prog.stitches.size());
    lay.nodes.resize(N);
    lay.rowPrev.assign(N, -1);
    lay.rowNext.assign(N, -1);
    lay.waleHint.resize(N);
    lay.nodeState.resize(N);
    out.faces.resize(N);
    for (int k = 0; k < N; ++k) {
        const StitchRecord& s = prog.stitches[k];
        const Node& nd = nodes[nodeOf[k]];
        if (k > 0 && prog.stitches[k - 1].turn == s.turn) lay.rowPrev[k] = k - 1;
        if (k + 1 < N && prog.stitches[k + 1].turn == s.turn) lay.rowNext[k] = k + 1;
        LayoutNode& ln = lay.nodes[k];
        ln.p = nd.p;
        ln.n = surf.normalAt(nd.tri, nd.bary);
        ln.a = nd.r;
        ln.b = nd.i;
        lay.nodeState[k] = surf.stateAt(nd.tri, nd.bary, nd.p);
        lay.waleHint[k] = nd.walked ? nd.v3 : A.pushForward(nd.tri, waleDirection(A.charts()[s.lattice.chart], p.angleDeg));
        StitchFaces& f = out.faces[k];
        f.own = s.face;
        f.child = s.nChildren ? prog.stitches[s.children[0]].face : s.face;
        f.left = lay.rowPrev[k] >= 0 ? prog.stitches[k - 1].face : s.face;
    }
    computeFrames(surf, prog, p.g, lay);
}

// The curves on the surface, with the projection binding; numerics stores them.
void placeYarn(const MeshSurface& surf, const TileParams& p, KnitNumerics& numerics, LatticeKnit& out) {
    out.bind.surface = &surf;
    out.bind.cp = surfaceWalkParams(p);
    out.curves = placeRows(surf, out.lay, out.prog, out.pp, {p.sh, CurveKind::Surface, &out.faces, &out.bind});
    numerics.store(out.curves, controlsOf(out.curves), {});
    out.curves.waleSpacing = p.g.w;
    out.curves.courseSpacing = p.g.h;
}

// Kept stitches per chart and the share of the chart they cover.
void recordChartCoverage(const UvAtlas& A, const GaugeParams& g, LatticeKnit& out) {
    std::vector<int> perChart(A.charts().size(), 0);
    for (const StitchRecord& s : out.prog.stitches) ++perChart[s.lattice.chart];
    for (ChartTiling& ct : out.charts) {
        ct.kept = perChart[ct.chart];
        const double area = A.charts()[ct.chart].area;
        ct.coverage = area > 0 ? ct.kept * g.w * g.h / area : 0;
    }
}

} // namespace

bool footprintClear(const UvAtlas& A, const LatticeChartQuery& q, const NodeExtent& node) {
    const std::vector<std::array<Vec2, 2>>& border = q.footprintCheck ? A.boundary(q.chart) : A.openBoundary(q.chart);
    const double x = node.x, y = node.y;
    double mU = q.margin;
    for (int it = 0; it < 2; ++it) {
        const double x0 = x - node.left - mU, x1 = x + node.right + mU, y0 = y + q.Ymin * q.h - mU, y1 = y + q.Ymax * q.h + mU;
        Vec2 lo{1e300, 1e300}, hi{-1e300, -1e300};
        for (double cx : {x0, x1})
            for (double cy : {y0, y1}) {
                const Vec2 c = q.c2 * cx + q.v2 * cy;
                lo = {std::min(lo.x, c.x), std::min(lo.y, c.y)};
                hi = {std::max(hi.x, c.x), std::max(hi.y, c.y)};
            }
        const double s2 = A.minS2In(q.chart, lo, hi);
        const double next = q.margin / std::max(1e-3, std::min(1.0, s2));
        if (!(next > mU)) break;
        mU = next;
    }
    const UvRect box{x - node.left - mU, x + node.right + mU, y + q.Ymin * q.h - mU, y + q.Ymax * q.h + mU};
    for (const auto& e : border) {
        const Vec2 a{dot2(e[0], q.c2), dot2(e[0], q.v2)}, b{dot2(e[1], q.c2), dot2(e[1], q.v2)};
        if (segmentMeetsRect(a, b, box)) return false;
    }
    return true;
}

bool tileMesh(const ChartedMesh& on, const TileParams& p, KnitNumerics& numerics, LatticeKnit& out, std::string* err) {
    const UvAtlas& A = on.atlas;
    out = LatticeKnit{};
    LatticeRunBuilder lattice(on, p, numerics);
    out.pp = lattice.pattern();
    if (!buildRuns(lattice, A, p, out.charts, err)) return false;
    const std::vector<int> nodeOf = buildStitchProgram(lattice.runs(), lattice.nodes(), out.pp.pattern, out.prog);
    layOutNodes(on, p, lattice.nodes(), nodeOf, out);
    placeYarn(on.surface, p, numerics, out);
    out.lattice = checkLatticeProgram(out.prog);
    out.layout = checkLayout(out.prog, out.lay, p.g);
    recordChartCoverage(A, p.g, out);
    char b[200];
    std::snprintf(b, sizeof b, "%zu stitches in %zu runs on %zu charts", out.prog.stitches.size(), lattice.runs().size(), A.charts().size());
    out.detail = b;
    return true;
}

} // namespace knit
