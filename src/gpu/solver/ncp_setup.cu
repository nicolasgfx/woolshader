// ncp-gpu: one relax call on the device: samples, rows, upload, the capacities, download, and the single
// operations the tests run (ncp.h). Everything here is host preparation or a test: the iteration
// itself is the graph of ncp_control.cu.
#include "gpu/solver/ncp_device.cuh"

#include <cstring>

#include "knitting/layout/stitch_layout.h"
#include "knitting/knit/piece_knit.h"
#include "knitting/layout/uv_atlas.h"

namespace knit::gpu {

namespace {

long long gIndexLimitForTests = 0x7fffffffLL;
double gHeadroomForTests = -1;
long long gHeadroomMinForTests = -1;
double gLiveHeadroomForTests = -1;
long long gLiveHeadroomMinForTests = -1;

// Headroom of the contact capacity over the starting state's rows and candidates. Contacts enter and
// leave while the yarn moves; the measured growth over a teddy settle call is reported per call
// (NcpStats::maxContactRows against initialContactRows).
constexpr double kContactHeadroom = 1.5;
constexpr long long kContactHeadroomMin = 65536;
// The least capacity a memory cap may leave: the rows plus 5 % (growth beyond it is met by regrowing,
// solveNcp). Below it the solve is refused by the cap instead of starting at the edge.
constexpr double kCapacityFloor = 1.05;
// The grid's fixed range: the starting samples' box grown by this much on every side [mm] (items that
// move beyond it are clamped into its border cells).
constexpr double kGridMargin = 8.0;
// ContactStorage::LiveRows: the capacities over the starting state's near rows and live rows
// (over 52 measured calls the largest growth within a
// call was 4.1 % for the near rows and 5 % for the live rows). No regrow: an overflow ends the solve.
constexpr double kLiveRowHeadroom = 1.10;
constexpr double kLiveListHeadroom = 1.15;
constexpr long long kLiveHeadroomMin = 65536;

}  // namespace

double NcpEvalSummary::maxConstraint() const {
    double m = 0;
    for (int f = 0; f < knit::kRowFamilies; ++f)
        if (knit::isHard(static_cast<knit::RowFamily>(f))) m = std::max(m, maxViol[f]);
    return m;
}

NcpDevice::State::~State() {
    if (stream) cudaStreamDestroy(stream);  // cleanup in a destructor: never throws
}

NcpDevice::~NcpDevice() = default;
int NcpDevice::vars() const { return state_->N; }
void NcpDevice::setIndexLimitForTests(long long limit) { gIndexLimitForTests = limit; }
void NcpDevice::setContactHeadroomForTests(double factor, long long minimum) {
    gHeadroomForTests = factor;
    gHeadroomMinForTests = minimum;
}
void NcpDevice::setLiveHeadroomForTests(double factor, long long minimum) {
    gLiveHeadroomForTests = factor;
    gLiveHeadroomMinForTests = minimum;
}
NcpStats& NcpDevice::stats() { return state_->st; }
NcpDevice::State& NcpDevice::state() { return *state_; }

long long NcpDevice::contactRows() const {
    const State& m = *state_;
    check(cudaStreamSynchronize(m.stream), "ncp-gpu");
    NcpCounts c;
    m.countsBuffer.download(&c, sizeof c);
    return c.contactRows;
}

void NcpDevice::State::uploadMesh(const knit::TriMesh& tm) {
    std::vector<double> pos, vn, fn;
    pos.reserve(tm.pos.size() * 3);
    for (const knit::Vec3& v : tm.pos) pos.insert(pos.end(), {v.x, v.y, v.z});
    for (const knit::Vec3& v : tm.vnormal) vn.insert(vn.end(), {v.x, v.y, v.z});
    for (size_t t = 0; t < tm.tri.size(); ++t) {
        const knit::Vec3 f = tm.faceNormal(static_cast<int>(t));
        fn.insert(fn.end(), {f.x, f.y, f.z});
    }
    if (vn.size() != pos.size()) throw DeviceError("ncp-gpu: the band mesh has no vertex normals");
    mPos.upload(pos);
    mVn.upload(vn);
    mFn.upload(fn);
    DoubleMesh& dm = bp.mesh;
    dm.pos = mPos.doubles();
    dm.vnormal = mVn.doubles();
    dm.faceN = mFn.doubles();
    dm.tri = mesh->view.tri;
    dm.grid = mesh->view.grid;
    dm.adj = mesh->adjacency;
    knit::Vec3 lo, hi;
    tm.bounds(lo, hi);
    dm.lo[0] = lo.x; dm.lo[1] = lo.y; dm.lo[2] = lo.z;
    dm.hi[0] = hi.x; dm.hi[1] = hi.y; dm.hi[2] = hi.z;
}

namespace {

// Host arrays of a call, on their way to the device.
struct ControlArrays {
    std::vector<d3> P;
    std::vector<int> varOf, ctrlOfVar;
};
struct SampleArrays {
    std::vector<int> base, yarn;
    std::vector<double> w, s;
    std::vector<unsigned char> active, segOk;
};
struct GeneralRowArrays {
    std::vector<int> fam, idx, sub, n, var;
    std::vector<double> rest;
};

// The parameters of KnitProblem.
void setParameters(NcpDevice::State& m, const knit::RelaxMoveParams& mp, const knit::KnitProblemParams& kp) {
    m.pr.D = mp.thr > 0 ? mp.thr : 1.0;
    m.pr.thr = mp.thr;
    m.pr.kappaMax = mp.kappaMax;
    m.pr.hardC = kp.hardContact > 0 ? kp.hardContact : mp.thr;
    m.pr.hardK = kp.hardKappa > 0 ? kp.hardKappa : mp.kappaMax;
    m.pr.targets = kp.hardContact > 0 || kp.hardKappa > 0;
    m.pr.bandHalf = mp.bandHalf;
    m.pr.cap = mp.maxDisplacement;
    m.pr.wL = kp.wL;
    m.pr.wT = kp.wT;
    m.pr.wA = kp.wA;
    m.pr.kDeg = kp.kDegenerate;
    m.pr.curvQ = std::clamp(kp.curvatureSamples, 2, 255);
    m.pr.contactModel = static_cast<int>(kp.contactModel);
    m.pr.parallelCos = kp.parallelCos;
    m.delta = kp.delta;
    m.window = mp.window;
}

void setBand(NcpDevice::State& m, const NcpCall& call) {
    const knit::BaseSurface* surface = call.surface;
    const knit::SurfaceBinding* bind = call.binding;
    const knit::Controls& Ph = call.controls;
    m.bp = BandParams{};
    m.hasBind = false;
    if (bind && bind->surface) {
        const knit::MeshSurface* ms = bind->surface;
        m.mesh = residentMesh(ms->mesh(), MeshWalk::Adjacency);
        m.uploadMesh(ms->mesh());
        m.bp.mode = kBandBind;
        m.bp.cont.lambdaW = bind->cp.lambdaW;
        m.bp.cont.minNormalCos = bind->cp.minNormalCos;
        m.bp.cont.band = bind->cp.band;
        m.bp.cont.maxSteps = bind->cp.maxSteps;
        if (bind->state.size() != Ph.size() || bind->result.size() != Ph.size()) throw DeviceError("ncp-gpu: binding not sized like the controls");
        m.hasBind = true;
    } else if (surface) {
        const knit::TriMesh* tm = nullptr;
        if (const auto* ps = dynamic_cast<const knit::PieceSurface*>(surface)) {
            tm = &ps->chart().mesh();
        } else if (const auto* ms = dynamic_cast<const knit::MeshSurface*>(surface)) {
            tm = &ms->mesh();
            m.bp.faceFallback = true;
        } else if (const auto* cone = dynamic_cast<const knit::ConeSurface*>(surface)) {
            m.bp.mode = kBandCone;
            m.bp.R0 = cone->radiusAt(0);
            m.bp.R1 = cone->radiusAt(1);
            m.bp.H = cone->height();
        } else if (dynamic_cast<const knit::PlaneSurface*>(surface)) {
            m.bp.mode = kBandPlane;
        } else {
            throw DeviceError("ncp-gpu: this BaseSurface type has no device projection");
        }
        if (m.bp.mode == kBandNone) {
            m.mesh = residentMesh(*tm, MeshWalk::None);
            m.uploadMesh(*tm);
            m.bp.mode = kBandMesh;
        }
    }
}

ControlArrays setControls(NcpDevice::State& m, const NcpCall& call) {
    const knit::CurveSet& cs = call.curves;
    const knit::Controls& Ph = call.controls;
    const knit::Controls& Prest = call.rest;
    const knit::ControlMask& freeCtrl = call.freeControls;
    const int Y = static_cast<int>(Ph.size());
    ControlArrays a;
    m.ctrlBase.assign(Y + 1, 0);
    for (int y = 0; y < Y; ++y) {
        if (Ph[y].size() != cs.yarns[y].curve.ctrl.size() || freeCtrl[y].size() != Ph[y].size() || Prest[y].size() != Ph[y].size())
            throw DeviceError("ncp-gpu: per-yarn control counts disagree");
        m.ctrlBase[y + 1] = m.ctrlBase[y] + static_cast<int>(Ph[y].size());
    }
    m.nCtrl = m.ctrlBase[Y];
    m.Pshape = Ph;
    a.P.resize(m.nCtrl);
    a.varOf.assign(m.nCtrl, -1);
    for (int y = 0; y < Y; ++y)
        for (size_t c = 0; c < Ph[y].size(); ++c) {
            const int g = m.ctrlBase[y] + static_cast<int>(c);
            a.P[g] = {Ph[y][c].x, Ph[y][c].y, Ph[y][c].z};
            if (freeCtrl[y][c]) {
                a.varOf[g] = static_cast<int>(a.ctrlOfVar.size());
                a.ctrlOfVar.push_back(g);
            }
        }
    m.N = static_cast<int>(a.ctrlOfVar.size());
    return a;
}

// The samples of KnitProblem's constructor.
SampleArrays sampleArrays(NcpDevice::State& m, const NcpCall& call) {
    const knit::CurveSet& cs = call.curves;
    const knit::ControlMask& freeCtrl = call.freeControls;
    const knit::RelaxMoveParams& mp = call.move;
    const knit::KnitProblemParams& kp = call.problem;
    const int Y = static_cast<int>(call.controls.size());
    SampleArrays a;
    const int sps = mp.samplesPerSpan;
    for (int y = 0; y < Y; ++y) {
        const knit::YarnPath& yp = cs.yarns[y];
        for (int k = 0; k < yp.curve.numSpans(); ++k) {
            bool act = false;
            for (int i = 0; i < 4; ++i) act = act || freeCtrl[y][k + i];
            int n = sps;  // KnitProblem's constructor: adaptive samples per span
            if (kp.maxSampleSpacing > 0)
                n = std::max(n, static_cast<int>(std::ceil((yp.arcLengthAt(k + 1.0) - yp.arcLengthAt(k)) / kp.maxSampleSpacing)));
            for (int j = 0; j < n; ++j) {
                const double t = double(j) / n, s1 = 1 - t;
                a.base.push_back(m.ctrlBase[y] + k);
                a.w.push_back(s1 * s1 * s1 / 6);
                a.w.push_back((3 * t * t * t - 6 * t * t + 4) / 6);
                a.w.push_back((-3 * t * t * t + 3 * t * t + 3 * t + 1) / 6);
                a.w.push_back(t * t * t / 6);
                a.yarn.push_back(y);
                a.s.push_back(yp.arcLengthAt(k + t));
                a.active.push_back(act ? 1 : 0);
            }
        }
    }
    m.nSamples = static_cast<int>(a.base.size());
    a.segOk.resize(m.nSamples);
    for (int i = 0; i < m.nSamples; ++i) a.segOk[i] = i + 1 < m.nSamples && a.yarn[i + 1] == a.yarn[i];
    return a;
}

// The fixed general rows: edges, curvature, band, cap, anchor.
GeneralRowArrays generalRowArrays(NcpDevice::State& m, const NcpCall& call, const ControlArrays& controls) {
    const knit::Controls& Ph = call.controls;
    const knit::Controls& Prest = call.rest;
    const knit::RelaxMoveParams& mp = call.move;
    const knit::KnitProblemParams& kp = call.problem;
    const int Y = static_cast<int>(Ph.size());
    GeneralRowArrays r;
    auto addRow = [&](int f, int g, int s, const std::vector<int>& vars, double rs) {
        r.fam.push_back(f);
        r.idx.push_back(g);
        r.sub.push_back(s);
        r.n.push_back(static_cast<int>(vars.size()));
        for (int e = 0; e < kGeneralRowEntries; ++e) r.var.push_back(e < static_cast<int>(vars.size()) ? vars[e] : -1);
        r.rest.push_back(rs);
    };
    for (int y = 0; y < Y; ++y)
        for (size_t c = 0; c + 1 < Ph[y].size(); ++c) {
            const int g = m.ctrlBase[y] + static_cast<int>(c);
            const int va = controls.varOf[g], vb = controls.varOf[g + 1];
            if (va < 0 && vb < 0) continue;
            std::vector<int> vars;
            int side = 0;  // bit e: entry e is control g (gradient -u)
            if (vb >= 0) vars.push_back(vb);
            if (va >= 0) {
                side |= 1 << static_cast<int>(vars.size());
                vars.push_back(va);
            }
            addRow(FEdge, g, side, vars, knit::norm(Prest[y][c + 1] - Prest[y][c]));
        }
    // Curvature rows q in [qBegin, qEnd) of every span with a free control.
    auto addCurvatureRows = [&](int qBegin, int qEnd) {
        for (int y = 0; y < Y; ++y) {
            const int nSp = static_cast<int>(Ph[y].size()) - 3;
            for (int k = 0; k < nSp; ++k) {
                const int g = m.ctrlBase[y] + k;
                std::vector<int> vars;
                int offs = 0;
                for (int i = 0; i < 4; ++i)
                    if (controls.varOf[g + i] >= 0) {
                        offs |= i << (8 + 2 * static_cast<int>(vars.size()));
                        vars.push_back(controls.varOf[g + i]);
                    }
                if (vars.empty()) continue;
                for (int q = qBegin; q < qEnd; ++q) addRow(FCurv, g, q | offs, vars, 0);
            }
        }
    };
    addCurvatureRows(0, m.pr.curvQ);
    const bool band = m.bp.mode != kBandNone;
    for (int a = 0; a < m.N; ++a) {
        const int g = controls.ctrlOfVar[a];
        if (band) {
            addRow(FBandU, g, 0, {a}, 0);
            addRow(FBandL, g, 0, {a}, 0);
        }
        if (mp.maxDisplacement > 0) addRow(FCap, g, 0, {a}, 0);
        if (kp.wA > 0)
            for (int ax = 0; ax < 3; ++ax) addRow(FAnchor, g, ax, {a}, 0);
    }
    // CurvatureRows::SamplesAndSpanMaximum: q = curvQ marks the row at the span's curvature maximum, which
    // evalGeneralRowsKernel locates with about 60 curvature evaluations. These rows come last, together, so the
    // warps that search hold no other row (mixed among the samples, every warp waited for its search: tube-600
    // took 405 s instead of 195 s).
    if (call.options.curvatureRows == CurvatureRows::SamplesAndSpanMaximum) addCurvatureRows(m.pr.curvQ, m.pr.curvQ + 1);
    return r;
}

// The general rows as uploaded (without their rest values), downloaded for the tests.
GeneralRowArrays downloadGeneralRows(const NcpDevice::State& m) {
    GeneralRowArrays r;
    m.general.fam.download(r.fam, m.nG);
    m.general.idx.download(r.idx, m.nG);
    m.general.sub.download(r.sub, m.nG);
    m.general.n.download(r.n, m.nG);
    m.general.var.download(r.var, size_t(m.nG) * kGeneralRowEntries);
    return r;
}

// Row keys (ncp_reference.cpp rowKey) for the tests' comparison with the CPU rows.
std::vector<uint64_t> generalRowKeys(const NcpDevice::State& m, const GeneralRowArrays& r) {
    std::vector<int> yOf(m.nCtrl), cOf(m.nCtrl);
    for (size_t y = 0; y + 1 < m.ctrlBase.size(); ++y)
        for (int g = m.ctrlBase[y]; g < m.ctrlBase[y + 1]; ++g) {
            yOf[g] = static_cast<int>(y);
            cOf[g] = g - m.ctrlBase[y];
        }
    auto key = [](int f, uint64_t a, uint64_t b) { return (uint64_t(f) << 60) ^ (a << 30) ^ b; };
    std::vector<uint64_t> keys(m.nG);
    for (int k = 0; k < m.nG; ++k) {
        const int g = r.idx[k], y = yOf[g], c = cOf[g];
        switch (r.fam[k]) {
        case FCurv: keys[k] = key(FCurv, uint64_t(y) << 20 | uint64_t(c), uint64_t(r.sub[k] & 0xff)); break;
        case FAnchor: keys[k] = key(FAnchor, uint64_t(y) << 2 | uint64_t(r.sub[k]), uint64_t(c)); break;
        default: keys[k] = key(r.fam[k], uint64_t(y), uint64_t(c)); break;
        }
    }
    return keys;
}

void uploadBinding(NcpDevice::State& m, const knit::Controls& Ph, const knit::SurfaceBinding& bind) {
    const int Y = static_cast<int>(Ph.size());
    std::vector<DoubleProjectionState> hs(m.nCtrl);
    std::vector<unsigned char> hr(m.nCtrl);
    for (int y = 0; y < Y; ++y)
        for (size_t c = 0; c < Ph[y].size(); ++c) {
            const knit::ProjState& s = bind.state[y][c];
            DoubleProjectionState& d = hs[m.ctrlBase[y] + c];
            d.tri = s.tri;
            for (int k = 0; k < 3; ++k) d.b[k] = s.b[k];
            d.foot = {s.foot.x, s.foot.y, s.foot.z};
            d.n = {s.n.x, s.n.y, s.n.z};
            d.query = {s.query.x, s.query.y, s.query.z};
            hr[m.ctrlBase[y] + c] = static_cast<unsigned char>(bind.result[y][c]);
        }
    m.bindState.upload(hs.data(), hs.size() * sizeof(hs[0]));
    m.bindResult.upload(hr);
}

// The grid's fixed shape: the cell of KnitProblem::refresh (pad + 1), over the starting samples' box
// grown by kGridMargin. The sample points come from the host controls (the same B-spline weights).
void fixGridShape(NcpDevice::State& m, const ControlArrays& controls, const SampleArrays& samples) {
    const float pad = static_cast<float>(m.pr.thr + 2 * m.delta) + 1e-3f;
    FixedGrid& g = m.grid;
    g.cell = pad + 1.0f;
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (int i = 0; i < m.nSamples; ++i) {
        const int b = samples.base[i];
        const double* w = &samples.w[4 * size_t(i)];
        const d3 p = controls.P[b] * w[0] + controls.P[b + 1] * w[1] + controls.P[b + 2] * w[2] + controls.P[b + 3] * w[3];
        const double c[3] = {p.x, p.y, p.z};
        for (int a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], c[a]);
            hi[a] = std::max(hi[a], c[a]);
        }
    }
    if (m.nSamples == 0) lo[0] = lo[1] = lo[2] = hi[0] = hi[1] = hi[2] = 0;
    for (;;) {
        long long cells = 1;
        for (int a = 0; a < 3; ++a) {
            g.lo[a] = static_cast<long long>(std::floor((lo[a] - kGridMargin) / g.cell));
            g.n[a] = static_cast<long long>(std::floor((hi[a] + kGridMargin) / g.cell)) - g.lo[a] + 1;
            cells *= g.n[a];
        }
        if (cells <= (1ll << 24)) break;
        g.cell *= 1.25f;
    }
    const long long nCells = g.n[0] * g.n[1] * g.n[2];
    g.endBit = 32;
    while (g.endBit < 64 && (1ull << (g.endBit - 32)) <= static_cast<unsigned long long>(nCells)) ++g.endBit;
    // A segment's box is shorter than a cell (samples are a fraction of a span apart), so it
    // registers in at most 2 x 2 x 2 cells; a longer one is an explicit overflow.
    m.capacity.registrations = 8LL * m.nSamples + 1024;
    g.boxLo.alloc(size_t(m.nSamples) * 12 + 12);
    g.boxHi.alloc(size_t(m.nSamples) * 12 + 12);
    g.counts.alloc(size_t(m.nSamples + 1) * 8);
    g.offsets.alloc(size_t(m.nSamples + 1) * 8);
    g.keys.alloc(size_t(m.capacity.registrations) * 8);
    g.sorted.alloc(size_t(m.capacity.registrations) * 8);
    g.entries.alloc(size_t(m.capacity.registrations) * 4);
    g.cellCounts.alloc(size_t(nCells + 1) * 4);
    g.cellStart.alloc(size_t(nCells + 1) * 4);
    g.view = GridView{};
    g.view.cell = g.cell;
    for (int a = 0; a < 3; ++a) {
        g.view.lo[a] = g.lo[a];
        g.view.hi[a] = g.lo[a] + g.n[a] - 1;
    }
    g.view.cellStart = g.cellStart.uints();
    g.view.entries = g.entries.uints();
    g.view.boxLo = g.boxLo.floats();
    g.view.boxHi = g.boxHi.floats();
    g.view.nItems = m.nSamples;
}

// One row value set: its general rows' values and terms, and (set 0 only) the contact rows' storage,
// which both sets share.
void allocateValueSet(NcpDevice::State& m, RowValueSet& b, const RowValueSet* shareContact) {
    const size_t nG = static_cast<size_t>(m.nG), C = static_cast<size_t>(m.contactListCapacity());
    b.gVal.alloc(nG * 8 + 8);
    b.gTval.alloc(nG * 8 + 8);
    b.gAux.alloc(nG * 8 + 8);
    b.gGrad.alloc(nG * kGeneralRowEntries * sizeof(d3) + 8);
    b.gFlag.alloc(nG + 8);
    b.wres.alloc(nG * 8 + 8);
    b.wgn.alloc(nG * 8 + 8);
    b.merit.alloc(nG * 8 + 8);
    b.obj.alloc(nG * 8 + 8);
    b.strain2.alloc(nG * 8 + 8);
    if (!shareContact) {
        if (m.storage == ContactStorage::LiveRows) {
            b.cDist.alloc(C * 8 + 8);  // the values follow from the distance (contactHardValue, contactTargetValue)
        } else {
            b.cVal.alloc(C * 8 + 8);
            b.cTval.alloc(C * 8 + 8);
        }
        b.cWres.alloc(C * 8 + 8);
        b.cWgn.alloc(C * 8 + 8);
        b.cGradOn.alloc(C + 8);
    }
    const RowValueSet& c = shareContact ? *shareContact : b;
    b.view = {b.gVal.doubles(), b.gTval.doubles(), b.gAux.doubles(), static_cast<d3*>(b.gGrad.data()), b.gFlag.chars(), b.wres.doubles(), b.wgn.doubles(),
              b.merit.doubles(), b.obj.doubles(), b.strain2.doubles(), c.cVal.doubles(), c.cTval.doubles(), c.cWres.doubles(), c.cWgn.doubles(), c.cGradOn.chars(),
              c.cDist.doubles()};
}

}  // namespace

void NcpDevice::State::allocateValueSets() {
    allocateValueSet(*this, values[0], nullptr);
    allocateValueSet(*this, values[1], &values[0]);
}

// Every buffer whose size follows the contact rows, at this capacity. The constructor sizes them (the
// starting state's rows with headroom); the iteration never allocates. The active slots are one per row
// (general and contact), so the active lists cannot overflow; the active entries' memory also holds the
// multiplier carry-over's sort during a refresh (carryBuffers).
void NcpDevice::State::allocateCapacity(long long contactRows) {
    capacity.contactRows = contactRows;
    const long long C = capacity.contactRows, R = nG + C;
    checkIndexLimit(R, "rows");
    checkIndexLimit(4LL * nG + kContactRowEntries * C, "row entries");
    if (R >= (1LL << 27)) throw DeviceError("ncp-gpu: more than 2^27 rows (active list entries hold slot << 4)");
    capacity.activeEntries = 4LL * nG + kContactRowEntries * C;
    contacts.pair.alloc(size_t(C) * sizeof(int2) + 8);
    contacts.type.alloc(size_t(C) + 8);
    contacts.lambda.alloc(size_t(C) * 8 + 8);
    allocateValueSets();
    active.slotOf.alloc(size_t(R + 1) * 4 + 8);
    active.rows.alloc(size_t(R) * 4 + 8);
    active.count.alloc(size_t(N + 1) * 4);
    active.start.alloc(size_t(N + 1) * 4);
    active.cursor.alloc(size_t(N + 1) * 4);
    active.entry.alloc(std::max(size_t(capacity.activeEntries) * 4, size_t(C) * 32) + 8);  // the carry-over: 4 arrays of C x 8 bytes
    active.normal.alloc(size_t(R) * sizeof(d3) + 8);
    active.st.alloc(size_t(R) * 16 + 8);
    system.rowT.alloc(size_t(R) * 8 + 8);
    system.rowS.alloc(size_t(nG) * sizeof(d3) + 8);
    const long long nCells = grid.n[0] * grid.n[1] * grid.n[2];
    prims.reserveScratch(std::max({capacity.registrations, C + 1, R + 1, nCells + 1, static_cast<long long>(N) + 1, static_cast<long long>(nSamples) + 1}));
}

// Mirrors allocateCapacity, buffer by buffer.
long long NcpDevice::State::capacityBytes(long long contactRows) const {
    const long long C = contactRows, R = nG + C, G = nG;
    const long long activeEntries = 4LL * nG + kContactRowEntries * C;
    long long b = (C * 8 + 8) + (C + 8) + (C * 8 + 8);                                       // pair, type, lambda
    b += 2 * (8 * (G * 8 + 8) + (G * kGeneralRowEntries * 24 + 8) + (G + 8));                // general values and terms, 2 sets
    b += 4 * (C * 8 + 8) + (C + 8);                                                          // contact values and terms
    b += ((R + 1) * 4 + 8) + (R * 4 + 8) + 3 * (static_cast<long long>(N) + 1) * 4;          // slots, rows, count, start, cursor
    b += std::max(activeEntries * 4, C * 32) + 8 + (R * 24 + 8) + (R * 16 + 8);              // entries, normals, s t
    b += (R * 8 + 8) + (G * 24 + 8);                                                         // rowT, rowS
    const long long nCells = grid.n[0] * grid.n[1] * grid.n[2];
    const long long items = std::max({capacity.registrations, C + 1, R + 1, nCells + 1, static_cast<long long>(N) + 1, static_cast<long long>(nSamples) + 1});
    const long long scratch = static_cast<long long>(DevicePrimitives::scratchBytesFor(items)) - static_cast<long long>(prims.scratchBytes());
    return b + std::max(0LL, scratch);
}

long long NcpDevice::State::capacityBytesNow() const {
    long long b = 0;
    for (const DeviceBuffer* x : {&contacts.pair, &contacts.type, &contacts.lambda, &active.slotOf, &active.rows, &active.count, &active.start, &active.cursor,
                                  &active.entry, &active.normal, &active.st, &system.rowT, &system.rowS})
        b += static_cast<long long>(x->bytes());
    for (const RowValueSet* v : {&values[0], &values[1]})
        for (const DeviceBuffer* x : {&v->gVal, &v->gTval, &v->gAux, &v->gGrad, &v->gFlag, &v->wres, &v->wgn, &v->merit, &v->obj, &v->strain2, &v->cVal, &v->cTval,
                                      &v->cWres, &v->cWgn, &v->cGradOn})
            b += static_cast<long long>(x->bytes());
    return b;
}

long long NcpDevice::State::chooseCapacity(long long rows) const {
    const double headroom = gHeadroomForTests >= 0 ? gHeadroomForTests : kContactHeadroom;
    const long long headroomMin = gHeadroomForTests >= 0 ? gHeadroomMinForTests : kContactHeadroomMin;
    const long long desired = static_cast<long long>(headroom * rows) + headroomMin;
    const size_t limit = AllocationCap::limit();
    if (limit == 0) return desired;
    // What the cap leaves once the capacity buffers in place now are given back.
    const long long available = static_cast<long long>(limit) - (static_cast<long long>(memStats().current) - capacityBytesNow());
    if (capacityBytes(desired) <= available) return desired;
    long long lo = std::min(desired, static_cast<long long>(std::ceil(kCapacityFloor * rows)));
    if (capacityBytes(lo) > available) return lo;  // refused by the cap when allocated
    long long hi = desired;
    while (hi - lo > 1) {  // the largest capacity that fits
        const long long mid = lo + (hi - lo) / 2;
        (capacityBytes(mid) <= available ? lo : hi) = mid;
    }
    return lo;
}

void NcpDevice::State::regrowContactRows(long long needed) {
    check(cudaStreamSynchronize(stream), "ncp-gpu regrow");
    NcpCounts c;
    countsBuffer.download(&c, sizeof c);
    const size_t rows = static_cast<size_t>(c.contactRows);
    std::vector<int2> pair(rows);
    std::vector<unsigned char> type(rows);
    std::vector<double> lam(rows);
    contacts.pair.download(pair.data(), rows * sizeof(int2));
    contacts.type.download(type.data(), rows);
    contacts.lambda.download(lam.data(), rows * 8);
    allocateCapacity(std::max(chooseCapacity(needed), needed));  // at least the rows the refresh needs
    check(cudaMemcpy(contacts.pair.data(), pair.data(), rows * sizeof(int2), cudaMemcpyHostToDevice), "ncp-gpu regrow");
    check(cudaMemcpy(contacts.type.data(), type.data(), rows, cudaMemcpyHostToDevice), "ncp-gpu regrow");
    check(cudaMemcpy(contacts.lambda.data(), lam.data(), rows * 8, cudaMemcpyHostToDevice), "ncp-gpu regrow");
    c.overflow = kCapNone;
    c.overflowNeeded = c.overflowCapacity = 0;
    countsBuffer.upload(&c, sizeof c);
    st.capacityContactRows = capacity.contactRows;
}

std::vector<NcpMemoryFamily> NcpDevice::State::memoryFamilies() const {
    const auto sum = [](std::initializer_list<const DeviceBuffer*> buffers) {
        long long b = 0;
        for (const DeviceBuffer* buffer : buffers) b += static_cast<long long>(buffer->bytes());
        return b;
    };
    const RowValueSet &v0 = values[0], &v1 = values[1];
    std::vector<NcpMemoryFamily> f;
    f.push_back({"controls", sum({&controls.P, &controls.Q, &controls.P0, &controls.Xref, &controls.varOf, &controls.ctrlOfVar, &controls.pin})});
    f.push_back({"samples", sum({&samples.base, &samples.w, &samples.yarn, &samples.s, &samples.active, &samples.segOk, &samplesAtP, &samplesAtQ,
                                 &contacts.sampleRows, &contacts.sampleRowOffset})});
    f.push_back({"grid", sum({&grid.boxLo, &grid.boxHi, &grid.counts, &grid.offsets, &grid.keys, &grid.sorted, &grid.cellCounts, &grid.cellStart, &grid.entries})});
    f.push_back({"general rows", sum({&general.fam, &general.idx, &general.sub, &general.n, &general.var, &general.rest, &general.lambda})});
    f.push_back({"general values and terms (2 sets)", sum({&v0.gVal, &v0.gTval, &v0.gAux, &v0.gGrad, &v0.gFlag, &v0.wres, &v0.wgn, &v0.merit, &v0.obj,
                                                           &v0.strain2, &v1.gVal, &v1.gTval, &v1.gAux, &v1.gGrad, &v1.gFlag, &v1.wres, &v1.wgn, &v1.merit,
                                                           &v1.obj, &v1.strain2})});
    f.push_back({"contact rows", sum({&contacts.pair, &contacts.type, &contacts.lambda})});
    if (storage == ContactStorage::LiveRows) f.push_back({"live list", sum({&liveList.slot, &liveList.live, &liveList.row})});
    f.push_back({"contact values and terms", sum({&v0.cVal, &v0.cTval, &v0.cWres, &v0.cWgn, &v0.cGradOn, &v0.cDist})});
    f.push_back({"active slots", sum({&active.slotOf, &active.rows, &active.normal, &active.st, &active.count, &active.start, &active.cursor})});
    f.push_back({"active entries and carry-over", sum({&active.entry})});
    f.push_back({"linear system", sum({&system.grad, &system.gnorm, &system.hdiag, &system.hblock, &system.rowT, &system.rowS, &system.x, &system.r, &system.z,
                                       &system.p, &system.q, &system.precond, &system.partials})});
    f.push_back({"reductions", sum({&reductions.evaluation, &reductions.acceptEvaluation, &reductions.blocks, &reductions.regained, &reductions.perVariable,
                                    &reductions.perVariable2, &countsBuffer, &controlBuffer})});
    f.push_back({"band", sum({&bandFootBuffer, &bandNrmBuffer, &bandOkBuffer, &bindState, &bindResult, &mPos, &mVn, &mFn})});
    f.push_back({"scratch", static_cast<long long>(prims.scratchBytes())});
    return f;
}

NcpDevice::NcpDevice(const NcpCall& call, NcpSizing sizing) : state_(std::make_unique<State>()) {
    State& m = *state_;
    m.capLimit = gIndexLimitForTests;
    m.storage = call.options.contactStorage;
    const DeviceProperties& device = deviceProperties();
    if (!device.ok) throw DeviceError(device.error);
    const size_t yarns = call.controls.size();
    if (yarns != call.curves.yarns.size() || call.freeControls.size() != yarns || call.rest.size() != yarns)
        throw DeviceError("ncp-gpu: controls, mask and curve set disagree");
    check(cudaStreamCreate(&m.stream), "ncp-gpu stream");
    setParameters(m, call.move, call.problem);
    setBand(m, call);
    const ControlArrays controls = setControls(m, call);
    const SampleArrays samples = sampleArrays(m, call);
    const GeneralRowArrays rows = generalRowArrays(m, call, controls);
    m.nG = static_cast<int>(rows.fam.size());
    std::vector<double> zeros(m.nG, 0.0);
    m.controls.P.upload(controls.P.data(), controls.P.size() * sizeof(d3));
    m.controls.Q.upload(controls.P.data(), controls.P.size() * sizeof(d3));
    m.controls.P0.upload(controls.P.data(), controls.P.size() * sizeof(d3));
    m.controls.Xref.upload(controls.P.data(), controls.P.size() * sizeof(d3));
    m.controls.varOf.upload(controls.varOf);
    m.controls.ctrlOfVar.upload(controls.ctrlOfVar);
    m.samples.base.upload(samples.base);
    m.samples.w.upload(samples.w);
    m.samples.yarn.upload(samples.yarn);
    m.samples.s.upload(samples.s);
    m.samples.active.upload(samples.active);
    m.samples.segOk.upload(samples.segOk);
    m.general.fam.upload(rows.fam);
    m.general.idx.upload(rows.idx);
    m.general.sub.upload(rows.sub);
    m.general.n.upload(rows.n);
    m.general.var.upload(rows.var);
    m.general.rest.upload(rows.rest);
    m.general.lambda.upload(zeros);
    if (m.hasBind) uploadBinding(m, call.controls, *call.binding);
    m.S = {m.samples.base.ints(), m.samples.w.doubles(), m.samples.yarn.ints(), m.samples.s.doubles(), m.samples.active.chars(), m.samples.segOk.chars(),
           m.nSamples};
    m.controls.pin.alloc(size_t(m.N) + 8);
    check(cudaMemset(m.controls.pin.data(), 0, size_t(m.N) + 8), "ncp-gpu setup");
    m.countsBuffer.alloc(sizeof(NcpCounts));
    check(cudaMemset(m.countsBuffer.data(), 0, sizeof(NcpCounts)), "ncp-gpu setup");
    m.controlBuffer.alloc(sizeof(NcpControl));
    check(cudaMemset(m.controlBuffer.data(), 0, sizeof(NcpControl)), "ncp-gpu setup");
    m.reductions.evaluation.alloc(sizeof(EvalReductions));
    m.reductions.acceptEvaluation.alloc(sizeof(EvalReductions));
    m.reductions.blocks.alloc(4 * kSumBlocks * sizeof(double));
    m.reductions.regained.alloc(size_t(m.nG) * 8 + 8);
    m.reductions.perVariable.alloc(size_t(m.N) * 8 + 8);
    m.reductions.perVariable2.alloc(size_t(m.N) * 8 + 8);
    m.bandFootBuffer.alloc(size_t(m.nCtrl) * sizeof(d3) + 8);
    m.bandNrmBuffer.alloc(size_t(m.nCtrl) * sizeof(d3) + 8);
    m.bandOkBuffer.alloc(size_t(m.nCtrl) + 8);
    LinearSystem& s = m.system;
    for (DeviceBuffer* b : {&s.grad, &s.hdiag, &s.x, &s.r, &s.z, &s.p, &s.q}) b->alloc(size_t(3 * m.N) * 8 + 8);
    s.gnorm.alloc(size_t(m.N) * 8 + 8);
    s.hblock.alloc(size_t(m.N) * 6 * 8 + 8);
    s.precond.alloc(size_t(m.N) * 6 * 8 + 8);
    s.partials.alloc(4 * kSumBlocks * sizeof(double));
    m.samplesAtP.alloc(size_t(m.nSamples) * sizeof(d3) + 8);
    m.samplesAtQ.alloc(size_t(m.nSamples) * sizeof(d3) + 8);
    m.contacts.sampleRows.alloc(size_t(m.nSamples + 1) * 8);
    m.contacts.sampleRowOffset.alloc(size_t(m.nSamples + 1) * 8);
    fixGridShape(m, controls, samples);
    checkLaunch("ncp-gpu: init");
    m.st.samples = m.nSamples;
    m.st.generalRows = m.nG;
    long long countedRows = 0;
    // The starting state's sizes for the memory report (host preparation: the counts are read).
    if (m.N > 0) {
        const float pad = static_cast<float>(m.pr.thr + 2 * m.delta) + 1e-3f;
        const long long nCells = m.grid.n[0] * m.grid.n[1] * m.grid.n[2];
        m.prims.reserveScratch(std::max({m.capacity.registrations, nCells + 1, static_cast<long long>(m.nSamples) + 1}));
        DeviceBuffer totals;
        totals.alloc(4 * sizeof(unsigned long long));
        m.enqueueSamplePoints(RowSet::Current);
        m.enqueueGrid();
        m.enqueueStartCounts(pad, totals.ulongs());
        check(cudaStreamSynchronize(m.stream), "ncp-gpu start counts");
        NcpCounts c;
        m.countsBuffer.download(&c, sizeof c);
        if (c.overflow == kCapRegistrations)
            throw DeviceError("ncp-gpu: a sample segment is longer than the broad phase's grid cell (" + std::to_string(c.overflowNeeded) + " registrations)");
        unsigned long long t[4];
        totals.download(t, sizeof t);
        m.st.initialCandidates = static_cast<long long>(t[0]);
        m.st.initialActiveContactRows = static_cast<long long>(t[2]);
        m.st.initialNearRows = static_cast<long long>(t[3]);
        m.st.registrations = static_cast<long long>(c.registrations);
        check(cudaMemset(m.countsBuffer.data(), 0, sizeof(NcpCounts)), "ncp-gpu setup");
        countedRows = static_cast<long long>(t[1]);
        if (sizing == NcpSizing::CountOnly) m.st.initialContactRows = countedRows;
    }
    if (sizing == NcpSizing::CountOnly) {
        m.st.memory = m.memoryFamilies();
        m.st.setupPeakMB = memStats().peak / 1048576.0;
        return;
    }
    // The capacity from the starting state's rows (counted above) with headroom. The first refresh, run
    // here once, must write exactly the counted rows: the count and the fill follow one rule.
    if (m.N > 0) {
        if (m.storage == ContactStorage::LiveRows) {
            // The rows are the near ones (no multiplier is carried at the start); the live rows at the start
            // are those with a term at multipliers 0.
            countedRows = m.st.initialNearRows;
            const double liveHeadroom = gLiveHeadroomForTests >= 0 ? gLiveHeadroomForTests : kLiveListHeadroom;
            const long long liveMin = gLiveHeadroomForTests >= 0 ? gLiveHeadroomMinForTests : kLiveHeadroomMin;
            m.allocateLiveCapacity(static_cast<long long>(std::ceil(kLiveRowHeadroom * countedRows)) + kLiveHeadroomMin,
                                   static_cast<long long>(std::ceil(liveHeadroom * m.st.initialActiveContactRows)) + liveMin);
            m.st.capacityLiveRows = m.capacity.liveRows;
        } else {
            m.allocateCapacity(m.chooseCapacity(countedRows));
        }
        m.st.initialContactRows = countedRows;
        const float pad = static_cast<float>(m.pr.thr + 2 * m.delta) + 1e-3f;
        m.enqueueContactRows(pad);
        check(cudaStreamSynchronize(m.stream), "ncp-gpu setup");
        NcpCounts c;
        m.countsBuffer.download(&c, sizeof c);
        // A capacity below the starting rows (tests) overflows here as it will in the solve's first refresh.
        const long long written = c.overflow == kCapContactRows ? c.overflowNeeded : c.overflow == kCapNone ? c.contactRows : -1;
        if (written != countedRows)
            throw DeviceError("ncp-gpu: the starting state's row count (" + std::to_string(countedRows) + ") disagrees with its contact rows (" +
                              std::to_string(written) + ", overflow " + std::to_string(c.overflow) + ")");
        check(cudaMemset(m.countsBuffer.data(), 0, sizeof(NcpCounts)), "ncp-gpu setup");
    } else {
        m.allocateCapacity(0);
    }
    m.st.capacityContactRows = m.capacity.contactRows;
    m.st.memory = m.memoryFamilies();
    m.st.setupPeakMB = memStats().peak / 1048576.0;
}

void NcpDevice::setContactCapacityForTests(long long rows) { state_->capacity.contactRows = rows; }

std::vector<uint64_t> NcpDevice::contactKeys() {
    State& m = *state_;
    const size_t n = static_cast<size_t>(contactRows());
    std::vector<int2> pair(n);
    std::vector<unsigned char> type(n);
    m.contacts.pair.download(pair.data(), n * sizeof(int2));
    m.contacts.type.download(type.data(), n);
    std::vector<uint64_t> keys(n);
    for (size_t c = 0; c < n; ++c)
        keys[c] = (2ull << 60) ^ (static_cast<uint64_t>(pair[c].x) << 30) ^ static_cast<uint64_t>(pair[c].y) ^ (static_cast<uint64_t>(type[c]) << 57);
    return keys;
}

void NcpDevice::overwriteContactRowsForTests(long long first, long long count) {
    State& m = *state_;
    check(cudaStreamSynchronize(m.stream), "ncp-gpu test rows");
    int2 pair;
    unsigned char type;
    m.contacts.pair.download(&pair, sizeof pair, size_t(first - 1) * sizeof(int2));
    m.contacts.type.download(&type, 1, size_t(first - 1));
    for (long long c = first; c < first + count; ++c) {
        check(cudaMemcpy(static_cast<int2*>(m.contacts.pair.data()) + c, &pair, sizeof pair, cudaMemcpyHostToDevice), "ncp-gpu test rows");
        check(cudaMemcpy(m.contacts.type.chars() + c, &type, 1, cudaMemcpyHostToDevice), "ncp-gpu test rows");
    }
}

void NcpDevice::downloadResult(knit::Controls& Ph, knit::SurfaceBinding* bind) {
    State& m = *state_;
    check(cudaStreamSynchronize(m.stream), "ncp-gpu download");
    std::vector<d3> hP(m.nCtrl);
    m.controls.P.download(hP.data(), hP.size() * sizeof(d3));
    Ph = m.Pshape;
    for (size_t y = 0; y < Ph.size(); ++y)
        for (size_t c = 0; c < Ph[y].size(); ++c) {
            const d3& v = hP[m.ctrlBase[y] + c];
            Ph[y][c] = {v.x, v.y, v.z};
        }
    if (bind && m.hasBind) {
        std::vector<DoubleProjectionState> hs(m.nCtrl);
        std::vector<unsigned char> hr(m.nCtrl);
        m.bindState.download(hs.data(), hs.size() * sizeof(DoubleProjectionState));
        m.bindResult.download(hr.data(), hr.size());
        for (size_t y = 0; y < Ph.size(); ++y)
            for (size_t c = 0; c < Ph[y].size(); ++c) {
                const int g = m.ctrlBase[y] + static_cast<int>(c);
                const DoubleProjectionState& d = hs[g];
                knit::ProjState& s = bind->state[y][c];
                s.tri = d.tri;
                for (int k = 0; k < 3; ++k) s.b[k] = d.b[k];
                s.foot = {d.foot.x, d.foot.y, d.foot.z};
                s.n = {d.n.x, d.n.y, d.n.z};
                s.query = {d.query.x, d.query.y, d.query.z};
                bind->result[y][c] = static_cast<knit::ProjResult>(hr[g]);
            }
    }
}

// --- single operations for the tests (each synchronizes and reads its result)

void NcpDevice::refresh() {
    State& m = *state_;
    m.enqueueRefresh();
    check(cudaStreamSynchronize(m.stream), "ncp-gpu refresh");
    NcpCounts c;
    m.countsBuffer.download(&c, sizeof c);
    if (c.overflow) throw DeviceError("ncp-gpu: contact capacity exceeded (needed " + std::to_string(c.overflowNeeded) + ")");
    if (c.contactEntryOverflow) throw DeviceError("ncp-gpu: a contact row touches more than 10 free controls");
}

double NcpDevice::motionSinceRefresh() {
    State& m = *state_;
    m.enqueueMotion();
    check(cudaStreamSynchronize(m.stream), "ncp-gpu motion");
    double h = 0;
    m.controlBuffer.download(&h, 8, offsetof(NcpControl, motion));
    return h;
}

NcpEvalSummary NcpDevice::linearize(ControlSet which, double rho) {
    State& m = *state_;
    check(cudaMemcpyAsync(reinterpret_cast<char*>(m.control()) + offsetof(NcpControl, rho), &rho, 8, cudaMemcpyHostToDevice, m.stream), "ncp-gpu evaluate");
    const RowSet set = which == ControlSet::Trial ? RowSet::Trial : RowSet::Current;
    m.enqueueLinearize(set);
    check(cudaStreamSynchronize(m.stream), "ncp-gpu evaluate");
    NcpEvalSummary out;
    m.controlBuffer.download(&out, sizeof out, set == RowSet::Trial ? offsetof(NcpControl, trialEval) : offsetof(NcpControl, currentEval));
    return out;
}

SweptResult NcpDevice::swept(double eps) {
    State& m = *state_;
    check(cudaMemcpyAsync(reinterpret_cast<char*>(m.control()) + offsetof(NcpControl, sweptEps), &eps, 8, cudaMemcpyHostToDevice, m.stream), "ncp-gpu swept");
    m.enqueueSwept();
    check(cudaStreamSynchronize(m.stream), "ncp-gpu swept");
    unsigned long long h[3];
    m.controlBuffer.download(h, sizeof h, offsetof(NcpControl, swept));
    return {static_cast<long long>(h[0]), static_cast<long long>(h[1])};
}

void NcpDevice::accept() {
    State& m = *state_;
    m.enqueueAccept();
    check(cudaStreamSynchronize(m.stream), "ncp-gpu accept");
}

void NcpDevice::setTrial(const knit::Controls& Qh) {
    State& m = *state_;
    std::vector<d3> h(m.nCtrl);
    for (size_t y = 0; y < Qh.size(); ++y)
        for (size_t c = 0; c < Qh[y].size(); ++c) h[m.ctrlBase[y] + c] = {Qh[y][c].x, Qh[y][c].y, Qh[y][c].z};
    check(cudaStreamSynchronize(m.stream), "ncp-gpu trial");
    m.controls.Q.upload(h.data(), h.size() * sizeof(h[0]));
    m.enqueueSamplePoints(RowSet::Trial);
}

namespace {

// Tests (debugRows): every contact row's controls, coefficients, normal and key as KnitProblem::linearize
// gives them; the normal 0 where the row had no term in the merit when evaluated.
struct DebugContacts {
    const int* rows;  // the rows to write (null: every row, in row order)
    int count;        // with rows: how many
    int* n;
    int* ctrl;
    double* K;
    d3* normal;
    unsigned long long* key;
    double* value;
    double* tval;
};
__global__ void debugContactRowsKernel(NcpView v, RowValuePair pair, const d3* points, DebugContacts out) {
    const NcpRows& R = v.rows;
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const int n = out.rows ? out.count : *R.nC;
    for (long long kk = loopStart(); kk < n; kk += loopStride()) {
        const int k = static_cast<int>(kk);
        const int c = out.rows ? out.rows[k] : k;
        const int row = contactRowOf(R, c);
        const int2 p = R.cPair[row];
        int* C = out.ctrl + kContactRowEntries * k;
        bool overflow = false;
        const int ne = contactControls(v.samples, v.varOf, p.x, p.y, C, overflow);
        d3 nrm{0, 0, 0};
        double s = 0, t = 0;
        if (V.cGradOn[c]) {
            d3 pa, pb;
            const double d = contactDistance(points, p.x, p.y, R.cType[row], pa, pb, s, t);
            nrm = d > 1e-9 ? (pa - pb) / d : d3{0, 0, 1};
        }
        for (int e = 0; e < kContactRowEntries; ++e) out.K[kContactRowEntries * k + e] = e < ne ? contactCoefficient(v.samples, p.x, p.y, s, t, C[e], v.params.D) : 0;
        out.n[k] = ne;
        out.normal[k] = nrm;
        out.key[k] = contactKey(p.x, p.y, R.cType[row]);
        out.value[k] = contactHardValue(V, contactValueScale(v.params), c);
        out.tval[k] = contactTargetValue(V, contactValueScale(v.params), c);
    }
}

// The contact rows touching variable a (any order), at most `capacity` of them; their number in *count.
__global__ void contactRowsTouchingKernel(NcpView v, int a, int capacity, int* rows, int* count) {
    const NcpRows& R = v.rows;
    const int nC = *R.nC;
    for (long long cc = loopStart(); cc < nC; cc += loopStride()) {
        const int c = static_cast<int>(cc);
        const int2 p = R.cPair[contactRowOf(R, c)];
        int controls[kContactRowEntries];
        bool overflow = false;
        const int n = contactControls(v.samples, v.varOf, p.x, p.y, controls, overflow);
        bool touches = false;
        for (int e = 0; e < n; ++e) touches = touches || v.varOf[controls[e]] == a;
        if (!touches) continue;
        const int at = atomicAdd(count, 1);
        if (at < capacity) rows[at] = c;
    }
}

}  // namespace

std::vector<NcpDevice::DebugRow> NcpDevice::debugRows(int onlyVariable) {
    State& m = *state_;
    check(cudaStreamSynchronize(m.stream), "ncp-gpu debug rows");
    NcpCounts counts;
    m.countsBuffer.download(&counts, sizeof counts);
    int cur = 0;
    m.controlBuffer.download(&cur, sizeof cur, offsetof(NcpControl, current));
    const RowValueSet& V = m.values[cur];
    // The contact rows to download: every row, or those touching onlyVariable.
    constexpr int kTouchingCapacity = 1 << 20;
    DeviceBuffer list, listCount;
    // LiveRows: the live rows only (the rows that are not live hold no values).
    const int contactListLength = m.storage == ContactStorage::LiveRows ? counts.liveRows : counts.contactRows;
    int nC = onlyVariable == kDebugGeneralRowsOnly ? 0 : contactListLength;
    if (onlyVariable >= 0) {
        list.alloc(size_t(kTouchingCapacity) * 4 + 8);
        listCount.alloc(8);
        check(cudaMemset(listCount.data(), 0, 8), "ncp-gpu debug rows");
        contactRowsTouchingKernel<<<kLoopBlocks, kThreadsPerBlock, 0, m.stream>>>(m.view(), onlyVariable, kTouchingCapacity, list.ints(), listCount.ints());
        checkLaunch("ncp-gpu: debug rows");
        check(cudaStreamSynchronize(m.stream), "ncp-gpu debug rows");
        listCount.download(&nC, sizeof nC);
        if (nC > kTouchingCapacity) throw DeviceError("ncp-gpu debug rows: more than 2^20 rows touch the variable");
    }
    DeviceBuffer dn, dctrl, dK, dnormal, dkey, dval, dtval;
    dn.alloc(size_t(nC) * 4 + 8);
    dctrl.alloc(size_t(nC) * kContactRowEntries * 4 + 8);
    dK.alloc(size_t(nC) * kContactRowEntries * 8 + 8);
    dnormal.alloc(size_t(nC) * sizeof(d3) + 8);
    dkey.alloc(size_t(nC) * 8 + 8);
    dval.alloc(size_t(nC) * 8 + 8);
    dtval.alloc(size_t(nC) * 8 + 8);
    if (nC) {
        const DebugContacts out{onlyVariable >= 0 ? list.ints() : nullptr, nC, dn.ints(), dctrl.ints(), dK.doubles(), static_cast<d3*>(dnormal.data()),
                                dkey.ulongs(), dval.doubles(), dtval.doubles()};
        debugContactRowsKernel<<<kLoopBlocks, kThreadsPerBlock, 0, m.stream>>>(m.view(), m.valuePair(), m.pointsAtP(), out);
        checkLaunch("ncp-gpu: debug rows");
        check(cudaStreamSynchronize(m.stream), "ncp-gpu debug rows");
    }
    std::vector<double> gVal(m.nG), gTval(m.nG), cVal(nC), cTval(nC), cK(size_t(nC) * kContactRowEntries);
    std::vector<d3> gGrad(size_t(m.nG) * kGeneralRowEntries), cN(nC);
    std::vector<unsigned char> gFlag(m.nG);
    std::vector<int> cNe(nC), cCtrl(size_t(nC) * kContactRowEntries), varOf(m.nCtrl);
    std::vector<unsigned long long> cKey(nC);
    V.gVal.download(gVal.data(), gVal.size() * 8);
    V.gTval.download(gTval.data(), gTval.size() * 8);
    V.gGrad.download(gGrad.data(), gGrad.size() * sizeof(d3));
    V.gFlag.download(gFlag.data(), gFlag.size());
    m.controls.varOf.download(varOf.data(), varOf.size() * 4);
    if (nC) {
        dval.download(cVal.data(), cVal.size() * 8);
        dtval.download(cTval.data(), cTval.size() * 8);
        dK.download(cK.data(), cK.size() * 8);
        dnormal.download(cN.data(), cN.size() * sizeof(d3));
        dn.download(cNe.data(), cNe.size() * 4);
        dctrl.download(cCtrl.data(), cCtrl.size() * 4);
        dkey.download(cKey.data(), cKey.size() * 8);
    }
    const GeneralRowArrays rows = downloadGeneralRows(m);
    const std::vector<uint64_t> keys = generalRowKeys(m, rows);
    std::vector<DebugRow> out;
    const double sc = m.pr.hardK / m.pr.kappaMax;
    for (int r = 0; r < m.nG; ++r) {
        if (gFlag[r]) continue;
        bool touches = onlyVariable < 0;  // every general row, or those touching onlyVariable
        for (int e = 0; e < rows.n[r]; ++e) touches = touches || rows.var[size_t(kGeneralRowEntries) * r + e] == onlyVariable;
        if (!touches) continue;
        DebugRow d;
        d.family = rows.fam[r];
        d.key = keys[r];
        d.value = gVal[r];
        d.n = rows.n[r];
        for (int e = 0; e < d.n; ++e) {
            d.var[e] = rows.var[size_t(kGeneralRowEntries) * r + e];
            const d3 g = gGrad[size_t(kGeneralRowEntries) * r + e];
            d.g[e][0] = g.x;
            d.g[e][1] = g.y;
            d.g[e][2] = g.z;
        }
        out.push_back(d);
        if (d.family == FCurv && gTval[r] > 0) {
            DebugRow t = d;
            t.family = FCurvT;
            t.key = (uint64_t(FCurvT) << 60) ^ (d.key & ((uint64_t(1) << 60) - 1));
            t.value = gTval[r];
            for (int e = 0; e < t.n; ++e)
                for (int k = 0; k < 3; ++k) t.g[e][k] *= sc;
            out.push_back(t);
        }
    }
    for (int c = 0; c < nC; ++c) {
        DebugRow d;
        d.family = FContact;
        d.key = cKey[c];
        d.value = cVal[c];
        d.n = cNe[c];
        d.hasGradient = cN[c].x != 0 || cN[c].y != 0 || cN[c].z != 0;
        for (int e = 0; e < d.n; ++e) {
            d.var[e] = varOf[cCtrl[size_t(kContactRowEntries) * c + e]];
            const double k = cK[size_t(kContactRowEntries) * c + e];
            d.g[e][0] = cN[c].x * k;
            d.g[e][1] = cN[c].y * k;
            d.g[e][2] = cN[c].z * k;
        }
        out.push_back(d);
        if (cTval[c] > 0) {
            DebugRow t = d;
            t.family = FContactT;
            t.key = (uint64_t(FContactT) << 60) ^ (d.key & ((uint64_t(1) << 60) - 1));
            t.value = cTval[c];
            out.push_back(t);
        }
    }
    return out;
}

std::vector<NcpDevice::FrozenPass> NcpDevice::frozenLinearizations(double rho, const std::vector<double>& p) {
    State& m = *state_;
    if (p.size() != size_t(3 * m.N)) throw DeviceError("ncp-gpu frozen: the product vector is not 3 per variable");
    std::vector<FrozenPass> out;
    for (int pass = 0; pass < 2; ++pass) {
        if (pass == 1) m.enqueueMultipliers();
        m.enqueueSamplePoints(RowSet::Current);
        refresh();
        FrozenPass f;
        f.eval = linearize(ControlSet::Current, rho);
        m.enqueueAssemble();
        const double zero = 0;
        check(cudaMemcpyAsync(reinterpret_cast<char*>(m.control()) + offsetof(NcpControl, mu), &zero, 8, cudaMemcpyHostToDevice, m.stream), "ncp-gpu frozen");
        check(cudaMemcpyAsync(reinterpret_cast<char*>(m.control()) + offsetof(NcpControl, eps), &zero, 8, cudaMemcpyHostToDevice, m.stream), "ncp-gpu frozen");
        check(cudaMemcpyAsync(m.system.p.data(), p.data(), p.size() * 8, cudaMemcpyHostToDevice, m.stream), "ncp-gpu frozen");
        m.enqueueOperatorProduct();
        check(cudaStreamSynchronize(m.stream), "ncp-gpu frozen");
        f.grad.resize(size_t(3 * m.N));
        f.hblock.resize(size_t(6 * m.N));
        f.product.resize(size_t(3 * m.N));
        m.system.grad.download(f.grad.data(), f.grad.size() * 8);
        m.system.hblock.download(f.hblock.data(), f.hblock.size() * 8);
        m.system.q.download(f.product.data(), f.product.size() * 8);
        NcpCounts c;
        m.countsBuffer.download(&c, sizeof c);
        f.contactRows = c.contactRows;
        f.activeRows = c.activeRows;
        f.activeEntries = c.activeEntries;
        out.push_back(std::move(f));
    }
    return out;
}

} // namespace knit::gpu
