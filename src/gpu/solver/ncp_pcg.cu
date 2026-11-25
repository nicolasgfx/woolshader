// ncp-gpu: the damped Gauss-Newton system (H + mu diag(H) + eps I) dx = -grad solved by conjugate
// gradients with a 3x3 block-Jacobi preconditioner, matrix-free and decided on the device (ncp.h).
//
// The product H p sums per row and then per variable: pass 1 writes t_r = w_r (J_r . p) for every
// active row (and a stretched edge's Hessian term kk (I - u u^T) (p_a - p_b)); pass 2 gathers, for each
// variable, g_{r,a} t_r over its sorted active entries (ActiveLists), so every sum has a fixed order.
// No matrix is assembled: the diagonal blocks of H (assembleKernel) feed the damping and the
// preconditioner. Pinned variables are held at zero (their right-hand side, preconditioned residual
// and product vanish), which solves the system of the others with them fixed.
//
// The scalars (rz, p.q, alpha, beta, r.r) stay in NcpControl; the stop rule |r| <= tol |b| or the
// iteration cap is tested on the device after every iteration and drives the graph's WHILE node.
#include "gpu/solver/ncp_device.cuh"

namespace knit::gpu {

namespace {

constexpr int kSumStride = kSumBlocks * kSumThreads;

__device__ __forceinline__ d3 vec3At(const double* x, int a) { return {x[3 * a], x[3 * a + 1], x[3 * a + 2]}; }

// Block sums of up to two per-thread values: partial[block], partial[kSumBlocks + block].
__device__ void blockPartials(double s0, double s1, double* partial) {
    __shared__ double a[kSumThreads], b[kSumThreads];
    a[threadIdx.x] = s0;
    b[threadIdx.x] = s1;
    __syncthreads();
    for (int k = kSumThreads / 2; k > 0; k >>= 1) {
        if (static_cast<int>(threadIdx.x) < k) {
            a[threadIdx.x] += a[threadIdx.x + k];
            b[threadIdx.x] += b[threadIdx.x + k];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        partial[blockIdx.x] = a[0];
        partial[kSumBlocks + blockIdx.x] = b[0];
    }
}

// The sum of kSumBlocks partials (one block of kSumBlocks threads); every thread gets it.
__device__ double sumOfPartials(const double* partial) {
    __shared__ double sh[kSumBlocks];
    sh[threadIdx.x] = partial[threadIdx.x];
    __syncthreads();
    for (int k = kSumBlocks / 2; k > 0; k >>= 1) {
        if (static_cast<int>(threadIdx.x) < k) sh[threadIdx.x] += sh[threadIdx.x + k];
        __syncthreads();
    }
    const double s = sh[0];
    __syncthreads();
    return s;
}

// z = M r with the block inverse (6 per variable); pinned: z = 0.
__device__ __forceinline__ d3 precondition(const double* inv, const unsigned char* pin, int a, d3 r) {
    if (pin[a]) return {0, 0, 0};
    const double* o = inv + 6 * a;
    return {o[0] * r.x + o[1] * r.y + o[2] * r.z, o[1] * r.x + o[3] * r.y + o[4] * r.z, o[2] * r.x + o[4] * r.y + o[5] * r.z};
}

// The damped diagonal blocks' inverses (Cholesky L L^T, then L^-T L^-1); a block that is not positive
// definite falls back to plain Jacobi on its positive diagonal entries (counted).
__global__ void preconditionerKernel(int N, const double* hblock, NcpControl* ctl, double* inv) {
    const int a = threadItem();
    if (a >= N) return;
    const double mu = ctl->mu, eps = ctl->eps;
    const double* h = hblock + 6 * a;
    const double xx = h[0] * (1 + mu) + eps, xy = h[1], xz = h[2], yy = h[3] * (1 + mu) + eps, yz = h[4], zz = h[5] * (1 + mu) + eps;
    const double l00 = xx > 0 ? sqrt(xx) : 0;
    const double l10 = l00 > 0 ? xy / l00 : 0, l20 = l00 > 0 ? xz / l00 : 0;
    const double d1 = yy - l10 * l10;
    const double l11 = d1 > 0 ? sqrt(d1) : 0;
    const double l21 = l11 > 0 ? (yz - l20 * l10) / l11 : 0;
    const double d2 = zz - l20 * l20 - l21 * l21;
    const double l22 = d2 > 0 ? sqrt(d2) : 0;
    double* o = inv + 6 * a;
    if (!(l00 > 0 && l11 > 0 && l22 > 0)) {
        atomicAdd(reinterpret_cast<unsigned long long*>(&ctl->nonSpdBlocks), 1ull);
        o[0] = xx > 0 ? 1 / xx : 1; o[1] = 0; o[2] = 0;
        o[3] = yy > 0 ? 1 / yy : 1; o[4] = 0;
        o[5] = zz > 0 ? 1 / zz : 1;
        return;
    }
    const double i00 = 1 / l00, i11 = 1 / l11, i22 = 1 / l22;
    const double i10 = -l10 * i00 * i11;
    const double i21 = -l21 * i11 * i22;
    const double i20 = -(l20 * i00 + l21 * i10) * i22;
    o[0] = i00 * i00 + i10 * i10 + i20 * i20;
    o[1] = i10 * i11 + i20 * i21;
    o[2] = i20 * i22;
    o[3] = i11 * i11 + i21 * i21;
    o[4] = i21 * i22;
    o[5] = i22 * i22;
}

struct PcgVectors {
    double* x;
    double* r;
    double* z;
    double* p;
    double* q;
    const double* inv;
};

// x = 0, r = b = -grad (pinned: 0), z = M r, p = z; partials r.z and b.b.
__global__ void pcgStartVectorsKernel(int N, const double* grad, const unsigned char* pin, PcgVectors w, double* partial) {
    double rz = 0, bb = 0;
    for (long long aa = static_cast<long long>(blockIdx.x) * kSumThreads + threadIdx.x; aa < N; aa += kSumStride) {
        const int a = static_cast<int>(aa);
        const d3 b = pin[a] ? d3{0, 0, 0} : vec3At(grad, a) * -1.0;
        const d3 z = precondition(w.inv, pin, a, b);
        for (int k = 0; k < 3; ++k) {
            w.x[3 * a + k] = 0;
            w.r[3 * a + k] = comp(b, k);
            w.z[3 * a + k] = comp(z, k);
            w.p[3 * a + k] = comp(z, k);
        }
        rz += dot(b, z);
        bb += dot(b, b);
    }
    blockPartials(rz, bb, partial);
}

__global__ void pcgStartScalarsKernel(const double* partial, NcpControl* ctl, cudaGraphConditionalHandle loop) {
    const double rz = sumOfPartials(partial), bb = sumOfPartials(partial + kSumBlocks);
    if (threadIdx.x != 0) return;
    ctl->pcgRz = rz;
    ctl->pcgBb = bb;
    ctl->pcgRr = bb;
    ctl->pcgIteration = 0;
    ctl->pcgBreakdown = 0;
    cudaGraphSetConditional(loop, bb > 0 && ctl->pcgMaxIterations > 0 ? 1 : 0);
}

// J_r . p of active slot k: a general row over its entries; a contact row over its controls (in entry
// order), each with its coefficient from the slot's closest points.
__device__ __forceinline__ double rowDotP(const NcpView& v, const NcpRowValues& V, const ActiveLists& L, int k, const double* p) {
    const NcpRows& R = v.rows;
    const int r = L.rows[k];
    double jp = 0;
    if (r < R.nG) {
        for (int e = 0; e < R.gN[r]; ++e) jp += dot(V.gGrad[kGeneralRowEntries * r + e], vec3At(p, R.gVar[kGeneralRowEntries * r + e]));
        return jp;
    }
    const int2 pair = R.cPair[contactRowOf(R, r - R.nG)];
    int controls[kContactRowEntries];
    bool overflow = false;
    const int n = contactControls(v.samples, v.varOf, pair.x, pair.y, controls, overflow);
    const d3 nrm = L.normal[k];
    const double s = L.st[2 * k], t = L.st[2 * k + 1];
    for (int e = 0; e < n; ++e) {
        const d3 g = nrm * contactCoefficient(v.samples, pair.x, pair.y, s, t, controls[e], v.params.D);
        jp += dot(g, vec3At(p, v.varOf[controls[e]]));
    }
    return jp;
}

// Pass 1 of the product: t_k = w_r (J_r . p) per active slot; a stretched edge also kk (I - u u^T) (p_0 - p_1).
__global__ void rowProductKernel(NcpView v, RowValuePair pair, ActiveLists L, const NcpCounts* counts, const double* p, double* rowT, d3* rowS) {
    const NcpRows& R = v.rows;
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const int nAct = counts->activeRows;
    for (long long kk = loopStart(); kk < nAct; kk += loopStride()) {
        const int k = static_cast<int>(kk);
        const int r = L.rows[k];
        const int n = r < R.nG ? R.gN[r] : 0;
        rowT[k] = rowWgn(V, R.nG, r) * rowDotP(v, V, L, k, p);
        if (r < R.nG && R.gFam[r] == FEdge) {
            const double val = V.gVal[r], l = V.gAux[r];
            const d3 g0 = V.gGrad[kGeneralRowEntries * r];
            const double gl = norm(g0);
            d3 s{0, 0, 0};
            if (val > 0 && l > 1e-12 && gl > 0) {
                const d3 u = g0 / gl;
                const double stiff = v.params.wL * val / (l * v.params.D);
                d3 dp = vec3At(p, R.gVar[kGeneralRowEntries * r]);
                if (n > 1) dp = dp - vec3At(p, R.gVar[kGeneralRowEntries * r + 1]);
                s = (dp - u * dot(u, dp)) * stiff;
            }
            rowS[r] = s;
        }
    }
}

// Pass 2: q_a = sum over a's active entries of g_{r,a} t_r (+- an edge's s_r) + (mu diag(H) + eps) p_a;
// pinned: 0. Partials of p.q.
__global__ void variableProductKernel(NcpView v, RowValuePair pair, ActiveLists L, const double* hdiag, const NcpControl* ctl, const double* rowT,
                                      const d3* rowS, PcgVectors w, double* partial) {
    const NcpRows& R = v.rows;
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const double mu = ctl->mu, eps = ctl->eps;
    double pq = 0;
    for (long long aa = static_cast<long long>(blockIdx.x) * kSumThreads + threadIdx.x; aa < v.variables; aa += kSumStride) {
        const int a = static_cast<int>(aa);
        d3 y{0, 0, 0};
        if (!v.pin[a]) {
            for (int q = L.start[a]; q < L.start[a + 1]; ++q) {
                const int k = L.entry[q] >> 4, e = L.entry[q] & 15;
                const int r = L.rows[k];
                y += entryGradient(v, V, L, k, e, a) * rowT[k];
                if (r < R.nG && R.gFam[r] == FEdge) y += rowS[r] * (e == 0 ? 1.0 : -1.0);
            }
            const d3 pa = vec3At(w.p, a);
            y += d3{(mu * hdiag[3 * a] + eps) * pa.x, (mu * hdiag[3 * a + 1] + eps) * pa.y, (mu * hdiag[3 * a + 2] + eps) * pa.z};
        }
        for (int k = 0; k < 3; ++k) w.q[3 * a + k] = comp(y, k);
        pq += dot(vec3At(w.p, a), y);
    }
    blockPartials(pq, 0, partial);
}

// Pass 2 with one warp per variable (LiveRows; the M1 profile measured the one-thread gather at 42 % of a
// tube call, and a lane-0 sum in entry order at 46 % of a live call):
// lane l sums the entries l, l + 32, l + 64, ... in order, then the 32 lane sums are added in a
// fixed tree (warpSum). Every run adds the same numbers in the same order (deterministic); against
// variableProductKernel's single sequential sum the grouping differs (rounding only).
constexpr int kWarp = 32;

// The sum of the 32 lanes' values in a fixed tree (offsets 16, 8, 4, 2, 1); lane 0 holds it.
__device__ __forceinline__ double warpSum(double x) {
    for (int offset = kWarp / 2; offset > 0; offset /= 2) x += __shfl_down_sync(0xffffffffu, x, offset);
    return x;
}

__global__ void variableProductByWarpKernel(NcpView v, RowValuePair pair, ActiveLists L, const double* hdiag, const NcpControl* ctl, const double* rowT,
                                            const d3* rowS, PcgVectors w) {
    const NcpRows& R = v.rows;
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const int lane = threadIdx.x % kWarp;
    const long long aa = (static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x) / kWarp;
    if (aa >= v.variables) return;  // the whole warp: one variable per warp
    const int a = static_cast<int>(aa);
    d3 y{0, 0, 0};
    if (!v.pin[a]) {
        const int end = L.start[a + 1];
        for (int q = L.start[a] + lane; q < end; q += kWarp) {
            const int k = L.entry[q] >> 4, e = L.entry[q] & 15;
            const int r = L.rows[k];
            y += entryGradient(v, V, L, k, e, a) * rowT[k];
            if (r < R.nG && R.gFam[r] == FEdge) y += rowS[r] * (e == 0 ? 1.0 : -1.0);
        }
    }
    y = d3{warpSum(y.x), warpSum(y.y), warpSum(y.z)};
    if (lane != 0) return;
    if (!v.pin[a]) {
        const double mu = ctl->mu, eps = ctl->eps;
        const d3 pa = vec3At(w.p, a);
        y += d3{(mu * hdiag[3 * a] + eps) * pa.x, (mu * hdiag[3 * a + 1] + eps) * pa.y, (mu * hdiag[3 * a + 2] + eps) * pa.z};
    }
    for (int k = 0; k < 3; ++k) w.q[3 * a + k] = comp(y, k);
}

// The partials of p.q over q, in variableProductKernel's order (thread t sums variables t, t + stride, ...).
__global__ void pqPartialsKernel(int variables, PcgVectors w, double* partial) {
    double pq = 0;
    for (long long aa = static_cast<long long>(blockIdx.x) * kSumThreads + threadIdx.x; aa < variables; aa += kSumStride) {
        const int a = static_cast<int>(aa);
        pq += dot(vec3At(w.p, a), vec3At(w.q, a));
    }
    blockPartials(pq, 0, partial);
}

// alpha = rz / p.q (0 and a breakdown when p.q <= 0: the iterate is kept).
__global__ void pcgAlphaKernel(const double* partial, NcpControl* ctl) {
    const double pq = sumOfPartials(partial);
    if (threadIdx.x != 0) return;
    ctl->pcgPq = pq;
    if (pq > 0) {
        ctl->pcgAlpha = ctl->pcgRz / pq;
    } else {
        ctl->pcgAlpha = 0;
        ctl->pcgBreakdown = 1;
    }
}

// x += alpha p, r -= alpha q, z = M r; partials r.r and r.z.
__global__ void pcgUpdateKernel(int N, const unsigned char* pin, const NcpControl* ctl, PcgVectors w, double* partial) {
    const double alpha = ctl->pcgAlpha;
    double rr = 0, rz = 0;
    for (long long aa = static_cast<long long>(blockIdx.x) * kSumThreads + threadIdx.x; aa < N; aa += kSumStride) {
        const int a = static_cast<int>(aa);
        d3 r = vec3At(w.r, a);
        r = r - vec3At(w.q, a) * alpha;
        const d3 z = precondition(w.inv, pin, a, r);
        for (int k = 0; k < 3; ++k) {
            w.x[3 * a + k] += alpha * w.p[3 * a + k];
            w.r[3 * a + k] = comp(r, k);
            w.z[3 * a + k] = comp(z, k);
        }
        rr += dot(r, r);
        rz += dot(r, z);
    }
    blockPartials(rr, rz, partial);
}

// The stop rule after an iteration, and beta = rz_new / rz for the next one.
__global__ void pcgDecideKernel(const double* partial, NcpControl* ctl, cudaGraphConditionalHandle loop) {
    const double rr = sumOfPartials(partial), rz = sumOfPartials(partial + kSumBlocks);
    if (threadIdx.x != 0) return;
    ++ctl->pcgIteration;
    ctl->pcgRr = rr;
    const double target = ctl->pcgTolerance * ctl->pcgTolerance * ctl->pcgBb;
    const bool done = ctl->pcgBreakdown || rr <= target || ctl->pcgIteration >= ctl->pcgMaxIterations;
    ctl->pcgBeta = ctl->pcgRz != 0 ? rz / ctl->pcgRz : 0;
    ctl->pcgRz = rz;
    cudaGraphSetConditional(loop, done ? 0 : 1);
}

// p = z + beta p.
__global__ void pcgDirectionKernel(int N, const NcpControl* ctl, PcgVectors w) {
    const long long i = threadItem64();
    if (i >= 3LL * N) return;
    w.p[i] = w.z[i] + ctl->pcgBeta * w.p[i];
}

// r . x (the model curvature of an inexact solve: dx.H.dx = -g.dx - r.dx - dDd) and the statistics.
__global__ void pcgResidualDotKernel(int N, PcgVectors w, double* partial) {
    double s = 0;
    for (long long aa = static_cast<long long>(blockIdx.x) * kSumThreads + threadIdx.x; aa < N; aa += kSumStride) {
        const int a = static_cast<int>(aa);
        s += dot(vec3At(w.r, a), vec3At(w.x, a));
    }
    blockPartials(s, 0, partial);
}

__global__ void pcgFinishKernel(const double* partial, NcpControl* ctl) {
    const double rx = sumOfPartials(partial);
    if (threadIdx.x != 0) return;
    ctl->rDotX = rx;
    ctl->pcgIterations += ctl->pcgIteration;
    ++ctl->solves;
    if (ctl->pcgBreakdown) ++ctl->pcgBreakdowns;
    if (ctl->pcgRr > ctl->pcgTolerance * ctl->pcgTolerance * ctl->pcgBb && ctl->pcgBb > 0) ++ctl->pcgCapped;
}

}  // namespace

void NcpDevice::State::enqueuePcgSetup(cudaGraphConditionalHandle loop) {
    LinearSystem& s = system;
    if (N) preconditionerKernel<<<blocksFor(N), kThreadsPerBlock, 0, stream>>>(N, s.hblock.doubles(), control(), s.precond.doubles());
    const PcgVectors w{s.x.doubles(), s.r.doubles(), s.z.doubles(), s.p.doubles(), s.q.doubles(), s.precond.doubles()};
    pcgStartVectorsKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(N, s.grad.doubles(), controls.pin.chars(), w, s.partials.doubles());
    pcgStartScalarsKernel<<<1, kSumBlocks, 0, stream>>>(s.partials.doubles(), control(), loop);
    checkLaunch("ncp-gpu: pcg setup");
}

// q := (H + mu diag(H) + eps I) p and the partials of p.q (pass 1, then pass 2 by storage: the warp
// gather for LiveRows, the thread gather for CompactRows; the same q to rounding).
void NcpDevice::State::enqueueProductPasses() {
    LinearSystem& s = system;
    const PcgVectors w{s.x.doubles(), s.r.doubles(), s.z.doubles(), s.p.doubles(), s.q.doubles(), s.precond.doubles()};
    rowProductKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(view(), valuePair(), activeLists(), counts(), s.p.doubles(), s.rowT.doubles(),
                                                                   static_cast<d3*>(s.rowS.data()));
    if (storage == ContactStorage::LiveRows) {
        const long long threads = static_cast<long long>(N) * kWarp;
        if (N)
            variableProductByWarpKernel<<<static_cast<unsigned>((threads + kThreadsPerBlock - 1) / kThreadsPerBlock), kThreadsPerBlock, 0, stream>>>(
                view(), valuePair(), activeLists(), s.hdiag.doubles(), control(), s.rowT.doubles(), static_cast<const d3*>(s.rowS.data()), w);
        pqPartialsKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(N, w, s.partials.doubles());
    } else {
        variableProductKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(view(), valuePair(), activeLists(), s.hdiag.doubles(), control(), s.rowT.doubles(),
                                                                      static_cast<const d3*>(s.rowS.data()), w, s.partials.doubles());
    }
}

void NcpDevice::State::enqueuePcgIteration(cudaGraphConditionalHandle loop) {
    LinearSystem& s = system;
    const PcgVectors w{s.x.doubles(), s.r.doubles(), s.z.doubles(), s.p.doubles(), s.q.doubles(), s.precond.doubles()};
    enqueueProductPasses();
    pcgAlphaKernel<<<1, kSumBlocks, 0, stream>>>(s.partials.doubles(), control());
    pcgUpdateKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(N, controls.pin.chars(), control(), w, s.partials.doubles());
    pcgDecideKernel<<<1, kSumBlocks, 0, stream>>>(s.partials.doubles(), control(), loop);
    if (N) pcgDirectionKernel<<<blocksFor(3LL * N), kThreadsPerBlock, 0, stream>>>(N, control(), w);
    checkLaunch("ncp-gpu: pcg iteration");
}

void NcpDevice::State::enqueuePcgFinish() {
    LinearSystem& s = system;
    const PcgVectors w{s.x.doubles(), s.r.doubles(), s.z.doubles(), s.p.doubles(), s.q.doubles(), s.precond.doubles()};
    pcgResidualDotKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(N, w, s.partials.doubles());
    pcgFinishKernel<<<1, kSumBlocks, 0, stream>>>(s.partials.doubles(), control());
    checkLaunch("ncp-gpu: pcg finish");
}

void NcpDevice::State::enqueueOperatorProduct() {
    enqueueProductPasses();
    checkLaunch("ncp-gpu: operator product");
}

}  // namespace knit::gpu
