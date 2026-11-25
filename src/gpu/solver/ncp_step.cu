// ncp-gpu: assembly (gradient and diagonal blocks), the step's reductions, the trial, the swept test,
// pinning, the regained merit, the KKT measures and the coverage audit of the tests (ncp.h). The
// damping, the step length and every other scalar come from the device control record.
#include "gpu/solver/ncp_device.cuh"

namespace knit::gpu {

namespace {

constexpr int kSweptMaxSteps = 256;  // conservative advancement steps of one pair
constexpr int kSumStride = kSumBlocks * kSumThreads;

// A stretched edge's Hessian beyond Gauss-Newton: kk (I - u u^T) between its controls (0 otherwise).
__device__ __forceinline__ double edgeStretchStiffness(const NcpRows& R, const NcpRowValues& V, const NcpParams& pr, int r, d3& u) {
    const double val = V.gVal[r], l = V.gAux[r];
    const d3 g0 = V.gGrad[kGeneralRowEntries * r];
    const double gl = norm(g0);
    if (!(val > 0 && l > 1e-12 && gl > 0)) return 0;
    u = g0 / gl;
    return pr.wL * val / (l * pr.D);
}

struct AssembledOutput {
    double* grad;   // 3 per variable
    double* gnorm;  // per variable
    double* hdiag;  // 3 per variable: H's diagonal
    double* hblock; // 6 per variable: H's diagonal block xx xy xz yy yz zz
};

// One thread per variable over its sorted active entries: the merit gradient and the 3x3 diagonal
// block of H = J^T W J (+ the stretched edges' Hessian). gradMax and diagMax by integer atomics.
__global__ void assembleKernel(NcpView v, RowValuePair pair, ActiveLists L, NcpControl* ctl, AssembledOutput out) {
    const NcpRows& R = v.rows;
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const int a = threadItem();
    if (a >= v.variables) return;
    double h[6] = {0, 0, 0, 0, 0, 0};
    d3 ga{0, 0, 0};
    for (int q = L.start[a]; q < L.start[a + 1]; ++q) {
        const int k = L.entry[q] >> 4, e = L.entry[q] & 15;
        const int r = L.rows[k];
        const double w = rowWgn(V, R.nG, r), wr = rowWres(V, R.nG, r);
        const d3 g = entryGradient(v, V, L, k, e, a);
        ga += g * wr;
        h[0] += w * g.x * g.x; h[1] += w * g.x * g.y; h[2] += w * g.x * g.z;
        h[3] += w * g.y * g.y; h[4] += w * g.y * g.z; h[5] += w * g.z * g.z;
        if (r < R.nG && R.gFam[r] == FEdge) {
            d3 u;
            const double kk = edgeStretchStiffness(R, V, v.params, r, u);
            h[0] += kk * (1 - u.x * u.x); h[1] -= kk * u.x * u.y; h[2] -= kk * u.x * u.z;
            h[3] += kk * (1 - u.y * u.y); h[4] -= kk * u.y * u.z; h[5] += kk * (1 - u.z * u.z);
        }
    }
    if (v.pin[a]) ga = {0, 0, 0};  // pinned variables carry no gradient: the stationarity is over the others
    out.grad[3 * a] = ga.x;
    out.grad[3 * a + 1] = ga.y;
    out.grad[3 * a + 2] = ga.z;
    const double gn = sqrt(dot(ga, ga));
    out.gnorm[a] = gn;
    out.hdiag[3 * a] = h[0];
    out.hdiag[3 * a + 1] = h[3];
    out.hdiag[3 * a + 2] = h[5];
    for (int k = 0; k < 6; ++k) out.hblock[6 * a + k] = h[k];
    atomicMaxPositive(&ctl->gradMax, gn);
    atomicMaxPositive(&ctl->diagMax, fmax(h[0], fmax(h[3], h[5])));
}

// The same assembly with one warp per variable (LiveRows; 10 % of a live tube call as one thread per
// variable): lane l takes the entries l, l + 32, ... in
// order, then the lanes' sums are added in a fixed tree. Deterministic; against assembleKernel's single
// sequential sum the grouping differs (rounding only).
constexpr int kWarp = 32;
__device__ __forceinline__ double warpSum(double x) {
    for (int offset = kWarp / 2; offset > 0; offset /= 2) x += __shfl_down_sync(0xffffffffu, x, offset);
    return x;
}

__global__ void assembleByWarpKernel(NcpView v, RowValuePair pair, ActiveLists L, NcpControl* ctl, AssembledOutput out) {
    const NcpRows& R = v.rows;
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const int lane = threadIdx.x % kWarp;
    const long long aa = (static_cast<long long>(blockIdx.x) * blockDim.x + threadIdx.x) / kWarp;
    if (aa >= v.variables) return;  // the whole warp: one variable per warp
    const int a = static_cast<int>(aa);
    double h[6] = {0, 0, 0, 0, 0, 0};
    d3 ga{0, 0, 0};
    for (int q = L.start[a] + lane; q < L.start[a + 1]; q += kWarp) {
        const int k = L.entry[q] >> 4, e = L.entry[q] & 15;
        const int r = L.rows[k];
        const double w = rowWgn(V, R.nG, r), wr = rowWres(V, R.nG, r);
        const d3 g = entryGradient(v, V, L, k, e, a);
        ga += g * wr;
        h[0] += w * g.x * g.x; h[1] += w * g.x * g.y; h[2] += w * g.x * g.z;
        h[3] += w * g.y * g.y; h[4] += w * g.y * g.z; h[5] += w * g.z * g.z;
        if (r < R.nG && R.gFam[r] == FEdge) {
            d3 u;
            const double kk = edgeStretchStiffness(R, V, v.params, r, u);
            h[0] += kk * (1 - u.x * u.x); h[1] -= kk * u.x * u.y; h[2] -= kk * u.x * u.z;
            h[3] += kk * (1 - u.y * u.y); h[4] -= kk * u.y * u.z; h[5] += kk * (1 - u.z * u.z);
        }
    }
    ga = d3{warpSum(ga.x), warpSum(ga.y), warpSum(ga.z)};
    for (int k = 0; k < 6; ++k) h[k] = warpSum(h[k]);
    if (lane != 0) return;
    if (v.pin[a]) ga = {0, 0, 0};
    out.grad[3 * a] = ga.x;
    out.grad[3 * a + 1] = ga.y;
    out.grad[3 * a + 2] = ga.z;
    const double gn = sqrt(dot(ga, ga));
    out.gnorm[a] = gn;
    out.hdiag[3 * a] = h[0];
    out.hdiag[3 * a + 1] = h[3];
    out.hdiag[3 * a + 2] = h[5];
    for (int k = 0; k < 6; ++k) out.hblock[6 * a + k] = h[k];
    atomicMaxPositive(&ctl->gradMax, gn);
    atomicMaxPositive(&ctl->diagMax, fmax(h[0], fmax(h[3], h[5])));
}

__device__ void sumBlocksInto(const double* partial, double* out) {
    __shared__ double sh[kSumBlocks];
    sh[threadIdx.x] = partial[threadIdx.x];
    __syncthreads();
    for (int k = kSumBlocks / 2; k > 0; k >>= 1) {
        if (static_cast<int>(threadIdx.x) < k) sh[threadIdx.x] += sh[threadIdx.x + k];
        __syncthreads();
    }
    if (threadIdx.x == 0) *out = sh[0];
    __syncthreads();
}

// The step's reductions from the linear solution x: g . dx and dx . (mu D + eps I) dx (block partials),
// the coverage step fraction (largest alpha keeping every control within delta of Xref) and max |dx|.
__global__ void stepKernel(NcpView v, const double* x, const double* grad, const double* hdiag, const d3* P, const d3* Xref, double delta, NcpControl* ctl,
                           double* partial) {
    __shared__ double sg[kSumThreads], sd[kSumThreads];
    const double mu = ctl->mu, eps = ctl->eps;
    double gdx = 0, dDd = 0;
    for (long long aa = static_cast<long long>(blockIdx.x) * kSumThreads + threadIdx.x; aa < v.variables; aa += kSumStride) {
        const int a = static_cast<int>(aa);
        d3 dx{x[3 * a], x[3 * a + 1], x[3 * a + 2]};
        if (v.pin[a]) dx = {0, 0, 0};
        gdx += grad[3 * a] * dx.x + grad[3 * a + 1] * dx.y + grad[3 * a + 2] * dx.z;
        dDd += (mu * hdiag[3 * a] + eps) * dx.x * dx.x + (mu * hdiag[3 * a + 1] + eps) * dx.y * dx.y + (mu * hdiag[3 * a + 2] + eps) * dx.z * dx.z;
        const double len = norm(dx);
        atomicMaxPositive(&ctl->dxMax, len);
        if (len > 0) {
            const int g = v.ctrlOfVar[a];
            atomicMinPositive(&ctl->alphaCov, (delta - norm(P[g] - Xref[g])) / len);
        }
    }
    sg[threadIdx.x] = gdx;
    sd[threadIdx.x] = dDd;
    __syncthreads();
    for (int k = kSumThreads / 2; k > 0; k >>= 1) {
        if (static_cast<int>(threadIdx.x) < k) {
            sg[threadIdx.x] += sg[threadIdx.x + k];
            sd[threadIdx.x] += sd[threadIdx.x + k];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        partial[blockIdx.x] = sg[0];
        partial[kSumBlocks + blockIdx.x] = sd[0];
    }
}

__global__ void stepSumsKernel(const double* partial, NcpControl* ctl) {
    sumBlocksInto(partial, &ctl->gdx);
    sumBlocksInto(partial + kSumBlocks, &ctl->dDd);
}

__global__ void stepResetKernel(NcpControl* ctl) {
    ctl->alphaCov = 1.0;
    ctl->dxMax = 0;
}

// Q = proj_cap(P + alpha dx), alpha = ctl->alpha.
__global__ void trialKernel(NcpView v, const double* x, const NcpControl* ctl, const d3* P, d3* Q) {
    const int a = threadItem();
    if (a >= v.variables) return;
    const int g = v.ctrlOfVar[a];
    if (v.pin[a]) {
        Q[g] = P[g];
        return;
    }
    const d3 dx{x[3 * a], x[3 * a + 1], x[3 * a + 2]};
    d3 q = P[g] + dx * ctl->alpha;
    if (v.params.cap > 0) {  // exact projection onto the displacement ball (non-expansive: coverage holds)
        const d3 d = q - v.P0[g];
        const double l = norm(d);
        if (l > v.params.cap) q = v.P0[g] + d * (v.params.cap / l);
    }
    Q[g] = q;
}

// Conservative advancement along the straight motion P -> Q of a candidate pair's two sample
// segments: the distance is L-Lipschitz in the motion parameter with L = max endpoint motion of
// segment a + that of segment b (a point of a segment moves by a convex combination of its end
// points' motions), so from alpha with distance d no contact closer than eps occurs before
// alpha + (d - eps) / L. The advancement aims at eps / 2 (each step gains at least eps / (2 L), so
// the iteration ends) while a distance <= eps counts as touching. Pairs closer than eps at the
// start are counted apart (not a new touch). counts: touching, unresolved, closer at the start.
__global__ void sweptKernel(NcpView v, const d3* P, const d3* Q, NcpControl* ctl) {
    const NcpRows& R = v.rows;
    const double eps = ctl->sweptEps;
    const int nC = *R.nRows;  // every row (a row without a term can still be touched)
    for (long long cc = loopStart(); cc < nC; cc += loopStride()) {
        const int c = static_cast<int>(cc);
        const int i = R.cPair[c].x, j = R.cPair[c].y;
        const d3 a0 = P[i], a1 = P[i + 1], b0 = P[j], b1 = P[j + 1];
        const d3 da0 = Q[i] - a0, da1 = Q[i + 1] - a1, db0 = Q[j] - b0, db1 = Q[j + 1] - b1;
        const double L = fmax(norm(da0), norm(da1)) + fmax(norm(db0), norm(db1));
        double s, t;
        double d = segmentDistance(a0, a1, b0, b1, s, t);
        if (d <= eps) {
            atomicAdd(&ctl->swept[2], 1ull);
            continue;
        }
        if (!(L > 0) || d - L > eps) continue;  // cannot come within eps
        double al = 0;
        bool decided = false;
        for (int it = 0; it < kSweptMaxSteps && !decided; ++it) {
            al += (d - 0.5 * eps) / L;
            if (al >= 1) {
                decided = true;
                break;
            }
            d = segmentDistance(a0 + da0 * al, a1 + da1 * al, b0 + db0 * al, b1 + db1 * al, s, t);
            if (d <= eps) {
                atomicAdd(&ctl->swept[0], 1ull);
                decided = true;
            }
        }
        if (!decided) atomicAdd(&ctl->swept[1], 1ull);
    }
}

// Band rows whose control is Owned at the current state and not at the trial: pin the control.
__global__ void pinLostKernel(NcpRows R, RowValuePair pair, unsigned char* pin, NcpControl* ctl) {
    const int r = threadItem();
    if (r >= R.nG || R.gFam[r] != FBandU || R.gN[r] < 1) return;
    const unsigned char* flagCur = rowValues(pair, RowSet::Current).gFlag;
    const unsigned char* flagTrial = rowValues(pair, RowSet::Trial).gFlag;
    if (!flagTrial[r] || flagCur[r]) return;
    const int v = R.gVar[kGeneralRowEntries * r];
    if (!pin[v]) {
        pin[v] = 1;
        atomicAdd(&ctl->scratchCount, 1ull);
    }
}

// General rows off at the current state (not Owned / degenerate) and on at the trial (a binding
// regained): their trial merit, as block partials, so the acceptance compares both states on the same rows.
__global__ void regainedMeritKernel(int nG, RowValuePair pair, double* partial) {
    __shared__ double sh[kSumThreads];
    const NcpRowValues& cur = rowValues(pair, RowSet::Current);
    const NcpRowValues& tri = rowValues(pair, RowSet::Trial);
    double s = 0;
    for (long long r = static_cast<long long>(blockIdx.x) * kSumThreads + threadIdx.x; r < nG; r += kSumStride)
        if (cur.gFlag[r] && !tri.gFlag[r]) s += tri.merit[r];
    sh[threadIdx.x] = s;
    __syncthreads();
    for (int k = kSumThreads / 2; k > 0; k >>= 1) {
        if (static_cast<int>(threadIdx.x) < k) sh[threadIdx.x] += sh[threadIdx.x + k];
        __syncthreads();
    }
    if (threadIdx.x == 0) partial[blockIdx.x] = sh[0];
}

__global__ void regainedSumKernel(const double* partial, NcpControl* ctl) { sumBlocksInto(partial, &ctl->regained); }

// Complementarity with lambda+ = max(0, lambda + rho g): max over hard rows of min(lambda+, -g) (g < 0).
__global__ void complementarityKernel(NcpRows R, ContactValueScale scale, RowValuePair pair, NcpControl* ctl) {
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const double rho = ctl->rho;
    const long long nR = R.nG + static_cast<long long>(*R.nC);
    for (long long rr = loopStart(); rr < nR; rr += loopStride()) {
        const int r = static_cast<int>(rr);
        double g, lam;
        if (r < R.nG) {
            const int fam = R.gFam[r];
            if (fam == FEdge || fam == FAnchor || V.gFlag[r]) continue;
            g = V.gVal[r];
            lam = R.gLam[r];
        } else {
            g = contactHardValue(V, scale, r - R.nG);
            lam = R.cLam[r - R.nG];
        }
        const double lp = fmax(0.0, lam + rho * g);
        if (lp > 0) atomicAdd(&ctl->activeMultipliers, 1ull);
        if (g < 0) atomicMaxPositive(&ctl->complementarity, fmin(lp, -g));
    }
}

// Tests (NcpOptions::auditCoverage): every pair of active sample segments closer than thr at the
// trial Q or at the midpoint of P -> Q must have been within the pad at Xref (a candidate with rows).
// Brute force over all pairs: small fixtures only.
__global__ void coverageAuditKernel(NcpView v, const d3* P, const d3* Q, const d3* Xref, double pad, double window, NcpControl* ctl) {
    const NcpSamples& S = v.samples;
    const int i = threadItem();
    if (i >= S.n || !S.segOk[i] || !S.active[i]) return;
    const double thr = v.params.thr;
    unsigned long long checks = 0, missing = 0;
    for (int j = 0; j < S.n; ++j) {
        if (j == i || !S.segOk[j]) continue;
        if (S.yarn[j] == S.yarn[i]) {
            const double gap = j > i ? S.s[j] - S.s[i + 1] : S.s[i] - S.s[j + 1];
            if (gap < window) continue;  // the same-yarn arc window: never a candidate
        }
        double s, t;
        const d3 qa0 = evalSample(Q, S, i), qa1 = evalSample(Q, S, i + 1), qb0 = evalSample(Q, S, j), qb1 = evalSample(Q, S, j + 1);
        const d3 pa0 = evalSample(P, S, i), pa1 = evalSample(P, S, i + 1), pb0 = evalSample(P, S, j), pb1 = evalSample(P, S, j + 1);
        const double dq = segmentDistance(qa0, qa1, qb0, qb1, s, t);
        const double dm = segmentDistance((pa0 + qa0) * 0.5, (pa1 + qa1) * 0.5, (pb0 + qb0) * 0.5, (pb1 + qb1) * 0.5, s, t);
        if (dq >= thr && dm >= thr) continue;
        ++checks;
        const double dx = segmentDistance(evalSample(Xref, S, i), evalSample(Xref, S, i + 1), evalSample(Xref, S, j), evalSample(Xref, S, j + 1), s, t);
        if (dx > pad) ++missing;
    }
    if (checks) atomicAdd(&ctl->coverageChecks, checks);
    if (missing) atomicAdd(&ctl->coverageViolations, missing);
}

}  // namespace

void NcpDevice::State::enqueueAssemble() {
    enqueueActiveLists(RowSet::Current);
    NcpControl* ctl = control();
    check(cudaMemsetAsync(&ctl->gradMax, 0, sizeof(double), stream), "ncp-gpu assemble");
    check(cudaMemsetAsync(&ctl->diagMax, 0, sizeof(double), stream), "ncp-gpu assemble");
    const AssembledOutput out{system.grad.doubles(), system.gnorm.doubles(), system.hdiag.doubles(), system.hblock.doubles()};
    if (N && storage == ContactStorage::LiveRows) {
        const long long threads = static_cast<long long>(N) * kWarp;
        assembleByWarpKernel<<<static_cast<unsigned>((threads + kThreadsPerBlock - 1) / kThreadsPerBlock), kThreadsPerBlock, 0, stream>>>(view(), valuePair(), activeLists(),
                                                                                                                                       ctl, out);
    } else if (N) {
        assembleKernel<<<blocksFor(N), kThreadsPerBlock, 0, stream>>>(view(), valuePair(), activeLists(), ctl, out);
    }
    checkLaunch("ncp-gpu: assemble");
}

void NcpDevice::State::enqueueStepSummary() {
    double* partial = system.partials.doubles();
    stepResetKernel<<<1, 1, 0, stream>>>(control());
    stepKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(view(), system.x.doubles(), system.grad.doubles(), system.hdiag.doubles(), controls.current(),
                                                       controls.reference(), delta, control(), partial);
    stepSumsKernel<<<1, kSumBlocks, 0, stream>>>(partial, control());
    checkLaunch("ncp-gpu: step");
}

void NcpDevice::State::enqueueTrial() {
    if (N) trialKernel<<<blocksFor(N), kThreadsPerBlock, 0, stream>>>(view(), system.x.doubles(), control(), controls.current(), controls.trial());
    enqueueSamplePoints(RowSet::Trial);
}

void NcpDevice::State::enqueueSwept() {
    check(cudaMemsetAsync(control()->swept, 0, sizeof(control()->swept), stream), "ncp-gpu swept");
    sweptKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(view(), pointsAtP(), pointsAtQ(), control());
}

void NcpDevice::State::enqueueCoverageAudit() {
    const double pad = pr.thr + 2 * delta;
    if (nSamples)
        coverageAuditKernel<<<blocksFor(nSamples), kThreadsPerBlock, 0, stream>>>(view(), controls.current(), controls.trial(), controls.reference(), pad, window, control());
}

void NcpDevice::State::enqueuePinLost() {
    check(cudaMemsetAsync(&control()->scratchCount, 0, sizeof(unsigned long long), stream), "ncp-gpu pin");
    if (hasBind && nG) pinLostKernel<<<blocksFor(nG), kThreadsPerBlock, 0, stream>>>(rows(), valuePair(), controls.pin.chars(), control());
}

void NcpDevice::State::enqueueRegained() {
    double* partial = system.partials.doubles();
    regainedMeritKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(nG, valuePair(), partial);
    regainedSumKernel<<<1, kSumBlocks, 0, stream>>>(partial, control());
}

// KKT measures at the current controls, BEFORE the multiplier update: the merit gradient at
// (lambda, rho) equals grad f + sum lambda+ grad g with lambda+ = max(0, lambda + rho g), so
// stationarity = max per-control |that gradient| (control()->gradMax); complementarity = max over hard
// rows with g < 0 of min(lambda+, -g); activeMultipliers = #(lambda+ > 0).
void NcpDevice::State::enqueueKkt() {
    enqueueTerms();
    enqueueAssemble();
    NcpControl* ctl = control();
    check(cudaMemsetAsync(&ctl->complementarity, 0, sizeof(double), stream), "ncp-gpu kkt");
    check(cudaMemsetAsync(&ctl->activeMultipliers, 0, sizeof(unsigned long long), stream), "ncp-gpu kkt");
    complementarityKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(rows(), contactValueScale(pr), valuePair(), ctl);
    checkLaunch("ncp-gpu: kkt");
}

}  // namespace knit::gpu
