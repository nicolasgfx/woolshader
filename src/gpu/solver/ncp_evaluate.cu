// ncp-gpu: row values, merit, multipliers (ncp.h), and the fixed-grid sums every part uses. Each
// sequence is enqueued on the call's stream; the penalty rho and the row set's index come from the
// device control record, so the captured kernels serve every iteration.
#include "gpu/solver/ncp_device.cuh"

namespace knit::gpu {

namespace {

constexpr int kSumStride = kSumBlocks * kSumThreads;

// What rowSumsPartialKernel sums: the general rows, and with nC the contact rows, whose terms come
// from their stored values at the multipliers and the penalty now (contactMeritTerms).
struct RowSums {
    int nG;
    const int* nC;  // null: the general rows only (a trial: its contact rows were reduced)
    const double* cLam;
    double wT;
    // LiveRows: the values from the stored distance (contactHardValue's expressions). Explicit fields, not
    // the whole NcpParams: with NcpParams in this struct, compact solves of ncp_gpu_scaling_cylinders took
    // another path (iteration cap at 24 turns; the summed values are the same; cause not found).
    double hardC, thr, D;
    bool targets;
};

// Block partial sums of the row set's merit, objective and edge strain^2 over the general rows and,
// with nC, the contact rows: partial[k * kSumBlocks + block], k = 0 merit, 1 objective, 2 strain^2.
// The rows are visited in one order (general, then contact), as when every row's terms were stored.
__global__ void rowSumsPartialKernel(RowValuePair pair, RowSet which, RowSums rs, const NcpControl* ctl, double* partial) {
    const NcpRowValues& V = rowValues(pair, which);
    __shared__ double sm[kSumThreads], so[kSumThreads], ss[kSumThreads];
    const int nG = rs.nG;
    const long long n = nG + (rs.nC ? static_cast<long long>(*rs.nC) : 0);
    const double rho = ctl->rho;
    double m = 0, o = 0, st = 0;
    for (long long i = static_cast<long long>(blockIdx.x) * kSumThreads + threadIdx.x; i < n; i += kSumStride) {
        if (i < nG) {
            m += V.merit[i];
            o += V.obj[i];
            st += V.strain2[i];
        } else {
            const long long c = i - nG;
            double cm, co;
            double val = 0, tval = 0;
            if (V.cDist) {  // LiveRows: from the stored distance, as evalContactRowsKernel computes them
                const double d = V.cDist[c];
                val = (rs.hardC - d) / rs.D;
                tval = rs.targets && d < rs.thr ? (rs.thr - d) / rs.D : 0.0;
            } else {
                val = V.cVal[c];
                tval = V.cTval[c];
            }
            contactMeritTerms(val, tval, rs.cLam[c], rho, rs.wT, cm, co);
            m += cm;
            o += co;
        }
    }
    sm[threadIdx.x] = m;
    so[threadIdx.x] = o;
    ss[threadIdx.x] = st;
    __syncthreads();
    for (int k = kSumThreads / 2; k > 0; k >>= 1) {
        if (static_cast<int>(threadIdx.x) < k) {
            sm[threadIdx.x] += sm[threadIdx.x + k];
            so[threadIdx.x] += so[threadIdx.x + k];
            ss[threadIdx.x] += ss[threadIdx.x + k];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        partial[blockIdx.x] = sm[0];
        partial[kSumBlocks + blockIdx.x] = so[0];
        partial[2 * kSumBlocks + blockIdx.x] = ss[0];
    }
}

// The sum of kSumBlocks partials into *out (one block of kSumBlocks threads, fixed order).
__device__ void sumBlocks(const double* partial, double* out) {
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

// The three sums of rowSumsPartialKernel into the given places (null: not needed).
__global__ void rowSumsFinalKernel(const double* partial, double* merit, double* objective, double* strain2) {
    if (merit) sumBlocks(partial, merit);
    if (objective) sumBlocks(partial + kSumBlocks, objective);
    if (strain2) sumBlocks(partial + 2 * kSumBlocks, strain2);
}

// A trial's contact merit and objective from reduceContactRowsKernel's partials.
__global__ void contactSumsFinalKernel(const double* blockSums, EvalReductions* sum) {
    sumBlocks(blockSums, &sum->contactMerit);
    sumBlocks(blockSums + kSumBlocks, &sum->contactObjective);
}

// The band's surface at each free control's position in X, once per evaluation (both band rows read it).
__global__ void bandFeetKernel(NcpView v, const d3* X, BandFeet feet) {
    const int a = threadItem();
    if (a >= v.variables) return;
    const int g = v.ctrlOfVar[a];
    d3 foot{0, 0, 0}, nrm{0, 0, 1};
    feet.ok[g] = bandSurface(v.band, v.bindState, g, X[g], foot, nrm) ? 1 : 0;
    feet.foot[g] = foot;
    feet.nrm[g] = nrm;
}

// The curvature |C' x C''| / |C'|^3 of span g (controls g .. g + 3) at parameter t; 1e300 where C' vanishes.
__device__ double spanCurvatureAt(const d3* P, int g, double t) {
    double w[4], b1[4], b2[4];
    bsplineBasis(t, w, b1, b2);
    d3 d1{0, 0, 0}, d2{0, 0, 0};
    for (int i = 0; i < 4; ++i) {
        d1 += P[g + i] * b1[i];
        d2 += P[g + i] * b2[i];
    }
    const double L = norm(d1);
    return L > 1e-150 ? norm(cross(d1, d2)) / (L * L * L) : 1e300;
}

// The parameter of span g's curvature maximum (CurvatureRows::SamplesAndSpanMaximum): the largest of 33 equally
// spaced values, refined by golden-section search between its two neighbours (24 steps shrink that bracket of
// 1/16 to 3e-7). The scan's own point stays when the refined one is lower (a maximum at a span end). A point
// where C' vanishes is the maximum, so the row reports the span as degenerate.
__device__ double spanMaximumParameter(const d3* P, int g) {
    constexpr int kScan = 32;
    constexpr int kGoldenSteps = 24;
    constexpr double kInvPhi = 0.6180339887498949;
    int best = 0;
    double kBest = -1;
    for (int j = 0; j <= kScan; ++j) {
        const double k = spanCurvatureAt(P, g, double(j) / kScan);
        if (k > kBest) {
            kBest = k;
            best = j;
        }
    }
    double a = fmax(0.0, double(best - 1) / kScan), b = fmin(1.0, double(best + 1) / kScan);
    double x1 = b - kInvPhi * (b - a), x2 = a + kInvPhi * (b - a);
    double f1 = spanCurvatureAt(P, g, x1), f2 = spanCurvatureAt(P, g, x2);
    for (int s = 0; s < kGoldenSteps; ++s) {
        if (f1 < f2) {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + kInvPhi * (b - a);
            f2 = spanCurvatureAt(P, g, x2);
        } else {
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = b - kInvPhi * (b - a);
            f1 = spanCurvatureAt(P, g, x1);
        }
    }
    const double tRefined = 0.5 * (a + b);
    return spanCurvatureAt(P, g, tRefined) >= kBest ? tRefined : double(best) / kScan;
}

__global__ void evalGeneralRowsKernel(NcpView v, const d3* P, RowValuePair pair, RowSet which, BandFeet feet, const NcpControl* ctl, EvalReductions* sum) {
    const NcpRows& R = v.rows;
    const NcpParams& pr = v.params;
    const NcpRowValues& V = rowValues(pair, which);
    const double rho = ctl->rho;
    const int r = threadItem();
    if (r >= R.nG) return;
    const int fam = R.gFam[r], g = R.gIdx[r], n = R.gN[r];
    d3* G = V.gGrad + kGeneralRowEntries * r;
    for (int e = 0; e < kGeneralRowEntries; ++e) G[e] = {0, 0, 0};
    double val = 0, tval = 0, aux = 0, wres = 0, wgn = 0, m = 0, o = 0, st2 = 0;
    unsigned char flag = 0;
    if (fam == FEdge) {
        const d3 d = P[g + 1] - P[g];
        const double l = norm(d);
        val = (l - R.gRest[r]) / pr.D;
        aux = l;
        if (l > 1e-12) {
            // Entries: control g + 1 (+u) first, then g (-u) (ncp_reference.cpp linearize), free ones
            // only; gSub bit e set: entry e is control g.
            const d3 u = d / (l * pr.D);
            for (int e = 0; e < n; ++e) G[e] = (R.gSub[r] >> e & 1) ? u * -1.0 : u;
        }
        wres = pr.wL * val;
        wgn = pr.wL;
        m = o = 0.5 * pr.wL * val * val;
        if (R.gRest[r] > 0) {
            const double s = fabs(l - R.gRest[r]) / R.gRest[r];
            st2 = s * s;
            atomicMaxPositive(&sum->strainMax, s);
            atomicAdd(&sum->edges, 1ull);
        }
        atomicMaxPositive(&sum->maxViol[FEdge], fabs(val));
    } else if (fam == FCurv) {
        double w[4], b1[4], b2[4];
        // q < curvQ: a fixed sample; q = curvQ: the span's curvature maximum at this state. Where that maximum is
        // unique, its gradient is the curvature's gradient at the maximizing parameter (Danskin's theorem), the
        // expression below.
        const int q = R.gSub[r] & 0xff;
        bsplineBasis(q == pr.curvQ ? spanMaximumParameter(P, g) : double(q) / (pr.curvQ - 1), w, b1, b2);
        d3 d1{0, 0, 0}, d2{0, 0, 0};
        for (int i = 0; i < 4; ++i) {
            d1 += P[g + i] * b1[i];
            d2 += P[g + i] * b2[i];
        }
        const double L = norm(d1);
        if (L < pr.kDeg) {
            flag = 1;
            atomicAdd(&sum->degenerate, 1ull);
        } else {
            const d3 c = cross(d1, d2);
            const double C = norm(c), L3 = L * L * L;
            const double kappa = C / L3;
            d3 g1 = d1 * (-3 * C / (L3 * L * L)), g2{0, 0, 0};
            if (C > 1e-15) {
                g1 += cross(d2, c) / (C * L3);
                g2 = cross(c, d1) / (C * L3);
            }
            val = kappa / pr.hardK - 1;
            // entries: the free controls among g .. g + 3 in order (gSub bits 8 + 2e: their offsets)
            for (int e = 0; e < n; ++e) {
                const int i = (R.gSub[r] >> (8 + 2 * e)) & 3;
                G[e] = (g1 * b1[i] + g2 * b2[i]) / pr.hardK;
            }
            hardTerm(val, R.gLam[r], rho, wres, wgn, m);
            const double vv = val > 0 ? val : 0;
            atomicMaxPositive(&sum->maxViol[FCurv], vv);
            if (val > kViolationTolerance) atomicAdd(&sum->nViol[FCurv], 1ull);
            const double sc = pr.hardK / pr.kappaMax;  // the target row's gradient is sc * the hard row's
            if (pr.targets && kappa > pr.kappaMax) {
                tval = kappa / pr.kappaMax - 1;
                wres += pr.wT * tval * sc;
                m += 0.5 * pr.wT * tval * tval;
                o += 0.5 * pr.wT * tval * tval;
                atomicMaxPositive(&sum->maxViol[FCurvT], tval);
                if (tval > kViolationTolerance) atomicAdd(&sum->nViol[FCurvT], 1ull);
            }
            if (pr.targets && kappa / pr.kappaMax - 1 > 0) wgn += pr.wT * sc * sc;
        }
    } else if (fam == FBandU || fam == FBandL) {
        const bool owned = feet.ok[g] != 0;
        const d3 foot = feet.foot[g], nrm = feet.nrm[g];
        if (n > 0 && v.pin[R.gVar[kGeneralRowEntries * r]]) {
            // A pinned control (pinLostOwnership) is held like a fixed one: its own band row is off
            // (its value is a constant of the call); counted and its violation reported.
            flag = 1;
            if (owned) {
                const d3 dv = P[g] - foot;
                const double z = (dot(dv, nrm) >= 0 ? 1.0 : -1.0) * norm(dv);
                atomicMaxPositive(&sum->pinnedBandViol, (fam == FBandU ? z : -z) / pr.bandHalf - 1);
            }
            if (fam == FBandU) atomicAdd(&sum->pinnedBand, 1ull);
        } else if (!owned) {
            flag = 1;
            if (fam == FBandU) atomicAdd(&sum->notOwned, 1ull);
        } else {
            // Signed distance to the foot (sign of the surface normal there) and its exact gradient.
            const d3 dv = P[g] - foot;
            const double l = norm(dv), s = dot(dv, nrm) >= 0 ? 1.0 : -1.0;
            const double z = s * l;
            const d3 u = l > 1e-12 ? dv * (s / l) : nrm;
            const double sgn = fam == FBandU ? 1.0 : -1.0;
            val = sgn * z / pr.bandHalf - 1;
            G[0] = u * (sgn / pr.bandHalf);
            hardTerm(val, R.gLam[r], rho, wres, wgn, m);
            atomicMaxPositive(&sum->maxViol[fam], val);
            if (val > kViolationTolerance) atomicAdd(&sum->nViol[fam], 1ull);
        }
    } else if (fam == FCap) {
        const d3 d = P[g] - v.P0[g];
        const double cap2 = pr.cap * pr.cap;
        val = dot(d, d) / cap2 - 1;
        G[0] = d * (2.0 / cap2);
        hardTerm(val, R.gLam[r], rho, wres, wgn, m);
        atomicMaxPositive(&sum->maxViol[FCap], val);
        if (val > kViolationTolerance) atomicAdd(&sum->nViol[FCap], 1ull);
    } else if (fam == FAnchor) {
        const int ax = R.gSub[r];
        val = comp(P[g] - v.P0[g], ax) / pr.D;
        G[0] = ax == 0 ? d3{1 / pr.D, 0, 0} : ax == 1 ? d3{0, 1 / pr.D, 0} : d3{0, 0, 1 / pr.D};
        wres = pr.wA * val;
        wgn = pr.wA;
        m = o = 0.5 * pr.wA * val * val;
    }
    V.gVal[r] = val;
    V.gTval[r] = tval;
    V.gAux[r] = aux;
    V.gFlag[r] = flag;
    V.wres[r] = wres;
    V.wgn[r] = wgn;
    V.merit[r] = m;
    V.obj[r] = o;
    V.strain2[r] = st2;
}

// atTrialRho: the accepted state's contacts at the rho its trial was evaluated with (accept).
// The contact rows' values and terms at the sample points (one storage for both sets). Their merit and
// objective are summed from these values (rowSumsPartialKernel); their gradients are formed for the
// active lists (fillActiveEntriesKernel) where cGradOn is set.
__global__ void evalContactRowsKernel(NcpView v, const d3* points, RowValuePair pair, RowSet which, const NcpControl* ctl, bool atTrialRho, EvalReductions* sum) {
    const NcpRows& R = v.rows;
    const NcpParams& pr = v.params;
    const NcpRowValues& V = rowValues(pair, which);
    const double rho = atTrialRho ? ctl->lastTrialRho : ctl->rho;
    const int nC = *R.nC;
    unsigned long long live = 0;
    for (long long cc = loopStart(); cc < nC; cc += loopStride()) {
        const int c = static_cast<int>(cc);
        d3 pa, pb;
        double s, t;
        const double d = contactDistance(points, R.cPair[c].x, R.cPair[c].y, R.cType[c], pa, pb, s, t);
        const double val = (pr.hardC - d) / pr.D;
        double wres = 0, wgn = 0, m = 0, tval = 0;
        hardTerm(val, R.cLam[c], rho, wres, wgn, m);
        if (wgn != 0 || R.cLam[c] > 0 || (pr.targets && d < pr.thr)) ++live;
        // The gradient only where the row has a term in the merit. A row without one keeps none until
        // it is evaluated again: a multiplier update leaves it at lambda = 0, a larger rho lowers
        // lambda / rho, and its target term depends on the distance alone.
        V.cGradOn[c] = wgn != 0 || (pr.targets && d < pr.thr) ? 1 : 0;
        if (pr.targets && (pr.thr - d) / pr.D > 0) wgn += pr.wT;
        if (pr.targets && d < pr.thr) {
            tval = (pr.thr - d) / pr.D;
            wres += pr.wT * tval;
            atomicMaxPositive(&sum->maxViol[FContactT], tval);
            if (tval > kViolationTolerance) atomicAdd(&sum->nViol[FContactT], 1ull);
        }
        atomicMaxPositive(&sum->maxViol[FContact], val);
        if (val > kViolationTolerance) atomicAdd(&sum->nViol[FContact], 1ull);
        V.cVal[c] = val;
        V.cTval[c] = tval;
        V.cWres[c] = wres;
        V.cWgn[c] = wgn;
    }
    if (live) atomicAdd(&sum->liveRows, live);
}

// The memory report's largest live set (EvalReductions::liveRows of a current-state evaluation).
__global__ void liveRowsMaxKernel(const EvalReductions* sum, NcpControl* ctl) {
    const long long live = static_cast<long long>(sum->liveRows);
    if (live > ctl->maxLiveRows) ctl->maxLiveRows = live;
}

// Lite trials: a trial's contact rows without per-row storage. The merit and objective sums are
// fixed-grid partials (rowSumsPartialKernel's pattern: deterministic); maxima and counts as evalContactRowsKernel.
__global__ void reduceContactRowsKernel(NcpView v, const d3* points, const NcpControl* ctl, EvalReductions* sum, double* blockSums) {
    const NcpRows& R = v.rows;
    const NcpParams& pr = v.params;
    const double rho = ctl->rho;
    __shared__ double shm[kSumThreads], sho[kSumThreads];
    double am = 0, ao = 0;
    const int nC = *R.nRows;  // every row: a row that is not live enters with multiplier 0
    for (long long cc = static_cast<long long>(blockIdx.x) * kSumThreads + threadIdx.x; cc < nC; cc += kSumStride) {
        const int c = static_cast<int>(cc);
        d3 pa, pb;
        double s, t;
        const double d = contactDistance(points, R.cPair[c].x, R.cPair[c].y, R.cType[c], pa, pb, s, t);
        const double val = (pr.hardC - d) / pr.D;
        double wres = 0, wgn = 0, m = 0, o = 0;
        hardTerm(val, rowMultiplier(R, c), rho, wres, wgn, m);
        if (pr.targets && d < pr.thr) {
            const double tval = (pr.thr - d) / pr.D;
            m += 0.5 * pr.wT * tval * tval;
            o = 0.5 * pr.wT * tval * tval;
            atomicMaxPositive(&sum->maxViol[FContactT], tval);
            if (tval > kViolationTolerance) atomicAdd(&sum->nViol[FContactT], 1ull);
        }
        atomicMaxPositive(&sum->maxViol[FContact], val);
        if (val > kViolationTolerance) atomicAdd(&sum->nViol[FContact], 1ull);
        am += m;
        ao += o;
    }
    shm[threadIdx.x] = am;
    sho[threadIdx.x] = ao;
    __syncthreads();
    for (int k = kSumThreads / 2; k > 0; k >>= 1) {
        if (static_cast<int>(threadIdx.x) < k) {
            shm[threadIdx.x] += shm[threadIdx.x + k];
            sho[threadIdx.x] += sho[threadIdx.x + k];
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        blockSums[blockIdx.x] = shm[0];
        blockSums[kSumBlocks + blockIdx.x] = sho[0];
    }
}

// The evaluation's reductions become the summary the controller compares (current or trial).
__global__ void summaryKernel(const EvalReductions* h, RowSet which, NcpControl* ctl) {
    NcpEvalSummary& out = which == RowSet::Trial ? ctl->trialEval : ctl->currentEval;
    out.merit = h->merit + (which == RowSet::Trial ? h->contactMerit : 0);
    out.objective = h->objective + (which == RowSet::Trial ? h->contactObjective : 0);
    for (int f = 0; f < knit::kRowFamilies; ++f) {
        out.maxViol[f] = h->maxViol[f];
        out.nViol[f] = static_cast<long long>(h->nViol[f]);
    }
    out.degenerate = static_cast<long long>(h->degenerate);
    out.notOwned = static_cast<long long>(h->notOwned);
    out.pinned = static_cast<long long>(h->pinnedBand);
    out.pinnedBandViol = h->pinnedBandViol;
    out.strainMax = h->strainMax;
    out.strainRms = h->edges ? sqrt(h->strain2 / static_cast<double>(h->edges)) : 0;
    if (which == RowSet::Trial) ctl->lastTrialRho = ctl->rho;
}

// Terms only (new rho or multipliers; values unchanged).
__global__ void termsKernel(NcpRows R, RowValuePair pair, NcpParams pr, const NcpControl* ctl) {
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const double rho = ctl->rho;
    const long long nR = R.nG + static_cast<long long>(*R.nC);
    for (long long rr = loopStart(); rr < nR; rr += loopStride()) {
        const int r = static_cast<int>(rr);
        double wres = 0, wgn = 0, m = 0;
        if (r < R.nG) {
            const int fam = R.gFam[r];
            if (V.gFlag[r]) {
            } else if (fam == FEdge) {
                wres = pr.wL * V.gVal[r];
                wgn = pr.wL;
                m = 0.5 * pr.wL * V.gVal[r] * V.gVal[r];
            } else if (fam == FAnchor) {
                wres = pr.wA * V.gVal[r];
                wgn = pr.wA;
                m = 0.5 * pr.wA * V.gVal[r] * V.gVal[r];
            } else {
                hardTerm(V.gVal[r], R.gLam[r], rho, wres, wgn, m);
                if (fam == FCurv && pr.targets) {
                    const double sc = pr.hardK / pr.kappaMax, tv = V.gTval[r];
                    const double raw = (V.gVal[r] + 1) * sc - 1;  // kappa / kappaMax - 1
                    if (tv > 0) {
                        wres += pr.wT * tv * sc;
                        m += 0.5 * pr.wT * tv * tv;
                    }
                    if (raw > 0) wgn += pr.wT * sc * sc;
                }
            }
        } else {
            const int c = r - R.nG;
            const double val = contactHardValue(V, contactValueScale(pr), c);
            hardTerm(val, R.cLam[c], rho, wres, wgn, m);
            const double tv = contactTargetValue(V, contactValueScale(pr), c);
            if (tv > 0) wres += pr.wT * tv;
            // (thr - d) / D from the hard value (hardC - d) / D
            if (pr.targets && val + (pr.thr - pr.hardC) / pr.D > 0) wgn += pr.wT;
            V.cWres[c] = wres;
            V.cWgn[c] = wgn;
            continue;  // its merit is summed from the values (rowSumsPartialKernel)
        }
        V.wres[r] = wres;
        V.wgn[r] = wgn;
        V.merit[r] = m;
    }
}

__global__ void multipliersKernel(NcpRows R, ContactValueScale scale, RowValuePair pair, const NcpControl* ctl) {
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const double rho = ctl->rho;
    const long long nR = R.nG + static_cast<long long>(*R.nC);
    for (long long rr = loopStart(); rr < nR; rr += loopStride()) {
        const int r = static_cast<int>(rr);
        if (r < R.nG) {
            const int fam = R.gFam[r];
            if (fam == FEdge || fam == FAnchor || V.gFlag[r]) continue;
            R.gLam[r] = fmax(0.0, R.gLam[r] + rho * V.gVal[r]);
        } else {
            const int c = r - R.nG;
            R.cLam[c] = fmax(0.0, R.cLam[c] + rho * contactHardValue(V, scale, c));
        }
    }
}

// Accept the trial: the trial row set becomes current, and its summary with it.
__global__ void acceptKernel(NcpControl* ctl) {
    ctl->current = 1 - ctl->current;
    ctl->currentEval = ctl->trialEval;
}

}  // namespace

void NcpDevice::State::enqueueLinearize(RowSet which) {
    const bool trial = which == RowSet::Trial;
    const d3* X = trial ? controls.trial() : controls.current();
    EvalReductions* sum = static_cast<EvalReductions*>(reductions.evaluation.data());
    double* partial = reductions.blocks.doubles();
    check(cudaMemsetAsync(sum, 0, sizeof(EvalReductions), stream), "ncp-gpu evaluate");
    const NcpView v = view();
    const RowValuePair pair = valuePair();
    if (bp.mode != kBandNone && N) bandFeetKernel<<<blocksFor(N), kThreadsPerBlock, 0, stream>>>(v, X, feet());
    if (nG) evalGeneralRowsKernel<<<blocksFor(nG), kThreadsPerBlock, 0, stream>>>(v, X, pair, which, feet(), control(), sum);
    if (trial) {
        reduceContactRowsKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(v, pointsAtQ(), control(), sum, partial);
        contactSumsFinalKernel<<<1, kSumBlocks, 0, stream>>>(partial, sum);
    } else if (storage == ContactStorage::LiveRows) {
        enqueueLiveEvaluation(EvaluationPenalty::CurrentRho, sum);
        liveRowsMaxKernel<<<1, 1, 0, stream>>>(sum, control());
    } else {
        evalContactRowsKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(v, pointsAtP(), pair, which, control(), false, sum);
        liveRowsMaxKernel<<<1, 1, 0, stream>>>(sum, control());
    }
    // A trial's per-row sums cover its general rows (its contacts were reduced above); the current state's
    // cover its contact list (LiveRows: the live rows, the only rows with a term).
    const RowSums rs{nG, trial ? nullptr : rows().nC, contacts.lambda.doubles(), pr.wT, pr.hardC, pr.thr, pr.D, pr.targets};
    rowSumsPartialKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(pair, which, rs, control(), partial);
    rowSumsFinalKernel<<<1, kSumBlocks, 0, stream>>>(partial, &sum->merit, &sum->objective, &sum->strain2);
    summaryKernel<<<1, 1, 0, stream>>>(sum, which, control());
    checkLaunch("ncp-gpu: evaluate");
}

void NcpDevice::State::enqueueAccept() {
    check(cudaMemcpyAsync(controls.P.data(), controls.Q.data(), size_t(nCtrl) * sizeof(d3), cudaMemcpyDeviceToDevice, stream), "ncp-gpu accept");
    check(cudaMemcpyAsync(samplesAtP.data(), samplesAtQ.data(), size_t(nSamples) * sizeof(d3), cudaMemcpyDeviceToDevice, stream), "ncp-gpu accept");
    acceptKernel<<<1, 1, 0, stream>>>(control());
    EvalReductions* sum = static_cast<EvalReductions*>(reductions.acceptEvaluation.data());
    check(cudaMemsetAsync(sum, 0, sizeof(EvalReductions), stream), "ncp-gpu accept");
    if (storage == ContactStorage::LiveRows) {
        enqueueLiveEvaluation(EvaluationPenalty::LastTrialRho, sum);
    } else {
        evalContactRowsKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(view(), pointsAtP(), valuePair(), RowSet::Current, control(), true, sum);
    }
    liveRowsMaxKernel<<<1, 1, 0, stream>>>(sum, control());
    checkLaunch("ncp-gpu: accept");
}

void NcpDevice::State::enqueueMultipliers() {
    multipliersKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(rows(), contactValueScale(pr), valuePair(), control());
}

void NcpDevice::State::enqueueTerms() {
    termsKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(rows(), valuePair(), pr, control());
    double* partial = reductions.blocks.doubles();
    const RowSums rs{nG, &counts()->contactRows, contacts.lambda.doubles(), pr.wT};
    rowSumsPartialKernel<<<kSumBlocks, kSumThreads, 0, stream>>>(valuePair(), RowSet::Current, rs, control(), partial);
    rowSumsFinalKernel<<<1, kSumBlocks, 0, stream>>>(partial, &control()->currentEval.merit, nullptr, nullptr);
}

}  // namespace knit::gpu
