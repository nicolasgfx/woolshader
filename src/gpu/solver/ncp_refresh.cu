// ncp-gpu: broad phase, contact rows and the active rows of the products (ncp.h). The kernels, then the
// sequences that enqueue them, in run order: the motion since the last refresh; the refresh (binding,
// grid, candidates, contact rows, multipliers carried over); the active lists of one assembly. Every
// sequence runs on the call's stream without a host read, so the graph of ncp_control.cu captures it.
#include "gpu/solver/ncp_device.cuh"

namespace knit::gpu {

namespace {

// --- motion since the last refresh
__global__ void motionKernel(int N, const int* ctrlOfVar, const d3* P, const d3* Xref, NcpControl* ctl) {
    const int a = threadItem();
    if (a >= N) return;
    const int g = ctrlOfVar[a];
    atomicMaxPositive(&ctl->motion, norm(P[g] - Xref[g]));
}

// --- binding continuation at refresh (KnitProblem::refresh): advanced only on Owned.
__global__ void refreshBindingKernel(NcpView v, const d3* P, const NcpControl* ctl, DoubleProjectionState* state, unsigned char* result) {
    if (ctl->resume) return;  // done by the refresh that overflowed (solveNcp's resumed launch)
    const int a = threadItem();
    if (a >= v.variables) return;
    const int g = v.ctrlOfVar[a];
    DoubleProjectionState next;
    const int r = continueFrom(v.band.mesh, state[g], P[g], next, v.band.cont);
    result[g] = static_cast<unsigned char>(r);
    if (r == 0) state[g] = next;
}

// The first capacity a refresh exceeded, and what it needed (kept: the iteration stops on it).
__device__ void recordOverflow(NcpCounts* counts, int kind, long long needed, long long capacity) {
    if (counts->overflow) return;
    counts->overflow = kind;
    counts->overflowNeeded = needed;
    counts->overflowCapacity = capacity;
}

__global__ void refreshBeginKernel(NcpCounts* counts) {
    counts->previousRows = counts->contactRows;
    counts->candidates = 0;
    counts->nearRows = 0;
    counts->carriedFar = 0;
    counts->carried = 0;
}

// The sample points of controls X (evalSample).
__global__ void samplePointsKernel(const d3* X, NcpSamples S, d3* points) {
    const int i = threadItem();
    if (i < S.n) points[i] = evalSample(X, S, i);
}

// A segment's box (from the sample points); a segment that takes no part gets an inverted box, which
// the grid leaves out.
__global__ void sampleBoxesKernel(const d3* points, NcpSamples S, float* boxLo, float* boxHi) {
    const int i = threadItem();
    if (i >= S.n) return;
    if (!S.segOk[i]) {
        boxLo[3 * i] = 1;
        boxHi[3 * i] = 0;
        return;
    }
    const f3 a = tof(points[i]), b = tof(points[i + 1]);
    const f3 lo = vmin(a, b), hi = vmax(a, b);
    boxLo[3 * i] = lo.x; boxLo[3 * i + 1] = lo.y; boxLo[3 * i + 2] = lo.z;
    boxHi[3 * i] = hi.x; boxHi[3 * i + 1] = hi.y; boxHi[3 * i + 2] = hi.z;
}

// --- the grid of fixed shape (FixedGrid): cells clamped into the range
struct GridShapeView {
    float inv;
    long long lo[3], n[3];
};

__device__ __forceinline__ long long clampedCell(float v, float inv, long long lo, long long n) {
    const long long c = static_cast<long long>(floorf(v * inv)) - lo;
    return c < 0 ? 0 : c >= n ? n - 1 : c;
}

__device__ __forceinline__ void itemCells(const float* boxLo, const float* boxHi, int i, const GridShapeView& g, long long c0[3], long long c1[3]) {
    for (int a = 0; a < 3; ++a) {
        c0[a] = clampedCell(boxLo[3 * i + a], g.inv, g.lo[a], g.n[a]);
        c1[a] = clampedCell(boxHi[3 * i + a], g.inv, g.lo[a], g.n[a]);
    }
}

__global__ void gridCountKernel(const float* boxLo, const float* boxHi, int items, GridShapeView g, unsigned long long* count) {
    const int i = threadItem();
    if (i == 0) count[items] = 0;
    if (i >= items) return;
    if (boxLo[3 * i] > boxHi[3 * i]) {
        count[i] = 0;
        return;
    }
    long long c0[3], c1[3];
    itemCells(boxLo, boxHi, i, g, c0, c1);
    count[i] = static_cast<unsigned long long>((c1[0] - c0[0] + 1) * (c1[1] - c0[1] + 1) * (c1[2] - c0[2] + 1));
}

__global__ void gridCheckKernel(const unsigned long long* offsets, int items, long long capacity, NcpCounts* counts) {
    const unsigned long long total = offsets[items];
    counts->registrations = total;
    if (total > static_cast<unsigned long long>(capacity)) recordOverflow(counts, kCapRegistrations, static_cast<long long>(total), capacity);
}

// One key per cell of the item: (cell << 32) | item, so sorting orders the items within a cell.
__global__ void gridEmitKernel(const float* boxLo, const float* boxHi, int items, GridShapeView g, const unsigned long long* offsets, long long capacity,
                               unsigned long long* keys, unsigned* cellCounts) {
    const int i = threadItem();
    if (i >= items || boxLo[3 * i] > boxHi[3 * i] || offsets[i + 1] > static_cast<unsigned long long>(capacity)) return;
    long long c0[3], c1[3];
    itemCells(boxLo, boxHi, i, g, c0, c1);
    unsigned long long at = offsets[i];
    for (long long x = c0[0]; x <= c1[0]; ++x)
        for (long long y = c0[1]; y <= c1[1]; ++y)
            for (long long z = c0[2]; z <= c1[2]; ++z) {
                const unsigned cell = static_cast<unsigned>((x * g.n[1] + y) * g.n[2] + z);
                keys[at++] = (static_cast<unsigned long long>(cell) << 32) | static_cast<unsigned>(i);
                atomicAdd(&cellCounts[cell], 1u);
            }
}

// A sorted key becomes an entry: the item with its "lowest cell" flags (GridCursor's de-duplication).
__global__ void gridEntriesKernel(const float* boxLo, const float* boxHi, GridShapeView g, const unsigned long long* keys, const NcpCounts* counts,
                                  long long capacity, unsigned* entries) {
    const long long n = static_cast<long long>(counts->registrations) < capacity ? static_cast<long long>(counts->registrations) : capacity;
    for (long long p = loopStart(); p < n; p += loopStride()) {
        const long long cell = static_cast<long long>(keys[p] >> 32);
        const int i = static_cast<int>(keys[p] & 0xffffffffu);
        const long long z = cell % g.n[2], y = (cell / g.n[2]) % g.n[1], x = cell / (g.n[2] * g.n[1]);
        long long c0[3], c1[3];
        itemCells(boxLo, boxHi, i, g, c0, c1);
        const unsigned flags = (x == c0[0] ? 1u : 0u) | (y == c0[1] ? 2u : 0u) | (z == c0[2] ? 4u : 0u);
        entries[p] = static_cast<unsigned>(i) | flags << 29;
    }
}

// --- candidates and contact rows
// Which candidates become rows: every one (CompactRows), or those whose segments are within the near
// cutoff and those carrying a positive multiplier (LiveRows, NcpOptions::contactStorage).
enum class ResidentRule { EveryCandidate, NearOrCarried };
struct CandidateQuery {
    GridView grid;
    NcpSamples samples;
    float pad;
    double window;
    ResidentRule rule;
    double nearCutoff;                     // NearOrCarried: nearCutoff()
    const unsigned long long* carriedKeys; // NearOrCarried: the refresh's sorted carry-over keys
    const int* carriedCount;               // their number (NcpCounts::carried)
};

// KnitProblem::refresh's candidate rules: active i, each active pair once, same yarn excluded
// when the arc gap between the two segments is below the window.
__device__ __forceinline__ bool isCandidatePair(const NcpSamples& S, int i, int j, double window) {
    if (j == i || (S.active[j] && j < i)) return false;
    if (S.yarn[j] != S.yarn[i]) return true;
    const double gap = j > i ? S.s[j] - S.s[i + 1] : S.s[i] - S.s[j + 1];
    return gap >= window;
}

// The query of item i in the fixed grid: its box grown by the pad (the cursor clamps to the range).
__device__ __forceinline__ void queryBox(const CandidateQuery& q, int i, f3& lo, f3& hi) {
    const f3 pv{q.pad, q.pad, q.pad};
    lo = fromF(q.grid.boxLo + 3 * i) - pv;
    hi = fromF(q.grid.boxHi + 3 * i) + pv;
}

// The contact type of candidate pair (i, j) by the contact model (KnitProblem::refresh): 0 segment -
// segment, 1 point - point, 2 point i - segment j (hybrid: near-parallel segments). A point-segment
// candidate gives two rows, (i, j) and (j, i).
__device__ __forceinline__ int contactTypeOf(const NcpParams& pr, const d3* points, int i, int j) {
    const int model = pr.contactModel;
    int ty = model == 0 ? 0 : model == 1 ? 1 : 2;
    if (model == 3) {
        const d3 a = points[i + 1] - points[i], b = points[j + 1] - points[j];
        const double la = norm(a), lb = norm(b);
        ty = la > 0 && lb > 0 && fabs(dot(a, b)) > pr.parallelCos * la * lb ? 2 : 0;
    }
    return ty;
}

// The near cutoff of the memory report: max(hardC, thr) + 2 delta + guard (the pad is thr + 2 delta +
// guard). Until the next refresh every control stays within delta of Xref (the step's coverage limit;
// the cap projection is non-expansive), a sample point is a convex combination of controls, so two
// sample segments farther apart than this at the refresh stay farther than max(hardC, thr) + guard at
// every trial and along every straight motion: no term in the merit and no swept touch.
__device__ __forceinline__ double nearCutoff(const CandidateQuery& q, const NcpParams& pr) {
    return static_cast<double>(q.pad) + fmax(0.0, pr.hardC - pr.thr);
}

// Whether the sample segments i and j are closer than the cutoff (the full segments, whatever the row's
// contact type: a row's own distance is never below its segments' distance, and the swept test uses them).
__device__ __forceinline__ bool segmentsNear(const d3* points, int i, int j, double cutoff) {
    double s, t;
    return segmentDistance(points[i], points[i + 1], points[j], points[j + 1], s, t) < cutoff;
}

// Whether row key k carries a positive multiplier into this refresh (binary search of the sorted table).
__device__ __forceinline__ bool keyCarried(const CandidateQuery& q, unsigned long long k) {
    int lo = 0, hi = *q.carriedCount;
    while (lo < hi) {
        const int mid = (lo + hi) >> 1;
        if (q.carriedKeys[mid] < k) lo = mid + 1;
        else hi = mid;
    }
    return lo < *q.carriedCount && q.carriedKeys[lo] == k;
}

// Whether candidate (i, j) of contact type ty becomes rows. A point-segment candidate's two rows are kept
// or left out together. A left-out candidate has no multiplier, and its segments stay farther apart than
// max(hardC, thr) until the next refresh (nearCutoff): its rows would have neither a term nor a touch.
__device__ __forceinline__ bool candidateResident(const CandidateQuery& q, const d3* points, int i, int j, int ty) {
    if (q.rule == ResidentRule::EveryCandidate) return true;
    if (segmentsNear(points, i, j, q.nearCutoff)) return true;
    return keyCarried(q, contactKey(i, j, ty)) || (ty == 2 && keyCarried(q, contactKey(j, i, ty)));
}

// Whether contact row (i, j, type) has a term in the merit at multiplier 0 (evalContactRowsKernel's
// test with lambda = 0): its hard row violated, or its distance below the target threshold.
__device__ __forceinline__ bool activeAtZeroMultiplier(const NcpParams& pr, const d3* X, int i, int j, int ty) {
    d3 pa, pb;
    double s, t;
    const double d = contactDistance(X, i, j, ty, pa, pb, s, t);
    return (pr.hardC - d) / pr.D > 0 || (pr.targets && d < pr.thr);
}

// The memory report's starting sizes, summed with integer atomics: totals[0] candidates, [1] contact
// rows, [2] contact rows with a term in the merit at multipliers 0, [3] contact rows whose segments are
// closer than the near cutoff.
__global__ void startCountsKernel(CandidateQuery q, NcpParams pr, const d3* points, unsigned long long* totals) {
    const NcpSamples& S = q.samples;
    const int i = threadItem();
    if (i >= S.n || !S.segOk[i] || !S.active[i]) return;
    const double cutoff = nearCutoff(q, pr);
    unsigned long long candidates = 0, rows = 0, active = 0, near = 0;
    f3 lo, hi;
    queryBox(q, i, lo, hi);
    for (GridCursor cursor(q.grid, lo, hi); cursor.next();) {
        const int j = cursor.item();
        if (!isCandidatePair(S, i, j, q.window)) continue;
        ++candidates;
        const int ty = contactTypeOf(pr, points, i, j);
        const int rowsOfPair = ty == 2 ? 2 : 1;
        rows += rowsOfPair;
        if (segmentsNear(points, i, j, cutoff)) near += rowsOfPair;
        active += activeAtZeroMultiplier(pr, points, i, j, ty) ? 1 : 0;
        if (ty == 2) active += activeAtZeroMultiplier(pr, points, j, i, ty) ? 1 : 0;
    }
    atomicAdd(&totals[0], candidates);
    atomicAdd(&totals[1], rows);
    atomicAdd(&totals[2], active);
    atomicAdd(&totals[3], near);
}

// The contact rows of sample i's candidates (two for a point-segment candidate); the candidates are
// summed for the statistics.
__global__ void countRowsKernel(CandidateQuery q, NcpParams pr, const d3* points, unsigned long long* rowsPerSample, NcpCounts* counts) {
    const NcpSamples& S = q.samples;
    const int i = threadItem();
    if (i == 0) rowsPerSample[S.n] = 0;
    if (i >= S.n) return;
    unsigned long long rows = 0, candidates = 0;
    if (S.segOk[i] && S.active[i]) {
        f3 lo, hi;
        queryBox(q, i, lo, hi);
        for (GridCursor cursor(q.grid, lo, hi); cursor.next();) {
            const int j = cursor.item();
            if (!isCandidatePair(S, i, j, q.window)) continue;
            ++candidates;
            const int ty = contactTypeOf(pr, points, i, j);
            if (!candidateResident(q, points, i, j, ty)) continue;
            rows += ty == 2 ? 2 : 1;
        }
    }
    rowsPerSample[i] = rows;
    if (candidates) atomicAdd(&counts->candidates, candidates);
}

// The new row count, or the overflow. On overflow the previous rows stay as they are (nothing is
// written over them): the controller ends the iteration at the last accepted state.
__global__ void rowCheckKernel(const unsigned long long* rowOffset, int samples, long long rowCapacity, NcpCounts* counts, NcpControl* ctl) {
    ++ctl->refreshes;
    if (!counts->overflow) {
        const long long rows = static_cast<long long>(rowOffset[samples]);
        if (rows > rowCapacity) recordOverflow(counts, kCapContactRows, rows, rowCapacity);
        else counts->contactRows = static_cast<int>(rows);
    }
    if (counts->contactRows > ctl->maxContactRows) ctl->maxContactRows = counts->contactRows;
    if (static_cast<long long>(counts->candidates) > ctl->maxCandidates) ctl->maxCandidates = static_cast<long long>(counts->candidates);
}

// The previous rows' keys and multipliers for the carry-over sort; the places past them get the
// largest key, so the sort over the capacity puts them last.
__global__ void oldKeysKernel(NcpRows R, const NcpCounts* counts, long long capacity, unsigned long long* keys, double* lam) {
    const int nOld = counts->previousRows;
    for (long long c = loopStart(); c < capacity; c += loopStride()) {
        if (c < nOld) {
            keys[c] = contactKey(R.cPair[c].x, R.cPair[c].y, R.cType[c]);
            lam[c] = R.cLam[c];
        } else {
            keys[c] = ~0ull;
            lam[c] = 0;
        }
    }
}

// Sample i's rows at its offset, in countRowsKernel's order: (i, j), and (j, i) right after it for a
// point-segment candidate. A row touching more than kContactRowEntries free controls is recorded (the
// solve refuses it).
__global__ void fillRowsKernel(CandidateQuery q, NcpParams pr, const d3* points, const int* varOf, const unsigned long long* offsets, NcpCounts* counts,
                               int2* pairs, unsigned char* type) {
    if (counts->overflow) return;
    const NcpSamples& S = q.samples;
    const int i = threadItem();
    if (i >= S.n || !S.segOk[i] || !S.active[i]) return;
    f3 lo, hi;
    queryBox(q, i, lo, hi);
    unsigned long long at = offsets[i];
    bool entryOverflow = false;
    int controls[kContactRowEntries];
    const double cutoff = nearCutoff(q, pr);
    unsigned long long near = 0;
    for (GridCursor cursor(q.grid, lo, hi); cursor.next();) {
        const int j = cursor.item();
        if (!isCandidatePair(S, i, j, q.window)) continue;
        const int ty = contactTypeOf(pr, points, i, j);
        if (!candidateResident(q, points, i, j, ty)) continue;
        pairs[at] = make_int2(i, j);
        type[at] = static_cast<unsigned char>(ty);
        ++at;
        if (ty == 2) {
            pairs[at] = make_int2(j, i);
            type[at] = 2;
            ++at;
        }
        if (segmentsNear(points, i, j, cutoff)) near += ty == 2 ? 2 : 1;
        contactControls(S, varOf, i, j, controls, entryOverflow);  // (j, i) touches the same controls
    }
    if (entryOverflow) counts->contactEntryOverflow = 1;
    if (near) atomicAdd(&counts->nearRows, near);
}

// Multipliers carried over by row key (the previous keys sorted). The memory report counts the rows
// that carried a positive multiplier although their segments are beyond the near cutoff.
__global__ void carryMultipliersKernel(const unsigned long long* oldKeys, const double* oldLambda, NcpRows R, const d3* points, double cutoff, NcpCounts* counts) {
    if (counts->overflow) return;
    const int nC = counts->contactRows, nOld = counts->previousRows;
    unsigned long long far = 0;
    for (long long cc = loopStart(); cc < nC; cc += loopStride()) {
        const int c = static_cast<int>(cc);
        const unsigned long long k = contactKey(R.cPair[c].x, R.cPair[c].y, R.cType[c]);
        int lo = 0, hi = nOld;
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (oldKeys[mid] < k) lo = mid + 1;
            else hi = mid;
        }
        const double lambda = lo < nOld && oldKeys[lo] == k ? oldLambda[lo] : 0.0;
        R.cLam[c] = lambda;
        if (lambda > 0 && !segmentsNear(points, R.cPair[c].x, R.cPair[c].y, cutoff)) ++far;
    }
    if (far) atomicAdd(&counts->carriedFar, far);
}

// The memory report's largest broad phase (after the carry-over: its counts are complete).
__global__ void refreshCountsKernel(const NcpCounts* counts, NcpControl* ctl) {
    if (counts->overflow) return;
    const long long near = static_cast<long long>(counts->nearRows), far = static_cast<long long>(counts->carriedFar);
    if (near > ctl->maxNearRows) ctl->maxNearRows = near;
    if (far > ctl->maxCarriedFar) ctl->maxCarriedFar = far;
}

// --- the active rows of one assembly
// Per row: 1 when it has a term in the merit (and every edge: its stretch Hessian); 0 past the rows up
// to the capacity. Their exclusive scan gives each active row its slot, in row order.
__global__ void activeFlagsKernel(NcpView v, RowValuePair pair, RowSet which, long long rowCapacity, int* flags) {
    const NcpRows& R = v.rows;
    const NcpRowValues& V = rowValues(pair, which);
    const long long nR = R.nG + static_cast<long long>(*R.nC);
    for (long long rr = loopStart(); rr < rowCapacity; rr += loopStride()) {
        const int r = static_cast<int>(rr);
        int flag = 0;
        if (rr < nR) {
            const bool edge = r < R.nG && R.gFam[r] == FEdge;
            flag = rowWgn(V, R.nG, r) != 0 || rowWres(V, R.nG, r) != 0 || edge ? 1 : 0;
        }
        flags[r] = flag;
    }
}

// Each active row into its slot; its entries counted per variable.
__global__ void activeSlotsKernel(NcpView v, const int* slotOf, long long rowCapacity, NcpCounts* counts, int* rows, int* perVariable) {
    const NcpRows& R = v.rows;
    const long long nR = R.nG + static_cast<long long>(*R.nC);
    if (loopStart() == 0) counts->activeRows = slotOf[rowCapacity];
    for (long long rr = loopStart(); rr < nR; rr += loopStride()) {
        const int k = slotOf[rr];
        if (slotOf[rr + 1] == k) continue;  // not active
        const int r = static_cast<int>(rr);
        rows[k] = r;
        if (r < R.nG) {
            for (int e = 0; e < R.gN[r]; ++e) atomicAdd(&perVariable[R.gVar[kGeneralRowEntries * r + e]], 1);
        } else {
            const int2 p = R.cPair[contactRowOf(R, r - R.nG)];
            int controls[kContactRowEntries];
            bool overflow = false;
            const int n = contactControls(v.samples, v.varOf, p.x, p.y, controls, overflow);
            for (int e = 0; e < n; ++e) atomicAdd(&perVariable[v.varOf[controls[e]]], 1);
        }
    }
}

// The entries (slot << 4) | e per variable, and each contact slot's gradient at the current sample
// points: the unit normal and the closest points' parameters of a row that had a term in the merit
// when it was evaluated (evalContactRowsKernel), else a zero normal.
__global__ void fillActiveEntriesKernel(NcpView v, RowValuePair pair, RowSet which, const d3* points, const NcpCounts* counts, const int* rows, const int* start,
                                        int* cursor, int* entry, d3* normal, double* st) {
    const NcpRows& R = v.rows;
    const NcpRowValues& V = rowValues(pair, which);
    const int nAct = counts->activeRows;
    for (long long kk = loopStart(); kk < nAct; kk += loopStride()) {
        const int k = static_cast<int>(kk);
        const int r = rows[k];
        if (r < R.nG) {
            for (int e = 0; e < R.gN[r]; ++e) {
                const int a = R.gVar[kGeneralRowEntries * r + e];
                entry[start[a] + atomicAdd(&cursor[a], 1)] = (k << 4) | e;
            }
            continue;
        }
        const int c = r - R.nG;
        const int row = contactRowOf(R, c);
        const int2 p = R.cPair[row];
        int controls[kContactRowEntries];
        bool overflow = false;
        const int n = contactControls(v.samples, v.varOf, p.x, p.y, controls, overflow);
        for (int e = 0; e < n; ++e) {
            const int a = v.varOf[controls[e]];
            entry[start[a] + atomicAdd(&cursor[a], 1)] = (k << 4) | e;
        }
        d3 nrm{0, 0, 0};
        double s = 0, t = 0;
        if (V.cGradOn[c]) {
            d3 pa, pb;
            const double d = contactDistance(points, p.x, p.y, R.cType[row], pa, pb, s, t);
            nrm = d > 1e-9 ? (pa - pb) / d : d3{0, 0, 1};
        }
        normal[k] = nrm;
        st[2 * k] = s;
        st[2 * k + 1] = t;
    }
}

// Each variable's entries in slot order, which is row order (they were appended in any order): every
// later sum over a variable's rows adds in the same order in every run.
__global__ void sortActiveEntriesKernel(int N, const int* start, int* entry) {
    const int a = threadItem();
    if (a >= N) return;
    const int b = start[a], e = start[a + 1];
    for (int i = b + 1; i < e; ++i) {
        const int x = entry[i];
        int j = i - 1;
        while (j >= b && entry[j] > x) {
            entry[j + 1] = entry[j];
            --j;
        }
        entry[j + 1] = x;
    }
}

// The same sort, one block per variable (LiveRows; the M1 profile measured the one-thread insertion sort at
// 20 % of a tube call): the variable's entries in shared
// memory, a bitonic sort padded to a power of two, written back. The entries are distinct, so the result is
// the insertion sort's. A variable with more entries than the shared array holds is sorted as before.
constexpr int kSortThreads = 256;
constexpr int kSortSharedEntries = 8192;  // 32 KB of shared memory per block
__global__ void sortActiveEntriesByBlockKernel(const int* start, int* entry) {
    __shared__ int sh[kSortSharedEntries];
    const int a = blockIdx.x;
    const int b = start[a], n = start[a + 1] - b;
    if (n <= 1) return;
    if (n > kSortSharedEntries) {
        if (threadIdx.x != 0) return;
        for (int i = b + 1; i < b + n; ++i) {
            const int x = entry[i];
            int j = i - 1;
            while (j >= b && entry[j] > x) {
                entry[j + 1] = entry[j];
                --j;
            }
            entry[j + 1] = x;
        }
        return;
    }
    int padded = 1;
    while (padded < n) padded <<= 1;
    for (int i = threadIdx.x; i < padded; i += kSortThreads) sh[i] = i < n ? entry[b + i] : 0x7fffffff;
    __syncthreads();
    for (int size = 2; size <= padded; size <<= 1) {
        for (int stride = size >> 1; stride > 0; stride >>= 1) {
            for (int i = threadIdx.x; i < padded; i += kSortThreads) {
                const int j = i ^ stride;
                if (j <= i) continue;
                const bool ascending = (i & size) == 0;
                const int x = sh[i], y = sh[j];
                if ((x > y) == ascending) {
                    sh[i] = y;
                    sh[j] = x;
                }
            }
            __syncthreads();
        }
    }
    for (int i = threadIdx.x; i < n; i += kSortThreads) entry[b + i] = sh[i];
}

__global__ void activeTotalKernel(const int* start, int N, NcpCounts* counts, NcpControl* ctl) {
    counts->activeEntries = start[N];
    if (counts->activeRows > ctl->maxActiveRows) ctl->maxActiveRows = counts->activeRows;
    if (counts->activeEntries > ctl->maxActiveEntries) ctl->maxActiveEntries = counts->activeEntries;
}

}  // namespace

void NcpDevice::State::enqueueSamplePoints(RowSet which) {
    const d3* X = which == RowSet::Trial ? controls.trial() : controls.current();
    d3* out = static_cast<d3*>(which == RowSet::Trial ? samplesAtQ.data() : samplesAtP.data());
    if (nSamples) samplePointsKernel<<<blocksFor(nSamples), kThreadsPerBlock, 0, stream>>>(X, S, out);
}

void NcpDevice::State::enqueueMotion() {
    check(cudaMemsetAsync(&control()->motion, 0, sizeof(double), stream), "ncp-gpu motion");
    if (N) motionKernel<<<blocksFor(N), kThreadsPerBlock, 0, stream>>>(N, controls.ctrlOfVar.ints(), controls.current(), controls.reference(), control());
}

void NcpDevice::State::enqueueRefresh() {
    if (hasBind && N > 0)
        refreshBindingKernel<<<blocksFor(N), kThreadsPerBlock, 0, stream>>>(view(), controls.current(), control(), static_cast<DoubleProjectionState*>(bindState.data()),
                                                                             bindResult.chars());
    check(cudaMemcpyAsync(controls.Xref.data(), controls.P.data(), size_t(nCtrl) * sizeof(d3), cudaMemcpyDeviceToDevice, stream), "ncp-gpu refresh");
    // Candidates at P: pad thr + 2 delta, plus a guard for the float boxes' rounding.
    enqueueContactRows(static_cast<float>(pr.thr + 2 * delta) + 1e-3f);
}

// The broad phase's grid of fixed shape over the segments at P (FixedGrid): boxes, registrations
// sorted by cell, cell starts and entries.
void NcpDevice::State::enqueueGrid() {
    const int Ns = nSamples;
    NcpCounts* cnt = counts();
    FixedGrid& g = grid;
    const GridShapeView shape{1.0f / g.cell, {g.lo[0], g.lo[1], g.lo[2]}, {g.n[0], g.n[1], g.n[2]}};
    const long long nCells = g.n[0] * g.n[1] * g.n[2];
    sampleBoxesKernel<<<blocksFor(Ns), kThreadsPerBlock, 0, stream>>>(pointsAtP(), S, g.boxLo.floats(), g.boxHi.floats());
    gridCountKernel<<<blocksFor(Ns + 1), kThreadsPerBlock, 0, stream>>>(g.boxLo.floats(), g.boxHi.floats(), Ns, shape, g.counts.ulongs());
    prims.countsToOffsetsOnStream(g.counts.ulongs(), g.offsets.ulongs(), Ns, stream);
    gridCheckKernel<<<1, 1, 0, stream>>>(g.offsets.ulongs(), Ns, capacity.registrations, cnt);
    check(cudaMemsetAsync(g.keys.data(), 0xff, size_t(capacity.registrations) * 8, stream), "ncp-gpu grid");
    check(cudaMemsetAsync(g.cellCounts.data(), 0, size_t(nCells + 1) * 4, stream), "ncp-gpu grid");
    gridEmitKernel<<<blocksFor(Ns), kThreadsPerBlock, 0, stream>>>(g.boxLo.floats(), g.boxHi.floats(), Ns, shape, g.offsets.ulongs(), capacity.registrations,
                                                                  g.keys.ulongs(), g.cellCounts.uints());
    prims.sortKeysOnStream({g.keys.ulongs(), g.sorted.ulongs(), capacity.registrations, g.endBit}, stream);
    prims.countsToOffsetsOnStream(g.cellCounts.uints(), g.cellStart.uints(), nCells, stream);
    gridEntriesKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(g.boxLo.floats(), g.boxHi.floats(), shape, g.sorted.ulongs(), cnt, capacity.registrations,
                                                                    g.entries.uints());
}

void NcpDevice::State::enqueueStartCounts(float pad, unsigned long long* totals) {
    check(cudaMemsetAsync(totals, 0, 4 * sizeof(unsigned long long), stream), "ncp-gpu start counts");
    const CandidateQuery query{grid.view, S, pad, window};
    if (nSamples) startCountsKernel<<<blocksFor(nSamples), kThreadsPerBlock, 0, stream>>>(query, pr, pointsAtP(), totals);
    checkLaunch("ncp-gpu: start counts");
}

// The previous rows' keys and multipliers, and their sorted copies, in the active entries' memory:
// nothing reads the active lists between a refresh and the next assembly, which builds them again.
CarryBuffers NcpDevice::State::carryBuffers() const {
    const long long C = contactListCapacity();  // CompactRows: the rows; LiveRows: the live list (its carry-over)
    unsigned long long* keys = active.entry.ulongs();
    double* lam = reinterpret_cast<double*>(keys + 2 * C);
    return {keys, keys + C, lam, lam + C, C};
}

// Contact rows of the candidate pairs at P with the given pad (multipliers carried by key): counted per
// sample, checked against the capacity, then written in place; an overflow leaves the previous rows.
void NcpDevice::State::enqueueContactRows(float pad) {
    const int Ns = nSamples;
    NcpCounts* cnt = counts();
    refreshBeginKernel<<<1, 1, 0, stream>>>(cnt);
    enqueueGrid();
    const double cutoff = static_cast<double>(pad) + std::max(0.0, pr.hardC - pr.thr);
    if (storage == ContactStorage::LiveRows) {
        // The live rows' carried multipliers first (they decide which far candidates stay), then the rows
        // near or carried, then the carried multipliers handed to the new rows.
        const CarryTable table = enqueueLiveCarryTable();
        const CandidateQuery query{grid.view, S, pad, window, ResidentRule::NearOrCarried, cutoff, table.keys, &cnt->carried};
        countRowsKernel<<<blocksFor(Ns + 1), kThreadsPerBlock, 0, stream>>>(query, pr, pointsAtP(), contacts.sampleRows.ulongs(), cnt);
        prims.countsToOffsetsOnStream(contacts.sampleRows.ulongs(), contacts.sampleRowOffset.ulongs(), Ns, stream);
        rowCheckKernel<<<1, 1, 0, stream>>>(contacts.sampleRowOffset.ulongs(), Ns, capacity.contactRows, cnt, control());
        fillRowsKernel<<<blocksFor(Ns), kThreadsPerBlock, 0, stream>>>(query, pr, pointsAtP(), controls.varOf.ints(), contacts.sampleRowOffset.ulongs(), cnt,
                                                                      static_cast<int2*>(contacts.pair.data()), contacts.type.chars());
        enqueueLiveCarryOver(table);
        refreshCountsKernel<<<1, 1, 0, stream>>>(cnt, control());
        checkLaunch("ncp-gpu: contact rows");
        return;
    }
    // ---- the previous rows' keys and multipliers, sorted by key, for the carry-over
    const CarryBuffers carry = carryBuffers();
    oldKeysKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(rows(), cnt, carry.n, carry.keys, carry.lambda);
    const bool sortedInAlt = prims.sortPairsAlternatingOnStream({carry.keys, carry.keysAlt, carry.lambda, carry.lambdaAlt, carry.n}, stream);
    const unsigned long long* oldKeys = sortedInAlt ? carry.keysAlt : carry.keys;
    const double* oldLambda = sortedInAlt ? carry.lambdaAlt : carry.lambda;
    // ---- the new rows: per sample, their offsets, the capacity, the rows themselves
    const CandidateQuery query{grid.view, S, pad, window};
    countRowsKernel<<<blocksFor(Ns + 1), kThreadsPerBlock, 0, stream>>>(query, pr, pointsAtP(), contacts.sampleRows.ulongs(), cnt);
    prims.countsToOffsetsOnStream(contacts.sampleRows.ulongs(), contacts.sampleRowOffset.ulongs(), Ns, stream);
    rowCheckKernel<<<1, 1, 0, stream>>>(contacts.sampleRowOffset.ulongs(), Ns, capacity.contactRows, cnt, control());
    fillRowsKernel<<<blocksFor(Ns), kThreadsPerBlock, 0, stream>>>(query, pr, pointsAtP(), controls.varOf.ints(), contacts.sampleRowOffset.ulongs(), cnt,
                                                                  static_cast<int2*>(contacts.pair.data()), contacts.type.chars());
    carryMultipliersKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(oldKeys, oldLambda, rows(), pointsAtP(), cutoff, cnt);
    refreshCountsKernel<<<1, 1, 0, stream>>>(cnt, control());
    checkLaunch("ncp-gpu: contact rows");
}

// Per variable, the active rows of the given row set in slots, sorted (ActiveLists).
void NcpDevice::State::enqueueActiveLists(RowSet which) {
    NcpCounts* cnt = counts();
    const long long rowCapacity = nG + contactListCapacity();
    check(cudaMemsetAsync(active.count.data(), 0, size_t(N + 1) * 4, stream), "ncp-gpu active rows");
    check(cudaMemsetAsync(active.cursor.data(), 0, size_t(N + 1) * 4, stream), "ncp-gpu active rows");
    activeFlagsKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(view(), valuePair(), which, rowCapacity, active.slotOf.ints());
    prims.countsToOffsetsOnStream(active.slotOf.ints(), active.slotOf.ints(), rowCapacity, stream);  // in place: flags -> slots
    activeSlotsKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(view(), active.slotOf.ints(), rowCapacity, cnt, active.rows.ints(), active.count.ints());
    prims.countsToOffsetsOnStream(active.count.ints(), active.start.ints(), N, stream);
    activeTotalKernel<<<1, 1, 0, stream>>>(active.start.ints(), N, cnt, control());
    fillActiveEntriesKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(view(), valuePair(), which, pointsAtP(), cnt, active.rows.ints(), active.start.ints(),
                                                                          active.cursor.ints(), active.entry.ints(), static_cast<d3*>(active.normal.data()),
                                                                          active.st.doubles());
    if (storage == ContactStorage::LiveRows) {
        if (N) sortActiveEntriesByBlockKernel<<<N, kSortThreads, 0, stream>>>(active.start.ints(), active.entry.ints());
    } else {
        sortActiveEntriesKernel<<<blocksFor(N), kThreadsPerBlock, 0, stream>>>(N, active.start.ints(), active.entry.ints());
    }
    checkLaunch("ncp-gpu: active rows");
}

} // namespace knit::gpu
