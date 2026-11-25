// ncp-gpu, ContactStorage::LiveRows (ncp.h): the contact
// rows held as two lists. The row list holds every row of the last refresh whose segments are within the
// near cutoff (or that carries a positive multiplier): its pair and type, 14 bytes with its live flag and
// slot. The live list holds the rows with a term in the merit or a positive multiplier: their row, multiplier,
// values and terms, and it is the contact list every term, multiplier and active-list kernel runs over.
//
// The live list is renumbered at every evaluation of the current state, in row order (an exclusive scan of
// the live flags), and the multipliers move along with their rows. Between two evaluations a row cannot
// become live: after the multiplier update lambda = max(0, rho * value) stays 0 for a row without a term
// (its value is <= 0), and a larger rho lowers nothing to a term. A trial evaluates every row (the rows that
// are not live with multiplier 0); the accepted trial's evaluation at the new state gives the rows that
// gained a term their place. A refresh carries the positive multipliers by row key, as CompactRows does.
//
// Capacities are fixed per call from the starting state (ncp_setup.cu): an overflow ends the solve with an
// explicit UNRESOLVED status at the last accepted state, never a truncation and never a regrow.
#include "gpu/solver/ncp_device.cuh"

namespace knit::gpu {

namespace {

__device__ __forceinline__ double penaltyOf(const NcpControl* ctl, EvaluationPenalty penalty) {
    return penalty == EvaluationPenalty::LastTrialRho ? ctl->lastTrialRho : ctl->rho;
}

// --- the evaluation of the current state
// Pass 1, every row: whether it is live at the current multipliers (a term in the merit, a target term, or a
// positive multiplier), as 1 / 0 in `flags` (the row slots, scanned next; 0 past the rows up to the capacity).
// The violation maxima and counts of evalContactRowsKernel: a row that is not live has none (its value is
// <= 0 and it has no target term), so they are complete even when the live list overflows.
__global__ void liveFlagsKernel(NcpRows R, NcpParams pr, const d3* points, const NcpControl* ctl, EvaluationPenalty penalty, long long rowCapacity, int* flags,
                                EvalReductions* sum) {
    const double rho = penaltyOf(ctl, penalty);
    const int nRows = *R.nRows;
    unsigned long long live = 0;
    for (long long rr = loopStart(); rr <= rowCapacity; rr += loopStride()) {
        const int r = static_cast<int>(rr);
        if (rr >= nRows) {
            flags[r] = 0;
            continue;
        }
        d3 pa, pb;
        double s, t;
        const double d = contactDistance(points, R.cPair[r].x, R.cPair[r].y, R.cType[r], pa, pb, s, t);
        const double val = (pr.hardC - d) / pr.D;
        const double lambda = rowMultiplier(R, r);  // reads this row's slot before the flag replaces it
        double wres = 0, wgn = 0, m = 0;
        hardTerm(val, lambda, rho, wres, wgn, m);
        const bool target = pr.targets && d < pr.thr;
        const bool isLive = wgn != 0 || lambda > 0 || target;
        flags[r] = isLive ? 1 : 0;
        live += isLive ? 1 : 0;
        if (target) {
            const double tval = (pr.thr - d) / pr.D;
            atomicMaxPositive(&sum->maxViol[FContactT], tval);
            if (tval > kViolationTolerance) atomicAdd(&sum->nViol[FContactT], 1ull);
        }
        atomicMaxPositive(&sum->maxViol[FContact], val);
        if (val > kViolationTolerance) atomicAdd(&sum->nViol[FContact], 1ull);
    }
    if (live) atomicAdd(&sum->liveRows, live);
}

// The scanned flags' total is the new live list's length; beyond the capacity it is an overflow (the
// live list is left as it was and the controller ends the solve, afterAssembleKernel / innerEndKernel).
__global__ void liveCountKernel(const int* slots, long long rowCapacity, long long liveCapacity, NcpCounts* counts) {
    const long long total = slots[rowCapacity];
    counts->previousLive = counts->liveRows;
    if (counts->overflow) return;
    if (total > liveCapacity) {
        counts->overflow = kCapLiveRows;
        counts->overflowNeeded = total;
        counts->overflowCapacity = liveCapacity;
        return;
    }
    counts->liveRows = static_cast<int>(total);
}

// Pass 2, every row: a live row's place in the new numbering (its slot), its row there, multiplier 0 until
// the move; the live flags of the new numbering.
__global__ void liveSlotsKernel(const int* nRowsOf, const int* slots, const NcpCounts* counts, unsigned char* live, int* rowNext, double* lambdaNext) {
    if (counts->overflow) return;
    const int nRows = *nRowsOf;
    for (long long rr = loopStart(); rr < nRows; rr += loopStride()) {
        const int r = static_cast<int>(rr);
        const bool isLive = slots[r + 1] != slots[r];
        live[r] = isLive ? 1 : 0;
        if (!isLive) continue;
        rowNext[slots[r]] = r;
        lambdaNext[slots[r]] = 0;
    }
}

// The positive multipliers of the previous numbering to their rows' new slots (such a row is live: its
// multiplier makes it so).
__global__ void liveMoveKernel(const int* row, const double* lambda, const int* slots, const unsigned char* live, const NcpCounts* counts, double* lambdaNext) {
    if (counts->overflow) return;
    const int nOld = counts->previousLive;
    for (long long kk = loopStart(); kk < nOld; kk += loopStride()) {
        const int k = static_cast<int>(kk);
        const int r = row[k];
        if (r < 0 || !(lambda[k] > 0) || !live[r]) continue;
        lambdaNext[slots[r]] = lambda[k];
    }
}

// The new numbering becomes the live list.
__global__ void liveAdoptKernel(const NcpCounts* counts, const int* rowNext, const double* lambdaNext, int* row, double* lambda) {
    if (counts->overflow) return;
    const int n = counts->liveRows;
    for (long long kk = loopStart(); kk < n; kk += loopStride()) {
        const int k = static_cast<int>(kk);
        row[k] = rowNext[k];
        lambda[k] = lambdaNext[k];
    }
}

// Pass 3, every live row: its values and terms as evalContactRowsKernel computes them (the same operations
// in the same order), stored at its slot.
__global__ void liveValuesKernel(NcpRows R, NcpParams pr, const d3* points, RowValuePair pair, const NcpControl* ctl, EvaluationPenalty penalty, const NcpCounts* counts) {
    if (counts->overflow) return;
    const NcpRowValues& V = rowValues(pair, RowSet::Current);
    const double rho = penaltyOf(ctl, penalty);
    const int n = counts->liveRows;
    for (long long kk = loopStart(); kk < n; kk += loopStride()) {
        const int k = static_cast<int>(kk);
        const int r = R.cRow[k];
        d3 pa, pb;
        double s, t;
        const double d = contactDistance(points, R.cPair[r].x, R.cPair[r].y, R.cType[r], pa, pb, s, t);
        const double val = (pr.hardC - d) / pr.D;
        double wres = 0, wgn = 0, m = 0, tval = 0;
        hardTerm(val, R.cLam[k], rho, wres, wgn, m);
        V.cGradOn[k] = wgn != 0 || (pr.targets && d < pr.thr) ? 1 : 0;
        if (pr.targets && (pr.thr - d) / pr.D > 0) wgn += pr.wT;
        if (pr.targets && d < pr.thr) {
            tval = (pr.thr - d) / pr.D;
            wres += pr.wT * tval;
        }
        V.cDist[k] = d;  // val and tval follow from it (contactHardValue, contactTargetValue)
        V.cWres[k] = wres;
        V.cWgn[k] = wgn;
    }
}

// --- the refresh's carry-over
// The live rows' keys and positive multipliers (the largest key and 0 elsewhere, up to the capacity, so the
// sort over the capacity puts them last); their number into NcpCounts::carried. Runs before the rows are
// rebuilt: the pairs are still those of the live rows.
__global__ void liveCarryKeysKernel(NcpRows R, long long liveCapacity, NcpCounts* counts, unsigned long long* keys, double* lambda) {
    const int nLive = counts->liveRows;
    unsigned long long carried = 0;
    for (long long kk = loopStart(); kk < liveCapacity; kk += loopStride()) {
        const int k = static_cast<int>(kk);
        if (k < nLive && R.cLam[k] > 0) {
            const int r = R.cRow[k];
            keys[k] = contactKey(R.cPair[r].x, R.cPair[r].y, R.cType[r]);
            lambda[k] = R.cLam[k];
            ++carried;
        } else {
            keys[k] = ~0ull;
            lambda[k] = 0;
        }
    }
    if (carried) atomicAdd(&counts->carried, static_cast<int>(carried));
}

// The sorted table becomes the live list: its multipliers in place, every row unknown until a new row
// claims it (liveCarryRowsKernel). Skipped when the refresh overflowed (the previous rows stay).
__global__ void liveCarryResetKernel(CarryTable table, const NcpCounts* counts, int* row, double* lambda) {
    if (counts->overflow) return;
    const int n = counts->carried;
    for (long long kk = loopStart(); kk < n; kk += loopStride()) {
        const int k = static_cast<int>(kk);
        row[k] = -1;
        lambda[k] = table.lambda[k];
    }
}

// Every new row: live with its carried multiplier when its key is in the table (keys are unique: a table
// entry is claimed at most once), otherwise not live.
__global__ void liveCarryRowsKernel(CarryTable table, const int2* pairs, const unsigned char* type, const NcpCounts* counts, unsigned char* live, int* slots, int* row) {
    if (counts->overflow) return;
    const int nRows = counts->contactRows, nTable = counts->carried;
    for (long long rr = loopStart(); rr < nRows; rr += loopStride()) {
        const int r = static_cast<int>(rr);
        const unsigned long long k = contactKey(pairs[r].x, pairs[r].y, type[r]);
        int lo = 0, hi = nTable;
        while (lo < hi) {
            const int mid = (lo + hi) >> 1;
            if (table.keys[mid] < k) lo = mid + 1;
            else hi = mid;
        }
        const bool found = lo < nTable && table.keys[lo] == k;
        live[r] = found ? 1 : 0;
        if (!found) continue;
        slots[r] = lo;
        row[lo] = r;
    }
}

// The carried table is the live list until the evaluation that follows the refresh renumbers it.
__global__ void liveCarryCountKernel(NcpCounts* counts) {
    if (counts->overflow) return;
    counts->liveRows = counts->carried;
}

}  // namespace

// Every buffer that follows the two capacities, sized once per call (no regrow). The active lists run over
// the general rows and the live list; the active entries' memory also holds the carry-over's sort.
void NcpDevice::State::allocateLiveCapacity(long long rows, long long live) {
    capacity.contactRows = rows;
    capacity.liveRows = live;
    const long long R = nG + live;
    checkIndexLimit(rows + 1, "contact rows");
    checkIndexLimit(R, "rows");
    checkIndexLimit(4LL * nG + kContactRowEntries * live, "row entries");
    if (R >= (1LL << 27)) throw DeviceError("ncp-gpu: more than 2^27 rows (active list entries hold slot << 4)");
    capacity.activeEntries = 4LL * nG + kContactRowEntries * live;
    contacts.pair.alloc(size_t(rows) * sizeof(int2) + 8);
    contacts.type.alloc(size_t(rows) + 8);
    contacts.lambda.alloc(size_t(live) * 8 + 8);
    liveList.slot.alloc(size_t(rows + 1) * 4 + 8);
    liveList.live.alloc(size_t(rows) + 8);
    liveList.row.alloc(size_t(live) * 4 + 8);
    allocateValueSets();
    active.slotOf.alloc(size_t(R + 1) * 4 + 8);
    active.rows.alloc(size_t(R) * 4 + 8);
    active.count.alloc(size_t(N + 1) * 4);
    active.start.alloc(size_t(N + 1) * 4);
    active.cursor.alloc(size_t(N + 1) * 4);
    active.entry.alloc(std::max(size_t(capacity.activeEntries) * 4, size_t(live) * 32) + 8);  // the carry-over: 4 arrays of live x 8 bytes
    active.normal.alloc(size_t(R) * sizeof(d3) + 8);
    active.st.alloc(size_t(R) * 16 + 8);
    system.rowT.alloc(size_t(R) * 8 + 8);  // also the next multipliers while the live list is renumbered
    system.rowS.alloc(size_t(nG) * sizeof(d3) + 8);
    const long long nCells = grid.n[0] * grid.n[1] * grid.n[2];
    ScratchPlan plan;
    plan.scanItems = std::max({rows + 1, R + 1, nCells + 1, static_cast<long long>(N) + 1, static_cast<long long>(nSamples) + 1});
    plan.keySortItems = capacity.registrations;
    plan.pairSortItems = live;
    prims.reserveScratch(plan);
}

CarryTable NcpDevice::State::enqueueLiveCarryTable() {
    const CarryBuffers carry = carryBuffers();
    liveCarryKeysKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(rows(), carry.n, counts(), carry.keys, carry.lambda);
    const bool sortedInAlt = prims.sortPairsAlternatingOnStream({carry.keys, carry.keysAlt, carry.lambda, carry.lambdaAlt, carry.n}, stream);
    checkLaunch("ncp-gpu: live carry table");
    return {sortedInAlt ? carry.keysAlt : carry.keys, sortedInAlt ? carry.lambdaAlt : carry.lambda};
}

void NcpDevice::State::enqueueLiveCarryOver(const CarryTable& table) {
    NcpCounts* cnt = counts();
    liveCarryResetKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(table, cnt, liveList.row.ints(), contacts.lambda.doubles());
    liveCarryRowsKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(table, static_cast<const int2*>(contacts.pair.data()), contacts.type.chars(), cnt,
                                                                      liveList.live.chars(), liveList.slot.ints(), liveList.row.ints());
    liveCarryCountKernel<<<1, 1, 0, stream>>>(cnt);
    checkLaunch("ncp-gpu: live carry-over");
}

void NcpDevice::State::enqueueLiveEvaluation(EvaluationPenalty penalty, EvalReductions* sum) {
    NcpCounts* cnt = counts();
    const NcpRows R = rows();
    const long long rowCapacity = capacity.contactRows;
    int* slots = liveList.slot.ints();
    // The renumbering's scratch lives in buffers that are free during an evaluation: the next rows in the
    // active slots' scan (used only while the active lists are built), the next multipliers in the product's
    // row scalars (ncp_pcg.cu writes them before it reads them).
    int* rowNext = active.slotOf.ints();
    double* lambdaNext = system.rowT.doubles();
    liveFlagsKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(R, pr, pointsAtP(), control(), penalty, rowCapacity, slots, sum);
    prims.countsToOffsetsOnStream(slots, slots, rowCapacity, stream);  // in place: flags -> slots
    liveCountKernel<<<1, 1, 0, stream>>>(slots, rowCapacity, capacity.liveRows, cnt);
    liveSlotsKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(R.nRows, slots, cnt, liveList.live.chars(), rowNext, lambdaNext);
    liveMoveKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(liveList.row.ints(), contacts.lambda.doubles(), slots, liveList.live.chars(), cnt, lambdaNext);
    liveAdoptKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(cnt, rowNext, lambdaNext, liveList.row.ints(), contacts.lambda.doubles());
    liveValuesKernel<<<kLoopBlocks, kThreadsPerBlock, 0, stream>>>(R, pr, pointsAtP(), valuePair(), control(), penalty, cnt);
    checkLaunch("ncp-gpu: live evaluation");
}

} // namespace knit::gpu
