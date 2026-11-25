#pragma once
// The exact segment budget over all parts, evaluating as few candidates as possible.
//
// The export keeps the accepted fibres with float(key) < q, q being the first key at which the
// summed segments of the key-sorted accepted fibres would exceed the budget
// (budgetThreshold). A candidate's key is a hash known before any geometry
// (CandidateSource::key), so only the candidates with small keys need evaluating:
//   1. count: every candidate's key, into buckets by the top 16 bits of its float bit pattern
//      (non-negative floats order like their bits; equal keys share a bucket; 1.0f is alone in
//      the top bucket);
//   2. windows: runs of whole buckets over all parts: a pilot of about 2 % of the candidates, then
//      the estimated remainder x 1.1, then at least doubling;
//   3. accumulate: the accepted segments per bucket, summed in bucket order. At the first bucket
//      whose sum would exceed the budget its (key, segments) are sorted and budgetThreshold's scan
//      continues there, which gives the same q as evaluating every candidate. A sum that never
//      exceeds the budget (exact exhaustion included: the condition is strict) evaluates
//      everything: q = 2.
// Output per part: the accepted fibres with key < keyBound (the first bucket not evaluated, 2 when
// complete), in (strand, index) order, with the counters of exactly those candidates. The fibres
// below q are therefore the full run's and the export is byte-identical; the counters cover the
// evaluated subset.
//
// BudgetSchedule is the plan as data: it says which window to evaluate next and collects what the caller
// evaluated. The caller (scene/fibres/fibre_stage.cpp) owns the evaluation, its device support, its
// concurrency, cancellation and the log.
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "fibres/fibre_generator.h"

namespace knit {

constexpr int kKeyBuckets = 0x3F81;  // keys in [0, 1]: buckets 0 .. 0x3F80 (1.0f)
int keyBucket(float key);            // top 16 bits of the bit pattern (key in [0, 1])
float bucketLowerKey(int b);         // smallest float of bucket b; 2.0 for b >= kKeyBuckets

// A part's candidates per key bucket, by CandidateSource::key: what a schedule is planned from.
std::vector<int64_t> candidatesPerBucket(const CandidateSource& source);

// A run of whole buckets that every part evaluates next.
struct KeyWindow {
    int firstBucket = 0, endBucket = 0;
    float keyLo = 0, keyHi = 0;        // bucketLowerKey of both ends: the candidates with key in [keyLo, keyHi)
    std::vector<int64_t> candidates;   // per part, in the window; 0: the part has nothing to evaluate
};

struct BudgetPartResult {
    FibreSet fs;              // accepted fibres with key < keyBound, (strand, index) order
    float keyBound = 0;       // every candidate with key < keyBound is covered (2.0: all)
    int64_t candidates = 0;   // of the source
    int64_t covered = 0;      // candidates with key < keyBound
    int64_t evaluated = 0;    // the same, as counted by the windows
};

// What a schedule decided. Its times and the support it rebuilt are the caller's.
struct BudgetRun {
    float q = 2.0f;               // export threshold: keep float(key) < q
    float keyBound = 2.0f;        // common coverage (every part's keyBound)
    bool complete = false;        // every candidate covered
    long long budget = 0;
    long long segmentsBelowQ = 0; // accepted segments with key < q (<= budget when q < 2)
    long long segmentsCovered = 0;
    int batches = 0;              // windows evaluated
    int crossingBucket = -1;
    int64_t candidates = 0, covered = 0, evaluated = 0;
    std::vector<BudgetPartResult> parts;
};

// The steps of the budget in order: nextWindow, the caller evaluates it, addWindow; until finished, then
// finish. One window is pending at a time, and every part answers it.
class BudgetSchedule {
public:
    // perBucket: per part, its candidatesPerBucket. budget <= 0: one window holds every candidate.
    BudgetSchedule(std::vector<std::vector<int64_t>> perBucket, long long budget);
    bool finished() const;  // the budget was crossed, or every bucket is evaluated
    // Plans the next window from the segments accepted so far and holds it as pending. Throws if the
    // schedule is finished or a window is already pending.
    const KeyWindow& nextWindow();
    // The pending window's fibres, one set per part (an empty set where the window holds none of its
    // candidates). Throws if no window is pending, the set count is not the part count, or a part's
    // candidate count differs from the window's. Sums the window's buckets in key order; at the bucket
    // that crosses the budget, finds q exactly.
    void addWindow(std::vector<FibreSet> fibres);
    std::string windowLine() const;  // the last added window, for the log
    // Once, when finished: every part's fibres in candidate order, what they cover, and q.
    BudgetRun finish();

private:
    // One part's accepted fibres: each window's appended in candidate order.
    struct CollectedFibres {
        FibreSet set;
        std::vector<size_t> windowStart;
    };
    int64_t windowTarget(int64_t remaining) const;
    float crossingKey(int b) const;
    bool accumulateBuckets(int from, int to);

    std::vector<std::vector<int64_t>> perBucket_;  // [part][bucket]
    std::vector<int64_t> allBuckets_;              // summed over the parts
    long long budget_ = 0;
    std::vector<long long> segmentsInBucket_;      // accepted, summed as windows arrive
    long long accumulated_ = 0;                    // accepted segments of the buckets already accumulated
    long long seen_ = 0;                           // accepted segments of everything evaluated
    std::vector<CollectedFibres> collected_;
    std::optional<KeyWindow> pending_;
    int64_t pendingTarget_ = 0, lastTarget_ = 0;   // the candidates the pending and the last window aimed at
    int next_ = 0;                                 // the first bucket not evaluated
    bool crossed_ = false, finished_ = false;
    std::string lastLine_;
    BudgetRun run_;
};

// The key threshold q such that the fibres with key < q have at most `budget` segments in total (2 if
// all fit or budget <= 0). Depends only on the multiset of (key, segments).
float budgetThreshold(std::vector<std::pair<float, int>> keySegments, long long budget);

} // namespace knit
