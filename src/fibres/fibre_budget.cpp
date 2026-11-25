// The budget schedule: windows of key buckets until the budget's threshold is known.
#include "fibres/fibre_budget.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>

namespace knit {

namespace {

constexpr int64_t kSmallestBatch = 4096;
constexpr double kPilotFraction = 0.02;  // the first window, of all candidates

uint32_t keyBits(float k) {
    uint32_t b;
    std::memcpy(&b, &k, 4);
    return b;
}

// A fibre's place in candidate order.
struct FibreOrder {
    int strand = 0, candidate = 0;
    size_t at = 0;
};
bool comesFirst(const FibreOrder& a, const FibreOrder& b) { return a.strand != b.strand ? a.strand < b.strand : a.candidate < b.candidate; }

// The bucket after a window of at least one bucket and about `target` candidates from `first` on.
// Empty buckets cost nothing, so they join the window.
int windowEnd(const std::vector<int64_t>& perBucket, int first, int64_t target) {
    int end = first;
    int64_t inWindow = 0;
    while (end < kKeyBuckets && (inWindow < target || end == first)) inWindow += perBucket[end++];
    while (end < kKeyBuckets && perBucket[end] == 0) ++end;
    return end;
}

// One part's fibres in (strand, index) order: every window is in that order already, so its runs are
// merged, then records and points are copied once into that order.
FibreSet inCandidateOrder(const FibreSet& set, const std::vector<size_t>& windowStart) {
    const std::vector<Fibre>& fibres = set.fibres;
    std::vector<FibreOrder> order;
    order.reserve(fibres.size());
    for (size_t k = 0; k < fibres.size(); ++k) {
        if (fibres[k].candidate < 0) throw std::runtime_error("BudgetSchedule: a fibre without a candidate identity");
        order.push_back({fibres[k].strand, fibres[k].candidate, k});
    }
    for (size_t k = 1; k < windowStart.size(); ++k) {
        const size_t next = k + 1 < windowStart.size() ? windowStart[k + 1] : order.size();
        std::inplace_merge(order.begin(), order.begin() + static_cast<std::ptrdiff_t>(windowStart[k]),
                           order.begin() + static_cast<std::ptrdiff_t>(next), comesFirst);
    }
    FibreSet merged;
    merged.rep = set.rep;
    merged.fibres.reserve(order.size());
    merged.points.reserve(set.points.size());
    for (size_t k = 0; k < order.size(); ++k) {
        if (k && !comesFirst(order[k - 1], order[k])) throw std::runtime_error("BudgetSchedule: a candidate accepted twice");
        const Fibre& f = fibres[order[k].at];
        merged.addFibre(f, set.pointsOf(f));
    }
    return merged;
}

} // namespace

int keyBucket(float key) {
    if (!(key >= 0.0f) || key > 1.0f) throw std::invalid_argument("keyBucket: key outside [0, 1]");
    return static_cast<int>(keyBits(key) >> 16);
}

float bucketLowerKey(int b) {
    if (b >= kKeyBuckets) return 2.0f;
    if (b <= 0) return 0.0f;
    const uint32_t bits = static_cast<uint32_t>(b) << 16;
    float k;
    std::memcpy(&k, &bits, 4);
    return k;
}

std::vector<int64_t> candidatesPerBucket(const CandidateSource& source) {
    std::vector<int64_t> count(kKeyBuckets, 0);
    for (int s = 0; s < source.strands(); ++s)
        for (int i = 0; i < std::max(0, source.count[s]); ++i) ++count[keyBucket(source.key(s, i))];
    return count;
}

BudgetSchedule::BudgetSchedule(std::vector<std::vector<int64_t>> perBucket, long long budget)
    : perBucket_(std::move(perBucket)), allBuckets_(kKeyBuckets, 0), budget_(budget), segmentsInBucket_(kKeyBuckets, 0),
      collected_(perBucket_.size()) {
    run_.budget = budget;
    run_.parts.resize(perBucket_.size());
    for (size_t p = 0; p < perBucket_.size(); ++p) {
        if (perBucket_[p].size() != static_cast<size_t>(kKeyBuckets)) throw std::invalid_argument("BudgetSchedule: a part's buckets are not kKeyBuckets");
        for (int k = 0; k < kKeyBuckets; ++k) {
            allBuckets_[k] += perBucket_[p][k];
            run_.parts[p].candidates += perBucket_[p][k];
        }
        run_.candidates += run_.parts[p].candidates;
    }
}

bool BudgetSchedule::finished() const { return crossed_ || next_ >= kKeyBuckets; }

// How many candidates the next window should hold: all of them without a budget, else a pilot, then
// the estimate of what the budget still needs, then at least double the last window.
int64_t BudgetSchedule::windowTarget(int64_t remaining) const {
    if (budget_ <= 0) return remaining;
    if (run_.batches == 0 || run_.evaluated == 0) return std::max(kSmallestBatch, static_cast<int64_t>(kPilotFraction * double(remaining)));
    const double segmentsPerCandidate = double(seen_) / double(run_.evaluated);
    const double need = segmentsPerCandidate > 0 ? double(budget_ - accumulated_) / segmentsPerCandidate * 1.1 : double(remaining);
    const int64_t target = run_.batches == 1 ? static_cast<int64_t>(need) : std::max(static_cast<int64_t>(need), 2 * lastTarget_);
    return std::max(target, kSmallestBatch);
}

const KeyWindow& BudgetSchedule::nextWindow() {
    if (finished()) throw std::logic_error("BudgetSchedule::nextWindow: the schedule is finished");
    if (pending_) throw std::logic_error("BudgetSchedule::nextWindow: a window is pending");
    int64_t remaining = 0;
    for (int k = next_; k < kKeyBuckets; ++k) remaining += allBuckets_[k];
    pendingTarget_ = windowTarget(remaining);
    KeyWindow w;
    w.firstBucket = next_;
    w.endBucket = windowEnd(allBuckets_, next_, pendingTarget_);
    w.keyLo = bucketLowerKey(w.firstBucket);
    w.keyHi = bucketLowerKey(w.endBucket);
    w.candidates.assign(perBucket_.size(), 0);
    for (size_t p = 0; p < perBucket_.size(); ++p)
        for (int k = w.firstBucket; k < w.endBucket; ++k) w.candidates[p] += perBucket_[p][k];
    pending_ = std::move(w);
    return *pending_;
}

// The key at which bucket b's accepted fibres cross the budget: budgetThreshold's scan over the
// bucket, from what the buckets before it hold.
float BudgetSchedule::crossingKey(int b) const {
    struct Entry {
        float key;
        int segments;
    };
    std::vector<Entry> entries;
    for (const CollectedFibres& p : collected_)
        for (const Fibre& f : p.set.fibres)
            if (keyBucket(float(f.key)) == b) entries.push_back({float(f.key), f.segments()});
    std::sort(entries.begin(), entries.end(), [](const Entry& x, const Entry& y) { return x.key != y.key ? x.key < y.key : x.segments < y.segments; });
    long long accumulated = accumulated_;
    for (const Entry& e : entries) {
        if (accumulated + e.segments > budget_) return e.key;
        accumulated += e.segments;
    }
    return 2.0f;  // not reached: the bucket's sum exceeds the budget
}

// Sums buckets [from, to) in order. True: the budget was crossed (q and the crossing bucket set).
bool BudgetSchedule::accumulateBuckets(int from, int to) {
    for (int b = from; b < to; ++b) {
        if (budget_ > 0 && accumulated_ + segmentsInBucket_[b] > budget_) {
            run_.q = crossingKey(b);
            run_.crossingBucket = b;
            return true;
        }
        accumulated_ += segmentsInBucket_[b];
    }
    return false;
}

void BudgetSchedule::addWindow(std::vector<FibreSet> fibres) {
    if (!pending_) throw std::logic_error("BudgetSchedule::addWindow: no window is pending");
    if (fibres.size() != collected_.size()) throw std::invalid_argument("BudgetSchedule::addWindow: one fibre set per part");
    const KeyWindow window = std::move(*pending_);
    pending_.reset();
    int64_t evaluated = 0;
    for (size_t p = 0; p < collected_.size(); ++p) {
        const FibreSet& batch = fibres[p];
        if (batch.rep.candidates != window.candidates[p])
            throw std::runtime_error("BudgetSchedule: a key window evaluated " + std::to_string(batch.rep.candidates) +
                                     " candidates, the source has " + std::to_string(window.candidates[p]));
        CollectedFibres& part = collected_[p];
        part.set.rep.addCounters(batch.rep);
        part.windowStart.push_back(part.set.fibres.size());
        part.set.appendSet(batch);
        for (const Fibre& f : batch.fibres) {
            segmentsInBucket_[keyBucket(float(f.key))] += f.segments();
            seen_ += f.segments();
        }
        run_.parts[p].evaluated += window.candidates[p];
        evaluated += window.candidates[p];
    }
    run_.evaluated += evaluated;
    ++run_.batches;
    lastTarget_ = pendingTarget_;
    crossed_ = accumulateBuckets(window.firstBucket, window.endBucket);
    char b[256];
    std::snprintf(b, sizeof b, "budget batch %d: buckets [%d, %d) keys [%.6f, %.6f): %lld candidates; S %lld of %lld%s", run_.batches,
                  window.firstBucket, window.endBucket, double(window.keyLo), double(window.keyHi), static_cast<long long>(evaluated),
                  accumulated_, budget_, crossed_ ? " -> crossed" : "");
    lastLine_ = b;
    next_ = window.endBucket;
}

std::string BudgetSchedule::windowLine() const { return lastLine_; }

BudgetRun BudgetSchedule::finish() {
    if (finished_) throw std::logic_error("BudgetSchedule::finish: already finished");
    if (pending_ || !finished()) throw std::logic_error("BudgetSchedule::finish: windows remain");
    finished_ = true;
    run_.complete = next_ >= kKeyBuckets;
    run_.keyBound = bucketLowerKey(next_);
    if (!crossed_) run_.q = 2.0f;
    for (size_t p = 0; p < collected_.size(); ++p) {
        BudgetPartResult& r = run_.parts[p];
        r.keyBound = run_.keyBound;
        for (int k = 0; k < next_; ++k) r.covered += perBucket_[p][k];
        run_.covered += r.covered;
        r.fs = inCandidateOrder(collected_[p].set, collected_[p].windowStart);
        resumFibreLengths(r.fs);
        for (const Fibre& f : r.fs.fibres) {
            run_.segmentsCovered += f.segments();
            if (float(f.key) < run_.q) run_.segmentsBelowQ += f.segments();
        }
    }
    return std::move(run_);
}

float budgetThreshold(std::vector<std::pair<float, int>> keySegments, long long budget) {
    std::sort(keySegments.begin(), keySegments.end());
    long long total = 0;
    if (budget > 0)
        for (const auto& [key, segments] : keySegments) {
            if (total + segments > budget) return key;  // every fibre with this key or above is dropped
            total += segments;
        }
    return 2.0f;
}

} // namespace knit
