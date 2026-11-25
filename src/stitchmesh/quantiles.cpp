#include "stitchmesh/quantiles.h"

#include <algorithm>
#include <cmath>

namespace knit {

namespace {

double atRank(const std::vector<double>& sorted, double p) {
    const double rank = p * (sorted.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(rank));
    const size_t hi = std::min(lo + 1, sorted.size() - 1);
    const double t = rank - lo;
    return sorted[lo] * (1 - t) + sorted[hi] * t;
}

} // namespace

Quantiles quantiles(std::vector<double> values) {
    Quantiles q;
    if (values.empty()) return q;
    std::sort(values.begin(), values.end());
    q.count = static_cast<int>(values.size());
    q.min = values.front();
    q.p5 = atRank(values, 0.05);
    q.p50 = atRank(values, 0.5);
    q.p95 = atRank(values, 0.95);
    q.max = values.back();
    return q;
}

} // namespace knit
