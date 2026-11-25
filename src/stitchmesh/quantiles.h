#pragma once
// The five-number summary of a sample: what a report says about a distribution (edge lengths, stitch widths, angles)
// instead of a mean. Shared by knitting (layout metrics) and stitchmesh (coarse-mesh quality), which do not see each other.
#include <vector>

namespace knit {

struct Quantiles {
    double min = 0, p5 = 0, p50 = 0, p95 = 0, max = 0;
    int count = 0;  // 0: no sample, every value 0
};

// The quantiles of `values` (a copy is sorted; the input stays). Quantile p sits at rank p (n - 1), linearly
// interpolated between the two neighbouring sorted values.
Quantiles quantiles(std::vector<double> values);

} // namespace knit
