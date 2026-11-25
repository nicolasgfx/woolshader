#include "stitchmesh/stage_report.h"

#include <cstdio>

namespace knit {

std::string hexDigest(uint64_t digest) {
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(digest));
    return b;
}

std::string StageReport::text() const {
    std::string t = "stage " + stage + " (algorithm " + std::to_string(algorithmVersion) + "): " + toString(status) + "\n";
    if (!detail.empty()) t += "  " + detail + "\n";
    t += "  input " + hexDigest(inputDigest) + ", output " + hexDigest(outputDigest) + ", seed " + std::to_string(seed) + "\n";
    char b[160];
    std::snprintf(b, sizeof b, "  gauge %.3f x %.3f mm, %.0f ms", gaugeW, gaugeH, elapsedMs);
    t += b;
    if (peakHostMb > 0) {
        std::snprintf(b, sizeof b, ", peak host %.0f MB", peakHostMb);
        t += b;
    }
    t += "\n";
    for (const std::pair<std::string, double>& c : counts) {
        if (c.second == static_cast<long long>(c.second)) std::snprintf(b, sizeof b, "  %s = %lld\n", c.first.c_str(), static_cast<long long>(c.second));
        else std::snprintf(b, sizeof b, "  %s = %.4f\n", c.first.c_str(), c.second);
        t += b;
    }
    if (!settings.empty()) t += "  settings:\n" + settings;
    return t;
}

} // namespace knit
