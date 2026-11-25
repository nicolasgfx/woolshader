#pragma once
// What every stage of the stitch-mesh layout returns beside its product: identities, counts, time and a verdict.
// The report is the evidence a failed or incomplete artifact is inspected by; UNRESOLVED never counts as PASS.
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "model/status.h"

namespace knit {

struct StageReport {
    std::string stage;                 // "remesh"
    int algorithmVersion = 1;          // bumped when the stage's rules change (a cache key input later)
    uint64_t inputDigest = 0, outputDigest = 0;
    uint64_t seed = 0;
    double gaugeW = 0, gaugeH = 0;     // the effective gauge [mm]
    double elapsedMs = 0;
    double peakHostMb = 0;             // filled by the caller that measures the process (the lab); 0: not measured
    std::string settings;              // "key = value" lines, in order
    std::vector<std::pair<std::string, double>> counts;
    Status status = Status::Unresolved;
    std::string detail;                // why FAIL or UNRESOLVED; what was kept on PASS
    std::string text() const;          // one block, the same words every writer of the report uses
};

// A decimal digest for the report lines.
std::string hexDigest(uint64_t digest);

} // namespace knit
