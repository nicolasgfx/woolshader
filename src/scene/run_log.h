#pragma once
// Stage timing and memory probes. Profiling only: nothing here changes a result. The log text and
// elapsed-time helpers are model/log.h.
#include <chrono>
#include <string>

#include "knitting/checks/checks.h"

namespace knit {

// One stage: wall time, process CPU time and validate() time, so a complete run shows where its time
// goes and how many threads each stage keeps busy. Started at construction, read with line().
class StageClock {
public:
    StageClock();
    double wallMs() const;
    // "<label> <wall> ms (threads busy <cpu / wall>)[; validate <n> calls <ms> ms (curvature, sampling,
    // signature, clearance), other <wall - validate> ms]"; the validate part only if it was called.
    std::string line(const char* label) const;

private:
    std::chrono::steady_clock::time_point start_;
    double cpuStartMs_;
    ValidateStats validateStart_;
};

double processPrivateMB();  // current private bytes; 0 where unknown
double peakWorkingSetMB();  // 0 where unknown

} // namespace knit
