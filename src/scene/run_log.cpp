// StageClock (wall, CPU and validate time per stage) and the process memory probes (run_log.h). Profiling only.
#include "scene/run_log.h"

#include "model/log.h"

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
#else
#include <ctime>
#endif


namespace knit {

namespace {

// Process CPU time (user + kernel, all threads) [ms].
double processCpuMs() {
#ifdef _WIN32
    FILETIME created, ended, kernel, user;
    if (!GetProcessTimes(GetCurrentProcess(), &created, &ended, &kernel, &user)) return 0;
    const auto ms = [](const FILETIME& f) { return ((static_cast<unsigned long long>(f.dwHighDateTime) << 32) | f.dwLowDateTime) / 1e4; };
    return ms(kernel) + ms(user);
#else
    return 1e3 * std::clock() / CLOCKS_PER_SEC;
#endif
}

} // namespace

StageClock::StageClock() : start_(std::chrono::steady_clock::now()), cpuStartMs_(processCpuMs()), validateStart_(validateStats()) {}

double StageClock::wallMs() const { return msSince(start_); }

std::string StageClock::line(const char* label) const {
    const double wall = wallMs(), cpu = processCpuMs() - cpuStartMs_;
    const ValidateStats v = validateStats().since(validateStart_);
    std::string text;
    appendFormat(text, "%s %.0f ms (threads busy %.1f)", label, wall, wall > 0 ? cpu / wall : 0.0);
    if (v.calls)
        appendFormat(text, "; validate %lld calls %.0f ms (curvature %.0f, sampling %.0f, signature %.0f, clearance %.0f), other %.0f ms",
                     v.calls, v.ms, v.msCurvature, v.msSampling, v.msSignature, v.msClearance, wall - v.ms);
    return text;
}

double processPrivateMB() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX counters;
    if (GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&counters), sizeof counters))
        return counters.PrivateUsage / 1048576.0;
#endif
    return 0;
}

double peakWorkingSetMB() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS counters;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof counters)) return counters.PeakWorkingSetSize / 1048576.0;
#endif
    return 0;
}

} // namespace knit
