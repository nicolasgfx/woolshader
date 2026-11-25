#pragma once
// Log text and elapsed time, for the stages' logs and timings. Profiling only: nothing here changes a result.
#include <chrono>
#include <string>

namespace knit {

// Appends printf-formatted text, of any length.
void appendFormat(std::string& log, const char* format, ...);
double msSince(std::chrono::steady_clock::time_point start);
double secondsSince(std::chrono::steady_clock::time_point start);

} // namespace knit
