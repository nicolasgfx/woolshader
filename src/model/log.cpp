// Log text and elapsed time for the stages' logs (log.h). Profiling only.
#include "model/log.h"

#include <cstdarg>
#include <cstdio>
#include <vector>

namespace knit {

void appendFormat(std::string& log, const char* format, ...) {
    va_list values, again;
    va_start(values, format);
    va_copy(again, values);
    const int length = std::vsnprintf(nullptr, 0, format, values);
    va_end(values);
    if (length > 0) {
        std::vector<char> text(static_cast<size_t>(length) + 1);
        std::vsnprintf(text.data(), text.size(), format, again);
        log.append(text.data(), static_cast<size_t>(length));
    }
    va_end(again);
}

double msSince(std::chrono::steady_clock::time_point start) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
}

double secondsSince(std::chrono::steady_clock::time_point start) { return msSince(start) / 1000; }

} // namespace knit
