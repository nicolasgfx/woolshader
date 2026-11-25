#pragma once
// Environment variables through _dupenv_s on MSVC (std::getenv is C4996 there).
#include <cstdlib>
#include <optional>
#include <string>

namespace knit {

// The variable's value, nullopt if it is unset.
inline std::optional<std::string> envVar(const char* name) {
#ifdef _MSC_VER
    char* v = nullptr;
    size_t n = 0;
    if (_dupenv_s(&v, &n, name) != 0 || !v) return std::nullopt;
    std::string r = v;
    std::free(v);
    return r;
#else
    const char* v = std::getenv(name);
    if (!v) return std::nullopt;
    return std::string(v);
#endif
}

}  // namespace knit
