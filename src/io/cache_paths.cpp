// Where the caches live: out/cache/<name>, and the import store's folder names (cache_paths.h).
#include "io/cache_paths.h"

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <optional>

#include "io/files.h"
#include "model/env.h"

namespace knit {

std::string outputCacheDir(const std::string& name) {
    std::error_code ec;
    if (std::filesystem::is_directory("scenes", ec)) return "out/cache/" + name;
    const std::string tree = sourceTreeDir();
    return tree.empty() ? "out/cache/" + name : tree + "/out/cache/" + name;
}

std::string importCacheRoot() {
    // KNIT_IMPORT_CACHE: another root (tests, portable installs).
    const std::optional<std::string> env = envVar("KNIT_IMPORT_CACHE");
    if (env && !env->empty()) return *env;
    return outputCacheDir("imports");
}

std::string markingsCacheRoot() {
    // KNIT_MARKINGS_CACHE: another root (tests that run side by side on the same mesh).
    const std::optional<std::string> env = envVar("KNIT_MARKINGS_CACHE");
    if (env && !env->empty()) return *env;
    return outputCacheDir("markings");
}

std::string identityHex(uint64_t identity) {
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(identity));
    return b;
}

bool parseIdentityHex(const std::string& s, uint64_t& identity) {
    if (s.size() != 16) return false;
    char* end = nullptr;
    identity = std::strtoull(s.c_str(), &end, 16);
    return end && *end == 0;
}

std::string importCacheDir(const std::string& root, uint64_t identity) { return root + "/" + identityHex(identity); }

} // namespace knit
