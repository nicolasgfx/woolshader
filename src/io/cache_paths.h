#pragma once
// Where the caches live: out/cache/<name> beside the checkout, and the folder names of the import store
// (out/cache/imports/<identity>), which saved settings and the scene registry refer to.
#include <cstdint>
#include <string>

namespace knit {

// out/cache/<name>: relative when the program runs in the checkout (a scenes/ folder here), else in the
// source tree (sourceTreeDir), else relative.
std::string outputCacheDir(const std::string& name);

std::string importCacheRoot();  // out/cache/imports (next to the scene cache), or $KNIT_IMPORT_CACHE
std::string markingsCacheRoot();  // out/cache/markings (the Mark step's automatic saves), or $KNIT_MARKINGS_CACHE
std::string importCacheDir(const std::string& root, uint64_t identity);
std::string identityHex(uint64_t identity);  // 16 hex digits, the folder name of an import
bool parseIdentityHex(const std::string& s, uint64_t& identity);

} // namespace knit
