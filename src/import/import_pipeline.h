#pragma once
// One preparation from start to end, in the order it happens: the file split by the user's selection, the surface of the
// knitted nodes (the file's own, or rebuilt as one closed surface by the import's own stages), the plain nodes as the
// file has them, both placed where the whole object stands, then one stored entry under the import's identity. Whichever way
// the surface was made, what comes out is an ordinary prepared import.
#include <cstdint>
#include <filesystem>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "import/prepare_import.h"

namespace knit {

// Prepares the request and stores it in the import cache. changes: what the preparation did, one line each. false with a
// message for the user. Cancellable.
[[nodiscard]] bool prepareAndStore(const ImportRequest& r, std::vector<std::string>& changes, std::string* err,
                                   const ImportProgress& progress);

// The rebuild ends here: `folder` holds the rebuilt surface (wrapped.obj in mm, Y up, with its charts as uv, and its MTL and
// texture), which goes back through the ordinary reader and repairs into `out`; `o` gives the source file's unit and up axis,
// `summary` what the rebuild did.
bool readRebuiltSurface(const std::filesystem::path& folder, const ImportOptions& o, const std::vector<std::string>& summary, PreparedMesh& out,
                        std::string* err, const ImportProgress& progress);

// The file's triangles by node, each side with only the vertices it uses. file.nodes must be finished (finishNodes).
struct NodeSplit {
    ImportedMesh knitted, plain;
};
NodeSplit splitByNodes(const ImportedMesh& file, const std::set<std::string>& plainNodes);

} // namespace knit
