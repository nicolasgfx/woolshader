// The architecture's tests (knit_tests_architecture): what the compiler cannot hold, checked on the
// sources and on build/components.json. Run one case: knit_tests_architecture <name>; list: --list.
#include <algorithm>
#include <cctype>
#include <atomic>
#include <chrono>
#include <cstring>
#include <mutex>
#include <thread>
#include <functional>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>

#include "harness.h"

using namespace knit;

// ---------------------------------------------------------- architecture ----

// STYLE.md (no new templates) pins the number of template
// declarations and generic lambdas in the production sources (src, tools, third_party;
// not the shim folders, which stand in for third-party template APIs): it may only fall. A change that removes one lowers the pinned number in the same commit; a new one
// fails here and prints where it is.
namespace {

constexpr int kPinnedTemplateDeclarations = 0;
constexpr int kPinnedGenericLambdas = 0;

struct PatternCount {
    int count = 0;
    std::vector<std::string> places;  // "file:line" of each match
};

// Lines of the production sources matching `pattern`; `cheapFilter` is a substring every match
// contains (skips the regex on most lines).
PatternCount countInProductionSources(const std::regex& pattern, const char* cheapFilter) {
    PatternCount found;
    for (const char* dir : {"src", "tools", "third_party"}) {
        if (!std::filesystem::exists(std::filesystem::path(KNIT_SOURCE_DIR) / dir)) continue;  // tools/ is not published
        for (const auto& entry : std::filesystem::recursive_directory_iterator(std::filesystem::path(KNIT_SOURCE_DIR) / dir)) {
            const std::string ext = entry.path().extension().string();
            if (ext != ".h" && ext != ".cpp" && ext != ".cu" && ext != ".cuh") continue;
            if (entry.path().generic_string().find("/shim/") != std::string::npos) continue;
            if (entry.path().generic_string().find("/third_party/slim/src/") != std::string::npos) continue;  // vendored MPL-2.0 code
            std::ifstream in(entry.path());
            std::string line;
            for (int number = 1; std::getline(in, line); ++number) {
                if (line.find(cheapFilter) == std::string::npos || !std::regex_search(line, pattern)) continue;
                ++found.count;
                found.places.push_back(entry.path().lexically_relative(KNIT_SOURCE_DIR).generic_string() + ":" +
                                       std::to_string(number));
            }
        }
    }
    return found;
}

void printPlaces(const PatternCount& found) {
    for (const std::string& place : found.places) std::printf("    %s\n", place.c_str());
}

} // namespace

// The declared graph (build/components.json, written by cmake/components.cmake; src/ARCHITECTURE.md): a file
// includes only its own folder and the folders its CMakeLists.txt declares as dependencies. The compiler
// holds what is on the include path, but a PUBLIC dependency hands its own dependencies' headers on, so
// only this test holds the direct includes. No include is spelled relative ("../") or absolute either:
// such a spelling would get around the forwarding headers.
TEST(architecture_module_dependencies) {
    namespace fs = std::filesystem;
    std::ifstream json(KNIT_COMPONENTS_JSON);
    const std::string text((std::istreambuf_iterator<char>(json)), std::istreambuf_iterator<char>());
    std::map<std::string, std::set<std::string>> allowed;  // folder -> itself and its declared dependencies
    const std::regex entryRe(R"re(\{"name": "(\w+)", "target": "\w+",\s*"public_deps": \[([^\]]*)\],\s*"private_deps": \[([^\]]*)\])re");
    const std::regex quotedRe(R"re("(\w+)")re");
    const std::sregex_iterator none;
    for (std::sregex_iterator it(text.begin(), text.end(), entryRe); it != none; ++it) {
        std::set<std::string>& may = allowed[(*it)[1].str()];
        may.insert((*it)[1].str());
        for (int group : {2, 3}) {
            const std::string deps = (*it)[group].str();
            for (std::sregex_iterator q(deps.begin(), deps.end(), quotedRe); q != none; ++q) may.insert((*q)[1].str());
        }
    }
    const std::regex include(R"re(^\s*#\s*include\s+"([^"]+)")re");
    int files = 0;
    std::vector<std::string> wrong;
    for (const auto& [folder, may] : allowed)
        for (const auto& entry : fs::recursive_directory_iterator(fs::path(KNIT_SOURCE_DIR) / "src" / folder)) {
            if (!entry.is_regular_file()) continue;
            ++files;
            std::ifstream in(entry.path());
            std::string line;
            std::smatch m;
            while (std::getline(in, line)) {
                if (!std::regex_search(line, m, include)) continue;
                const std::string spelled = m[1].str();
                const std::string where = folder + "/" + entry.path().filename().string() + ": \"" + spelled + "\"";
                if (spelled.rfind("..", 0) == 0 || spelled.rfind("/", 0) == 0 || spelled.find(':') != std::string::npos) {
                    wrong.push_back(where + " is relative or absolute");
                    continue;
                }
                const std::string target = spelled.substr(0, spelled.find('/'));
                if (spelled.find('/') != std::string::npos && allowed.count(target) && !may.count(target))
                    wrong.push_back(where + ": " + target + " is not a declared dependency");
            }
        }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(allowed.size() >= 11 && files > 250 && wrong.empty());
}

// The old layout's folders and targets do not come back in the sources, the docs beside them or
// the build, and the production scans cover src (a scan of an empty tree would pass). That goes for
// src/knitting's flat paths too, in comments as well: only numerics.h sits at the folder's root.
TEST(architecture_no_stale_paths) {
    namespace fs = std::filesystem;
    const std::regex stale(R"re((^|[^\w/.-])(geom/(src|include)|gpu/(src|include)|knitgeom|knitscene|knitgpu|third_party/remesh)\b)re");
    const std::regex flatKnitting(R"re(knitting/(?!numerics\.h|ARCHITECTURE\.md)\w+\.\w+)re");
    int files = 0;
    std::vector<std::string> wrong;
    const std::set<std::string> scanned = {".h", ".cpp", ".cu", ".cuh", ".txt", ".md", ".cmake"};
    for (const char* dir : {"src", "third_party", "cmake"})
        for (const auto& entry : fs::recursive_directory_iterator(fs::path(KNIT_SOURCE_DIR) / dir)) {
            if (!entry.is_regular_file() || !scanned.count(entry.path().extension().string())) continue;
            ++files;
            std::ifstream in(entry.path());
            std::string line;
            for (int number = 1; std::getline(in, line); ++number)
                if (std::regex_search(line, stale) || std::regex_search(line, flatKnitting))
                    wrong.push_back(entry.path().filename().string() + ":" + std::to_string(number));
        }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files > 200 && wrong.empty());
}

// Every case defined in tests/*.cpp is named in tests/CMakeLists.txt: registered, or named in a
// comment that says why not. An unnamed case would never run and nothing would notice. The GPU
// cases are exempt: CMake registers them from knit_tests_gpu --list.
TEST(architecture_every_test_registered) {
    namespace fs = std::filesystem;
    std::ifstream cmake(fs::path(KNIT_SOURCE_DIR) / "tests" / "CMakeLists.txt");
    const std::string lists((std::istreambuf_iterator<char>(cmake)), std::istreambuf_iterator<char>());
    const std::regex definition(R"(^TEST\((\w+)\))");
    std::vector<std::string> missing;
    int cases = 0;
    for (const auto& entry : fs::directory_iterator(fs::path(KNIT_SOURCE_DIR) / "tests")) {
        if (entry.path().extension() != ".cpp") continue;
        std::ifstream source(entry.path());
        std::string line;
        std::smatch name;
        while (std::getline(source, line)) {
            if (!std::regex_search(line, name, definition)) continue;
            ++cases;
            if (lists.find(name[1].str()) == std::string::npos) missing.push_back(name[1].str());
        }
    }
    for (const std::string& m : missing) std::printf("    not in tests/CMakeLists.txt: %s\n", m.c_str());
    CHECK(cases > 300 && missing.empty());
}

// src/model is plain values: no file I/O and no GPU qualifiers (src/model/ARCHITECTURE.md).
TEST(architecture_model_is_plain) {
    namespace fs = std::filesystem;
    const std::regex banned(R"(#\s*include\s*<(fstream|filesystem)>|__CUDACC__|__host__|__device__)");
    int files = 0;
    std::vector<std::string> wrong;
    for (const auto& entry : fs::directory_iterator(fs::path(KNIT_SOURCE_DIR) / "src" / "model")) {
        ++files;
        std::ifstream in(entry.path());
        std::string line;
        for (int number = 1; std::getline(in, line); ++number)
            if (std::regex_search(line, banned)) wrong.push_back(entry.path().filename().string() + ":" + std::to_string(number));
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files > 10 && wrong.empty());
}

TEST(architecture_no_new_templates) {
    const PatternCount found = countInProductionSources(std::regex(R"(^\s*template\s*<)"), "template");
    std::printf("  %d template declarations (pinned %d)\n", found.count, kPinnedTemplateDeclarations);
    if (found.count > kPinnedTemplateDeclarations) printPlaces(found);
    if (found.count < kPinnedTemplateDeclarations) std::printf("  lower kPinnedTemplateDeclarations to %d\n", found.count);
    CHECK(found.count == kPinnedTemplateDeclarations);
}

PatternCount genericLambdas() {
    // a lambda whose parameter list has an `auto` parameter: [captures](... auto ...)
    return countInProductionSources(std::regex(R"(\[[^\[\]]*\]\s*\([^)]*\bauto\b)"), "auto");
}

TEST(architecture_no_new_generic_lambdas) {
    const PatternCount found = genericLambdas();
    std::printf("  %d generic lambdas (pinned %d)\n", found.count, kPinnedGenericLambdas);
    if (found.count > kPinnedGenericLambdas) printPlaces(found);
    if (found.count < kPinnedGenericLambdas) std::printf("  lower kPinnedGenericLambdas to %d\n", found.count);
    CHECK(found.count == kPinnedGenericLambdas);
}

// src/export has one door (src/export/ARCHITECTURE.md): the program includes export.h, and --strand-demo
// the two layers it builds a package from; a format reads an ExportScene and never sees how it was
// made; and what is said about a format is in its table row, not in the viewer or the jobs.
TEST(architecture_export_door) {
    namespace fs = std::filesystem;
    const std::regex include(R"re(^\s*#\s*include\s+"([^"]+)")re"), formatLiteral(R"re("(obj|blend|ppt)")re");
    const std::vector<std::string> notForAFormat = {"settings/scene_registry.h", "export/collect_scene.h", "export/export.h"};
    int files = 0, formats = 0;
    std::vector<std::string> wrong;
    for (const auto& entry : fs::recursive_directory_iterator(fs::path(KNIT_SOURCE_DIR) / "src")) {
        if (!entry.is_regular_file()) continue;
        const std::string name = entry.path().filename().string(), folder = entry.path().parent_path().filename().string();
        const bool inExport = folder == "export";
        const bool neutral = inExport && (name.rfind("format_", 0) == 0 || name.rfind("export_scene", 0) == 0);
        const bool interfaceSide = entry.path().generic_string().find("/src/viewer/") != std::string::npos || entry.path().generic_string().find("/src/jobs/") != std::string::npos;
        ++files;
        formats += neutral;
        std::ifstream in(entry.path());
        std::string line;
        std::smatch m;
        while (std::getline(in, line)) {
            if (interfaceSide && std::regex_search(line, m, formatLiteral)) wrong.push_back(name + " names the format " + m[1].str());
            if (!std::regex_search(line, m, include)) continue;
            const std::string header = m[1].str();
            if (neutral && (header.rfind("gpu/", 0) == 0 || header.rfind("scene/", 0) == 0 ||
                            std::find(notForAFormat.begin(), notForAFormat.end(), header) != notForAFormat.end()))
                wrong.push_back(name + " includes " + header);
            if (inExport || header.rfind("export/", 0) != 0 || header == "export/export.h") continue;
            if (name != "pieces_cli.cpp" || header != "export/collect_scene.h") wrong.push_back(name + " goes round the door: " + header);
        }
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files > 200 && formats >= 5 && wrong.empty());
}

// src/gpu's folders use each other in one direction (src/gpu/ARCHITECTURE.md): each folder includes itself,
// the folders listed and the public headers listed; the port, at the top, composes the capabilities.
// device/ and device.h, which everything includes, stay free of knit types apart from the cancel token.
// Every folder of src/gpu has a layer here.
TEST(architecture_gpu_layers) {
    namespace fs = std::filesystem;
    struct Layer {
        std::set<std::string> folders, publicHeaders;
    };
    const std::map<std::string, Layer> allowed = {
        {"device", {{"device"}, {"device.h"}}},
        {"geometry", {{"device", "geometry"}, {"device.h"}}},
        {"solver", {{"device", "geometry", "solver"}, {"device.h"}}},
        {"validation", {{"device", "geometry", "validation"}, {"device.h"}}},
        {"surface", {{"device", "geometry", "surface"}, {"device.h", "surfaces.h"}}},
        {"fibres", {{"device", "geometry", "fibres"}, {"device.h", "fibre_generator.h"}}},
        {"port", {{"device", "geometry", "solver", "validation", "surface", "port"}, {"device.h", "knitting.h", "surfaces.h"}}},
    };
    const std::regex folderInclude(R"re(^\s*#\s*include\s+"gpu/(\w+)/)re"), publicInclude(R"re(^\s*#\s*include\s+"gpu/(\w+\.h)")re"),
        projectInclude(R"re(^\s*#\s*include\s+"(\w+)/([\w/]+\.\w+)")re");
    const fs::path gpu = fs::path(KNIT_SOURCE_DIR) / "src" / "gpu";
    int files = 0;
    std::vector<std::string> wrong;
    for (const auto& entry : fs::directory_iterator(gpu))
        if (entry.is_directory() && !allowed.count(entry.path().filename().string()))
            wrong.push_back("a folder without a layer: " + entry.path().filename().string());
    for (const auto& [folder, may] : allowed)
        for (const auto& entry : fs::recursive_directory_iterator(gpu / folder)) {
            if (!entry.is_regular_file()) continue;
            ++files;
            std::ifstream in(entry.path());
            std::string line;
            std::smatch m;
            while (std::getline(in, line)) {
                if (std::regex_search(line, m, folderInclude) && !may.folders.count(m[1].str()))
                    wrong.push_back(folder + " includes " + m[1].str() + "/: " + entry.path().filename().string());
                if (std::regex_search(line, m, publicInclude) && !may.publicHeaders.count(m[1].str()))
                    wrong.push_back(folder + " includes gpu/" + m[1].str() + ": " + entry.path().filename().string());
                if (folder == "device" && std::regex_search(line, m, projectInclude) && m[1].str() != "gpu" && m[0].str().find("model/cancel.h") == std::string::npos)
                    wrong.push_back("device includes " + m[1].str() + "/" + m[2].str() + ": " + entry.path().filename().string());
            }
        }
    std::ifstream deviceHeader(gpu / "device.h");
    std::string line;
    std::smatch m;
    while (std::getline(deviceHeader, line))
        if (std::regex_search(line, m, projectInclude)) wrong.push_back("device.h includes " + m[1].str() + "/" + m[2].str());
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files >= 25 && wrong.empty());
}

// src/gpu is written without templates: concrete accessors, overloads, cursors, two short kernels.
TEST(architecture_gpu_no_templates) {
    std::vector<std::string> inGpu;
    for (const PatternCount& found : {countInProductionSources(std::regex(R"(^\s*template\s*<)"), "template"), genericLambdas()})
        for (const std::string& place : found.places)
            if (place.rfind("src/gpu/", 0) == 0) inGpu.push_back(place);
    for (const std::string& place : inGpu) std::printf("    %s\n", place.c_str());
    CHECK(inGpu.empty());
}

// Every CUDA call in src/gpu is checked (an argument of check), handled explicitly (its result assigned, returned, or tested inside a condition), a destroy or
// free (cleanup, which must not throw), an error-text function, a device function without a status (the
// solver's graph conditions), or discarded on purpose with (void). A call whose result is dropped otherwise fails. Statements, not lines: comments and literals are blanked first, and a
// statement may span lines.
TEST(architecture_gpu_cuda_calls_checked) {
    namespace fs = std::filesystem;
    const std::regex call(R"(cuda[A-Z]\w*\s*\()"), assignment(R"([^=!<>]=[^=]|^=[^=])");
    const std::set<std::string> errorText = {"cudaGetErrorString", "cudaGetLastError"};
    const std::set<std::string> noStatus = {"cudaGraphSetConditional"};  // device code, returns void
    int files = 0, calls = 0;
    std::vector<std::string> wrong;
    for (const auto& entry : fs::recursive_directory_iterator(fs::path(KNIT_SOURCE_DIR) / "src" / "gpu")) {
        const std::string ext = entry.path().extension().string();
        if (ext != ".cu" && ext != ".cpp" && ext != ".cuh" && ext != ".h") continue;
        ++files;
        std::ifstream in(entry.path(), std::ios::binary);
        std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        // Blank comments, string and character literals (newlines kept for the line numbers).
        for (size_t i = 0; i < text.size(); ++i) {
            if (text.compare(i, 2, "//") == 0) {
                while (i < text.size() && text[i] != '\n') text[i++] = ' ';
            } else if (text.compare(i, 2, "/*") == 0) {
                for (; i < text.size() && text.compare(i, 2, "*/") != 0; ++i)
                    if (text[i] != '\n') text[i] = ' ';
                if (i < text.size()) text[i] = text[i + 1] = ' ';
            } else if (text[i] == '"' || text[i] == '\'') {
                const char quote = text[i];
                for (++i; i < text.size() && text[i] != quote; ++i) {
                    if (text[i] == '\\') text[i++] = ' ';
                    if (text[i] != '\n') text[i] = ' ';
                }
            }
        }
        size_t start = 0;
        for (size_t i = 0; i <= text.size(); ++i) {
            if (i < text.size() && text[i] != ';' && text[i] != '{' && text[i] != '}') continue;
            const std::string statement = text.substr(start, i - start);
            for (std::sregex_iterator it(statement.begin(), statement.end(), call), end; it != end; ++it) {
                const std::string name = it->str().substr(0, it->str().find_first_of(" \t\r\n("));
                const std::string before = statement.substr(0, it->position());
                ++calls;
                if (name.find("Destroy") != std::string::npos || name.find("Free") != std::string::npos || errorText.count(name) || noStatus.count(name)) continue;
                std::string trimmed = before;
                while (!trimmed.empty() && std::isspace(static_cast<unsigned char>(trimmed.back()))) trimmed.pop_back();
                const bool checked = trimmed.size() >= 6 && trimmed.ends_with("check(");
                const bool discarded = trimmed.ends_with("(void)");
                int depth = 0;
                for (char ch : before) depth += ch == '(' ? 1 : ch == ')' ? -1 : 0;
                std::string head = before;
                head.erase(0, head.find_first_not_of(" \t\r\n"));
                const bool inCondition = depth > 0 && (head.rfind("if", 0) == 0 || head.rfind("while", 0) == 0);
                const bool used = std::regex_search(before, assignment) || head.rfind("return", 0) == 0;
                if (checked || discarded || inCondition || used) continue;
                const long long line = 1 + std::count(text.begin(), text.begin() + static_cast<std::ptrdiff_t>(start + it->position()), '\n');
                wrong.push_back(fs::relative(entry.path(), fs::path(KNIT_SOURCE_DIR)).generic_string() + ":" + std::to_string(line) + " " + name);
            }
            start = i + 1;
        }
    }
    for (const std::string& w : wrong) std::printf("    result dropped: %s\n", w.c_str());
    CHECK(files >= 25 && calls >= 100 && wrong.empty());
}

namespace {

// A folder's ARCHITECTURE.md is the way into it: its file tree (from "## File tree" to "**Reading order")
// names every file there, and only those; the introduction above it may name other folders' files.
// Returns the files in question, after counting them.
std::vector<std::string> architectureDocGaps(const std::filesystem::path& root, size_t* fileCount) {
    namespace fs = std::filesystem;
    std::ifstream in(root / "ARCHITECTURE.md");
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const size_t treeStart = text.find("## File tree");
    const size_t from = treeStart == std::string::npos ? 0 : treeStart;
    const std::string tree = text.substr(from, text.find("**Reading order") - from);
    std::set<std::string> existing;
    for (const auto& entry : fs::recursive_directory_iterator(root))
        if (entry.is_regular_file()) existing.insert(entry.path().filename().string());
    // In the tree "knit_numerics.h, .cpp" names two files: a bare extension repeats the stem before it.
    std::set<std::string> named;
    const std::regex token(R"((\w*)(\.(?:h|cu|cuh|cpp|md|txt))\b)");
    std::string stem;
    for (std::sregex_iterator it(tree.begin(), tree.end(), token), end; it != end; ++it) {
        if (!(*it)[1].str().empty()) stem = (*it)[1].str();
        named.insert(stem + (*it)[2].str());
    }
    std::vector<std::string> wrong;
    for (const std::string& f : existing)
        if (!named.count(f)) wrong.push_back("not in ARCHITECTURE.md: " + f);
    for (const std::string& f : named)
        if (!existing.count(f)) wrong.push_back("no such file: " + f);
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    *fileCount = existing.size();
    return wrong;
}

} // namespace

TEST(architecture_gpu_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "gpu", &files);
    CHECK(files >= 25 && wrong.empty());
}

TEST(architecture_import_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "import", &files);
    CHECK(files >= 20 && wrong.empty());
}

// The import's layers (src/import/ARCHITECTURE.md): each file includes only the import files its row names, so it can
// be read knowing only those. model and io may be included anywhere; knitting only by import_analysis.cpp, which asks the
// knit's own predicates; each third-party header only by the one file that wraps it (angle-bracket includes, such as
// SLIM's Eigen, are not scanned). A public header includes no private one, so jobs and scene see neither the rules nor
// the built-in engine.
TEST(architecture_import_layers) {
    namespace fs = std::filesystem;
    const std::map<std::string, std::set<std::string>> allowed = {
        {"repair_rules", {}},
        {"wrap_stages", {}},
        {"wrap_grid", {"wrap_stages"}},
        {"wrap_nets", {"wrap_stages"}},
        {"wrap_remesh", {"wrap_stages"}},
        {"wrap_tidy", {"wrap_stages"}},
        {"wrap_charts", {"wrap_stages"}},
        {"wrap_flatten", {"wrap_stages"}},
        {"wrap_quad_split", {"wrap_stages"}},
        {"wrap_pack", {"wrap_stages"}},
        {"wrap_bake", {"wrap_stages"}},
        {"wrap_surface", {"wrap_stages", "prepare_import"}},
        {"prepare_import", {"repair_rules"}},
        {"import_analysis", {"prepare_import", "repair_rules"}},
        {"import_store", {"prepare_import"}},
        {"import_pipeline", {"prepare_import", "import_store", "wrap_surface"}},
    };
    const std::set<std::string> publicStems = {"import_analysis", "prepare_import", "import_pipeline", "import_store"};
    const std::map<std::string, std::string> thirdParty = {{"slim.h", "wrap_flatten.cpp"}, {"level_set_remesh.h", "wrap_remesh.cpp"}};
    const fs::path root = fs::path(KNIT_SOURCE_DIR) / "src" / "import";
    const std::regex include(R"re(^\s*#\s*include\s+"([^"]+)")re"), ours(R"(import/(\w+)\.h)");
    int files = 0;
    std::vector<std::string> wrong;
    for (const fs::directory_entry& entry : fs::directory_iterator(root)) {
        if (!entry.is_regular_file() || entry.path().extension() == ".md" || entry.path().filename() == "CMakeLists.txt") continue;
        ++files;
        const std::string name = entry.path().filename().string(), stem = entry.path().stem().string();
        const bool publicHeader = entry.path().extension() == ".h" && publicStems.count(stem);
        if (!allowed.count(stem)) {
            wrong.push_back(name + " sits outside the layers");
            continue;
        }
        std::ifstream in(entry.path());
        std::string line;
        std::smatch m, to;
        while (std::getline(in, line)) {
            if (!std::regex_search(line, m, include)) continue;
            const std::string target = m[1].str();
            if (target.rfind("model/", 0) == 0 || target.rfind("io/", 0) == 0) continue;
            if (target.rfind("knitting/", 0) == 0) {
                wrong.push_back(name + " includes " + target + " (the import does not depend on knitting)");
                continue;
            }
            if (thirdParty.count(target)) {
                if (thirdParty.at(target) != name) wrong.push_back(name + " includes " + target + " (only " + thirdParty.at(target) + " may)");
                continue;
            }
            if (!std::regex_match(target, to, ours)) {
                wrong.push_back(name + " includes " + target);
                continue;
            }
            const std::string other = to[1].str();
            if (other != stem && !allowed.at(stem).count(other)) wrong.push_back(name + " includes " + target);
            if (publicHeader && !publicStems.count(other)) wrong.push_back(name + " (a public header) includes the private " + target);
        }
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files >= 20 && wrong.empty());
}

TEST(architecture_jobs_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "jobs", &files);
    CHECK(files >= 20 && wrong.empty());
}

// The jobs' layers (src/jobs/ARCHITECTURE.md): each file includes only the jobs files its row names, so the jobs
// never see the queue, the display builders never see a job, and a public header includes no private one (the
// viewer and knit_check see neither the display builders nor the jobs' shared steps). Other components are
// checked by architecture_module_dependencies.
TEST(architecture_jobs_layers) {
    namespace fs = std::filesystem;
    const std::set<std::string> data = {"built_scene", "results"};
    const std::set<std::string> jobSteps = {"job_steps", "jobs", "display_scene", "scene_cache", "built_scene", "results"};
    const std::map<std::string, std::set<std::string>> allowed = {
        {"results", {}},
        {"crash_report", {}},
        {"built_scene", {"results"}},
        {"jobs", data},
        {"display_scene", {"built_scene", "results"}},
        {"scene_cache", {"built_scene", "results"}},
        {"job_steps", {"jobs", "display_scene", "scene_cache", "built_scene", "results"}},
        {"previews", {"display_scene", "built_scene", "results"}},
        {"fur_preview", {"display_scene", "built_scene", "results"}},
        {"import_session", {"display_scene", "built_scene", "results"}},
        {"marking_session", {"built_scene", "results"}},
        {"knit_job", jobSteps},
        {"strand_mesh_job", jobSteps},
        {"fibre_job", jobSteps},
        {"export_job", jobSteps},
        {"render_preview", jobSteps},
        {"scene_builder", {"jobs", "built_scene", "results"}},
    };
    const std::set<std::string> publicStems = {"built_scene", "crash_report", "fur_preview", "import_session", "jobs", "marking_session", "previews", "render_preview",
                                               "results", "scene_builder", "scene_cache"};
    const fs::path root = fs::path(KNIT_SOURCE_DIR) / "src" / "jobs";
    const std::regex include(R"re(^\s*#\s*include\s+"([^"]+)")re"), ours(R"(jobs/(\w+)\.h)");
    int files = 0;
    std::vector<std::string> wrong;
    for (const fs::directory_entry& entry : fs::directory_iterator(root)) {
        if (!entry.is_regular_file() || entry.path().extension() == ".md" || entry.path().filename() == "CMakeLists.txt") continue;
        ++files;
        const std::string name = entry.path().filename().string(), stem = entry.path().stem().string();
        const bool publicHeader = entry.path().extension() == ".h" && publicStems.count(stem);
        if (!allowed.count(stem)) {
            wrong.push_back(name + " sits outside the layers");
            continue;
        }
        std::ifstream in(entry.path());
        std::string line;
        std::smatch m, to;
        while (std::getline(in, line)) {
            if (!std::regex_search(line, m, include)) continue;
            const std::string target = m[1].str();
            if (!std::regex_match(target, to, ours)) continue;  // another component: architecture_module_dependencies
            const std::string other = to[1].str();
            if (other != stem && !allowed.at(stem).count(other)) wrong.push_back(name + " includes " + target);
            if (publicHeader && !publicStems.count(other)) wrong.push_back(name + " (a public header) includes the private " + target);
        }
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files >= 20 && wrong.empty());
}

TEST(architecture_viewer_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "viewer", &files);
    CHECK(files >= 20 && wrong.empty());
}

TEST(architecture_fibres_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "fibres", &files);
    CHECK(files >= 20 && wrong.empty());
}

TEST(architecture_knitting_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "knitting", &files);
    CHECK(files >= 45 && wrong.empty());
}

TEST(architecture_stitchmesh_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "stitchmesh", &files);
    CHECK(files >= 10 && wrong.empty());
}

TEST(architecture_model_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "model", &files);
    CHECK(files >= 25 && wrong.empty());
}

TEST(architecture_io_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "io", &files);
    CHECK(files >= 15 && wrong.empty());
}

TEST(architecture_settings_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "settings", &files);
    CHECK(files >= 12 && wrong.empty());
}

TEST(architecture_export_doc_complete) {
    size_t files = 0;
    const std::vector<std::string> wrong = architectureDocGaps(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "export", &files);
    CHECK(files >= 10 && wrong.empty());
}

namespace {

// One component as build/components.json declares it, or as the table in src/ARCHITECTURE.md states it.
struct ComponentDeclaration {
    std::set<std::string> publicDeps, privateDeps;
    int publicHeaders = 0, privateHeaders = 0;
};

// The quoted names in a JSON list: ["model", "io"] -> {model, io}.
std::set<std::string> quotedNames(const std::string& list) {
    const std::regex quoted(R"re("([^"]+)")re");
    std::set<std::string> names;
    for (std::sregex_iterator it(list.begin(), list.end(), quoted), end; it != end; ++it) names.insert((*it)[1].str());
    return names;
}

std::map<std::string, ComponentDeclaration> componentsFromJson() {
    std::ifstream json(KNIT_COMPONENTS_JSON);
    const std::string text((std::istreambuf_iterator<char>(json)), std::istreambuf_iterator<char>());
    const std::regex entry(R"re(\{"name": "(\w+)", "target": "\w+",\s*"public_deps": \[([^\]]*)\],\s*"private_deps": \[([^\]]*)\],)re"
                           R"re(\s*"public_headers": \[([^\]]*)\],\s*"private_headers": \[([^\]]*)\])re");
    std::map<std::string, ComponentDeclaration> components;
    for (std::sregex_iterator it(text.begin(), text.end(), entry), end; it != end; ++it) {
        ComponentDeclaration& c = components[(*it)[1].str()];
        c.publicDeps = quotedNames((*it)[2].str());
        c.privateDeps = quotedNames((*it)[3].str());
        c.publicHeaders = static_cast<int>(quotedNames((*it)[4].str()).size());
        c.privateHeaders = static_cast<int>(quotedNames((*it)[5].str()).size());
    }
    return components;
}

// A dependency cell of the table: "model, io", "none", "(CUDA)" (an outside library, not a component) or
// "jobs; export's private headers" (access to a component's private headers counts as depending on it).
std::set<std::string> componentsInCell(std::string cell) {
    cell = std::regex_replace(cell, std::regex(R"(\([^)]*\))"), "");
    std::replace(cell.begin(), cell.end(), ';', ',');
    std::set<std::string> names;
    std::stringstream items(cell);
    std::string item;
    while (std::getline(items, item, ',')) {
        std::smatch name;
        if (!std::regex_search(item, name, std::regex(R"((\w+))"))) continue;
        if (name[1].str() != "none") names.insert(name[1].str());
    }
    return names;
}

// The dependency graph table of src/ARCHITECTURE.md: | `folder` | public | private | headers public / private |.
std::map<std::string, ComponentDeclaration> componentsFromTable() {
    std::ifstream in(std::filesystem::path(KNIT_SOURCE_DIR) / "src" / "ARCHITECTURE.md");
    const std::regex row(R"(^\| `(\w+)` \| ([^|]+)\| ([^|]+)\| ([^|]+)\|\s*$)");
    const std::regex counts(R"((\d+) / (\d+))");
    std::map<std::string, ComponentDeclaration> components;
    std::string line;
    while (std::getline(in, line)) {
        std::smatch m, n;
        if (!std::regex_match(line, m, row)) continue;
        const std::string headerCell = m[4].str();
        const bool hasCounts = std::regex_search(headerCell, n, counts);
        if (!hasCounts && headerCell.find("none") == std::string::npos) continue;  // the folder table above it
        ComponentDeclaration& c = components[m[1].str()];
        c.publicDeps = componentsInCell(m[2].str());
        c.privateDeps = componentsInCell(m[3].str());
        c.publicHeaders = hasCounts ? std::stoi(n[1].str()) : 0;
        c.privateHeaders = hasCounts ? std::stoi(n[2].str()) : 0;
    }
    return components;
}

std::string joinedNames(const std::set<std::string>& names) {
    std::string text;
    for (const std::string& name : names) text += (text.empty() ? "" : ", ") + name;
    return text.empty() ? "none" : text;
}

} // namespace

// The dependency graph table of src/ARCHITECTURE.md states what the CMakeLists.txt files declare: every
// component's public and private dependencies and its number of public and private headers.
TEST(architecture_graph_table_current) {
    const std::map<std::string, ComponentDeclaration> declared = componentsFromJson();
    const std::map<std::string, ComponentDeclaration> table = componentsFromTable();
    std::vector<std::string> wrong;
    for (const auto& [name, d] : declared) {
        if (!table.count(name)) {
            wrong.push_back(name + ": no row in the table");
            continue;
        }
        const ComponentDeclaration& t = table.at(name);
        if (t.publicDeps != d.publicDeps)
            wrong.push_back(name + ": public deps " + joinedNames(t.publicDeps) + " in the table, " + joinedNames(d.publicDeps) + " declared");
        if (t.privateDeps != d.privateDeps)
            wrong.push_back(name + ": private deps " + joinedNames(t.privateDeps) + " in the table, " + joinedNames(d.privateDeps) + " declared");
        if (t.publicHeaders != d.publicHeaders || t.privateHeaders != d.privateHeaders)
            wrong.push_back(name + ": headers " + std::to_string(t.publicHeaders) + " / " + std::to_string(t.privateHeaders) + " in the table, " +
                            std::to_string(d.publicHeaders) + " / " + std::to_string(d.privateHeaders) + " declared");
    }
    for (const auto& [name, t] : table)
        if (!declared.count(name)) wrong.push_back(name + ": a table row without a component");
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(declared.size() >= 12 && wrong.empty());
}

// src/model's files sit in four layers, and their includes run one way (src/model/ARCHITECTURE.md): the
// foundation includes only the foundation; the vocabulary, the curves and the surfaces each include the
// foundation and their own layer. Only a forbidden include fails: removing one never does.
TEST(architecture_model_layers) {
    namespace fs = std::filesystem;
    const std::map<std::string, std::string> layerOf = {
        {"vec", "foundation"},        {"coordinates", "foundation"}, {"colour", "foundation"},
        {"status", "foundation"},     {"hash", "foundation"},        {"env", "foundation"},
        {"log", "foundation"},        {"compute_caps", "foundation"}, {"cancel", "foundation"},
        {"parallel", "foundation"},   {"noise", "foundation"},       {"params", "vocabulary"},      {"choices", "vocabulary"},
        {"bspline", "curves"},        {"topology", "curves"},        {"yarn_path", "curves"},
        {"spatial_hash", "surfaces"}, {"mesh", "surfaces"},          {"mesh_query", "surfaces"},
        {"mesh_check", "surfaces"},   {"surface_marking", "markings"}, {"marking_review", "markings"},
        {"surface_path", "markings"}, {"surface_pieces", "markings"}, {"stitch_mesh_data", "markings"},
    };
    const std::map<std::string, std::set<std::string>> allowed = {
        {"foundation", {"foundation"}},
        {"vocabulary", {"foundation", "vocabulary"}},
        {"curves", {"foundation", "curves"}},
        {"surfaces", {"foundation", "surfaces"}},
        {"markings", {"foundation", "surfaces", "markings"}},
    };
    const fs::path root = fs::path(KNIT_SOURCE_DIR) / "src" / "model";
    const std::regex include(R"re(^\s*#\s*include\s+"model/(\w+)\.h")re");
    int files = 0;
    std::vector<std::string> wrong;
    for (const fs::directory_entry& entry : fs::directory_iterator(root)) {
        const std::string extension = entry.path().extension().string();
        if (extension != ".h" && extension != ".cpp") continue;  // ARCHITECTURE.md, CMakeLists.txt
        ++files;
        const std::string name = entry.path().filename().string(), stem = entry.path().stem().string();
        if (!layerOf.count(stem)) {
            wrong.push_back(name + " sits outside the layers");
            continue;
        }
        const std::set<std::string>& may = allowed.at(layerOf.at(stem));
        std::ifstream in(entry.path());
        std::string line;
        std::smatch m;
        while (std::getline(in, line)) {
            if (!std::regex_search(line, m, include)) continue;
            const std::string target = m[1].str();
            if (!layerOf.count(target) || !may.count(layerOf.at(target))) wrong.push_back(name + " includes model/" + target + ".h");
        }
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files >= 25 && wrong.empty());
}

// The viewer's layers (src/viewer/ARCHITECTURE.md). The compiler holds the state: its target has no
// ImGui on its include path. This holds the rest: the interface has no GL and no GLFW, so the UI tests
// draw it without a window; the 3D view has no ImGui; presentation.cpp and view_content.cpp only read.
TEST(architecture_viewer_layers) {
    namespace fs = std::filesystem;
    const fs::path root = fs::path(KNIT_SOURCE_DIR) / "src" / "viewer";
    const std::regex gl(R"re(#\s*include\s+[<"](GL/|GLFW/|gl_functions\.h|renderer\.h|mesh_lighting\.h|fibre_lines\.h|fur_discs\.h|imgui_impl))re");
    const std::regex imgui(R"re(#\s*include\s+"(imgui|interface\.h))re");
    const std::regex guiCalls(R"(\bImGui::)");
    const std::map<std::string, std::vector<const std::regex*>> banned = {
        {"workflow.h", {&gl, &imgui, &guiCalls}},        {"workflow.cpp", {&gl, &imgui, &guiCalls}},
        {"presentation.cpp", {&gl, &imgui, &guiCalls}},  {"viewer_script.h", {&gl, &imgui, &guiCalls}},
        {"viewer_script.cpp", {&gl, &imgui, &guiCalls}}, {"camera.h", {&gl, &imgui, &guiCalls}},
        {"camera.cpp", {&gl, &imgui, &guiCalls}},        {"viewer_prefs.cpp", {&gl, &imgui, &guiCalls}},
        {"about.h", {&gl, &imgui, &guiCalls}},           {"about.cpp", {&gl, &imgui, &guiCalls}},
        {"about_texts.h", {&gl, &imgui, &guiCalls}},     {"about_window.cpp", {&gl}},
        {"view_content.h", {&gl, &imgui, &guiCalls}},    {"view_content.cpp", {&gl, &imgui, &guiCalls}},
        {"interface.h", {&gl}},                          {"interface.cpp", {&gl}},
        {"widgets.cpp", {&gl}},                          {"panel_import.cpp", {&gl}},
        {"panel_knit.cpp", {&gl}},                       {"panel_fibres.cpp", {&gl}},  {"panel_fur.cpp", {&gl}},
        {"panel_yarn_fibres.cpp", {&gl}},                {"panel_mark.cpp", {&gl}},
        {"marking_flow.h", {&gl, &imgui, &guiCalls}},    {"marking_flow.cpp", {&gl, &imgui, &guiCalls}},
        {"overlay_lines.h", {&imgui, &guiCalls}},        {"overlay_lines.cpp", {&imgui, &guiCalls}},
        {"panel_export.cpp", {&gl}},                     {"renderer.h", {&imgui, &guiCalls}},
        {"renderer.cpp", {&imgui, &guiCalls}},           {"mesh_lighting.h", {&imgui, &guiCalls}},
        {"mesh_lighting.cpp", {&imgui, &guiCalls}},      {"fibre_lines.h", {&imgui, &guiCalls}},
        {"fibre_lines.cpp", {&imgui, &guiCalls}},        {"fur_discs.h", {&imgui, &guiCalls}}, {"fur_discs.cpp", {&imgui, &guiCalls}},        {"gl_functions.h", {&imgui, &guiCalls}},
    };
    std::vector<std::string> wrong;
    for (const auto& [file, patterns] : banned) {
        std::ifstream in(root / file);
        if (!in) wrong.push_back("no such file: " + file);
        std::string line;
        for (int number = 1; std::getline(in, line); ++number)
            for (const std::regex* p : patterns)
                if (std::regex_search(line, *p)) wrong.push_back(file + ":" + std::to_string(number));
    }
    // Every SceneUi in presentation.cpp and view_content.cpp is const: what changes it is in workflow.cpp.
    const std::regex mutableUi(R"((^|[^\w])SceneUi&)");
    const std::regex constUi(R"(const SceneUi&)");
    int readers = 0;
    for (const char* reader : {"presentation.cpp", "view_content.cpp"}) {
        std::ifstream in(root / reader);
        std::string line;
        for (int number = 1; std::getline(in, line); ++number) {
            const std::string rest = std::regex_replace(line, constUi, "");
            if (std::regex_search(rest, mutableUi)) wrong.push_back(std::string(reader) + ":" + std::to_string(number) + " takes a mutable SceneUi");
            readers += rest != line;
        }
    }
    // Each layer includes only the components it uses (a direct #include "<component>/..."); every file is in a layer.
    const std::vector<std::string> state = {"model", "io", "settings", "export", "jobs"};
    const std::vector<std::string> interface = {"model", "settings", "export", "jobs"};
    const std::vector<std::string> view3d = {"model", "io", "jobs"};
    const std::vector<std::string> program = {"model", "jobs"};
    const std::map<std::string, const std::vector<std::string>*> layerOf = {
        {"workflow.h", &state},       {"workflow.cpp", &state},       {"presentation.cpp", &state},   {"view_content.h", &state},
        {"view_content.cpp", &state}, {"viewer_script.h", &state},    {"viewer_script.cpp", &state},  {"viewer_prefs.h", &state},
        {"viewer_prefs.cpp", &state}, {"camera.h", &state},           {"camera.cpp", &state},         {"interface.h", &interface},
        {"about.h", &state},          {"about.cpp", &state},          {"about_texts.h", &state},      {"about_window.cpp", &interface},
        {"interface.cpp", &interface}, {"widgets.cpp", &interface},   {"panel_import.cpp", &interface}, {"panel_knit.cpp", &interface},
        {"panel_fibres.cpp", &interface}, {"panel_fur.cpp", &interface}, {"panel_export.cpp", &interface}, {"renderer.h", &view3d},   {"renderer.cpp", &view3d},
        {"panel_yarn_fibres.cpp", &interface}, {"panel_mark.cpp", &interface}, {"marking_flow.h", &state}, {"marking_flow.cpp", &state},
        {"overlay_lines.h", &view3d}, {"overlay_lines.cpp", &view3d},
        {"mesh_lighting.h", &view3d}, {"mesh_lighting.cpp", &view3d}, {"fibre_lines.h", &view3d},     {"fibre_lines.cpp", &view3d}, {"fur_discs.h", &view3d}, {"fur_discs.cpp", &view3d},
        {"gl_functions.h", &view3d},  {"main.cpp", &program},         {"native_dialogs.h", &program}, {"native_dialogs.cpp", &program},
    };
    const std::regex projectInclude(R"re(#\s*include\s+"(\w+)/)re");
    const std::vector<std::string> components = {"model", "knitting", "fibres", "gpu", "io", "settings", "import", "scene", "export", "jobs", "cli"};
    for (const fs::directory_entry& entry : fs::directory_iterator(root)) {
        const std::string file = entry.path().filename().string();
        const std::string ext = entry.path().extension().string();
        if (ext != ".h" && ext != ".cpp") continue;
        if (layerOf.count(file) == 0) {
            wrong.push_back(file + " is in no layer");
            continue;
        }
        const std::vector<std::string>& allowed = *layerOf.at(file);
        std::ifstream in(entry.path());
        std::string line;
        std::smatch m;
        for (int number = 1; std::getline(in, line); ++number) {
            if (!std::regex_search(line, m, projectInclude)) continue;
            const std::string component = m[1].str();
            const bool known = std::find(components.begin(), components.end(), component) != components.end();
            if (known && std::find(allowed.begin(), allowed.end(), component) == allowed.end())
                wrong.push_back(file + ":" + std::to_string(number) + " includes " + component + ", which its layer does not use");
        }
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(readers > 30 && wrong.empty());
}

// src/scene's folders follow the pipeline, and their includes run one way (src/scene/ARCHITECTURE.md).
TEST(architecture_scene_layers) {
    namespace fs = std::filesystem;
    const std::map<std::string, std::set<std::string>> allowed = {
        {"", {""}},
        {"source", {"", "source"}},
        {"appearance", {"", "appearance"}},
        {"knit", {"", "source", "appearance", "knit"}},
        {"fibres", {"", "source", "appearance", "fibres"}},
    };
    const fs::path root = fs::path(KNIT_SOURCE_DIR) / "src" / "scene";
    const std::regex include(R"re(^\s*#\s*include\s+"scene/([^"]+)")re");
    int files = 0;
    std::vector<std::string> wrong;
    for (const fs::directory_entry& entry : fs::recursive_directory_iterator(root)) {
        if (!entry.is_regular_file() || entry.path().extension() == ".md" || entry.path().filename() == "CMakeLists.txt") continue;
        ++files;
        const std::string name = entry.path().lexically_relative(root).generic_string();
        const size_t slash = name.find('/');
        const std::string folder = slash == std::string::npos ? "" : name.substr(0, slash);
        if (!allowed.count(folder)) {
            wrong.push_back(name + " sits outside the layers");
            continue;
        }
        std::ifstream in(entry.path());
        std::string line;
        std::smatch m;
        while (std::getline(in, line)) {
            if (!std::regex_search(line, m, include)) continue;
            const std::string target = m[1].str();
            const size_t cut = target.find('/');
            if (!allowed.at(folder).count(cut == std::string::npos ? "" : target.substr(0, cut))) wrong.push_back(name + " includes " + target);
        }
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files >= 45 && wrong.empty());
}

// ARCHITECTURE.md's tree names every file of src/scene, and no file that is gone.
TEST(architecture_scene_doc_complete) {
    namespace fs = std::filesystem;
    const fs::path root = fs::path(KNIT_SOURCE_DIR) / "src" / "scene";
    std::ifstream doc(root / "ARCHITECTURE.md");
    std::string line, tree;
    bool inTree = false, seenTree = false;
    while (std::getline(doc, line)) {
        if (line.rfind("```", 0) == 0) {
            if (seenTree) break;
            inTree = seenTree = true;
            continue;
        }
        if (inTree) tree += line + "\n";
    }
    CHECK(!tree.empty());
    std::set<std::string> named;
    // A file is named by its header or, without one, by its source.
    const std::regex file(R"re(\b([a-z_]+)\.(?:h|cpp)\b)re");
    for (auto it = std::sregex_iterator(tree.begin(), tree.end(), file); it != std::sregex_iterator(); ++it) named.insert((*it)[1].str());
    std::set<std::string> onDisk;
    for (const fs::directory_entry& entry : fs::recursive_directory_iterator(root))
        if (entry.is_regular_file() && entry.path().extension() != ".md" && entry.path().filename() != "CMakeLists.txt")
            onDisk.insert(entry.path().stem().string());
    int wrong = 0;
    for (const std::string& name : onDisk)
        if (!named.count(name)) {
            std::printf("    %s is not in the tree\n", name.c_str());
            ++wrong;
        }
    for (const std::string& name : named)
        if (!onDisk.count(name)) {
            std::printf("    %s is gone\n", name.c_str());
            ++wrong;
        }
    CHECK(named.size() >= 25 && wrong == 0);
}

// src/knitting's folders are the stages of knitting, and their includes run one way
// (src/knitting/ARCHITECTURE.md): a folder includes itself, the ones listed and src/model. numerics.h,
// the port, sits at the root: any folder may include it, and it includes only src/model, so it adds
// no edge between them.
TEST(architecture_knitting_layers) {
    namespace fs = std::filesystem;
    const std::map<std::string, std::set<std::string>> allowed = {
        {"stitch", {"stitch"}},
        {"checks", {"stitch", "checks"}},
        {"layout", {"stitch", "layout"}},
        {"tube", {"tube"}},
        {"relax", {"stitch", "checks", "layout", "tube", "relax"}},
        {"knit", {"stitch", "checks", "layout", "tube", "relax", "knit"}},
    };
    const fs::path root = fs::path(KNIT_SOURCE_DIR) / "src" / "knitting";
    const std::regex include(R"re(^\s*#\s*include\s+"([^"]+)")re"), layered(R"(knitting/(\w+)/[\w.]+)");
    int files = 0;
    std::vector<std::string> wrong;
    for (const fs::directory_entry& entry : fs::recursive_directory_iterator(root)) {
        if (!entry.is_regular_file() || entry.path().extension() == ".md" || entry.path().filename() == "CMakeLists.txt") continue;
        ++files;
        const std::string name = entry.path().lexically_relative(root).generic_string();
        const std::string folder = name.substr(0, name.find('/'));
        const bool port = name == "numerics.h";
        if (!port && !allowed.count(folder)) {
            wrong.push_back(name + " sits outside the layers");
            continue;
        }
        std::ifstream in(entry.path());
        std::string line;
        std::smatch m, to;
        while (std::getline(in, line)) {
            if (!std::regex_search(line, m, include)) continue;
            const std::string target = m[1].str();
            if (target.rfind("model/", 0) == 0) continue;
            const bool inward = !port && (target == "knitting/numerics.h" ||
                                          (std::regex_match(target, to, layered) && allowed.at(folder).count(to[1].str())));
            if (!inward) wrong.push_back(name + " includes " + target);
        }
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files >= 40 && wrong.empty());
}

// The fibres' layers (src/fibres/ARCHITECTURE.md): each file includes only the fibres files below it, so
// it can be read knowing only those. fibre_record.h, strands_file.h and yarn_body include nothing of the folder;
// fur.cpp only strands_file.h, and the Fur workflow's files only fur and each other (look -> field -> loops). The headers include nothing of knitting (they name its types through pointers, so the
// component's users do not compile knitting); only the .cpp files may.
TEST(architecture_fibres_layers) {
    namespace fs = std::filesystem;
    const std::map<std::string, std::set<std::string>> allowed = {
        {"fibre_record", {}},
        {"strands_file", {}},
        {"yarn_body", {}},
        {"yarn_bake", {"yarn_body"}},
        {"fur", {"strands_file"}},
        {"fur_look", {"fur"}},
        {"fur_field", {"fur_look", "fur", "fibre_islands"}},
        {"fur_loops", {"fur", "fur_field", "fur_look", "strands_file"}},
        {"fibre_islands", {"fibre_record"}},
        {"fibre", {"fibre_record", "fibre_islands", "strands_file"}},
        {"fibre_shape", {"fibre", "fibre_record"}},
        {"fibre_checks", {"fibre", "fibre_record"}},
        {"fibre_strands", {"fibre", "fibre_record", "strands_file"}},
        {"fibre_generator", {"fibre", "fibre_record", "fibre_islands", "fibre_shape", "strands_file"}},
        {"fibre_budget", {"fibre_generator", "fibre", "fibre_record"}},
    };
    const fs::path root = fs::path(KNIT_SOURCE_DIR) / "src" / "fibres";
    const std::regex include(R"re(^\s*#\s*include\s+"([^"]+)")re"), ours(R"(fibres/(\w+)\.h)");
    int files = 0;
    std::vector<std::string> wrong;
    for (const fs::directory_entry& entry : fs::directory_iterator(root)) {
        if (!entry.is_regular_file() || entry.path().extension() == ".md" || entry.path().filename() == "CMakeLists.txt") continue;
        ++files;
        const std::string name = entry.path().filename().string(), stem = entry.path().stem().string();
        const bool header = entry.path().extension() == ".h";
        if (!allowed.count(stem)) {
            wrong.push_back(name + " sits outside the layers");
            continue;
        }
        std::ifstream in(entry.path());
        std::string line;
        std::smatch m, to;
        while (std::getline(in, line)) {
            if (!std::regex_search(line, m, include)) continue;
            const std::string target = m[1].str();
            if (target.rfind("model/", 0) == 0) continue;
            if (target.rfind("knitting/", 0) == 0) {
                if (header) wrong.push_back(name + " (a public header) includes " + target);
                continue;
            }
            const bool below = std::regex_match(target, to, ours) && (to[1].str() == stem || allowed.at(stem).count(to[1].str()));
            if (!below) wrong.push_back(name + " includes " + target);
        }
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files >= 20 && wrong.empty());
}

// Knitting computes; the modules above it read and write. The one exception is knit/pieces.cpp, which
// reads the pieces manifest and the pieces' meshes.
TEST(architecture_knitting_no_file_io) {
    namespace fs = std::filesystem;
    const fs::path root = fs::path(KNIT_SOURCE_DIR) / "src" / "knitting";
    const std::regex banned(R"(^\s*#\s*include\s*<(fstream|filesystem)>)");
    int files = 0;
    std::vector<std::string> wrong;
    for (const fs::directory_entry& entry : fs::recursive_directory_iterator(root)) {
        const std::string name = entry.path().lexically_relative(root).generic_string();
        if (!entry.is_regular_file() || name == "knit/pieces.cpp") continue;
        ++files;
        std::ifstream in(entry.path());
        std::string line;
        for (int number = 1; std::getline(in, line); ++number)
            if (std::regex_search(line, banned)) wrong.push_back(name + ":" + std::to_string(number));
    }
    for (const std::string& w : wrong) std::printf("    %s\n", w.c_str());
    CHECK(files >= 40 && wrong.empty());
}

