// Starting blender.exe in the background, stopping it on Cancel, and failing it after its time limit.
#include "io/blender_process.h"
#include "io/child_process.h"

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "model/cancel.h"
#include "io/files.h"
#include "model/env.h"

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace knit {

namespace fs = std::filesystem;

namespace {

constexpr int kPollMilliseconds = 100;
constexpr size_t kMessageTailBytes = 600;

#ifdef _WIN32

std::string outputTail(const std::string& logPath) {
    std::ifstream log(logPath, std::ios::binary);
    std::stringstream text;
    text << log.rdbuf();
    const std::string all = text.str();
    return all.size() > kMessageTailBytes ? all.substr(all.size() - kMessageTailBytes) : all;
}

// Hands on the lines Blender finished since `consumed`, the log's bytes already handed on.
void forwardNewLines(const BlenderScriptRun& run, std::streamoff& consumed) {
    if (!run.onLine) return;
    std::ifstream log(run.logPath, std::ios::binary);
    log.seekg(consumed);
    std::string line;
    while (std::getline(log, line) && !log.eof()) {  // at the end of the file the line is not finished yet
        consumed += static_cast<std::streamoff>(line.size()) + 1;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        run.onLine(line);
    }
}

#endif

} // namespace

#ifdef _WIN32

std::string findBlender() {
    std::error_code ec;
    if (const std::optional<std::string> chosen = envVar("KNIT_BLENDER")) return fs::exists(*chosen, ec) ? *chosen : std::string();
    const std::optional<std::string> programFiles = envVar("ProgramFiles");
    if (!programFiles) return "";
    std::string newest;
    int newestMinor = 0;
    for (const fs::directory_entry& folder : fs::directory_iterator(fs::path(*programFiles) / "Blender Foundation", ec)) {
        int major = 0, minor = 0;
        if (sscanf_s(folder.path().filename().string().c_str(), "Blender %d.%d", &major, &minor) != 2) continue;
        const fs::path blender = folder.path() / "blender.exe";
        if (major == 5 && minor >= 1 && minor > newestMinor && fs::exists(blender, ec)) {
            newest = blender.string();
            newestMinor = minor;
        }
    }
    return newest;
}

bool runBlenderScript(const BlenderScriptRun& run, BlenderScriptResult& result, std::string* err) {
    ChildProcessStart what;
    what.program = run.blender;
    what.arguments = {"-b", "--factory-startup", "--python-exit-code", "1", "--python", fs::absolute(run.script).string(), "--"};
    what.arguments.insert(what.arguments.end(), run.arguments.begin(), run.arguments.end());
    what.logPath = run.logPath;
    // Its destructor stops Blender, so leaving the job in any way (Cancelled included) frees the staging
    // directory before the publisher removes it.
    ChildProcess process;
    if (!process.start(what, err)) return false;
    const auto started = std::chrono::steady_clock::now();
    std::streamoff consumed = 0;
    while (!process.exitedWithin(kPollMilliseconds)) {
        cancelCheckpoint();
        forwardNewLines(run, consumed);
        const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - started);
        if (elapsed.count() > run.timeCapSeconds) return fail(err, "Blender did not finish within " + std::to_string(run.timeCapSeconds / 60) + " minutes");
    }
    forwardNewLines(run, consumed);
    result.exitCode = process.exitCode();
    result.outputTail = outputTail(run.logPath);
    return true;
}

#else

std::string findBlender() { return ""; }

bool runBlenderScript(const BlenderScriptRun&, BlenderScriptResult&, std::string* err) { return fail(err, "Blender runs on Windows only"); }

#endif

} // namespace knit
