// The Blender format. Blender itself writes the .blend: it runs in the background on
// blender/ppt_package_to_blend.py, which reads the photon_path_tracer package, so this format
// starts from a package folder and not from an ExportScene. A .blend is always metres, Z up,
// Blender's convention: the OBJ alone follows an import's units and up axis.
#include <filesystem>

#include "export/export_scene.h"
#include "io/blender_process.h"
#include "io/files.h"

namespace knit {

namespace fs = std::filesystem;

namespace {

constexpr const char* kConversionScript = "blender/ppt_package_to_blend.py";
// A 14 M segment export converts in about a minute; this only catches a Blender that hangs.
constexpr int kTimeCapSeconds = 1800;

} // namespace

bool blenderConversionAvailable() { return !findBlender().empty() && !findProjectFile(kConversionScript).empty(); }

bool writeBlendFile(const std::string& packageFolder, const FormatTarget& target, std::string& log, std::string* err) {
    const std::string blender = findBlender(), script = findProjectFile(kConversionScript);
    if (blender.empty()) return fail(err, "Blender 5.1 or later was not found (KNIT_BLENDER names its blender.exe)");
    if (script.empty()) return fail(err, std::string("the conversion script is missing: ") + kConversionScript);
    const fs::path outFolder = fs::absolute(target.folder), blendPath = outFolder / (target.baseName + ".blend");
    const fs::path logPath = outFolder / ".blender_output.txt";
    BlenderScriptRun run;
    run.blender = blender;
    run.script = script;
    run.arguments = {fs::absolute(packageFolder).string(), blendPath.string()};
    run.logPath = logPath.string();
    run.timeCapSeconds = kTimeCapSeconds;
    BlenderScriptResult result;
    if (!runBlenderScript(run, result, err)) return false;
    std::error_code ec;
    const bool written = result.exitCode == 0 && fs::exists(blendPath, ec);
    fs::remove(logPath, ec);
    if (!written) return fail(err, "Blender could not write the .blend:\n" + result.outputTail);
    log += "Blender: " + std::to_string(fs::file_size(blendPath, ec) / (1024 * 1024)) + " MB -> " + target.baseName + ".blend\n";
    return true;
}

} // namespace knit
