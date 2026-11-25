#pragma once
// Running Blender in the background on a script: the Export step's Blender format runs
// blender/ppt_package_to_blend.py (export/format_blend.cpp), the only use of Blender. This file finds Blender
// and supervises that process: a cancelled job stops it, and a run over its time cap fails, so nothing
// waits on Blender for ever. Windows only: elsewhere no Blender is found.
#include <functional>
#include <string>
#include <vector>

namespace knit {

// The blender.exe to run, "" if there is none. KNIT_BLENDER in the environment names it; otherwise the
// newest Blender 5.x (from 5.1, the version the scripts are written for) under Program Files.
std::string findBlender();

struct BlenderScriptRun {
    std::string blender;                             // findBlender()
    std::string script;                              // the .py file
    std::vector<std::string> arguments;              // what the script receives after "--"
    std::string logPath;                             // takes Blender's output
    int timeCapSeconds = 1800;
    std::function<void(const std::string&)> onLine;  // optional: each output line, as it appears
};

struct BlenderScriptResult {
    int exitCode = 1;
    std::string outputTail;  // the end of the output: a failing script prints its message there
};

// Runs the script to its end. false with *err if Blender did not start or overran its time cap. Throws
// Cancelled when the job is cancelled (model/cancel.h); Blender is stopped first.
bool runBlenderScript(const BlenderScriptRun& run, BlenderScriptResult& result, std::string* err);

} // namespace knit
