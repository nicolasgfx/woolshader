#pragma once
// The viewer's own preferences, kept between sessions: %APPDATA%/knit_viewer/viewer.ini (one
// "key = value" per line; unknown keys are ignored, so older and newer viewers share the file).
// Only explicit user choices are written: the UI scale, the window size and the folder of the last
// import. Settings files (--settings) and the export destination are not preferences.
#include <string>

namespace viewer {

struct ViewerPrefs {
    float uiScale = 1.0f;
    int windowW = 0, windowH = 0;                // 0: the default size
    std::string lastImportDir;
    std::string path;                            // where they are kept ("" : not kept)
};

ViewerPrefs loadViewerPrefs();
bool saveViewerPrefs(const ViewerPrefs& p);

// ImGui's layout file (window sizes, the workflow panel's width), kept next to the preferences.
// The first run that finds no layout file there copies an imgui.ini from the working directory
// (where earlier viewers kept it); an existing layout file is never overwritten. "" when the
// preferences are not kept (ImGui then uses its default, imgui.ini in the working directory).
std::string imguiIniPath(const ViewerPrefs& p);

} // namespace viewer
