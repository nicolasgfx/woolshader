// The viewer's own preferences, viewer.ini, and the ImGui layout file (viewer_prefs.h).
#include "viewer/viewer_prefs.h"

#include "io/files.h"
#include "model/env.h"
#include "settings/file_presets.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace viewer {

namespace {

std::string prefsPath() {
    if (const std::optional<std::string> v = knit::envVar("KNIT_VIEWER_PREFS"); v && !v->empty()) return *v;  // tests, portable installs
#ifdef _MSC_VER
    const std::optional<std::string> home = knit::envVar("APPDATA");
    return home && !home->empty() ? (std::filesystem::path(*home) / "knit_viewer" / "viewer.ini").string() : "";
#else
    const std::optional<std::string> home = knit::envVar("HOME");
    return home ? *home + "/.config/knit_viewer/viewer.ini" : "";
#endif
}

} // namespace

ViewerPrefs loadViewerPrefs() {
    ViewerPrefs p;
    p.path = prefsPath();
    if (p.path.empty()) return p;
    std::ifstream f(p.path);
    std::string line;
    while (std::getline(f, line)) {
        const size_t eq = line.find('=');
        if (line.empty() || line[0] == '#' || eq == std::string::npos) continue;
        const std::string k = knit::trimmed(line.substr(0, eq)), v = knit::trimmed(line.substr(eq + 1));
        if (k == "ui_scale") p.uiScale = static_cast<float>(std::atof(v.c_str()));
        else if (k == "window_w") p.windowW = std::atoi(v.c_str());
        else if (k == "window_h") p.windowH = std::atoi(v.c_str());
        else if (k == "last_import_dir") p.lastImportDir = v;
    }
    if (!(p.uiScale >= 0.5f && p.uiScale <= 3.0f)) p.uiScale = 1.0f;
    return p;
}

bool saveViewerPrefs(const ViewerPrefs& p) {
    if (p.path.empty()) return false;
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(p.path).parent_path(), ec);
    std::ostringstream o;
    o << "# knit_viewer preferences (explicit choices only)\n";
    o << "ui_scale = " << p.uiScale << "\n";
    if (p.windowW > 0 && p.windowH > 0) o << "window_w = " << p.windowW << "\nwindow_h = " << p.windowH << "\n";
    if (!p.lastImportDir.empty()) o << "last_import_dir = " << p.lastImportDir << "\n";
    return knit::replaceFile(p.path, o.str(), nullptr);
}

std::string imguiIniPath(const ViewerPrefs& p) {
    namespace fs = std::filesystem;
    if (p.path.empty()) return "";
    const fs::path ini = fs::path(p.path).parent_path() / "imgui.ini";
    std::error_code ec;
    if (!fs::exists(ini, ec) && fs::is_regular_file("imgui.ini", ec)) {
        fs::create_directories(ini.parent_path(), ec);
        fs::copy_file("imgui.ini", ini, fs::copy_options::skip_existing, ec);
    }
    return ini.string();
}

} // namespace viewer
