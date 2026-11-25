// knit_viewer: the scene viewer. This file is the program: it opens the window and runs the frame
// loop, and frame() below is the order of everything that happens. src/viewer/ARCHITECTURE.md is
// the way into the folder.
#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#endif
#include <GL/gl.h>
#include <GLFW/glfw3.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl2.h"
#include "viewer/interface.h"
#include "jobs/scene_builder.h"
#include "jobs/scene_cache.h"
#include "model/coordinates.h"
#include "model/env.h"
#include "viewer/native_dialogs.h"
#include "viewer/renderer.h"
#include "viewer/view_content.h"
#include "jobs/crash_report.h"
#include "viewer/viewer_prefs.h"
#include "viewer/viewer_script.h"

using namespace viewer;

namespace {

// Everything a frame needs. Declared in runViewer, so its GL objects go before the window does.
struct Viewer {
    explicit Viewer(const ViewerOptions& o) : options(o), script(o) {}
    ViewerOptions options;
    SceneUi ui;
    ViewerPrefs prefs;
    CameraRig rig;
    CameraDrag drag;
    Renderer renderer;
    ScriptedRun script;
    // The frame readout (F3): the time between frames, smoothed.
    bool showFrameStats = false, f3Held = false;
    double frameMs = 0, lastFrameAt = 0;
};

// The display's content scale: the monitor's, or KNIT_CONTENT_SCALE (layout checks at another DPI
// without another monitor).
float contentScale(GLFWwindow* win) {
    static const std::optional<std::string> forced = knit::envVar("KNIT_CONTENT_SCALE");
    if (forced && std::atof(forced->c_str()) > 0) return static_cast<float>(std::atof(forced->c_str()));
    float sx = 1, sy = 1;
    glfwGetWindowContentScale(win, &sx, &sy);
    return std::max(sx, sy);
}

// The interface's request to the window, after it is drawn: a dialog blocks, then its answer goes to the action.
void answerPlatformRequest(Viewer& v, GLFWwindow* win) {
    const PlatformRequest request = v.ui.platform.request;
    v.ui.platform.request = PlatformRequest::None;
    switch (request) {
    case PlatformRequest::None: return;
    case PlatformRequest::OpenMesh:
        if (const std::optional<std::string> path = openMeshDialog(win, v.prefs.lastImportDir)) openMeshFile(v.ui, *path);
        return;
    case PlatformRequest::OpenMaskImage:
        if (const std::optional<std::string> path = openImageDialog(win, v.prefs.lastImportDir)) loadMaskImage(v.ui, *path);
        return;
    case PlatformRequest::ChooseExportFolder:
        if (const std::optional<std::string> dir = pickFolderDialog(win, exportFolderStart(v.ui))) chooseExportFolder(v.ui, *dir);
        return;
    case PlatformRequest::ShowExportFolder:
        if (v.ui.results.exported) openInFileManager(v.ui.results.exported->exportedTo);
        return;
    case PlatformRequest::SaveMarkings:
        if (const std::optional<std::string> path = saveMarkingsDialog(win, v.prefs.lastImportDir)) saveMarkingsAs(v.ui, *path);
        return;
    case PlatformRequest::OpenMarkings:
        if (const std::optional<std::string> path = openMarkingsDialog(win, v.prefs.lastImportDir)) loadMarkingsFrom(v.ui, *path);
        return;
    case PlatformRequest::SaveStitchMesh:
        if (const std::optional<std::string> path = saveStitchMeshDialog(win, v.prefs.lastImportDir, stitchMeshFileName(v.ui))) exportStitchMesh(v.ui, *path);
        return;
    }
}

// An interactive window opens maximized on the primary monitor's work area (the taskbar excluded).
// Its restored size (the kept one, else the created one) is first clamped to that area, frame
// included, and centred, so un-maximizing never leaves the window partly off screen.
void openMaximized(GLFWwindow* win, int restoredW, int restoredH) {
    if (GLFWmonitor* monitor = glfwGetPrimaryMonitor()) {
        int areaX = 0, areaY = 0, areaW = 0, areaH = 0;
        glfwGetMonitorWorkarea(monitor, &areaX, &areaY, &areaW, &areaH);
        int left = 0, top = 0, right = 0, bottom = 0;
        glfwGetWindowFrameSize(win, &left, &top, &right, &bottom);
        const int fitW = areaW - left - right, fitH = areaH - top - bottom;
        const int w = std::min(restoredW, fitW), h = std::min(restoredH, fitH);
        if (w > 0 && h > 0) {
            glfwSetWindowSize(win, w, h);
            glfwSetWindowPos(win, areaX + left + (fitW - w) / 2, areaY + top + (fitH - h) / 2);
        }
    }
    glfwMaximizeWindow(win);
}

// Lighting, preferences, the window's size and everything the window provides.
void setUpWindow(Viewer& v, GLFWwindow* win) {
    const ViewerOptions& o = v.options;
    if (o.plain >= 0) v.ui.view.unknittedParts = o.plain != 0;
    std::string studioError;
    v.ui.view.studioAvailable = v.renderer.init(&studioError);
    if (!v.ui.view.studioAvailable) {
        std::fprintf(stderr, "Studio lighting unavailable (Flat lighting used): %s\n", studioError.c_str());
        v.ui.view.studioNote = studioError;
    }
    v.ui.view.lighting = v.ui.view.studioAvailable && o.lighting != "flat" ? Lighting::Studio : Lighting::Flat;
    if (o.occlusion >= 0) v.ui.view.occlusion = static_cast<float>(std::min(o.occlusion, 1.0));
    v.prefs = loadViewerPrefs();
    v.ui.prefs = &v.prefs;
    static const std::string imguiIni = imguiIniPath(v.prefs);  // ImGui keeps the pointer past this viewer
    if (!imguiIni.empty()) ImGui::GetIO().IniFilename = imguiIni.c_str();
    uiScaleSetUser(o.uiScale > 0 ? o.uiScale : v.prefs.uiScale);
    // The window was created hidden: a shot keeps it so, an explicit size opens at that size, and
    // otherwise it opens maximized.
    if (o.width > 0 && o.height > 0) glfwSetWindowSize(win, o.width, o.height);
    if (o.shotPath.empty()) {
        if (o.width > 0 && o.height > 0) glfwShowWindow(win);
        else {
            int w = v.prefs.windowW, h = v.prefs.windowH;
            if (w <= 0 || h <= 0) glfwGetWindowSize(win, &w, &h);
            openMaximized(win, w, h);
            glfwShowWindow(win);
        }
    }
    v.ui.platform.dialogs = nativeDialogsAvailable();  // else the controls say how instead
    v.ui.platform.fileManager = true;
    glfwSetWindowUserPointer(win, &v);
    glfwSetDropCallback(win, [](GLFWwindow* w, int n, const char** paths) {
        Viewer* dropped = static_cast<Viewer*>(glfwGetWindowUserPointer(w));
        if (n > 0 && dropped && !dropped->ui.job.active) openMeshFile(dropped->ui, paths[0]);
    });
    glfwSetScrollCallback(win, ImGui_ImplGlfw_ScrollCallback);  // zoom reads ImGui's wheel
}

// The mouse in framebuffer pixels (high-DPI displays), as the viewport the pan scales to. A drag or
// the wheel acts on the camera when it starts inside the viewport (the interface's units = window
// coordinates) and not over the interface.
MouseState readMouse(GLFWwindow* win, const SceneUi& ui) {
    int width, height, winW, winH;
    glfwGetFramebufferSize(win, &width, &height);
    glfwGetWindowSize(win, &winW, &winH);
    const ImGuiIO& io = ImGui::GetIO();
    double wx, wy;
    glfwGetCursorPos(win, &wx, &wy);
    MouseState m;
    m.x = wx * (winW > 0 ? double(width) / winW : 1);
    m.y = wy * (winH > 0 ? double(height) / winH : 1);
    m.wheel = io.MouseWheel;
    for (int b = 0; b < 3; ++b) m.down[b] = glfwGetMouseButton(win, b) == GLFW_PRESS;
    m.shift = glfwGetKey(win, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS || glfwGetKey(win, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
    const bool inside = viewportRect(ui, io.DisplaySize.x, io.DisplaySize.y).contains(float(wx * io.DisplaySize.x / std::max(1, winW)), float(wy * io.DisplaySize.y / std::max(1, winH)));
    m.overViewport = inside && !io.WantCaptureMouse;
    return m;
}

// A click in framebuffer pixels as a pixel of the viewport: the click's y runs down, the viewport's origin is its lower corner.
ViewPixel viewPixelOf(const knit::Vec2& click, const FrameTarget& target) {
    const ViewportPixels& v = target.viewport;
    return {click.x - v.x, click.y - (target.height - v.y - v.height), double(v.width), double(v.height)};
}

FrameTarget frameTarget(const Viewer& v, GLFWwindow* win) {
    FrameTarget t;
    glfwGetFramebufferSize(win, &t.width, &t.height);
    const ImVec2 display = ImGui::GetIO().DisplaySize;
    t.viewport = viewportPixels(v.ui, t.width, t.height, display.x, display.y);
    t.cameraMoving = v.rig.moving(glfwGetTime());
    return t;
}

void beginImGuiFrame(GLFWwindow* win) {
    if (uiScaleBeginFrame(contentScale(win))) {
        ImGui_ImplOpenGL2_DestroyFontsTexture();
        ImGui_ImplOpenGL2_CreateFontsTexture();
    }
    ImGui_ImplOpenGL2_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
}

// F3 shows the frame time and the fibre segments drawn.
void updateFrameStats(Viewer& v, GLFWwindow* win) {
    const bool f3 = glfwGetKey(win, GLFW_KEY_F3) == GLFW_PRESS;
    if (f3 && !v.f3Held) v.showFrameStats = !v.showFrameStats;
    v.f3Held = f3;
    const double now = glfwGetTime();
    if (v.lastFrameAt > 0) v.frameMs = 0.9 * v.frameMs + 0.1 * 1000 * (now - v.lastFrameAt);
    v.lastFrameAt = now;
}

// --bench: the 3D view's cost per frame (median, each frame finished with glFinish, no vsync):
// frames with the camera at rest, then frames orbiting once around the target. Then what is resident on
// the GPU and how long its uploads took (once, when a view first showed it: no frame above pays for it).
void runBench(Viewer& v, FrameTarget target, int frames) {
    glfwSwapInterval(0);
    const ViewportFrame f = viewportFrame(v.ui);
    std::vector<double> restMs, orbitMs;
    long long restSegments = 0, orbitSegments = 0;
    for (int i = 0; i < frames; ++i) {
        const double t0 = glfwGetTime();
        v.renderer.draw(v.ui, f, v.rig.camera(), target);
        glFinish();
        restMs.push_back(1000 * (glfwGetTime() - t0));
        restSegments = v.renderer.fibreSegmentsDrawn();
    }
    const Camera kept = v.rig.camera();
    target.cameraMoving = true;
    for (int i = 0; i < frames; ++i) {
        v.rig.camera().orbit(2 * knit::kPi / frames / 0.006, 0);  // Camera::orbit turns 0.006 rad per pixel
        const double t0 = glfwGetTime();
        v.renderer.draw(v.ui, f, v.rig.camera(), target);
        glFinish();
        orbitMs.push_back(1000 * (glfwGetTime() - t0));
        orbitSegments = std::max(orbitSegments, v.renderer.fibreSegmentsDrawn());
    }
    v.rig.camera() = kept;
    std::sort(restMs.begin(), restMs.end());
    std::sort(orbitMs.begin(), orbitMs.end());
    std::printf("bench: %d frames %d x %d: rest %.2f ms (fibre segments %lld), orbit %.2f ms (max %lld), orbit worst %.2f ms; resident GPU "
                "%.1f MB, uploaded in %.1f ms%s\n",
                frames, target.width, target.height, restMs[restMs.size() / 2], restSegments, orbitMs[orbitMs.size() / 2], orbitSegments, orbitMs.back(),
                v.renderer.residentBytes() / 1048576.0, v.renderer.uploadMs(), f.source ? " (compared with the source)" : "");
    glfwSwapInterval(1);
}

// The finished frame, still in the back buffer, as a PNG. Returns the process's exit code.
int takeShot(Viewer& v, const FrameTarget& target) {
    if (v.options.benchFrames > 0) {
        runBench(v, target, v.options.benchFrames);
        FrameTarget atRest = target;
        atRest.cameraMoving = false;
        v.renderer.draw(v.ui, viewportFrame(v.ui), v.rig.camera(), atRest);  // the shot's frame again (the bench drew over it)
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
    }
    const bool written = writeBackBufferPng(v.options.shotPath, target.width, target.height);
    std::printf("%s %s (%d x %d)\n", written ? "wrote" : "FAILED to write", v.options.shotPath.c_str(), target.width, target.height);
    return v.script.reportShot(v.ui, written);
}

// One frame, in the order it happens. A value: the run ends with this exit code.
std::optional<int> frame(Viewer& v, GLFWwindow* win) {
    glfwPollEvents();
    pollSceneUi(v.ui);
    if (!v.script.beforeFrame(v.ui, glfwGetTime())) return 1;
    beginImGuiFrame(win);
    // The mouse acts with the layout of the last frame: this one's is not drawn yet.
    if (v.drag.update(v.rig.camera(), readMouse(win, v.ui), frameTarget(v, win).viewport.height)) v.rig.noteMoved(glfwGetTime());
    if (const std::optional<knit::Vec2> click = v.drag.takeClick()) {
        const ViewPixel pixel = viewPixelOf(*click, frameTarget(v, win));
        if (v.ui.stage == Stage::Mark) markClick(v.ui, v.rig.camera(), pixel, ImGui::GetIO().KeyCtrl ? ClickIntent::Select : ClickIntent::Place);
        else pickNode(v.ui, v.rig.camera(), pixel);
    }
    updateFrameStats(v, win);
    updateRenderPreview(v.ui, v.rig.camera(), v.rig.moving(glfwGetTime()), glfwGetTime());
    v.renderer.uploadRenderPreview(v.ui);
    const ImVec2 display = ImGui::GetIO().DisplaySize;  // the interface lays out in window coordinates
    drawSceneUi(v.ui, display.x, display.y);
    answerPlatformRequest(v, win);  // a dialog the interface asked for, before the viewport shows what it opened
    const ViewportFrame f = viewportFrame(v.ui);  // the shown content, and the source while comparing
    const bool hasContent = !f.shown.meshes.empty();
    v.rig.show(f.shown.viewKey, hasContent, f.bounds, v.options.startView);
    v.rig.apply(v.ui.camera, f.bounds, hasContent);  // after the interface: the content and its bounds are current
    v.ui.view.cameraPreset = presetOf(v.rig.camera());  // the View button's word next frame
    if (hasContent) v.script.frameShot(v.rig.camera(), f.bounds);
    const FrameTarget target = frameTarget(v, win);
    v.renderer.draw(v.ui, f, v.rig.camera(), target);
    if (v.showFrameStats) frameStatsLabel(v.ui, display.x, display.y, v.frameMs, v.renderer.fibreSegmentsDrawn());
    ImGui::Render();
    ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
    v.renderer.releaseUnused(v.ui);
    if (v.script.shotDue()) return takeShot(v, target);  // before the swap: the shot reads the back buffer
    glfwSwapBuffers(win);
    markFrameShown(v.ui);  // the quick check drawn in this frame answers its request (the edit latency)
    return std::nullopt;
}

int runViewer(GLFWwindow* win, const ViewerOptions& options) {
    Viewer v(options);
    setUpWindow(v, win);
    if (!loadSettings(v.ui, options)) return 2;
    v.ui.cacheRoot = knit::sceneCacheRoot();
    v.ui.loadKnit = options.loadKnit;
    v.ui.builder = std::make_unique<knit::SceneBuilder>();
    openObject(v.ui, options);
    while (!glfwWindowShouldClose(win))
        if (const std::optional<int> exitCode = frame(v, win)) return *exitCode;
    if (options.shotPath.empty()) {  // only a closed window keeps its size: a shot's is not the user's
        v.prefs.uiScale = uiScaleUser();
        // The restored size is kept; a maximized or minimized window's size is not that one.
        if (!glfwGetWindowAttrib(win, GLFW_MAXIMIZED) && !glfwGetWindowAttrib(win, GLFW_ICONIFIED))
            glfwGetWindowSize(win, &v.prefs.windowW, &v.prefs.windowH);
        saveViewerPrefs(v.prefs);
    }
    v.ui.builder.reset();
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    // A crash is reported to stderr and to crash.txt next to the viewer's preferences.
    const std::string prefs = loadViewerPrefs().path;
    knit::installCrashReport(prefs.empty() ? std::string() : (std::filesystem::path(prefs).parent_path() / "crash.txt").string());
    if (const std::string retired = retiredOption(argc, argv); !retired.empty()) {
        std::fprintf(stderr, "knit_viewer: %s is no longer available (the lab viewers were removed)\n", retired.c_str());
        return 2;
    }
    const ViewerOptions options = parseViewerOptions(argc, argv);
    if (!glfwInit()) { std::fprintf(stderr, "glfwInit failed\n"); return 1; }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);  // shown by setUpWindow once sized; --shot keeps it hidden
    GLFWwindow* win = glfwCreateWindow(1600, 1000, "Knit viewer", nullptr, nullptr);
    if (!win) { std::fprintf(stderr, "window creation failed\n"); glfwTerminate(); return 1; }
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(win, true);
    ImGui_ImplOpenGL2_Init();

    const int rc = runViewer(win, options);
    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(win);
    glfwTerminate();
    return rc;
}
