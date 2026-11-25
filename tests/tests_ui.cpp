// UI tests of the scene viewer's panel (knit_tests_ui): no window. ImGui is built with
// IMGUI_ENABLE_TEST_ENGINE (tests/imgui_test_config.h), so every item registration reaches the
// hooks below. Two independent checks of widget IDs:
//   (a) registration: an ID registered twice in one frame is a collision (both items share
//       open state, activation and hover), whatever the mouse does;
//   (b) ImGui's own detector (ConfigDebugHighlightIdConflicts) only sees a collision through
//       the previous frame's HOVERED id (imgui.cpp: HoveredIdPreviousFrameItemCount), so the
//       mouse is moved over every registered item for a few frames.
// Every collapsible item is opened first (the collisions of sections that are closed by
// default only exist when they are expanded together).
#include <algorithm>
#include <atomic>
#include <functional>
#include <chrono>
#include <cmath>
#include <cstdarg>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <regex>
#include <set>
#include <stdexcept>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "model/coordinates.h"
#include "viewer/camera.h"
#include "viewer/fibre_lines.h"
#include "harness.h"
#include "io/cache_paths.h"
#include "io/files.h"
#include "io/marking_file.h"
#include "io/stitch_mesh_obj.h"
#include "imgui.h"
#include "imgui_internal.h"
#include "viewer/about.h"
#include "viewer/interface.h"
#include "viewer/mesh_lighting.h"
#include "viewer/view_content.h"
#include "export/export.h"
#include "settings/scene_registry.h"
#include "settings/scene_settings.h"
#include "settings/file_presets.h"
#include "settings/settings_profiles.h"
#include "settings/thickness_presets.h"
#include "settings/settings_mapping.h"
#include "viewer/viewer_prefs.h"
#include "viewer/viewer_script.h"

namespace {

struct ItemRec {
    ImGuiID id = 0;
    ImRect bb;
    ImGuiWindow* window = nullptr;
    std::string label, windowName;
    ImGuiItemStatusFlags flags = 0;
    bool disabled = false;  // drawn inside BeginDisabled (ImGuiItemFlags_Disabled)
};

std::vector<ItemRec> g_items;

} // namespace

// ---- test engine hooks (imgui_internal.h) ----
void ImGuiTestEngineHook_ItemAdd(ImGuiContext* ctx, ImGuiID id, const ImRect& bb, const ImGuiLastItemData*) {
    ItemRec r;
    r.id = id;
    r.bb = bb;
    r.window = ctx->CurrentWindow;
    if (r.window) r.windowName = r.window->Name;
    r.disabled = (ctx->CurrentItemFlags & ImGuiItemFlags_Disabled) != 0;
    g_items.push_back(r);
}
void ImGuiTestEngineHook_ItemInfo(ImGuiContext*, ImGuiID id, const char* label, ImGuiItemStatusFlags flags) {
    for (auto it = g_items.rbegin(); it != g_items.rend(); ++it)
        if (it->id == id) {
            if (label) it->label = label;
            it->flags = flags;
            return;
        }
}
void ImGuiTestEngineHook_Log(ImGuiContext*, const char*, ...) {}
const char* ImGuiTestEngine_FindItemDebugLabel(ImGuiContext*, ImGuiID) { return nullptr; }

namespace {

// The whole interface of a scene, without a window or GL, on a builder that runs no job: it records
// each job it receives and fails it (nothing here polls by itself, so the interface keeps showing
// the job as running). Tests read and set `ui` directly; this holds only what takes more than a
// line: completed results made without a job, a running job without a job, settings by key.
struct UiFixture {
    struct Jobs {
        std::mutex m;
        std::vector<std::string> received;
    };
    std::shared_ptr<Jobs> jobs = std::make_shared<Jobs>();
    viewer::SceneUi ui;  // after jobs: its builder (and thread) goes first

    explicit UiFixture(const std::string& scene) {
        const knit::SceneEntry* e = knit::findScene(scene);
        if (!e) throw std::runtime_error("unknown scene " + scene);
        ui.s = knit::viewerDefaults(*e);
        ui.workflowChosen = true;  // as a script or a settings file opens a scene: Knit, the settings' default
        ui.builder = std::make_unique<knit::SceneBuilder>([jobs = jobs](const knit::JobRequest& j, const knit::JobControl&, std::string* err) {
            std::string what = "knit";
            if (const knit::FibreRequest* f = std::get_if<knit::FibreRequest>(&j)) what = "fibres knit=" + std::to_string(f->knit ? f->knit->id : 0);
            if (const knit::ExportRequest* x = std::get_if<knit::ExportRequest>(&j)) what = "export dir=" + x->dir + " format=" + x->format;
            {
                std::lock_guard<std::mutex> lock(jobs->m);
                jobs->received.push_back(what);
            }
            if (err) *err = "test harness: no jobs";
            return std::shared_ptr<knit::BuiltScene>();
        });
        ui.exportTo.name = viewer::defaultExportName(ui);
        ui.exportTo.dir = (std::filesystem::temp_directory_path() / ui.exportTo.name).make_preferred().string();
        // The Export step's tests read the same on a machine with Blender and on one without.
        ui.exportTo.format = "ppt";
    }

    // Stage 2 of the Import step with an analysed file.
    bool openImport(const std::string& path) {
        viewer::openMeshFile(ui, path);
        for (int i = 0; i < 2000; ++i) {
            viewer::pollSceneUi(ui);
            if (ui.import.session.snapshot().phase == knit::ImportSession::Phase::Analysed) return true;
            if (ui.import.session.snapshot().phase == knit::ImportSession::Phase::Failed) return false;
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        return false;
    }

    // Runs the import to its end with stage 2's options, as the footer's action does (startImport).
    bool finishImport() {
        viewer::startImport(ui);
        for (int i = 0; i < 12000; ++i) {
            viewer::pollSceneUi(ui);
            if (!ui.import.fileOpen) return true;
            if (ui.import.session.snapshot().phase == knit::ImportSession::Phase::Failed) return false;
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        return false;
    }

    // A completed knit of the current settings, at Final quality (the full profile applied first)
    // or as they are, with `failingParts` of 4 parts failing.
    void injectKnit(bool final, int failingParts = 0) {
        const knit::SceneEntry* e = viewer::entryOf(ui);
        if (final && e) knit::applyFullProfile(ui.s, *e);
        auto k = std::make_shared<knit::KnitResult>();
        k->id = 1;
        k->settings = ui.s;
        k->final = e && knit::isFinalProfile(k->settings, *e);
        k->parts = 4;
        k->fail = failingParts;
        k->pass = k->parts - failingParts;
        for (int i = 0; i < failingParts; ++i) k->failing.push_back("part " + std::to_string(i));
        auto b = std::make_shared<knit::BuiltScene>();
        b->settings = k->settings;
        b->knit = k;
        knit::DisplayMesh yarn;  // 3 triangles of yarn
        yarn.kind = knit::DisplayMesh::Kind::Yarn;
        yarn.pos = {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0};
        yarn.idx = {0, 1, 2, 1, 3, 2, 0, 2, 3};
        b->meshes.push_back(yarn);
        ui.results.knit = b;
    }

    // Wool fibres of the shown knit for the current settings.
    std::shared_ptr<knit::BuiltScene> fakeFibres(std::vector<knit::ObjectSegments> objects = {{4321, 4321, 4321}}) const {
        auto f = std::make_shared<knit::FibreResult>();
        f->id = 1;
        f->knitId = viewer::knitResult(ui) ? viewer::knitResult(ui)->id : 0;
        f->settings = ui.s;
        f->key = knit::fibreResultKeyOf(f->settings);
        f->fibres = f->displayedFibres = 1234;
        f->segments = 4321;
        f->objectSegments = std::move(objects);
        auto b = std::make_shared<knit::BuiltScene>();
        b->kind = knit::BuiltScene::Kind::Fibres;
        b->settings = f->settings;
        b->fibreResult = f;
        if (ui.results.knit) b->meshes = ui.results.knit->meshes;
        b->fibres.pos.assign(6 * 10, 0.0f);  // 10 segments
        return b;
    }
    void injectFibres() { ui.results.fibres = fakeFibres(); }
    void injectFibres(const std::vector<knit::ObjectSegments>& objects) { ui.results.fibres = fakeFibres(objects); }

    void setKnitYarnLength(double mm) {
        if (!ui.results.knit) return;
        auto b = std::make_shared<knit::BuiltScene>(*ui.results.knit);
        b->yarnLength = mm;
        ui.results.knit = b;
    }

    // The Export step's detail levels by label ("Full", "Reduced", "Light").
    static const viewer::ExportDetail* exportDetailNamed(const std::string& label) {
        for (const viewer::ExportDetail& d : viewer::kExportDetails)
            if (label == d.label) return &d;
        return nullptr;
    }
    std::string exportDetailForecast(const std::string& label) const {
        const viewer::ExportDetail* d = exportDetailNamed(label);
        return d ? viewer::exportDetailForecast(ui, *d) : "no such level";
    }
    std::string exportDetailMemoryWarning(const std::string& label) const {
        const viewer::ExportDetail* d = exportDetailNamed(label);
        return d ? viewer::exportDetailMemoryWarning(ui, *d) : "no such level";
    }
    bool chooseExportDetail(const std::string& label) {
        const viewer::ExportDetail* d = exportDetailNamed(label);
        return d && viewer::applyExportDetail(ui, *d);
    }
    std::string exportDetailChosen() const {
        for (const viewer::ExportDetail& d : viewer::kExportDetails)
            if (viewer::exportDetailSelected(ui, d)) return d.label;
        return "";
    }
    // The settings a level owns: "<sides>/<ring>/<segment budget>".
    std::string exportDetailSettings() const {
        char b[96];
        std::snprintf(b, sizeof b, "%d/%g/%lld", ui.s.tubeSides, ui.s.tubeRing, ui.s.segmentBudget);
        return b;
    }

    // A running job without a job: the builder stays idle and nothing polls, so the state stays.
    // kind: 0 knit, 1 fibres, 2 export.
    void injectRunning(int kind, const std::string& stage, int done, int total, bool partial, bool cancelRequested = false) {
        ui.job.active = true;
        ui.job.kind = kind == 0 ? knit::BuiltScene::Kind::Knit : kind == 1 ? knit::BuiltScene::Kind::Fibres : knit::BuiltScene::Kind::Export;
        ui.job.meshOnly = false;
        ui.job.cancelRequested = cancelRequested;
        ui.job.progress = knit::JobProgress{};
        ui.job.progress.revision = 1;
        ui.job.progress.stage = stage;
        ui.job.progress.done = done;
        ui.job.progress.total = total;
        ui.job.progress.seconds = 400;
        ui.results.partial.reset();
        if (partial) {
            auto b = std::make_shared<knit::BuiltScene>();
            b->partial = true;
            if (ui.results.knit) b->meshes = ui.results.knit->meshes;
            ui.results.partial = b;
        }
    }

    // Settings by key (scene_settings.h), set as a settings file sets them; read back as written.
    bool setSetting(const std::string& key, const std::string& value) {
        std::string err;
        return knit::setSetting(ui.s, key, value, &err);
    }
    std::string setting(const std::string& key) const {
        const std::string text = "\n" + knit::writeSettings(ui.s), want = "\n" + key + " = ";
        const size_t at = text.find(want);
        if (at == std::string::npos) return "";
        const size_t from = at + want.size();
        return text.substr(from, text.find('\n', from) - from);
    }

    // A gauge preset matches the values (the panel shows it, not Custom).
    bool gaugePresetMatches() const {
        const knit::SceneEntry* e = viewer::entryOf(ui);
        return e && knit::matchingPreset(ui.s, *e) != nullptr;
    }
    bool finalProfile() const {
        const knit::SceneEntry* e = viewer::entryOf(ui);
        return e && knit::isFinalProfile(ui.s, *e);
    }
    // The matching mesh detail level's label ("" = custom).
    std::string meshDetail() const {
        const knit::FilePreset* l = knit::matchingFilePreset(knit::PresetGroup::StrandMesh, ui.s);
        return l ? l->label : "";
    }

    // The key the preview worker builds for after one poll (as each viewer frame polls).
    uint64_t fibrePreviewRequested() {
        viewer::pollSceneUi(ui);
        return ui.previews.fibres.requestedKey();
    }

    void setExportName(const std::string& name) {
        ui.exportTo.name = name;
        ui.exportTo.nameEdited = true;
    }
    // The folder the next "Choose folder ..." click picks.
    // The folder dialog's answer, as main.cpp gives it (answerPlatformRequest): the click's request, then the action.
    bool answerFolderDialog(const std::string& dir) {
        if (ui.platform.request != viewer::PlatformRequest::ChooseExportFolder) return false;
        ui.platform.request = viewer::PlatformRequest::None;
        viewer::chooseExportFolder(ui, dir);
        return true;
    }
    // A completed export job that published `path`, without a job.
    void injectExportDone(const std::string& path) {
        auto b = std::make_shared<knit::BuiltScene>();
        b->kind = knit::BuiltScene::Kind::Export;
        b->knit = ui.results.knit ? ui.results.knit->knit : nullptr;
        b->exportedTo = path;
        ui.results.exported = b;
        ui.results.exportAttempt = viewer::Attempt{};
        ui.results.exportAttempt.outcome = knit::JobOutcome::Done;
    }
    std::string exportSummary() const {
        const knit::ExportFormat* f = knit::findExportFormat(ui.exportTo.format);
        return f ? viewer::formatSummary(ui, *f) : "";
    }

    // The jobs the builder received, in order, waiting up to 2 s for at least `atLeast` of them.
    std::vector<std::string> submittedJobs(size_t atLeast) const {
        for (int i = 0; i < 200; ++i) {
            {
                std::lock_guard<std::mutex> lock(jobs->m);
                if (jobs->received.size() >= atLeast) return jobs->received;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        std::lock_guard<std::mutex> lock(jobs->m);
        return jobs->received;
    }
};

void setAttempt(viewer::Attempt& a, const std::string& message, bool failed) {
    a.outcome = failed ? knit::JobOutcome::Failed : knit::JobOutcome::Cancelled;
    a.message = message;
}

struct UiContext {
    ImGuiContext* ctx = nullptr;
    UiContext(float W, float H) {
        ctx = ImGui::CreateContext();
        ImGuiIO& io = ImGui::GetIO();
        io.IniFilename = nullptr;
        io.DisplaySize = ImVec2(W, H);
        io.DeltaTime = 1.0f / 60.0f;
        io.ConfigDebugHighlightIdConflicts = true;
        unsigned char* px = nullptr;
        int w = 0, h = 0;
        io.Fonts->AddFontDefault();
        io.Fonts->GetTexDataAsRGBA32(&px, &w, &h);
        ctx->TestEngineHookItems = true;
    }
    ~UiContext() { ImGui::DestroyContext(ctx); }
};

// One frame of the panel; the registered items of that frame are returned.
std::vector<ItemRec> frame(UiFixture& p, float W, float H) {
    g_items.clear();
    ImGui::NewFrame();
    viewer::drawSceneUi(p.ui, W, H);
    ImGui::Render();
    return g_items;
}

// Opens every openable item (collapsing headers, tree nodes) through the window state storage
// (where TreeNode keeps its open flag), frame after frame until nothing new opens.
void openEverything(UiFixture& p, float W, float H) {
    for (int round = 0; round < 32; ++round) {
        int opened = 0;
        for (const ItemRec& r : frame(p, W, H))
            if ((r.flags & ImGuiItemStatusFlags_Openable) && !(r.flags & ImGuiItemStatusFlags_Opened) && r.window) {
                r.window->StateStorage.SetInt(r.id, 1);
                ++opened;
            }
        if (opened == 0) return;
    }
    throw std::runtime_error("openEverything: items keep appearing after 32 rounds");
}

std::string describe(const ItemRec& r) { return "'" + (r.label.empty() ? std::string("?") : r.label) + "' in " + r.windowName; }

// (a) IDs registered more than once in one frame.
std::vector<std::string> duplicateIds(const std::vector<ItemRec>& items) {
    std::map<ImGuiID, std::vector<const ItemRec*>> byId;
    for (const ItemRec& r : items) byId[r.id].push_back(&r);
    std::vector<std::string> out;
    for (const auto& [id, rs] : byId)
        if (rs.size() > 1) {
            std::string s = "id " + std::to_string(id) + ":";
            for (const ItemRec* r : rs) s += " " + describe(*r);
            out.push_back(s);
        }
    return out;
}

// (b) ImGui's detector with the mouse over each registered item (a few frames each).
std::vector<std::string> hoverConflicts(UiFixture& p, float W, float H) {
    const std::vector<ItemRec> items = frame(p, W, H);
    std::set<ImGuiID> seen;
    std::vector<std::string> out;
    for (const ItemRec& r : items) {
        if (!seen.insert(r.id).second) continue;
        const ImVec2 c((r.bb.Min.x + r.bb.Max.x) * 0.5f, (r.bb.Min.y + r.bb.Max.y) * 0.5f);
        if (c.x < 0 || c.y < 0 || c.x >= W || c.y >= H) continue;  // off screen (the display is tall: none expected)
        ImGui::GetIO().MousePos = c;
        ImGuiID conflict = 0;
        for (int k = 0; k < 3; ++k) {
            frame(p, W, H);
            if (GImGui->DebugDrawIdConflicts != 0) conflict = GImGui->DebugDrawIdConflicts;
        }
        if (conflict) out.push_back("hovering " + describe(r) + ": ImGui reports conflicting id " + std::to_string(conflict));
    }
    ImGui::GetIO().MousePos = ImVec2(-FLT_MAX, -FLT_MAX);
    return out;
}

void report(const char* what, const std::vector<std::string>& lines) {
    for (const std::string& l : lines) std::printf("    %s: %s\n", what, l.c_str());
}

// One state of the interface: every section open, both checks.
bool clickSuffix(UiFixture& p, float W, float H, const std::string& suffix);  // below, with the other input helpers
void pressKey(UiFixture& p, float W, float H, ImGuiKey key);

void checkState(UiFixture& p, float W, float H, const char* what) {
    openEverything(p, W, H);
    const std::vector<ItemRec> items = frame(p, W, H);
    CHECK(items.size() > 20);  // the interface was drawn
    int openable = 0;
    for (const ItemRec& r : items) openable += (r.flags & ImGuiItemStatusFlags_Openable) != 0;
    std::printf("    %s at %.0f px: %zu items registered, %d collapsible (all open)\n", what, W, items.size(), openable);
    const std::vector<std::string> dup = duplicateIds(items), hov = hoverConflicts(p, W, H);
    if (!dup.empty() || !hov.empty()) std::printf("    in %s (%.0f x %.0f):\n", what, W, H);
    report("duplicate", dup);
    report("hover", hov);
    CHECK(dup.empty());
    CHECK(hov.empty());
}

// The whole interface of one scene: every stage, then the dialogs, every section open, at two
// window widths. A tall display keeps every item on screen (no scrolling needed for the hover sweep).
void checkScene(const char* scene) {
    namespace fs = std::filesystem;
    const fs::path obj = fs::temp_directory_path() / "knit_ui_test_import.obj";
    {
        std::ofstream f(obj);
        f << "v 0 0 0\nv 40 0 0\nv 40 30 0\nv 0 30 0\nv 0 0 0\nf 1 2 3\nf 1 3 4\nf 1 2 2\n";  // a duplicate vertex, a degenerate face
    }
    for (const float W : {1280.0f, 1920.0f}) {
        const float H = 12000;
        UiContext ui(W, H);
        UiFixture p(scene);
        const char* stages[] = {"Import", "Mark", "Knit", "Fibres", "Export"};
        for (int st = 0; st < viewer::kStages; ++st) {
            p.ui.stage = static_cast<viewer::Stage>(st);
            checkState(p, W, H, stages[st]);
        }
        // The Fibres step in both profiles and both views (wool / fur x reference / model).
        for (const bool fur : {false, true})
            for (const bool model : {false, true}) {
                viewer::chooseWorkflow(p.ui, fur ? knit::Workflow::Fur : knit::Workflow::Knit);
                p.ui.stage = viewer::Stage::Fibres;
                p.ui.view.fibres = model ? viewer::FibreView::Model : viewer::FibreView::Reference;
                const std::string what = std::string(fur ? "Fur" : "Wool") + (model ? ", model view" : ", reference view");
                checkState(p, W, H, what.c_str());
                if (!fur) {  // the population editor's other tab
                    p.ui.panels.fibreTabRequest = 1;
                    frame(p, W, H);  // the selection takes effect on the next frame
                    checkState(p, W, H, (what + ", stray fibres tab").c_str());
                    p.ui.panels.fibreTabRequest = 0;
                }
                p.ui.stage = viewer::Stage::Export;
                checkState(p, W, H, (std::string("Export, ") + (fur ? "fur" : "wool")).c_str());
            }
        viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
        // The toolbar's popovers, open one at a time, with a knit on show (the yarn's rows).
        p.injectKnit(false);
        p.injectFibres();
        p.ui.stage = viewer::Stage::Knit;
        p.ui.view.mode = viewer::ViewMode::Fibres;
        for (const char* menu : {"###viewMenu", "###displayMenu"}) {
            CHECK(clickSuffix(p, W, H, menu));
            checkState(p, W, H, menu);
            pressKey(p, W, H, ImGuiKey_Escape);
        }
        p.ui.results.knit.reset();
        p.ui.results.fibres.reset();
        p.ui.stage = viewer::Stage::Export;
        CHECK(p.openImport(obj.string()));
        checkState(p, W, H, "stage 2 of Import");
        viewer::closeImport(p.ui);
    }
}

} // namespace

// Knit or Fur (owner, 2026-11-19: no model is tied to a workflow): an example opens unchosen (the Import step asks,
// the later steps say "Choose Knit or Fur", the export is blocked with that fix); Knit shows Import -> Mark -> Knit ->
// Fibres -> Export, Fur Import -> Fur -> Export; the default folder says _knit / _fur; switching parks the other
// workflow's results and brings them back on the same source, never after a source change.
TEST(ui_unified_workflow_steps) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    auto stageButtons = [&](UiFixture& p) {
        int n = 0;
        for (const ItemRec& r : frame(p, W, H)) n += r.label == "##stage";
        return n;
    };
    UiFixture teddy("teddy");
    viewer::openExample(teddy.ui, "teddy-brown");
    std::printf("    opened: Import \"%s\", Knit \"%s\", Fibres \"%s\", Export \"%s\"\n", viewer::importState(teddy.ui).word.c_str(),
                viewer::knitState(teddy.ui).word.c_str(), viewer::fibreState(teddy.ui).word.c_str(), viewer::exportState(teddy.ui).word.c_str());
    CHECK(!teddy.ui.workflowChosen && teddy.ui.stage == viewer::Stage::Import && viewer::importState(teddy.ui).word == "Choose");
    CHECK(viewer::knitState(teddy.ui).word == viewer::kChooseWorkflow && viewer::fibreState(teddy.ui).word == viewer::kChooseWorkflow);
    CHECK(viewer::eligibility(teddy.ui).fix == viewer::kChooseWorkflow && !viewer::fibresInputReady(teddy.ui));
    CHECK(stageButtons(teddy) == viewer::kStages);
    {
        std::lock_guard<std::mutex> lock(teddy.jobs->m);
        CHECK(teddy.jobs->received.empty());  // nothing starts before the choice
    }
    viewer::chooseWorkflow(teddy.ui, knit::Workflow::Knit);
    CHECK(teddy.ui.workflowChosen && knit::knitWorkflow(teddy.ui.s) && !viewer::furShown(teddy.ui) && teddy.ui.stage == viewer::Stage::Knit);
    CHECK(teddy.ui.exportTo.name == "teddy-brown_knit" && stageButtons(teddy) == viewer::kStages);
    for (int i = 0; i < 400 && teddy.ui.job.active; ++i) {  // the first choice's knit (the harness refuses it at once)
        viewer::pollSceneUi(teddy.ui);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    teddy.injectKnit(true);
    const std::shared_ptr<const knit::BuiltScene> knitted = teddy.ui.results.knit;
    viewer::chooseWorkflow(teddy.ui, knit::Workflow::Fur);
    std::printf("    fur: strip %d steps, Fibres label \"%s\", Knit \"%s\", name %s\n", stageButtons(teddy),
                viewer::stageLabel(teddy.ui, viewer::Stage::Fibres), viewer::knitState(teddy.ui).word.c_str(), teddy.ui.exportTo.name.c_str());
    CHECK(viewer::furShown(teddy.ui) && teddy.ui.stage == viewer::Stage::Fibres && !teddy.ui.results.knit && teddy.ui.results.parked.knit == knitted);
    CHECK(stageButtons(teddy) == 3 && std::string(viewer::stageLabel(teddy.ui, viewer::Stage::Fibres)) == "1  Fur");
    CHECK(!viewer::stageInStrip(teddy.ui, viewer::Stage::Knit) && !viewer::stageInStrip(teddy.ui, viewer::Stage::Mark));
    CHECK(teddy.ui.exportTo.name == "teddy-brown_fur" && viewer::fibresInputReady(teddy.ui));
    viewer::chooseWorkflow(teddy.ui, knit::Workflow::Knit);  // the same source: the knit comes back, nothing starts
    CHECK(teddy.ui.results.knit == knitted && teddy.ui.exportTo.name == "teddy-brown_knit");
    {
        std::lock_guard<std::mutex> lock(teddy.jobs->m);
        CHECK(teddy.jobs->received.size() == 1);  // the first choice's knit only
    }
    // Another source while Fur is shown: the parked knit is not brought back.
    viewer::chooseWorkflow(teddy.ui, knit::Workflow::Fur);
    teddy.ui.s.importSurface = "wrapped";
    viewer::chooseWorkflow(teddy.ui, knit::Workflow::Knit);
    CHECK(!teddy.ui.results.knit);
}

// The export's blocker is visible where the user looks (viewer plan P3b): a Draft knit says in the
// Knit step that the export needs Final (state Draft, not Ready) and offers "Knit final" there
// (the Export step's fix is "Knit at Final quality"); the stage bar's Export word names what is missing; the Export step's fixed
// footer (below the settings, outside their scrolling child) holds the blocker's fix, and the
// export action once it is ready. Final without fibres is Ready (fibres are optional), Final with
// failing parts is Diagnostic, and a ticked fibre box without fibres asks for them.
TEST(ui_export_blocker_visible) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    auto itemY = [&](UiFixture& p, const std::string& label) {
        for (const ItemRec& r : frame(p, W, H))
            if (r.label == label) return r.bb.Min.y;
        return -1.0f;
    };
    auto inFooter = [&](UiFixture& p, const std::string& label) {  // the panel itself, not its scrolling child
        for (const ItemRec& r : frame(p, W, H))
            if (r.label == label) return r.windowName == "Workflow";
        return false;
    };
    for (const char* scene : {"teddy"}) {
        UiFixture p(scene);
        p.injectKnit(false);  // the viewer's default: Draft
        std::printf("    %s Draft: Knit \"%s\" (%s), Export \"%s\" (%s)\n", scene, viewer::knitState(p.ui).word.c_str(), viewer::knitState(p.ui).detail.c_str(),
                    viewer::exportState(p.ui).word.c_str(), viewer::exportState(p.ui).detail.c_str());
        CHECK(viewer::knitState(p.ui).word == "Draft" && viewer::knitState(p.ui).detail.find("export") != std::string::npos);
        CHECK((viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked) && viewer::exportState(p.ui).word == "Needs Final knit");
        p.ui.stage = viewer::Stage::Knit;
        CHECK(itemY(p, "Knit final") > 0);
        p.ui.stage = viewer::Stage::Export;
        // In the footer, below the Contents (the fibre box) and the folder button.
        const float fix = itemY(p, "Knit at Final quality"), fibres = itemY(p, "Include fibres"), folder = itemY(p, "Choose folder ...");
        std::printf("    Export panel: fix at y %.0f, Include fibres at %.0f, folder button at %.0f\n", fix, fibres, folder);
        CHECK(fix > 0 && fibres > 0 && folder > 0 && fix > fibres && fix > folder);
        CHECK(inFooter(p, "Knit at Final quality"));
        p.injectKnit(true);
        CHECK(viewer::knitState(p.ui).word == "Final" && !(viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked) && viewer::exportState(p.ui).word == "Ready");  // fibres are optional
        CHECK(itemY(p, "Knit at Final quality") < 0 && itemY(p, "Export package") > folder && inFooter(p, "Export package"));
        p.ui.exportTo.format = "obj";
        CHECK(itemY(p, "Export package") < 0 && inFooter(p, "Export OBJ"));
        p.ui.exportTo.format = "ppt";
        p.ui.s.fibres = true;
        CHECK((viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked) && viewer::exportState(p.ui).word == "Needs fibres");
        p.injectFibres();
        CHECK(!(viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked) && viewer::exportState(p.ui).word == "Ready");
        p.ui.s.fibres = false;
        p.injectKnit(true, 1);
        CHECK(!(viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked) && viewer::exportState(p.ui).word == "Diagnostic");
    }
}

// The canvas overlay's counts (viewer plan P4): the triangles shown, the fibres generated with
// their segments, and how many are shown (none in the Knit view); the Fibres step's reference
// patch is labelled as the patch, never as the model.
TEST(ui_overlay_counts) {
    UiContext ctx(1280, 12000);
    UiFixture p("teddy");
    p.injectKnit(true);
    p.injectFibres();
    p.ui.stage = viewer::Stage::Knit;
    p.ui.view.mode = viewer::ViewMode::Fibres;
    const std::string fib = viewer::viewCounts(p.ui, viewer::viewContent(p.ui));
    p.ui.view.mode = viewer::ViewMode::Knit;
    const std::string knit = viewer::viewCounts(p.ui, viewer::viewContent(p.ui));
    p.ui.stage = viewer::Stage::Fibres;
    p.ui.view.fibres = viewer::FibreView::Reference;
    const std::string ref = viewer::viewCounts(p.ui, viewer::viewContent(p.ui));
    std::printf("    Fibres view: %s\n    Knit view: %s\n    reference: %s\n", fib.c_str(), knit.c_str(), ref.c_str());
    CHECK(fib.rfind("Shown: 3 triangles", 0) == 0);
    CHECK(fib.find("1,234 fibres generated (4,321 segments), 1,234 shown") != std::string::npos);
    CHECK(knit.rfind("Shown: 3 triangles", 0) == 0 && knit.find("not shown in this view") != std::string::npos);
    CHECK(ref.rfind("Reference patch: ", 0) == 0 && ref.find("generated") == std::string::npos);
}

// The Export step's detail levels (GPU memory in the renderer): each sets the strand mesh level and
// the segment budget together, and forecasts the triangles (2 x sides x yarn length / ring spacing)
// and the fibre segments (the budget, or the generated count if fewer).
TEST(ui_export_detail) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(true);
    p.ui.s.fibres = true;
    p.ui.stage = viewer::Stage::Export;
    CHECK(p.exportDetailForecast("Full") == "triangles: knit first, <= 14.0 M fibre segments per object");
    p.setKnitYarnLength(390000);  // 390 m of tubed curves
    const std::string full = p.exportDetailForecast("Full"), reduced = p.exportDetailForecast("Reduced"),
                      light = p.exportDetailForecast("Light");
    std::printf("    Full: %s; Reduced: %s; Light: %s\n", full.c_str(), reduced.c_str(), light.c_str());
    CHECK(full == "~15.6 M triangles, <= 14.0 M fibre segments per object");    // 2 x 8 x 390000 / 0.4
    CHECK(reduced == "~7.8 M triangles, <= 7.0 M fibre segments per object");    // 2 x 6 x 390000 / 0.6
    CHECK(light.rfind("~5.", 0) == 0 && light.find("<= 5.0 M fibre segments per object") != std::string::npos);  // 5.85 M
    CHECK(p.exportDetailChosen() == "Reduced" && p.exportDetailSettings() == "6/0.6/7000000");  // the defaults
    CHECK(!p.chooseExportDetail("Reduced"));  // already chosen: nothing changes
    CHECK(p.chooseExportDetail("Full") && p.exportDetailChosen() == "Full" && p.exportDetailSettings() == "8/0.4/14000000");
    CHECK(p.chooseExportDetail("Light") && p.exportDetailSettings() == "6/0.8/5000000");
    // Generated fibres give the segments and the renderer's memory (the harness's fibres: 4,321 segments,
    // grown at Light's budget; at Full's they are scaled, "~"). 15.6 M triangles: 1256 + 63.7 x 15.6 MiB.
    p.injectFibres();
    CHECK(p.exportDetailForecast("Full") == "~15.6 M triangles, ~4 k fibre segments, ~2.2 GiB GPU while rendering");
    p.ui.s.fibres = false;
    CHECK(p.exportDetailForecast("Full") == "~15.6 M triangles, no fibres, ~2.2 GiB GPU while rendering");
    p.ui.s.fibres = true;
    // The white teddy's bear at Reduced: the exported count at the result's own budget, exact; the late
    // taper's vertices (12.88 M for 7 M kept) scale the forecast at Full's budget.
    CHECK(p.chooseExportDetail("Reduced"));
    p.injectFibres({{18464844, 7000000, 12880000}});
    CHECK(p.exportDetailForecast("Reduced") == "~7.8 M triangles, 12.9 M fibre segments, ~4.1 GiB GPU while rendering");
    CHECK(p.exportDetailForecast("Full").find(", ~25.8 M fibre segments") != std::string::npos);  // 14 M x 12.88 / 7
    CHECK(p.exportDetailMemoryWarning("Reduced").empty());
    // Two bears of 34 M exported segments each: about 14.5 GiB, above the 9 GiB cap; one bear fits.
    p.injectFibres({{18464844, 18464844, 33970516}, {18464844, 18464844, 33970516}});
    std::printf("    two bears: %s; %s\n", p.exportDetailForecast("Reduced").c_str(), p.exportDetailMemoryWarning("Reduced").c_str());
    CHECK(p.exportDetailMemoryWarning("Reduced").find("above the 9 GiB") != std::string::npos);
    p.injectFibres({{18464844, 18464844, 33970516}});
    CHECK(p.exportDetailMemoryWarning("Reduced").empty());
    // A changed late taper invalidates the exported count.
    p.ui.s.tailTaperTip += 0.1;
    CHECK(p.exportDetailForecast("Reduced").find("export to count (the late taper changed)") != std::string::npos);
    CHECK(p.exportDetailMemoryWarning("Reduced").empty());
    p.ui.s.tailTaperTip -= 0.1;
    // Plies: the tubed length is the plies' already (their helix included): not multiplied by the plies.
    p.ui.s.plies = 3;
    const knit::FilePreset* reducedLevel = viewer::meshLevelOf(*p.exportDetailNamed("Reduced"));
    CHECK(reducedLevel != nullptr);
    if (reducedLevel) {
        char want[48];
        std::snprintf(want, sizeof want, "~%.1f M triangles", 2.0 * reducedLevel->values.plyTubeSides / reducedLevel->values.plyTubeRing * 390000 / 1e6);
        CHECK(p.exportDetailForecast("Reduced").rfind(want, 0) == 0);
    }
    p.ui.s.plies = 1;
    CHECK(p.chooseExportDetail("Full") && p.exportDetailSettings() == "8/0.4/14000000");
}

// ImGui's layout file lives next to the viewer preferences. Its first use copies an imgui.ini of
// the working directory (the old location; the copy keeps the panel width), the original stays,
// and an existing layout file is never overwritten.
TEST(ui_imgui_ini_migration) {
    namespace fs = std::filesystem;
    const fs::path root = fs::temp_directory_path() / "knit_ui_imgui_ini";
    std::error_code ec;
    fs::remove_all(root, ec);
    fs::create_directories(root / "work");
    auto write = [](const fs::path& p, const std::string& text) { std::ofstream(p, std::ios::binary) << text; };
    auto read = [](const fs::path& p) {
        std::ifstream f(p, std::ios::binary);
        return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
    };
    write(root / "work" / "imgui.ini", "[Window][Workflow]\nSize=512,900\n");
    const fs::path cwd = fs::current_path();
    fs::current_path(root / "work");
    viewer::ViewerPrefs prefs;
    prefs.path = (root / "prefs" / "viewer.ini").string();
    const std::string ini = viewer::imguiIniPath(prefs);
    CHECK(fs::path(ini) == root / "prefs" / "imgui.ini");
    CHECK(read(ini) == "[Window][Workflow]\nSize=512,900\n");                // copied
    CHECK(read(root / "work" / "imgui.ini") == "[Window][Workflow]\nSize=512,900\n");  // not moved
    write(root / "work" / "imgui.ini", "[Window][Workflow]\nSize=300,900\n");
    CHECK(viewer::imguiIniPath(prefs) == ini);
    CHECK(read(ini) == "[Window][Workflow]\nSize=512,900\n");                // never overwritten
    viewer::ViewerPrefs unkept;
    CHECK(viewer::imguiIniPath(unkept).empty());
    fs::current_path(cwd);
    fs::remove_all(root, ec);
}


TEST(ui_panel_unique_ids_teddy) { checkScene("teddy"); }
TEST(ui_panel_unique_ids_wool_swatch) { checkScene("wool-swatch"); }

namespace {

// The display at another interface scale, as the viewer builds it (widgets.cpp): the theme's
// logical sizes, then every size and the font scaled; uiPx follows. Restored on destruction.
struct ScaledUi {
    UiContext ctx;
    explicit ScaledUi(float W, float H, float scale) : ctx(W, H) {
        ImGuiStyle& st = ImGui::GetStyle();
        st = ImGuiStyle();
        ImGui::StyleColorsDark(&st);
        viewer::applySceneTheme(st);
        st.ScaleAllSizes(scale);
        ImGuiIO& io = ImGui::GetIO();
        io.Fonts->Clear();
        ImFontConfig cfg;
        cfg.SizePixels = 13.0f * scale;
        io.Fonts->AddFontDefault(&cfg);
        unsigned char* px = nullptr;
        int w = 0, h = 0;
        io.Fonts->GetTexDataAsRGBA32(&px, &w, &h);
        viewer::uiScaleSetUser(scale);
    }
    ~ScaledUi() { viewer::uiScaleSetUser(1.0f); }
};

const ItemRec* findItem(const std::vector<ItemRec>& items, const std::string& label) {
    for (const ItemRec& r : items)
        if (r.label == label) return &r;
    return nullptr;
}

// An item whose label starts with `prefix` (labels that carry a value, e.g. "Stray fibres 20%###stray").
const ItemRec* findPrefix(const std::vector<ItemRec>& items, const std::string& prefix) {
    for (const ItemRec& r : items)
        if (r.label.rfind(prefix, 0) == 0) return &r;
    return nullptr;
}

// A message of n words (it wraps over many lines).
std::string longMessage(int words) {
    std::string m;
    for (int i = 0; i < words; ++i) m += "word" + std::to_string(i) + " ";
    return m;
}

int countItems(const std::vector<ItemRec>& items, const std::string& label) {
    int n = 0;
    for (const ItemRec& r : items) n += r.label == label;
    return n;
}

// A few frames, so that measured sizes (auto-resizing windows, the Knit footer) settle.
std::vector<ItemRec> settle(UiFixture& p, float W, float H) {
    for (int i = 0; i < 4; ++i) frame(p, W, H);
    return frame(p, W, H);
}

std::string renderedText(UiFixture& p, float W, float H);  // below, with the snapshot: every node open, the frame's text

// A click at a point of the display, through ImGui's event queue (one transition per frame).
void clickAt(UiFixture& p, float W, float H, ImVec2 c) {
    ImGuiIO& io = ImGui::GetIO();
    io.AddMousePosEvent(c.x, c.y);
    frame(p, W, H);
    io.AddMouseButtonEvent(0, true);
    frame(p, W, H);
    io.AddMouseButtonEvent(0, false);
    frame(p, W, H);
    io.AddMousePosEvent(-FLT_MAX, -FLT_MAX);
    frame(p, W, H);
}

ImVec2 centreOf(const ItemRec& r) { return ImVec2((r.bb.Min.x + r.bb.Max.x) * 0.5f, (r.bb.Min.y + r.bb.Max.y) * 0.5f); }

// A click on an item's centre (the item labelled `label`, or the first whose label starts with it).
bool click(UiFixture& p, float W, float H, const std::string& label, bool prefix = false) {
    const std::vector<ItemRec> items = settle(p, W, H);
    const ItemRec* r = prefix ? findPrefix(items, label) : findItem(items, label);
    if (!r) return false;
    clickAt(p, W, H, centreOf(*r));
    return true;
}

// The item whose label ends with `suffix`: a toolbar menu button, "<shown>###viewMenu".
const ItemRec* findSuffix(const std::vector<ItemRec>& items, const std::string& suffix) {
    for (const ItemRec& r : items)
        if (r.label.size() >= suffix.size() && r.label.compare(r.label.size() - suffix.size(), suffix.size(), suffix) == 0) return &r;
    return nullptr;
}

bool clickSuffix(UiFixture& p, float W, float H, const std::string& suffix) {
    const std::vector<ItemRec> items = settle(p, W, H);
    const ItemRec* r = findSuffix(items, suffix);
    if (!r) return false;
    clickAt(p, W, H, centreOf(*r));
    return true;
}

// A key pressed and released, one transition per frame.
void pressKey(UiFixture& p, float W, float H, ImGuiKey key) {
    ImGuiIO& io = ImGui::GetIO();
    io.AddKeyEvent(key, true);
    frame(p, W, H);
    io.AddKeyEvent(key, false);
    frame(p, W, H);
}

// The item labelled `label` inside an open popover (ImGui names a popup's window "##Popup_<id>").
const ItemRec* popoverItem(const std::vector<ItemRec>& items, const std::string& label) {
    for (const ItemRec& r : items)
        if (r.windowName.rfind("##Popup", 0) == 0 && r.label == label) return &r;
    return nullptr;
}

// The labels of both toolbar popovers, each opened by its button and closed with Escape.
std::set<std::string> popoverLabels(UiFixture& p, float W, float H) {
    std::set<std::string> l;
    for (const char* menu : {"###viewMenu", "###displayMenu"}) {
        if (!clickSuffix(p, W, H, menu)) continue;
        for (const ItemRec& r : settle(p, W, H))
            if (r.windowName.rfind("##Popup", 0) == 0) l.insert(r.label);
        pressKey(p, W, H, ImGuiKey_Escape);
    }
    return l;
}

// A mouse drag of an item by dx pixels.
bool drag(UiFixture& p, float W, float H, const std::string& label, float dx) {
    const std::vector<ItemRec> items = settle(p, W, H);
    const ItemRec* r = findItem(items, label);
    if (!r) return false;
    const ImVec2 c((r->bb.Min.x + r->bb.Max.x) * 0.5f, (r->bb.Min.y + r->bb.Max.y) * 0.5f);
    ImGuiIO& io = ImGui::GetIO();
    io.AddMousePosEvent(c.x, c.y);
    frame(p, W, H);
    io.AddMouseButtonEvent(0, true);
    frame(p, W, H);
    for (int i = 1; i <= 10; ++i) {
        io.AddMousePosEvent(c.x + dx * i / 10, c.y);
        frame(p, W, H);
    }
    io.AddMouseButtonEvent(0, false);
    frame(p, W, H);
    io.AddMousePosEvent(-FLT_MAX, -FLT_MAX);
    frame(p, W, H);
    return true;
}

// The item labelled `label` lies inside the panel ("Workflow", the root of its window) and the display.
bool insidePanel(UiFixture& p, float W, float H, const std::string& label) {
    const std::vector<ItemRec> items = settle(p, W, H);
    const ItemRec* r = findItem(items, label);
    if (!r) {
        std::printf("    \"%s\": not drawn\n", label.c_str());
        return false;
    }
    const ImGuiWindow* win = r->window;
    while (win && win->ParentWindow) win = win->ParentWindow;
    const ImRect panel(win->Pos, ImVec2(win->Pos.x + win->Size.x, win->Pos.y + win->Size.y));
    std::printf("    \"%s\" at y %.0f..%.0f in %s y %.0f..%.0f\n", label.c_str(), r->bb.Min.y, r->bb.Max.y, win->Name, panel.Min.y, panel.Max.y);
    return r->bb.Min.y >= panel.Min.y && r->bb.Max.y <= panel.Max.y && r->bb.Max.y <= H;
}

// Every item's label in one frame with every section open.
std::set<std::string> allLabels(UiFixture& p, float W, float H) {
    openEverything(p, W, H);
    std::set<std::string> l;
    for (const ItemRec& r : frame(p, W, H)) l.insert(r.label);
    return l;
}

std::string missingLabels(const std::set<std::string>& have, std::initializer_list<const char*> want) {
    std::string m;
    for (const char* w : want)
        if (!have.count(w)) m += std::string(" '") + w + "'";
    return m;
}

} // namespace

// The one progress rule (workflow.h progressView) on the knit job's real step names
// (scene_builder.cpp): the counters are parts only while knitting; other steps get a bar without
// unit words, or an indeterminate one without a total; nothing is a time estimate.
TEST(ui_progress_classification) {
    UiContext ctx(1280, 720);
    auto view = [](const char* stage, int done, int total, bool partial = false) {
        knit::JobProgress p;
        p.stage = stage;
        p.done = done;
        p.total = total;
        p.seconds = 400;
        viewer::RunningJob job;
        job.progress = p;
        return viewer::progressView(job, partial);
    };
    const viewer::JobProgressView knitting = view("knitting white_bear_fabric_c4", 7, 18, true);
    std::printf("    \"%s\" | \"%s\" | \"%s\" | %.3f | \"%s\"\n", knitting.title.c_str(), knitting.count.c_str(), knitting.qualifier.c_str(),
                knitting.fraction, knitting.status.c_str());
    CHECK(knitting.title == "Knitting ..." && knitting.count == "7 of 18 parts");
    CHECK(std::fabs(knitting.fraction - 7.0f / 18.0f) < 1e-6f);
    CHECK(knitting.qualifier.rfind("Partial build shown", 0) == 0);
    CHECK(knitting.status.rfind("Knitting white_bear_fabric_c4", 0) == 0 && knitting.status.find("Elapsed 6 min 40 s") != std::string::npos);
    CHECK(view("knitting", 0, 18).count == "0 of 18 parts");
    CHECK(view("knitted brown_bear_fabric_c8", 18, 18).count == "18 of 18 parts");
    // Not parts: the strand mesh rebuild's counter, the display, an unknown step with a counter.
    const viewer::JobProgressView rebuild = view("rebuilding the strand mesh", 2, 5);
    CHECK(rebuild.count.empty() && std::fabs(rebuild.fraction - 0.4f) < 1e-6f && rebuild.qualifier == "Rebuilding the strand mesh");
    const viewer::JobProgressView display = view("building the display", 0, 0);
    CHECK(display.count.empty() && display.fraction < 0 && display.qualifier == "Building the display");
    const viewer::JobProgressView unknown = view("sorting the rows", 3, 9);
    CHECK(unknown.count.empty() && std::fabs(unknown.fraction - 1.0f / 3.0f) < 1e-6f && unknown.qualifier == "Sorting the rows");
    // A mesh-only rebuild never counts parts; a cancel request says what it waits for.
    knit::JobProgress p;
    p.stage = "knitting x";
    p.done = 1;
    p.total = 2;
    viewer::RunningJob job;
    job.progress = p;
    job.meshOnly = true;
    CHECK(viewer::progressView(job, false).count.empty());
    CHECK(viewer::progressView(job, false).title == "Meshing ...");
    job.meshOnly = false;
    job.cancelRequested = true;
    const viewer::JobProgressView stopping = viewer::progressView(job, false);
    CHECK(stopping.title == "Stopping ..." && stopping.qualifier.find("stopping after knitting x") != std::string::npos);
}

// The knit job's state has one home (the Knit footer): there exactly one Cancel and no second one
// in the status bar; on another stage the status bar carries the Cancel (the job stays
// cancellable from everywhere).
TEST(ui_knit_progress_single_home) {
    const float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(false);
    p.injectRunning(0, "knitting white_bear_fabric_c4", 7, 18, true);
    p.ui.stage = viewer::Stage::Knit;
    std::vector<ItemRec> items = settle(p, W, H);
    int inPanel = 0, inStatus = 0;
    for (const ItemRec& r : items)
        if (r.label == "Cancel") (r.windowName.rfind("Workflow", 0) == 0 ? inPanel : inStatus) += 1;
    std::printf("    Knit stage: Cancel in the panel %d, elsewhere %d\n", inPanel, inStatus);
    CHECK(inPanel == 1 && inStatus == 0);
    p.ui.stage = viewer::Stage::Export;
    items = settle(p, W, H);
    inPanel = inStatus = 0;
    for (const ItemRec& r : items)
        if (r.label == "Cancel") (r.windowName == "##status" ? inStatus : inPanel) += 1;
    std::printf("    Export stage: Cancel in the status bar %d, elsewhere %d\n", inStatus, inPanel);
    CHECK(inStatus == 1 && inPanel == 0);
}

namespace {

// The Knit footer's two actions (or Cancel) on a window of this size at 1.5 x interface scale, whatever
// the status block holds: true when every one of them is inside the panel.
bool knitFooterActionsReachable(float W, float H) {
    ScaledUi ctx(W, H, 1.5f);
    bool all = true;
    auto inside = [&](UiFixture& p, const std::string& label) {
        const std::vector<ItemRec> items = settle(p, W, H);
        const ItemRec* r = findItem(items, label);
        if (!r) {
            std::printf("    \"%s\": not drawn\n", label.c_str());
            return false;
        }
        const ImGuiWindow* win = r->window;
        while (win && win->ParentWindow) win = win->ParentWindow;  // the panel ("Workflow")
        const ImRect panel(win->Pos, ImVec2(win->Pos.x + win->Size.x, win->Pos.y + win->Size.y));
        std::printf("    \"%s\" at y %.0f..%.0f in %s y %.0f..%.0f\n", label.c_str(), r->bb.Min.y, r->bb.Max.y, win->Name, panel.Min.y, panel.Max.y);
        return r->bb.Min.y >= panel.Min.y && r->bb.Max.y <= panel.Max.y && r->bb.Max.y <= H;
    };
    std::printf("    window %.0f x %.0f at 1.5 x\n", W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Knit;
    p.injectKnit(false);
    p.injectRunning(0, "knitting white_bear_fabric_c4", 7, 18, true);
    all &= inside(p, "Cancel");
    UiFixture q("teddy");
    q.ui.stage = viewer::Stage::Knit;
    all &= q.setSetting("strand.diameter", "1.46");  // before the knit: a current Draft result of thick yarn
    q.injectKnit(false);
    setAttempt(q.ui.results.knitAttempt, longMessage(150), true);
    all &= inside(q, "Knit draft") && inside(q, "Knit final");
    all &= q.setSetting("strand.wale", "5.4");  // the knit is now outdated: the status says previous result
    all &= inside(q, "Knit draft") && inside(q, "Knit final");
    UiFixture r("teddy");
    r.ui.stage = viewer::Stage::Knit;
    all &= r.setSetting("strand.diameter", "1.46");
    setAttempt(r.ui.results.knitAttempt, longMessage(100), false);
    all &= inside(r, "Knit draft") && inside(r, "Knit final");
    return all;
}

} // namespace

// The Knit footer's actions stay on screen on a short window and on a narrow one at 1.5 x interface
// scale, whatever the status block holds (running with a partial build; a long failure, thick yarn and
// a Draft result; an outdated result): the status block scrolls, the buttons keep their room.
TEST(ui_knit_footer_reachable) {
    CHECK(knitFooterActionsReachable(1280, 720));
    CHECK(knitFooterActionsReachable(1100, 720));
}

// Every former Knit-step control is still drawn (presence only; the bindings are the interaction
// tests below), with every section open, in both yarn constructions and with the custom mesh.
TEST(ui_knit_controls_inventory) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    auto labels = [&](UiFixture& p) {
        openEverything(p, W, H);
        std::set<std::string> l;
        for (const ItemRec& r : frame(p, W, H)) l.insert(r.label);
        return l;
    };
    auto missing = [](const std::set<std::string>& have, std::initializer_list<const char*> want) {
        std::string m;
        for (const char* w : want)
            if (!have.count(w)) m += std::string(" '") + w + "'";
        return m;
    };
    // The refresh (2026-11-11): the Draft / Final toggle became the footer's Knit draft / Knit final, the
    // hand-made preset list the character tiles, and the detail sections are collapsible nodes.
    const std::initializer_list<const char*> common = {"Stockinette", "Reverse stockinette", "Garter", "Rib 1 x 1",
                                                       "##Gauge", "##Yarn diameter", "##Stitch width",
                                                       "##Row height", "Single tube", "Ply tubes", "Fibre strands", "Fibre maps", "##Depth",
                                                       "##Clearance", "Contact flattening",
                                                       "Flatten yarn that is too thick for the stitch", "##Ply texture", "Swept plies",
                                                       "Circles", "Knit draft", "Knit final", "##Yarn tessellation", "Off", "Natural",
                                                       "Hand-made", "Rustic", "##Uneven stitches", "##Column shift", "##Column run",
                                                       "##Raised loops", "##Loop height", "##Thick and thin spots", "##Spot strength",
                                                       "##Spot length", "##Fabric wave", "##Variation", "##Seams", "Gauge & yarn geometry",
                                                       "Fine-tune character", "Stitch layout", "Mesh detail",
                                                       "Source", "Knit", "Fibres", "Fit"};
    for (const char* scene : {"teddy"}) {
        UiFixture p(scene);
        p.ui.stage = viewer::Stage::Knit;
        CHECK(p.setSetting("strand.flatten", "true"));
        CHECK(p.setSetting("strand.tube_sides", "5"));  // no level matches: the custom fields
        std::set<std::string> have = labels(p);
        std::string m = missing(have, common) + missing(have, {"##Inflation", "##Sides around", "##Ring spacing"});
        for (const char* gone : {"Draft", "Final", "##Hand-made preset", "Solver & tessellation details"})
            if (have.count(gone)) m += std::string(" (replaced control ") + gone + ")";
        // The camera presets and the display rows, in their popovers (the yarn's rows need a knit on show).
        p.injectKnit(false);
        m += missing(popoverLabels(p, W, H), {"3/4 view", "Front", "Side", "Top", "Back", "Upright", "Studio", "Flat", "Occlusion", "Wireframe",
                                              "Ply texture", "Compare with source"});
        p.ui.results.knit.reset();
        CHECK(p.setSetting("material.yarn_colour_source", "single"));  // the teddy opens with its colours from the object
        have = labels(p);
        m += missing(have, {"##Row direction", "##Yarn colour", "##Colour"});
        for (const char* gone : {"##Body stitch count", "##White bear", "##Brown bear"})
            if (have.count(gone)) m += std::string(" (former teddy control ") + gone + ")";
        CHECK(p.setSetting("strand.plies", "3"));  // separate ply tubes
        have = labels(p);
        m += missing(have, {"##Ply tubes", "##Ply sides", "##Ply ring spacing"});
        std::printf("    %s: missing%s\n", scene, m.empty() ? " none" : m.c_str());
        CHECK(m.empty());
    }
}

// A yarn of fibres (yarn_look.h): chosen on the Knit step's Construction tiles, it turns the yarn body on
// and disables the controls it makes meaningless (the ply texture); contact flattening and its Max
// flattening stay usable. Its rows are in the Fibres step (Yarn fibres), not in the Knit step
// (2026-11-13). Back to the single tube, the rows go and the ply texture returns.
TEST(ui_knit_fibres_disables) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Knit;
    CHECK(p.setSetting("strand.flatten", "true"));
    openEverything(p, W, H);
    auto state = [&](const std::string& label) {  // 1 disabled, 0 enabled, -1 not drawn
        const std::vector<ItemRec> items = settle(p, W, H);
        const ItemRec* r = findItem(items, label);
        return r ? int(r->disabled) : -1;
    };
    CHECK(state("##Ply texture") == 0 && state("##Fibres per yarn") == -1);
    const std::string profile = p.setting("strand.texture_profile");
    CHECK(click(p, W, H, "Fibre strands"));
    openEverything(p, W, H);
    // The ply texture keeps its value (unused): the knit stays current.
    CHECK(p.setting("yarn_body.enabled") == "true" && p.setting("yarn_body.render") == "fibres" && p.setting("strand.texture_profile") == profile);
    const int plyTexture = state("##Ply texture"), flatten = state("Contact flattening"), maxFlat = state("##Max flattening");
    std::printf("    fibre strands: ply texture %d, contact flattening %d, max flattening %d\n", plyTexture, flatten, maxFlat);
    CHECK(plyTexture == -1 && flatten == 0 && maxFlat == 0);  // the ply texture is the Single tube's sub-option: gone with it
    CHECK(state("##Fibres per yarn") == -1 && state("Set up in the Fibres step") == 0);  // the rows live in the Fibres step
    p.ui.stage = viewer::Stage::Fibres;
    openEverything(p, W, H);  // the fibre structure and appearance disclosures
    for (const char* rowId : {"##Look", "##Plies", "##Fibres per yarn", "##Interior fibres", "##Fibre diameter \xC2\xB7 um",
                              "##Shell depth", "##Fibre angle", "##Irregularity", "##Core", "##Fibre seed", "##Fibre IOR"})
        if (state(rowId) != 0) {
            std::printf("    missing or disabled: %s\n", rowId);
            CHECK(false);
        }
    CHECK(state("##Rendered as") == -1 && state("Render preview") == -1);  // both went with the tiles: the rendering is a tile, Preview is in the status bar
    p.ui.stage = viewer::Stage::Knit;
    CHECK(click(p, W, H, "Single tube"));
    CHECK(p.setting("yarn_body.enabled") == "false" && p.setting("strand.texture_profile") == profile);
    CHECK(state("##Ply texture") == 0);
    p.ui.stage = viewer::Stage::Fibres;
    CHECK(state("##Fibres per yarn") == -1 && state("##Look") == -1);
}

namespace {
int disabledItem(UiFixture& p, float W, float H, const std::string& label);  // below, with the Export tests

// The stage panel's scrolling child (interface.cpp drawSceneUi), whatever it shows: rows outside its view are
// clipped and register no item, so it is found by its name, not through an item.
ImGuiWindow* stagePanelWindow() {
    for (ImGuiWindow* w : GImGui->Windows)
        if (std::string(w->Name).find("##stagePanel") != std::string::npos) return w;
    return nullptr;
}

// The frame's text with every run of white space as one space: wrapped lines read as one.
std::string flatText(const std::string& text) {
    std::string out;
    for (char c : text) {
        const bool space = c == ' ' || c == '\n' || c == '\t' || c == '\r';
        if (space && !out.empty() && out.back() == ' ') continue;
        out += space ? ' ' : c;
    }
    return out;
}

size_t occurrences(const std::string& text, const std::string& what) {
    size_t n = 0;
    for (size_t at = text.find(what); at != std::string::npos; at = text.find(what, at + what.size())) ++n;
    return n;
}
} // namespace

// The review of the yarn made of fibres (2026-11-11, moved to the Fibres step 2026-11-13), on the blue swatch:
//   - Colourwork is offered disabled, and a settings file that chose it is warned about;
//   - in the Fibres step the three choices that make the yarn come first, the structure and the shading in
//     their own disclosures, and the footer's action is the fuzz's; Ply angle is disabled with one ply; the
//     Knit step has no fibre rows and no second Preview button;
//   - a refused setting is said under its control and below the rows, with its remedy, once in the Knit step
//     with a link, and blocks the Export step (its fix and heading), which names the remedy too;
//   - the Export switch is the fuzz's, the yarn fibres are listed apart, and the ready line tells them apart.
TEST(ui_fibre_yarn_feedback) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("wool-swatch-blue");
    p.injectKnit(true);
    p.ui.stage = viewer::Stage::Fibres;
    openEverything(p, W, H);
    const std::string missing = missingLabels(allLabels(p, W, H), {"##Look", "##Plies", "##Fibres per yarn", "##Fibre diameter \xC2\xB7 um",
                                                                   "Fibre structure", "##Shell depth", "##Core", "Fibre appearance", "##Fibre IOR",
                                                                   "Generate fuzz"});
    if (!missing.empty()) std::printf("    missing:%s\n", missing.c_str());
    CHECK(missing.empty());
    p.ui.stage = viewer::Stage::Knit;
    openEverything(p, W, H);
    const std::set<std::string> knitLabels = allLabels(p, W, H);
    CHECK(knitLabels.count("Render preview") == 0 && knitLabels.count("##Plies") == 0 && knitLabels.count("Set up in the Fibres step") == 1);
    // Colourwork: in the list, disabled.
    CHECK(click(p, W, H, "##Yarn colour"));
    int colourwork = -1;
    for (const ItemRec& r : settle(p, W, H))
        if (r.windowName.rfind("##Combo", 0) == 0 && r.label == "Colourwork (per stitch)") colourwork = int(r.disabled);
    pressKey(p, W, H, ImGuiKey_Escape);
    std::printf("    colourwork item: %d (1 disabled)\n", colourwork);
    CHECK(colourwork == 1);
    CHECK(p.setSetting("material.yarn_colour_source", "texture"));
    CHECK(flatText(renderedText(p, W, H)).find("Colourwork is not exported with a yarn made of fibres") != std::string::npos);
    CHECK(p.setSetting("material.yarn_colour_source", "single"));
    // Ply angle: nothing to twist with one ply.
    p.ui.stage = viewer::Stage::Fibres;
    CHECK(disabledItem(p, W, H, "##Ply angle") == 0);
    CHECK(p.setSetting("yarn_body.plies", "1"));
    CHECK(disabledItem(p, W, H, "##Ply angle") == 1);
    CHECK(p.setSetting("yarn_body.plies", "3"));
    // A refused setting: an opaque core at shell depth 100 %.
    const std::string remedy = "Reduce Shell depth or choose Core: None.";
    CHECK(flatText(renderedText(p, W, H)).find(remedy) == std::string::npos);
    CHECK(p.setSetting("yarn_body.shell_depth", "1"));
    const size_t fibresSays = occurrences(flatText(renderedText(p, W, H)), remedy);
    std::printf("    remedy shown %zu times in the Fibres step (under Core, below the rows)\n", fibresSays);
    CHECK(fibresSays >= 2);
    p.ui.stage = viewer::Stage::Knit;
    const std::string knitText = flatText(renderedText(p, W, H));
    CHECK(occurrences(knitText, remedy) == 1 && knitText.find("Correct it in the Fibres step") != std::string::npos);
    p.ui.stage = viewer::Stage::Export;
    const knit::ExportEligibility el = viewer::eligibility(p.ui);
    CHECK(el.state == knit::ExportEligibility::State::Blocked && el.fix == "Correct the yarn's fibres");
    const std::string exportText = flatText(renderedText(p, W, H));
    CHECK(exportText.find("Yarn fibres refused") != std::string::npos && exportText.find(remedy) != std::string::npos);
    CHECK(p.setSetting("yarn_body.shell_depth", "0.6"));
    // The Export step's contents: the fuzz's switch, the yarn fibres apart.
    const std::set<std::string> have = allLabels(p, W, H);
    CHECK(have.count("Include fuzz") == 1 && have.count("Include fibres") == 0);
    CHECK(flatText(renderedText(p, W, H)).find("| Yarn fibres (?) | Included through Construction |") != std::string::npos);
    CHECK(p.setSetting("fibres.enabled", "false"));
    const std::string ready = viewer::exportReadyLine(p.ui);
    std::printf("    ready line: %s\n", ready.c_str());
    CHECK(ready.find("yarn made of fibres") != std::string::npos && ready.find("no fuzz") != std::string::npos);
}

// The yarn of fibres rendered as maps (yarn_body.render): the Knit step's Fibre maps tile chooses it, the
// Fibres step's map rows appear with it (texel, cavity, the three appearance rows) and the core proxy's row
// goes, the cost line counts texels, the summaries say so; the Export step's contents and ready line follow;
// the Look preset row names the blue swatch's look and turns Custom on an edit.
TEST(ui_knit_fibres_rendered_as) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("wool-swatch-blue");
    p.injectKnit(true);
    p.setKnitYarnLength(12000);  // 12 m of yarn: the cost lines count
    p.ui.stage = viewer::Stage::Fibres;
    openEverything(p, W, H);
    auto drawn = [&](const std::string& label) { return findItem(settle(p, W, H), label) != nullptr; };
    CHECK(p.setting("yarn_body.render") == "fibres" && drawn("##Core") && !drawn("##Map texel \xC2\xB7 um"));
    std::string text = flatText(renderedText(p, W, H));
    CHECK(text.find("Yarn fibres: ~") != std::string::npos && text.find("| Look (?) | { Smooth } |") != std::string::npos);
    p.ui.stage = viewer::Stage::Knit;
    CHECK(click(p, W, H, "Fibre maps"));
    CHECK(p.setting("yarn_body.render") == "maps" && p.setting("yarn_body.enabled") == "true");
    openEverything(p, W, H);
    text = flatText(renderedText(p, W, H));
    CHECK(text.find("maps 40 um") != std::string::npos && text.find("Fibre maps") != std::string::npos);  // the summary line, the Yarn finish summary
    p.ui.stage = viewer::Stage::Fibres;
    openEverything(p, W, H);
    const std::string missing = missingLabels(allLabels(p, W, H), {"##Map texel \xC2\xB7 um", "##Cavity shading", "##Map roughness", "##Map specular",
                                                                   "##Highlight anisotropy", "##Look"});
    if (!missing.empty()) std::printf("    missing:%s\n", missing.c_str());
    CHECK(missing.empty() && !drawn("##Core"));
    text = flatText(renderedText(p, W, H));
    CHECK(text.find("Yarn maps: ~") != std::string::npos);
    // An edited look value: the preset row says Custom; the look put back by its id.
    CHECK(p.setSetting("yarn_body.map_cavity", "0.2"));
    CHECK(flatText(renderedText(p, W, H)).find("| Look (?) | { Custom } |") != std::string::npos);
    CHECK(p.setSetting("yarn_body.map_cavity", "0.6"));
    p.ui.stage = viewer::Stage::Export;
    openEverything(p, W, H);
    const std::string exportText = flatText(renderedText(p, W, H));
    CHECK(exportText.find("| Yarn fibres (?) | Baked into maps (Construction) |") != std::string::npos);
    CHECK(exportText.find("yarn-map texels") != std::string::npos);
    const std::string ready = viewer::exportReadyLine(p.ui);
    std::printf("    ready line: %s\n", ready.c_str());
    CHECK(ready.find("yarn of fibres as maps") != std::string::npos);
    CHECK(p.setSetting("fibres.enabled", "false"));  // no fuzz requested: nothing else blocks the export
    CHECK(viewer::eligibility(p.ui).state != knit::ExportEligibility::State::Blocked);
    // A texel too fine for a 500 mm yarn's atlas is refused with the texel's remedy (the fixture's knit holds no
    // curves, so the door is asked directly; the eligibility asks it with the knit's longest yarn).
    CHECK(knit::checkYarnMaps(p.ui.s, 500).accepted);
    CHECK(p.setSetting("yarn_body.map_texel_um", "10"));
    const knit::YarnBodyCheck tooFine = knit::checkYarnMaps(p.ui.s, 500);
    std::printf("    at 10 um: %s %s\n", tooFine.reason.c_str(), tooFine.remedy.c_str());
    CHECK(!tooFine.accepted && tooFine.field == knit::YarnBodyField::MapTexel && tooFine.remedy == "Choose a coarser Map texel.");
    CHECK(p.setSetting("yarn_body.map_texel_um", "40"));
    p.ui.stage = viewer::Stage::Knit;
    CHECK(click(p, W, H, "Fibre strands"));
    CHECK(p.setting("yarn_body.render") == "fibres");
}

// The Display popover's "Yarn maps" toggle appears once the knit on show carries baked cavity maps on its tubes
// (one mesh per yarn with an albedo texture), and not before; the Fibres step (Yarn fibres) offers "Bake the
// yarn maps" under maps until then. The Fibres step's viewport shows the model (the Display rows judge the
// shown scene; the reference patch has no maps).
TEST(ui_display_yarn_maps_toggle) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("wool-swatch-blue");
    p.injectKnit(true);
    p.ui.stage = viewer::Stage::Fibres;
    p.ui.view.fibres = viewer::FibreView::Model;
    CHECK(p.setSetting("yarn_body.render", "maps"));
    openEverything(p, W, H);
    CHECK(!viewer::displayRows(p.ui).yarnMaps);
    CHECK(allLabels(p, W, H).count("Bake the yarn maps") == 1);
    // The shown scene with its tubes under an atlas (what the strand mesh job builds).
    auto b = std::make_shared<knit::BuiltScene>(*p.ui.results.knit);
    b->settings = p.ui.s;
    auto tex = std::make_shared<knit::DisplayTexture>();
    tex->width = tex->height = 2;
    tex->rgba.assign(16, 255);
    b->albedo.push_back(tex);
    for (knit::DisplayMesh& m : b->meshes)
        if (m.kind == knit::DisplayMesh::Kind::Yarn) {
            m.name += "#y000";
            m.albedo = 0;
            m.albedoUv.assign(m.pos.size() / 3 * 2, 0.5f);
        }
    p.ui.results.knit = b;
    CHECK(viewer::displayRows(p.ui).yarnMaps && knit::carriesYarnMaps(*p.ui.results.knit));
    openEverything(p, W, H);
    CHECK(allLabels(p, W, H).count("Bake the yarn maps") == 0);
    CHECK(flatText(renderedText(p, W, H)).find("The yarn maps are baked") != std::string::npos);
    // A changed texel: the shown maps are stale, the bake is offered again.
    CHECK(p.setSetting("yarn_body.map_texel_um", "60"));
    CHECK(allLabels(p, W, H).count("Bake the yarn maps again") == 1);
    // The Display popover (a popup, not a disclosure) offers the toggle while the maps are on show.
    p.ui.view.mode = viewer::ViewMode::Knit;
    CHECK(clickSuffix(p, W, H, "###displayMenu") && popoverItem(settle(p, W, H), "Yarn maps"));
}

// The Construction tiles (2026-11-13, yarn_look.h): four, each with its detail line; the one the settings
// describe is selected; a click applies the construction (Fibre maps writes yarn_body.render) and, between a
// yarn of fibres and the single tube, the knit stays current; no job runs.
TEST(ui_construction_tiles) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(true);
    p.ui.stage = viewer::Stage::Knit;
    const std::string missing = missingLabels(allLabels(p, W, H), {"Single tube", "Ply tubes", "Fibre strands", "Fibre maps"});
    if (!missing.empty()) std::printf("    missing:%s\n", missing.c_str());
    CHECK(missing.empty());
    const std::string text = flatText(renderedText(p, W, H));
    for (const char* detail : {"One smooth tube", "Real ply geometry", "Individual fibres", "lower render cost"}) CHECK(text.find(detail) != std::string::npos);
    CHECK(knit::yarnConstructionOf(p.ui.s) == knit::YarnConstruction::SingleTube && viewer::knitState(p.ui).word == "Final");
    CHECK(click(p, W, H, "Fibre maps"));
    CHECK(p.setting("yarn_body.enabled") == "true" && p.setting("yarn_body.render") == "maps" && p.setting("strand.plies") == "1");
    CHECK(knit::yarnConstructionOf(p.ui.s) == knit::YarnConstruction::FibreMaps && viewer::knitState(p.ui).word == "Final");
    CHECK(flatText(renderedText(p, W, H)).find("Set up in the Fibres step") != std::string::npos);
    CHECK(click(p, W, H, "Fibre strands"));
    CHECK(p.setting("yarn_body.render") == "fibres" && knit::yarnConstructionOf(p.ui.s) == knit::YarnConstruction::FibreStrands);
    CHECK(click(p, W, H, "Ply tubes"));
    CHECK(p.setting("yarn_body.enabled") == "false" && p.setting("strand.plies") == "3" && findItem(settle(p, W, H), "##Ply tubes") != nullptr);
    CHECK(viewer::knitState(p.ui).word == "Outdated");  // ply tubes change the knit's surface
    CHECK(click(p, W, H, "Single tube"));
    CHECK(p.setting("strand.plies") == "1" && viewer::knitState(p.ui).word == "Final");
    CHECK(!p.ui.job.active && p.submittedJobs(0).empty());
}

// The Fibres step's Yarn fibres section (panel_yarn_fibres.cpp): shown for a yarn of fibres only (not for a
// tube, ply tubes or under Hair); its keys are Export-stage keys, so an edit leaves the knit Final and the
// fuzz Current; the strip reads "Yarn fibres" without fuzz and the fuzz's own words with it; the footer's
// action is the fuzz's by name (Generate fuzz / Update fuzz).
TEST(ui_fibres_yarn_section) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(true);
    p.ui.stage = viewer::Stage::Fibres;
    std::set<std::string> have = allLabels(p, W, H);
    CHECK(have.count("##Look") == 0 && have.count("Apply to model") == 1);
    CHECK(p.setSetting("yarn_body.enabled", "true"));
    have = allLabels(p, W, H);
    CHECK(have.count("##Look") == 1 && have.count("Fibre structure") == 1 && have.count("Generate fuzz") == 1 && have.count("Apply to model") == 0);
    CHECK(flatText(renderedText(p, W, H)).find("Edits here need no action") != std::string::npos);
    // Fuzz off: the strip says the yarn's fibres; on: the fuzz's own state.
    CHECK(p.setSetting("fibres.enabled", "false"));
    CHECK(viewer::fibreState(p.ui).word == "Yarn fibres");
    CHECK(p.setSetting("fibres.enabled", "true"));
    CHECK(viewer::fibreState(p.ui).word == "Not generated");
    p.injectFibres();
    CHECK(viewer::fibreState(p.ui).word == "Ready" && viewer::fibreOutputState(p.ui).word == "Current");
    CHECK(allLabels(p, W, H).count("Generate fuzz again") == 1);
    // A yarn fibre edit: the knit and the fuzz stay; a fuzz edit asks for Update fuzz.
    CHECK(p.setSetting("yarn_body.fibres", "500") && p.setSetting("yarn_body.shell_depth", "0.5"));
    CHECK(viewer::knitState(p.ui).word == "Final" && viewer::fibreOutputState(p.ui).word == "Current");
    CHECK(p.setSetting("fibres.density", "77"));
    CHECK(viewer::fibreOutputState(p.ui).word == "Outdated" && allLabels(p, W, H).count("Update fuzz") == 1);
    // Ply tubes: no section. Hair: no section (knitting skipped).
    CHECK(p.setSetting("yarn_body.enabled", "false") && p.setSetting("strand.plies", "3"));
    CHECK(allLabels(p, W, H).count("##Look") == 0);
    CHECK(p.setSetting("strand.plies", "1") && p.setSetting("yarn_body.enabled", "true"));
    viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);
    CHECK(viewer::furShown(p.ui) && allLabels(p, W, H).count("##Look") == 0);
}

// A link to the yarn's fibres reveals them (2026-11-13): the stage panel is one scrolling child shared by the
// steps, so a step left scrolled far down lands the next step there too. In a short window, with Fibre
// structure collapsed and Shell depth refused (the Core row is named): the Export fix clicked with the Export
// step scrolled to its end, and the Knit step's link clicked at the top of its view, both open the Fibres step
// with the node open and the Core row in view; "Set up in the Fibres step" brings the Look row into view.
TEST(ui_yarn_fibres_reveal) {
    const float W = 1600, H = 360;
    UiContext ctx(W, H);
    UiFixture p("wool-swatch-blue");
    p.injectKnit(true);
    CHECK(p.setSetting("yarn_body.shell_depth", "1"));  // refused: a full shell with an opaque core (the Core row is named)
    CHECK(!knit::checkYarnBody(p.ui.s).accepted && knit::checkYarnBody(p.ui.s).field == knit::YarnBodyField::Core);
    auto collapseStructure = [&] {
        p.ui.stage = viewer::Stage::Fibres;
        for (const ItemRec& r : settle(p, W, H))
            if (r.label == "Fibre structure" && r.window) r.window->StateStorage.SetInt(r.id, 0);
        settle(p, W, H);
    };
    // The step's panel scrolled so that the item `label` sits at the top of its view ("": to the end).
    auto scrollPanel = [&](viewer::Stage stage, const std::string& label) {
        p.ui.stage = stage;
        const std::vector<ItemRec> items = settle(p, W, H);
        ImGuiWindow* panel = stagePanelWindow();
        CHECK(panel != nullptr);
        if (!panel) return;
        if (label.empty()) {
            panel->ScrollTarget.y = 1e6f;
        } else {
            const ItemRec* r = findItem(items, label);
            CHECK(r != nullptr);
            if (!r) return;
            panel->ScrollTarget.y = panel->Scroll.y + (r->bb.Min.y - panel->InnerRect.Min.y);
        }
        panel->ScrollTargetCenterRatio.y = 0;
        settle(p, W, H);
        std::printf("    %s step scrolled to %.0f px\n", stage == viewer::Stage::Knit ? "Knit" : "Export", panel->Scroll.y);
        CHECK(panel->Scroll.y > 60);  // far enough that the Fibres step's heading would be above the view
    };
    collapseStructure();
    scrollPanel(viewer::Stage::Export, "");
    CHECK(click(p, W, H, "Correct the yarn's fibres"));
    CHECK(p.ui.stage == viewer::Stage::Fibres);
    CHECK(insidePanel(p, W, H, "##Core"));
    CHECK(p.ui.panels.revealYarnFibres == viewer::Reveal::None);  // used up
    // The Knit step's link does the same.
    collapseStructure();
    scrollPanel(viewer::Stage::Knit, "Correct it in the Fibres step");
    CHECK(click(p, W, H, "Correct it in the Fibres step"));
    CHECK(p.ui.stage == viewer::Stage::Fibres && insidePanel(p, W, H, "##Core"));
    // "Set up in the Fibres step": the section's rows come into view.
    CHECK(p.setSetting("yarn_body.shell_depth", "0.6"));
    collapseStructure();
    scrollPanel(viewer::Stage::Knit, "Set up in the Fibres step");
    CHECK(click(p, W, H, "Set up in the Fibres step"));
    CHECK(p.ui.stage == viewer::Stage::Fibres && insidePanel(p, W, H, "##Look"));
}

// The Knit step's summary line counts the generated fuzz as the Fibres step's cost line does (one owner,
// yarnFibresCost): the same GiB in both, above the yarn's own; without fuzz generated it says so.
TEST(ui_yarn_fibres_cost_counts_fuzz) {
    UiFixture p("wool-swatch-blue");
    p.injectKnit(true);
    p.setKnitYarnLength(12000);
    CHECK(p.setSetting("fibres.enabled", "true"));
    std::string summary = viewer::yarnFibresSummary(p.ui);
    std::printf("    %s\n", summary.c_str());
    CHECK(summary.find("fuzz not generated yet") != std::string::npos);
    p.injectFibres({{5000000, 5000000, 4000000}});
    const viewer::YarnFibresCost c = viewer::yarnFibresCost(p.ui);
    summary = viewer::yarnFibresSummary(p.ui);
    const std::string line = viewer::yarnFibresCostLine(c);
    std::printf("    %s\n    %s\n", summary.c_str(), line.c_str());
    CHECK(c.knitted && c.fuzzSegments == 4000000 && c.yarnSegments > 0);
    char gib[32];
    std::snprintf(gib, sizeof gib, "~%.1f GiB", c.gib);
    CHECK(summary.find("+ 4.0 M fuzz") != std::string::npos && summary.find(gib) != std::string::npos && line.find(gib) != std::string::npos);
    CHECK(c.gib > knit::rendererMemoryGiB(std::llround(c.yarnSegments), 0));
    CHECK(p.setSetting("fibres.enabled", "false"));
    CHECK(viewer::yarnFibresSummary(p.ui).find("fuzz") == std::string::npos);
}

// The yarn map bake started from the Fibres step (startYarnMapBake, a strand mesh job of the Knit kind): the
// section shows its progress and Cancel, the Knit footer the same job by name, and the Fibres footer still
// names the fuzz.
TEST(ui_yarn_map_bake_progress) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("wool-swatch-blue");
    p.injectKnit(true);
    CHECK(p.setSetting("yarn_body.render", "maps"));
    p.ui.stage = viewer::Stage::Fibres;
    CHECK(allLabels(p, W, H).count("Bake the yarn maps") == 1);
    CHECK(click(p, W, H, "Bake the yarn maps"));
    CHECK(p.ui.job.active && p.ui.job.meshOnly && p.ui.job.bakeMaps);
    const std::string text = flatText(renderedText(p, W, H));
    CHECK(text.find("Baking yarn maps ...") != std::string::npos && text.find("Fuzz") != std::string::npos);
    CHECK(allLabels(p, W, H).count("Cancel the bake") == 1);
    CHECK(viewer::knitState(p.ui).word == "Baking" && viewer::progressView(p.ui.job, false).title == "Baking yarn maps ...");
    p.ui.stage = viewer::Stage::Knit;
    CHECK(flatText(renderedText(p, W, H)).find("Baking yarn maps ...") != std::string::npos);
    for (int i = 0; i < 2000 && p.ui.job.active; ++i) {  // the harness fails the job at once
        viewer::pollSceneUi(p.ui);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    CHECK(!p.ui.job.active);
}

// The baked maps' material (yarn_body.map_material, Fibres step, Fibre appearance): Fabric by default; the three
// highlight rows are enabled under both materials (since the renderer's R1 the fabric carries the highlight too);
// choosing Shiny in the combo writes the key; the summary says which, with the highlight's values.
namespace {
bool pickItem(UiFixture& p, float W, float H, const std::string& label);  // below, with the Export tests
}

TEST(ui_yarn_map_material) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("wool-swatch-blue");
    p.injectKnit(true);
    CHECK(p.setSetting("yarn_body.render", "maps"));
    p.ui.stage = viewer::Stage::Fibres;
    openEverything(p, W, H);
    const auto disabled = [&](const std::string& label) {
        const ItemRec* r = findItem(settle(p, W, H), label);
        CHECK(r != nullptr);
        return r && r->disabled;
    };
    const char* const highlight[] = {"##Map roughness", "##Map specular", "##Highlight anisotropy"};
    CHECK(p.setting("yarn_body.map_material") == "fabric" && !disabled("##Map material"));
    for (const char* l : highlight) CHECK(!disabled(l));
    std::string text = flatText(renderedText(p, W, H));
    CHECK(text.find("| Map material (?) | { Fabric } |") != std::string::npos && text.find("maps: fabric, roughness") != std::string::npos);
    CHECK(click(p, W, H, "##Map material") && pickItem(p, W, H, "Shiny"));
    CHECK(p.setting("yarn_body.map_material") == "shiny");
    for (const char* l : highlight) CHECK(!disabled(l));
    text = flatText(renderedText(p, W, H));
    CHECK(text.find("| Map material (?) | { Shiny } |") != std::string::npos && text.find("maps: shiny, roughness") != std::string::npos);
    CHECK(click(p, W, H, "##Map material") && pickItem(p, W, H, "Fabric"));
    CHECK(p.setting("yarn_body.map_material") == "fabric");
}

// One workflow (user, 2026-10-28): every knitted example shows the same controls in the Knit, Fibres and
// Export steps as the teddy, with knitting on, one yarn colour and no fibres in the export (the examples' own
// colour sources, fibre defaults and listed settings set aside: a row such as Inflation shows for a setting's
// value, the same in every scene). Imports differ only in the Import step (units, axis, charts), never here (owner, 2026-11-11: one
// workflow for every preset and every import, no scene of its own kind).
TEST(ui_same_workflow_all_scenes) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    auto controls = [&](const std::string& scene, int stage) {
        UiFixture p(scene);
        viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
        CHECK(p.setSetting("material.yarn_colour_source", "single"));
        CHECK(p.setSetting("fibres.enabled", "false"));
        const knit::SceneSettings reference = knit::viewerDefaults(*knit::findScene("teddy"));
        for (const auto& [key, value] : knit::findScene(scene)->defaultSettings) knit::copySetting(p.ui.s, reference, key);
        p.injectKnit(true);
        p.ui.stage = static_cast<viewer::Stage>(stage);
        openEverything(p, W, H);
        std::set<std::string> l;
        frame(p, W, H);  // the footer keeps last frame's height: one frame to settle after it grows
        for (const ItemRec& r : frame(p, W, H))
            if (r.windowName.rfind("Workflow", 0) == 0 && !r.label.empty()) l.insert(r.label);
        return l;
    };
    std::string m;
    for (int stage = 1; stage <= 3; ++stage) {
        const std::set<std::string> r = controls("teddy", stage);
        for (const knit::SceneEntry& e : knit::sceneRegistry()) {
            if (e.id == "teddy") continue;
            const std::set<std::string> have = controls(e.id, stage);
            std::string only, missing;
            for (const std::string& l : have)
                if (!r.count(l)) only += " '" + l + "'";
            for (const std::string& l : r)
                if (!have.count(l)) missing += " '" + l + "'";
            std::printf("    stage %d %s: extra%s, missing%s\n", stage, e.id.c_str(), only.empty() ? " none" : only.c_str(),
                        missing.empty() ? " none" : missing.c_str());
            if (!(only.empty() && missing.empty())) m += " " + e.id;        }
    }
    CHECK(m.empty());
}

// The teddy knits as spiral pieces behind the same controls (user, 2026-10-28): opening it gives its Dense
// (1.86 at 5.5 x 3.4, Final fallback from 1.80), stockinette, one colour per bear and its fibre defaults,
// whatever the scene before had; a round trip through the finished wool swatch gives the same again.
TEST(ui_teddy_pieces_scene_switch) {
    knit_test::requireScene(knit::findProjectFile(knit::findScene("teddy")->piecesManifest));  // the Knit step's controls follow the loaded pieces
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    auto state = [&] {
        return p.setting("strand.thickness_preset") + " " + p.setting("strand.diameter") + " " + p.setting("strand.wale") + " " +
               p.setting("strand.course") + " fallback [" + p.setting("strand.diameter_fallback") + "] " + p.setting("strand.pattern") +
               " " + p.setting("material.yarn_colour_source");
    };
    auto num = [&](const char* key) { return std::atof(p.setting(key).c_str()); };
    CHECK(p.setSetting("strand.pattern", "rib1x1"));
    CHECK(p.setSetting("material.yarn_colour_source", "single"));
    viewer::openExample(p.ui, "teddy");
    std::printf("    teddy (Draft): %s\n", state().c_str());
    CHECK(knit::isPreviewProfile(p.ui.s) && num("strand.diameter") == 1.86 && num("strand.wale") == 5.5 && num("strand.course") == 3.4);
    CHECK(p.setting("strand.thickness_preset") == "thick-dense");
    CHECK(p.setting("strand.pattern") == "stockinette" && p.setting("material.yarn_colour_source") == "object");
    CHECK(p.gaugePresetMatches());
    // The teddy's own fibre defaults (owner, 2026-11-05): density 135 and Natural islands.
    CHECK(num("fibres.density") == 135 && num("fibres.islands.rate") == 0.06 && num("fibres.islands.matted") == 1);
    CHECK(p.setting("knit.coverage") == "component");  // spiral pieces have no charts to continue
    p.injectKnit(true);  // the full profile applied first
    std::printf("    teddy Final: %s\n", state().c_str());
    CHECK(p.finalProfile() && p.setting("strand.diameter_fallback").find("1.73") != std::string::npos);
    viewer::openExample(p.ui, "wool-swatch");
    viewer::openExample(p.ui, "teddy");
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);  // an opened example asks first
    std::printf("    wool-swatch -> teddy (Final): %s\n", state().c_str());
    CHECK(num("strand.diameter") == 1.86 && p.setting("material.yarn_colour_source") == "object");
    CHECK(p.finalProfile() && p.setting("strand.diameter_fallback").find("1.73") != std::string::npos);
    CHECK(p.setting("knit.coverage") == "component");
    // The Knit step draws the same controls; the pieces' fixed ones are disabled, not removed.
    p.ui.stage = viewer::Stage::Knit;
    openEverything(p, W, H);  // Stitch layout
    bool rowDirection = false;
    for (const ItemRec& r : frame(p, W, H)) rowDirection |= r.label == "##Row direction";
    CHECK(rowDirection);
}

// Bindings through real input: dragging Stitch width changes the settings' wale and makes the
// gauge Custom (no preset matches).
TEST(ui_knit_drag_dimension) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Knit;
    const std::string before = p.setting("strand.wale");
    CHECK(p.gaugePresetMatches());
    CHECK(drag(p, W, H, "##Stitch width", 80));
    std::printf("    stitch width %s -> %s mm, preset %s\n", before.c_str(), p.setting("strand.wale").c_str(), p.setting("strand.thickness_preset").c_str());
    CHECK(std::atof(p.setting("strand.wale").c_str()) > std::atof(before.c_str()));
    CHECK(p.setting("strand.thickness_preset") == "custom" && !p.gaugePresetMatches());
}

// The one Gauge list, through real input: Snug / Regular / Loose derive the stitch from the yarn (a
// round yarn fits: the flattening checkbox is inert), a diameter drag moves the stitch with it, the
// stitch is a readout (no field to drag), Custom brings the fields back with the values kept, a preset
// returns to its fixed stitch.
TEST(ui_knit_gauge_from_yarn) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    auto num = [](UiFixture& p, const char* k) { return std::atof(p.setting(k).c_str()); };
    auto fits = [&](UiFixture& p) {  // yarnFit's need >= 1 (stockinette's limits)
        const double need = (std::min(num(p, "strand.wale") / knit::kFitW, num(p, "strand.course") / knit::kFitH) - num(p, "strand.clearance")) /
                            num(p, "strand.diameter");
        return need >= 1 - 1e-9;
    };
    // A click on the open combo's item whose label starts with `prefix` (the popup window "##Combo_...").
    auto pick = [&](UiFixture& p, const std::string& prefix) {
        const std::vector<ItemRec> items = settle(p, W, H);
        for (const ItemRec& r : items)
            if (r.windowName.rfind("##Combo", 0) == 0 && r.label.rfind(prefix, 0) == 0) return click(p, W, H, r.label);
        return false;
    };
    // 1 / 0: the item labelled `label` is drawn disabled / enabled; -1: not drawn.
    auto disabled = [&](UiFixture& p, const std::string& label) {
        const std::vector<ItemRec> items = settle(p, W, H);
        const ItemRec* r = findItem(items, label);
        return r ? int(r->disabled) : -1;
    };
    // The readout of a stitch dimension in the panel's text ("11.60 mm" on its row).
    auto readout = [&](UiFixture& p, double mm) {
        char b[32];
        std::snprintf(b, sizeof b, "%.2f mm", mm);
        return renderedText(p, W, H).find(b) != std::string::npos;
    };
    const std::string flatten = "Flatten yarn that is too thick for the stitch";
    for (const char* scene : {"teddy"}) {
        UiFixture p(scene);
        p.ui.stage = viewer::Stage::Knit;
        openEverything(p, W, H);  // Gauge & yarn geometry: the flattening checkbox
        // Dense: a fixed stitch too small for a round 1.86 mm yarn, so the flattening applies.
        CHECK(p.setting("strand.stitch_size") == "custom" && p.gaugePresetMatches());
        CHECK(disabled(p, "##Stitch width") == 0 && disabled(p, flatten) == 0);
        CHECK(click(p, W, H, "##Gauge") && pick(p, "Regular"));
        CHECK(p.setting("strand.stitch_size") == "yarn" && p.setting("strand.tightness") == "standard" &&
              p.setting("strand.thickness_preset") == "custom");
        // The stitch follows the yarn: readouts instead of fields, and a round yarn fits (nothing to flatten).
        CHECK(disabled(p, "##Stitch width") == -1 && disabled(p, "##Row height") == -1 && disabled(p, "##Yarn diameter") == 0 &&
              disabled(p, flatten) == 1);
        CHECK(readout(p, num(p, "strand.wale")) && readout(p, num(p, "strand.course")));
        const double w0 = num(p, "strand.wale"), h0 = num(p, "strand.course"), D0 = num(p, "strand.diameter");
        CHECK(fits(p));
        CHECK(drag(p, W, H, "##Yarn diameter", 80));  // thicker yarn: a larger stitch, still Regular
        const double w1 = num(p, "strand.wale"), D1 = num(p, "strand.diameter");
        CHECK(D1 > D0 && w1 > w0 && num(p, "strand.course") > h0 && fits(p) && p.setting("strand.stitch_size") == "yarn");
        CHECK(readout(p, w1));
        CHECK(!drag(p, W, H, "##Stitch width", 80));  // no field in yarn mode: nothing to drag
        CHECK(num(p, "strand.wale") == w1);
        CHECK(click(p, W, H, "##Gauge") && pick(p, "Loose"));
        CHECK(p.setting("strand.tightness") == "loose" && num(p, "strand.wale") > w1 && fits(p));
        CHECK(click(p, W, H, "##Gauge") && pick(p, "Snug"));
        CHECK(p.setting("strand.tightness") == "tight" && num(p, "strand.wale") < w1 && fits(p));
        std::printf("    %s: D %.2f -> %.2f, w %.2f -> %.2f (regular), snug %.2f x %.2f\n", scene, D0, D1, w0, w1, num(p, "strand.wale"),
                    num(p, "strand.course"));
        // Custom keeps the values and makes the stitch editable.
        const double wSnug = num(p, "strand.wale");
        CHECK(click(p, W, H, "##Gauge") && pick(p, "Custom"));
        CHECK(p.setting("strand.stitch_size") == "custom" && num(p, "strand.wale") == wSnug && disabled(p, "##Stitch width") == 0);
        CHECK(click(p, W, H, "##Gauge") && pick(p, "Dense"));
        CHECK(p.setting("strand.stitch_size") == "custom" && p.gaugePresetMatches() && disabled(p, flatten) == 0);
    }
}

// The footer's two actions: each sets its real solve profile (the export rule's predicate) and submits
// one knit job (not a mesh rebuild); from a Final profile Knit draft goes back to Draft.
TEST(ui_knit_draft_final_actions) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Knit;
    knit::applyFullProfile(p.ui.s, *viewer::entryOf(p.ui));
    CHECK(click(p, W, H, "Knit draft"));
    std::vector<std::string> jobs = p.submittedJobs(1);
    CHECK(knit::isPreviewProfile(p.ui.s) && !p.finalProfile());
    CHECK(p.ui.job.active && p.ui.job.kind == knit::BuiltScene::Kind::Knit && !p.ui.job.meshOnly && jobs.size() == 1 && jobs[0] == "knit");
    UiFixture q("teddy");
    q.ui.stage = viewer::Stage::Knit;
    CHECK(knit::isPreviewProfile(q.ui.s));
    CHECK(click(q, W, H, "Knit final"));
    jobs = q.submittedJobs(1);
    CHECK(q.finalProfile() && !knit::isPreviewProfile(q.ui.s));
    CHECK(q.ui.job.active && q.ui.job.kind == knit::BuiltScene::Kind::Knit && !q.ui.job.meshOnly && jobs.size() == 1 && jobs[0] == "knit");
}

// Yarn tessellation stays mesh-only: choosing another level with a current knit submits the mesh
// rebuild (StrandMeshRequest, no knitting) and the knit stays current (not Outdated). The fibre warning
// follows real invalidation: none while the generated fibres match the tubes (Mesh detail open, level
// unchanged), drawn once another level makes them stale (their state Outdated).
TEST(ui_knit_mesh_detail_mesh_only) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    const std::string stale = "The generated fibres no longer match these yarn tubes";
    auto warned = [&](UiFixture& f) { return renderedText(f, W, H).find(stale) != std::string::npos; };
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Knit;
    p.injectKnit(true);
    p.injectFibres();
    openEverything(p, W, H);  // Mesh detail
    const std::string level = p.meshDetail(), knitWord = viewer::knitState(p.ui).word;
    CHECK(!level.empty() && knitWord == "Final");
    CHECK(!viewer::fibresNeedNewTubes(p.ui) && viewer::fibreOutputState(p.ui).word == "Current" && !warned(p));
    CHECK(click(p, W, H, "##Yarn tessellation"));  // opens the list
    std::string other;
    for (const ItemRec& r : frame(p, W, H))
        if (r.windowName.rfind("##Combo", 0) == 0 && !r.label.empty() && r.label != level && r.label != "Custom") other = r.label;
    CHECK(!other.empty());
    CHECK(click(p, W, H, other));
    std::printf("    mesh detail %s -> %s; job submitted %d, mesh only %d; Knit \"%s\"\n", level.c_str(), p.meshDetail().c_str(), p.ui.job.active,
                (p.ui.job.active && p.ui.job.meshOnly), viewer::knitState(p.ui).word.c_str());
    CHECK(p.meshDetail() == other);
    CHECK(p.ui.job.active && (p.ui.job.active && p.ui.job.meshOnly));
    CHECK(viewer::knitCurrent(p.ui));  // the knitting key is unchanged: no knitting needed
    std::printf("    fibres after the change: need new tubes %d, state \"%s\"\n", viewer::fibresNeedNewTubes(p.ui),
                viewer::fibreOutputState(p.ui).word.c_str());
    CHECK(viewer::fibresNeedNewTubes(p.ui) && viewer::fibreOutputState(p.ui).word == "Outdated" && warned(p));
}

// ---- the Fibres step ----

// Every former Fibres-step control is still drawn (presence only; bindings: the tests below): both
// population tabs, every node open, the fuzz's own colour, both shading models, the toolbar's
// preview target and the footer; the hair profile on its own.
TEST(ui_fibre_controls_inventory) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    p.injectKnit(true);
    CHECK(p.setSetting("material.fibre_inherits_yarn", "false"));
    std::set<std::string> have = allLabels(p, W, H);
    std::string m = missingLabels(have, {"##Wool preset", "##Stray preset", "##Density preset", "##Density", "##Stray fibres", "##Shortest", "##Typical", "##Longest",
                                         "##Length spread", "##Thinnest", "##Thickest", "##Bow", "##Max lift angle", "##Lie on the stitches", "Match yarn colour",
                                         "##Colour", "##Shading", "##Roughness along", "##Roughness around", "##Segment budget",
                                         "##More at curved yarn", "##Segments per fibre", "##Generation",
                                         "Reference patch", "Model", "Source", "Knit", "Fibres", "Apply to model"});
    // Standard wool's fuzz is hair-shaded; the other shading model has its own row.
    CHECK(p.setSetting("material.fibre_directed.model", "diffuse-sheen"));
    m += missingLabels(allLabels(p, W, H), {"##Sheen"});
    // The export's fibre switch is the Export step's only (Contents, Include fibres).
    if (have.count("Include fibres in export") || have.count("Include fibres")) m += " (the export's fibre switch in the Fibres step)";
    CHECK(findPrefix(frame(p, W, H), "Yarn fuzz ") && findPrefix(frame(p, W, H), "Stray fibres "));
    p.ui.panels.fibreTabRequest = 1;
    settle(p, W, H);  // the selection takes effect on the next frame
    CHECK(p.setting("material.fibre_random.model") == "fibre" && p.setting("material.fibre_random_inherits_yarn") == "true");
    have = allLabels(p, W, H);
    if (have.count("##Colour")) m += " (a stray colour row while the strays match the yarn)";
    CHECK(p.setSetting("material.fibre_random_inherits_yarn", "false"));
    have = allLabels(p, W, H);
    m += missingLabels(have, {"##Shortest", "##Typical", "##Longest", "##Length spread", "##Thinnest", "##Thickest", "##Bend", "##Max angle",
                              "##Direction spread", "Match yarn colour", "##Colour", "##Shading", "##Roughness along", "##Roughness around"});
    if (have.count("##Bow") || have.count("##Lie on the stitches")) m += " (the fuzz's controls on the stray tab)";
    UiFixture hair("teddy");
    viewer::chooseWorkflow(hair.ui, knit::Workflow::Fur);
    hair.ui.stage = viewer::Stage::Fibres;
    have = allLabels(hair, W, H);
    m += missingLabels(have, {"Reference sphere", "Model", "Grow fur", "Entire surface", "##Density", "##Length", "##Thickness", "##Taper",
                              "##Bend", "##Curl", "##Gravity", "##Length Var.", "##Thickness Var.", "##Direction Var.", "##Bend Var.",
                              "##Curl Var.", "Fine fibres over the loops", "##Density##fibres", "##Lift", "##Roughness along", "Seed",
                              "Segment budget"});
    if (have.count("##Wool preset") || have.count("##Stray preset") || have.count("##Density preset") || have.count("##Profile"))
        m += " (wool-only controls in the Fur step)";
    std::printf("    missing%s\n", m.empty() ? " none" : m.c_str());
    CHECK(m.empty());
}

// The population tabs only choose what is edited: switching changes no setting, view, preview
// target or job; an edit in one tab changes that population only, and both keep their values.
TEST(ui_fibre_tabs_preserve) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    const std::string before = knit::writeSettings(p.ui.s);
    const int view = static_cast<int>(p.ui.view.mode);
    const bool model = (p.ui.view.fibres == viewer::FibreView::Model);
    CHECK(click(p, W, H, "Stray fibres ", true));
    CHECK(click(p, W, H, "Yarn fuzz ", true));
    CHECK(knit::writeSettings(p.ui.s) == before && static_cast<int>(p.ui.view.mode) == view && (p.ui.view.fibres == viewer::FibreView::Model) == model && !p.ui.job.active);
    const char* dirKey = "fibres.directed.length_median";
    const char* rndKey = "fibres.random.length_median";
    const std::string d0 = p.setting(dirKey), r0 = p.setting(rndKey);
    CHECK(drag(p, W, H, "##Typical", 60));  // the fuzz tab
    const std::string d1 = p.setting(dirKey);
    CHECK(d1 != d0 && p.setting(rndKey) == r0);
    CHECK(click(p, W, H, "Stray fibres ", true));
    CHECK(drag(p, W, H, "##Typical", 60));
    const std::string r1 = p.setting(rndKey);
    std::printf("    typical length: fuzz %s -> %s, stray %s -> %s; preset %s\n", d0.c_str(), d1.c_str(), r0.c_str(), r1.c_str(),
                p.setting("fibres.preset").c_str());
    CHECK(p.setting(dirKey) == d1 && r1 != r0);
    CHECK(p.setting("fibres.preset") == "custom");
    CHECK(click(p, W, H, "Yarn fuzz ", true));
    CHECK(p.setting(dirKey) == d1 && p.setting(rndKey) == r1 && !p.ui.job.active);
}

// The stray share is a percentage on screen and a fraction in the settings (exact at 0, 20 and
// 100 %); the tabs show both shares; the candidate line is density x yarn length.
TEST(ui_fibre_share_percent) {
    CHECK(viewer::sharePercent(0) == 0 && viewer::shareFraction(0) == 0);
    CHECK(viewer::sharePercent(0.2) == 20 && viewer::shareFraction(20) == 0.2);
    CHECK(viewer::sharePercent(1) == 100 && viewer::shareFraction(100) == 1);
    CHECK(viewer::shareFraction(viewer::sharePercent(0.2)) == 0.2);
    const std::string estimate = viewer::candidateEstimate(125, 129500);
    std::printf("    %s\n", estimate.c_str());
    CHECK(estimate == "Approx. 16.2M candidates \xC2\xB7 129.5 m yarn");
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    CHECK(p.setSetting("fibres.random_fraction", "0.2"));
    std::vector<ItemRec> items = settle(p, W, H);
    CHECK(findPrefix(items, "Yarn fuzz 80%###") && findPrefix(items, "Stray fibres 20%###"));
    CHECK(drag(p, W, H, "##Stray fibres", 2000));
    const double high = std::atof(p.setting("fibres.random_fraction").c_str());
    items = settle(p, W, H);
    CHECK(high == 1.0 && findPrefix(items, "Stray fibres 100%###") && findPrefix(items, "Yarn fuzz 0%###"));
    CHECK(drag(p, W, H, "##Stray fibres", -4000));
    const double low = std::atof(p.setting("fibres.random_fraction").c_str());
    std::printf("    stray share after dragging right %g, left %g\n", high, low);
    CHECK(low == 0.0);
}

// The wool preset combo applies a preset's values (no job); editing a value makes it Custom.
TEST(ui_fibre_preset_custom) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    CHECK(click(p, W, H, "##Wool preset"));
    CHECK(click(p, W, H, "Hard wool"));
    const std::string hard = p.setting("fibres.preset");
    CHECK(click(p, W, H, "##Wool preset"));
    CHECK(click(p, W, H, "Standard wool"));
    const std::string standard = p.setting("fibres.preset"), fraction = p.setting("fibres.random_fraction");
    CHECK(drag(p, W, H, "##Bow", 40));
    std::printf("    presets: %s, %s (stray share %s); after an edit: %s; jobs %zu\n", hard.c_str(), standard.c_str(), fraction.c_str(),
                p.setting("fibres.preset").c_str(), p.submittedJobs(0).size());
    CHECK(hard == "hard" && standard == "standard" && std::atof(fraction.c_str()) == 0.03);
    CHECK(p.setting("fibres.preset") == "custom");
    CHECK(!p.ui.job.active && p.submittedJobs(0).empty());
}

// The density preset combo sets fibres.density only (no job); the teddy opens with its own density (135,
// scene_registry.cpp), and dragging Coverage's Density leaves every preset's value.
TEST(ui_fibre_density_preset) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    const std::string initial = p.setting("fibres.density"), preset = p.setting("fibres.preset");
    CHECK(click(p, W, H, "##Density preset"));
    CHECK(click(p, W, H, "Light"));
    const std::string light = p.setting("fibres.density");
    CHECK(click(p, W, H, "##Density preset"));
    CHECK(click(p, W, H, "Dense"));
    const std::string dense = p.setting("fibres.density");
    CHECK(drag(p, W, H, "##Density", 40));
    const double dragged = std::atof(p.setting("fibres.density").c_str());
    std::printf("    density: initial %s, Light %s, Dense %s, dragged %g; wool preset %s -> %s; jobs %zu\n", initial.c_str(), light.c_str(),
                dense.c_str(), dragged, preset.c_str(), p.setting("fibres.preset").c_str(), p.submittedJobs(0).size());
    CHECK(initial == "135" && light == "30" && dense == "125");
    CHECK(dragged != 30 && dragged != 60 && dragged != 125);
    CHECK(p.setting("fibres.preset") == preset);  // the wool preset does not own the density
    CHECK(!p.ui.job.active && p.submittedJobs(0).empty());
}

// The Stitch character tiles are the hand-knit presets: a click applies one (and selects its tile), a
// Fine-tune edit makes the character Custom (no tile selected, the heading says so); Off is a preset
// like the others; the per-stitch noise and the columns are separate values. No job.
TEST(ui_hand_preset) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Knit;
    // 1 / 0: the tile labelled `label` is drawn selected / not selected; -1: not drawn.
    auto selected = [&](const std::string& label) {
        const std::vector<ItemRec> items = settle(p, W, H);
        const ItemRec* r = findItem(items, label);
        if (!r) return -1;
        const knit::FilePreset* m = knit::matchingFilePreset(knit::PresetGroup::HandKnit, p.ui.s);
        return m && m->label == label ? 1 : 0;
    };
    std::vector<std::string> tiles;
    for (const knit::FilePreset& f : knit::presets(knit::PresetGroup::HandKnit)) tiles.push_back(f.label);
    std::printf("    tiles:");
    for (const std::string& t : tiles) std::printf(" %s (%d)", t.c_str(), selected(t));
    std::printf("\n");
    CHECK(tiles.size() == 4 && selected("Hand-made") == 1 && selected("Rustic") == 0 && selected("Off") == 0 && selected("Natural") == 0);
    CHECK(click(p, W, H, "Rustic"));
    CHECK(p.setting("knit.hand.unevenness") == "32" && p.setting("knit.hand.column_shift") == "1");
    CHECK(knit::matchingFilePreset(knit::PresetGroup::HandKnit, p.ui.s)->id == "rustic" && selected("Rustic") == 1);
    openEverything(p, W, H);  // Fine-tune character
    CHECK(drag(p, W, H, "##Column shift", 20));
    CHECK(p.setting("knit.hand.column_shift") != "1" && p.setting("knit.hand.unevenness") == "32");
    CHECK(knit::matchingFilePreset(knit::PresetGroup::HandKnit, p.ui.s) == nullptr);  // Custom
    for (const std::string& t : tiles) CHECK(selected(t) == 0);
    CHECK(renderedText(p, W, H).find("Your own values (Fine-tune character).") != std::string::npos);
    CHECK(click(p, W, H, "Off"));
    CHECK(!p.ui.s.hand.inUse() && selected("Off") == 1);
    CHECK(!p.ui.job.active && p.submittedJobs(0).empty());
}

// The stray preset combo sets the strays' deviation only (no job): the wool preset keeps its name,
// and editing Direction spread on the stray tab makes the stray preset Custom, not the wool preset.
TEST(ui_stray_preset) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    const std::string wool = p.setting("fibres.preset");
    CHECK(click(p, W, H, "##Stray preset"));
    CHECK(click(p, W, H, "Wild"));
    const knit::FilePreset* wild = knit::matchingFilePreset(knit::PresetGroup::Strays, p.ui.s);
    CHECK(wild && wild->id == "wild" && p.setting("fibres.random.spread_deg") == "180" && p.setting("fibres.random.max_angle_deg") == "60");
    settle(p, W, H);
    CHECK(p.setting("fibres.preset") == wool && wool != "custom");
    p.ui.panels.fibreTabRequest = 1;
    settle(p, W, H);
    CHECK(drag(p, W, H, "##Direction spread", -40));
    settle(p, W, H);
    std::printf("    wool preset %s; Wild spread 180 -> %s after a drag; jobs %zu\n", wool.c_str(), p.setting("fibres.random.spread_deg").c_str(),
                p.submittedJobs(0).size());
    CHECK(p.setting("fibres.random.spread_deg") != "180" && knit::matchingFilePreset(knit::PresetGroup::Strays, p.ui.s) == nullptr);
    CHECK(p.setting("fibres.preset") == wool);
    CHECK(!p.ui.job.active && p.submittedJobs(0).empty());
}

TEST(ui_island_preset) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    CHECK(click(p, W, H, "##Island preset"));
    CHECK(click(p, W, H, "Well-loved"));
    CHECK(p.ui.s.islands.inUse() && p.setting("fibres.islands.strength") == "100");
    CHECK(click(p, W, H, "##Island preset"));
    CHECK(click(p, W, H, "Off"));
    CHECK(!p.ui.s.islands.inUse());
    CHECK(!p.ui.job.active && p.submittedJobs(0).empty());
}

// The footer keeps the generated fibres and the export switch apart: "Off" in the stage bar is
// the export switch; the output is Not generated / Current / Outdated; an appearance edit (an
// Export-stage key) leaves the fibres current; Generate submits one fibres job on the shown knit
// and does not touch the export switch.
TEST(ui_fibre_footer_states) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    p.injectKnit(true);
    auto drawn = [&](const char* label) { return findItem(settle(p, W, H), label) != nullptr; };
    CHECK(viewer::fibreOutputState(p.ui).word == "Not generated" && viewer::fibreState(p.ui).word == "Off" && drawn("Apply to model"));
    p.ui.s.fibres = true;
    CHECK(viewer::fibreOutputState(p.ui).word == "Not generated" && viewer::fibreState(p.ui).word == "Not generated" && (viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked));
    p.injectFibres();
    CHECK(viewer::fibreOutputState(p.ui).word == "Current" && viewer::fibreState(p.ui).word == "Ready" && drawn("Apply again") && !(viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked));
    p.ui.s.fibres = false;
    std::printf("    generated, not included: output \"%s\", stage bar \"%s\" (%s)\n", viewer::fibreOutputState(p.ui).word.c_str(), viewer::fibreState(p.ui).word.c_str(),
                viewer::fibreState(p.ui).detail.c_str());
    CHECK(viewer::fibreOutputState(p.ui).word == "Current" && viewer::fibreState(p.ui).word == "Off" && viewer::fibreState(p.ui).detail == "fibres generated; not in the export");
    CHECK(p.setSetting("material.fibre_random.sheen", "0.7"));
    CHECK(viewer::fibreOutputState(p.ui).word == "Current");
    CHECK(p.setSetting("fibres.density", "90"));
    CHECK(viewer::fibreOutputState(p.ui).word == "Outdated" && viewer::fibreState(p.ui).word == "Outdated" && drawn("Apply to model"));
    const std::string included = p.setting("fibres.enabled");
    CHECK(!(p.ui.view.fibres == viewer::FibreView::Model));  // the step opens on the reference patch
    CHECK(click(p, W, H, "Apply to model"));
    const std::vector<std::string> jobs = p.submittedJobs(1);
    std::printf("    jobs: %s; fibres.enabled %s -> %s; model view %d, view %d\n", jobs.empty() ? "none" : jobs[0].c_str(), included.c_str(),
                p.setting("fibres.enabled").c_str(), (p.ui.view.fibres == viewer::FibreView::Model), static_cast<int>(p.ui.view.mode));
    CHECK(jobs.size() == 1 && jobs[0] == "fibres knit=1");
    CHECK(p.setting("fibres.enabled") == included && !(p.ui.view.fibres == viewer::FibreView::Model) && p.ui.view.mode == viewer::ViewMode::Fibres);  // the patch stays shown
}

namespace {
bool pickItem(UiFixture& p, float W, float H, const std::string& label);  // below, with the Export tests
}

// On the reference patch every Fibres-step control rebuilds the preview by itself: each edit
// changes the preview's key, and the next poll (each viewer frame) requests a build for it - no
// button. A control that changes no key, or whose key is not requested, is listed.
TEST(ui_fibre_preview_follows_edits) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    p.injectKnit(true);
    CHECK(p.setSetting("material.fibre_directed.model", "diffuse-sheen"));  // the matte fuzz: its Sheen row drawn
    openEverything(p, W, H);  // every node open: all controls drawn
    CHECK(!(p.ui.view.fibres == viewer::FibreView::Model));
    CHECK(p.fibrePreviewRequested() == knit::fibrePreviewKey(p.ui.s));
    std::string missed;
    int followed = 0;
    auto verdict = [&](const std::string& what, uint64_t before) {
        const uint64_t now = knit::fibrePreviewKey(p.ui.s);
        if (now == before) missed += " " + what + " (no change)";
        else if (p.fibrePreviewRequested() != now) missed += " " + what + " (not requested)";
        else ++followed;
    };
    auto slider = [&](const std::string& label) {
        const uint64_t before = knit::fibrePreviewKey(p.ui.s);
        if (!drag(p, W, H, label, 60)) {
            missed += " " + label + " (not drawn)";
            return;
        }
        if (knit::fibrePreviewKey(p.ui.s) == before) drag(p, W, H, label, -120);  // at its maximum
        verdict(label, before);
    };
    auto combo = [&](const std::string& label) {
        const uint64_t before = knit::fibrePreviewKey(p.ui.s);
        if (!click(p, W, H, label)) {
            missed += " " + label + " (not drawn)";
            return;
        }
        std::vector<std::string> entries;
        for (const ItemRec& r : settle(p, W, H))
            if (r.windowName.rfind("##Combo", 0) == 0 && !r.label.empty()) entries.push_back(r.label);
        for (size_t i = 0; i < entries.size() && knit::fibrePreviewKey(p.ui.s) == before; ++i) {
            if (i > 0) click(p, W, H, label);
            pickItem(p, W, H, entries[i]);
        }
        verdict(label, before);
    };
    for (const char* l : {"##Density", "##Stray fibres", "##Shortest", "##Typical", "##Longest", "##Length spread", "##Thinnest",
                          "##Thickest", "##Bow", "##Max lift angle", "##Sheen", "##More at curved yarn", "##Segments per fibre"})
        slider(l);
    for (const char* l : {"##Wool preset", "##Stray preset", "##Density preset", "##Shading", "##Generation"}) combo(l);
    {
        CHECK(p.setSetting("material.fibre_inherits_yarn", "false"));  // as the inventory: the checkbox drawn
        const uint64_t before = knit::fibrePreviewKey(p.ui.s);
        if (click(p, W, H, "Match yarn colour")) verdict("Match yarn colour", before);
        else missed += " Match yarn colour (not drawn)";
    }
    p.ui.panels.fibreTabRequest = 1;
    settle(p, W, H);
    openEverything(p, W, H);
    for (const char* l : {"##Shortest", "##Typical", "##Longest", "##Thinnest", "##Thickest", "##Bend", "##Max angle"}) slider(std::string(l));
    std::printf("    %d controls followed; missed%s\n", followed, missed.empty() ? " none" : missed.c_str());
    CHECK(missed.empty());
}

// The preview target (toolbar) switches the view only: no job.
TEST(ui_fibre_preview_target) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    p.ui.view.fibres = viewer::FibreView::Reference;
    CHECK(click(p, W, H, "Model"));
    CHECK((p.ui.view.fibres == viewer::FibreView::Model) && !p.ui.job.active);
    CHECK(click(p, W, H, "Reference patch"));
    CHECK(!(p.ui.view.fibres == viewer::FibreView::Model) && !p.ui.job.active && p.submittedJobs(0).empty());
    CHECK(click(p, W, H, "Model") && (p.ui.view.fibres == viewer::FibreView::Model));
    p.ui.stage = viewer::Stage::Knit;
    CHECK(!findItem(settle(p, W, H), "Reference patch"));  // the Fibres step's own control
    p.ui.stage = viewer::Stage::Fibres;
    settle(p, W, H);
    CHECK(!(p.ui.view.fibres == viewer::FibreView::Model));  // the step always opens on the reference patch
}

// The fibre job's state has one home (the Fibres footer); on another stage the status bar has Cancel.
TEST(ui_fibre_progress_single_home) {
    const float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(true);
    p.injectRunning(1, "generating fibres", 3, 10, false);
    p.ui.stage = viewer::Stage::Fibres;
    std::vector<ItemRec> items = settle(p, W, H);
    int panel = 0, status = 0;
    for (const ItemRec& r : items) {
        if (r.label == "Cancel fibres" && r.windowName.rfind("Workflow", 0) == 0) ++panel;
        if (r.label == "Cancel" && r.windowName == "##status") ++status;
    }
    std::printf("    Fibres stage: Cancel fibres in the panel %d, Cancel in the status bar %d\n", panel, status);
    CHECK(panel == 1 && status == 0);
    p.ui.stage = viewer::Stage::Knit;
    items = settle(p, W, H);
    panel = countItems(items, "Cancel fibres");
    status = 0;
    for (const ItemRec& r : items) status += r.label == "Cancel" && r.windowName == "##status";
    std::printf("    Knit stage: Cancel fibres %d, Cancel in the status bar %d\n", panel, status);
    CHECK(panel == 0 && status == 1);
}

// The Fibres footer's action stays on screen on a short window at 1.5 x interface scale: running,
// after a long failure with the fibres included in the export, and in the hair profile.
TEST(ui_fibre_footer_reachable) {
    const float W = 1280, H = 720;
    ScaledUi ctx(W, H, 1.5f);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    p.injectKnit(true);
    p.injectRunning(1, "generating fibres", 3, 10, false);
    CHECK(insidePanel(p, W, H, "Cancel fibres"));
    UiFixture q("teddy");
    q.ui.stage = viewer::Stage::Fibres;
    q.injectKnit(true);
    q.ui.s.fibres = true;
    setAttempt(q.ui.results.fibreAttempt, longMessage(150), true);
    CHECK(insidePanel(q, W, H, "Apply to model"));
    UiFixture r("teddy");
    viewer::chooseWorkflow(r.ui, knit::Workflow::Fur);
    r.ui.stage = viewer::Stage::Fibres;
    setAttempt(r.ui.results.fibreAttempt, longMessage(150), false);
    CHECK(insidePanel(r, W, H, "Grow fur"));
}

// ---- the panel's layout is at rest ----

namespace {

// The items of the panel ("Workflow" and its children) in one frame.
std::vector<ItemRec> panelItems(const std::vector<ItemRec>& items) {
    std::vector<ItemRec> out;
    for (const ItemRec& r : items)
        if (r.windowName.rfind("Workflow", 0) == 0) out.push_back(r);
    return out;
}

// "" when both frames draw the same items at the same places, else the first difference.
std::string layoutDifference(const std::vector<ItemRec>& a, const std::vector<ItemRec>& b) {
    char t[256];
    if (a.size() != b.size()) {
        std::snprintf(t, sizeof t, "%zu items, then %zu", a.size(), b.size());
        return t;
    }
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i].id != b[i].id || a[i].bb.Min.x != b[i].bb.Min.x || a[i].bb.Min.y != b[i].bb.Min.y || a[i].bb.Max.x != b[i].bb.Max.x ||
            a[i].bb.Max.y != b[i].bb.Max.y) {
            std::snprintf(t, sizeof t, "%s at (%.1f, %.1f)-(%.1f, %.1f), then (%.1f, %.1f)-(%.1f, %.1f)", describe(a[i]).c_str(), a[i].bb.Min.x,
                          a[i].bb.Min.y, a[i].bb.Max.x, a[i].bb.Max.y, b[i].bb.Min.x, b[i].bb.Min.y, b[i].bb.Max.x, b[i].bb.Max.y);
            return t;
        }
    return "";
}

// The panel widths in [w0, w1], 1 px apart, at which the layout keeps moving: after it settles,
// six more frames must draw every item where the first did. Every width is printed with its
// first difference.
int restlessWidths(UiFixture& p, float W, float H, float w0, float w1) {
    frame(p, W, H);  // the panel exists
    ImGuiWindow* panel = ImGui::FindWindowByName("Workflow");
    if (!panel) throw std::runtime_error("restlessWidths: no Workflow window");
    int restless = 0;
    for (float pw = w0; pw <= w1; pw += 1) {
        panel->SizeFull.x = pw;  // as a drag of its edge (the size constraints still apply)
        for (int i = 0; i < 6; ++i) frame(p, W, H);
        const std::vector<ItemRec> ref = panelItems(frame(p, W, H));
        for (int i = 0; i < 6; ++i) {
            const std::string d = layoutDifference(ref, panelItems(frame(p, W, H)));
            if (d.empty()) continue;
            if (restless < 8) std::printf("    width %.0f, frame +%d: %s\n", panel->Size.x, i + 1, d.c_str());
            ++restless;
            break;
        }
    }
    return restless;
}

// A running job whose progress line changes every frame (as a knit reports its parts): the
// frames, of `n`, in which the panel ("Workflow") itself shows a scrollbar - each such frame
// narrows every control under it for that frame (the stage panel scrolls; the panel never does).
int panelScrollbarFrames(UiFixture& p, float W, float H, int kind, int n) {
    int shown = 0;
    for (int i = 0; i < n; ++i) {
        const std::string stage = i % 3 == 0 ? "knitting part " + std::to_string(i) : "knitting " + longMessage(4 + i % 7);
        p.injectRunning(kind, stage, i % 18, 18, true);
        frame(p, W, H);
        const ImGuiWindow* panel = ImGui::FindWindowByName("Workflow");
        if (panel && panel->ScrollbarY) {
            if (shown < 4) std::printf("    frame %d: the panel shows a scrollbar (progress \"%s\")\n", i, stage.c_str());
            ++shown;
        }
    }
    return shown;
}

} // namespace

// User, 2026-10-29: the panel's controls flickered between two positions. Each step with a
// footer, with and without a long status message, at every panel width from the narrowest up
// (so the footers' wrapped lines cross every wrap edge): the layout comes to rest.
TEST(ui_layout_at_rest) {
    const float W = 1600, H = 900;
    ScaledUi ctx(W, H, 1.5f);
    const float w0 = viewer::uiPx(440), w1 = viewer::uiPx(600);  // from the narrowest panel (interface.cpp)
    int restless = 0;
    {
        UiFixture p("teddy");
        p.ui.stage = viewer::Stage::Knit;
        p.injectKnit(false);
        std::printf("  Knit, Draft result:\n");
        restless += restlessWidths(p, W, H, w0, w1);
        setAttempt(p.ui.results.knitAttempt, longMessage(40), true);
        std::printf("  Knit, failed attempt:\n");
        restless += restlessWidths(p, W, H, w0, w1);
    }
    {
        UiFixture p("teddy");
        p.ui.stage = viewer::Stage::Fibres;
        p.injectKnit(true);
        std::printf("  Fibres:\n");
        restless += restlessWidths(p, W, H, w0, w1);
        setAttempt(p.ui.results.fibreAttempt, longMessage(40), false);
        std::printf("  Fibres, attempt message:\n");
        restless += restlessWidths(p, W, H, w0, w1);
        p.ui.stage = viewer::Stage::Export;
        std::printf("  Export:\n");
        restless += restlessWidths(p, W, H, w0, w1);
    }
    int flashes = 0;
    {
        UiFixture p("teddy");
        p.ui.stage = viewer::Stage::Knit;
        p.injectKnit(false);
        std::printf("  Knit, running:\n");
        flashes += panelScrollbarFrames(p, W, H, 0, 120);
    }
    {
        UiFixture p("teddy");
        p.ui.stage = viewer::Stage::Fibres;
        p.injectKnit(true);
        std::printf("  Fibres, running:\n");
        flashes += panelScrollbarFrames(p, W, H, 1, 120);
    }
    std::printf("  %d restless widths, %d frames with a panel scrollbar\n", restless, flashes);
    CHECK(restless == 0);
    CHECK(flashes == 0);
}

// ---- the Export step ----

namespace {

// The labels of the panel's items ("Workflow" and its children: not the toolbar or the status
// bar) with every disclosure open.
std::set<std::string> panelLabels(UiFixture& p, float W, float H) {
    openEverything(p, W, H);
    std::set<std::string> l;
    for (const ItemRec& r : frame(p, W, H))
        if (r.windowName.rfind("Workflow", 0) == 0) l.insert(r.label);
    return l;
}

// 1 / 0: the item labelled `label` is drawn disabled / enabled; -1: not drawn.
int disabledItem(UiFixture& p, float W, float H, const std::string& label) {
    const std::vector<ItemRec> items = settle(p, W, H);
    const ItemRec* r = findItem(items, label);
    return r ? int(r->disabled) : -1;
}

// A click on the open combo's item labelled `label` (the popup window "##Combo_...").
bool pickItem(UiFixture& p, float W, float H, const std::string& label) {
    for (const ItemRec& r : settle(p, W, H))
        if (r.windowName.rfind("##Combo", 0) == 0 && r.label == label) return click(p, W, H, label);
    return false;
}

} // namespace

// Every former Export-step control is still drawn (presence; bindings: the tests below), every
// disclosure open: the package and OBJ, a Draft knit (its fix), the Fur workflow (nothing knitted:
// no fibre, filling or detail choice; the body colour stays and still sets the colour source).
// The OBJ's own description makes no package-only claim.
TEST(ui_export_controls_inventory) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(false);
    p.ui.stage = viewer::Stage::Export;
    std::set<std::string> have = panelLabels(p, W, H);
    std::string m = missingLabels(have, {"##Format", "What the package holds", "Include fibres", "Filling", "Colour", "##Unknitted parts",
                                         "##Level", "Choose folder ...", "Existing files", "Knit at Final quality"});
    if (have.count("##File name")) m += " (a file name for the package, whose names are fixed)";
    p.ui.exportTo.format = "obj";
    have = panelLabels(p, W, H);
    m += missingLabels(have, {"##Format", "What the files hold", "Include fibres", "Filling", "##Unknitted parts", "##Level",
                              "Choose folder ...", "##File name", "Existing files"});
    if (have.count("What the package holds")) m += " (the package's disclosure with OBJ)";
    std::string err;
    CHECK(knit::manifestScene("own/manifest.txt", "piece.obj", 24, &err) != nullptr);  // not a preset: its fur export keeps the body
    UiFixture hair("pieces");
    viewer::chooseWorkflow(hair.ui, knit::Workflow::Fur);
    hair.ui.stage = viewer::Stage::Export;
    have = panelLabels(hair, W, H);
    m += missingLabels(have, {"##Format", "What the package holds", "##Body colour", "Choose folder ...", "Existing files"});
    if (have.count("Include fibres") || have.count("Filling") || have.count("##Level") || have.count("##Unknitted parts"))
        m += " (knit-only controls with knitting skipped)";
    std::printf("    missing%s\n", m.empty() ? " none" : m.c_str());
    CHECK(m.empty());
    CHECK(click(hair, W, H, "##Body colour") && pickItem(hair, W, H, "Average colour"));
    CHECK(hair.setting("material.plain_colour_source") == "average");
    const knit::ExportFormat* obj = knit::findExportFormat("obj");
    CHECK(obj && obj->fibres && !obj->summary.empty());  // fibres and hair as OBJ lines
    for (const std::string& t : {obj->summary, obj->contents})
        CHECK(t.find("pptstrands") == std::string::npos && t.find("scene.json") == std::string::npos && t.find("shading maps (") == std::string::npos);
}

// A preset with knitting skipped exports like an import: its mesh with the hair. Only the missing
// hair blocks it, and the body colour is offered.
TEST(ui_preset_with_knitting_skipped_exports_its_mesh) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);
    p.ui.stage = viewer::Stage::Export;
    const knit::ExportEligibility el = viewer::eligibility(p.ui);
    std::printf("    blocked %d, fix \"%s\"\n", (viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked), el.fix.c_str());
    CHECK((viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked) && el.fix == "Grow fur" && el.reasons.size() == 1);
    std::set<std::string> have = panelLabels(p, W, H);
    CHECK(have.count("##Body colour") && !have.count("Knit at Final quality") && !have.count("Filling"));
}

// The Contents and Tessellation controls write the settings they always wrote: Include fibres is
// the export's only fibre switch (the Fibres footer only reports it); the OBJ format omits requested fibres without dropping the request; the
// filling colour only edits with the filling on; the colour source; a detail level sets the strand
// mesh and budget, rebuilds the mesh only and leaves Draft / Final alone.
TEST(ui_export_bindings) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(true);
    p.ui.stage = viewer::Stage::Export;
    CHECK(p.setting("fibres.enabled") == "false" && viewer::fibreExport(p.ui) == viewer::FibreExport::Excluded && viewer::fibreExportWord(p.ui) == "Excluded");
    CHECK(click(p, W, H, "Include fibres"));
    CHECK(p.setting("fibres.enabled") == "true" && viewer::fibreExport(p.ui) == viewer::FibreExport::Generate && (viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked));
    // The Fibres footer has no switch of its own: it says what the export does with the fibres.
    p.ui.stage = viewer::Stage::Fibres;
    CHECK(findItem(settle(p, W, H), "Include fibres") == nullptr && findItem(settle(p, W, H), "Include fibres in export") == nullptr);
    CHECK(std::string(viewer::fibreFooterExportNote(p.ui)) == std::string("Export: included, generate the fibres first"));
    p.ui.stage = viewer::Stage::Export;
    CHECK(click(p, W, H, "##Format") && pickItem(p, W, H, "Wavefront OBJ"));
    std::printf("    OBJ with fibres requested: \"%s\", ready line \"%s\"\n", viewer::fibreExportWord(p.ui).c_str(), viewer::exportReadyLine(p.ui).c_str());
    // OBJ carries the fibres too (as lines): the same request, the same need to generate them.
    CHECK(p.ui.exportTo.format == "obj" && p.setting("fibres.enabled") == "true" && viewer::fibreExport(p.ui) == viewer::FibreExport::Generate);
    CHECK(disabledItem(p, W, H, "Include fibres") == 0 && (viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked));
    CHECK(click(p, W, H, "##Format") && pickItem(p, W, H, "Photon Path Tracer package"));
    CHECK(p.ui.exportTo.format == "ppt" && viewer::fibreExport(p.ui) == viewer::FibreExport::Generate && (viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked));
    p.injectFibres();
    CHECK(viewer::fibreExport(p.ui) == viewer::FibreExport::Included && viewer::fibreExportWord(p.ui) == "Included" && !(viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked));
    CHECK(viewer::exportReadyLine(p.ui) == "Photon Path Tracer package \xC2\xB7 fibres included");
    CHECK(click(p, W, H, "##Format") && pickItem(p, W, H, "Wavefront OBJ"));
    std::printf("    OBJ with fibres: \"%s\", ready line \"%s\"\n", viewer::fibreExportWord(p.ui).c_str(), viewer::exportReadyLine(p.ui).c_str());
    CHECK(viewer::fibreExport(p.ui) == viewer::FibreExport::Included && viewer::fibreExportWord(p.ui) == "Included" && !(viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked));
    CHECK(viewer::exportReadyLine(p.ui) == "Wavefront OBJ \xC2\xB7 fibres included");
    CHECK(click(p, W, H, "##Format") && pickItem(p, W, H, "Photon Path Tracer package"));
    CHECK(click(p, W, H, "Include fibres"));
    CHECK(p.setting("fibres.enabled") == "false" && viewer::fibreExport(p.ui) == viewer::FibreExport::Excluded);
    CHECK(std::string(viewer::fibreFooterExportNote(p.ui)) == std::string("Export: not included (Export step, Contents)"));
    // The filling: on by default, its colour only edits with it on.
    CHECK(p.setting("filling.enabled") == "true" && disabledItem(p, W, H, "Colour") == 0);
    CHECK(click(p, W, H, "Filling"));
    CHECK(p.setting("filling.enabled") == "false" && disabledItem(p, W, H, "Colour") == 1);
    CHECK(click(p, W, H, "Filling") && p.setting("filling.enabled") == "true");
    const bool draft = knit::isPreviewProfile(p.ui.s), final = p.finalProfile();
    CHECK(p.exportDetailChosen() == "Reduced" && !p.ui.job.active);  // the default
    CHECK(click(p, W, H, "##Level"));
    int levels = 0;
    for (const ItemRec& r : settle(p, W, H))
        if (r.windowName.rfind("##Combo", 0) == 0 && (r.label == "Full" || r.label == "Reduced" || r.label == "Light")) ++levels;
    CHECK(levels == 3);
    CHECK(click(p, W, H, "Light"));
    std::printf("    level %s (%s); job %d, mesh only %d; knit current %d\n", p.exportDetailChosen().c_str(), p.exportDetailSettings().c_str(),
                p.ui.job.active, (p.ui.job.active && p.ui.job.meshOnly), viewer::knitCurrent(p.ui));
    CHECK(p.exportDetailChosen() == "Light" && p.exportDetailSettings() == "6/0.8/5000000");
    CHECK(knit::isPreviewProfile(p.ui.s) == draft && p.finalProfile() == final);
    CHECK(p.ui.job.active && (p.ui.job.active && p.ui.job.meshOnly) && viewer::knitCurrent(p.ui));
}

// Generated wool fibres tick Export > Contents > Include fibres by themselves; unticking it is
// kept until the next fibre result.
TEST(ui_export_fibres_auto_include) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(true);
    p.ui.stage = viewer::Stage::Export;
    CHECK(p.setting("fibres.enabled") == "false" && viewer::fibreExportWord(p.ui) == "Excluded");
    viewer::takeFibreResult(p.ui, p.fakeFibres());
    std::printf("    after the fibre job: fibres.enabled %s, \"%s\"\n", p.setting("fibres.enabled").c_str(), viewer::fibreExportWord(p.ui).c_str());
    CHECK(p.setting("fibres.enabled") == "true" && viewer::fibreExportWord(p.ui) == "Included");
    CHECK(click(p, W, H, "Include fibres") && p.setting("fibres.enabled") == "false");
    settle(p, W, H);
    CHECK(p.setting("fibres.enabled") == "false" && viewer::fibreExportWord(p.ui) == "Excluded: fibres are generated, tick to include");
    viewer::takeFibreResult(p.ui, p.fakeFibres());  // a new result
    CHECK(p.setting("fibres.enabled") == "true");
}

// The Blender format: first in the list, a new panel starts on it where Blender is installed, and
// the panel describes its contents. The fibres go with it, as curves.
TEST(ui_export_blender_format) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    CHECK(knit::exportFormats().front().id == "blend");
    UiFixture p("teddy");
    CHECK(viewer::SceneUi{}.exportTo.format == knit::defaultExportFormat());
    p.injectKnit(true);
    p.ui.exportTo.format = "blend";
    p.ui.stage = viewer::Stage::Export;
    const std::set<std::string> have = panelLabels(p, W, H);
    CHECK(missingLabels(have, {"##Format", "What the files hold", "Include fibres", "##Level", "Export Blender file"}).empty());
    CHECK(p.exportSummary() == knit::findExportFormat("blend")->summary);
    p.ui.s.fibres = true;
    CHECK((viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked) && viewer::fibreExportWord(p.ui) == "Requested: generate them (Fibres)");
    p.injectFibres();
    CHECK(!(viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked) && viewer::fibreExportWord(p.ui) == "Included");
    CHECK(viewer::exportReadyLine(p.ui) == "Blender file \xC2\xB7 fibres included");
}

// The Output folder line is the folder the export job is given (every format); an invalid name
// shows the destination check's reason and disables both export actions, but never a fix that
// knits or generates fibres.
TEST(ui_export_target_matches) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    for (const char* f : {"ppt", "obj", "blend"}) {
        UiFixture p("teddy");
        p.injectKnit(true);
        p.ui.exportTo.format = f;
        p.ui.stage = viewer::Stage::Export;
        const std::string shown = viewer::exportOutputText(p.ui);
        CHECK(shown == p.ui.exportTo.dir);
        const std::string format = f;
        CHECK(click(p, W, H, format == "ppt" ? "Export package" : format == "obj" ? "Export OBJ" : "Export Blender file"));
        const std::vector<std::string> jobs = p.submittedJobs(1);
        std::printf("    %s: shown %s; job %s\n", f, shown.c_str(), jobs.empty() ? "none" : jobs[0].c_str());
        CHECK(jobs.size() == 1 && jobs[0] == "export dir=" + shown + " format=" + f);
    }
    for (const char* bad : {"a/b", "CON", ""}) {
        UiFixture q("teddy");
        q.ui.stage = viewer::Stage::Export;
        q.ui.exportTo.format = "obj";
        q.setExportName(bad);
        std::printf("    name \"%s\": %s\n", bad, viewer::exportOutputText(q.ui).c_str());
        CHECK(viewer::exportOutputText(q.ui).rfind("Cannot export there: ", 0) == 0);
        q.injectKnit(false);  // Draft
        CHECK(disabledItem(q, W, H, "Knit at Final quality") == 0);
        q.injectKnit(true);
        CHECK(disabledItem(q, W, H, "Export OBJ") == 1);
        q.injectKnit(true, 1);
        CHECK(disabledItem(q, W, H, "Export OBJ") == 1);
        q.ui.s.fibres = true;
        CHECK(disabledItem(q, W, H, "Generate fibres") == 0);
        q.ui.exportTo.format = "ppt";  // the package's file names are fixed: the typed name does not stop it
        q.injectKnit(true);
        CHECK(viewer::exportOutputText(q.ui) == q.ui.exportTo.dir);
    }
}

// The folder the user chooses is the folder the job writes: no subfolder, for every format, with
// the file name passed on its own. It stays while the knit is switched off and goes with the scene.
TEST(ui_export_writes_into_the_chosen_folder) {
    namespace fs = std::filesystem;
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    const std::string chosen = (fs::temp_directory_path() / "knit_ui_test_chosen_folder").make_preferred().string();
    std::error_code ec;
    fs::create_directories(chosen, ec);
    for (const char* f : {"ppt", "obj", "blend"}) {
        UiFixture p("teddy");
        p.injectKnit(true);
        p.ui.exportTo.format = f;
        p.ui.stage = viewer::Stage::Export;
        p.ui.platform.dialogs = true;
        CHECK(click(p, W, H, "Choose folder ...") && p.answerFolderDialog(chosen));
        settle(p, W, H);
        std::printf("    %s: target %s\n", f, p.ui.exportTo.dir.c_str());
        CHECK(p.ui.exportTo.dir == chosen && viewer::exportOutputText(p.ui) == chosen);
        const std::string format = f;
        CHECK(click(p, W, H, format == "ppt" ? "Export package" : format == "obj" ? "Export OBJ" : "Export Blender file"));
        const std::vector<std::string> jobs = p.submittedJobs(1);
        CHECK(jobs.size() == 1 && jobs[0] == "export dir=" + chosen + " format=" + f);
    }
    UiFixture s("teddy");
    s.ui.stage = viewer::Stage::Export;
    s.ui.platform.dialogs = true;
    CHECK(click(s, W, H, "Choose folder ...") && s.answerFolderDialog(chosen));
    settle(s, W, H);
    viewer::chooseWorkflow(s.ui, knit::Workflow::Fur);
    CHECK(s.ui.exportTo.dir == chosen && s.ui.exportTo.name == "teddy_fur");
    viewer::chooseWorkflow(s.ui, knit::Workflow::Knit);
    viewer::openExample(s.ui, "teddy");
    std::printf("    after another scene: %s\n", s.ui.exportTo.dir.c_str());
    CHECK(s.ui.exportTo.dir != chosen && fs::path(s.ui.exportTo.dir).filename() == s.ui.exportTo.name);
    fs::remove_all(chosen, ec);
}

// The footer's actions do what they did: the Draft fix knits at Final quality in the Knit step
// (one knit job, no export); a missing fibre result is generated in the Fibres step; a diagnostic
// export starts at once, as any other. A completed export names the folder that job published,
// after the name and format have changed (the Output folder line follows the new settings).
TEST(ui_export_footer_states) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(false);
    p.ui.stage = viewer::Stage::Export;
    CHECK(knit::isPreviewProfile(p.ui.s) && (viewer::eligibility(p.ui).state == knit::ExportEligibility::State::Blocked));
    CHECK(click(p, W, H, "Knit at Final quality"));
    std::vector<std::string> jobs = p.submittedJobs(1);
    CHECK(jobs.size() == 1 && jobs[0] == "knit" && p.ui.stage == viewer::Stage::Knit && p.finalProfile());
    UiFixture f("teddy");
    f.injectKnit(true);
    f.ui.s.fibres = true;
    f.ui.stage = viewer::Stage::Export;
    CHECK(click(f, W, H, "Generate fibres"));
    jobs = f.submittedJobs(1);
    CHECK(jobs.size() == 1 && jobs[0] == "fibres knit=1" && f.ui.stage == viewer::Stage::Fibres);
    UiFixture q("teddy");
    q.injectKnit(true, 1);
    q.ui.stage = viewer::Stage::Export;
    CHECK(click(q, W, H, "Export package"));
    jobs = q.submittedJobs(1);
    CHECK(jobs.size() == 1 && jobs[0] == "export dir=" + q.ui.exportTo.dir + " format=ppt");
    UiFixture d("teddy");
    d.injectKnit(true);
    d.ui.stage = viewer::Stage::Export;
    const std::string published = d.ui.exportTo.dir;
    d.injectExportDone(published);
    CHECK(viewer::exportDoneText(d.ui) == "Exported to " + published && viewer::exportState(d.ui).word == "Exported");
    const std::string another = (std::filesystem::temp_directory_path() / "knit_ui_test_another_folder").make_preferred().string();
    d.ui.platform.dialogs = true;
    CHECK(click(d, W, H, "Choose folder ...") && d.answerFolderDialog(another));
    d.ui.exportTo.format = "obj";
    settle(d, W, H);
    std::printf("    done: \"%s\"; output now %s\n", viewer::exportDoneText(d.ui).c_str(), viewer::exportOutputText(d.ui).c_str());
    CHECK(viewer::exportDoneText(d.ui) == "Exported to " + published);
    CHECK(viewer::exportOutputText(d.ui) == d.ui.exportTo.dir && viewer::exportOutputText(d.ui) != published);
}

// The Export footer's action stays on screen on a short window at 1.5 x interface scale: blocked
// with several reasons, ready and diagnostic after a long failure, and running.
TEST(ui_export_footer_reachable) {
    const float W = 1280, H = 720;
    ScaledUi ctx(W, H, 1.5f);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Export;
    p.injectKnit(false);
    p.ui.s.fibres = true;
    CHECK(insidePanel(p, W, H, "Knit at Final quality"));
    UiFixture q("teddy");
    q.ui.stage = viewer::Stage::Export;
    q.injectKnit(true);
    setAttempt(q.ui.results.exportAttempt, longMessage(150), true);
    CHECK(insidePanel(q, W, H, "Export package"));
    q.injectKnit(true, 3);
    CHECK(insidePanel(q, W, H, "Export package"));
    UiFixture r("teddy");
    r.ui.stage = viewer::Stage::Export;
    r.injectKnit(true);
    r.injectRunning(2, "writing the photon_path_tracer package", 0, 0, false);
    CHECK(insidePanel(r, W, H, "Cancel export"));
}

// One home per job: an export running in the Export step has its Cancel in the footer (none in the
// status bar); a knit started from the Export step shows no export progress or Cancel there, its
// Cancel is the status bar's; an export running while another step is shown is cancelled from the
// status bar.
TEST(ui_export_progress_single_home) {
    const float W = 1600, H = 1000;
    UiContext ctx(W, H);
    auto count = [&](const std::vector<ItemRec>& items, const std::string& label, bool panel) {
        int n = 0;
        for (const ItemRec& r : items) n += r.label == label && (r.windowName.rfind("Workflow", 0) == 0) == panel;
        return n;
    };
    UiFixture p("teddy");
    p.injectKnit(true);
    p.injectRunning(2, "writing the photon_path_tracer package", 0, 0, false);
    p.ui.stage = viewer::Stage::Export;
    std::vector<ItemRec> items = settle(p, W, H);
    CHECK(count(items, "Cancel export", true) == 1 && count(items, "Cancel", false) == 0);
    p.ui.stage = viewer::Stage::Knit;
    items = settle(p, W, H);
    CHECK(count(items, "Cancel export", true) == 0 && count(items, "Cancel", false) == 1);
    UiFixture q("teddy");
    q.injectKnit(false);
    q.injectRunning(0, "knitting body", 2, 5, false);
    q.ui.stage = viewer::Stage::Export;
    items = settle(q, W, H);
    std::printf("    knit running in the Export step: Cancel export %d, status bar Cancel %d, waiting %d\n", count(items, "Cancel export", true),
                count(items, "Cancel", false), disabledItem(q, W, H, "Waiting for the running job"));
    CHECK(count(items, "Cancel export", true) == 0 && count(items, "Cancel", false) == 1);
    CHECK(disabledItem(q, W, H, "Waiting for the running job") == 1 && disabledItem(q, W, H, "Knit at Final quality") == -1);
}

// The workflow review's fixes: generated fibres left out are said in the warning colour; the
// detail note names the Knit step's mesh level and what makes the settings Custom; knitting
// Fur: the summary names the mesh and its fur and the fur's fix grows it in the Fur step;
// without a workflow the fix opens the Import step's choice (no job); an existing
// output folder is announced.
TEST(ui_export_workflow_fixes) {
    namespace fs = std::filesystem;
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(true);
    p.ui.stage = viewer::Stage::Export;
    CHECK(viewer::fibreExportWord(p.ui) == "Excluded");
    p.injectFibres();
    std::printf("    generated, excluded: \"%s\"; ready line \"%s\"\n", viewer::fibreExportWord(p.ui).c_str(), viewer::exportReadyLine(p.ui).c_str());
    CHECK(viewer::fibreExportWord(p.ui) == "Excluded: fibres are generated, tick to include");
    CHECK(viewer::exportReadyLine(p.ui) == "Photon Path Tracer package \xC2\xB7 no fibres (generated, not included)");
    std::printf("    detail: %s\n", viewer::exportDetailNote(p.ui).c_str());
    CHECK(viewer::exportDetailNote(p.ui) == "Sets the Knit step's Yarn tessellation to Reduced and at most 7.0 M fibre segments.");
    CHECK(p.setSetting("fibres.segment_budget", "3000000"));
    std::printf("    detail: %s\n", viewer::exportDetailNote(p.ui).c_str());
    CHECK(viewer::exportDetailNote(p.ui).rfind("Custom: the tessellation is Reduced's (Reduced), but the fibre segment budget is 3.0 M instead of 7.0 M", 0) == 0);
    CHECK(p.setSetting("strand.tube_sides", "16") && p.setSetting("strand.tube_ring", "0.2") && p.setSetting("strand.ply_tube_sides", "10") &&
          p.setSetting("strand.ply_tube_ring", "0.15"));
    std::printf("    detail: %s\n", viewer::exportDetailNote(p.ui).c_str());
    CHECK(viewer::exportDetailNote(p.ui).rfind("Custom: the Knit step's Yarn tessellation is Fine, which no level uses", 0) == 0);
    std::string why;
    CHECK(knit::manifestScene("own/manifest.txt", "piece.obj", 24, &why) != nullptr);  // not a preset: fur with its body
    UiFixture hair("pieces");
    viewer::chooseWorkflow(hair.ui, knit::Workflow::Fur);
    hair.ui.stage = viewer::Stage::Export;
    CHECK(!knit::knitWorkflow(hair.ui.s) && (viewer::eligibility(hair.ui).state == knit::ExportEligibility::State::Blocked));
    std::printf("    fur summary: %s\n", hair.exportSummary().c_str());
    CHECK(hair.exportSummary() == "The mesh with its fur, and report.json.");
    hair.ui.exportTo.format = "obj";
    CHECK(hair.exportSummary() == "One .obj and .mtl with the mesh and its fur (as lines), and report.json.");
    hair.ui.exportTo.format = "ppt";
    CHECK(click(hair, W, H, "Grow fur"));
    std::vector<std::string> jobs = hair.submittedJobs(1);
    std::printf("    Grow fur: jobs %s, stage %d\n", jobs.empty() ? "none" : jobs[0].c_str(), static_cast<int>(hair.ui.stage));
    CHECK(jobs.size() == 1 && jobs[0] == "fibres knit=0" && hair.ui.stage == viewer::Stage::Fibres);
    UiFixture c("teddy");
    viewer::openExample(c.ui, "teddy");  // a newly opened example: no workflow yet
    c.ui.stage = viewer::Stage::Export;
    CHECK(click(c, W, H, viewer::kChooseWorkflow));
    CHECK(c.ui.stage == viewer::Stage::Import && !c.ui.job.active && c.submittedJobs(0).empty());
    // The warning: only a folder with something in it.
    const fs::path dir = fs::temp_directory_path() / "knit_ui_test_existing_export";
    std::error_code ec;
    fs::remove_all(dir, ec);
    fs::create_directories(dir, ec);
    p.ui.platform.dialogs = true;
    CHECK(click(p, W, H, "Choose folder ...") && p.answerFolderDialog(dir.string()));
    settle(p, W, H);
    CHECK(viewer::exportExistingNote(p.ui).empty());
    std::ofstream(dir / "mine.txt") << "mine";
    std::printf("    existing: %s\n", viewer::exportExistingNote(p.ui).c_str());
    CHECK(viewer::exportExistingNote(p.ui).rfind("This folder is not empty.", 0) == 0);
    fs::remove_all(dir, ec);
    CHECK(viewer::exportExistingNote(p.ui).empty());
}

// ---- UI text snapshot (simplification parity, 2026-10-28) ----
// Every text the interface renders (ImGui's log captures each RenderText call) and every widget
// with its enabled state, for each example scene, in a few result states, at every step, with
// every section open. Compared exactly with tests/data/ui_text_snapshot.txt, so a refactor that
// changes a label, a state word, a disabled state or a line of help text shows up by name.
// KNIT_UI_SNAPSHOT_WRITE=1 writes the reference instead (an intended UI change, reviewed in the diff).

namespace {

// Replaces every occurrence of `from` in `text` with `to`.
void replaceAll(std::string& text, const std::string& from, const std::string& to) {
    if (from.empty()) return;
    for (size_t at = text.find(from); at != std::string::npos; at = text.find(from, at + to.size()))
        text.replace(at, from.size(), to);
}

// Machine paths in the text (source tree, working folder, temp folder) become fixed tokens, so
// two checkouts compare equal.
std::string maskPaths(std::string text) {
    namespace fs = std::filesystem;
    const std::string paths[][2] = {{fs::path(KNIT_SOURCE_DIR).make_preferred().string(), "<src>"},
                                    {fs::path(KNIT_SOURCE_DIR).generic_string(), "<src>"},
                                    {fs::current_path().string(), "<cwd>"},
                                    {fs::current_path().generic_string(), "<cwd>"},
                                    {fs::temp_directory_path().string(), "<tmp>"},
                                    {fs::temp_directory_path().generic_string(), "<tmp>"}};
    for (const auto& p : paths) replaceAll(text, p[0], p[1]);
    return text;
}

// One state of the interface: its rendered text and its widgets (label, disabled), all sections open.
// A toolbar popover's widgets (label, disabled), opened by its button and closed with Escape.
std::string popoverSnapshot(UiFixture& p, float W, float H, const char* menu, const char* title) {
    if (!clickSuffix(p, W, H, menu)) return std::string("--- ") + title + ": no button\n";
    std::string out = std::string("--- ") + title + "\n";
    for (const ItemRec& r : settle(p, W, H))
        if (r.windowName.rfind("##Popup", 0) == 0 && !r.label.empty()) out += r.label + (r.disabled ? "  [disabled]\n" : "\n");
    pressKey(p, W, H, ImGuiKey_Escape);
    return out;
}

std::string snapshotState(UiFixture& p, float W, float H, const std::string& title) {
    openEverything(p, W, H);
    g_items.clear();
    ImGui::NewFrame();
    ImGui::LogToBuffer();
    viewer::drawSceneUi(p.ui, W, H);
    std::string text = GImGui->LogBuffer.c_str();
    ImGui::LogFinish();
    ImGui::Render();
    std::string out = "### " + title + "\n" + text + "\n--- widgets\n";
    for (const ItemRec& r : g_items)
        if (!r.label.empty()) out += r.label + (r.disabled ? "  [disabled]\n" : "\n");
    return maskPaths(out + popoverSnapshot(p, W, H, "###displayMenu", "Display popover"));
}

// The states of one scene: as opened, then a Final knit, then its wool fibres, then the Fur workflow.
std::string snapshotScene(const char* scene, float W, float H) {
    struct Step {
        viewer::Stage stage;
        const char* name;
    };
    const Step steps[] = {{viewer::Stage::Import, "Import"}, {viewer::Stage::Mark, "Mark"}, {viewer::Stage::Knit, "Knit"},
                          {viewer::Stage::Fibres, "Fibres"}, {viewer::Stage::Export, "Export"}};
    std::string out;
    {
        UiFixture p(scene);
        for (const Step& st : steps) {
            p.ui.stage = st.stage;
            out += snapshotState(p, W, H, std::string(scene) + " / opened / " + st.name);
        }
        out += popoverSnapshot(p, W, H, "###viewMenu", "View popover");
    }
    {
        UiFixture p(scene);
        p.injectKnit(true);
        for (const Step& st : steps) {
            if (st.stage == viewer::Stage::Import || st.stage == viewer::Stage::Mark) continue;
            p.ui.stage = st.stage;
            out += snapshotState(p, W, H, std::string(scene) + " / Final knit / " + st.name);
        }
        p.injectFibres();
        for (const Step& st : steps) {
            if (st.stage != viewer::Stage::Fibres && st.stage != viewer::Stage::Export) continue;
            p.ui.stage = st.stage;
            out += snapshotState(p, W, H, std::string(scene) + " / Final knit + fibres / " + st.name);
        }
    }
    {
        UiFixture p(scene);
        viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);
        for (const Step& st : steps) {
            if (st.stage == viewer::Stage::Import || st.stage == viewer::Stage::Mark) continue;
            p.ui.stage = st.stage;
            out += snapshotState(p, W, H, std::string(scene) + " / fur workflow / " + st.name);
        }
    }
    return out;
}

// The index of the first line where `a` and `b` differ, with both lines.
std::string firstDifference(const std::string& a, const std::string& b) {
    std::istringstream sa(a), sb(b);
    std::string la, lb, section;
    for (int n = 1;; ++n) {
        const bool ga = static_cast<bool>(std::getline(sa, la)), gb = static_cast<bool>(std::getline(sb, lb));
        if (!ga && !gb) return "";
        if (ga && la.rfind("### ", 0) == 0) section = la;
        if (!ga || !gb || la != lb)
            return "line " + std::to_string(n) + " (" + section + "):\n      reference: " + (ga ? la : "<end>") +
                   "\n      now:       " + (gb ? lb : "<end>");
    }
}

} // namespace

namespace {

// The text the interface renders in this state.
std::string renderedText(UiFixture& p, float W, float H) { return snapshotState(p, W, H, "state"); }

// A closed box as an OBJ with a colour, in a fresh folder.
std::string writeBoxObj(const char* folder) {
    namespace fs = std::filesystem;
    const fs::path d = fs::temp_directory_path() / "knit_ui_import" / folder;
    std::error_code ec;
    fs::remove_all(d, ec);
    fs::create_directories(d);
    std::ofstream(d / "box.mtl") << "newmtl paint\nKd 0.6 0.2 0.1\n";
    std::ofstream(d / "box.obj")
        << "mtllib box.mtl\nv 0 0 0\nv 90 0 0\nv 90 90 0\nv 0 90 0\nv 0 0 90\nv 90 0 90\nv 90 90 90\nv 0 90 90\nusemtl paint\n"
           "f 1 4 3\nf 1 3 2\nf 1 2 6\nf 1 6 5\nf 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\nf 2 3 7\nf 2 7 6\nf 5 6 7\nf 5 7 8\nf 1 2 2\n";
    return (d / "box.obj").string();
}

} // namespace

// A knitted import is always rebuilt (2026-11-11): the dialog says how and offers no other surface, no
// repairs, no charts and no knitting in the round.
TEST(ui_import_default_rebuilds) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("default")));
    CHECK((p.ui.import.options.surface == knit::ImportOptions::Surface::Wrapped));
    const std::string text = renderedText(p, W, H);
    CHECK(text.find("rebuilt as one closed surface at stitch scale") != std::string::npos && text.find("Surface preparation") != std::string::npos);
    CHECK(text.find("Keep the file's mesh") == std::string::npos && text.find("in the round") == std::string::npos);
    CHECK(text.find("Repairs change the prepared copy") == std::string::npos && text.find("Knitting charts") == std::string::npos);
    CHECK(text.find("Up to about") != std::string::npos && text.find("stitches at the current stitch size") != std::string::npos);
    viewer::closeImport(p.ui);
}

namespace {

// Settings naming an import, adopted as the viewer does at start, then polled until the session has answered
// (read, then kept, prepared again, or left open for the user).
void adoptUntilSettled(UiFixture& q, const knit::SceneSettings& s) {
    using Ph = knit::ImportSession::Phase;
    q.ui.s = s;
    viewer::adoptSettings(q.ui);
    for (int i = 0; i < 6000; ++i) {
        viewer::pollSceneUi(q.ui);
        const Ph phase = q.ui.import.session.snapshot().phase;
        if (phase != Ph::Reading && phase != Ph::Preparing && !q.ui.import.session.busy()) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    viewer::pollSceneUi(q.ui);  // the last phase integrated
}

} // namespace

// An import named in saved settings is read again at start (adoptSettings): the same file keeps its stored
// entry; a changed one is prepared again and says so; one that can no longer be imported opens the Import
// step with its surface preparation; a missing one uses the prepared copy and says so.
TEST(ui_import_revalidation) {
    namespace fs = std::filesystem;
    using Ph = knit::ImportSession::Phase;
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_revalidation";  // its own: CTest runs cases in parallel
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    const std::string box = writeBoxObj("revalidation");
    UiFixture p("teddy");
    CHECK(p.openImport(box));
    p.ui.import.options.surface = knit::ImportOptions::Surface::KeepFile;  // the file's mesh: no rebuild
    CHECK(p.finishImport());
    const knit::SceneSettings saved = p.ui.s;
    CHECK(saved.scene == "import" && !saved.importIdentity.empty());
    UiFixture same("teddy");
    adoptUntilSettled(same, saved);
    CHECK(same.ui.import.session.snapshot().phase == Ph::Analysed && !same.ui.import.revalidating && !same.ui.import.fileOpen);
    CHECK(same.ui.import.note.empty() && same.ui.s.importIdentity == saved.importIdentity);
    std::ofstream(box, std::ios::app) << "# edited\n";
    UiFixture changed("teddy");
    adoptUntilSettled(changed, saved);
    CHECK(changed.ui.import.note.find("changed since it was imported") != std::string::npos);
    CHECK(changed.ui.import.session.snapshot().phase == Ph::Prepared && changed.ui.s.importIdentity != saved.importIdentity);
    std::ofstream(box) << "v 1 1 1\nv 1 1 1\nv 1 1 1\nf 1 2 3\n";  // all at one point: analysed, but fatal
    UiFixture fatal("teddy");
    adoptUntilSettled(fatal, saved);
    CHECK(fatal.ui.import.fileOpen && fatal.ui.import.session.snapshot().analysis && fatal.ui.import.session.snapshot().analysis->fatal());
    fatal.ui.stage = viewer::Stage::Import;
    fatal.ui.import.options.surface = knit::ImportOptions::Surface::Wrapped;
    const std::string text = renderedText(fatal, W, H);
    CHECK(text.find("Surface preparation") != std::string::npos);
    fs::remove(box, ec);
    UiFixture missing("teddy");
    adoptUntilSettled(missing, saved);
    CHECK(missing.ui.import.note.find("is missing: using the prepared copy") != std::string::npos && !missing.ui.import.fileOpen);
    _putenv_s("KNIT_IMPORT_CACHE", "");
}

// A Fur import keeps the file's mesh and prepares no charts; this box has no texture coordinates, so it cannot be
// knitted, and choosing Knit afterwards shows no Knit footer. pollSceneUi reads it once per import, not per frame.
TEST(ui_import_lacks_charts) {
    namespace fs = std::filesystem;
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_lacks_charts";
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("lacks_charts")));
    CHECK(!p.ui.workflowChosen);  // a newly opened file asks
    viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);
    CHECK(p.ui.import.options.surface == knit::ImportOptions::Surface::KeepFile);
    renderedText(p, W, H);  // the Import step as the user sees it, before Import
    const bool finished = p.finishImport();
    viewer::pollSceneUi(p.ui);
    _putenv_s("KNIT_IMPORT_CACHE", "");
    CHECK(finished && p.ui.s.scene == "import");
    CHECK(p.ui.import.lacksCharts);
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
    CHECK(!viewer::knitFooterShown(p.ui));
}

// Knit -> Fur -> Knit on an import (handover rev. 3, section 4): Knit imports the rebuilt surface; switching to Fur
// prepares the file's own mesh (every setting and the component choices kept) and parks the knit; switching back
// selects the rebuilt surface again from the import cache and brings the parked knit back. While a job runs the
// switch is refused, so its result cannot land in the other workflow.
TEST(ui_import_workflow_switch_selects_source) {
    namespace fs = std::filesystem;
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_workflow_switch";
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("workflow_switch")));
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Other);
    viewer::setCustomImportUnit(p.ui, 0.5);  // a 45 mm box: a quick rebuild
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
    const bool finished = p.finishImport();
    std::printf("    knit import: finished %d, scene %s, surface %s, phase %d, error %s\n", finished, p.ui.s.scene.c_str(), p.ui.s.importSurface.c_str(),
                static_cast<int>(p.ui.import.session.snapshot().phase), p.ui.import.session.snapshot().error.c_str());
    CHECK(finished && p.ui.s.scene == "import" && p.ui.s.importSurface == "wrapped");
    const std::string rebuilt = p.ui.s.importIdentity, plain = p.ui.s.importPlainNodes;
    p.injectKnit(false);
    const std::shared_ptr<const knit::BuiltScene> knitted = p.ui.results.knit;
    p.injectRunning(0, "knitting", 1, 4, false);
    viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);
    CHECK(knit::knitWorkflow(p.ui.s) && p.ui.results.knit == knitted);  // refused while the job runs
    p.ui.job.active = false;
    viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);
    bool prepared = false;
    for (int i = 0; i < 12000 && !prepared; ++i) {
        viewer::pollSceneUi(p.ui);
        prepared = p.ui.s.importSurface == "file";
        if (!prepared) std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    std::printf("    rebuilt %s -> own mesh %s (prepared %d), stage %d\n", rebuilt.c_str(), p.ui.s.importIdentity.c_str(), prepared,
                static_cast<int>(p.ui.stage));
    CHECK(prepared && p.ui.s.importIdentity != rebuilt && viewer::furShown(p.ui) && p.ui.stage == viewer::Stage::Fibres);
    CHECK(p.ui.s.importPlainNodes == plain && !p.ui.results.knit && p.ui.results.parked.knit == knitted);
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
    _putenv_s("KNIT_IMPORT_CACHE", "");
    CHECK(p.ui.s.importIdentity == rebuilt && p.ui.s.importSurface == "wrapped" && p.ui.results.knit == knitted);
}

// The Fur step's disc preview (jobs/fur_preview.h, handover section 12): choosing Fur opens on Discs, the preview
// target offers Discs / Reference sphere / Model, the view shows a published preview's discs over its body with the
// counts (a draft says it is a quarter), the View section's colouring is a viewer choice (never a setting), and the
// knit's preview targets have no Discs.
TEST(ui_fur_disc_preview) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy-brown");
    viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);
    CHECK(p.ui.stage == viewer::Stage::Fibres && p.ui.view.fibres == viewer::FibreView::Discs);
    std::vector<ItemRec> items = frame(p, W, H);
    CHECK(findItem(items, "Discs") && findItem(items, "Reference sphere") && findItem(items, "Model") && findItem(items, "##Colour by"));
    CHECK(click(p, W, H, "Reference sphere") && p.ui.view.fibres == viewer::FibreView::Reference);
    CHECK(click(p, W, H, "Discs") && p.ui.view.fibres == viewer::FibreView::Discs);
    auto preview = std::make_shared<knit::FurPreview>();
    preview->loops = 2;
    preview->fibres = 1;
    preview->keyBelow = 0.25f;
    preview->discPos.assign(18, 0.0f);
    preview->discNormal.assign(18, 0.0f);
    preview->discLocal.assign(18, 0.0f);
    preview->discColour.assign(6, 0xFFFFFFFFu);
    auto body = std::make_shared<knit::BuiltScene>();
    knit::DisplayMesh mesh;
    mesh.pos = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    mesh.idx = {0, 1, 2};
    body->meshes.push_back(mesh);
    preview->body = body;
    p.ui.previews.discShown = preview;
    const viewer::ViewContent c = viewer::viewContent(p.ui);
    std::printf("    caption: %s\n", c.caption.c_str());
    CHECK(c.furDiscs == preview && c.meshes.size() == 1 && c.caption.find("2 loops, 1 fine fibres (a quarter shown while editing)") != std::string::npos);
    bool owned = false;
    for (const auto& s : viewer::meshOwners(p.ui).scenes) owned |= s == body;
    CHECK(owned);  // the body stays uploaded while it is shown
    const std::string before = knit::writeSettings(p.ui.s);
    p.ui.view.discColour = knit::FurDiscColour::Tilt;
    CHECK(knit::writeSettings(p.ui.s) == before);
    UiFixture wool("teddy");
    wool.ui.stage = viewer::Stage::Fibres;
    CHECK(!findItem(frame(wool, W, H), "Discs"));
}

// Fur keeps the file's mesh: the fur grows on the file's surface. The surface preparation goes.
TEST(ui_import_skip_forces_keep) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("skip")));
    CHECK(click(p, W, H, "Fur"));
    const std::string text = renderedText(p, W, H);
    CHECK(!(p.ui.import.options.surface == knit::ImportOptions::Surface::Wrapped));
    CHECK(text.find("Surface preparation") == std::string::npos);
    viewer::closeImport(p.ui);
}

// With the Fur switch off (knit::furEnabled) Knit is the only workflow: an example opens on its knit (asked for at
// once), an opened file asks nothing and shows no Knit / Fur buttons, and the script's choose-fur changes nothing.
TEST(ui_fur_switched_off) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    knit::setFurEnabled(false);
    UiFixture teddy("teddy");
    viewer::openExample(teddy.ui, "teddy-brown");
    const bool exampleOnKnit = teddy.ui.workflowChosen && knit::knitWorkflow(teddy.ui.s) && teddy.ui.stage == viewer::Stage::Knit &&
                               teddy.ui.job.active;
    UiFixture p("teddy");
    const bool opened = p.openImport(writeBoxObj("fur_off"));
    const bool fileChosen = p.ui.workflowChosen && knit::knitWorkflow(p.ui.s);
    const std::string text = renderedText(p, W, H);
    const bool buttonsShown = text.find("What to make") != std::string::npos || text.find(viewer::kChooseWorkflow) != std::string::npos;
    viewer::applyScriptStep(p.ui, "choose-fur");
    const bool stillKnit = knit::knitWorkflow(p.ui.s);
    viewer::closeImport(p.ui);
    knit::setFurEnabled(true);
    CHECK(exampleOnKnit);
    CHECK(opened && fileChosen);
    CHECK(!buttonsShown);
    CHECK(stillKnit);
}

namespace {

// A box 10 wide, 80 tall along y and 20 deep along z, in a fresh folder: its height differs per up axis.
std::string writeTallBoxObj(const char* folder) {
    namespace fs = std::filesystem;
    const fs::path d = fs::temp_directory_path() / "knit_ui_import" / folder;
    std::error_code ec;
    fs::remove_all(d, ec);
    fs::create_directories(d);
    std::ofstream(d / "tall.obj") << "v 0 0 0\nv 10 0 0\nv 10 80 0\nv 0 80 0\nv 0 0 20\nv 10 0 20\nv 10 80 20\nv 0 80 20\n"
                                     "f 1 4 3\nf 1 3 2\nf 1 2 6\nf 1 6 5\nf 3 4 8\nf 3 8 7\nf 4 1 5\nf 4 5 8\nf 2 3 7\nf 2 7 6\nf 5 6 7\nf 5 7 8\n";
    return (d / "tall.obj").string();
}

// Holds the import session's worker at "start" until released; released on scope exit (declare it after the
// fixture, so it goes first and the session's destructor never waits on it).
struct HoldImportWorker {
    std::atomic<bool> reached{false}, released{false};
    std::function<void(const char*)> hook() {
        return [this](const char* at) {
            if (std::string(at) != "start") return;
            reached = true;
            while (!released) std::this_thread::sleep_for(std::chrono::milliseconds(1));
        };
    }
    bool waitReached() {
        for (int i = 0; i < 5000 && !reached; ++i) std::this_thread::sleep_for(std::chrono::milliseconds(1));
        return reached;
    }
    ~HoldImportWorker() { released = true; }
};

} // namespace

// The size choices: a file unit; Custom starts from the current size; "Set height" (a click) starts from the
// height on screen, so the size stays, follows a change of the up axis, and File units brings the unit back.
TEST(ui_import_unit_choices) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(p.openImport(writeTallBoxObj("unit_choices")));
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Cm);
    CHECK(p.ui.import.unit == viewer::UnitChoice::Cm && p.ui.import.options.unitMm == 10 && p.ui.import.unitTouched);
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Inch);
    CHECK(p.ui.import.options.unitMm == 25.4);
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Other);
    CHECK(p.ui.import.unit == viewer::UnitChoice::Other && p.ui.import.customUnitMm == 25.4 && p.ui.import.options.unitMm == 25.4);
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Cm);
    CHECK(click(p, W, H, "Set height") && p.ui.import.sizeByHeight && p.ui.import.unit == viewer::UnitChoice::Cm);
    CHECK(std::abs(p.ui.import.targetHeightMm - 800) < 1e-9 && std::abs(p.ui.import.options.unitMm - 10) < 1e-12);  // the size stays
    viewer::chooseImportUp(p.ui, knit::UpAxis::Z);
    CHECK(std::abs(p.ui.import.options.unitMm - 800.0 / 20) < 1e-12);
    viewer::chooseImportUp(p.ui, knit::UpAxis::Y);
    CHECK(std::abs(p.ui.import.options.unitMm - 800.0 / 80) < 1e-12);
    viewer::setImportHeight(p.ui, 400);
    CHECK(std::abs(p.ui.import.options.unitMm - 5) < 1e-12);
    CHECK(click(p, W, H, "File units") && !p.ui.import.sizeByHeight && p.ui.import.options.unitMm == 10);
    viewer::closeImport(p.ui);
}

// Fur keeps the file's mesh; Knit rebuilds the surface again, and so does a click on Knit while it is chosen. Import
// waits for the choice.
TEST(ui_import_skip_rules) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    using O = knit::ImportOptions;
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("skip_rules")));
    CHECK(!p.ui.workflowChosen && viewer::importFooterView(p.ui).whyDisabled.find(viewer::kChooseWorkflow) != std::string::npos);
    CHECK(click(p, W, H, "Fur") && knit::furWorkflow(p.ui.s) && p.ui.import.options.surface == O::Surface::KeepFile);
    CHECK(click(p, W, H, "Knit") && knit::knitWorkflow(p.ui.s) && p.ui.import.options.surface == O::Surface::Wrapped);
    CHECK(click(p, W, H, "Knit") && p.ui.import.options.surface == O::Surface::Wrapped);
    viewer::closeImport(p.ui);
}

// Cancel while the file is read closes it; Cancel while it is prepared goes back to its choices.
TEST(ui_import_cancel_paths) {
    namespace fs = std::filesystem;
    using Ph = knit::ImportSession::Phase;
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_cancel_paths";
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    const std::string box = writeBoxObj("cancel_paths");
    UiFixture p("teddy");
    HoldImportWorker reading;
    p.ui.import.session.setTestHook(reading.hook());
    viewer::openMeshFile(p.ui, box);
    CHECK(reading.waitReached());
    CHECK(click(p, W, H, "Cancel"));
    CHECK(!p.ui.import.fileOpen && p.ui.import.session.snapshot().phase == Ph::Idle);
    reading.released = true;
    p.ui.import.session.setTestHook(nullptr);
    CHECK(p.openImport(box));
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);  // Import waits for the choice
    p.ui.import.options.surface = knit::ImportOptions::Surface::KeepFile;  // quick if it ran
    HoldImportWorker preparing;
    p.ui.import.session.setTestHook(preparing.hook());
    CHECK(click(p, W, H, "Rebuild & import"));
    CHECK(preparing.waitReached() && p.ui.import.session.snapshot().phase == Ph::Preparing);
    CHECK(click(p, W, H, "Cancel preparing"));
    CHECK(p.ui.import.session.snapshot().phase == Ph::Analysed && p.ui.import.fileOpen);
    preparing.released = true;
    p.ui.import.session.setTestHook(nullptr);
    _putenv_s("KNIT_IMPORT_CACHE", "");
    viewer::closeImport(p.ui);
}

namespace {

} // namespace

// A wrap that fails falls back to nothing: the dialog stays, still rebuilding, and the user decides.
TEST(ui_import_failed_wrap_keeps_the_choice) {
    namespace fs = std::filesystem;
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_failed";
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("failed")));
    p.ui.import.options.unitMm = 0.01;  // under a mm: the blur leaves nothing of it
    CHECK((p.ui.import.options.surface == knit::ImportOptions::Surface::Wrapped));
    const bool finished = p.finishImport();
    _putenv_s("KNIT_IMPORT_CACHE", "");
    CHECK(!finished);
    CHECK((p.ui.import.options.surface == knit::ImportOptions::Surface::Wrapped));  // nothing falls back by itself: the user chooses
    viewer::closeImport(p.ui);
}

// After a failed rebuild the user may keep the file's own surface: offered (never pressed by itself), it imports the file
// as it is, once. Defect: an automatic fallback (ui_import_failed_wrap_keeps_the_choice fails), or no way past the failure
// for an open tube (finished stays false here).
TEST(ui_import_failed_wrap_keep_file_on_request) {
    namespace fs = std::filesystem;
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_failed_keep";
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("failed_keep")));
    p.ui.import.options.unitMm = 0.01;  // under a mm: the blur leaves nothing of it
    CHECK(!p.finishImport());
    CHECK(viewer::canKeepFileAfterFailedRebuild(p.ui));
    CHECK(renderedText(p, W, 12000).find("Keep the file's own surface") != std::string::npos);
    viewer::keepFileAfterFailedRebuild(p.ui);
    CHECK((p.ui.import.options.surface == knit::ImportOptions::Surface::KeepFile));
    CHECK(!viewer::canKeepFileAfterFailedRebuild(p.ui));
    bool finished = false;
    for (int i = 0; i < 12000 && !finished; ++i) {
        viewer::pollSceneUi(p.ui);
        finished = !p.ui.import.fileOpen;
        if (p.ui.import.session.snapshot().phase == knit::ImportSession::Phase::Failed) break;
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
    _putenv_s("KNIT_IMPORT_CACHE", "");
    CHECK(finished);
    CHECK(p.ui.s.importSurface == "file");
}

// A rebuilt import opens with its baked colours, per stitch.
TEST(ui_import_opens_with_colourwork) {
    namespace fs = std::filesystem;
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache";
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("colourwork")));
    p.ui.import.options.unitMm = 1.0;  // 90 mm: the unit the dialog guesses makes this a 900 mm box of 2 million quads
    const bool finished = p.finishImport();
    _putenv_s("KNIT_IMPORT_CACHE", "");
    CHECK(finished);
    CHECK(p.ui.s.yarnColourSource == "texture");
    p.ui.stage = viewer::Stage::Import;
    const std::string text = renderedText(p, W, 12000);
    CHECK(text.find("rebuilt as a closed surface") != std::string::npos && text.find("rebuilt by") == std::string::npos);
}

namespace {

// Two named 40 mm boxes side by side, 60 mm apart, with the file's Y up.
std::string writeTwoBoxObj(const char* folder) {
    namespace fs = std::filesystem;
    const fs::path d = fs::temp_directory_path() / "knit_ui_import" / folder;
    std::error_code ec;
    fs::remove_all(d, ec);
    fs::create_directories(d);
    std::string obj;
    for (const int x : {0, 100})
        for (const int z : {0, 40})
            for (const int y : {0, 40})
                for (const int dx : {0, 40}) obj += "v " + std::to_string(x + dx) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
    for (int box = 0; box < 2; ++box) {
        obj += box == 0 ? "o left\n" : "o right\n";
        const int faces[12][3] = {{1, 3, 4}, {1, 4, 2}, {5, 6, 8}, {5, 8, 7}, {1, 2, 6}, {1, 6, 5}, {3, 7, 8}, {3, 8, 4}, {1, 5, 7}, {1, 7, 3}, {2, 4, 8}, {2, 8, 6}};
        for (const int (&f)[3] : faces) obj += "f " + std::to_string(f[0] + 8 * box) + " " + std::to_string(f[1] + 8 * box) + " " + std::to_string(f[2] + 8 * box) + "\n";
    }
    std::ofstream(d / "boxes.obj") << obj;
    return (d / "boxes.obj").string();
}

} // namespace

namespace {

// The nodes along the viewport's middle row, without repeats: -1 for the background.
std::vector<int> nodesAcross(const viewer::SceneUi& ui, const viewer::Camera& camera, const viewer::ViewPixel& viewport) {
    std::vector<int> seen;
    for (double x = 0; x < viewport.width; x += 2) {
        const std::optional<int> node = viewer::nodeUnderPixel(ui, camera, {x, viewport.height / 2, viewport.width, viewport.height});
        const int n = node.value_or(-1);
        if (seen.empty() || seen.back() != n) seen.push_back(n);
    }
    return seen;
}

// The first pixel of the middle row that shows this node (-1: the background).
viewer::ViewPixel pixelOfNode(const viewer::SceneUi& ui, const viewer::Camera& camera, const viewer::ViewPixel& viewport, int wanted) {
    for (double x = 0; x < viewport.width; x += 2) {
        const viewer::ViewPixel pixel{x, viewport.height / 2, viewport.width, viewport.height};
        if (viewer::nodeUnderPixel(ui, camera, pixel).value_or(-1) == wanted) return pixel;
    }
    return {};
}

// A camera on the -y side of the library's axes, framing what the viewport shows: x runs to the right of the picture.
viewer::Camera frontCamera(const viewer::SceneUi& ui) {
    viewer::Camera camera;
    camera.look({0, 1, 0}, {0, 0, 1});
    viewer::frameBox(camera, viewer::viewContent(ui).bounds);
    return camera;
}

} // namespace

// Stage 2: the object's components, in the list, with the quick selections; the other steps wait; Cancel is the way back.
TEST(ui_import_two_stages) {
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Import;
    std::vector<ItemRec> items = settle(p, W, H);
    CHECK(findItem(items, "Open mesh (OBJ, PLY) ...") != nullptr);
    CHECK(p.openImport(writeTwoBoxObj("stages")));
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
    items = settle(p, W, H);
    for (const char* label : {"left", "right", "Knit left", "Knit right", "Knit all", "Keep all as mesh", "Rebuild & import", "Cancel"})
        CHECK(findItem(items, label) != nullptr);
    CHECK(findItem(items, "Open mesh (OBJ, PLY) ...") == nullptr);
    int locked = 0;
    for (const ItemRec& r : items) locked += r.label == "##stage" && r.disabled;
    CHECK(locked == viewer::kStages - 1);

    CHECK(click(p, W, H, "Knit left"));
    CHECK((p.ui.import.options.plainNodes == std::set<std::string>{"left"}));
    CHECK(viewer::importFooterView(p.ui).summary.rfind(std::string("1 to knit") + viewer::kDot + "1 kept as mesh", 0) == 0);
    CHECK(click(p, W, H, "Keep all as mesh"));
    CHECK((p.ui.import.options.plainNodes == std::set<std::string>{"left", "right"}));
    CHECK(findItem(settle(p, W, H), "Rebuild & import")->disabled);  // nothing to knit: the pipeline refuses it
    CHECK(renderedText(p, W, H).find("Nothing to knit: tick a component.") != std::string::npos);
    CHECK(click(p, W, H, "Knit all"));
    CHECK(p.ui.import.options.plainNodes.empty());
    CHECK(!findItem(settle(p, W, H), "Rebuild & import")->disabled);
    CHECK(click(p, W, H, "Cancel"));
    CHECK(!p.ui.import.fileOpen);
    items = settle(p, W, H);
    CHECK(findItem(items, "Open mesh (OBJ, PLY) ...") != nullptr);
    locked = 0;
    for (const ItemRec& r : items) locked += r.label == "##stage" && r.disabled;
    CHECK(locked == 0);
}

// A click on the object focuses its component (in either up axis, the nearer box first) and asks for its row; the
// ticks are the checkboxes' only. The background clears the focus.
TEST(ui_import_pick_focuses) {
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(p.openImport(writeTwoBoxObj("click")));
    p.ui.import.options.unitMm = 1;
    const viewer::ViewPixel viewport{0, 0, 800, 600};
    for (const knit::UpAxis up : {knit::UpAxis::Y, knit::UpAxis::Z}) {
        p.ui.import.options.up = up;
        CHECK((nodesAcross(p.ui, frontCamera(p.ui), viewport) == std::vector<int>{-1, 0, -1, 1, -1}));
    }
    p.ui.import.options.up = knit::UpAxis::Y;
    const viewer::Camera camera = frontCamera(p.ui);
    viewer::pickNode(p.ui, camera, pixelOfNode(p.ui, camera, viewport, 1));
    CHECK(p.ui.import.focusedNode == "right" && p.ui.import.revealFocus && p.ui.import.options.plainNodes.empty());
    int focused = 0;
    for (const knit::DisplayMesh& d : p.ui.import.preview->meshes) focused += viewer::previewNodeIsFocused(p.ui, d) ? 1 : 0;
    CHECK(focused >= 1);
    p.ui.import.options.up = knit::UpAxis::Z;
    const viewer::Camera zCamera = frontCamera(p.ui);
    viewer::pickNode(p.ui, zCamera, pixelOfNode(p.ui, zCamera, viewport, 0));
    CHECK(p.ui.import.focusedNode == "left" && p.ui.import.options.plainNodes.empty());
    p.ui.import.options.up = knit::UpAxis::Y;
    viewer::pickNode(p.ui, camera, pixelOfNode(p.ui, camera, viewport, -1));  // the background clears the focus
    CHECK(p.ui.import.focusedNode.empty() && !p.ui.import.revealFocus && p.ui.import.options.plainNodes.empty());

    // Looking along the row of boxes, the nearer one answers.
    viewer::Camera along = camera;
    along.look({1, 0, 0}, {0, 0, 1});
    CHECK(viewer::nodeUnderPixel(p.ui, along, {400, 300, 800, 600}).value_or(-1) == 0);
    along.look({-1, 0, 0}, {0, 0, 1});
    CHECK(viewer::nodeUnderPixel(p.ui, along, {400, 300, 800, 600}).value_or(-1) == 1);
    viewer::closeImport(p.ui);
}

TEST(ui_import_drag_is_no_click) {
    viewer::Camera c;
    viewer::CameraDrag drag;
    viewer::MouseState m;
    m.overViewport = true;
    m.x = 100;
    m.y = 50;
    drag.update(c, m, 500);
    m.down[0] = true;
    drag.update(c, m, 500);
    CHECK(!drag.takeClick());  // pressed, not yet released
    m.x = 102;                 // a hand's tremor
    drag.update(c, m, 500);
    m.down[0] = false;
    drag.update(c, m, 500);
    const std::optional<knit::Vec2> click = drag.takeClick();
    CHECK(click && click->x == 102 && click->y == 50);
    CHECK(!drag.takeClick());  // once

    m.down[0] = true;  // an orbit that returns to where it started
    drag.update(c, m, 500);
    m.x = 160;
    drag.update(c, m, 500);
    m.x = 102;
    drag.update(c, m, 500);
    m.down[0] = false;
    drag.update(c, m, 500);
    CHECK(!drag.takeClick());

    m.overViewport = false;  // a press over the interface
    m.down[0] = true;
    drag.update(c, m, 500);
    m.down[0] = false;
    drag.update(c, m, 500);
    CHECK(!drag.takeClick());
}

// The ray of a pixel passes through what the matrices put on that pixel.
TEST(ui_view_ray_matches_the_camera) {
    viewer::Camera c;
    c.target = {3, 4, 5};
    c.dist = 50;
    c.look({1, 2, -1}, {0, 0, 1});
    const double width = 800, height = 600;
    const viewer::Mat4 m = viewer::multiply(c.proj(width / height), c.view());
    const knit::Vec3 point{7, -2, 9};
    const double p[3] = {point.x, point.y, point.z};
    double clip[4];
    for (int row = 0; row < 4; ++row) clip[row] = m.m[row] * p[0] + m.m[row + 4] * p[1] + m.m[row + 8] * p[2] + m.m[row + 12];
    const viewer::ViewPixel pixel{(clip[0] / clip[3] + 1) / 2 * width, (1 - clip[1] / clip[3]) / 2 * height, width, height};
    const viewer::Ray ray = viewer::viewRay(c, pixel);
    const knit::Vec3 away = point - ray.origin;
    CHECK(knit::norm(knit::cross(away, ray.direction)) / knit::norm(ray.direction) < 1e-9);
}

// Everything stage 2 draws has an id of its own, whatever the list holds.
TEST(ui_import_inventory) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(p.openImport(writeTwoBoxObj("inventory")));
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
    const std::vector<ItemRec> items = settle(p, W, H);
    for (const char* label : {"Knit", "Fur", "Knit all", "Keep all as mesh", "left", "right", "Knit left", "Knit right", "File units",
                              "Set height", "##fileUnit", "##upAxis", "Preparation details", "File details & compatibility",
                              "Rebuild & import", "Cancel"})
        CHECK(findItem(items, label) != nullptr);
    const std::string text = renderedText(p, W, H);
    CHECK(text.find("Components") != std::string::npos && text.find("2 components") != std::string::npos);
    viewer::closeImport(p.ui);
}

TEST(ui_import_footer_reachable) {
    const float W = 1280, H = 720;
    ScaledUi ctx(W, H, 1.5f);
    UiFixture p("teddy");
    CHECK(p.openImport(writeTwoBoxObj("footer")));
    CHECK(insidePanel(p, W, H, "Rebuild & import"));
    CHECK(insidePanel(p, W, H, "Cancel"));
    viewer::closeImport(p.ui);
}

namespace {

// `parts` separate 10 mm cubes in a row along x, 20 mm apart, without names: the file's components are "Part 1" ...
std::string writeManyPartsObj(const char* folder, int parts) {
    namespace fs = std::filesystem;
    const fs::path d = fs::temp_directory_path() / "knit_ui_import" / folder;
    std::error_code ec;
    fs::remove_all(d, ec);
    fs::create_directories(d);
    std::string obj;
    const int faces[12][3] = {{1, 3, 4}, {1, 4, 2}, {5, 6, 8}, {5, 8, 7}, {1, 2, 6}, {1, 6, 5}, {3, 7, 8}, {3, 8, 4}, {1, 5, 7}, {1, 7, 3}, {2, 4, 8}, {2, 8, 6}};
    for (int k = 0; k < parts; ++k) {
        for (const int z : {0, 10})
            for (const int y : {0, 10})
                for (const int dx : {0, 10}) obj += "v " + std::to_string(20 * k + dx) + " " + std::to_string(y) + " " + std::to_string(z) + "\n";
        for (const int (&f)[3] : faces) obj += "f " + std::to_string(f[0] + 8 * k) + " " + std::to_string(f[1] + 8 * k) + " " + std::to_string(f[2] + 8 * k) + "\n";
    }
    std::ofstream(d / "parts.obj") << obj;
    return (d / "parts.obj").string();
}

// A flat 40 x 40 square in the file's xz plane: nothing along Y, the OBJ's up axis.
std::string writeFlatSquareObj(const char* folder) {
    namespace fs = std::filesystem;
    const fs::path d = fs::temp_directory_path() / "knit_ui_import" / folder;
    std::error_code ec;
    fs::remove_all(d, ec);
    fs::create_directories(d);
    std::ofstream(d / "flat.obj") << "v 0 0 0\nv 40 0 0\nv 40 0 40\nv 0 0 40\nf 1 2 3\nf 1 3 4\n";
    return (d / "flat.obj").string();
}

// The stage panel's window (the form that scrolls), from an item drawn in it.
ImGuiWindow* stagePanelOf(const std::vector<ItemRec>& items, const std::string& label) {
    const ItemRec* r = findItem(items, label);
    ImGuiWindow* w = r ? r->window : nullptr;
    while (w && std::string(w->Name).find("##stagePanel") == std::string::npos) w = w->ParentWindow;
    return w;
}

} // namespace

// A long list in a small window: the list scrolls by itself, and the count, the bulk buttons and the footer's
// actions stay in reach.
TEST(ui_import_long_list_footer_reachable) {
    const float W = 900, H = 720;
    ScaledUi ctx(W, H, 1.5f);
    UiFixture p("teddy");
    CHECK(p.openImport(writeManyPartsObj("long_list", 40)));
    CHECK(static_cast<int>(p.ui.import.session.snapshot().mesh->nodes.size()) == 40);
    CHECK(insidePanel(p, W, H, "Knit all"));
    CHECK(insidePanel(p, W, H, "Rebuild & import"));
    CHECK(insidePanel(p, W, H, "Cancel"));
    CHECK(click(p, W, H, "Keep all as mesh") && p.ui.import.options.plainNodes.size() == 40);  // the rows below the list's view too
    viewer::closeImport(p.ui);
}

// A pick reveals its row through both scroll containers: the form scrolled to its end, Components above the
// view, the row of part 37 (below the list's view) comes into sight, and the request is used up.
TEST(ui_import_pick_reveals_row) {
    const float W = 1280, H = 720;
    ScaledUi ctx(W, H, 1.5f);
    UiFixture p("teddy");
    CHECK(p.openImport(writeManyPartsObj("reveal", 40)));
    ImGuiWindow* form = stagePanelOf(settle(p, W, H), "Knit all");
    CHECK(form != nullptr);
    if (!form) return;
    form->ScrollTarget.y = 1e6f;  // the end of the form (FLT_MAX means "no target" to ImGui)
    form->ScrollTargetCenterRatio.y = 0;
    settle(p, W, H);
    CHECK(!insidePanel(p, W, H, "Knit all"));  // Components is above the view
    p.ui.import.focusedNode = "Part 37";       // what pickNode sets for a hit on that cube
    p.ui.import.revealFocus = true;
    frame(p, W, H);
    CHECK(p.ui.import.revealFocus);  // the list was off screen: the form scrolls first, the request stays
    settle(p, W, H);
    CHECK(!p.ui.import.revealFocus);
    const std::vector<ItemRec> items = settle(p, W, H);
    const ItemRec* row = findItem(items, "Part 37");
    CHECK(row != nullptr);
    if (!row) return;
    const ImRect list = row->window->InnerClipRect, shown = form->InnerClipRect;
    std::printf("    row y %.0f..%.0f, list y %.0f..%.0f, form y %.0f..%.0f\n", row->bb.Min.y, row->bb.Max.y, list.Min.y, list.Max.y, shown.Min.y, shown.Max.y);
    CHECK(row->bb.Min.y >= list.Min.y - 1 && row->bb.Max.y <= list.Max.y + 1);   // inside the list's view
    CHECK(row->bb.Min.y >= shown.Min.y - 1 && row->bb.Max.y <= shown.Max.y + 1); // and the form's
    viewer::closeImport(p.ui);
}

// The readouts are the import's own mapping: x the width; Y up: z the depth, y the height; Z up: y the depth.
TEST(ui_import_size_readouts) {
    UiFixture p("teddy");
    CHECK(p.openImport(writeTallBoxObj("readouts")));
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Mm);
    viewer::ImportSize s = viewer::importSize(p.ui);
    CHECK(s.valid && s.width == 10 && s.depth == 20 && s.height == 80);
    viewer::chooseImportUp(p.ui, knit::UpAxis::Z);
    s = viewer::importSize(p.ui);
    CHECK(s.valid && s.width == 10 && s.depth == 80 && s.height == 20);
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Cm);
    s = viewer::importSize(p.ui);
    CHECK(s.width == 100 && s.depth == 800 && s.height == 200);
    CHECK(viewer::importFooterView(p.ui).summary.find("200.0 mm high") != std::string::npos);
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Other);
    viewer::setCustomImportUnit(p.ui, 0);
    CHECK(!viewer::importSize(p.ui).valid);
    viewer::closeImport(p.ui);
}

// A click on a name focuses it, a click on its checkbox changes only its tick.
TEST(ui_import_focus_and_tick_are_separate) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(p.openImport(writeTwoBoxObj("focus_tick")));
    CHECK(click(p, W, H, "right") && p.ui.import.focusedNode == "right" && p.ui.import.options.plainNodes.empty());
    CHECK(click(p, W, H, "Knit left") && p.ui.import.focusedNode == "right");
    CHECK((p.ui.import.options.plainNodes == std::set<std::string>{"left"}));
    CHECK(renderedText(p, W, H).find(std::string("Shown in the preview") + viewer::kDot + "right") != std::string::npos);
    viewer::closeImport(p.ui);
}

// Opening another file starts a fresh draft: no focus, no reveal, the size by file units and suggested again.
TEST(ui_import_second_file_resets_draft) {
    UiFixture p("teddy");
    CHECK(p.openImport(writeTwoBoxObj("reset_first")));
    viewer::focusNode(p.ui, "left");
    p.ui.import.revealFocus = true;
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Other);
    viewer::setCustomImportUnit(p.ui, 3);
    viewer::chooseSizeByHeight(p.ui, true);
    CHECK(p.ui.import.sizeByHeight);
    CHECK(p.openImport(writeTallBoxObj("reset_second")));
    CHECK(!p.ui.import.sizeByHeight && p.ui.import.focusedNode.empty() && !p.ui.import.revealFocus);
    CHECK(!p.ui.import.unitTouched && p.ui.import.unit != viewer::UnitChoice::Other && p.ui.import.customUnitMm == 1);
    viewer::closeImport(p.ui);
}

// A size that cannot be used disables the import and says why: a Custom unit of 0 or less, a height on a file flat
// along the up axis. File units still take that flat file (height 0 mm).
TEST(ui_import_invalid_size) {
    UiFixture p("teddy");
    CHECK(p.openImport(writeTallBoxObj("invalid_custom")));
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Other);
    for (const double bad : {0.0, -1.0}) {
        viewer::setCustomImportUnit(p.ui, bad);
        const viewer::ImportFooterView v = viewer::importFooterView(p.ui);
        CHECK(p.ui.import.options.unitMm == 0 && v.whyDisabled == "Enter a positive size for one file unit.");
    }
    viewer::closeImport(p.ui);
    UiFixture flat("teddy");
    CHECK(flat.openImport(writeFlatSquareObj("invalid_flat")));
    viewer::chooseWorkflow(flat.ui, knit::Workflow::Knit);
    viewer::chooseImportUnit(flat.ui, viewer::UnitChoice::Mm);
    viewer::ImportSize s = viewer::importSize(flat.ui);
    CHECK(s.valid && s.height == 0 && s.width == 40 && viewer::importFooterView(flat.ui).whyDisabled.empty());
    viewer::chooseSizeByHeight(flat.ui, true);
    CHECK(flat.ui.import.options.unitMm == 0);
    CHECK(viewer::importFooterView(flat.ui).whyDisabled == "The object is flat along the up axis: set its size with File units.");
    viewer::chooseSizeByHeight(flat.ui, false);
    CHECK(flat.ui.import.options.unitMm == 1 && viewer::importFooterView(flat.ui).whyDisabled.empty());
    viewer::closeImport(flat.ui);
}

// Custom survives Set height: Custom 2 mm, Set height, a target of 50 mm, File units: 2 mm again.
TEST(ui_import_custom_survives_height) {
    UiFixture p("teddy");
    CHECK(p.openImport(writeTallBoxObj("custom_height")));
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Other);
    viewer::setCustomImportUnit(p.ui, 2);
    viewer::chooseSizeByHeight(p.ui, true);
    CHECK(std::abs(p.ui.import.targetHeightMm - 160) < 1e-9);
    viewer::setImportHeight(p.ui, 50);
    CHECK(std::abs(p.ui.import.options.unitMm - 50.0 / 80) < 1e-12);
    viewer::chooseSizeByHeight(p.ui, false);
    CHECK(p.ui.import.options.unitMm == 2 && p.ui.import.unit == viewer::UnitChoice::Other && p.ui.import.customUnitMm == 2);
    viewer::closeImport(p.ui);
}

// The suggestion names the unit the bounds suggest (this 90-unit box: cm), until a unit is chosen.
TEST(ui_import_suggestion_names_unit) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("suggestion")));
    CHECK(p.ui.import.unit == viewer::UnitChoice::Cm);
    CHECK(renderedText(p, W, H).find("Centimetres suggested from the model bounds") != std::string::npos);
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Mm);
    CHECK(renderedText(p, W, H).find("suggested from the model bounds") == std::string::npos);
    viewer::closeImport(p.ui);
}

// A mixed import: the unticked component is kept as the file has it, beside the rebuilt one.
TEST(ui_import_mixed_keeps_plain_component) {
    namespace fs = std::filesystem;
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_mixed";
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    UiFixture p("teddy");
    const bool opened = p.openImport(writeTwoBoxObj("mixed"));
    viewer::chooseImportUnit(p.ui, viewer::UnitChoice::Other);
    viewer::setCustomImportUnit(p.ui, 2);  // 80 mm boxes: as quick to rebuild as the colourwork test's 90 mm box
    viewer::toggleNode(p.ui, "left");
    const bool finished = opened && p.finishImport();
    _putenv_s("KNIT_IMPORT_CACHE", "");
    CHECK(finished);
    CHECK(p.ui.s.scene == "import" && p.ui.s.importPlainNodes == "left" && p.ui.s.importSurface == "wrapped");
}

// The window's start with nothing named (user, 2026-11-01): the Import step and no object, no job,
// and the steps ask for an object. Each first object gets the settings it gets from the default scene.
TEST(ui_start_without_object) {
    namespace fs = std::filesystem;
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    viewer::startWithoutObject(p.ui);
    CHECK(p.ui.stage == viewer::Stage::Import && viewer::importState(p.ui).word == "No object");
    CHECK(viewer::viewContent(p.ui).caption == "No object: open a mesh (Import)" && static_cast<int>(viewer::viewContent(p.ui).meshes.size()) == 0);
    CHECK(renderedText(p, W, 12000).find("No object yet.") != std::string::npos);
    for (int stage = 1; stage < viewer::kStages; ++stage) {
        p.ui.stage = static_cast<viewer::Stage>(stage);
        CHECK(renderedText(p, W, 12000).find("Import an object first (Import).") != std::string::npos);
    }
    CHECK(!p.ui.job.active);
    for (const knit::SceneEntry& e : knit::sceneRegistry()) {
        const char* id = e.id.c_str();
        UiFixture fromMan("teddy"), fromNothing("teddy");
        viewer::startWithoutObject(fromNothing.ui);
        viewer::openExample(fromMan.ui, id);
        viewer::openExample(fromNothing.ui, id);
        CHECK(knit::writeSettings(fromNothing.ui.s) == knit::writeSettings(fromMan.ui.s));
    }
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_no_object";  // its own: CTest runs cases in parallel
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    const std::string box = writeBoxObj("no_object");
    UiFixture fromMan("teddy"), fromNothing("teddy");
    viewer::startWithoutObject(fromNothing.ui);
    CHECK(fromMan.openImport(box) && fromNothing.openImport(box));
    fromMan.ui.import.options.surface = knit::ImportOptions::Surface::KeepFile;  // the file's mesh: no rebuild, which takes minutes here
    fromNothing.ui.import.options.surface = knit::ImportOptions::Surface::KeepFile;
    const bool finished = fromMan.finishImport() && fromNothing.finishImport();
    _putenv_s("KNIT_IMPORT_CACHE", "");
    CHECK(finished);
    CHECK(knit::writeSettings(fromNothing.ui.s) == knit::writeSettings(fromMan.ui.s));
}

// Every preset opens with its own settings (SceneEntry::defaultSettings), whatever opens it and whatever was open
// before: the command line's order, viewerDefaults, Import > Examples from the defaults and from a scene changed on
// purpose (a yarn-sized stitch, another diameter, Worn-like islands, other fibres), in Draft and after Knit at Final
// quality. The effective value is compared with the listed one, parsed the same way. Leaving the scene, its own
// settings stay with it: the next scene starts from its own values, except for a listed key the user changed. No
// preset lists a key the profile owns (those Draft and Final set differently).
TEST(ui_scene_settings_defaults_intended) {
    UiContext ctx(1280, 1000);
    auto wrong = [](const knit::SceneSettings& s, const knit::SceneEntry& e, const char* path) {
        std::string m;
        for (const auto& [key, value] : e.defaultSettings) {
            knit::SceneSettings want = s;
            std::string err;
            if (!knit::setSetting(want, key, value, &err)) m += std::string(" ") + key + " (refused: " + err + ")";
            else if (!knit::sameSetting(s, want, key)) m += std::string(" ") + key;
        }
        if (!m.empty()) std::printf("    %s, %s: not as listed:%s\n", e.id.c_str(), path, m.c_str());
        return m;
    };
    int listed = 0;
    for (const knit::SceneEntry& e : knit::sceneRegistry()) {
        listed += static_cast<int>(e.defaultSettings.size());
        knit::SceneSettings draft = knit::viewerDefaults(e), final = draft;
        knit::applyPreviewProfile(draft);
        knit::applyFullProfile(final, e);
        const std::vector<std::string> owned = knit::diffSettings(draft, final);
        for (const auto& [key, value] : e.defaultSettings)
            if (std::find(owned.begin(), owned.end(), key) != owned.end()) {
                std::printf("    %s lists %s, which the profile owns\n", e.id.c_str(), key.c_str());
                CHECK(!"a preset lists a profile-owned key");
            }
        knit::SceneSettings cli;
        cli.scene = e.id;
        knit::applySceneDefaults(cli, e);
        CHECK(wrong(cli, e, "command line").empty());
        CHECK(wrong(knit::viewerDefaults(e), e, "viewerDefaults").empty());
        UiFixture fresh("teddy");
        viewer::openExample(fresh.ui, e.id);
        CHECK(wrong(fresh.ui.s, e, "opened from the defaults").empty());
        UiFixture changed("teddy");
        CHECK(changed.setSetting("strand.stitch_size", "yarn") && changed.setSetting("strand.diameter", "1.5"));
        CHECK(changed.setSetting("fibres.islands.rate", "0.06") && changed.setSetting("fibres.density", "60"));
        CHECK(changed.setSetting("fibres.random.spread_deg", "180") && changed.setSetting("fibres.directed.length_max", "6"));
        viewer::openExample(changed.ui, e.id);
        CHECK(wrong(changed.ui.s, e, "opened from changed settings (Draft)").empty());
        viewer::knitAtFinalQuality(changed.ui);
        CHECK(knit::isFinalProfile(changed.ui.s, e));
        CHECK(wrong(changed.ui.s, e, "after Knit at Final quality").empty());
        // Leaving it for the other preset: the one listed key the user changed goes along, the others do not.
        if (e.defaultSettings.empty()) continue;
        const std::string other = e.id == "teddy" ? "wool-swatch" : "teddy";
        UiFixture leaving("teddy");
        viewer::openExample(leaving.ui, e.id);
        std::string edited;  // the first listed key that takes a number
        for (const auto& [key, value] : e.defaultSettings)
            if (edited.empty() && leaving.setSetting(key, "0.123")) edited = key;
        CHECK(!edited.empty());
        const std::string mine = leaving.setting(edited);
        viewer::openExample(leaving.ui, other);
        const knit::SceneSettings start = knit::viewerDefaults(*knit::findScene(other));
        std::string carried;
        for (const auto& [key, value] : e.defaultSettings)
            if (key != edited && !knit::sameSetting(leaving.ui.s, start, key)) carried += " " + key;
        if (!carried.empty()) std::printf("    %s -> %s: kept the scene's own%s\n", e.id.c_str(), other.c_str(), carried.c_str());
        CHECK(carried.empty() && leaving.setting(edited) == mine);
    }
    std::printf("    %d listed settings checked on 5 paths\n", listed);
    CHECK(listed > 0);
}

namespace {

// The examples whose models ship with the source (scenes/wool-swatch); the snapshot's Examples list shows only
// these in every checkout, so its text does not depend on which other models are on disk.
const std::vector<std::string> kShippedExamples = {"wool-swatch", "wool-swatch-blue"};

// The snapshot of `scenes` compared with (or, with KNIT_UI_SNAPSHOT_WRITE=1, written to) tests/data/<file>. The
// Fur workflow is switched on for it (the snapshot covers that workflow's text too) and back afterwards.
void compareSnapshot(const std::vector<std::string>& scenes, const char* file) {
    namespace fs = std::filesystem;
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    const bool furWas = knit::furEnabled();
    knit::setFurEnabled(true);
    knit::restrictExamplesForTests(kShippedExamples);
    std::string now;
    for (const std::string& id : scenes) now += snapshotScene(id.c_str(), W, H);
    knit::restrictExamplesForTests({});
    knit::setFurEnabled(furWas);
    const fs::path reference = fs::path(KNIT_SOURCE_DIR) / "tests" / "data" / file;
    if (const auto write = knit::envVar("KNIT_UI_SNAPSHOT_WRITE"); write && *write == "1") {
        std::ofstream(reference, std::ios::binary) << now;
        std::printf("    wrote %s (%zu bytes)\n", reference.string().c_str(), now.size());
        return;
    }
    std::ifstream in(reference, std::ios::binary);
    if (!in) SKIP("no reference (write it with KNIT_UI_SNAPSHOT_WRITE=1)");
    const std::string want((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const std::string diff = firstDifference(want, now);
    if (!diff.empty()) {
        const std::string actual = fs::path(file).stem().string() + ".actual.txt";
        std::ofstream(actual, std::ios::binary) << now;
        std::printf("    first difference at %s\n    the whole text: %s\n", diff.c_str(), actual.c_str());
    }
    CHECK(diff.empty());
}

} // namespace

// The examples whose models ship with the source: the same text in every checkout.
TEST(ui_text_snapshot) { compareSnapshot(kShippedExamples, "ui_text_snapshot.txt"); }

// The other registered examples, compared only where their models are on disk (they are not distributed).
TEST(ui_text_snapshot_local_models) {
    std::vector<std::string> local;
    for (const knit::SceneEntry& e : knit::sceneRegistry()) {
        if (std::find(kShippedExamples.begin(), kShippedExamples.end(), e.id) != kShippedExamples.end()) continue;
        if (!e.isPreset()) continue;
        if (knit::findProjectFile(e.path).empty() || (!e.piecesManifest.empty() && knit::findProjectFile(e.piecesManifest).empty()))
            SKIP("an example's model is not in this checkout (only the wool swatch's model ships with the source)");
        local.push_back(e.id);
    }
    compareSnapshot(local, "ui_text_snapshot_local_models.txt");
}

namespace {

// The Knit step's text once the flat swatch has its verdict (a worker computes it). "" = no verdict in time.
std::string knitTextWithSwatchVerdict(UiFixture& p, float W, float H) {
    p.ui.previews.flat.request(p.ui.s);
    for (int attempt = 0; attempt < 600; ++attempt) {
        const std::string text = snapshotState(p, W, H, "Knit");
        if (text.find("Checking the pattern on a flat swatch") == std::string::npos) return text;
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
    return {};
}

} // namespace

// A flat swatch that does not pass is a quiet note on every scene: every gauge preset fails it and
// still knits (thickness_presets_flat_matrix), so the panel neither warns nor names presets as the
// remedy. An UNRESOLVED swatch (stockinette at 1.12 mm on 5.5 x 3.4) reads as unresolved, not as a failure.
TEST(ui_flat_swatch_is_a_note) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Knit;
    const std::string failed = knitTextWithSwatchVerdict(p, W, H);
    CHECK(failed.find("Flat swatch fails (clearance); the teddy's pieces are checked by their knit.") != std::string::npos);
    CHECK(failed.find("Does not fit a flat swatch") == std::string::npos);
    CHECK(failed.find("Medium / Thin gauge") == std::string::npos);
    CHECK(p.setSetting("strand.diameter", "1.12"));
    const std::string unresolved = knitTextWithSwatchVerdict(p, W, H);
    CHECK(unresolved.find("Flat swatch is unresolved (") != std::string::npos);
    CHECK(unresolved.find("Flat swatch fails") == std::string::npos);
}

// The pattern tiles (2026-11-13): on a chart scene they offer stockinette, reverse stockinette, garter and
// rib 1 x 1, every one selectable (a flat swatch that fails is a note on the tile, never a lock; rib 2 x 2
// is not offered); a click sets the pattern and outdates the knit; each tile's verdict at the current gauge
// is the flat swatch's (thickness_presets_flat_matrix: rib 1 x 1 passes, garter fails); on the teddy's
// spiral pieces every tile but Stockinette is disabled, as the seams.
TEST(ui_knit_patterns) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("wool-swatch");
    p.injectKnit(true);
    p.ui.stage = viewer::Stage::Knit;
    const std::vector<ItemRec> items = settle(p, W, H);
    for (const char* label : {"Stockinette", "Reverse stockinette", "Garter", "Rib 1 x 1"}) {
        const ItemRec* r = findItem(items, label);
        CHECK(r && !r->disabled);
    }
    CHECK(findItem(items, "Rib 2 x 2") == nullptr && findItem(items, "##Pattern") == nullptr);
    CHECK(viewer::knitState(p.ui).word == "Final");
    CHECK(click(p, W, H, "Garter") && p.ui.s.pattern == "garter");
    CHECK(viewer::knitState(p.ui).word == "Outdated");  // knit.pattern is a knitting key
    // The verdicts, once the flat worker is done.
    const std::string text = knitTextWithSwatchVerdict(p, W, H);
    CHECK(!text.empty());
    const std::shared_ptr<const std::vector<knit::PatternReference>> refs = p.ui.previews.flat.get(p.ui.s);
    CHECK(refs != nullptr);
    if (refs) {
        const viewer::FlatSwatchView rib = viewer::flatSwatchOf(p.ui, refs.get(), knit::StitchPattern::Rib1x1);
        const viewer::FlatSwatchView garter = viewer::flatSwatchOf(p.ui, refs.get(), knit::StitchPattern::Garter);
        std::printf("    rib 1 x 1: %s\n    garter: %s\n", rib.text.c_str(), garter.text.c_str());
        CHECK(rib.status == viewer::StageStatus::Ready && garter.status == viewer::StageStatus::Optional && garter.outcome.rfind("fails", 0) == 0);
        CHECK(text.find("Flat swatch " + garter.outcome + "; the parts are checked after knitting.") != std::string::npos);
    }
    UiFixture t("teddy");
    t.ui.stage = viewer::Stage::Knit;
    const std::vector<ItemRec> teddy = settle(t, W, H);
    const ItemRec* stockinette = findItem(teddy, "Stockinette");
    const ItemRec* garterTile = findItem(teddy, "Garter");
    CHECK(stockinette && !stockinette->disabled && garterTile && garterTile->disabled);
    click(t, W, H, "Garter");  // disabled: nothing happens
    CHECK(t.ui.s.pattern == "stockinette");
}

// CameraDrag: a press over the interface is no drag (ImGui keeps the mouse captured for the whole
// press, so overViewport stays false); a Shift-drag's first 4 pixels choose the axis, then only it moves.
TEST(ui_camera_drag_lock) {
    viewer::Camera c;
    c.look({0, 1, 0}, {0, 0, 1});
    viewer::CameraDrag drag;
    auto same = [&](const knit::Vec3& a, const knit::Vec3& b) { return knit::norm(a - b) < 1e-12; };
    viewer::MouseState m;
    m.down[0] = true;  // pressed outside
    CHECK(!drag.update(c, m, 500));
    const knit::Vec3 start = c.forward;
    m.x = 100;
    CHECK(!drag.update(c, m, 500));
    m.x = 200;
    CHECK(!drag.update(c, m, 500));
    CHECK(same(c.forward, start));
    m.down[0] = false;
    drag.update(c, m, 500);
    m.overViewport = true;
    m.down[0] = true;  // pressed inside: the same mouse now drags
    CHECK(!drag.update(c, m, 500));
    m.x += 50;
    CHECK(drag.update(c, m, 500));
    CHECK(!same(c.forward, start));
    m.down[0] = false;
    drag.update(c, m, 500);

    // Shift: 3 pixels are no choice yet, the fourth chooses Yaw, a vertical move then does nothing.
    viewer::Camera locked;
    locked.look({0, 1, 0}, {0, 0, 1});
    m.shift = true;
    m.down[0] = true;
    drag.update(locked, m, 500);
    const knit::Vec3 before = locked.forward;
    m.x += 3;
    CHECK(!drag.update(locked, m, 500));
    CHECK(same(locked.forward, before));
    m.x += 1;
    CHECK(drag.update(locked, m, 500));
    const knit::Vec3 yawed = locked.forward;
    CHECK(!same(yawed, before));
    m.y += 40;
    drag.update(locked, m, 500);
    CHECK(same(locked.forward, yawed));

    // The wheel zooms over the viewport only.
    m = viewer::MouseState{};
    m.wheel = 1;
    const double dist = locked.dist;
    CHECK(!drag.update(locked, m, 500) && locked.dist == dist);
    m.overViewport = true;
    CHECK(drag.update(locked, m, 500) && locked.dist < dist);
}

namespace {

// Every mesh the viewport shows belongs to a scene the renderer keeps, and so do its textures: the
// one lifetime rule (renderer.h). Empty on success, else what is not held.
std::string notHeldContent(const viewer::ViewContent& c, const viewer::MeshOwners& o) {
    auto holds = [](const std::vector<knit::DisplayMesh>& meshes, const knit::DisplayMesh* d) {
        return !meshes.empty() && d >= meshes.data() && d < meshes.data() + meshes.size();
    };
    for (const knit::DisplayMesh* d : c.meshes) {
        bool held = o.preview && holds(o.preview->meshes, d);
        for (const std::shared_ptr<const knit::BuiltScene>& s : o.scenes) held = held || holds(s->meshes, d);
        if (!held) return "mesh of " + c.caption;
    }
    if (!c.meshOwner) return "";
    for (const std::shared_ptr<const knit::BuiltScene>& s : o.scenes)
        if (s.get() == c.meshOwner) return "";
    return "owner of " + c.caption;
}

// Both halves while comparing.
std::string notHeld(const UiFixture& p) {
    const viewer::ViewportFrame f = viewer::viewportFrame(p.ui);
    const viewer::MeshOwners o = viewer::meshOwners(p.ui);
    const std::string shown = notHeldContent(f.shown, o);
    if (!shown.empty() || !f.source) return shown;
    return notHeldContent(*f.source, o);
}

// A hair result on one textured base mesh.
std::shared_ptr<knit::BuiltScene> hairScene(const UiFixture& p) {
    auto hair = std::make_shared<knit::FibreResult>(*p.fakeFibres()->fibreResult);
    hair->workflow = knit::Workflow::Fur;
    auto b = p.fakeFibres();
    b->fibreResult = hair;
    knit::DisplayMesh base;
    base.kind = knit::DisplayMesh::Kind::Base;
    base.pos = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    base.idx = {0, 1, 2};
    base.albedo = 0;
    b->meshes.push_back(base);
    b->albedo.push_back(std::make_shared<knit::DisplayTexture>());
    return b;
}

} // namespace

// The viewer's textures of a hair scene were freed every 120 frames, because its owner was missing
// from the list of live scenes. Every scene a view can show must be listed.
TEST(ui_shown_meshes_cover_view) {
    float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(false);
    for (const viewer::ViewMode mode : {viewer::ViewMode::Mesh, viewer::ViewMode::Knit, viewer::ViewMode::Fibres}) {
        p.ui.view.mode = mode;
        CHECK(notHeld(p).empty());
    }
    p.injectFibres();
    p.ui.view.mode = viewer::ViewMode::Fibres;
    CHECK(notHeld(p).empty());
    p.injectRunning(0, "knitting", 1, 4, true);  // a partial build shows its own scene
    CHECK(viewer::viewContent(p.ui).meshOwner == p.ui.results.partial.get());
    CHECK(notHeld(p).empty());

    UiFixture hair("teddy");
    viewer::chooseWorkflow(hair.ui, knit::Workflow::Fur);  // hair on the bears' mesh
    viewer::chooseWorkflow(hair.ui, knit::Workflow::Fur);
    hair.ui.results.fibres = hairScene(hair);
    hair.ui.stage = viewer::Stage::Fibres;
    hair.ui.view.fibres = viewer::FibreView::Model;
    CHECK(viewer::viewContent(hair.ui).meshOwner == hair.ui.results.fibres.get() && !viewer::viewContent(hair.ui).meshes.empty());
    CHECK(notHeld(hair).empty());
    hair.ui.view.coverage = true;  // the coverage overlay is a scene of its own
    hair.ui.coverage.scene = hairScene(hair);
    CHECK(viewer::viewContent(hair.ui).meshOwner == hair.ui.coverage.scene.get());
    CHECK(notHeld(hair).empty());
    // A parked result of the other workflow is not on screen, so it holds no GPU memory.
    hair.ui.results.parked.fibres = hairScene(hair);
    for (const auto& s : viewer::meshOwners(hair.ui).scenes) CHECK(s != hair.ui.results.parked.fibres);

    UiFixture reference("teddy");  // the reference sample has no scene, only its preview
    reference.ui.stage = viewer::Stage::Fibres;
    auto sample = std::make_shared<knit::FibrePreview>();
    knit::DisplayMesh patch;
    patch.pos = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    patch.idx = {0, 1, 2};
    sample->meshes.push_back(patch);
    reference.ui.previews.shown = sample;
    CHECK(!viewer::viewContent(reference.ui).meshes.empty());
    CHECK(notHeld(reference).empty());
}

// Every flag of the command line lands in its field (parseViewerOptions, viewer_script.cpp).
TEST(ui_viewer_options) {
    std::vector<std::string> words = {"knit_viewer", "--scene", "teddy", "--shot", "a.png", "--view", "top", "--settings", "s.ks", "--import", "f.obj",
                                      "--units", "cm", "--up", "z", "--plain-nodes", "bow|eyes",
                                      "--stage", "export", "--do", "knit,final", "--size", "800x600", "--ui-scale", "1.5", "--zoom", "4", "--plain", "off",
                                      "--lighting", "flat", "--occlusion", "0.5", "--bench", "20", "--target-frac", "0.1,0.2,0.3"};
    std::vector<char*> argv;
    for (std::string& w : words) argv.push_back(w.data());
    const viewer::ViewerOptions o = viewer::parseViewerOptions(static_cast<int>(argv.size()), argv.data());
    CHECK(o.scene == "teddy" && o.shotPath == "a.png" && o.startView == viewer::CameraView::Top && o.settingsFile == "s.ks");
    CHECK(o.importFile == "f.obj" && o.units == "cm" && o.up == "z");
    CHECK((o.importChoices.plainNodes == std::set<std::string>{"bow", "eyes"}));
    CHECK(o.stage == "export" && o.steps == "knit,final" && o.width == 800 && o.height == 600 && o.uiScale == 1.5f);
    CHECK(o.zoom == 4 && o.plain == 0 && o.lighting == "flat" && o.occlusion == 0.5 && o.benchFrames == 20);
    CHECK(o.targetFraction[0] == 0.1 && o.targetFraction[1] == 0.2 && o.targetFraction[2] == 0.3);
    CHECK(viewer::retiredOption(static_cast<int>(argv.size()), argv.data()).empty());
    std::string lab = "--lab";
    char* retired[] = {argv[0], lab.data()};
    CHECK(viewer::retiredOption(2, retired) == "--lab");
    CHECK(viewer::parseUnits("mm") == 1 && viewer::parseUnits("cm") == 10 && viewer::parseUnits("m") == 1000 && viewer::parseUnits("inch") == 25.4);
    CHECK(viewer::parseUnits("2.5") == 2.5);
    CHECK(viewer::stageNamed("import") == viewer::Stage::Import && viewer::stageNamed("anything") == viewer::Stage::Knit);
    CHECK(viewer::cameraViewNamed("back") == viewer::CameraView::Back && !viewer::cameraViewNamed("sideways"));
    CHECK(viewer::cameraViewNamed("threequarter") == viewer::CameraView::ThreeQuarter);
}

// Every --do step (applyScriptStep, viewer_script.cpp) does what it says, the older names included.
TEST(ui_script_steps) {
    float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    viewer::applyScriptStep(p.ui, "final");
    CHECK(p.finalProfile());
    viewer::applyScriptStep(p.ui, "draft");
    CHECK(knit::isPreviewProfile(p.ui.s));
    viewer::applyScriptStep(p.ui, "mesh");
    CHECK(p.ui.view.mode == viewer::ViewMode::Mesh);
    viewer::applyScriptStep(p.ui, "knitview");
    CHECK(p.ui.view.mode == viewer::ViewMode::Knit);
    viewer::applyScriptStep(p.ui, "fibreview");
    CHECK(p.ui.view.mode == viewer::ViewMode::Fibres);
    viewer::applyScriptStep(p.ui, "model");
    CHECK(p.ui.view.fibres == viewer::FibreView::Model);
    viewer::applyScriptStep(p.ui, "reference");
    CHECK(p.ui.view.fibres == viewer::FibreView::Reference);
    viewer::applyScriptStep(p.ui, "coverage");
    CHECK(p.ui.view.coverage);
    viewer::applyScriptStep(p.ui, "compare");
    CHECK(p.ui.view.compare);
    viewer::applyScriptStep(p.ui, "wire");
    CHECK(p.ui.view.wire);
    viewer::applyScriptStep(p.ui, "populations");
    CHECK(p.ui.view.fibrePopulations);
    viewer::applyScriptStep(p.ui, "about:3");
    CHECK(p.ui.panels.aboutOpen && p.ui.panels.aboutTabRequest == 3);
    for (const char* name : {"stage:import", "stage:fibres", "stage:export", "stage:knit"}) {
        viewer::applyScriptStep(p.ui, name);
        CHECK(p.ui.stage == viewer::stageNamed(std::string(name).substr(6)));
    }
    viewer::applyScriptStep(p.ui, "choose-fur");
    CHECK(viewer::furShown(p.ui));
    viewer::applyScriptStep(p.ui, "choose-knit");
    CHECK(knit::knitWorkflow(p.ui.s) && !viewer::furShown(p.ui));
    p.injectKnit(false);
    viewer::applyScriptStep(p.ui, "knit");
    CHECK(p.ui.job.active && p.ui.job.kind == knit::BuiltScene::Kind::Knit);
    viewer::applyScriptStep(p.ui, "cancel");
    CHECK(p.ui.job.cancelRequested);
    p.ui.job.active = false;
    viewer::applyScriptStep(p.ui, "fibres");
    CHECK(p.ui.s.fibres && p.ui.job.active && p.ui.job.kind == knit::BuiltScene::Kind::Fibres);
    p.ui.job.active = false;
    viewer::applyScriptStep(p.ui, "obj");
    CHECK(p.ui.exportTo.format == "obj" && p.ui.job.active && p.ui.job.kind == knit::BuiltScene::Kind::Export);
    p.ui.job.active = false;
    viewer::applyScriptStep(p.ui, "export");
    CHECK(p.ui.exportTo.format == "ppt");
}

// The SSAO kernel (mesh_lighting.h): samples in the unit hemisphere around +z, off the tangent
// plane, near samples first (the scale grows with the index), the same for a seed.
TEST(ui_ssao_kernel) {
    const std::vector<knit::Vec3> k = viewer::ssaoKernel(16, 7);
    CHECK(k.size() == 16);
    double firstHalf = 0, secondHalf = 0;
    for (size_t i = 0; i < k.size(); ++i) {
        const double len = knit::norm(k[i]);
        CHECK(len <= 1.0 + 1e-12);
        CHECK(k[i].z >= 0.15 * len - 1e-12);
        (i < k.size() / 2 ? firstHalf : secondHalf) += len;
    }
    CHECK(firstHalf < 0.5 * secondHalf);
    const std::vector<knit::Vec3> again = viewer::ssaoKernel(16, 7);
    for (size_t i = 0; i < k.size(); ++i) CHECK(knit::norm(k[i] - again[i]) == 0);
    const std::vector<knit::Vec3> r = viewer::ssaoRotations(11);
    CHECK(r.size() == 16);
    for (const knit::Vec3& v : r) CHECK(std::abs(knit::norm(v) - 1) < 1e-12 && v.z == 0);
}

// The turntable camera (camera.cpp): any drag keeps the view unrolled (the right axis stays
// horizontal) and the pitch within straight down / straight up; Shift's locks move one axis
// only; Upright levels the view facing the same way.
TEST(ui_camera_turntable) {
    viewer::Camera c;
    c.look({0, 1, 0}, {0, 0, 1});
    const double drags[][2] = {{40, 15}, {-120, 300}, {75, -900}, {13, 7}, {-2000, 33}, {5, 5000}};
    for (const auto& d : drags) {
        c.orbit(d[0], d[1]);
        CHECK(std::abs(c.right().z) < 1e-9);
        CHECK(c.up.z >= -1e-9);
        CHECK(std::abs(knit::dot(c.forward, c.up)) < 1e-9);
    }
    c.look({0, 1, 0}, {0, 0, 1});
    c.orbit(0, 100);
    const double pitch = std::asin(c.forward.z);
    c.orbit(250, 80, viewer::OrbitLock::Yaw);  // turn only
    CHECK(std::abs(std::asin(c.forward.z) - pitch) < 1e-9);
    const knit::Vec3 h = c.heading();
    c.orbit(250, 80, viewer::OrbitLock::Pitch);  // tilt only
    CHECK(knit::norm(c.heading() - h) < 1e-9);
    c.level();
    CHECK(std::abs(c.forward.z) < 1e-12 && c.up.z == 1 && knit::norm(c.heading() - h) < 1e-9);
    c.look({0, 0, -1}, {0, 1, 0});  // Top: straight down, +y at the top of the screen
    c.orbit(0, 0);
    CHECK(knit::norm(c.forward - knit::Vec3{0, 0, -1}) < 1e-9 && knit::norm(c.up - knit::Vec3{0, 1, 0}) < 1e-9);
    c.orbit(0, -50);  // tilting up from Top looks along +y
    CHECK(c.forward.y > 0 && std::abs(c.right().z) < 1e-9);
}

// The fibres on screen follow the colour sources as the export writes them: the texture their own
// colour, "object" their part's. One yarn colour keeps the populations told apart, as does a stray or
// oriented fibre whose match with the yarn is off.
TEST(ui_fibre_colour_sources) {
    knit::DisplayLines lines;
    lines.population = {0, 1, 2};
    lines.partColor.assign(9, 0.2f);
    lines.ownColor.assign(9, 0.8f);
    lines.pos.assign(9, 0);
    const float part = viewer::toSrgb(0.2f), own = viewer::toSrgb(0.8f);
    auto shows = [&](const knit::SceneSettings& s, size_t v, float colour) -> bool {
        return std::fabs(viewer::fibreVertexColour(lines, v, viewer::fibreColoursOf(s)).x - colour) < 1e-5;
    };
    auto isFixed = [&](const knit::SceneSettings& s, size_t v, const knit::Vec3& colour) -> bool {
        const knit::Vec3 c = viewer::fibreVertexColour(lines, v, viewer::fibreColoursOf(s));
        return c.x == colour.x && c.y == colour.y && c.z == colour.z;
    };
    knit::SceneSettings s;
    s.pileMaterial.color = {0.1, 0.2, 0.3};
    s.yarnColourSource = "texture";
    s.pileColourSource = "texture";
    CHECK(shows(s, 0, own) && shows(s, 1, own) && shows(s, 2, own));
    s.yarnColourSource = "object";
    s.pileColourSource = "object";
    CHECK(shows(s, 0, part) && shows(s, 1, part) && shows(s, 2, part));
    s.yarnColourSource = "single";
    s.pileColourSource = "single";
    CHECK(isFixed(s, 0, knit::kOrientedFibreColour) && isFixed(s, 1, knit::kStrayFibreColour) && isFixed(s, 2, s.pileMaterial.color));
    s.yarnColourSource = "texture";
    s.fibreRandomInheritsYarn = false;
    CHECK(shows(s, 0, own) && isFixed(s, 1, knit::kStrayFibreColour));
    s.fibreInheritsYarn = false;
    CHECK(isFixed(s, 0, knit::kOrientedFibreColour));
    // The display switch: wool yellow and blue whatever the source, hair and the settings as they were.
    for (const char* source : {"object", "texture"}) {
        knit::SceneSettings m;
        m.yarnColourSource = m.pileColourSource = source;
        const knit::SceneSettings before = m;
        const viewer::FibreColours exported = viewer::fibreColoursOf(m), told = viewer::populationColours(exported);
        const auto colourOf = [&](size_t v, const viewer::FibreColours& c) { return viewer::fibreVertexColour(lines, v, c); };
        const knit::Vec3 oriented = colourOf(0, told), stray = colourOf(1, told);
        CHECK(oriented.x == knit::kOrientedFibreColour.x && oriented.y == knit::kOrientedFibreColour.y && oriented.z == knit::kOrientedFibreColour.z);
        CHECK(stray.x == knit::kStrayFibreColour.x && stray.y == knit::kStrayFibreColour.y && stray.z == knit::kStrayFibreColour.z);
        CHECK(colourOf(2, told).x == colourOf(2, exported).x);
        CHECK(knit::writeSettings(m) == knit::writeSettings(before));
    }
    // Lines without their own colours: the part's.
    lines.ownColor.clear();
    s = knit::SceneSettings{};
    s.pileColourSource = "texture";
    CHECK(shows(s, 2, part));
}

// "Show fibre populations" is offered for wool on the model only: the reference patch tells them
// apart anyway, the Knit view has no fibres and hair is one population.
TEST(ui_fibre_populations_checkbox) {
    const float W = 1600, H = 1200;
    UiContext ctx(W, H);
    auto drawn = [&](UiFixture& p) { return popoverLabels(p, W, H).count("Fibre populations") == 1; };
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Fibres;
    p.ui.view.mode = viewer::ViewMode::Fibres;
    p.ui.view.fibres = viewer::FibreView::Model;
    CHECK(drawn(p));
    p.ui.view.fibres = viewer::FibreView::Reference;
    CHECK(!drawn(p));
    p.ui.view.fibres = viewer::FibreView::Model;
    p.ui.view.mode = viewer::ViewMode::Knit;
    CHECK(!drawn(p));
    p.ui.view.mode = viewer::ViewMode::Fibres;
    CHECK(drawn(p));
    CHECK(p.setSetting("workflow", "fur"));
    CHECK(!drawn(p));
}

// Fibre display chunks (arrangeFibreLines, fibre_lines.h): fur on a sphere of radius 20 mm, each
// hair three 1 mm segments along its normal. The chunks cover every vertex in whole segments,
// each a compact patch; the vertices are only reordered. A close-up camera skips chunks, and no
// vertex of a skipped chunk projects inside the view; turned away, it skips every chunk.
// The import's rules hold where the choice is made, with no frame drawn: an analysed file with Fur chosen keeps the
// file's mesh; "Change the import" keeps the object's workflow (Fur: the file's mesh) and rebuilds an import that kept
// the file's mesh under Knit; a height of 0 gives no unit, so Import stays disabled.
TEST(ui_import_rules_without_a_frame) {
    using O = knit::ImportOptions;
    UiFixture fur("teddy");
    CHECK(fur.openImport(writeBoxObj("rules_skip")));
    viewer::chooseWorkflow(fur.ui, knit::Workflow::Fur);
    CHECK(fur.ui.import.options.surface == O::Surface::KeepFile);
    viewer::closeImport(fur.ui);
    UiFixture change("teddy");
    change.ui.s.importSource = writeBoxObj("rules_change");
    change.ui.s.importSurface = "wrapped";
    change.ui.s.workflow = "fur";
    viewer::changeImport(change.ui);
    CHECK(change.ui.workflowChosen && change.ui.import.options.surface == O::Surface::KeepFile);
    viewer::closeImport(change.ui);
    // A draft of the old direct-mesh option (import.surface = file) is rebuilt under Knit.
    UiFixture rebuild("teddy");
    rebuild.ui.s.importSource = writeBoxObj("rules_rebuild");
    rebuild.ui.s.importSurface = "file";
    viewer::changeImport(rebuild.ui);
    CHECK(rebuild.ui.import.options.surface == O::Surface::Wrapped);
    viewer::closeImport(rebuild.ui);
    UiFixture height("teddy");
    CHECK(height.openImport(writeTallBoxObj("rules_height")));
    viewer::chooseSizeByHeight(height.ui, true);
    viewer::setImportHeight(height.ui, 400);
    CHECK(std::abs(height.ui.import.options.unitMm - 5) < 1e-12);
    viewer::setImportHeight(height.ui, 0);
    CHECK(height.ui.import.options.unitMm == 0);
    viewer::setCustomImportUnit(height.ui, -1);  // not the size while it is set by height
    CHECK(height.ui.import.customUnitMm == -1 && height.ui.import.options.unitMm == 0);
    viewer::closeImport(height.ui);
}

// A scripted import follows the panel's rule: Fur keeps the file's mesh (and its plain components: plain in both
// workflows), and a knitted import is rebuilt.
TEST(ui_import_scripted_skip_keeps_charts) {
    using O = knit::ImportOptions;
    UiFixture p("teddy");
    viewer::openMeshFile(p.ui, writeBoxObj("scripted_skip"));
    viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);  // as the command line's settings choose it
    viewer::presetImportUnits(p.ui, 1, knit::UpAxis::Y);  // before the first poll, as --units
    for (int i = 0; i < 2000 && p.ui.import.session.snapshot().phase != knit::ImportSession::Phase::Analysed; ++i) {
        viewer::pollSceneUi(p.ui);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    viewer::ImportChoices furChoices;
    furChoices.plainNodes = {"Part 1"};  // --plain-nodes: a plain component under Fur too
    viewer::confirmImport(p.ui, furChoices);
    CHECK(p.ui.import.options.surface == O::Surface::KeepFile && viewer::importRequest(p.ui.import).plainNodes == std::set<std::string>{"Part 1"});
    p.ui.import.session.cancel();
    viewer::closeImport(p.ui);
    viewer::openMeshFile(p.ui, writeBoxObj("scripted_knit"));
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
    viewer::presetImportUnits(p.ui, 1, knit::UpAxis::Y);
    for (int i = 0; i < 2000 && p.ui.import.session.snapshot().phase != knit::ImportSession::Phase::Analysed; ++i) {
        viewer::pollSceneUi(p.ui);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    p.ui.import.options.surface = O::Surface::KeepFile;  // whatever stands before: the confirmation decides
    viewer::confirmImport(p.ui, viewer::ImportChoices{});
    CHECK(p.ui.import.options.surface == O::Surface::Wrapped);
    p.ui.import.session.cancel();
    viewer::closeImport(p.ui);
}

// "Knit final" in the Knit footer on a Draft result: the Final profile, then the knit.
TEST(ui_knit_final_button) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.injectKnit(false);  // a Draft result
    p.ui.stage = viewer::Stage::Knit;
    CHECK(!knit::isFinalProfile(p.ui.s, *viewer::entryOf(p.ui)));
    CHECK(click(p, W, H, "Knit final"));
    CHECK(knit::isFinalProfile(p.ui.s, *viewer::entryOf(p.ui)));
    CHECK(p.ui.job.active && p.ui.job.kind == knit::BuiltScene::Kind::Knit && !p.ui.job.meshOnly);
}

// Each footer action's notes and words are about the job it submits (its profile applied to a copy),
// not about the profile the settings hold: with Draft settings on a machine without the knitting
// solver, Knit final carries the solver note (computeNotes on the settings themselves says nothing),
// Knit draft does not; the Final tooltip names the fallback the Final profile installs, also after the
// gauge became Custom while an older preset's fallback is still in the settings.
TEST(ui_knit_action_notes_per_profile) {
    UiFixture p("teddy");
    const knit::SceneEntry& e = *viewer::entryOf(p.ui);
    knit::ComputeCaps noSolver;
    noSolver.device = true;
    noSolver.ncp = false;
    noSolver.whyNoNcp = "test: no solver";
    CHECK(knit::isPreviewProfile(p.ui.s) && knit::computeNotes(p.ui.s, noSolver).empty());
    viewer::KnitActionView draft = viewer::knitActionView(p.ui, viewer::KnitQuality::Draft, noSolver);
    viewer::KnitActionView final = viewer::knitActionView(p.ui, viewer::KnitQuality::Final, noSolver);
    std::printf("    Draft settings: draft notes %zu, final notes %zu (%s)\n", draft.notes.size(), final.notes.size(),
                final.notes.empty() ? "-" : final.notes[0].c_str());
    CHECK(draft.notes.empty() && final.notes.size() == 1 && final.notes[0].find("Final knit cannot run") != std::string::npos);
    auto listed = [](const std::vector<double>& fallback) {
        std::string t;
        for (double d : fallback) {
            char b[16];
            std::snprintf(b, sizeof b, " %.2f", d);
            t += b;
        }
        return "Diameter fallback:" + t + " mm";
    };
    knit::SceneSettings full = p.ui.s;
    knit::applyFullProfile(full, e);
    CHECK(!full.diameterFallback.empty() && p.ui.s.diameterFallback.empty());
    CHECK(final.tooltip.find(listed(full.diameterFallback)) != std::string::npos && draft.tooltip.find("Diameter fallback") == std::string::npos);
    // Final with the Medium preset's fallback (the teddy's Dense one equals the scene default), then a
    // Custom diameter: the settings keep the preset's list, the Final copy takes the scene default's.
    knit::applyGaugeChoice(p.ui.s, e, "medium");
    knit::applyFullProfile(p.ui.s, e);
    const std::vector<double> stale = p.ui.s.diameterFallback;
    p.ui.s.diameter = 1.40;
    p.ui.s.thicknessPreset = "custom";
    full = p.ui.s;
    knit::applyFullProfile(full, e);
    final = viewer::knitActionView(p.ui, viewer::KnitQuality::Final, noSolver);
    std::printf("    stale fallback%s, Final installs%s\n", listed(stale).c_str() + 18, listed(full.diameterFallback).c_str() + 18);
    CHECK(full.diameterFallback != stale);
    CHECK(final.tooltip.find(listed(full.diameterFallback)) != std::string::npos && final.tooltip.find(listed(stale)) == std::string::npos);
}

// The Custom stitch fields hold every stitch the Gauge list derives (both ends of the diameter and
// clearance fields, Snug and Loose): switching to Custom and dragging moves the value by the drag,
// never jumps to a clamp.
TEST(ui_knit_custom_range_ends) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    struct Corner {
        double diameter, clearance;
        const char* gauge;
        float dx;  // towards the inside of the range
    };
    const Corner corners[] = {{viewer::kYarnDiameterRange.lo, viewer::kClearanceRange.lo, "yarn-tight", 10},
                              {viewer::kYarnDiameterRange.hi, viewer::kClearanceRange.hi, "yarn-loose", -10}};
    for (const Corner& c : corners) {
        UiFixture p("teddy");
        p.ui.stage = viewer::Stage::Knit;
        const knit::SceneEntry& e = *viewer::entryOf(p.ui);
        p.ui.s.diameter = c.diameter;
        p.ui.s.clearance = c.clearance;
        knit::applyGaugeChoice(p.ui.s, e, c.gauge);
        const double w = p.ui.s.wale, h = p.ui.s.course;
        const bool inside = w >= viewer::kStitchWidthRange.lo && w <= viewer::kStitchWidthRange.hi && h >= viewer::kRowHeightRange.lo &&
                            h <= viewer::kRowHeightRange.hi;
        CHECK(click(p, W, H, "##Gauge") && pickItem(p, W, H, "Custom"));
        CHECK(p.ui.s.stitchSize == "custom" && p.ui.s.wale == w && p.ui.s.course == h);  // Custom keeps the derived stitch
        CHECK(drag(p, W, H, "##Stitch width", c.dx));
        const double moved = p.ui.s.wale - w;
        std::printf("    D %.2f, c %.2f, %s: %.2f x %.2f mm (inside %d); first drag %+.2f mm\n", c.diameter, c.clearance, c.gauge, w, h, inside, moved);
        CHECK(inside && moved != 0 && std::fabs(moved) < 1.0 && (moved > 0) == (c.dx > 0));
    }
}

namespace {

// Polls the import until its file is analysed (stage 2) or failed.
bool waitForAnalysis(UiFixture& p) {
    for (int i = 0; i < 2000; ++i) {
        viewer::pollSceneUi(p.ui);
        const knit::ImportSession::Phase phase = p.ui.import.session.snapshot().phase;
        if (phase == knit::ImportSession::Phase::Analysed) return true;
        if (phase == knit::ImportSession::Phase::Failed) return false;
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    return false;
}

} // namespace

// The two fit checks are independent and shown together (a surface prepared for 5.5 mm under a Regular
// 1.86 mm stitch warns, the flat swatch still reports its own verdict), and Review Import is a round
// trip: the same file is prepared again for the current stitch width, and back in the Knit step every
// Knit choice is as it was (row direction, One colour, the character, the gauge) while the knit shows
// as outdated. Before the fix the completed import reset the row direction and the colour source.
TEST(ui_knit_review_import_round_trip) {
    namespace fs = std::filesystem;
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_review";
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    UiFixture p("teddy");
    const bool opened = p.openImport(writeBoxObj("review"));
    p.ui.import.options.unitMm = 1.0;  // a 90 mm box, as the colourwork test's
    const bool imported = opened && p.finishImport();
    CHECK(imported && p.ui.s.scene == "import" && p.ui.s.importSurface == "wrapped");
    if (imported) {
        const knit::SceneEntry& e = *viewer::entryOf(p.ui);
        const double prepared = p.ui.s.importWrapStitchMm;
        p.ui.stage = viewer::Stage::Knit;
        p.ui.s.diameter = 1.86;
        knit::applyGaugeChoice(p.ui.s, e, "yarn-standard");
        p.ui.s.angleDeg = 30;
        p.ui.s.yarnColourSource = "single";
        knit::applyPreset(knit::PresetGroup::HandKnit, p.ui.s, *knit::presetWithId(knit::PresetGroup::HandKnit, "rustic"));
        p.injectKnit(false);  // a knit of these settings: the new surface must outdate it
        const knit::SceneSettings before = p.ui.s;
        const viewer::PreparedSurface check = viewer::preparedSurfaceCheck(p.ui.s);
        std::printf("    prepared for %.2f mm, stitch %.2f mm: too wide %d\n", check.preparedMm, p.ui.s.wale, check.tooWide);
        CHECK(check.applies && check.tooWide && check.preparedMm == prepared);
        // Both checks on screen together: the warning with its action, and the swatch's own line.
        std::string text;
        for (int i = 0; i < 600; ++i) {
            text = renderedText(p, W, H);
            if (text.find("Checking the pattern on a flat swatch") == std::string::npos) break;
            viewer::pollSceneUi(p.ui);
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        const size_t swatch = text.find("Flat swatch");
        std::printf("    flat swatch: %s\n", swatch == std::string::npos ? "(no line)" : text.substr(swatch, text.find('\n', swatch) - swatch).c_str());
        CHECK(text.find("wider and may tangle") != std::string::npos && swatch != std::string::npos);
        CHECK(text.find("Flat swatch fits at these values.") != std::string::npos);  // the fixture: passes flat, warns on the surface
        CHECK(click(p, W, H, "Review Import"));
        CHECK(p.ui.stage == viewer::Stage::Import && p.ui.import.keepKnitChoices && p.ui.import.options.wrapStitchMm == before.wale);
        const bool reviewed = waitForAnalysis(p) && p.finishImport();
        CHECK(reviewed);
        std::printf("    after the review: stage %d, angle %.0f, colour %s, hand %s, gauge %s, prepared %.2f, knit \"%s\"\n",
                    static_cast<int>(p.ui.stage), p.ui.s.angleDeg, p.ui.s.yarnColourSource.c_str(),
                    knit::matchingFilePreset(knit::PresetGroup::HandKnit, p.ui.s) ? knit::matchingFilePreset(knit::PresetGroup::HandKnit, p.ui.s)->id.c_str() : "custom",
                    knit::currentGaugeChoice(p.ui.s, e).c_str(), p.ui.s.importWrapStitchMm, viewer::knitState(p.ui).word.c_str());
        CHECK(p.ui.stage == viewer::Stage::Knit && !p.ui.import.keepKnitChoices);
        CHECK(p.ui.s.angleDeg == 30 && p.ui.s.yarnColourSource == "single");
        CHECK(knit::matchingFilePreset(knit::PresetGroup::HandKnit, p.ui.s) && knit::matchingFilePreset(knit::PresetGroup::HandKnit, p.ui.s)->id == "rustic");
        CHECK(knit::currentGaugeChoice(p.ui.s, e) == "yarn-standard" && p.ui.s.wale == before.wale && p.ui.s.diameter == before.diameter);
        CHECK(p.ui.s.importWrapStitchMm == before.wale && !viewer::preparedSurfaceCheck(p.ui.s).tooWide);
        CHECK(viewer::knitState(p.ui).word == "Outdated");
    }
    _putenv_s("KNIT_IMPORT_CACHE", "");
}

// The collapsed sections' summaries say the settings they hold.
TEST(ui_knit_detail_summaries) {
    UiFixture p("teddy");
    const knit::SceneEntry& teddy = *viewer::entryOf(p.ui);
    knit::SceneSettings s = p.ui.s;
    s.plies = 1;
    s.depthAmplitude = 0.8;
    s.clearance = 0.1;
    s.flatten = false;
    CHECK(viewer::yarnGeometrySummary(s) == std::string("0.80 D") + viewer::kDot + "0.10 mm clearance");  // the construction has its own section
    s.flatten = true;
    CHECK(viewer::yarnGeometrySummary(s) == std::string("0.80 D") + viewer::kDot + "0.10 mm clearance" + viewer::kDot + "contact flattening");
    s.flatten = false;
    s.hand = knit::HandKnit{};
    knit::applyPreset(knit::PresetGroup::HandKnit, s, *knit::presetWithId(knit::PresetGroup::HandKnit, "off"));
    CHECK(viewer::characterSummary(s) == "Every stitch alike");
    knit::applyPreset(knit::PresetGroup::HandKnit, s, *knit::presetWithId(knit::PresetGroup::HandKnit, "handmade"));
    CHECK(viewer::characterSummary(s).rfind("20 % uneven", 0) == 0);
    CHECK(viewer::stitchLayoutSummary(s, teddy) == "Sewn pieces");  // the teddy's spiral pieces; the pattern has its own section
    UiFixture q("wool-swatch");
    const knit::SceneEntry& chart = *viewer::entryOf(q.ui);
    knit::SceneSettings c = q.ui.s;
    c.workflow = "knit";
    c.pattern = "rib1x1";
    c.coverage = "overlap";
    c.angleDeg = 30;
    std::printf("    %s: %s\n", chart.id.c_str(), viewer::stitchLayoutSummary(c, chart).c_str());
    if (!chart.knitsAsPieces())
        CHECK(viewer::stitchLayoutSummary(c, chart) == std::string("Overlapped seams") + viewer::kDot + "30 deg");
    const knit::FilePreset* level = knit::matchingFilePreset(knit::PresetGroup::StrandMesh, p.ui.s);
    CHECK(level && viewer::meshDetailSummary(p.ui.s, false) == level->label);
    knit::SceneSettings m = p.ui.s;
    m.tubeSides = 5;
    m.tubeRing = 0.4;
    CHECK(viewer::meshDetailSummary(m, false) == "Custom (5 sides, 0.40 mm)");
}

// An export level whose tessellation differs from the shown tubes rebuilds them at once (no knitting), only
// while the shown knit is current.
TEST(ui_tubes_follow_level) {
    UiFixture p("teddy");
    p.injectKnit(true);
    CHECK(p.chooseExportDetail("Light"));
    CHECK(p.ui.job.active && p.ui.job.meshOnly);
    for (int i = 0; i < 2000 && p.ui.job.active; ++i) {
        viewer::pollSceneUi(p.ui);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    CHECK(!p.ui.job.active);
    p.ui.s.wale += 0.5;  // the knit is outdated now
    CHECK(!viewer::knitCurrent(p.ui));
    CHECK(p.chooseExportDetail("Full"));
    CHECK(!p.ui.job.active);
}

namespace {

// The stitch width and row height are what the yarn gives (applyStitchSize changes nothing).
bool followsTheYarn(const knit::SceneSettings& s) {
    knit::SceneSettings derived = s;
    knit::applyStitchSize(derived);
    return derived.wale == s.wale && derived.course == s.course;
}

} // namespace

// In yarn mode the stitch size follows its inputs where they change (Q3, 2026-11-06), not once per frame of
// the Knit step: a scene switch that sets the pattern, a clearance drag; a frame with no edit writes nothing.
TEST(ui_stitch_size_on_change) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    p.ui.s.stitchSize = "yarn";
    p.ui.s.pattern = knit::toString(knit::StitchPattern::Garter);
    knit::applyStitchSize(p.ui.s);
    viewer::openExample(p.ui, "teddy");  // spiral pieces: stockinette (the pattern leaves w and h as they are)
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
    CHECK(p.ui.s.stitchSize == "yarn" && p.ui.s.pattern == knit::toString(knit::StitchPattern::Stockinette));
    CHECK(followsTheYarn(p.ui.s));
    p.ui.stage = viewer::Stage::Knit;
    openEverything(p, W, H);
    const double before = p.ui.s.wale;
    CHECK(drag(p, W, H, "##Clearance", 60));
    CHECK(followsTheYarn(p.ui.s) && p.ui.s.wale > before);
    p.ui.s.wale = 9.99;  // not what the yarn gives: only an edit of an input brings it back
    settle(p, W, H);
    CHECK(p.ui.s.wale == 9.99);
}

// The buttons that need the window ask for it (PlatformRequest); main.cpp answers after the frame. Without dialogs
// or a file manager they ask nothing, and the controls say so instead.
TEST(ui_platform_requests) {
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    using Rq = viewer::PlatformRequest;
    UiFixture p("teddy");
    p.ui.stage = viewer::Stage::Export;
    CHECK(click(p, W, H, "Choose folder ..."));
    CHECK(p.ui.platform.request == Rq::None);
    p.ui.platform.dialogs = true;
    const std::string chosen = (std::filesystem::temp_directory_path() / "knit_ui_requests").make_preferred().string();
    CHECK(click(p, W, H, "Choose folder ...") && p.ui.platform.request == Rq::ChooseExportFolder);
    CHECK(p.answerFolderDialog(chosen) && p.ui.exportTo.dir == chosen && p.ui.exportTo.dirChosen);
    p.ui.stage = viewer::Stage::Import;
    CHECK(click(p, W, H, "Open mesh (OBJ, PLY) ...") && p.ui.platform.request == Rq::OpenMesh);
    p.ui.platform.request = Rq::None;
    p.ui.stage = viewer::Stage::Export;
    p.injectExportDone(chosen);
    CHECK(renderedText(p, W, H).find("Open folder") == std::string::npos);  // no file manager: no button
    p.ui.platform.fileManager = true;
    CHECK(click(p, W, H, "Open folder") && p.ui.platform.request == Rq::ShowExportFolder);
    UiFixture hair("teddy");
    viewer::chooseWorkflow(hair.ui, knit::Workflow::Fur);  // hair on the bears' mesh
    hair.ui.stage = viewer::Stage::Fibres;
    hair.ui.s.pile.coverage = "mask";
    CHECK(click(hair, W, H, "Load ...") && hair.ui.platform.request == Rq::None);
    hair.ui.platform.dialogs = true;
    CHECK(click(hair, W, H, "Load ...") && hair.ui.platform.request == Rq::OpenMaskImage);
}

// The Fibres step's mask image check (checkMaskImage): the file's size and write time are read on every call,
// the file is hashed again only when they change, and Reload accepts the file as it is.
TEST(ui_fibre_mask_image_changed) {
    namespace fs = std::filesystem;
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    const fs::path d = fs::temp_directory_path() / "knit_ui_mask";
    std::error_code ec;
    fs::remove_all(d, ec);
    fs::create_directories(d);
    const fs::path mask = d / "mask.png", other = d / "other.png";
    std::ofstream(mask, std::ios::binary) << "first image";
    std::ofstream(other, std::ios::binary) << "the other, longer image";
    UiFixture p("teddy");
    viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);  // hair on the bears' mesh
    viewer::loadMaskImage(p.ui, mask.string());
    CHECK(!p.ui.s.pile.maskHash.empty() && !viewer::checkMaskImage(p.ui) && p.ui.panels.maskOnDisk.path == mask.string());
    // Edited: the same size, other bytes, a later write time.
    const fs::file_time_type loaded = fs::last_write_time(mask);
    std::ofstream(mask, std::ios::binary) << "frist image";
    fs::last_write_time(mask, loaded + std::chrono::seconds(2));
    CHECK(viewer::checkMaskImage(p.ui));
    p.ui.stage = viewer::Stage::Fibres;
    p.ui.s.pile.coverage = "mask";
    CHECK(renderedText(p, W, H).find("The image changed since it was loaded: Reload it.") != std::string::npos);
    viewer::reloadMaskImage(p.ui);
    CHECK(!viewer::checkMaskImage(p.ui));
    CHECK(renderedText(p, W, H).find("The image changed since it was loaded") == std::string::npos);
    // An edit that keeps the size and the write time is not seen: the check reads metadata, and hashes on a change.
    const fs::file_time_type edited = fs::last_write_time(mask);
    std::ofstream(mask, std::ios::binary) << "FRIST image";
    fs::last_write_time(mask, edited);
    CHECK(!viewer::checkMaskImage(p.ui));
    viewer::reloadMaskImage(p.ui);
    // Replaced by another file (renamed over it).
    fs::copy_file(other, d / "replacement.png");
    fs::rename(d / "replacement.png", mask);
    CHECK(viewer::checkMaskImage(p.ui));
    viewer::reloadMaskImage(p.ui);
    CHECK(!viewer::checkMaskImage(p.ui));
    // Deleted, then unreadable (a folder in its place).
    fs::remove(mask);
    CHECK(viewer::checkMaskImage(p.ui));
    fs::create_directory(mask);
    CHECK(viewer::checkMaskImage(p.ui));
    fs::remove(mask);
    // Another file: loaded with its own hash, and checked against it.
    viewer::loadMaskImage(p.ui, other.string());
    CHECK(!viewer::checkMaskImage(p.ui) && p.ui.panels.maskOnDisk.path == other.string());
    // No mask: nothing to check.
    p.ui.s.pile.maskPath.clear();
    p.ui.s.pile.maskHash.clear();
    CHECK(!viewer::checkMaskImage(p.ui));
}

TEST(ui_fibre_chunks) {
    knit::DisplayLines lines;
    const int hairs = 60000;
    knit::Vec3 sum;
    for (int i = 0; i < hairs; ++i) {
        // Even points on the sphere (Fibonacci spiral).
        const double z = 1 - 2 * (i + 0.5) / hairs, ring = std::sqrt(1 - z * z), a = i * knit::kPi * (3 - std::sqrt(5.0));
        const knit::Vec3 n{ring * std::cos(a), ring * std::sin(a), z};
        lines.fibreStarts.push_back(static_cast<uint32_t>(lines.pos.size() / 3));
        for (int k = 0; k < 3; ++k) {
            const knit::Vec3 p = n * (20 + k), q = n * (21 + k);
            lines.pos.insert(lines.pos.end(), {float(p.x), float(p.y), float(p.z), float(q.x), float(q.y), float(q.z)});
            sum = sum + knit::Vec3{float(p.x), float(p.y), float(p.z)} + knit::Vec3{float(q.x), float(q.y), float(q.z)};
            for (int e = 0; e < 2; ++e) {
                lines.population.push_back(2);
                lines.partColor.insert(lines.partColor.end(), {0.5f, 0.4f, 0.3f});
            }
        }
    }
    const size_t vertices = lines.pos.size() / 3;
    knit::arrangeFibreLines(lines);
    CHECK(lines.fibreStarts.empty() && lines.chunks.size() > 10);
    CHECK(lines.pos.size() == vertices * 3 && lines.population.size() == vertices && lines.partColor.size() == vertices * 3);
    knit::Vec3 after;
    for (size_t v = 0; v < vertices; ++v) after = after + knit::Vec3{lines.pos[3 * v], lines.pos[3 * v + 1], lines.pos[3 * v + 2]};
    CHECK(knit::norm(after - sum) < 1e-3 * double(vertices));
    auto point = [&](uint32_t v) { return knit::Vec3{lines.pos[3 * v], lines.pos[3 * v + 1], lines.pos[3 * v + 2]}; };
    uint32_t next = 0;
    for (const knit::FibreChunk& c : lines.chunks) {
        CHECK(c.firstVertex == next && c.vertexCount % 2 == 0 && c.vertexCount > 0);
        CHECK(knit::norm(c.hi - c.lo) < 25);  // a patch, not two far-apart regions (the hair spans 46 mm)
        for (uint32_t v = c.firstVertex; v < c.firstVertex + c.vertexCount; v += 2) CHECK(std::abs(knit::norm(point(v + 1) - point(v)) - 1) < 1e-4);
        next += c.vertexCount;
    }
    CHECK(next == vertices);

    viewer::Camera cam;
    cam.target = {0, -23, 0};
    cam.dist = 10;
    cam.look({0, 1, 0}, {0, 0, 1});  // a close-up of the -y side
    const viewer::FibreCamera fc = viewer::fibreCamera(cam.proj(1.6), cam.view());
    int skipped = 0;
    for (const knit::FibreChunk& c : lines.chunks) {
        if (viewer::chunkInFrustum(c, fc.clip)) continue;
        ++skipped;
        for (uint32_t v = c.firstVertex; v < c.firstVertex + c.vertexCount; ++v) {
            const knit::Vec3 p = point(v);
            double x[4];
            for (int r = 0; r < 4; ++r) x[r] = fc.clip.m[r] * p.x + fc.clip.m[4 + r] * p.y + fc.clip.m[8 + r] * p.z + fc.clip.m[12 + r];
            CHECK(!(std::abs(x[0]) <= x[3] && std::abs(x[1]) <= x[3] && std::abs(x[2]) <= x[3]));
        }
    }
    CHECK(skipped > 0);
    cam.target = {0, -60, 0};
    cam.look({0, -1, 0}, {0, 0, 1});  // the eye at y = -50, turned away: nothing in view
    const viewer::FibreCamera away = viewer::fibreCamera(cam.proj(1.6), cam.view());
    for (const knit::FibreChunk& c : lines.chunks) CHECK(!viewer::chunkInFrustum(c, away.clip));
}

// ---- the viewport toolbar (2026-11-11) ----

namespace {

// A completed knit whose scene holds the object's source surface and one unknitted part besides the yarn.
void injectSourceKnit(UiFixture& p) {
    p.injectKnit(false);
    auto b = std::make_shared<knit::BuiltScene>(*p.ui.results.knit);
    knit::DisplayMesh base;
    base.kind = knit::DisplayMesh::Kind::Base;
    base.pos = {-20, -20, 0, 20, -20, 0, 0, 20, 60};
    base.idx = {0, 1, 2};
    knit::DisplayMesh plain = base;
    plain.kind = knit::DisplayMesh::Kind::Plain;
    b->meshes.push_back(base);
    b->meshes.push_back(plain);
    b->lo = {-20, -20, 0};
    b->hi = {20, 20, 60};
    p.ui.results.knit = b;
}

bool hasKind(const viewer::ViewContent& c, knit::DisplayMesh::Kind kind) {
    for (const knit::DisplayMesh* d : c.meshes)
        if (d->kind == kind) return true;
    return false;
}

bool startsWith(const std::string& text, const std::string& prefix) { return text.rfind(prefix, 0) == 0; }

// One frame with the panel dragged as wide as it goes (its constraints keep it at 60 % of the window).
void widenPanel(UiFixture& p, float W, float H) {
    ImGui::NewFrame();
    ImGui::SetWindowSize("Workflow", ImVec2(W, H));
    viewer::drawSceneUi(p.ui, W, H);
    ImGui::Render();
}

// The open popover's window lies inside the viewport on both axes (half a pixel of rounding).
bool popoverInsideViewport(UiFixture& p, float W, float H, const char* menu, const std::string& itemInside) {
    if (!clickSuffix(p, W, H, menu)) return false;
    const std::vector<ItemRec> items = settle(p, W, H);
    const ItemRec* r = popoverItem(items, itemInside);
    pressKey(p, W, H, ImGuiKey_Escape);
    const viewer::ViewportRect v = viewer::viewportRect(p.ui, W, H);
    if (!r || !r->window) {
        std::printf("    %s: '%s' not drawn in the popover (viewport %.0f x %.0f)\n", menu, itemInside.c_str(), v.w, v.h);
        return false;
    }
    const ImVec2 lo = r->window->Pos, hi(r->window->Pos.x + r->window->Size.x, r->window->Pos.y + r->window->Size.y);
    std::printf("    %s: popover %.0f..%.0f x %.0f..%.0f, viewport %.0f..%.0f x %.0f..%.0f\n", menu, lo.x, hi.x, lo.y, hi.y, v.x, v.x + v.w, v.y, v.y + v.h);
    return lo.x >= v.x - 0.5f && hi.x <= v.x + v.w + 0.5f && lo.y >= v.y - 0.5f && hi.y <= v.y + v.h + 0.5f;
}

} // namespace

// The presets are orientations: each is recognised again, pan and zoom keep it, an orbit leaves it.
// The 3/4 view looks a little down from the front and the Side preset's side, with world +z up.
TEST(ui_camera_presets) {
    using viewer::CameraView;
    for (const CameraView v : {CameraView::Front, CameraView::Side, CameraView::Top, CameraView::Back, CameraView::ThreeQuarter}) {
        viewer::Camera c;
        viewer::lookFrom(c, v);
        CHECK(viewer::presetOf(c) == v);
        c.pan(30, -12, 500);
        c.zoom(2);
        CHECK(viewer::presetOf(c) == v);
    }
    viewer::Camera c;
    viewer::lookFrom(c, CameraView::Front);
    c.orbit(0, 40);  // tilted only
    CHECK(!viewer::presetOf(c));
    c.level();  // Upright: facing the same way, so Front again
    CHECK(viewer::presetOf(c) == CameraView::Front);
    c.orbit(25, 0);
    CHECK(!viewer::presetOf(c));
    viewer::Camera q;
    viewer::lookFrom(q, CameraView::ThreeQuarter);
    CHECK(q.forward.z < 0 && q.forward.y > 0 && q.forward.x < 0 && std::abs(q.right().z) < 1e-12);
    CHECK(viewer::cameraViewNamed("threequarter") == CameraView::ThreeQuarter);
}

// Every toolbar control is drawn, the camera presets and the display rows in their popovers; the old
// labels are gone and the status bar's interface size is "Interface".
TEST(ui_view_toolbar_inventory) {
    const float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    injectSourceKnit(p);
    p.injectFibres();
    p.ui.view.mode = viewer::ViewMode::Fibres;
    const std::vector<ItemRec> items = settle(p, W, H);
    std::set<std::string> have;
    for (const ItemRec& r : items) have.insert(r.label);
    std::string m = missingLabels(have, {"Source", "Knit", "Fibres", "Fit", "Interface"});
    if (!findSuffix(items, "###viewMenu")) m += " 'View button'";
    if (!findSuffix(items, "###displayMenu")) m += " 'Display button'";
    m += missingLabels(popoverLabels(p, W, H), {"3/4 view", "Front", "Side", "Top", "Back", "Upright", "Studio", "Flat", "Occlusion", "Wireframe",
                                                "Ply texture", "Unknitted parts", "Fibre populations", "Compare with source"});
    for (const char* gone : {"Source mesh", "Show ply texture", "Show unknitted parts", "Show fibre populations", "Frame", "Display"})
        if (have.count(gone)) m += std::string(" (old label ") + gone + ")";
    std::printf("    missing%s\n", m.empty() ? " none" : m.c_str());
    CHECK(m.empty());
}

// Each stage is available when there is something to draw, else says why; choosing one keeps the
// model's camera (one view key). The hair profile: Source is the bare mesh, Knit is skipped, Fibres
// stays available under the coverage overlay, and a Knit chosen before shows the hair.
TEST(ui_view_mode_availability) {
    const float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    using viewer::ViewMode;
    CHECK(!viewer::viewModeState(p.ui, ViewMode::Mesh).available);
    CHECK(viewer::viewModeState(p.ui, ViewMode::Knit).reason == "Not knitted yet: Knit (step 1).");
    CHECK(startsWith(viewer::viewModeState(p.ui, ViewMode::Fibres).reason, "Not knitted yet"));
    injectSourceKnit(p);
    CHECK(viewer::viewModeState(p.ui, ViewMode::Mesh).available && viewer::viewModeState(p.ui, ViewMode::Knit).available);
    CHECK(viewer::viewModeState(p.ui, ViewMode::Fibres).reason == "No fibres for this knit yet: Apply to model (step 2).");
    p.injectFibres();
    CHECK(viewer::viewModeState(p.ui, ViewMode::Fibres).available);
    for (const ViewMode m : {ViewMode::Mesh, ViewMode::Knit, ViewMode::Fibres}) {
        p.ui.view.mode = m;
        CHECK(viewer::viewContent(p.ui).viewKey == "model" && viewer::markedViewMode(p.ui) == m);
    }
    CHECK(click(p, W, H, "Source") && p.ui.view.mode == ViewMode::Mesh);
    const viewer::ViewContent source = viewer::viewContent(p.ui);
    CHECK(hasKind(source, knit::DisplayMesh::Kind::Base) && !hasKind(source, knit::DisplayMesh::Kind::Yarn) && !source.fibres);
    CHECK(click(p, W, H, "Fibres") && p.ui.view.mode == ViewMode::Fibres && viewer::viewContent(p.ui).fibres);

    UiFixture hair("teddy");
    viewer::chooseWorkflow(hair.ui, knit::Workflow::Fur);
    hair.ui.view.fibres = viewer::FibreView::Model;  // the Fur step opens on the reference sphere
    hair.ui.results.fibres = hairScene(hair);
    hair.ui.view.mode = ViewMode::Mesh;
    CHECK(!viewer::viewContent(hair.ui).meshes.empty() && !viewer::viewContent(hair.ui).fibres);
    CHECK(viewer::viewModeState(hair.ui, ViewMode::Knit).reason == "The Fur workflow knits nothing (step 1).");
    CHECK(viewer::viewModeState(hair.ui, ViewMode::Fibres).available);
    hair.ui.view.mode = ViewMode::Knit;
    CHECK(viewer::markedViewMode(hair.ui) == ViewMode::Fibres && viewer::viewContent(hair.ui).fibres);
    hair.ui.view.coverage = true;
    hair.ui.coverage.scene = hairScene(hair);
    CHECK(viewer::viewModeState(hair.ui, ViewMode::Fibres).available);
}

// An open Import file shows its preview whatever stage was chosen: the segment is disabled with the
// reason and marks Source, the stored stage stays, the Display rows of the knit and the comparison are
// off; once the file is closed the stored stage and the comparison are back.
TEST(ui_view_mode_open_file) {
    const float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    injectSourceKnit(p);
    p.injectFibres();
    p.ui.view.mode = viewer::ViewMode::Fibres;
    p.ui.view.compare = true;
    CHECK(viewer::viewportFrame(p.ui).source.has_value());
    CHECK(p.openImport(writeTwoBoxObj("openfile")));
    for (const viewer::ViewMode m : {viewer::ViewMode::Mesh, viewer::ViewMode::Knit, viewer::ViewMode::Fibres})
        CHECK(viewer::viewModeState(p.ui, m).reason == "The open file shows its preview until it is imported (Import step).");
    CHECK(viewer::markedViewMode(p.ui) == viewer::ViewMode::Mesh && p.ui.view.mode == viewer::ViewMode::Fibres);
    const viewer::DisplayRows rows = viewer::displayRows(p.ui);
    CHECK(!rows.plyTexture && !rows.unknittedParts && !rows.fibrePopulations && !rows.compare);
    CHECK(!viewer::viewportFrame(p.ui).source);
    CHECK(click(p, W, H, "Knit") && click(p, W, H, "Source"));  // disabled: nothing changes
    CHECK(p.ui.view.mode == viewer::ViewMode::Fibres);
    viewer::closeImport(p.ui);
    CHECK(viewer::viewModeState(p.ui, viewer::ViewMode::Fibres).available && viewer::markedViewMode(p.ui) == viewer::ViewMode::Fibres);
    CHECK(viewer::viewportFrame(p.ui).source.has_value());
}

// A knit running again keeps the previous fibres, which belong to the previous knit: the Fibres stage
// and the comparison are unavailable while the partial scene shows, running or cancelled; the Knit
// stage shows the partial build. Once the partial scene goes, both are back.
TEST(ui_view_mode_partial_rebuild) {
    const float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    injectSourceKnit(p);
    p.injectFibres();
    p.ui.view.mode = viewer::ViewMode::Fibres;
    p.ui.view.compare = true;
    p.injectRunning(0, "knitting", 1, 4, true);
    CHECK(startsWith(viewer::viewModeState(p.ui, viewer::ViewMode::Fibres).reason, "Knitting again: the fibres belong to the previous knit"));
    CHECK(!viewer::viewContent(p.ui).fibres && !viewer::shownFibres(p.ui));
    CHECK(!viewer::viewportFrame(p.ui).source);
    p.ui.view.mode = viewer::ViewMode::Knit;
    CHECK(viewer::viewModeState(p.ui, viewer::ViewMode::Knit).available);
    CHECK(viewer::displayRows(p.ui).compareWhyNot == "A partial build is not a result." && !viewer::viewportFrame(p.ui).source);
    p.ui.job.active = false;  // cancelled: the partial scene stays on show
    CHECK(startsWith(viewer::viewModeState(p.ui, viewer::ViewMode::Fibres).reason, "Incomplete knit: the fibres belong to the previous knit"));
    CHECK(!viewer::viewportFrame(p.ui).source);
    p.ui.results.partial.reset();
    p.ui.view.mode = viewer::ViewMode::Fibres;
    CHECK(viewer::viewModeState(p.ui, viewer::ViewMode::Fibres).available && viewer::viewportFrame(p.ui).source.has_value());
}

// The Display rows follow the stage on show and keep their values when hidden; Unknitted parts is
// judged on the knit's scene, so switching it off does not hide it.
TEST(ui_display_rows) {
    const float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    injectSourceKnit(p);
    p.ui.view.mode = viewer::ViewMode::Mesh;
    viewer::DisplayRows rows = viewer::displayRows(p.ui);
    CHECK(rows.occlusion && !rows.plyTexture && !rows.unknittedParts && !rows.compare);
    p.ui.view.lighting = viewer::Lighting::Flat;
    CHECK(!viewer::displayRows(p.ui).occlusion);
    p.ui.view.mode = viewer::ViewMode::Knit;
    rows = viewer::displayRows(p.ui);
    CHECK(rows.plyTexture && rows.plyTextureWhyNot.empty() && rows.unknittedParts && rows.compare && rows.compareWhyNot.empty());
    p.ui.view.unknittedParts = false;
    CHECK(viewer::displayRows(p.ui).unknittedParts);
    p.ui.view.plyTexture = true;
    p.ui.view.mode = viewer::ViewMode::Mesh;
    p.ui.view.mode = viewer::ViewMode::Knit;
    CHECK(p.ui.view.plyTexture && !p.ui.view.unknittedParts);
    CHECK(p.setSetting("strand.plies", "2"));
    CHECK(!viewer::displayRows(p.ui).plyTextureWhyNot.empty());
}

// Compare with source: the source's surface on the left, the stage on the right, bounds over both; it
// waits (suspended, still ticked) where it cannot be drawn and never changes the chosen stage. Both
// halves' meshes belong to kept scenes.
TEST(ui_compare_with_source) {
    const float W = 1600, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    injectSourceKnit(p);
    p.ui.view.mode = viewer::ViewMode::Knit;
    CHECK(!viewer::viewportFrame(p.ui).source);
    p.ui.view.compare = true;
    const viewer::ViewportFrame f = viewer::viewportFrame(p.ui);
    CHECK(f.source.has_value());
    if (f.source) {
        CHECK(hasKind(*f.source, knit::DisplayMesh::Kind::Base) && !hasKind(*f.source, knit::DisplayMesh::Kind::Yarn));
        CHECK(hasKind(f.shown, knit::DisplayMesh::Kind::Yarn));
        for (const viewer::Bounds* b : {&f.shown.bounds, &f.source->bounds}) {
            CHECK(f.bounds.lo.x <= b->lo.x && f.bounds.lo.y <= b->lo.y && f.bounds.lo.z <= b->lo.z);
            CHECK(f.bounds.hi.x >= b->hi.x && f.bounds.hi.y >= b->hi.y && f.bounds.hi.z >= b->hi.z);
        }
    }
    CHECK(notHeld(p).empty());
    p.ui.view.mode = viewer::ViewMode::Mesh;
    CHECK(!viewer::viewportFrame(p.ui).source && p.ui.view.compare);
    p.ui.view.mode = viewer::ViewMode::Knit;
    p.ui.results.knit.reset();
    CHECK(!viewer::viewportFrame(p.ui).source && p.ui.view.compare && !viewer::displayRows(p.ui).compareWhyNot.empty());
    injectSourceKnit(p);
    p.ui.stage = viewer::Stage::Fibres;
    p.ui.view.fibres = viewer::FibreView::Reference;
    CHECK(!viewer::viewportFrame(p.ui).source);
    p.ui.stage = viewer::Stage::Knit;
    CHECK(viewer::viewportFrame(p.ui).source.has_value());
    p.ui.view.compare = false;
    CHECK(p.ui.view.mode == viewer::ViewMode::Knit && !viewer::viewportFrame(p.ui).source);
}

// The popovers: a camera preset is requested and closes its popover; Display stays open while its
// rows change (Unknitted parts off and on again in the same popover); Escape, a click in the viewport
// and a stage closes it. F frames the model, except while a field is being typed in.
TEST(ui_view_popovers) {
    const float W = 1600, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    injectSourceKnit(p);
    p.ui.stage = viewer::Stage::Knit;
    p.ui.view.mode = viewer::ViewMode::Knit;
    CHECK(clickSuffix(p, W, H, "###viewMenu") && popoverItem(settle(p, W, H), "3/4 view"));
    std::vector<ItemRec> items = settle(p, W, H);
    const ItemRec* side = popoverItem(items, "Side");
    CHECK(side != nullptr);
    if (side) clickAt(p, W, H, centreOf(*side));
    CHECK(p.ui.camera.view == viewer::CameraView::Side && !popoverItem(settle(p, W, H), "3/4 view"));

    CHECK(clickSuffix(p, W, H, "###displayMenu"));
    for (const char* row : {"Wireframe", "Unknitted parts", "Unknitted parts"}) {
        items = settle(p, W, H);
        const ItemRec* r = popoverItem(items, row);
        CHECK(r != nullptr);
        if (r) clickAt(p, W, H, centreOf(*r));
    }
    CHECK(p.ui.view.wire && p.ui.view.unknittedParts && popoverItem(settle(p, W, H), "Unknitted parts"));
    pressKey(p, W, H, ImGuiKey_Escape);
    CHECK(!popoverItem(settle(p, W, H), "Wireframe"));
    CHECK(clickSuffix(p, W, H, "###displayMenu") && popoverItem(settle(p, W, H), "Wireframe"));
    const viewer::ViewportRect v = viewer::viewportRect(p.ui, W, H);
    clickAt(p, W, H, ImVec2(v.x + v.w * 0.5f, v.y + v.h * 0.5f));
    CHECK(!popoverItem(settle(p, W, H), "Wireframe"));
    CHECK(clickSuffix(p, W, H, "###displayMenu") && popoverItem(settle(p, W, H), "Wireframe"));
    // A click outside a popover closes it and goes no further (ImGui's popup rule): the stage waits for
    // the next click.
    CHECK(click(p, W, H, "Source"));
    CHECK(!popoverItem(settle(p, W, H), "Wireframe") && p.ui.view.mode == viewer::ViewMode::Knit);
    CHECK(click(p, W, H, "Source") && p.ui.view.mode == viewer::ViewMode::Mesh);

    p.ui.camera.frame = false;
    pressKey(p, W, H, ImGuiKey_F);
    CHECK(p.ui.camera.frame);
    p.ui.camera.frame = false;
    // Ctrl+click turns a drag field into text input; then F is a letter typed, not Fit.
    ImGui::GetIO().AddKeyEvent(ImGuiMod_Ctrl, true);
    CHECK(click(p, W, H, "##Yarn diameter"));
    ImGui::GetIO().AddKeyEvent(ImGuiMod_Ctrl, false);
    frame(p, W, H);
    CHECK(ImGui::GetIO().WantTextInput);
    pressKey(p, W, H, ImGuiKey_F);
    CHECK(!p.ui.camera.frame);
    pressKey(p, W, H, ImGuiKey_Escape);
}

// A narrow viewport at 1.5 x: View, Fit and Display wrap as one group, every toolbar item stays inside the
// toolbar, and each popover's whole rectangle stays inside the viewport; in a viewport narrower than the
// groups themselves (240 units) they break into rows rather than leave the toolbar (popoverRect clamps).
TEST(ui_toolbar_narrow) {
    const viewer::ViewportRect unit{100, 50, 300, 200};
    const viewer::PopoverRect left = viewer::popoverRect(unit, 120, 246, 5), right = viewer::popoverRect(unit, 900, 246, 5);
    const viewer::PopoverRect narrow = viewer::popoverRect(viewer::ViewportRect{100, 50, 150, 200}, 250, 246, 5);
    CHECK(left.x == 105 && left.y == 55 && left.maxHeight == 190);
    CHECK(right.x + right.width == 395);
    CHECK(narrow.width == 140 && narrow.x == 105);
    for (const ImVec2 size : {ImVec2(1100, 720), ImVec2(900, 560)}) {
        const float W = size.x, H = size.y;
        ScaledUi ctx(W, H, 1.5f);
        UiFixture p("teddy");
        injectSourceKnit(p);
        p.injectFibres();
        p.ui.stage = viewer::Stage::Knit;
        p.ui.view.mode = viewer::ViewMode::Fibres;
        settle(p, W, H);
        widenPanel(p, W, H);
        const std::vector<ItemRec> items = settle(p, W, H);
        const ItemRec* view = findSuffix(items, "###viewMenu");
        const ItemRec* fit = findItem(items, "Fit");
        const ItemRec* display = findSuffix(items, "###displayMenu");
        CHECK(view && fit && display);
        if (!view || !fit || !display) continue;
        std::printf("  %.0f x %.0f: View y %.0f, Fit y %.0f, Display y %.0f; viewport width %.0f\n", W, H, view->bb.Min.y, fit->bb.Min.y, display->bb.Min.y,
                    viewer::viewportRect(p.ui, W, H).w);
        if (W > 1000) CHECK(std::abs(view->bb.Min.y - fit->bb.Min.y) < 1 && std::abs(fit->bb.Min.y - display->bb.Min.y) < 1);  // one unit
        for (const ItemRec& r : items) {
            if (r.windowName != "##viewToolbar" || !r.window) continue;
            const float lo = r.window->Pos.x, hi = r.window->Pos.x + r.window->Size.x;
            if (r.bb.Min.x < lo - 0.5f || r.bb.Max.x > hi + 0.5f)
                std::printf("    outside the toolbar: %s at %.0f..%.0f y %.0f, window %.0f..%.0f\n", r.label.c_str(), r.bb.Min.x, r.bb.Max.x, r.bb.Min.y, lo, hi);
            CHECK(r.bb.Min.x >= lo - 0.5f && r.bb.Max.x <= hi + 0.5f);
        }
        CHECK(popoverInsideViewport(p, W, H, "###viewMenu", "Front"));
        CHECK(popoverInsideViewport(p, W, H, "###displayMenu", "Wireframe"));
    }
}

int main(int argc, char** argv) {
    int failed = 0, run = 0, skipped = 0;
    const std::string only = argc > 1 ? argv[1] : "";
    if (only == "--list") {
        for (const auto& c : knit_test::registry()) std::printf("%s\n", c.name);
        return 0;
    }
    for (const auto& c : knit_test::registry()) {
        if (!only.empty() && only != c.name) continue;
        ++run;
        try {
            c.fn();
            std::printf("[PASS] %s\n", c.name);
        } catch (const knit_test::Skipped& s) {
            ++skipped;
            std::printf("[SKIP] %s: %s\n", c.name, s.why);
        } catch (const std::exception& e) {
            ++failed;
            std::printf("[FAIL] %s: %s\n", c.name, e.what());
        }
    }
    if (run == 0) {
        std::printf("no test named %s\n", only.c_str());
        return 1;
    }
    std::printf("%d run, %d failed, %d skipped\n", run, failed, skipped);
    if (failed) return 1;
    return skipped == run ? 77 : 0;  // 77: every case run was skipped (CTest SKIP_RETURN_CODE)
}

// The Render preview window: opened from the status bar, it shows the preview's state and the newest
// frame's status line, has Pause and Resolution, keeps the last image's place while the renderer
// refuses a scene, and closes again. The renderer process itself is not started here (main.cpp's
// updateRenderPreview drives the worker; gpu_render_preview_lease_drill runs it).
TEST(ui_render_preview_window) {
    const float W = 1280, H = 900;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(renderedText(p, W, H).find("Render preview") == std::string::npos);
    CHECK(click(p, W, H, "Preview") && p.ui.panels.previewOpen);
    p.ui.renderPreview.state.phase = "refused";
    p.ui.renderPreview.state.message = "insufficient GPU budget: the preview needs 9000 MiB, 3000 MiB available";
    const std::string text = renderedText(p, W, H);
    CHECK(text.find("Render preview") != std::string::npos);
    CHECK(text.find("refused: insufficient GPU budget") != std::string::npos);
    CHECK(text.find("Pause") != std::string::npos && text.find("Resolution") != std::string::npos);
    CHECK(click(p, W, H, "Pause") && p.ui.renderPreview.paused);
    p.ui.renderPreview.state.phase = "rendering";
    p.ui.renderPreview.state.spp = 37;
    p.ui.renderPreview.state.msPerSpp = 21;
    CHECK(renderedText(p, W, H).find("rendering - 37 spp, 21.0 ms per sample") != std::string::npos);
    // What the preview waits for takes the image area (here without an image yet).
    p.ui.renderPreview.busy = "Loading the scene into the renderer...";
    CHECK(renderedText(p, W, H).find("Loading the scene into the renderer...") != std::string::npos);
    p.ui.renderPreview.busy.clear();
    CHECK(click(p, W, H, "Preview") && !p.ui.panels.previewOpen);
    CHECK(renderedText(p, W, H).find("Render preview") == std::string::npos);
}

// ---- the About window ----

// About in the status bar opens the About window with its five tabs, centred on the screen each time it
// opens (a move while open stays), and closes it again; the button
// ends inside the window, right-aligned, at the default size and at 1.5 x scale (the status bar's room
// for its tools counts every button).
TEST(ui_about_window) {
    {
        const float W = 1280, H = 900;
        UiContext ctx(W, H);
        UiFixture p("teddy");
        CHECK(!p.ui.panels.aboutOpen);
        CHECK(click(p, W, H, "About") && p.ui.panels.aboutOpen);
        std::set<std::string> have;
        for (const ItemRec& r : settle(p, W, H)) have.insert(r.label);
        const std::string m = missingLabels(have, {"Wool Shader", "Licence", "Third-party software", "Example models", "Blender and privacy"});
        std::printf("    tabs missing%s\n", m.empty() ? " none" : m.c_str());
        CHECK(m.empty());
        CHECK(click(p, W, H, "Example models") && renderedText(p, W, H).find("hectopod") != std::string::npos);
        // Centred on the screen when it opens; moved while open, it stays where it was put; opened again, centred again.
        const auto centreOffset = [&]() {
            const ImGuiWindow* w = ImGui::FindWindowByName("About Wool Shader");
            if (!w) return 1e9f;
            return std::max(std::abs(w->Pos.x + w->Size.x * 0.5f - W * 0.5f), std::abs(w->Pos.y + w->Size.y * 0.5f - H * 0.5f));
        };
        settle(p, W, H);
        std::printf("    opened: centre %.1f px from the screen's\n", centreOffset());
        CHECK(centreOffset() < 1.0f);
        ImGui::FindWindowByName("About Wool Shader")->Pos = ImVec2(20, 30);
        settle(p, W, H);
        CHECK(ImGui::FindWindowByName("About Wool Shader")->Pos.x == 20);
        CHECK(click(p, W, H, "About") && !p.ui.panels.aboutOpen);
        CHECK(renderedText(p, W, H).find("About Wool Shader") == std::string::npos);
        CHECK(click(p, W, H, "About") && p.ui.panels.aboutOpen);
        settle(p, W, H);
        std::printf("    opened again: centre %.1f px from the screen's\n", centreOffset());
        CHECK(centreOffset() < 1.0f);
    }
    for (const float scale : {1.0f, 1.5f}) {
        const float W = 1600, H = 900;
        ScaledUi ctx(W, H, scale);
        UiFixture p("teddy");
        const std::vector<ItemRec> items = settle(p, W, H);
        const ItemRec* about = findItem(items, "About");
        const ItemRec* report = findItem(items, "Report");
        CHECK(about && report);
        if (!about || !report) continue;
        std::printf("    scale %.1f: Report at x %.0f, About ends at x %.0f of %.0f\n", scale, report->bb.Min.x, about->bb.Max.x, W);
        CHECK(about->bb.Max.x <= W && about->bb.Max.x > W - viewer::uiPx(40));
        CHECK(about->bb.Min.y == report->bb.Min.y);  // one row
    }
}

// The licence texts in the program are the files they come from: the project's LICENSE word for word
// (a stale build would show an old one), and every library cmake/dependencies.cmake fetches has a notice
// with its licence text and source.
TEST(ui_about_licence_texts) {
    namespace fs = std::filesystem;
    std::ifstream in(fs::path(KNIT_SOURCE_DIR) / "LICENSE", std::ios::binary);
    std::string file((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    file.erase(std::remove(file.begin(), file.end(), '\r'), file.end());
    const std::string embedded = viewer::woolShaderLicenceText();
    CHECK(!file.empty() && embedded == file);
    CHECK(embedded.find("MIT License") != std::string::npos && embedded.find("Nicolas Menzel") != std::string::npos);

    std::ifstream deps(fs::path(KNIT_SOURCE_DIR) / "cmake" / "dependencies.cmake");
    const std::string text((std::istreambuf_iterator<char>(deps)), std::istreambuf_iterator<char>());
    const std::regex declare(R"(FetchContent_Declare\((\w+))");
    const std::vector<viewer::ThirdPartyNotice> notices = viewer::thirdPartyNotices();
    std::string missing;
    int declared = 0;
    for (std::sregex_iterator it(text.begin(), text.end(), declare), end; it != end; ++it) {
        ++declared;
        const std::string name = (*it)[1].str();
        bool found = false;
        for (const viewer::ThirdPartyNotice& n : notices) found = found || n.fetchName == name;
        if (!found) missing += " " + name;
    }
    std::printf("    %d fetched libraries, %zu notices; without a notice:%s\n", declared, notices.size(), missing.empty() ? " none" : missing.c_str());
    CHECK(declared >= 7 && missing.empty());
    for (const viewer::ThirdPartyNotice& n : notices) {
        const bool complete = n.text.size() > 200 && n.url.rfind("https://", 0) == 0 && !n.licence.empty() && !n.version.empty() && !n.usedFor.empty();
        if (!complete) std::printf("    incomplete notice: %s (%zu characters of text)\n", n.name.c_str(), n.text.size());
        CHECK(complete);
    }
}

// A source file as text, without carriage returns (the embedded texts have none).
static std::string sourceFileText(const std::string& relativePath) {
    std::ifstream in(std::filesystem::path(KNIT_SOURCE_DIR) / relativePath, std::ios::binary);
    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    text.erase(std::remove(text.begin(), text.end(), '\r'), text.end());
    return text;
}

// The About window names only what the published source holds. No text says tools/ (the author's own scripts, which
// are not published); every blender/ or third_party/ path a text names exists; the one script the Blender export runs is
// in blender/, is the path the export uses, and is GPL-3.0-or-later beside the GPL's text; the SLIM notice carries
// third_party/slim/CHANGES.md word for word, so a new change to those MPL files reaches the window; RIM's notice carries
// the Instant Meshes BSD text its compiled units need, and HiGHS's the notice of its bundled filereaderlp.
TEST(ui_about_public_paths) {
    std::string all = viewer::otherLicencesStatement() + "\n" + viewer::blenderStatement() + "\n" + viewer::privacyStatement() + "\n" +
                      viewer::exampleModelsStatement() + "\n" + viewer::teddyStatement() + "\n";
    std::string slimText, rimText, highsText;
    for (const viewer::ThirdPartyNotice& n : viewer::thirdPartyNotices()) {
        all += n.usedFor + "\n" + n.version + "\n" + n.text + "\n";
        if (n.name.find("SLIM") != std::string::npos) slimText = n.text;
        if (n.fetchName == "rim") rimText = n.text;
        if (n.fetchName == "highs") highsText = n.text;
    }
    for (const knit::ModelCredit& c : viewer::exampleModelCredits()) all += c.title + "\n" + c.source + "\n";
    const size_t privatePath = all.find("tools/");
    if (privatePath != std::string::npos) std::printf("    names tools/: ...%s...\n", all.substr(privatePath, 60).c_str());
    CHECK(privatePath == std::string::npos);

    const std::regex path(R"((?:^|[\s(])((?:blender|third_party)/[\w./-]*\w))");
    int named = 0;
    std::string missing;
    for (std::sregex_iterator it(all.begin(), all.end(), path), end; it != end; ++it) {
        ++named;
        const std::string p = (*it)[1].str();
        if (!std::filesystem::exists(std::filesystem::path(KNIT_SOURCE_DIR) / p)) missing += " " + p;
    }
    std::printf("    %d paths named; missing:%s\n", named, missing.empty() ? " none" : missing.c_str());
    CHECK(named >= 4 && missing.empty());

    const std::string script = "blender/ppt_package_to_blend.py";
    CHECK(viewer::blenderStatement().find(script) != std::string::npos);
    CHECK(sourceFileText("src/export/format_blend.cpp").find("\"" + script + "\"") != std::string::npos);
    CHECK(sourceFileText(script).rfind("# SPDX-License-Identifier: GPL-3.0-or-later\n", 0) == 0);
    const std::string gpl = sourceFileText("blender/LICENSE");
    CHECK(gpl.find("GNU GENERAL PUBLIC LICENSE") != std::string::npos && gpl.find("Version 3, 29 June 2007") != std::string::npos);

    const std::string changes = sourceFileText("third_party/slim/CHANGES.md");
    CHECK(!changes.empty() && slimText.find(changes) != std::string::npos);
    const std::string instantMeshes = sourceFileText("third_party/rim/INSTANT_MESHES_LICENSE.txt");
    CHECK(instantMeshes.find("Wenzel Jakob") != std::string::npos && rimText.find(instantMeshes) != std::string::npos);
    CHECK(highsText.find("Michael Feldmeier") != std::string::npos);
}

// The author's details, and the example models' credits as the registry has them: the teddy bears
// (hectopod, CC BY 4.0) and the wool swatch (CC BY-NC 4.0); the teddy's statement names its author,
// licence and CREDITS.txt.
TEST(ui_about_credits) {
    const viewer::AboutIdentity a = viewer::aboutIdentity();
    CHECK(a.product == "Wool Shader" && a.copyright.find("Nicolas Menzel") != std::string::npos);
    CHECK(a.email.find('@') != std::string::npos && a.portfolioUrl.rfind("https://", 0) == 0 && a.repositoryUrl.rfind("https://", 0) == 0);
    CHECK(!a.build.empty());
    bool teddy = false, swatch = false;
    for (const knit::ModelCredit& c : viewer::exampleModelCredits()) {
        teddy = teddy || (c.author == "hectopod" && c.licence == "CC BY 4.0" && c.source.rfind("https://sketchfab.com/", 0) == 0);
        swatch = swatch || (c.title == "Wool swatch on a needle" && c.licence == "CC BY-NC 4.0");
    }
    CHECK(teddy && swatch);
    const std::string t = viewer::teddyStatement();
    CHECK(t.find("hectopod") != std::string::npos && t.find("https://creativecommons.org/licenses/by/4.0/") != std::string::npos &&
          t.find("CREDITS.txt") != std::string::npos);
    CHECK(viewer::blenderStatement().find("not part of Wool Shader") != std::string::npos);
    CHECK(viewer::privacyStatement().find("no network connections") != std::string::npos);
}

// --- the Mark step (marking_flow.h, panel_mark.cpp)

namespace {

// The two boxes as an imported scene in the Mark step, its pieces read (false: not within 10 s).
bool openMarkStep(UiFixture& p, const char* folder) {
    // The boxes are the same mesh in every folder, so the automatic save goes to this test's own folder: tests that
    // run side by side (ctest -j) would otherwise share one file and remove or overwrite each other's markings.
    const std::filesystem::path saves = std::filesystem::temp_directory_path() / "knit_ui_import" / folder / "markings";
    std::error_code ec;
    std::filesystem::remove_all(saves, ec);
    _putenv_s("KNIT_MARKINGS_CACHE", saves.string().c_str());
    p.ui.s.scene = "import";
    p.ui.s.importSource = writeTwoBoxObj(folder);
    p.ui.s.importUnitMm = 1;
    p.ui.s.importUp = "y";
    p.ui.stage = viewer::Stage::Mark;
    for (int i = 0; i < 2000; ++i) {
        viewer::integrateMarking(p.ui);
        if (p.ui.mark.session.snapshot().pieces && !p.ui.mark.pieceLabels.empty()) return p.ui.mark.set.markings.empty();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    return false;
}

// The first pixel along the viewport's middle row whose surface belongs to the piece (-1: background).
std::optional<viewer::ViewPixel> pixelOfPiece(const viewer::SceneUi& ui, const viewer::Camera& camera, const viewer::ViewPixel& viewport, int wanted) {
    for (double x = 0; x < viewport.width; x += 2) {
        const viewer::ViewPixel pixel{x, viewport.height / 2, viewport.width, viewport.height};
        const std::optional<viewer::SurfaceHit> hit = viewer::surfaceUnderPixel(ui, camera, pixel);
        if ((hit ? hit->piece : -1) == wanted) return pixel;
    }
    return std::nullopt;
}

} // namespace

// The Mark step is a step of the strip: optional for an import, "Imports only" for a preset; its panel lists the
// file's pieces (the two boxes: "left, piece 1" and "right, piece 1") once the file is read.
TEST(ui_mark_stage_listed) {
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture teddy("teddy");
    CHECK(viewer::markState(teddy.ui).word == "Imports only");
    teddy.ui.stage = viewer::Stage::Mark;
    CHECK(findItem(settle(teddy, W, H), "Build stitch mesh") == nullptr);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "mark_listed"));
    CHECK(viewer::markState(p.ui).word == "Optional");
    const std::vector<ItemRec> items = settle(p, W, H);
    CHECK(findItem(items, "left, piece 1") != nullptr && findItem(items, "right, piece 1") != nullptr);
    CHECK(findItem(items, "Build stitch mesh") != nullptr);
    CHECK((p.ui.mark.pieceLabels == std::vector<std::string>{"left, piece 1", "right, piece 1"}));
}

// A click on a piece isolates it; with a marking active, the next click on it places a point there (the file
// triangle under the pixel, weights summing to 1), a click on the other piece isolates that one instead. Ctrl+click
// selects the point and the next click moves it. Control: the background changes nothing.
TEST(ui_mark_pick_adds_point) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "mark_pick"));
    const viewer::ViewPixel viewport{0, 0, 800, 600};
    const viewer::Camera camera = frontCamera(p.ui);
    const std::optional<viewer::ViewPixel> left = pixelOfPiece(p.ui, camera, viewport, 0), right = pixelOfPiece(p.ui, camera, viewport, 1);
    CHECK(left && right);
    viewer::markClick(p.ui, camera, *left, viewer::ClickIntent::Place);
    CHECK(p.ui.mark.piece == 0 && p.ui.mark.set.markings.empty());
    viewer::addMarking(p.ui, knit::MarkingMeaning::Start);
    CHECK(p.ui.mark.set.markings.size() == 1 && p.ui.mark.activeMarking == p.ui.mark.set.markings[0].id);
    const viewer::Camera onLeft = frontCamera(p.ui);  // framed on the isolated piece now
    const std::optional<viewer::ViewPixel> onPiece = pixelOfPiece(p.ui, onLeft, viewport, 0);
    CHECK(onPiece.has_value());
    viewer::markClick(p.ui, onLeft, *onPiece, viewer::ClickIntent::Place);
    const knit::SurfaceMarking& m = p.ui.mark.set.markings[0];
    CHECK(m.points.size() == 1 && m.points[0].tri >= 0 && m.points[0].tri < 12);
    CHECK(std::fabs(m.points[0].bary[0] + m.points[0].bary[1] + m.points[0].bary[2] - 1) < 1e-12);
    CHECK(!p.ui.mark.overlay.pointPos.empty());
    viewer::markClick(p.ui, onLeft, *onPiece, viewer::ClickIntent::Select);
    CHECK(p.ui.mark.selectedPoint == 0);
    const viewer::ViewPixel moved{onPiece->x + 6, onPiece->y, onPiece->width, onPiece->height};
    viewer::markClick(p.ui, onLeft, moved, viewer::ClickIntent::Place);
    CHECK(p.ui.mark.set.markings[0].points.size() == 1 && p.ui.mark.selectedPoint == -1);
    viewer::markClick(p.ui, onLeft, viewer::ViewPixel{0, 0, 800, 600}, viewer::ClickIntent::Place);  // background
    CHECK(p.ui.mark.set.markings[0].points.size() == 1 && p.ui.mark.piece == 0);
    viewer::setNeighbourLook(p.ui, viewer::NeighbourLook::Shown);
    if (const std::optional<viewer::ViewPixel> other = pixelOfPiece(p.ui, onLeft, viewport, 1)) {
        viewer::markClick(p.ui, onLeft, *other, viewer::ClickIntent::Place);
        CHECK(p.ui.mark.piece == 1 && p.ui.mark.set.markings[0].points.size() == 1);
    }
}

// Every edit can be undone and redone, in order; a new edit after an undo drops the redo; Ctrl+Z works from the keys.
TEST(ui_mark_undo_redo) {
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "mark_undo"));
    viewer::isolatePiece(p.ui, 0);
    viewer::addMarking(p.ui, knit::MarkingMeaning::Cut);
    const int id = p.ui.mark.activeMarking;
    viewer::setMarkingMeaning(p.ui, id, knit::MarkingMeaning::RowStartPath);
    CHECK(p.ui.mark.set.markings[0].meaning == knit::MarkingMeaning::RowStartPath && p.ui.mark.undo.size() == 2);
    viewer::undoMarking(p.ui);
    CHECK(p.ui.mark.set.markings[0].meaning == knit::MarkingMeaning::Cut);
    viewer::undoMarking(p.ui);
    CHECK(p.ui.mark.set.markings.empty() && p.ui.mark.redo.size() == 2);
    viewer::redoMarking(p.ui);
    viewer::redoMarking(p.ui);
    CHECK(p.ui.mark.set.markings.size() == 1 && p.ui.mark.set.markings[0].meaning == knit::MarkingMeaning::RowStartPath);
    viewer::undoMarking(p.ui);
    viewer::addMarking(p.ui, knit::MarkingMeaning::Finish);
    CHECK(p.ui.mark.redo.empty() && p.ui.mark.set.markings.size() == 2);
    settle(p, W, H);
    ImGui::GetIO().AddKeyEvent(ImGuiMod_Ctrl, true);
    pressKey(p, W, H, ImGuiKey_Z);
    ImGui::GetIO().AddKeyEvent(ImGuiMod_Ctrl, false);
    frame(p, W, H);
    CHECK(p.ui.mark.set.markings.size() == 1);
}

// A proposal is shown with Accept and Reject; accepting keeps it as input, rejecting keeps it listed but not input;
// editing a proposal makes it the user's marking and keeps where it came from.
TEST(ui_mark_proposal_accept_reject) {
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "mark_proposal"));
    knit::MarkingSet set;
    set.source = p.ui.s.importSource;
    set.digest = p.ui.mark.session.snapshot().digest;
    for (int id = 1; id <= 3; ++id) {
        knit::SurfaceMarking m;
        m.id = id;
        m.meaning = id == 1 ? knit::MarkingMeaning::Start : knit::MarkingMeaning::RowStartPath;
        m.shape = id == 1 ? knit::MarkingShape::ClosedCurve : knit::MarkingShape::Path;
        m.origin = knit::MarkingOrigin::Proposal;
        m.proposal = knit::ProposalState::Offered;
        m.reason = "test proposal";
        m.points.push_back(knit::SurfacePoint{2, {0.2, 0.3, 0.5}});
        set.markings.push_back(m);
    }
    set.nextId = 4;
    const std::string path = (std::filesystem::temp_directory_path() / "knit_ui_import" / "mark_proposal" / "proposal.markings").string();
    std::string err;
    CHECK(knit::saveMarkingFile(path, set, &err));
    viewer::loadMarkingsFrom(p.ui, path);
    CHECK(!p.ui.mark.stale && p.ui.mark.set.markings.size() == 3);
    viewer::isolatePiece(p.ui, 0);
    std::vector<ItemRec> items = settle(p, W, H);
    int accepts = 0;
    for (const ItemRec& r : items) accepts += r.label == "Accept";
    CHECK(accepts == 3);
    viewer::acceptProposal(p.ui, 1);
    viewer::rejectProposal(p.ui, 2);
    CHECK(p.ui.mark.set.markings[0].proposal == knit::ProposalState::Accepted && p.ui.mark.set.markings[1].proposal == knit::ProposalState::Rejected);
    viewer::setMarkingMeaning(p.ui, 3, knit::MarkingMeaning::Cut);
    const knit::SurfaceMarking& edited = p.ui.mark.set.markings[2];
    CHECK(edited.origin == knit::MarkingOrigin::User && edited.proposalId == 3 && edited.proposal == knit::ProposalState::None);
    items = settle(p, W, H);
    accepts = 0;
    for (const ItemRec& r : items) accepts += r.label == "Accept";
    CHECK(accepts == 0);
}

// Save as writes the set; Open reads it back. Every edit is also saved automatically and restored when the file is
// read again. Markings of another mesh are stale: listed, never applied, and no edit or build takes them.
TEST(ui_mark_save_load) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "mark_save"));
    viewer::isolatePiece(p.ui, 1);
    viewer::addMarking(p.ui, knit::MarkingMeaning::Boundary);
    knit::findMarking(p.ui.mark.set, p.ui.mark.activeMarking)->points.push_back(knit::SurfacePoint{14, {0.25, 0.25, 0.5}});
    const std::string path = (std::filesystem::temp_directory_path() / "knit_ui_import" / "mark_save" / "saved.markings").string();
    viewer::saveMarkingsAs(p.ui, path);
    knit::MarkingSet back;
    std::string err;
    CHECK(knit::loadMarkingFile(path, back, &err));
    CHECK(back.markings.size() == 1 && back.markings[0].points.size() == 1 && back.digest == p.ui.mark.session.snapshot().digest);
    CHECK(std::filesystem::exists(viewer::autosavePath(p.ui)));  // the edit before was saved automatically
    back.digest ^= 1;  // another mesh
    CHECK(knit::saveMarkingFile(path, back, &err));
    viewer::loadMarkingsFrom(p.ui, path);
    CHECK(p.ui.mark.stale && viewer::markState(p.ui).word == "Stale");
    const size_t before = p.ui.mark.set.markings.size();
    viewer::addMarking(p.ui, knit::MarkingMeaning::Start);
    CHECK(p.ui.mark.set.markings.size() == before);
    viewer::undoMarking(p.ui);  // back to the set before the load
    CHECK(p.ui.mark.set.markings.size() == 1);
}

// The scripted edit (mark-move): moving a proposal's control point to a file position snaps it to the nearest surface
// point of its piece, makes it the user's marking (its proposal id kept), is undoable and saved automatically.
// Control: an unknown marking or point index changes nothing.
TEST(ui_mark_move_point_edits_proposal) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "mark_move"));
    knit::MarkingSet set;
    set.source = p.ui.s.importSource;
    set.digest = p.ui.mark.session.snapshot().digest;
    knit::SurfaceMarking m;
    m.id = 1;
    m.meaning = knit::MarkingMeaning::Finish;
    m.shape = knit::MarkingShape::Point;
    m.origin = knit::MarkingOrigin::Proposal;
    m.proposal = knit::ProposalState::Offered;
    m.points.push_back(knit::SurfacePoint{2, {0.2, 0.3, 0.5}});
    set.markings.push_back(m);
    set.nextId = 2;
    const std::string path = (std::filesystem::temp_directory_path() / "knit_ui_import" / "mark_move" / "proposal.markings").string();
    std::string err;
    CHECK(knit::saveMarkingFile(path, set, &err));
    viewer::loadMarkingsFrom(p.ui, path);
    viewer::isolatePiece(p.ui, 0);
    const knit::MarkingSession::Snapshot snap = p.ui.mark.session.snapshot();
    const knit::SurfacePiece& piece = snap.pieces->pieces[knit::pieceOfMarking(*snap.pieces, p.ui.mark.set.markings[0])];
    const knit::Vec3 target = piece.mesh.pos[piece.mesh.tri[0][0]] * 0.6 + piece.mesh.pos[piece.mesh.tri[0][1]] * 0.4;
    CHECK(!viewer::moveMarkingPoint(p.ui, 7, 0, target) && !viewer::moveMarkingPoint(p.ui, 1, 3, target));
    CHECK(p.ui.mark.set.markings[0].origin == knit::MarkingOrigin::Proposal);
    CHECK(viewer::moveMarkingPoint(p.ui, 1, 0, target));
    const knit::SurfaceMarking& moved = p.ui.mark.set.markings[0];
    CHECK(moved.origin == knit::MarkingOrigin::User && moved.proposalId == 1 && moved.proposal == knit::ProposalState::None);
    knit::Vec3 at{0, 0, 0};
    for (int k = 0; k < 3; ++k) at += snap.mesh->pos[snap.mesh->tri[moved.points[0].tri][k]] * moved.points[0].bary[k];
    CHECK(knit::norm(at - target) < 1e-9);   // the file point is where the edit put it
    knit::MarkingSet back;
    CHECK(knit::loadMarkingFile(viewer::autosavePath(p.ui), back, &err) && back.markings[0].origin == knit::MarkingOrigin::User);
    viewer::undoMarking(p.ui);
    CHECK(p.ui.mark.set.markings[0].origin == knit::MarkingOrigin::Proposal);
}

// The Mark step's reference choice (smooth-reference report): an edit of the set like any other (undo, autosave), written
// as the v3 file's preparation line and read back by the session; an unknown value is ignored.
TEST(ui_mark_preparation_saved) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "mark_prep"));
    CHECK(p.ui.mark.set.preparation == "facets");
    viewer::setMarkingPreparation(p.ui, "catmull");
    CHECK(p.ui.mark.set.preparation == "facets");
    viewer::setMarkingPreparation(p.ui, "loop_limit");
    CHECK(p.ui.mark.set.preparation == "loop_limit");
    knit::MarkingSet back;
    std::string err;
    CHECK(knit::loadMarkingFile(viewer::autosavePath(p.ui), back, &err) && back.preparation == "loop_limit" && back.version == 5);  // v5: knitting components (v4: assisted stockinette)
    viewer::undoMarking(p.ui);
    CHECK(p.ui.mark.set.preparation == "facets");
}

// The guided route from the Mark step: Knit this piece makes the isolated piece the route's with the saved markings and
// their digest (Knit on, wool fibres), then opens Knit. An edit of the piece's markings changes the digest (the knit
// goes stale); a marking on the other piece does not. (Defect: a route that keeps the old digest.)
TEST(ui_mark_knit_this_piece) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "mark_knit_piece"));
    CHECK(!viewer::knitThisPiece(p.ui));  // nothing isolated
    viewer::isolatePiece(p.ui, 0);
    viewer::addMarking(p.ui, knit::MarkingMeaning::Start);
    CHECK(viewer::placeMarkingPoint(p.ui, knit::Vec3{20, 0, 20}));  // the left box's bottom face
    CHECK(viewer::knitThisPiece(p.ui));
    const knit::SceneSettings& s = p.ui.s;
    CHECK(knit::knitWorkflow(s) && s.knitRoute == "guided" && s.guidedPiece == 0 && p.ui.stage == viewer::Stage::Knit);
    CHECK(knit::knitsGuidedPiece(s));
    CHECK(std::filesystem::path(s.guidedMarkings).is_absolute() && std::filesystem::exists(s.guidedMarkings));
    const knit::MarkingSession::Snapshot snap = p.ui.mark.session.snapshot();
    CHECK(s.guidedMarkingsDigest == knit::identityHex(knit::effectiveMarkingsDigest(p.ui.mark.set, *snap.pieces, 0)));
    const std::string chosen = s.guidedMarkingsDigest;
    viewer::addMarking(p.ui, knit::MarkingMeaning::Finish);  // on piece 0: an empty marking still counts (its meaning)
    CHECK(viewer::placeMarkingPoint(p.ui, knit::Vec3{20, 40, 20}));
    CHECK(s.guidedMarkingsDigest != chosen);
    const std::string edited = s.guidedMarkingsDigest;
    viewer::isolatePiece(p.ui, 1);
    viewer::addMarking(p.ui, knit::MarkingMeaning::Start);
    CHECK(viewer::placeMarkingPoint(p.ui, knit::Vec3{120, 0, 20}));
    CHECK(s.guidedMarkingsDigest == edited);
}

// A Mark build shows as outdated once its markings or the gauge (the Knit settings') differ from what it was made of, and
// current again when they are back. Relax is not an input since assisted stockinette: a layout of today's markings at
// either depth is current and says whether it is relaxed (the automatic row update is unrelaxed). (Defect: the build kept
// on screen without a label after an edit.)
TEST(ui_mark_build_outdated_after_edit) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "mark_outdated"));
    viewer::isolatePiece(p.ui, 0);
    viewer::addMarking(p.ui, knit::MarkingMeaning::Start);
    CHECK(viewer::placeMarkingPoint(p.ui, knit::Vec3{20, 0, 20}));
    viewer::buildMarkedPiece(p.ui);
    for (int i = 0; i < 4000 && !viewer::shownBuild(p.ui); ++i) {
        viewer::integrateMarking(p.ui);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    CHECK(viewer::shownBuild(p.ui) != nullptr);
    CHECK(!viewer::shownBuildOutdated(p.ui));
    p.ui.s.course += 0.4;
    CHECK(viewer::shownBuildOutdated(p.ui));
    p.ui.s.course -= 0.4;
    CHECK(!viewer::shownBuildOutdated(p.ui));
    p.ui.mark.relax = !p.ui.mark.relax;
    CHECK(!viewer::shownBuildOutdated(p.ui));
    p.ui.mark.relax = !p.ui.mark.relax;
    CHECK(viewer::placeMarkingPoint(p.ui, knit::Vec3{10, 0, 10}));  // the start point moved
    CHECK(viewer::shownBuildOutdated(p.ui));
}

// The Knit step of an import: the Route section offers Surface tiling and Guided stitch mesh; the guided route without a
// piece points to the Mark step, and an import that kept the file's surface knits on the guided route without the
// "Rebuild the surface" detour. (Defect: the charts gate blocking the guided route.)
TEST(ui_knit_route_section) {
    namespace fs = std::filesystem;
    const fs::path cache = fs::temp_directory_path() / "knit_ui_import" / "cache_knit_route";
    std::error_code ec;
    fs::remove_all(cache, ec);
    _putenv_s("KNIT_IMPORT_CACHE", cache.string().c_str());
    const float W = 1280, H = 12000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(p.openImport(writeBoxObj("knit_route")));
    viewer::chooseWorkflow(p.ui, knit::Workflow::Fur);  // in stage 2: the file's surface kept
    const bool finished = p.finishImport();  // Fur: the file's surface kept
    viewer::pollSceneUi(p.ui);
    CHECK(finished && p.ui.s.scene == "import" && p.ui.import.lacksCharts);
    viewer::chooseWorkflow(p.ui, knit::Workflow::Knit);
    p.ui.stage = viewer::Stage::Knit;
    std::vector<ItemRec> items = settle(p, W, H);
    CHECK(findItem(items, "Surface tiling") != nullptr && findItem(items, "Guided stitch mesh") != nullptr);
    CHECK(findItem(items, "Rebuild the surface ...") != nullptr && !viewer::knitFooterShown(p.ui));
    p.ui.s.knitRoute = "guided";
    items = settle(p, W, H);
    _putenv_s("KNIT_IMPORT_CACHE", "");
    CHECK(findItem(items, "Open the Mark step") != nullptr);
    CHECK(findItem(items, "Rebuild the surface ...") == nullptr && viewer::knitFooterShown(p.ui));
}

// --- assisted stockinette in the Mark step (suggestions, previews, the quick check)

namespace {

// Two open tubes (r 10.5 mm, 40 mm long, 64 segments, 30 rings: rings - 1 prime), the second 60 mm along x: two pieces.
std::string writeTwoTubesObj(const char* folder) {
    namespace fs = std::filesystem;
    const fs::path d = fs::temp_directory_path() / "knit_ui_import" / folder;
    std::error_code ec;
    fs::remove_all(d, ec);
    fs::create_directories(d);
    std::ostringstream obj;
    const int segments = 64, rings = 30;
    for (int t = 0; t < 2; ++t)
        for (int j = 0; j < rings; ++j)
            for (int i = 0; i < segments; ++i) {
                const double a = 6.283185307179586 * i / segments;
                obj << "v " << 60.0 * t + 10.5 * std::cos(a) << " " << 40.0 * j / (rings - 1) << " " << -10.5 * std::sin(a) << "\n";
            }
    for (int t = 0; t < 2; ++t) {
        obj << (t == 0 ? "o left\n" : "o right\n");
        const int base = t * segments * rings + 1;
        for (int j = 0; j + 1 < rings; ++j)
            for (int i = 0; i < segments; ++i) {
                const int a = base + j * segments + i, b = base + j * segments + (i + 1) % segments, c = a + segments, e = b + segments;
                obj << "f " << a << " " << b << " " << e << "\nf " << a << " " << e << " " << c << "\n";
            }
    }
    std::ofstream(d / "tubes.obj") << obj.str();
    return (d / "tubes.obj").string();
}

bool openMarkTubes(UiFixture& p, const char* folder) {
    const std::filesystem::path saves = std::filesystem::temp_directory_path() / "knit_ui_import" / folder / "markings";
    std::error_code ec;
    std::filesystem::remove_all(saves, ec);
    _putenv_s("KNIT_MARKINGS_CACHE", saves.string().c_str());
    p.ui.s.scene = "import";
    p.ui.s.importSource = writeTwoTubesObj(folder);
    p.ui.s.importUnitMm = 1;
    p.ui.s.importUp = "y";
    p.ui.stage = viewer::Stage::Mark;
    for (int i = 0; i < 2000; ++i) {
        viewer::integrateMarking(p.ui);
        if (p.ui.mark.session.snapshot().pieces && p.ui.mark.pieceLabels.size() == 2) return true;
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    return false;
}

// Frames of the Mark step's background work until `done` (at most `seconds`).
bool markUntil(UiFixture& p, double seconds, bool (*done)(const viewer::SceneUi&)) {
    const auto end = std::chrono::steady_clock::now() + std::chrono::milliseconds(static_cast<int>(seconds * 1000));
    while (std::chrono::steady_clock::now() < end) {
        viewer::integrateMarking(p.ui);
        if (done(p.ui)) return true;
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    return done(p.ui);
}

bool hasCards(const viewer::SceneUi& ui) { return ui.mark.proposals && ui.mark.proposals->piece == ui.mark.piece && ui.mark.cardSetup.size() == ui.mark.proposals->setups.size(); }
bool quickNow(const viewer::SceneUi& ui) { return viewer::currentQuick(ui) != nullptr; }
bool deepIdle(const viewer::SceneUi& ui) { return ui.mark.session.snapshot().phase == knit::MarkingSession::Phase::Ready; }

// Suggestions for piece `piece` of the open tubes, arrived and merged.
bool suggestFor(UiFixture& p, int piece) {
    viewer::isolatePiece(p.ui, piece);
    const std::shared_ptr<const knit::ProposalSet> before = p.ui.mark.proposals;
    viewer::suggestMarkings(p.ui);
    for (int i = 0; i < 5000; ++i) {
        viewer::integrateMarking(p.ui);
        if (p.ui.mark.proposals != before && hasCards(p.ui)) return true;
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    return false;
}

int setupCount(const knit::MarkingSet& s, int setup, knit::ProposalState state) {
    int n = 0;
    for (const knit::SurfaceMarking& m : s.markings) n += m.setup == setup && m.proposal == state ? 1 : 0;
    return n;
}

bool findingWith(const knit::MarkingReview& r, const std::string& text) {
    for (const knit::MarkingFinding& f : r.findings)
        if (f.message.find(text) != std::string::npos) return true;
    return false;
}

} // namespace

// Suggestions arrive as offered markings (one undo step); a preview changes nothing saved (the set, the automatic save, the
// undo depth) and is reviewed on its own set; accepting a setup is one undo step that accepts exactly its three markings;
// undo and redo restore and repeat it. (Defect: accepting marking by marking, or a preview written into the saved set.)
TEST(ui_mark_setup_preview_accept_undo) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "mark_setup_accept"));
    CHECK(suggestFor(p, 0));
    CHECK(p.ui.mark.proposals->shape == knit::PieceShape::OpenTube && p.ui.mark.proposals->setups.size() == 2);
    const size_t undoAfterSuggest = p.ui.mark.undo.size();
    CHECK(undoAfterSuggest == 1);
    CHECK(viewer::suggestionState(p.ui, 1) == "offered" && viewer::suggestionState(p.ui, 2) == "offered");
    const std::string saved = knit::writeMarkings(p.ui.mark.set);
    std::string file;
    {
        std::ifstream f(viewer::autosavePath(p.ui), std::ios::binary);
        file.assign(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
    }
    viewer::previewSuggestion(p.ui, 1);
    CHECK(p.ui.mark.previewSetup == p.ui.mark.cardSetup[0]);
    std::vector<int> ids;
    const knit::MarkingSet eff = viewer::effectiveMarkingSet(p.ui, ids);
    CHECK(ids.size() == 3 && knit::writeMarkings(p.ui.mark.set) == saved && p.ui.mark.undo.size() == undoAfterSuggest);
    {
        std::ifstream f(viewer::autosavePath(p.ui), std::ios::binary);
        CHECK(std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()) == file);
    }
    CHECK(markUntil(p, 10, quickNow));
    CHECK(viewer::currentQuick(p.ui)->review.depth == knit::ReviewDepth::Quick && findingWith(viewer::currentQuick(p.ui)->review, "Research preview"));
    CHECK(viewer::acceptSuggestion(p.ui, 1, false));
    CHECK(p.ui.mark.previewSetup == 0 && p.ui.mark.undo.size() == undoAfterSuggest + 1);
    CHECK(viewer::suggestionState(p.ui, 1) == "accepted" && setupCount(p.ui.mark.set, p.ui.mark.cardSetup[0], knit::ProposalState::Accepted) == 3);
    CHECK(viewer::suggestionState(p.ui, 2) == "offered");
    viewer::undoMarking(p.ui);
    CHECK(viewer::suggestionState(p.ui, 1) == "offered" && setupCount(p.ui.mark.set, p.ui.mark.cardSetup[0], knit::ProposalState::Accepted) == 0);
    viewer::redoMarking(p.ui);
    CHECK(viewer::suggestionState(p.ui, 1) == "accepted");
}

// With setup A accepted, new suggestions keep A accepted; a preview of the other direction B is reviewed and built with
// B's start, termination and path alone (no "only one start"), and leaving it restores everything: the saved set, the undo
// depth, the guided route's digest, A in use. (Defect: a preview that adds B beside A.)
TEST(ui_mark_preview_after_accept) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "mark_preview_after_accept"));
    CHECK(suggestFor(p, 0));
    CHECK(viewer::acceptSuggestion(p.ui, 1, false));
    const int setupA = p.ui.mark.cardSetup[0];
    CHECK(viewer::knitThisPiece(p.ui));
    p.ui.stage = viewer::Stage::Mark;
    const std::string digest = p.ui.s.guidedMarkingsDigest;
    CHECK(suggestFor(p, 0));   // regenerate
    CHECK(setupCount(p.ui.mark.set, setupA, knit::ProposalState::Accepted) == 3 && p.ui.s.guidedMarkingsDigest == digest);
    const std::string saved = knit::writeMarkings(p.ui.mark.set);
    const size_t undo = p.ui.mark.undo.size();
    viewer::previewSuggestion(p.ui, 2);   // the other direction
    CHECK(p.ui.mark.previewSetup == p.ui.mark.cardSetup[1]);
    CHECK(markUntil(p, 10, quickNow));
    CHECK(!findingWith(viewer::currentQuick(p.ui)->review, "Only one"));
    viewer::buildMarkedPiece(p.ui);
    CHECK(markUntil(p, 60, [](const viewer::SceneUi& ui) { return deepIdle(ui) && viewer::shownBuild(ui) && !viewer::shownBuild(ui)->previewIds.empty(); }));
    CHECK(!findingWith(viewer::shownBuild(p.ui)->review, "Only one") && viewer::shownBuild(p.ui)->review.valid);
    viewer::previewSuggestion(p.ui, 0);
    CHECK(p.ui.mark.previewSetup == 0 && knit::writeMarkings(p.ui.mark.set) == saved && p.ui.mark.undo.size() == undo && p.ui.s.guidedMarkingsDigest == digest);
    CHECK(setupCount(p.ui.mark.set, setupA, knit::ProposalState::Accepted) == 3);
}

// After Knit this piece, suggestions arriving, regenerating them and previewing one keep the guided route's digest (so the
// Knit key); accepting another setup on the piece or moving a point changes it; suggestions and edits on the other piece do
// not. (Defect: offered markings or a preview counted as the Knit input.)
TEST(ui_mark_proposals_keep_knit_route) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "mark_keep_route"));
    CHECK(suggestFor(p, 0));
    CHECK(viewer::acceptSuggestion(p.ui, 1, false));
    CHECK(viewer::knitThisPiece(p.ui));
    p.ui.stage = viewer::Stage::Mark;
    const std::string digest = p.ui.s.guidedMarkingsDigest;   // the Knit key's markings input (guided_stage_boundaries)
    CHECK(suggestFor(p, 0));
    viewer::previewSuggestion(p.ui, 2);
    viewer::previewSuggestion(p.ui, 0);
    CHECK(p.ui.s.guidedMarkingsDigest == digest);
    CHECK(suggestFor(p, 1));   // the other piece
    CHECK(viewer::acceptSuggestion(p.ui, 1, false));
    CHECK(p.ui.s.guidedMarkingsDigest == digest);
    viewer::isolatePiece(p.ui, 0);
    CHECK(suggestFor(p, 0));
    CHECK(viewer::acceptSuggestion(p.ui, 2, false));   // the other direction on the knitted piece
    CHECK(p.ui.s.guidedMarkingsDigest != digest);
}

// A gauge change outdates the offered suggestions: they say so, cannot be accepted (no Accept button), and come back when
// the gauge does. (Defect: a suggestion computed for another gauge accepted as if current.)
TEST(ui_mark_setup_outdated_basis) {
    const float W = 1280, H = 1400;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "mark_outdated_basis"));
    CHECK(suggestFor(p, 0));
    CHECK(findItem(settle(p, W, H), "Accept") != nullptr);
    p.ui.s.course += 0.4;
    CHECK(viewer::suggestionState(p.ui, 1) == "outdated" && !viewer::acceptSuggestion(p.ui, 1, false));
    CHECK(findItem(settle(p, W, H), "Accept") == nullptr);
    p.ui.s.course -= 0.4;
    CHECK(viewer::suggestionState(p.ui, 1) == "offered");
}

// The quick check says what it knows: before any layout the rows are "not computed" (never 0); once the automatic layout
// ran they are current; after the start moves (with the automatic update off) they are outdated and the rows "pending".
// The cells, coverage and checks are never shown by the quick check. (Defect: unknown measures shown as 0 or PASS, or a
// path judged against an old start's rows.)
TEST(ui_mark_quick_feedback_states) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "mark_quick_states"));
    viewer::setLayoutRefresh(p.ui, false);
    CHECK(suggestFor(p, 0));
    CHECK(viewer::acceptSuggestion(p.ui, 1, false));
    CHECK(markUntil(p, 10, quickNow));
    std::string text = viewer::quickReviewText(p.ui);
    CHECK(text.find("rows: not computed") != std::string::npos && text.find("cells without a construction: not computed") != std::string::npos &&
          text.find("layout checks: not computed") != std::string::npos);
    viewer::setLayoutRefresh(p.ui, true);
    CHECK(markUntil(p, 60, [](const viewer::SceneUi& ui) {
        const std::shared_ptr<const knit::MarkingSession::QuickView> q = viewer::currentQuick(ui);
        return q && q->review.measures.rowsKnown == knit::Known::Current;
    }));
    text = viewer::quickReviewText(p.ui);
    CHECK(text.find("rows: 12") != std::string::npos && text.find("coverage: not computed") != std::string::npos);
    viewer::setLayoutRefresh(p.ui, false);
    int start = 0;
    for (const knit::SurfaceMarking& m : p.ui.mark.set.markings)
        if (m.meaning == knit::MarkingMeaning::Start && m.proposal == knit::ProposalState::Accepted) start = m.id;
    // the start's first point moved a few mm up into the tube: the start is no longer the opening
    const knit::SurfaceMarking* sm = knit::findMarking(p.ui.mark.set, start);
    CHECK(sm != nullptr);
    knit::Vec3 at{};
    for (size_t i = 0; i < p.ui.mark.set.markings.size(); ++i)
        if (p.ui.mark.set.markings[i].id == start) at = p.ui.mark.traces[i].points[p.ui.mark.traces[i].controlAt[0]].pos;
    const double inward = at.y < 1 ? 4 : -4;
    CHECK(viewer::moveMarkingPoint(p.ui, start, 0, knit::Vec3{at.x, at.y + inward, at.z}));
    CHECK(markUntil(p, 10, quickNow));
    text = viewer::quickReviewText(p.ui);
    CHECK(text.find("(outdated)") != std::string::npos || text.find("rows: not computed") != std::string::npos);
    CHECK(text.find("pending") != std::string::npos || text.find("need a layout") != std::string::npos || text.find("Row checks") != std::string::npos);
}

// A set with setups, an accepted, a rejected and an edited suggestion is saved (v4), opened again and equal; a v3 file
// opens unchanged. (Defect: setups or bases lost on the way through the file.)
TEST(ui_mark_save_load_setups) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "mark_save_setups"));
    CHECK(suggestFor(p, 0));
    CHECK(viewer::acceptSuggestion(p.ui, 1, false));
    viewer::rejectSuggestion(p.ui, 2);
    int path = 0;
    for (const knit::SurfaceMarking& m : p.ui.mark.set.markings)
        if (m.meaning == knit::MarkingMeaning::RowStartPath && m.proposal == knit::ProposalState::Accepted) path = m.id;
    const knit::SurfaceMarking* pm = knit::findMarking(p.ui.mark.set, path);
    CHECK(pm != nullptr && pm->points.size() >= 2);
    CHECK(viewer::moveMarkingPoint(p.ui, path, 1, knit::Vec3{10.5 * std::cos(0.4), 20, -10.5 * std::sin(0.4)}));
    CHECK(knit::findMarking(p.ui.mark.set, path)->origin == knit::MarkingOrigin::User && knit::findMarking(p.ui.mark.set, path)->setup > 0);
    const std::string file = (std::filesystem::temp_directory_path() / "knit_ui_import" / "mark_save_setups" / "saved.markings").string();
    viewer::saveMarkingsAs(p.ui, file);
    const std::string before = knit::writeMarkings(p.ui.mark.set);
    viewer::loadMarkingsFrom(p.ui, file);
    CHECK(knit::writeMarkings(p.ui.mark.set) == before && p.ui.mark.set.version == 5);
    std::string v3 = before;
    v3.replace(0, std::string("knit-markings v4").size(), "knit-markings v3");
    std::string stripped;
    std::istringstream lines(v3);
    for (std::string l; std::getline(lines, l);)
        if (l.rfind("setup = ", 0) != 0 && l.rfind("basis = ", 0) != 0) stripped += l + "\n";
    std::ofstream(file, std::ios::binary) << stripped;
    viewer::loadMarkingsFrom(p.ui, file);
    CHECK(p.ui.mark.set.version == 3 && p.ui.mark.set.markings.size() > 0);
}

// A finding's highlight lands on the review's own trace: the quick check of a path crossing a row twice names the row, and
// the highlighted overlay holds exactly the reported places and segments of that review's trace. (Defect: trace ranges
// applied to the viewer's own coarser trace, highlighting elsewhere.)
TEST(ui_mark_highlight_on_review_traces) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "mark_highlight"));
    CHECK(suggestFor(p, 0));
    CHECK(viewer::acceptSuggestion(p.ui, 1, false));
    CHECK(markUntil(p, 60, [](const viewer::SceneUi& ui) {
        const std::shared_ptr<const knit::MarkingSession::QuickView> q = viewer::currentQuick(ui);
        return q && q->review.measures.rowsKnown == knit::Known::Current;
    }));
    viewer::setLayoutRefresh(p.ui, false);
    int path = 0;
    for (const knit::SurfaceMarking& m : p.ui.mark.set.markings)
        if (m.meaning == knit::MarkingMeaning::RowStartPath && m.proposal == knit::ProposalState::Accepted) path = m.id;
    // three control points: up past a row, back below it, up again
    knit::MarkingSet edited = p.ui.mark.set;
    const knit::MarkingSession::Snapshot snap = p.ui.mark.session.snapshot();
    const knit::SurfacePiece& piece = snap.pieces->pieces[0];
    knit::SurfaceMarking* m = knit::findMarking(edited, path);
    m->points.clear();
    for (const double y : {1.0, 21.5, 18.5, 39.0}) {
        knit::SurfacePoint lp = knit::closestSurfacePoint(piece.mesh, knit::Vec3{10.5 * std::cos(0.3), y, -10.5 * std::sin(0.3)});
        lp.tri = piece.fileTriangle[lp.tri];
        m->points.push_back(lp);
    }
    p.ui.mark.undo.push_back(p.ui.mark.set);
    p.ui.mark.set = edited;
    viewer::undoMarking(p.ui);
    viewer::redoMarking(p.ui);   // the edited set, through the ordinary change path
    CHECK(markUntil(p, 10, quickNow));
    const std::shared_ptr<const knit::MarkingSession::QuickView> q = viewer::currentQuick(p.ui);
    int index = -1;
    for (size_t i = 0; i < q->review.findings.size(); ++i)
        if (q->review.findings[i].message.find("crosses row") != std::string::npos && q->review.findings[i].firstTracePoint >= 0) index = static_cast<int>(i);
    CHECK(index >= 0);
    if (index < 0) return;
    const knit::MarkingFinding& f = q->review.findings[index];
    viewer::highlightFinding(p.ui, index, true);
    const viewer::OverlayGeometry& g = p.ui.mark.overlay;
    // the places: every reported position is among the overlay's points (several crossings may lie close together)
    int places = 0;
    for (const knit::Vec3& x : f.at) {
        bool shown = false;
        for (size_t k = 0; k + 2 < g.pointPos.size() && !shown; k += 3)
            shown = std::fabs(g.pointPos[k] - x.x) < 1e-4 && std::fabs(g.pointPos[k + 1] - x.y) < 1e-4 && std::fabs(g.pointPos[k + 2] - x.z) < 1e-4;
        places += shown ? 1 : 0;
    }
    std::printf("  highlighted places %d of %zu (row %d, trace points %d-%d)\n", places, f.at.size(), f.row, f.firstTracePoint, f.lastTracePoint);
    CHECK(places == static_cast<int>(f.at.size()) && places >= 2);
    // the segments: every trace point of the range appears among the highlighted line vertices
    const knit::MarkingTrace* tr = nullptr;
    for (size_t i = 0; i < q->markingIds.size(); ++i)
        if (q->markingIds[i] == f.markingId) tr = &q->traces[i];
    CHECK(tr != nullptr);
    if (!tr) return;
    for (int k = f.firstTracePoint + 1; k <= f.lastTracePoint && k < static_cast<int>(tr->points.size()); ++k) {
        bool found = false;
        for (size_t v = 0; v + 2 < g.linePos.size() && !found; v += 3)
            found = std::fabs(g.linePos[v] - tr->points[k].pos.x) < 1e-4 && std::fabs(g.linePos[v + 1] - tr->points[k].pos.y) < 1e-4 &&
                    std::fabs(g.linePos[v + 2] - tr->points[k].pos.z) < 1e-4;
        CHECK(found);
    }
}

// A quick view that does not answer the newest request is not drawn as current ("checking"): with the reviews held back,
// an edit after a request leaves the older answer unshown, the newer one is shown when it arrives, and the latency record
// counts the older request as coalesced, not as answered. (Defect: an older answer drawn over a newer edit.)
TEST(ui_mark_quick_display_matches_request) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "mark_quick_latest"));
    CHECK(suggestFor(p, 0));
    CHECK(markUntil(p, 10, quickNow));
    viewer::markFrameShown(p.ui);
    const size_t answered = p.ui.mark.latencyMs.size();
    const int coalesced = p.ui.mark.coalesced;
    p.ui.mark.session.delayQuickForTest(250);
    viewer::previewSuggestion(p.ui, 1);    // request A
    viewer::previewSuggestion(p.ui, 2);    // request B: A is now older
    for (int i = 0; i < 100; ++i) {
        viewer::integrateMarking(p.ui);
        const std::shared_ptr<const knit::MarkingSession::QuickView> q = viewer::currentQuick(p.ui);
        CHECK(!q || q->serial == p.ui.mark.quickSerial);
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    CHECK(markUntil(p, 10, quickNow));
    viewer::markFrameShown(p.ui);
    CHECK(viewer::currentQuick(p.ui)->serial == p.ui.mark.quickSerial);
    CHECK(p.ui.mark.latencyMs.size() == answered + 1 && p.ui.mark.coalesced >= coalesced + 1);
    p.ui.mark.session.delayQuickForTest(0);
}

// --- the Stitch Mesh view

namespace {

// Piece `piece` of the open tubes with its first suggestion accepted; automatic layouts are off, so only explicit builds
// run.
bool acceptedTube(UiFixture& p, int piece) {
    viewer::setLayoutRefresh(p.ui, false);
    return suggestFor(p, piece) && viewer::acceptSuggestion(p.ui, 1, false);
}

// "Build stitch mesh" on the isolated piece; returns its request number.
uint64_t pressBuild(UiFixture& p) {
    viewer::buildMarkedPiece(p.ui);
    return p.ui.mark.explicitBuild;
}

// Frames until the session has published the build that answers `request` (at most `seconds`).
bool buildArrived(UiFixture& p, uint64_t request, double seconds) {
    const auto end = std::chrono::steady_clock::now() + std::chrono::milliseconds(static_cast<int>(seconds * 1000));
    while (std::chrono::steady_clock::now() < end) {
        viewer::integrateMarking(p.ui);
        for (const viewer::MarkingFlow::ArrivedBuild& a : p.ui.mark.builds)
            if (a.view->request == request) return true;
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
    return false;
}

bool contentHas(const viewer::ViewContent& c, const knit::DisplayMesh* mesh) { return std::find(c.meshes.begin(), c.meshes.end(), mesh) != c.meshes.end(); }

bool ownersHave(const viewer::SceneUi& ui, const std::shared_ptr<const knit::BuiltScene>& scene) {
    const viewer::MeshOwners o = viewer::meshOwners(ui);
    return std::find(o.scenes.begin(), o.scenes.end(), scene) != o.scenes.end();
}

} // namespace

// In the Mark step the toolbar offers Markings, Rows and Stitch Mesh in place of Source, Knit, Fibres: Stitch Mesh is open
// as soon as the file is read and says when there is nothing to show yet; Rows waits for a build. (Defect: the view
// hidden at the bottom of the panel and disabled before a build.)
TEST(ui_mark_stitch_view_in_toolbar) {
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "stitch_toolbar"));
    viewer::isolatePiece(p.ui, 0);
    const std::vector<ItemRec> items = settle(p, W, H);
    const ItemRec* stitch = findItem(items, "Stitch Mesh");
    const ItemRec* rows = findItem(items, "Rows");
    CHECK(stitch && !stitch->disabled && rows && rows->disabled && findItem(items, "Markings") && !findItem(items, "Source"));
    CHECK(click(p, W, H, "Stitch Mesh") && p.ui.mark.view == viewer::MarkView::StitchMesh);
    CHECK(renderedText(p, W, H).find("No stitch mesh for this piece yet") != std::string::npos);
    CHECK(viewer::viewContent(p.ui).caption.find("No stitch mesh for this piece yet") != std::string::npos);
}

// Build stitch mesh selects the Stitch Mesh view when its result arrives for the piece still shown: the stitch faces by
// kind with their polygon outlines, the piece's surface under them, the headline in the caption, and not one knit job.
// An automatic layout after an edit never selects the view. (Defect: the view chosen by hand after every build, or a
// view that needs yarn.)
TEST(ui_mark_stitch_view_after_build) {
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "stitch_after_build"));
    CHECK(acceptedTube(p, 0));
    // an automatic layout: arrives, the view stays Markings
    viewer::setLayoutRefresh(p.ui, true);
    CHECK(markUntil(p, 120, [](const viewer::SceneUi& ui) { return viewer::shownBuild(ui) != nullptr && deepIdle(ui); }));
    CHECK(p.ui.mark.view == viewer::MarkView::Markings);
    viewer::setLayoutRefresh(p.ui, false);
    // the explicit build
    const uint64_t request = pressBuild(p);
    CHECK(request > 0 && buildArrived(p, request, 120));
    const knit::MarkedPieceView* b = viewer::shownBuild(p.ui);
    CHECK(b && b->request == request && b->use.drawable && p.ui.mark.view == viewer::MarkView::StitchMesh);
    if (!b || !b->use.drawable) return;
    const viewer::ViewContent c = viewer::viewContent(p.ui);
    CHECK(!contentHas(c, &b->scene->meshes[0]) && contentHas(c, &b->scene->meshes[1]) && !contentHas(c, &b->scene->meshes[2]));
    CHECK(p.ui.mark.overlay.linePos.size() >= 6 * b->stitch.corner.size());   // every polygon edge drawn
    CHECK(c.caption.find(knit::stitchMeshHeadline(b->use, b->quality)) != std::string::npos);
    std::printf("  caption: %s\n", c.caption.c_str());
    // the switches
    p.ui.mark.stitchHighlight = true;
    p.ui.mark.stitchSource = true;
    const viewer::ViewContent h = viewer::viewContent(p.ui);
    CHECK(contentHas(h, &b->scene->meshes[2]) && !contentHas(h, &b->scene->meshes[1]) && contentHas(h, &b->scene->meshes[0]));
    const std::string text = renderedText(p, W, H);
    CHECK(text.find("Export Stitch Mesh ...") != std::string::npos && text.find("Highlight findings") != std::string::npos);
    CHECK(p.jobs->received.empty() && !p.ui.results.knit);
}

// A build whose result arrives after the user went to another piece is kept for its own piece and leaves the view alone.
// (Defect: one flag that switches the view on whichever result comes next.)
TEST(ui_mark_stitch_view_not_stolen) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "stitch_not_stolen"));
    CHECK(acceptedTube(p, 0));
    const uint64_t request = pressBuild(p);
    viewer::isolatePiece(p.ui, 1);   // before the result is taken over
    CHECK(buildArrived(p, request, 120));
    CHECK(p.ui.mark.view == viewer::MarkView::Markings && viewer::shownBuild(p.ui) == nullptr && p.ui.mark.explicitBuild == 0);
    viewer::isolatePiece(p.ui, 0);
    CHECK(viewer::shownBuild(p.ui) && viewer::shownBuild(p.ui)->request == request && p.ui.mark.view == viewer::MarkView::Markings);
}

// Build piece 0, build piece 1, go back to piece 0: the same mesh of piece 0 is shown, its scene stays among the uploaded
// owners, its drawn geometry maps back to the saved vertices through the import's conversion, and Export Stitch Mesh
// writes piece 0. (Defect: only the newest build drawn or kept.)
TEST(ui_mark_stitch_view_switch_back) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "stitch_switch_back"));
    CHECK(acceptedTube(p, 0));
    const uint64_t first = pressBuild(p);
    CHECK(buildArrived(p, first, 120));
    const std::shared_ptr<const knit::MarkedPieceView> body = p.ui.mark.builds.front().view;
    CHECK(body->piece == 0 && body->use.drawable);
    CHECK(acceptedTube(p, 1));
    const uint64_t second = pressBuild(p);
    CHECK(buildArrived(p, second, 120) && viewer::shownBuild(p.ui) && viewer::shownBuild(p.ui)->piece == 1);
    viewer::setLayoutRefresh(p.ui, true);   // the automatic layout after edits: never for a piece whose build is current
    viewer::isolatePiece(p.ui, 0);
    const auto settled = std::chrono::steady_clock::now() + std::chrono::milliseconds(1500);
    while (std::chrono::steady_clock::now() < settled) {
        viewer::integrateMarking(p.ui);
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    const knit::MarkedPieceView* shown = viewer::shownBuild(p.ui);
    CHECK(shown == body.get() && p.ui.mark.builds.size() == 2 && p.ui.mark.layoutDueAt == 0);
    viewer::setMarkView(p.ui, viewer::MarkView::StitchMesh);
    const viewer::ViewContent c = viewer::viewContent(p.ui);
    CHECK(contentHas(c, &body->scene->meshes[1]) && ownersHave(p.ui, body->scene));
    // the drawn file coordinates go to the library through the view's transform and back with the shared conversion
    CHECK(c.preview.has_value());
    if (!c.preview) return;
    const viewer::PreviewTransform& placed = *c.preview;
    double worst = 0;
    for (const knit::Vec3& v : body->stitch.pos) {
        const knit::Vec3 lib = viewer::toLibrary(placed, v);
        const knit::Vec3 back = knit::toFileCoordinates(placed.scale, placed.offset, placed.yUp ? knit::UpAxis::Y : knit::UpAxis::Z, lib);
        worst = std::max(worst, knit::norm(back - v));
    }
    CHECK(worst < 1e-9);
    const std::string path = (std::filesystem::temp_directory_path() / "knit_ui_import" / "stitch_switch_back" / "piece0.obj").string();
    CHECK(viewer::exportStitchMesh(p.ui, path));
    knit::StitchMeshData back;
    std::string err;
    CHECK(knit::readPolygonObj(path, back, &err) && back.faces() == body->stitch.faces() && back.pos == body->stitch.pos && back.corner == body->stitch.corner);
}

// After a good build, a rebuild without the row-start path returns without faces: the piece's earlier mesh is gone, the
// piece is drawn as imported, the reason is the backend's own finding, Export Stitch Mesh is refused with that reason, and
// the failed build never becomes the comparison baseline. (Defect: the old mesh left on screen.)
TEST(ui_mark_stitch_view_empty_rebuild) {
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "stitch_empty_rebuild"));
    CHECK(acceptedTube(p, 0));
    const uint64_t good = pressBuild(p);
    CHECK(buildArrived(p, good, 120));
    const std::shared_ptr<const knit::BuiltScene> goodScene = p.ui.mark.builds.front().view->scene;
    CHECK(p.ui.mark.view == viewer::MarkView::StitchMesh && p.ui.mark.baseline.get() == viewer::shownBuild(p.ui));
    for (const knit::SurfaceMarking& m : p.ui.mark.set.markings)
        if (m.meaning == knit::MarkingMeaning::RowStartPath && m.proposal == knit::ProposalState::Accepted) {
            viewer::deleteMarking(p.ui, m.id);
            break;
        }
    const uint64_t empty = pressBuild(p);
    CHECK(buildArrived(p, empty, 60));
    const knit::MarkedPieceView* b = viewer::shownBuild(p.ui);
    CHECK(b && b->request == empty && !b->use.drawable && b->scene->meshes.size() == 1 && p.ui.mark.builds.size() == 1);
    if (!b) return;
    bool named = false;
    for (const knit::MarkingFinding& f : b->review.findings) named = named || f.message == b->use.reason;
    CHECK(named);
    const viewer::ViewContent c = viewer::viewContent(p.ui);
    const knit::MarkingSession::Snapshot snap = p.ui.mark.session.snapshot();
    CHECK(contentHas(c, &snap.preview->meshes[0]) && !ownersHave(p.ui, goodScene));
    CHECK(c.caption.find("No stitch mesh: " + b->use.reason) != std::string::npos);
    std::printf("  caption: %s\n", c.caption.c_str());
    CHECK(!viewer::exportStitchMesh(p.ui, (std::filesystem::temp_directory_path() / "knit_ui_import" / "stitch_empty_rebuild" / "none.obj").string()));
    CHECK(p.ui.mark.note.find(b->use.reason) != std::string::npos);
    CHECK(p.ui.mark.baseline.get() != b);
    const std::string text = renderedText(p, W, H);
    CHECK(text.find("Export Stitch Mesh ...") == std::string::npos);
}

// Export Stitch Mesh in the view: the button asks for the save dialog (none without dialogs), its answer writes the file,
// whatever the checks found; the offered name names the file and the piece. (Defect: an export behind the Final knit or
// behind a passing fit.)
TEST(ui_mark_export_stitch_mesh) {
    const float W = 1280, H = 1000;
    UiContext ctx(W, H);
    UiFixture p("teddy");
    CHECK(openMarkTubes(p, "stitch_export"));
    CHECK(acceptedTube(p, 0));
    CHECK(buildArrived(p, pressBuild(p), 120));
    const knit::MarkedPieceView* b = viewer::shownBuild(p.ui);
    CHECK(b && b->use.drawable);
    if (!b) return;
    std::printf("  %s (clean: %s)\n", knit::stitchMeshHeadline(b->use, b->quality).c_str(), b->quality.clean ? "yes" : "no");
    CHECK(click(p, W, H, "Export Stitch Mesh ...") && p.ui.platform.request == viewer::PlatformRequest::None);
    p.ui.platform.dialogs = true;
    CHECK(click(p, W, H, "Export Stitch Mesh ...") && p.ui.platform.request == viewer::PlatformRequest::SaveStitchMesh);
    CHECK(viewer::stitchMeshFileName(p.ui) == "tubes-left_piece_1-stitch-mesh.obj");   // the piece label "left, piece 1" as one word
    const std::string path = (std::filesystem::temp_directory_path() / "knit_ui_import" / "stitch_export" / viewer::stitchMeshFileName(p.ui)).string();
    CHECK(viewer::exportStitchMesh(p.ui, path));
    knit::StitchMeshData back;
    std::string err;
    CHECK(knit::readPolygonObj(path, back, &err) && back.faces() == b->stitch.faces());
    CHECK(p.ui.mark.note == "Saved " + std::to_string(b->stitch.faces()) + " stitch faces to " + path);
    CHECK(p.jobs->received.empty() && !p.ui.results.knit);
}

// In the Stitch Mesh view a click on the isolated piece edits nothing (no point placed with a marking active); a click on
// another piece still isolates it. (Defect: inspecting the mesh moves a marking and makes the mesh outdated.)
TEST(ui_mark_stitch_view_click_selects_only) {
    UiContext ctx(1280, 1000);
    UiFixture p("teddy");
    CHECK(openMarkStep(p, "stitch_click"));
    viewer::isolatePiece(p.ui, 0);
    viewer::addMarking(p.ui, knit::MarkingMeaning::Start);
    viewer::setMarkView(p.ui, viewer::MarkView::StitchMesh);
    const viewer::ViewPixel viewport{0, 0, 800, 600};
    const viewer::Camera camera = frontCamera(p.ui);
    const std::optional<viewer::ViewPixel> onPiece = pixelOfPiece(p.ui, camera, viewport, 0);
    CHECK(onPiece.has_value());
    if (!onPiece) return;
    const size_t undo = p.ui.mark.undo.size();
    viewer::markClick(p.ui, camera, *onPiece, viewer::ClickIntent::Place);
    CHECK(p.ui.mark.set.markings[0].points.empty() && p.ui.mark.undo.size() == undo);
    viewer::setNeighbourLook(p.ui, viewer::NeighbourLook::Shown);
    if (const std::optional<viewer::ViewPixel> other = pixelOfPiece(p.ui, camera, viewport, 1)) {
        viewer::markClick(p.ui, camera, *other, viewer::ClickIntent::Place);
        CHECK(p.ui.mark.piece == 1);
    }
}
