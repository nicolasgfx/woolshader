// The Export job: checks the rule and what it writes from; exportScene does the rest.
#include "jobs/job_steps.h"
#include "jobs/results.h"
#include "model/log.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include "export/export.h"

namespace knit {

// --- Export job

// The export rule (exportEligibility, as the UI shows it) on the job's own snapshot: a blocked
// export is refused.
static bool exportAllowed(const ExportEligibility& el, std::string* err) {
    if (el.state != ExportEligibility::State::Blocked) return true;
    std::string m = "export refused:";
    for (const std::string& r : el.reasons) m += " " + r;
    if (err) *err = m;
    return false;
}

// What the export writes from: under Fur the fur result's base scene (nothing
// knitted), else the job's knit, which must be usable.
static bool exportInputsUsable(const ExportRequest& r, const JobControl& control, std::string* err) {
    if (knitWorkflow(*r.settings)) return usableKnit(r.knit.get(), *r.settings, control, err);
    if (!r.fibres || !r.fibres->base) {
        if (err) *err = "export refused: under Fur the export needs the fur grown on the mesh";
        return false;
    }
    return true;
}

// The scene the export writes: the hair result's base scene (knit skipped) or the knit.
static const Scene& exportedScene(const ExportRequest& r) { return knitWorkflow(*r.settings) ? *r.knit->scene : *r.fibres->base; }

// The Export job: the rule on the job's own snapshot, what it writes from, then exportScene, which
// moves the finished files into the export folder (publish.h, IntoFolder).
std::shared_ptr<BuiltScene> exportResults(const ExportRequest& r, const JobControl& control, std::string* err) {
    const SceneEntry* found = jobScene(*r.settings, BuiltScene::Kind::Export, err);
    if (!found) return nullptr;
    const SceneEntry& e = *found;
    const auto t0 = std::chrono::steady_clock::now();
    const SceneSettings& s = *r.settings;
    // The same rule as the UI (results.h), on the job's own snapshot.
    const FibreChoice choice = r.includeFibres ? FibreChoice::Included : FibreChoice::LeftOut;
    const ExportEligibility el = exportEligibility(s, r.knit.get(), r.fibres.get(), choice);
    if (!exportAllowed(el, err)) return nullptr;
    if (!exportInputsUsable(r, control, err)) return nullptr;
    const SceneFibres* fibres = r.includeFibres && r.fibres ? r.fibres->scene.get() : nullptr;
    if (r.includeFibres && !fibres) {
        if (err) *err = "export: fibres requested without a fibre result (growSceneFibres)";
        return nullptr;
    }
    std::shared_ptr<BuiltScene> out = newResultScene(BuiltScene::Kind::Export, control, s);
    out->knit = r.knit;
    out->fibreResult = r.fibres;
    ExportOptions o;
    o.format = r.format;
    o.folder = r.dir;
    o.name = r.name;
    o.policy = PublishPolicy::IntoFolder;
    o.status = el.state;
    o.stageProfile = true;  // per-part stage times in the export log
    o.onStage = [&](const std::string& stage, int done, int total) { reportProgress(control, stage, done, total); };
    std::string log;
    if (!exportScene({e, exportedScene(r), s, fibres}, o, &log, err)) return nullptr;
    out->exportedTo = std::filesystem::absolute(r.dir).string();
    const std::string diagnostic = el.state == ExportEligibility::State::Diagnostic ? "export status: diagnostic (parts failed validation)\n" : "";
    out->exportLog = log + diagnostic + findExportFormat(r.format)->label + " exported to " + out->exportedTo + "\n";
    out->ms = msSince(t0);
    return out;
}

bool exportPreviewPackage(const ExportRequest& r, const SceneEntry& e, const std::function<void(const std::string&)>& progress, std::string* err) {
    if (!exportInputsUsable(r, JobControl{}, err)) return false;
    const SceneFibres* fibres = r.includeFibres && r.fibres ? r.fibres->scene.get() : nullptr;
    ExportOptions o;
    o.format = kPackageFormat;
    o.folder = r.dir;
    o.name = r.name;
    o.policy = PublishPolicy::IntoFolder;
    o.status = ExportEligibility::State::Diagnostic;  // not a deliverable: a preview of whatever is there
    o.stageProfile = true;  // per-part stage times, in export_log.txt beside the package
    o.finishBakes = true;   // a first preview of new yarn maps bakes them once; a key change waits for that
    if (progress) o.onStage = [&](const std::string& stage, int done, int total) {
        progress(total > 0 ? stage + " (" + std::to_string(done) + " of " + std::to_string(total) + ")" : stage);
    };
    // What the preview shows follows the import's selection alone: a knitted part is its yarn and fibres
    // (never a filling behind them), a part left unknitted its own mesh.
    SceneSettings settings = *r.settings;
    settings.filling = false;
    std::string log;
    const bool written = exportScene({e, exportedScene(r), settings, fibres}, o, &log, err);
    std::ofstream(std::filesystem::path(r.dir).parent_path() / (std::filesystem::path(r.dir).filename().string() + "_export_log.txt")) << log;
    return written;
}

} // namespace knit
