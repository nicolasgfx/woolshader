#pragma once
// A run with nobody at the keyboard: the command line, the --do steps
// and when --shot is due. No GL and no GLFW, so all of it is tested without a window; main.cpp reads
// the pixels.
#include <string>
#include <vector>

#include "viewer/workflow.h"

namespace viewer {

struct ViewerOptions {
    std::string scene, settingsFile, stage, steps, lighting, shotPath;
    std::string importFile, units, up = "y";
    std::string loadKnit;        // --load-knit dir: the spiral pieces from a saved knit (knit_check --save-knit)
    ImportChoices importChoices;
    bool keepFileIfRebuildFails = false;  // --keep-file-if-rebuild-fails: press "Keep the file's own surface" after a failed rebuild
    CameraView startView = CameraView::Front;
    int width = 0, height = 0;   // --size (0: the saved or the default size)
    int plain = -1;              // --plain on|off: "Unknitted parts" (-1: as it is)
    int benchFrames = 0;         // --bench n (with --shot): time n frames at rest and n orbiting
    float uiScale = 0;           // --ui-scale (0: the saved one)
    double occlusion = -1;       // --occlusion k: the Studio occlusion strength (-1: the default)
    // --zoom k --target-frac fx,fy,fz (shots): the framed view aimed at that fraction of the scene's
    // box, its distance divided by k.
    double zoom = 0, targetFraction[3] = {0.5, 0.5, 0.5};
};

// A switch of the removed lab viewers, for main's refusal ("" none).
std::string retiredOption(int argc, char** argv);
ViewerOptions parseViewerOptions(int argc, char** argv);
double parseUnits(const std::string& unit);          // "mm", "cm", "m", "in", or mm per unit
Stage stageNamed(const std::string& name);           // anything else than import, fibres, export: Knit
void applyScriptStep(SceneUi& ui, const std::string& step);  // one --do step
// The settings: this machine's compute defaults for the scene (or no object), then --settings. False:
// refused (an unknown scene, a settings file not read), said on stderr.
[[nodiscard]] bool loadSettings(SceneUi& ui, const ViewerOptions& o);
// Opens what the command line or the settings name; a scripted import's unit is set before the first poll.
void openObject(SceneUi& ui, const ViewerOptions& o);

// With no --do, --shot or --units it does nothing.
class ScriptedRun {
public:
    explicit ScriptedRun(const ViewerOptions& options);
    // Once per frame before the interface: confirms a scripted import, runs the next step. False: a
    // scripted import failed (reported on stderr); nobody can answer its dialog, so the run ends.
    [[nodiscard]] bool beforeFrame(SceneUi& ui, double now);
    void frameShot(Camera& camera, const Bounds& content) const;  // --zoom, --target-frac: every frame
    bool shotDue();                                               // after the frame
    int reportShot(const SceneUi& ui, bool written) const;        // prints the outcome, returns the exit code

private:
    ViewerOptions options_;
    std::vector<std::string> todo_;
    bool importConfirmed_ = false, idle_ = false, shotNow_ = false;
    double runningSince_ = 0;
    double backgroundSince_ = 0;   // when the Mark step's deep work or a job began to run (mark-burst waits 0.2 s into it)
    int quietFrames_ = 0;
};

} // namespace viewer
