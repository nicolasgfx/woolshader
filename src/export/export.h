#pragma once
// The door to the export: which formats exist, the verdict a result is exported under, and exportScene,
// which writes one format into a folder that appears complete or not at all. The program includes
// this header only; ARCHITECTURE.md says what stands behind it.
#include <functional>
#include <string>
#include <vector>

#include "io/publish.h"

namespace knit {

struct Scene;
struct SceneEntry;
struct SceneFibres;
struct SceneSettings;

// --- the formats

enum class ExportFormatKind { Blend, Ppt, Obj };

inline constexpr const char* kPackageFormat = "ppt";

// Everything the program says about a format, in one row. A format is listed only when it is
// implemented; its description says what it holds and what it leaves out.
struct ExportFormat {
    ExportFormatKind kind = ExportFormatKind::Ppt;
    std::string id;                       // requests and scripts name a format by it
    std::string label;                    // in job messages and the log
    std::string name;                     // in the interface
    std::string shortName;                // "Export <shortName>", "Omitted by <shortName>"
    std::string contentsHeading;          // the disclosure over `contents`
    std::string contents;                 // what the folder holds
    std::string summary;                  // the contents in one line
    std::string meshSummary;              // the same with knitting skipped
    std::string omits;                    // what this format cannot carry ("": nothing)
    std::string coordinates;              // the sentence naming the files' coordinates ("": none)
    std::string importCoordinates;        // the same for an imported scene
    std::vector<std::string> extensions;  // of the files that take the export's name (none: names are fixed)
    bool fibres = false;                  // can carry the fibres
    bool available = true;
    std::string unavailable;              // why not, when it is not available
};

const std::vector<ExportFormat>& exportFormats();
const ExportFormat* findExportFormat(const std::string& id);

// The format a new Export step starts on: the Blender file, or the package where Blender is missing.
std::string defaultExportFormat();

// --- may it be exported? The rule over the reviewed results is the jobs' (exportEligibility,
// jobs/results.h); the export only writes its verdict into report.json.

struct ExportEligibility {
    enum class State { Ready, Diagnostic, Blocked } state = State::Blocked;
    std::vector<std::string> reasons;  // Blocked: why; Diagnostic: the failed checks
    std::string fix;                   // the next action (Blocked)
};
// A user-chosen name of the exported files: one plain path component (no separators, no "." / "..",
// no reserved Windows device name, no trailing dot or space, at most 120 characters).
bool validateExportName(const std::string& name, std::string* why);

// --- the renderer's GPU memory

// The most the renderer should take: the editor and the renderer share a 12 GB card (owner, 2026-11-04).
inline constexpr double kRendererMemoryCapGiB = 9.0;
// About how much GPU memory photon_path_tracer takes WHILE IT RENDERS a package of these counts, at the
// package's camera size (1600 x 1000) [GiB]: a fit to its render-phase process peak.
// It excludes the acceleration build, which peaks higher (one white bear at density 350: 9.8 against
// 7.6 GiB), so it is an estimate, not a capacity guarantee.
// textureTexels: 8-bit RGB texels of every map uploaded (the yarn maps), 2.7 B each, measured on the blue swatch's
// 24 x 4 atlases (884 MB).
double rendererMemoryGiB(long long strandSegments, long long triangles, long long textureTexels = 0);

// --- the yarn made of fibres (yarn_body.*, fibres/yarn_body.h)

// The yarn body setting a refusal concerns: the Knit step marks that control.
enum class YarnBodyField { None, Plies, Fibres, InteriorFibres, FibreDiameter, ShellDepth, FibreAngle, Irregularity, Core, MapTexel };

// The yarn body's settings checked at the yarn the export draws (the settings' rendered diameter).
struct YarnBodyCheck {
    bool accepted = true;
    YarnBodyField field = YarnBodyField::None;  // the control to change (None: a settings-file value)
    std::string reason;                         // the generator's words, naming the largest accepted count or depth
    std::string remedy;                         // what to change, in the Knit step's control names
};
YarnBodyCheck checkYarnBody(const SceneSettings& s);
// About how many fibre segments the yarn body writes for yarnLengthMm of yarn (a measured rate, for
// the forecasts; the export counts them exactly).
long long yarnBodySegments(const SceneSettings& s, double yarnLengthMm);
// The yarn maps (yarn_body.render = maps) checked at the longest yarn: its atlas must fit the renderer's texture
// side (yarn_bake.h planAtlas); refused with field MapTexel and the remedy. Accepted while the body is not rendered as maps.
YarnBodyCheck checkYarnMaps(const SceneSettings& s, double longestYarnMm);
// About how many texels the yarn maps of yarnLengthMm of yarn hold, over every map, and how many triangles
// their tubes have (the mapped tube's own tessellation, not the Export level's): the forecasts.
long long yarnMapTexels(const SceneSettings& s, double yarnLengthMm);
long long yarnMapTriangles(double yarnLengthMm);
// The settings that decide a yarn's maps (the body's structure and seed, the texel), as text: the viewport
// compares the shown maps' settings with the current ones by it.
std::string yarnMapInputsText(const SceneSettings& s);

// --- the export

// What is exported: reviewed results, never computed again.
struct ExportSource {
    const SceneEntry& entry;
    const Scene& scene;
    const SceneSettings& settings;
    const SceneFibres* fibres = nullptr;  // null: no fibres; else written as they are, below their key threshold
};

struct ExportOptions {
    std::string format = kPackageFormat;  // an ExportFormat::id
    std::string folder;
    std::string name;                     // stem of the files; "": the folder's name, else "knit"
    PublishPolicy policy = PublishPolicy::Strict;
    ExportEligibility::State status = ExportEligibility::State::Ready;  // Diagnostic: said in report.json
    int closeupPart = -1;                 // the part of the close-up cameras (-1: the largest knitted part; closeupPartOf)
    bool stageProfile = false;            // stage times in the log
    // Developer (knit_check --taper-points, --taper-tolerance): where the late taper's radii are
    // carried, a TaperPoints name (fibres/fibre_strands.h); "" is "dense", the export's own.
    std::string taperPoints;
    double taperTolerance = 0.05;
    // Developer (knit_check --body-file-segments, the tests): the most segments per yarn-body strand file.
    long long bodyFileSegments = 2'500'000;  // ExportScene::kBodyFileSegments
    // The stage about to run, for a progress line; done of total when the stage counts (a bake's yarns), else 0, 0.
    std::function<void(const std::string& stage, int done, int total)> onStage;
    // The render preview's exporter: a bake it starts finishes before a cancel acts (baked_maps.h BakeInterrupt),
    // the cache-bound minutes are never thrown away by a package key change.
    bool finishBakes = false;
};

// False with the reason; throws Cancelled at its checkpoints. A failed or cancelled export leaves the
// folder as it was. Eligibility is the caller's (jobs/results.h): the command line and the tests have no KnitResult.
[[nodiscard]] bool exportScene(const ExportSource& source, const ExportOptions& options, std::string* log,
                               std::string* err);

} // namespace knit
