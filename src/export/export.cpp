// The door: the format table, the export rule and exportScene, which runs an export from the first
// refusal to the published folder.
#include "export/export.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>

#include "export/collect_scene.h"
#include "export/export_scene.h"
#include "fibres/yarn_bake.h"
#include "fibres/yarn_body.h"
#include "io/files.h"
#include "model/log.h"
#include "scene/knit/knit_report.h"
#include "scene/result_keys.h"
#include "scene/run_log.h"
#include "scene/scene.h"
#include "scene/stage_params.h"
#include "settings/settings_mapping.h"
#include "settings/yarn_look.h"

namespace knit {

namespace fs = std::filesystem;

// --- the formats

const std::vector<ExportFormat>& exportFormats() {
    static const std::vector<ExportFormat> f = [] {
        std::vector<ExportFormat> v;
        ExportFormat b;
        b.kind = ExportFormatKind::Blend;
        b.id = "blend";
        b.label = "Blender file";
        b.name = "Blender file";
        b.shortName = "Blender file";
        b.contentsHeading = "What the files hold";
        b.contents = "One .blend for Blender 5.1 in metres, Z up: yarn tubes and an import's unknitted parts as meshes with their normals and "
                     "uvs, fibres and fur as native curves with their radii, materials with the textures packed in, the cameras, the ground and "
                     "the sky, and Cycles set up like the package's render (samples, bounces, tone curve); report.json with the settings and every "
                     "part's verdict.";
        b.omits = "The renderer's own sheen lobe (Cycles has none: yarn and wool fibres use a shader fitted to its look), and the ply maps no material "
                  "uses (height, AO)."
                  " A yarn of fibres rendered as maps keeps its cavity, normal, roughness and specular maps on a Principled BSDF; the highlight's "
                  "anisotropy is approximated (Cycles stretches highlights differently) and the height map goes unused.";
        b.summary = "One .blend with every part, fibres as curves, the camera and report.json.";
        b.meshSummary = "One .blend with the mesh, its fur as curves and the camera, and report.json.";
        b.coordinates = "Coordinates: metres, Z up.";
        b.importCoordinates = b.coordinates;
        b.extensions = {".blend"};
        b.fibres = true;
        b.available = blenderConversionAvailable();
        b.unavailable = "Needs Blender 5.1 or later and blender/ppt_package_to_blend.py. KNIT_BLENDER names the blender.exe to use.";
        v.push_back(b);
        ExportFormat p;
        p.kind = ExportFormatKind::Ppt;
        p.id = kPackageFormat;
        p.label = "photon_path_tracer package";
        p.name = "Photon Path Tracer package";
        p.shortName = "package";
        p.contentsHeading = "What the package holds";
        p.contents = "scene.json with materials, a camera, a grey ground plane and a uniform sky; yarn tubes and an import's unknitted "
                     "parts as PLY meshes (with their base colour texture); fibres as .pptstrands; ply shading maps (normal, height); report.json with the settings and every part's verdict.";
        p.summary = "Meshes, optional fibres, ply maps and report.";
        p.meshSummary = "The mesh with its fur, and report.json.";
        p.fibres = true;
        v.push_back(p);
        ExportFormat o;
        o.kind = ExportFormatKind::Obj;
        o.id = "obj";
        o.label = "Wavefront OBJ";
        o.name = "Wavefront OBJ";
        o.shortName = "OBJ";
        o.contentsHeading = "What the files hold";
        o.contents = "One .obj with every yarn tube and an import's unknitted parts as named objects, fibres and fur as "
                     "polylines (OBJ line elements, one object per fibre group), and an .mtl with one base colour per material "
                     "(unknitted parts: their base colour texture, beside the .obj); report.json with the settings and every part's "
                     "verdict.";
        o.omits = "The fibres' radius (OBJ lines have none: each fibre object's comment gives its mean radius), the ply shading maps "
                  "(they need the package's shading model) and the ground. Sheen, coat and fibre shading are reduced to the base colour."
                  " A yarn of fibres rendered as maps: its specular and roughness maps (the cavity map goes as map_Kd, the normal map as map_bump).";
        o.summary = "One .obj and .mtl with every part, fibres as lines, and report.json.";
        o.meshSummary = "One .obj and .mtl with the mesh and its fur (as lines), and report.json.";
        o.coordinates = "Coordinates: metres, Y up.";
        o.importCoordinates = "Coordinates: the imported file's units and up axis.";
        o.extensions = {".obj", ".mtl"};
        o.fibres = true;
        v.push_back(o);
        return v;
    }();
    return f;
}

const ExportFormat* findExportFormat(const std::string& id) {
    for (const ExportFormat& f : exportFormats())
        if (f.id == id) return &f;
    return nullptr;
}

std::string defaultExportFormat() { return findExportFormat("blend")->available ? "blend" : kPackageFormat; }

bool validateExportName(const std::string& name, std::string* why) {
    auto no = [&](const std::string& w) {
        if (why) *why = w;
        return false;
    };
    if (name.empty()) return no("enter a file name");
    if (name.size() > 120) return no("the file name is longer than 120 characters");
    if (name == "." || name == "..") return no("\"" + name + "\" is not a file name");
    for (unsigned char c : name) {
        if (c < 32) return no("the file name contains a control character");
        if (std::string("<>:\"/\\|?*").find(char(c)) != std::string::npos)
            return no(std::string("the file name may not contain ") + char(c) + " (a name, not a path)");
    }
    if (name.back() == '.' || name.back() == ' ') return no("the file name may not end with a dot or a space");
    // Reserved device names (with or without an extension).
    std::string stem = name.substr(0, name.find('.'));
    while (!stem.empty() && stem.back() == ' ') stem.pop_back();
    std::string up;
    for (char c : stem) up += char(std::toupper(static_cast<unsigned char>(c)));
    static const char* reserved[] = {"CON", "PRN", "AUX", "NUL", "COM1", "COM2", "COM3", "COM4", "COM5", "COM6", "COM7", "COM8", "COM9",
                                     "LPT1", "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8", "LPT9"};
    for (const char* r : reserved)
        if (up == r) return no("\"" + stem + "\" is a reserved device name on Windows");
    return true;
}

// --- the export

namespace {

// Why this export cannot start, "" if it can: refused before anything is written.
std::string exportRefusal(const ExportSource& source, const ExportOptions& options) {
    const Scene& scene = source.scene;
    if (options.closeupPart >= static_cast<int>(scene.parts.size()))
        return "close-up of part " + std::to_string(options.closeupPart) + ": the scene has " +
               std::to_string(scene.parts.size()) + " parts";
    return "";
}

// report.json: the settings and the per-part verdicts. The same from the viewer and knit_check.
bool writeReport(const ExportSource& source, ExportEligibility::State status, const std::string& folder,
                 std::string* err) {
    const SceneEntry& entry = source.entry;
    const Scene& scene = source.scene;
    std::ofstream f(folder + "/report.json");
    f << "{\n  \"generator\": \"wool_shader scene " << entry.id << "\",\n";
    if (status == ExportEligibility::State::Diagnostic) f << "  \"export_status\": \"diagnostic\",\n";
    if (entry.isPreset()) f << "  \"credit\": \"" << jsonEscaped(creditText(entry.credit, workflowOf(source.settings))) << "\",\n";
    f << "  \"settings\": \"" << jsonEscaped(writeSettings(source.settings));
    f << "\",\n  \"cross_component\": \"" << toString(scene.cross.status) << "\",\n  \"parts\": [\n";
    for (size_t i = 0; i < scene.parts.size(); ++i) {
        const ScenePart& r = scene.parts[i];
        f << "    {\"name\": \"" << r.name << "\", \"knitted\": " << (r.knit ? "true" : "false");
        if (r.knit)
            f << ", \"verdict\": \"" << toString(r.knit->verdict) << "\", \"gate_a\": \"" << toString(r.knit->after.verdict)
              << "\", \"stitches\": " << r.knit->t.prog.stitches.size() << ", \"yarn_diameter_mm\": " << r.yarnD;
        if (r.guided) {  // the guided route's provenance and outcome (only on that route: other reports stay as they were)
            const GuidedPart& g = *r.guided;
            f << ", \"route\": \"guided stitch mesh (experimental)\", \"guided\": {\"outcome\": \"" << toString(g.outcome) << "\", \"reason\": \""
              << jsonEscaped(g.reason) << "\", \"target_diameter_mm\": " << g.targetD << ", \"reached_diameter_mm\": " << g.reachedD << ", \"calls\": " << g.calls
              << ", \"calls_certified\": " << g.callsCertified << ", \"start_verified\": " << (g.startVerified ? "true" : "false") << ", \"lines\": [";
            const std::vector<std::string> lines = guidedPartLines(g);
            for (size_t k = 0; k < lines.size(); ++k) f << (k ? ", " : "") << "\"" << jsonEscaped(lines[k]) << "\"";
            f << "]}";
        }
        f << "}" << (i + 1 < scene.parts.size() ? "," : "") << "\n";
    }
    f << "  ]\n}\n";
    return f ? true : fail(err, "cannot write report.json");
}

// The verdicts and a preset's credit, which go with every format. Written beside the parts too,
// because Blender reads the credit from there.
bool writeSidecars(const ExportSource& source, ExportEligibility::State status, const std::string& folder,
                   std::string* err) {
    if (!writeReport(source, status, folder, err)) return false;
    if (!source.entry.isPreset()) return true;
    std::ofstream credits(folder + "/CREDITS.txt");
    credits << creditText(source.entry.credit, workflowOf(source.settings));
    return credits ? true : fail(err, "cannot write CREDITS.txt");
}

// One case per format. MSVC does not warn about a missing case: export_every_format_writes does.
bool writeFormat(const ExportFormat& format, const ExportScene& scene, const FormatTarget& target, std::string& log,
                 std::string* err) {
    switch (format.kind) {
    case ExportFormatKind::Ppt: return writeSceneJson(scene, err);
    case ExportFormatKind::Obj: return writeObjFile(scene, target, log, err);
    case ExportFormatKind::Blend: return writeSceneJson(scene, err) && writeBlendFile(scene.folder(), target, log, err);
    }
    return fail(err, "no writer for " + format.label);
}

std::string baseNameOf(const ExportOptions& options) {
    if (!options.name.empty()) return options.name;
    const std::string folderName = fs::path(options.folder).filename().string();
    return folderName.empty() ? "knit" : folderName;
}

} // namespace

double rendererMemoryGiB(long long strandSegments, long long triangles, long long textureTexels) {
    // photon_path_tracer 70ab4eeba's process peak on a white toy of a pieces scene, 2026-11-05, at
    // 1600 x 1000: 1512 MiB for 4.02 M triangles without strands, 1758 MiB for 7.89 M; 7832 MiB with
    // 33.97 M strand segments, 10636 MiB with 49.01 M. Within 3 MiB of all four. Textures: the blue
    // swatch's yarn maps, 2026-11-11, 24 x 4 maps of 4096 x 820 texels took 884 MB of texture memory.
    constexpr double kFixedMiB = 1256, kMiBPerMillionSegments = 186.1, kMiBPerMillionTriangles = 63.7, kBytesPerTexel = 2.7;
    const double mib = kFixedMiB + kMiBPerMillionSegments * double(strandSegments) / 1e6 + kMiBPerMillionTriangles * double(triangles) / 1e6 +
                       kBytesPerTexel * double(textureTexels) / (1024.0 * 1024.0);
    return mib / 1024;
}

YarnBodyCheck checkYarnBody(const SceneSettings& s) {
    YarnBodyCheck c;
    YarnBodyParams p;
    std::string why;
    if (!yarnBodyParams(s, p, &why)) {
        c.accepted = false;
        c.reason = why;
        c.remedy = "Correct the value in the settings file.";
        return c;
    }
    const double drawn = finalSurface(s, knitDiameter(s, s.diameter)).renderedDiameter();
    const BodyRefusal refusal = yarnBodyRefusal(p, 0.5 * drawn, &why);
    if (refusal == BodyRefusal::None) return c;
    c.accepted = false;
    const std::string prefix = "yarn body: ";
    c.reason = why.rfind(prefix, 0) == 0 ? why.substr(prefix.size()) : why;
    switch (refusal) {
    case BodyRefusal::None: break;
    case BodyRefusal::Plies: c.field = YarnBodyField::Plies; c.remedy = "Choose 1 to 3 plies."; break;
    case BodyRefusal::FibreCount: c.field = YarnBodyField::Fibres; c.remedy = "Choose at least one fibre per yarn."; break;
    case BodyRefusal::FibreTooThin: c.field = YarnBodyField::FibreDiameter; c.remedy = "Choose a fibre diameter of at least 10 um."; break;
    case BodyRefusal::FibreTooLarge: c.field = YarnBodyField::FibreDiameter; c.remedy = "Choose thinner fibres or fewer plies."; break;
    case BodyRefusal::ShellDepth: c.field = YarnBodyField::ShellDepth; c.remedy = "Choose a shell depth above 0 % and at most 100 %."; break;
    case BodyRefusal::Angles: c.field = YarnBodyField::FibreAngle; c.remedy = "Choose fibre and ply angles of 0 to 60 deg."; break;
    case BodyRefusal::Irregularity: c.field = YarnBodyField::Irregularity; c.remedy = "Choose an irregularity of 0 to 2."; break;
    case BodyRefusal::ShellFull:
        c.field = YarnBodyField::Fibres;
        c.remedy = "Lower Fibres per yarn, increase Shell depth or choose thinner fibres.";
        break;
    case BodyRefusal::InteriorFull:
        c.field = YarnBodyField::InteriorFibres;
        c.remedy = "Lower Interior fibres, or reduce Shell depth below 100 % to leave room inside.";
        break;
    case BodyRefusal::NoCoreRoom: c.field = YarnBodyField::Core; c.remedy = "Reduce Shell depth or choose Core: None."; break;
    }
    return c;
}

YarnBodyCheck checkYarnMaps(const SceneSettings& s, double longestYarnMm) {
    YarnBodyCheck c;
    if (yarnBodyRenderOf(s) != YarnBodyRender::Maps || longestYarnMm <= 0) return c;
    const double drawn = finalSurface(s, knitDiameter(s, s.diameter)).renderedDiameter();
    BakedYarnMaps plan;
    std::string why;
    if (planAtlas(longestYarnMm, 0.5 * drawn, bakeParams(s), plan, &why)) return c;
    c.accepted = false;
    c.field = YarnBodyField::MapTexel;
    const std::string prefix = "yarn bake: ";
    c.reason = why.rfind(prefix, 0) == 0 ? why.substr(prefix.size()) : why;
    c.remedy = "Choose a coarser Map texel.";
    return c;
}

long long yarnMapTexels(const SceneSettings& s, double yarnLengthMm) {
    // Every yarn's atlas is its length along by its circumference around, in texels, four maps each (cavity,
    // normal, specular, roughness); the unused tails of the strips are not counted.
    constexpr int kMapsPerYarn = 4;
    const double texel = bakeParams(s).texel;
    const double r = 0.5 * finalSurface(s, knitDiameter(s, s.diameter)).renderedDiameter();
    if (!(texel > 0) || !(r > 0)) return 0;
    const double along = yarnLengthMm / texel, around = std::max(8.0, std::round(2 * kPi * r / texel));
    return std::llround(kMapsPerYarn * along * around);
}

long long yarnMapTriangles(double yarnLengthMm) { return std::llround(2.0 * kBakeSides * yarnLengthMm / kBakeRing); }

std::string yarnMapInputsText(const SceneSettings& s) {
    const SceneSettings::YarnBody& b = s.yarnBody;
    char t[256];
    std::snprintf(t, sizeof t, "plies %d fibres %d interior %d diameter %.17g shell %.17g angles %.17g %.17g irregularity %.17g seed %d texel %.17g",
                  b.plies, b.fibres, b.interiorFibres, b.fibreDiameterUm, b.shellDepth, b.fibreAngleDeg, b.plyAngleDeg, b.irregularity, b.seed,
                  b.mapTexelUm);
    return t;
}

long long yarnBodySegments(const SceneSettings& s, double yarnLengthMm) {
    // Measured on the swatch (2026-11-07): 400 fibres of 36 um on 12.35 m of
    // yarn wrote 35.7 M segments, 7.2 per fibre and mm of yarn; at 24 um about 7 (the 6 degree turn
    // target sets the step, so the rate hardly depends on the fibre's size).
    constexpr double kSegmentsPerFibreMm = 7.2;
    const double fibres = double(s.yarnBody.fibres) + double(s.yarnBody.interiorFibres);
    return std::llround(fibres * yarnLengthMm * kSegmentsPerFibreMm);
}

bool exportScene(const ExportSource& source, const ExportOptions& options, std::string* log, std::string* err) {
    const ExportFormat* format = findExportFormat(options.format);
    if (!format || !format->available) return fail(err, "unknown or unavailable export format " + options.format);
    if (const std::string why = exportRefusal(source, options); !why.empty()) return fail(err, why);
    const auto stage = [&](const std::string& text) {
        if (options.onStage) options.onStage(text, 0, 0);
    };
    std::string ownLog;
    std::string& logText = log ? *log : ownLog;
    stage("preparing the export folder");
    recoverPublished(options.folder);
    Publisher publisher;
    if (!publisher.begin(options.folder, err, options.policy)) return false;
    const std::string staging = publisher.staging();
    // The parts are the package itself; for another format, a working folder removed at the end.
    const std::string partsFolder = format->kind == ExportFormatKind::Ppt ? staging : staging + "/.package";
    std::error_code ec;
    fs::create_directories(partsFolder + "/meshes", ec);

    stage("writing the " + format->label);
    ExportScene scene(partsFolder);
    if (!collectScene(source, options, scene, logText, err)) return false;
    const double gib = rendererMemoryGiB(static_cast<long long>(scene.strandSegmentCount()), static_cast<long long>(scene.triangleCount()));
    appendFormat(logText, "renderer GPU memory, render phase: about %.1f GiB (%.1f M strand segments, %.1f M triangles); "
                 "excludes the acceleration build's peak, which is higher\n", gib,
                 double(scene.strandSegmentCount()) / 1e6, double(scene.triangleCount()) / 1e6);
    if (gib > kRendererMemoryCapGiB)
        appendFormat(logText, "warning: above the %.0f GiB the editor and the renderer can share on a 12 GB card: %s\n", kRendererMemoryCapGiB,
                     source.settings.yarnBody.enabled ? "lower the yarn's Fibres per yarn (Knit step), or the fuzz's density"
                                                      : "choose a lighter Export level or a lower fibre density");
    const StageClock clock;
    if (!frameScene(scene, err)) return false;
    if (!writeSidecars(source, options.status, partsFolder, err)) return false;
    if (partsFolder != staging) stage("converting to " + format->label);
    if (!writeFormat(*format, scene, {staging, baseNameOf(options)}, logText, err)) return false;
    if (options.stageProfile) appendFormat(logText, "profile: %s\n", clock.line("export: the format's files").c_str());
    if (partsFolder != staging) {
        if (!writeSidecars(source, options.status, staging, err)) return false;
        fs::remove_all(partsFolder, ec);
    }

    stage("publishing");
    const uint64_t key = stageKey(source.settings, Stage::Export, source.scene.source.fileHash);
    if (!publisher.commit("scene " + source.entry.id + " " + format->id + " " + std::to_string(key), err)) return false;
    appendFormat(logText, "export: %zu meshes, %zu strand parts -> %s (%s)\n", scene.meshParts().size(),
                 scene.strandParts().size(), options.folder.c_str(), format->label.c_str());
    return true;
}

} // namespace knit
