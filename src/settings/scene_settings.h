#pragma once
// The scene settings shared by knit_viewer and knit_check --scene. One serializable struct; a table
// drives write, read and diff and the stage keys (result_keys.h), so every key round-trips and every
// key belongs to exactly one build stage.
//
// Text format ("knit-settings v1"): one `key = value` per line, dotted keys, '#'
// comments. Doubles are written with 17 significant digits (exact round trip). An
// unknown key or a missing / different version line is an error, never ignored.
#include <cstdint>
#include <string>
#include <vector>

#include "model/choices.h"
#include "model/vec.h"
#include "settings/material_spec.h"

namespace knit {

// Which build stage a key invalidates (a stage is rebuilt iff its key, which
// includes every upstream key and the source asset hash, changes).
enum class Stage : uint8_t { Load, Knit, Surface, Fibres, Export, Display };
const char* toString(Stage s);

// What is made of the model (key workflow). Knit: Import -> Knit -> Fibres -> Export, wool fibres on the
// knitted yarn. Fur: Import -> Fur -> Export, a loop pile and fine fibres grown on the model's own mesh
// (fibres/fur_field.h; no knit). Each workflow has its own keys; switching keeps both.
enum class Workflow : uint8_t { Knit, Fur };
const char* toString(Workflow w);

struct SceneSettings {
    std::string scene = "wool-swatch";  // registry id (scene_registry.h), or "import" (the import.* keys)

    // ---- an imported mesh (Load stage; scene = "import"; prepare_import.h, import_store.h). The identity
    // names the prepared cache entry (out/cache/imports/<identity>): it covers the source
    // file, its MTL / images, these options and the importer version. Inputs of
    // scene = "import" only.
    std::string importSource;      // the file as the user chose it
    std::string importIdentity;    // 16 hex digits
    double importUnitMm = 1;       // mm per file unit
    std::string importUp = "y";    // y | z
    std::string importRepairs = "weld,degenerate,duplicates,unreferenced,orientation";
    std::string importSurface = "file";  // file (not rebuilt: Fur, or scripted runs on the file's uv) | wrapped
    double importWrapStitchMm = 5.5;     // the stitch width a wrapped surface is prepared for
    std::string importPlainNodes;        // the nodes kept as plain meshes, | between the names ("" none: all knitted)

    // ---- the workflow (Load stage): knit | fur (Workflow). It is no line of a stage's text: the Knit
    // stage keeps the line knit.enabled=true it always held (result_keys.cpp), so saved knits keep their keys.
    std::string workflow = "knit";
    // ---- the guided route (experimental; scene = import, Knit stage): one imported piece knitted as a stockinette stitch
    // mesh from its Mark step markings (scene/knit/guided_part.cpp) instead of the surface tiling; the file's other pieces
    // stay plain meshes. The digest is the effective markings' (guidedMarkingsDigest, scene/source/guided_source.h): the Knit
    // stage refuses a markings file that no longer matches it, and an edit of the markings changes it, so the knit goes stale.
    std::string knitRoute = "surface";  // surface (the default: tiling or sewn pieces) | guided
    int guidedPiece = -1;               // the piece of import.source (SurfacePieces index: the Mark step's numbering)
    int guidedKnitComponent = 0;        // the piece's knitting component knitted (SurfaceMarking::knitComponent; 0: the first)
    std::string guidedMarkings;         // the markings file (knit-markings v5)
    std::string guidedMarkingsDigest;   // 16 hex digits

    // ---- strand (Knit stage unless noted)
    std::string thicknessPreset = "thick-dense";  // thick-dense (Dense) | medium | thin | custom (thickness_presets.h)
    double diameter = 1.16;        // D [mm]
    double clearance = 0.1;        // c [mm]
    // Visible coverage (tiled scenes): "component" (the default: one part per component, charts culled at their boundaries) or "overlap" (one part
    // per chart, continued knit.overlap_stitches stitches / rows beyond the chart on the surface so
    // the seams are covered by overlapping, partly hidden yarn; knit.gap_patches: local knitted
    // patches, independent parts, over the uncovered regions left). Patches are validated on their
    // own: no clearance between patches, no connection of their yarns. The viewer and the CLI's
    // scene defaults (applySceneDefaults) choose "overlap" for every uv-lattice scene: with
    // "component" every chart border is an open seam (user, 2026-10-29).
    std::string coverage = "component";
    int overlapStitches = 2;
    bool gapPatches = true;
    double depthAmplitude = 0.8;   // A [D]
    double wale = 5.5, course = 3.4;  // gauge w, h [mm]
    // Stitch size: "custom" (w and h as set) or "yarn": w and h follow from D, c and the tightness
    // (applyStitchSize). Both are labels of the derived w and h (never hashed). w x h is the
    // stockinette stitch; rib and garter knit their relaxed size of it (patternStitch).
    std::string stitchSize = "custom";  // custom | yarn
    std::string tightness = "standard"; // tight | standard | loose
    double angleDeg = 0;           // orientation on the mesh (tiled scenes)
    std::string pattern = "stockinette";
    std::vector<double> diameterFallback{1.12, 1.08, 1.04, 1.0};
    bool settle = true, repair = true;
    HandKnit hand;                 // uneven knitting (hand_tension.h)
    // Surface stage:
    int plies = 1;                 // substrands (1: single tube; > 1: explicit plies)
    double plyPitch = 2.5, plyFill = 0.95;
    int plyHandedness = 1;         // +1 Z, -1 S
    // Ply map set (plies == 1 only): "none" or "<n>-swept" / "<n>-circles", n >= 2. The default
    // is the default look: 3 plies, swept model.
    std::string textureProfile = "3-swept";
    bool flatten = false;          // contact flattening (FinalSurfaceSpec)
    double flattenInflate = 0.15;  // e: rendered D_r = D (1 + e)
    double maxFlattening = 0.4;    // the most a flattened yarn is squashed: b >= (1 - maxFlattening) r (0.05 to 0.6)
    // A gauge too tight for a round yarn of diameter D (yarnFit): knit at the contact diameter
    // k D that fits and draw the yarn at D (contact flattening). Off: knit at D (expect FAIL).
    bool compressToFit = true;
    bool surfaceGate = true, surfaceHalf = true;
    // The strand's polygonal tube (Fibres stage): sides around and ring spacing along [mm] of the
    // single yarn tube and of explicit plies. The export writes it, the fibres' roots and checks are
    // measured against it, the viewer draws it. Never a knitting input (changing it never re-knits).
    // The default is the "Reduced" level (presets/strand_mesh/reduced.preset; 2026-10-29).
    int tubeSides = 6;
    double tubeRing = 0.6;
    int plyTubeSides = 4;
    double plyTubeRing = 0.45;
    // Explicit computation budget for the diameter search: at most this many diameters (the
    // requested one first, then strand.diameter_fallback) are knitted per piece; 0 = all (full search,
    // the default). A search the budget stops is reported with the
    // diameters not attempted; it does not claim that they would fail.
    int diameterAttempts = 0;

    // ---- fibres (Fibres stage)
    bool fibres = false;
    double fibreDensity = 125;     // candidates per mm of yarn
    FibreIslands islands;          // patches of different fuzz (fibre_islands.h)
    double randomFraction = 0.5;
    double curvatureBias = 0.0;    // FibreParams::curvatureWeight
    // Wool preset (fibre_presets.h): a label of the values (never part of a key); "" or "custom".
    std::string fibrePreset;
    int fibreMaxSegments = 3;      // FibreParams::maxSegments (2..48)
    double fibreSegmentLength = 0.45;  // FibreParams::segLenTarget [mm]: the length of one straight piece
    long long segmentBudget = 7000000;  // Export "Reduced" (presentation.cpp kExportDetails)
    // fibres.tail_taper*: the late taper of the exported strands (StrandTaper, fibre_strands.h).
    // tail 0: none, every fibre tapers along its whole length (Fibre::vertexRadius).
    double tailTaper = 0, tailTaperMinDiameters = 8, tailTaperTip = 0.1;
    // generation: "budget" (only the
    // candidates below the global strand-budget key threshold are evaluated; the export is
    // byte-identical to "full", per-piece diagnostics cover the evaluated subset) or "full"
    // (every candidate).
    std::string fibreGeneration = "budget";
    // fibres.checks: "light" (stay out of the own strand, at most 20 % of the radius into another one;
    // GPU) or "full" (every clearance, approach and exposure test).
    std::string fibreChecks = "light";
    double sideBandCos = -0.1;     // FibreParams::frontCos (pieceFibres' value)
    // fibres.crimp_*: the directed fibres' crimp (FibreParams::crimpAmplitude, crimpWavelength) [mm]; 0: none.
    // A crimp needs chords well below its wavelength (fibres.segment_length).
    double crimpAmplitude = 0, crimpWavelength = 1.5;
    // strand.ply_*_jitter (TubeParams::uvPhaseJitter, uvPitchJitter) and strand.ply_crown_tilt_deg
    // (PlyProfileParams::blendTiltDeg: a lower tilt starts the rounded bridge earlier, a rounder crown).
    double plyPhaseJitter = 0, plyPitchJitter = 0, plyCrownTiltDeg = 60;
    // strand.normal_strength: the ply map's normalmap_scale in the export (the relief's depth; 1 = as built).
    double normalStrength = 1;
    struct Population {
        double lengthMin = 0.4, lengthMax = 1.6, lengthMedian = 0.8, lengthSpread = 0.35;  // [mm]
        // Thickness (diameter) [mm]: authoritative median t and sigma_ln (FibreParams); the
        // log-normal factor is clamped at +-2 sigma (fibre_shape.cpp), so the range is exactly
        // [t e^-2s, t e^2s] (thicknessMin/Max); inside it the distribution is the clamped
        // log-normal, not uniform.
        double thickness = 0.07, thicknessSigma = 0.04;
        double thicknessMin() const;
        double thicknessMax() const;
        void setThicknessRange(double mn, double mx);  // t = sqrt(mn mx), s = ln(mx / mn) / 4
        double curvature = 1.0;    // directed: bow scale; random: randomBend
        double maxAngleDeg = 24;   // directed: max nominal lift; random: max start elevation
        // Random only (fibres.random.curl_*): the helix of a curled stray; radius 0: none.
        double curlRadius = 0, curlLength = 8;  // [mm], [mm of fibre per turn]
        bool operator==(const Population&) const = default;
    } directed, random;
    // Directed only, not preset-owned: the rise levels off this many yarn radii above the yarn
    // (0: linear, FibreParams::hugHeight), and the fibres wind with the ply helix (strand.ply_*).
    double directedHugHeight = 0;
    bool directedFollowPlies = false;
    // Directed only, preset-owned (fibres.directed.drape, FibreParams::drape): a fuzz fibre that would sink into a
    // neighbouring stitch rests on it instead of being cut there (Mohair). Light checks only.
    bool directedDrape = false;
    // Random only, owned by the stray presets (fibres.random.spread_deg, FibreParams::randomSpreadDeg):
    // seen from above the yarn's surface, a stray starts at most this many degrees away from the yarn,
    // either way along it. 180: any direction around the yarn. 45: the Loose stray preset.
    double randomSpreadDeg = 45;

    // ---- the Fur workflow (Fibres stage; pile.*, inputs only under workflow = fur; fibres/fur_look.h has
    // the meanings). The names follow the Redshift "RS Procedural Fur" panel the owner works from; the
    // semantics are ours. The defaults are the approved look (frozen baseline gate/top15, 2026-11-19).
    struct Pile {
        // Coverage: where the pile grows. entire = the whole mesh; mask = maskPath (a greyscale image in
        // colour-texture space, white = full density). A mask is used only when chosen here.
        std::string coverage = "entire";
        std::string maskPath;
        std::string maskHash;          // the mask file's content hash (16 hex; Reload re-hashes it)
        std::string excludedParts;     // part names without fur, comma separated
        double density = 400;          // loop fibres per cm^2
        double length = 1;             // nominal loop lengths ("turns") of each fibre's own loop
        double straightLength = 4;     // [mm] a fibre's length at curl 0
        double thickness = 0.6;        // root diameter [mm]
        double taper = 0.8;            // tip diameter / root diameter
        double bend = 1;               // the curvature reaches 1 / curl after this many quarter turns' arc
        double curl = 0.9;             // loop radius [mm]; 0: straight fibres
        double gravity = 9.81;         // [m/s^2]
        double gravityFullLength = 30; // [mm] shorter fibres take less gravity, as (L / it)^3
        // Variation (shares in [0, 1]; density scale in mm)
        double varDensity = 0, varDensityScale = 10, varLength = 0.5, varThickness = 0.33, varDirection = 1, varBend = 0.2, varCurl = 1;
        // The fine fibres over the loops (straight; a log-normal length clamped to [median / 2, 2 median])
        bool fibresEnabled = true;
        double fibresDensity = 150;    // per cm^2
        double fibresLength = 4;       // [mm] median
        double fibresLengthSpread = 0.3;  // sigma of ln(length)
        double fibresThickness = 0.06;    // root diameter [mm]
        double fibresTaper = 0.3;
        double fibresLiftDeg = 85;
        double fibresDirection = 1;    // lift +- 40 v deg, comb turn +- 180 v deg
        // Islands (per cm^2, diameter mm, strength 0..1, relative kind shares) and the irregularities (clump reach mm,
        // strength, patch amplitudes on length and lift, wear). Off at 0.
        double islandsRate = 0, islandsSize = 15, islandsStrength = 0.6, islandsMatted = 1, islandsFluffy = 1, islandsThin = 1;
        double clumpSize = 2, clumpStrength = 0, patchLength = 0, patchLift = 0, wear = 0;
        int seed = 0;
        bool operator==(const Pile&) const = default;
    } pile;
    MaterialSpec pileMaterial;           // material.pile.* (Export): the hair material
    // object (one colour per part: its source's stored base colour) | single (material.pile.color) | texture
    // (every fibre the texture colour at its root; a part without a texture: as object)
    std::string pileColourSource = "object";

    // ---- materials (Export stage; Display previews colours)
    MaterialSpec yarnWhite, yarnBrown, fibreDirected, fibreRandom, stuffing;
    bool fibreInheritsYarn = true; // directed fibres use their yarn's colour
    bool fibreRandomInheritsYarn = true; // stray fibres take their yarn's colour as it stands
    // The fuzz's hair colour from the yarn: "layer" (a thin fuzz layer over the yarn) or "yarn" (a
    // dense coat that keeps the yarn's colour). Strays do not depend on it.
    std::string fibreHairColour = "layer";
    // Yarn colour of tiled scenes: "single" (material.yarn_white on every part), "object" (one
    // colour per knitted part: the area-weighted mean base colour of the part's source surface;
    // no colourwork inside a part) or "texture" (every stitch the mean texture colour over its
    // cell, stitch_colours.h; the default of an imported object with a texture).
    std::string yarnColourSource = "single";
    // Colour of the exported surfaces that are not knitted (plain parts, parts too small to knit,
    // the fur's base meshes): "texture" (the source's base colour texture on its uvs; the default),
    // "average" (one colour per part, the area-weighted mean base
    // colour) or "material" (the source material's base colour: glTF's factor without its texture,
    // the pre-texture export). Without a texture,
    // "texture" falls back to "average".
    std::string plainColourSource = "texture";

    // ---- filling (Export stage): a surface behind each knitted part's yarn (scene/appearance/filling.h)
    bool filling = true;
    double fillingMargin = 0.3;    // [mm] beyond the fabric's inner extent (fillingOffset)

    // ---- yarn body of structural fibres (Export stage; a prototype): each
    // knitted yarn drawn as 1 to 3 plies of long fibres in place of its tube (fibres/yarn_body.h). The
    // fuzz is separate (fibres.enabled). While off, none of these rows counts in a key.
    struct YarnBody {
        bool enabled = false;
        int plies = 3;                     // 1 to 3; needs strand.plies = 1 (the body makes its own plies); 3 is the standard yarn
        int fibres = 160;                  // shell fibres per yarn
        int interiorFibres = 0;            // fibres inside the shell per yarn (the dense reference)
        double fibreDiameterUm = 24;       // [um]
        double shellDepth = 0.35;          // fraction of the ply's centreline radius; 1 = filled
        double fibreAngleDeg = 20;         // surface fibres against the ply axis
        double plyAngleDeg = 20;           // plies against the yarn axis
        double irregularity = 1;           // fibre-scale variation: 0 regular .. 2 (knit.hand.* is the stitch scale)
        std::string core = "none";         // none | opaque | translucent: the core proxy (README)
        int seed = 1;
        // How the yarn of fibres is rendered (yarn_look.h YarnBodyRender): "fibres", the body's strands, or
        // "maps", one tube per yarn under the maps baked from those fibres (fibres/yarn_bake.h).
        std::string render = "fibres";
        double mapTexelUm = 40;            // [um] the maps' texel along and around the yarn (a bake input)
        double mapCavity = 0.6;            // albedo darkening by the baked ambient occlusion (material)
        double mapRoughness = 0.3;         // the smooth fibre's roughness floor under Toksvig's spread (material)
        double mapSpecular = 0.35;         // the highlight's strength, the renderer's `specular` (material)
        double mapAnisotropy = 0.5;        // the highlight's roughness along the yarn over across it (material)
        // The baked tube's material: "fabric" (matte with a sheen; the yarn's translucency passes) or "shiny"
        // (an anisotropic highlight from the roughness and specular maps; opaque). Material, never baked.
        std::string mapMaterial = "fabric";
        bool operator==(const YarnBody&) const = default;
    } yarnBody;

    // ---- display only (no build stage)
    // Run control (no build stage, never in a key): the knitting run's time cap
    // [s]. A part past it is UNRESOLVED ("time cap"), later parts are not started; 0 = none.
    double timeCapS = 3600;

    SceneSettings();               // material and population defaults
    bool operator==(const SceneSettings&) const = default;
};

std::string writeSettings(const SceneSettings& s);
// A retired key (one the program no longer uses) is ignored and listed in retiredKeys, unless it holds the
// default every saved file carries; an unknown key refuses the file. The former workflow selectors
// knit.enabled and fibres.profile are read before any key is applied (legacyWorkflow):
//   - without workflow: knitted wool (knit.enabled true or absent, fibres.profile wool or absent) reads as
//     workflow = knit; the Hair mode (fibres.profile = fur) is refused, its settings cannot be carried over
//     (owner, 2026-11-19); knit.enabled = false with wool, or true with fur, is refused naming both keys;
//   - with workflow: every selector present must agree with it;
//   - a selector given twice, or an unknown workflow, is refused.
bool readSettings(const std::string& text, SceneSettings& s, std::string* err, std::vector<std::string>* retiredKeys = nullptr);
// Keys whose values differ (in table order).
std::vector<std::string> diffSettings(const SceneSettings& a, const SceneSettings& b);
// Assign one key from text ("key=value" CLI overrides); false: unknown key or bad value. A retired
// key is accepted and changes nothing.
bool setSetting(SceneSettings& s, const std::string& key, const std::string& value, std::string* err);
// One key's value copied from another settings object, and compared (numbers within 1e-9).
void copySetting(SceneSettings& to, const SceneSettings& from, const std::string& key);
bool sameSetting(const SceneSettings& a, const SceneSettings& b, const std::string& key);

// One row of the settings table: a key, the stage it invalidates, and where its value lies inside
// SceneSettings. The table's order is part of the stage keys.
enum class SettingKind { Real, Integer, LongInteger, Flag, Text, Choice, Vector3, RealList, MaterialModel };
struct SettingField {
    std::string key;
    Stage stage;
    SettingKind kind;
    size_t offset;
    std::vector<std::string> choices;  // Choice: the accepted words
};
const std::vector<SettingField>& settingFields();
std::string settingValue(const SceneSettings& s, const SettingField& f);
std::string settingNumber(double v);  // a double as settings text: 17 digits, an exact round trip

Workflow workflowOf(const SceneSettings& s);  // unknown text: Knit (the reader refuses it)
inline bool knitWorkflow(const SceneSettings& s) { return workflowOf(s) == Workflow::Knit; }
inline bool furWorkflow(const SceneSettings& s) { return workflowOf(s) == Workflow::Fur; }

// The Fur workflow's switch: off by default (kFurByDefault; the Fur workflow is not released yet); the environment's
// KNIT_FUR=1 turns it on and KNIT_FUR=0 off, read once at the first call. While it is off nothing chooses Fur: readSettings and setSetting refuse
// workflow = fur, and the viewer asks no Knit-or-Fur question (Knit is chosen, workflow.cpp).
constexpr bool kFurByDefault = false;
bool furEnabled();
void setFurEnabled(bool on);  // tests: the switch without the environment

// The guided route is in use: an import, knitted, knit.route = guided. Its rows are stage inputs only then.
bool knitsGuidedPiece(const SceneSettings& s);

} // namespace knit
