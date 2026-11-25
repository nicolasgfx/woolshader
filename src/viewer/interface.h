#pragma once
// The viewer's interface: what the ImGui side shares. drawSceneUi draws one frame of it; the look,
// the widgets and the step panels are here so the panels need only this file. There is no GL here
// (the 3D view is renderer.h) and no GLFW, so the UI tests draw all of it without a window. The
// logic the panels show or call is in workflow.h.
#include <string>

#include "imgui.h"
#include "viewer/workflow.h"

namespace viewer {

// The whole interface for a window of this size (strip, panel, toolbar, labels, status bar, dialogs).
void drawSceneUi(SceneUi& ui, float width, float height);
void frameStatsLabel(const SceneUi& ui, float width, float height, double frameMs, long long segments);

// --- the steps (panel_*.cpp)

void importPanel(SceneUi& ui);
void importFooter(SceneUi& ui);
void markPanel(SceneUi& ui);  // panel_mark.cpp
bool markFooterShown(const SceneUi& ui);
void markFooter(SceneUi& ui);
void knitPanel(SceneUi& ui, const knit::SceneEntry& e);
void knitFooter(SceneUi& ui, const knit::SceneEntry& e);
void partsTable(SceneUi& ui, const knit::BuiltScene& b);  // in the Report window
void fibrePanel(SceneUi& ui, const knit::SceneEntry& e);
void yarnFibresSection(SceneUi& ui);  // the Fibres step's first section for a yarn made of fibres (panel_yarn_fibres.cpp)
void fibreFooter(SceneUi& ui);
void furPanel(SceneUi& ui);  // the Fur step, in the Fibres step's place under Fur (panel_fur.cpp)
void exportPanel(SceneUi& ui, const knit::SceneEntry& e);
void exportFooter(SceneUi& ui);

// --- the About window (about_window.cpp): the author, the licences, the credits

void aboutWindow(SceneUi& ui);

// --- the look and its scale (widgets.cpp)

const ImVec4 kAccent{0.26f, 0.52f, 0.90f, 1.0f}, kAccentHover{0.33f, 0.59f, 0.96f, 1.0f}, kAccentActive{0.20f, 0.44f, 0.80f, 1.0f};
const ImVec4 kOk{0.40f, 0.78f, 0.46f, 1.0f}, kWarn{0.96f, 0.70f, 0.28f, 1.0f}, kBad{0.94f, 0.40f, 0.36f, 1.0f};
const ImVec4 kInfo{0.55f, 0.75f, 1.0f, 1.0f};

ImVec4 muted();
void applySceneTheme(ImGuiStyle& st);  // logical sizes and colours

// Interface scale = the display's content scale x a user factor. Takes effect between frames: call
// uiScaleBeginFrame before NewFrame. True: the style and font atlas were rebuilt, so the caller
// makes a new font texture.
bool uiScaleBeginFrame(float contentScale);
void uiScaleSetUser(float userScale);  // clamped to [0.5, 3]
float uiScaleUser();
float uiScale();
inline float uiPx(float px) { return px * uiScale(); }  // a hard-coded layout size
void uiScaleSlider();

// --- words and states

void glyph(StageStatus s, float size = 0);  // its shape differs per state: colour is never the only cue
void drawGlyph(ImDrawList* dl, ImVec2 topLeft, float size, StageStatus s);  // the glyph alone, into a draw list (a tile's corner)
void help(const char* text);
void tooltip(const char* text);
void wrapped(const ImVec4& col, const std::string& s);
void mutedWrapped(const std::string& s);
void noteWithHelp(const ImVec4& col, const std::string& s, const char* helpText);
void stateCard(const StageState& st, const std::string& text);
void stepHeading(const char* title, const std::string& subtitle);
// A section's heading: its title, a rule, and a short muted note at the right edge (e.g. "Custom").
void sectionHeading(const char* title, const std::string& note);
// A collapsed section with its current settings in one muted line after the title (clipped to the
// panel). True when open: the caller draws the contents, then ImGui::TreePop.
bool detailsNode(const char* title, const std::string& summary);
void attemptNote(const Attempt& a, const std::string& alreadySaid);  // its message, unless already shown

// --- buttons and footers

bool primaryButton(const char* label, bool enabled, const char* whyDisabled = nullptr, float width = -FLT_MIN);  // full width by default
bool plainButton(const char* label, bool enabled, const char* whyDisabled, float width);  // an action beside a primary one
bool cancelButton(const char* label, bool requested);
bool choiceButton(const char* label, bool on, float width = 0);
// A toolbar button that opens a popover: the text, a drawn down-arrow, and an id that stays when the text
// changes ("Front" -> "Orbit"). menuButtonWidth: its width, for laying a row out before drawing it.
bool menuButton(const char* shown, const char* id, bool open);
float menuButtonWidth(const char* shown);
void beginFooterStatus(const char* id, float reserve);
void jobProgressBlock(const JobProgressView& v, const std::string& rightText);

// --- fields and rows

struct NumberRange {
    double lo, hi;
    const char* format = "%.2f";
};
// The Knit step's yarn and stitch fields. The stitch ranges hold every stitch the Gauge list's
// Snug / Regular / Loose derive over the diameter and clearance ranges (stitchSizeFromYarn: from
// 2.20 x 1.35 mm at Snug, 0.50 mm, no clearance, to 20.00 x 14.70 mm at Loose, 2.00 mm, 0.40 mm),
// so a switch to Custom never clamps on the first drag.
constexpr NumberRange kYarnDiameterRange{0.5, 2.0, "%.2f mm"};  // the Dense preset draws 1.86 (its stitch takes up to 1.94)
constexpr NumberRange kClearanceRange{0.0, 0.4, "%.2f mm"};
constexpr NumberRange kStitchWidthRange{2.0, 20.0, "%.2f mm"};
constexpr NumberRange kRowHeightRange{1.2, 15.0, "%.2f mm"};
bool sliderField(const char* label, double* v, const NumberRange& range);
bool logSliderField(const char* label, double* v, const NumberRange& range);
bool dragField(const char* id, double* v, const NumberRange& range);
bool dragField(const char* id, int* v, int lo, int hi);

// Label / value rows, one setting one statement (inside beginRows ... ImGui::EndTable). The field's
// id is "##" + the label; `row` + a field is for the few whose id differs.
bool beginRows(const char* id);
void row(const char* label, const char* helpText = nullptr);
bool dragRow(const char* label, double* v, const NumberRange& range, const char* helpText = nullptr);
bool dragRow(const char* label, int* v, int lo, int hi, const char* helpText = nullptr);
bool comboRow(const char* label, int* choice, const char* const items[], int count, const char* helpText = nullptr);
// A share in [0, 1] edited as a percentage ("60 %" for 0.6), within lo..hi (shares).
bool percentRow(const char* label, double* share, double lo, double hi, const char* helpText = nullptr);
void colourRow(const char* label, knit::Vec3& colour);

void testLabel(ImGuiID id, const char* label);
// ImGui::TextLink with its label reported to the UI tests (ImGui reports none for a link).
bool textLink(const char* label);
bool combo(const char* id, int* v, const char* const items[], int n);
bool colorEdit(const char* label, knit::Vec3& c);
// A combo the caller fills with a plain loop, so a closed list builds nothing (ends with EndCombo).
// flags: e.g. ImGuiComboFlags_HeightLargest for a list with captions that must show every entry.
bool beginList(const char* id, const char* shown, ImGuiComboFlags flags = 0);
// One row, one list: the group's presets and Custom. Choosing one sets every value it owns; the
// preset the settings match is shown, Custom when a value was edited. Returns the preset now shown.
const knit::FilePreset* presetRow(knit::SceneSettings& s, knit::PresetGroup group, const char* label, const char* customNote);

// --- tiles: a few exclusive choices as a row of pictures with a label (Stitch character, Construction)

// The row's layout from the panel's width: 4 across, 2 x 2 or 1 per row when narrow. withDetail adds two
// muted lines of detail under each label (wider tiles: 2 x 2 at the default panel width).
struct TileRow {
    int perRow = 4;
    float tileW = 0, tileH = 0, pad = 0, pictureH = 0;
    bool withDetail = false;
};
TileRow beginTiles(int count, bool withDetail);
struct Tile {
    const char* id;      // the ImGui id (a preset id, a construction's name)
    const char* label;
    const char* detail;  // the detail under the label, wrapped onto two lines ("" or null: none)
    const char* note;    // the tooltip ("" or null: none)
    bool selected;
};
// The picture area of a tile, for the caller to draw into after tile() returned.
struct TilePicture {
    ImVec2 lo, hi;
};
// Draws the tile `index` of the row: the selectable, its border, its label (and detail), the tooltip and
// the test label; returns true when clicked.
bool tile(const TileRow& r, int index, const Tile& t, TilePicture* picture);

} // namespace viewer
