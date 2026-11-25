// The scene viewer's look and its widgets (interface.h): colours and scale, state glyphs, buttons,
// settings rows and lists. No GL and no GLFW: main.cpp reads the display's scale from the window and
// keeps the font texture.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>

#include "imgui.h"
#include "imgui_internal.h"  // ImGui::RenderTextEllipsis (detailsNode); the UI tests' item hooks (testLabel)
#include "viewer/interface.h"
#include "settings/file_presets.h"

using knit::Vec3;

namespace viewer {

// --- scale: the display's content scale x the user's factor (1.0: the operating system's scaling)

namespace {

float g_user = 1.0f, g_content = 1.0f, g_applied = 0.0f;
bool g_rebuild = true;

std::string fontFile() {
    namespace fs = std::filesystem;
    std::error_code ec;
    for (const std::string& dir : {std::string("fonts"),
#ifdef KNIT_FONT_DIR
                                   std::string(KNIT_FONT_DIR),
#endif
         })
        if (fs::exists(dir + "/Roboto-Medium.ttf", ec)) return dir + "/Roboto-Medium.ttf";
    return "";
}

// The style and the font at the current scale. ScaleAllSizes is cumulative, hence the reset; the
// theme sets logical sizes, so it comes first and they are scaled once with the rest.
void rebuildLook() {
    const float s = g_content * g_user;
    ImGuiStyle& st = ImGui::GetStyle();
    st = ImGuiStyle();
    ImGui::StyleColorsDark(&st);
    applySceneTheme(st);
    st.ScaleAllSizes(s);
    ImGuiIO& io = ImGui::GetIO();
    io.Fonts->Clear();
    const std::string f = fontFile();
    ImFont* font = nullptr;
    if (!f.empty()) font = io.Fonts->AddFontFromFileTTF(f.c_str(), std::round(16.0f * s));
    if (!font) {
        ImFontConfig cfg;
        cfg.SizePixels = std::round(13.0f * s);
        io.Fonts->AddFontDefault(&cfg);
    }
    io.Fonts->Build();
    g_applied = s;
    g_rebuild = false;
}

} // namespace

bool uiScaleBeginFrame(float contentScale) {
    const float content = std::max(0.5f, contentScale);
    if (content != g_content) {
        g_content = content;
        g_rebuild = true;
    }
    if (!g_rebuild) return false;
    rebuildLook();
    return true;
}

void uiScaleSetUser(float userScale) {
    const float u = std::clamp(userScale, 0.5f, 3.0f);
    if (u != g_user) {
        g_user = u;
        g_rebuild = true;
    }
}

float uiScaleUser() { return g_user; }

float uiScale() { return g_applied > 0 ? g_applied : g_content * g_user; }

void uiScaleSlider() {
    float u = g_user;
    // Applied on release: rebuilding the font atlas every drag frame would reflow the slider.
    if (ImGui::SliderFloat("UI scale", &u, 0.5f, 3.0f, "%.2f")) g_user = std::clamp(u, 0.5f, 3.0f);
    if (ImGui::IsItemDeactivatedAfterEdit()) g_rebuild = true;
    ImGui::SameLine();
    ImGui::TextDisabled("(x %.2f display)", g_content);
}

void applySceneTheme(ImGuiStyle& st) {
    // Sizes in logical pixels (rebuildLook scales them afterwards).
    st.WindowRounding = 0;
    st.ChildRounding = 0;
    st.FrameRounding = 3;
    st.GrabRounding = 3;
    st.PopupRounding = 3;
    st.ScrollbarRounding = 3;
    st.TabRounding = 3;
    st.WindowPadding = ImVec2(12, 10);
    st.FramePadding = ImVec2(7, 4);
    st.ItemSpacing = ImVec2(8, 6);
    st.SeparatorTextPadding = ImVec2(0, 8);
    st.SeparatorTextBorderSize = 1;
    // Colours: the panels a little lighter than the viewport (renderer.cpp clears to 0.11 0.12
    // 0.14), inputs distinct from the panel, a restrained blue for choices and the main action,
    // neutral separators, readable secondary text.
    ImVec4* c = st.Colors;
    const ImVec4 panel(0.145f, 0.155f, 0.175f, 1), line(0.27f, 0.29f, 0.33f, 1);
    c[ImGuiCol_Text] = ImVec4(0.88f, 0.90f, 0.93f, 1);
    c[ImGuiCol_TextDisabled] = ImVec4(0.62f, 0.66f, 0.72f, 1);
    c[ImGuiCol_WindowBg] = panel;
    c[ImGuiCol_ChildBg] = ImVec4(0, 0, 0, 0);
    c[ImGuiCol_PopupBg] = ImVec4(0.12f, 0.13f, 0.15f, 0.98f);
    c[ImGuiCol_Border] = line;
    c[ImGuiCol_Separator] = line;
    c[ImGuiCol_FrameBg] = ImVec4(0.20f, 0.22f, 0.26f, 1);
    c[ImGuiCol_FrameBgHovered] = ImVec4(0.24f, 0.27f, 0.32f, 1);
    c[ImGuiCol_FrameBgActive] = ImVec4(0.27f, 0.31f, 0.37f, 1);
    c[ImGuiCol_Button] = ImVec4(0.21f, 0.23f, 0.27f, 1);
    c[ImGuiCol_ButtonHovered] = ImVec4(0.26f, 0.30f, 0.36f, 1);
    c[ImGuiCol_ButtonActive] = ImVec4(0.29f, 0.35f, 0.44f, 1);
    c[ImGuiCol_Header] = ImVec4(0.19f, 0.31f, 0.44f, 1);
    c[ImGuiCol_HeaderHovered] = ImVec4(0.23f, 0.36f, 0.50f, 1);
    c[ImGuiCol_HeaderActive] = ImVec4(0.26f, 0.41f, 0.57f, 1);
    c[ImGuiCol_CheckMark] = ImVec4(0.55f, 0.74f, 0.95f, 1);
    c[ImGuiCol_SliderGrab] = ImVec4(0.42f, 0.60f, 0.85f, 1);
    c[ImGuiCol_SliderGrabActive] = ImVec4(0.50f, 0.68f, 0.92f, 1);
    c[ImGuiCol_TitleBg] = c[ImGuiCol_TitleBgCollapsed] = ImVec4(0.12f, 0.13f, 0.15f, 1);
    c[ImGuiCol_TitleBgActive] = ImVec4(0.16f, 0.20f, 0.26f, 1);
    c[ImGuiCol_PlotHistogram] = ImVec4(0.42f, 0.60f, 0.85f, 1);  // progress bars
    c[ImGuiCol_TableBorderLight] = line;
    // Tabs (the Fibres step's population editor): neutral, the chosen one in the selection colour
    // with an accent line, the same whether or not the panel has focus.
    c[ImGuiCol_Tab] = c[ImGuiCol_TabDimmed] = c[ImGuiCol_Button];
    c[ImGuiCol_TabHovered] = c[ImGuiCol_ButtonHovered];
    c[ImGuiCol_TabSelected] = c[ImGuiCol_TabDimmedSelected] = c[ImGuiCol_Header];
    c[ImGuiCol_TabSelectedOverline] = c[ImGuiCol_TabDimmedSelectedOverline] = c[ImGuiCol_CheckMark];
}

// --- words and states

ImVec4 muted() { return ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled); }

namespace {

ImU32 stColor(StageStatus s) {
    switch (s) {
    case StageStatus::Ready: return ImGui::GetColorU32(kOk);
    case StageStatus::Diagnostic:
    case StageStatus::NeedsUpdate:
    case StageStatus::ViewOnly: return ImGui::GetColorU32(kWarn);
    case StageStatus::Failed: return ImGui::GetColorU32(kBad);
    case StageStatus::Running:
    case StageStatus::Cancelling: return ImGui::GetColorU32(kInfo);
    default: return ImGui::GetColorU32(muted());
    }
}

} // namespace

// A state glyph (its shape differs per state, so colour is never the only cue), at the cursor.
void drawGlyph(ImDrawList* dl, ImVec2 topLeft, float size, StageStatus s) {
    const float r = size * 0.36f;
    const ImVec2 c{topLeft.x + size * 0.5f, topLeft.y + size * 0.5f};
    const ImU32 col = stColor(s);
    const float th = std::max(1.0f, size * 0.09f);
    switch (s) {
    case StageStatus::Ready:
        dl->AddCircleFilled(c, r, col);
        dl->AddLine({c.x - r * 0.5f, c.y}, {c.x - r * 0.1f, c.y + r * 0.42f}, IM_COL32(20, 30, 20, 255), th * 1.3f);
        dl->AddLine({c.x - r * 0.1f, c.y + r * 0.42f}, {c.x + r * 0.55f, c.y - r * 0.45f}, IM_COL32(20, 30, 20, 255), th * 1.3f);
        break;
    case StageStatus::Diagnostic:
    case StageStatus::NeedsUpdate:
    case StageStatus::ViewOnly:
        dl->AddTriangleFilled({c.x, c.y - r}, {c.x + r, c.y + r * 0.8f}, {c.x - r, c.y + r * 0.8f}, col);
        dl->AddLine({c.x, c.y - r * 0.35f}, {c.x, c.y + r * 0.25f}, IM_COL32(30, 25, 10, 255), th * 1.2f);
        dl->AddCircleFilled({c.x, c.y + r * 0.55f}, th * 0.8f, IM_COL32(30, 25, 10, 255));
        break;
    case StageStatus::Failed:
        dl->AddCircleFilled(c, r, col);
        dl->AddLine({c.x - r * 0.45f, c.y - r * 0.45f}, {c.x + r * 0.45f, c.y + r * 0.45f}, IM_COL32(40, 10, 10, 255), th * 1.3f);
        dl->AddLine({c.x + r * 0.45f, c.y - r * 0.45f}, {c.x - r * 0.45f, c.y + r * 0.45f}, IM_COL32(40, 10, 10, 255), th * 1.3f);
        break;
    case StageStatus::Running:
    case StageStatus::Cancelling: {
        const float t = static_cast<float>(ImGui::GetTime() * 4.0);
        dl->PathArcTo(c, r, t, t + 4.4f, 16);
        dl->PathStroke(col, 0, th * 1.4f);
        break;
    }
    case StageStatus::Cancelled:
        dl->AddCircle(c, r, col, 16, th);
        dl->AddLine({c.x - r * 0.7f, c.y + r * 0.7f}, {c.x + r * 0.7f, c.y - r * 0.7f}, col, th);
        break;
    case StageStatus::Blocked:
        dl->AddRect({c.x - r * 0.8f, c.y - r * 0.8f}, {c.x + r * 0.8f, c.y + r * 0.8f}, col, 0, 0, th);
        break;
    case StageStatus::Optional:
        dl->AddCircle(c, r * 0.6f, col, 12, th);
        break;
    case StageStatus::NotStarted:
    default: dl->AddCircle(c, r, col, 16, th);
    }
}

// The glyph as an item on the text line (centred on it), with its room.
void glyph(StageStatus s, float size) {
    if (size <= 0) size = ImGui::GetFontSize();
    const ImVec2 p = ImGui::GetCursorScreenPos();
    drawGlyph(ImGui::GetWindowDrawList(), {p.x, p.y + (ImGui::GetTextLineHeight() - size) * 0.5f}, size, s);
    ImGui::Dummy({size, ImGui::GetTextLineHeight()});
}

void tooltip(const char* text) {
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled | ImGuiHoveredFlags_DelayShort)) {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 28);
        ImGui::TextUnformatted(text);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}

void help(const char* text) {
    ImGui::SameLine();
    ImGui::TextDisabled("(?)");
    tooltip(text);
}

void wrapped(const ImVec4& col, const std::string& s) {
    ImGui::PushStyleColor(ImGuiCol_Text, col);
    ImGui::TextWrapped("%s", s.c_str());
    ImGui::PopStyleColor();
}

void mutedWrapped(const std::string& s) { wrapped(muted(), s); }

// A wrapped note with its help marker at the end of the first line: the text wraps short of the
// panel's edge so the marker stays inside.
void noteWithHelp(const ImVec4& col, const std::string& s, const char* helpText) {
    ImGui::PushTextWrapPos(ImGui::GetContentRegionMax().x - ImGui::CalcTextSize("(?)").x - ImGui::GetStyle().ItemSpacing.x);
    wrapped(col, s);
    ImGui::PopTextWrapPos();
    help(helpText);
}

// A full-width primary action (the accent colour is reserved for these).
bool primaryButton(const char* label, bool enabled, const char* whyDisabled, float width) {
    ImGui::BeginDisabled(!enabled);
    ImGui::PushStyleColor(ImGuiCol_Button, kAccent);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, kAccentHover);
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, kAccentActive);
    ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1, 1, 1, 1));
    const bool hit = ImGui::Button(label, ImVec2(width, ImGui::GetFrameHeight() * 1.35f));
    ImGui::PopStyleColor(4);
    ImGui::EndDisabled();
    if (!enabled && whyDisabled) tooltip(whyDisabled);
    return hit && enabled;
}

// A normal button for an action beside a primary one (the Knit footer's "Knit draft"), at the
// primary buttons' height so the row lines up.
bool plainButton(const char* label, bool enabled, const char* whyDisabled, float width) {
    ImGui::BeginDisabled(!enabled);
    const bool hit = ImGui::Button(label, ImVec2(width, ImGui::GetFrameHeight() * 1.35f));
    ImGui::EndDisabled();
    if (!enabled && whyDisabled) tooltip(whyDisabled);
    return hit && enabled;
}

// A job's Cancel: a normal secondary button (red is kept for errors).
bool cancelButton(const char* label, bool requested) {
    ImGui::BeginDisabled(requested);
    const bool hit = ImGui::Button(label, ImVec2(-FLT_MIN, ImGui::GetFrameHeight() * 1.35f));
    ImGui::EndDisabled();
    return hit;
}

void attemptNote(const Attempt& a, const std::string& alreadySaid) {
    if (!a.message.empty() && a.message != alreadySaid) wrapped(a.outcome == knit::JobOutcome::Failed ? kBad : muted(), a.message);
}

// A step's heading: what it is for and what it applies to.
void stepHeading(const char* title, const std::string& subtitle) {
    ImGui::TextUnformatted(title);
    mutedWrapped(subtitle);
}

// As ImGui::SeparatorText, with the note right-aligned after the rule. A note that does not fit
// beside the title is left out of the line (the section's rows say the same).
void sectionHeading(const char* title, const std::string& note) {
    ImGui::Spacing();
    const ImVec2 start = ImGui::GetCursorScreenPos();
    const float right = start.x + ImGui::GetContentRegionAvail().x;
    const float gap = ImGui::GetStyle().ItemSpacing.x;
    ImGui::TextUnformatted(title);
    float ruleEnd = right;
    const float noteWidth = ImGui::CalcTextSize(note.c_str()).x;
    const float titleEnd = start.x + ImGui::CalcTextSize(title).x;
    if (!note.empty() && titleEnd + gap * 4 + noteWidth <= right) {
        ImGui::SameLine();
        ImGui::SetCursorScreenPos({right - noteWidth, start.y});
        ImGui::PushStyleColor(ImGuiCol_Text, muted());
        ImGui::TextUnformatted(note.c_str());
        ImGui::PopStyleColor();
        ruleEnd = right - noteWidth - gap;
    }
    const float y = start.y + ImGui::GetTextLineHeight() * 0.5f;
    if (ruleEnd > titleEnd + gap * 2)
        ImGui::GetWindowDrawList()->AddLine({titleEnd + gap, y}, {ruleEnd, y}, ImGui::GetColorU32(ImGuiCol_Separator));
}

bool detailsNode(const char* title, const std::string& summary) {
    const bool open = ImGui::TreeNodeEx(title, ImGuiTreeNodeFlags_SpanAvailWidth);  // the summary is not in the id: it stays
    if (summary.empty()) return open;
    // The summary sits on the node's row, after its title, clipped with "..." at the panel's edge.
    const ImRect node(ImGui::GetItemRectMin(), ImGui::GetItemRectMax());
    const ImGuiStyle& st = ImGui::GetStyle();
    const float x0 = node.Min.x + ImGui::GetTreeNodeToLabelSpacing() + ImGui::CalcTextSize(title).x + st.ItemSpacing.x * 2;
    const float x1 = node.Max.x - st.FramePadding.x;
    if (x1 > x0 + ImGui::GetFontSize() * 2) {
        const float y = node.Min.y + (node.GetHeight() - ImGui::GetTextLineHeight()) * 0.5f;
        ImGui::PushStyleColor(ImGuiCol_Text, muted());
        ImGui::RenderTextEllipsis(ImGui::GetWindowDrawList(), {x0, y}, {x1, node.Max.y}, x1, x1, summary.c_str(), nullptr, nullptr);
        ImGui::PopStyleColor();
    }
    return open;
}

// A stage's state card: glyph, state, details.
void stateCard(const StageState& st, const std::string& text) {
    glyph(st.status);
    ImGui::SameLine();
    ImGui::TextUnformatted(st.word.c_str());
    if (!text.empty()) {
        ImGui::SameLine();
        ImGui::PushStyleColor(ImGuiCol_Text, muted());
        ImGui::TextWrapped("%s", text.c_str());
        ImGui::PopStyleColor();
    }
}

// --- fields and rows

bool sliderField(const char* label, double* v, const NumberRange& r) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &r.lo, &r.hi, r.format);
}

bool logSliderField(const char* label, double* v, const NumberRange& r) {
    return ImGui::SliderScalar(label, ImGuiDataType_Double, v, &r.lo, &r.hi, r.format, ImGuiSliderFlags_Logarithmic);
}

// A number field: drag it, or double-click / ctrl-click to type a value; clamped to the range. The
// unit is part of the format.
bool dragField(const char* id, double* v, const NumberRange& r) {
    return ImGui::DragScalar(id, ImGuiDataType_Double, v, float((r.hi - r.lo) / 400), &r.lo, &r.hi, r.format, ImGuiSliderFlags_AlwaysClamp);
}

bool dragField(const char* id, int* v, int lo, int hi) { return ImGui::DragInt(id, v, 0.1f, lo, hi, "%d", ImGuiSliderFlags_AlwaysClamp); }

// Aligned label / value rows: a two-column table, the label left, the control filling the right
// column (its id is "##" + the label). Sections end the table before notes that span the panel.
bool beginRows(const char* id) {
    if (!ImGui::BeginTable(id, 2, ImGuiTableFlags_SizingStretchProp)) return false;
    ImGui::TableSetupColumn("##label", ImGuiTableColumnFlags_WidthStretch, 0.45f);
    ImGui::TableSetupColumn("##value", ImGuiTableColumnFlags_WidthStretch, 0.55f);
    return true;
}

void row(const char* label, const char* helpText) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::AlignTextToFramePadding();
    ImGui::TextUnformatted(label);
    if (helpText && *helpText) help(helpText);
    ImGui::TableNextColumn();
    ImGui::SetNextItemWidth(-FLT_MIN);
}

// One of a few exclusive choices as a normal button; the chosen one in the selection colour.
bool choiceButton(const char* label, bool on, float width) {
    if (on) {
        ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_Header));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImGui::GetStyleColorVec4(ImGuiCol_HeaderHovered));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImGui::GetStyleColorVec4(ImGuiCol_HeaderActive));
    }
    const bool hit = ImGui::Button(label, ImVec2(width, 0));
    if (on) ImGui::PopStyleColor(3);
    return hit;
}

// The down-arrow's room after a menu button's text.
float menuArrowWidth() { return ImGui::GetFontSize() * 0.9f; }

float menuButtonWidth(const char* shown) {
    return ImGui::CalcTextSize(shown).x + ImGui::GetStyle().ItemInnerSpacing.x + menuArrowWidth() + ImGui::GetStyle().FramePadding.x * 2;
}

// The text, then a small down-pointing triangle drawn in the text colour (the font atlas has no glyph for it).
// "###id" keeps the id when the text changes.
bool menuButton(const char* shown, const char* id, bool open) {
    const std::string label = std::string(shown) + "###" + id;
    ImGui::PushStyleVar(ImGuiStyleVar_ButtonTextAlign, ImVec2(0.0f, 0.5f));  // the text left, the arrow right
    const bool hit = choiceButton(label.c_str(), open, menuButtonWidth(shown));
    ImGui::PopStyleVar();
    const ImVec2 lo = ImGui::GetItemRectMin(), hi = ImGui::GetItemRectMax();
    const float s = ImGui::GetFontSize() * 0.32f, cx = hi.x - ImGui::GetStyle().FramePadding.x - menuArrowWidth() * 0.5f, cy = (lo.y + hi.y) * 0.5f;
    ImGui::GetWindowDrawList()->AddTriangleFilled({cx - s, cy - s * 0.5f}, {cx + s, cy - s * 0.5f}, {cx, cy + s * 0.6f}, ImGui::GetColorU32(ImGuiCol_Text));
    return hit;
}

// The UI tests' item hooks (IMGUI_ENABLE_TEST_ENGINE, tests/tests_ui.cpp) learn the label of a
// widget ImGui does not report one for (combos, colour buttons), by the id taken before the
// widget. Nothing in the viewer.
void testLabel(ImGuiID id, const char* label) {
#ifdef IMGUI_ENABLE_TEST_ENGINE
    ImGuiContext& g = *GImGui;
    IMGUI_TEST_ENGINE_ITEM_INFO(id, label, 0);
#else
    (void)id;
    (void)label;
#endif
}

bool combo(const char* id, int* v, const char* const items[], int n) {
    const ImGuiID iid = ImGui::GetID(id);
    const bool changed = ImGui::Combo(id, v, items, n);
    testLabel(iid, id);
    return changed;
}

bool colorEdit(const char* label, Vec3& c) {
    ImGui::PushID(label);
    const ImGuiID button = ImGui::GetID("##ColorButton");  // ColorEdit3's swatch (NoInputs)
    ImGui::PopID();
    float f[3] = {float(c.x), float(c.y), float(c.z)};
    const bool changed = ImGui::ColorEdit3(label, f, ImGuiColorEditFlags_NoInputs);
    testLabel(button, label);
    if (!changed) return false;
    c = {f[0], f[1], f[2]};
    return true;
}

namespace {

std::string fieldId(const char* label) { return std::string("##") + label; }

} // namespace

bool dragRow(const char* label, double* v, const NumberRange& range, const char* helpText) {
    row(label, helpText);
    return dragField(fieldId(label).c_str(), v, range);
}

bool dragRow(const char* label, int* v, int lo, int hi, const char* helpText) {
    row(label, helpText);
    return dragField(fieldId(label).c_str(), v, lo, hi);
}

bool comboRow(const char* label, int* choice, const char* const items[], int count, const char* helpText) {
    row(label, helpText);
    return combo(fieldId(label).c_str(), choice, items, count);
}

void colourRow(const char* label, knit::Vec3& colour) {
    row(label);
    colorEdit(fieldId(label).c_str(), colour);
}

bool percentRow(const char* label, double* share, double lo, double hi, const char* helpText) {
    double percent = *share * 100;
    if (!dragRow(label, &percent, {lo * 100, hi * 100, "%.0f %%"}, helpText)) return false;
    *share = percent / 100;
    return true;
}

bool textLink(const char* label) {
    const ImGuiID id = ImGui::GetID(label);
    const bool hit = ImGui::TextLink(label);
    testLabel(id, label);
    return hit;
}

// A combo the caller fills: ImGui reports no label for it, so the tests' hook gets one here.
bool beginList(const char* id, const char* shown, ImGuiComboFlags flags) {
    const ImGuiID item = ImGui::GetID(id);
    const bool open = ImGui::BeginCombo(id, shown, flags);
    testLabel(item, id);
    return open;
}

const knit::FilePreset* presetRow(knit::SceneSettings& s, knit::PresetGroup group, const char* label, const char* customNote) {
    const knit::FilePreset* match = knit::matchingFilePreset(group, s);
    row(label, match ? match->note.c_str() : customNote);
    if (beginList((std::string("##") + label).c_str(), match ? match->label.c_str() : "Custom")) {
        for (const knit::FilePreset& p : knit::presets(group)) {
            if (ImGui::Selectable(p.label.c_str(), match == &p)) {
                knit::applyPreset(group, s, p);
                match = &p;
            }
            tooltip((p.note + " (" + p.path + ")").c_str());
        }
        ImGui::Selectable("Custom", match == nullptr);
        tooltip(customNote);
        ImGui::EndCombo();
    }
    return match;
}

// --- footers

// The footers' status block (Knit, Fibres, Export): it grows with its content up to 40 % of the
// panel less `reserve` (what the footer draws after it), then scrolls; the caller ends it with
// EndChild.
void beginFooterStatus(const char* id, float reserve) {
    ImGui::Separator();
    const float budget = std::max(ImGui::GetFrameHeightWithSpacing() * 2, ImGui::GetWindowHeight() * 0.4f - reserve);
    ImGui::SetNextWindowSizeConstraints(ImVec2(0, 0), ImVec2(FLT_MAX, budget));
    ImGui::BeginChild(id, ImVec2(0, 0), ImGuiChildFlags_AutoResizeY);
}

// A running job in its footer: the state and title, a count on the right, a thin bar (indeterminate
// without a total) and what the job is doing. Parts done or work done, never a time estimate.
void jobProgressBlock(const JobProgressView& v, const std::string& rightText) {
    glyph(v.stopping ? StageStatus::Cancelling : StageStatus::Running);
    ImGui::SameLine();
    ImGui::TextUnformatted(v.title.c_str());
    if (!rightText.empty()) {
        const float width = ImGui::CalcTextSize(rightText.c_str()).x;
        ImGui::SameLine(std::max(ImGui::GetCursorPosX(), ImGui::GetContentRegionMax().x - width));
        ImGui::TextUnformatted(rightText.c_str());
    }
    ImGui::ProgressBar(v.fraction >= 0 ? v.fraction : -1.0f * float(ImGui::GetTime()), ImVec2(-FLT_MIN, std::max(4.0f, ImGui::GetFontSize() * 0.4f)), "");
    if (!v.qualifier.empty()) mutedWrapped(v.qualifier);
}

// --- tiles

TileRow beginTiles(int count, bool withDetail) {
    const ImGuiStyle& st = ImGui::GetStyle();
    const float avail = ImGui::GetContentRegionAvail().x;
    // A tile with a detail is wide enough for its two wrapped lines: 2 x 2 at the default panel width.
    const float minTile = ImGui::GetFontSize() * (withDetail ? 13 : 6);
    TileRow r;
    r.perRow = std::max(1, count);
    while (r.perRow > 1 && r.perRow * minTile + (r.perRow - 1) * st.ItemSpacing.x > avail) r.perRow = (r.perRow + 1) / 2;
    r.tileW = (avail - (r.perRow - 1) * st.ItemSpacing.x) / r.perRow;
    r.pad = st.FramePadding.y + 1;
    r.pictureH = ImGui::GetFontSize() * 1.8f;
    r.withDetail = withDetail;
    r.tileH = r.pad * 3 + r.pictureH + ImGui::GetTextLineHeight() * (withDetail ? 3 : 1);  // the detail wraps onto two lines
    return r;
}

// The tile's label is drawn with RenderTextClipped, which ImGui's log captures for the UI tests' text snapshot.
bool tile(const TileRow& r, int index, const Tile& t, TilePicture* picture) {
    const ImGuiStyle& st = ImGui::GetStyle();
    if (index % r.perRow != 0) ImGui::SameLine();
    ImGui::PushID(t.id);
    const ImVec2 lo = ImGui::GetCursorScreenPos();
    const ImVec2 hi{lo.x + r.tileW, lo.y + r.tileH};
    const ImGuiID id = ImGui::GetID("##tile");
    const bool hit = ImGui::Selectable("##tile", t.selected, ImGuiSelectableFlags_None, {r.tileW, r.tileH});
    testLabel(id, t.label);
    if (t.note && *t.note) tooltip(t.note);
    ImGui::GetWindowDrawList()->AddRect(lo, hi, ImGui::GetColorU32(t.selected ? kAccent : st.Colors[ImGuiCol_Border]), st.FrameRounding, 0,
                                        t.selected ? 2.0f : 1.0f);
    picture->lo = {lo.x + r.pad, lo.y + r.pad};
    picture->hi = {hi.x - r.pad, lo.y + r.pad + r.pictureH};
    const float textY = lo.y + r.pad * 2 + r.pictureH;
    ImGui::RenderTextClipped({lo.x + r.pad * 2, textY}, {hi.x - r.pad, hi.y}, t.label, nullptr, nullptr);
    if (r.withDetail && t.detail && *t.detail) {
        ImGui::PushStyleColor(ImGuiCol_Text, muted());
        ImGui::RenderTextWrapped({lo.x + r.pad * 2, textY + ImGui::GetTextLineHeight()}, t.detail, nullptr, r.tileW - r.pad * 3);
        ImGui::PopStyleColor();
    }
    ImGui::PopID();
    return hit;
}

} // namespace viewer
