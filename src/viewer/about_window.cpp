// The About window (interface.h), opened by About in the status bar: what about.h says, in five tabs.
// Wool Shader (the author, contact, links, build), Licence, Third-party software (each library with its
// licence text word for word, selectable for copying), Example models (credits and what applies to
// them), and Blender and privacy. Links open in the browser or the mail program (ImGui's shell call).
#include <string>
#include <utility>
#include <vector>

#include "imgui.h"
#include "viewer/about.h"
#include "viewer/interface.h"

namespace viewer {

namespace {

// A label and a value on one line; the value is a link when `url` is given.
void identityRow(const char* label, const std::string& value, const std::string& url) {
    ImGui::TextDisabled("%s", label);
    ImGui::SameLine(ImGui::GetFontSize() * 7);
    if (url.empty()) ImGui::TextUnformatted(value.c_str());
    else ImGui::TextLinkOpenURL(value.c_str(), url.c_str());
}

void identityTab(const AboutIdentity& a) {
    ImGui::TextUnformatted(a.product.c_str());
    mutedWrapped(a.summary);
    ImGui::Spacing();
    ImGui::TextUnformatted(a.copyright.c_str());
    ImGui::Spacing();
    identityRow("Contact", a.email, "mailto:" + a.email);
    identityRow("Portfolio", a.portfolioUrl, a.portfolioUrl);
    identityRow("Source code", a.repositoryUrl, a.repositoryUrl);
    identityRow("Build", a.build, "");
    ImGui::Spacing();
    mutedWrapped("Free software under the MIT License, without any warranty (Licence). The libraries it is built with and the "
                 "example models keep their own licences.");
}

// A licence text, word for word, in a read-only box whose text can be selected and copied.
void licenceBox(const char* id, const std::string& text, float height) {
    // ReadOnly: ImGui never writes to the buffer, so the string's own characters can be shown in place.
    ImGui::InputTextMultiline(id, const_cast<char*>(text.c_str()), text.size() + 1, ImVec2(-FLT_MIN, height), ImGuiInputTextFlags_ReadOnly);
}

void licenceTab(const std::string& licence, const std::string& others) {
    licenceBox("##woolShaderLicence", licence, ImGui::GetTextLineHeightWithSpacing() * 14);
    ImGui::Spacing();
    ImGui::TextWrapped("%s", others.c_str());
}

void thirdPartyTab(const std::vector<ThirdPartyNotice>& notices) {
    mutedWrapped("Wool Shader is built with these libraries and contains them. Each keeps its own licence; the texts are word for "
                 "word as their authors publish them.");
    for (size_t i = 0; i < notices.size(); ++i) {
        const ThirdPartyNotice& n = notices[i];
        ImGui::PushID(static_cast<int>(i));
        if (detailsNode(n.name.c_str(), n.licence + ", " + n.version)) {
            identityRow("Used for", n.usedFor, "");
            identityRow("Licence", n.licence, "");
            identityRow("Version", n.version, "");
            identityRow("Source", n.url, n.url);
            licenceBox("##text", n.text, ImGui::GetTextLineHeightWithSpacing() * 12);
            ImGui::TreePop();
        }
        ImGui::PopID();
    }
}

void exampleModelsTab(const std::vector<knit::ModelCredit>& credits, const std::string& teddy, const std::string& models) {
    for (size_t i = 0; i < credits.size(); ++i) {
        const knit::ModelCredit& c = credits[i];
        ImGui::PushID(static_cast<int>(i));
        ImGui::Bullet();
        ImGui::TextWrapped("%s", knit::creditLine(c).c_str());
        ImGui::Indent();
        if (c.source.rfind("http", 0) == 0) ImGui::TextLinkOpenURL(c.source.c_str(), c.source.c_str());
        else ImGui::TextDisabled("%s", c.source.c_str());
        ImGui::TextLinkOpenURL(c.licenceUrl.c_str(), c.licenceUrl.c_str());
        ImGui::Unindent();
        ImGui::PopID();
    }
    if (!teddy.empty()) {
        ImGui::SeparatorText("Teddy bears");
        ImGui::TextWrapped("%s", teddy.c_str());
    }
    ImGui::SeparatorText("Every model");
    ImGui::TextWrapped("%s", models.c_str());
}

// A tab of the window, selected when it is the one asked for (`request`). True when it is open: the caller
// draws it, then ImGui::EndTabItem. (A tab reports its label before it is added: testLabel names it afterwards.)
bool beginTab(const char* label, int request, int tab) {
    const ImGuiID id = ImGui::GetID(label);
    const bool open = ImGui::BeginTabItem(label, nullptr, request == tab ? ImGuiTabItemFlags_SetSelected : ImGuiTabItemFlags_None);
    testLabel(id, label);
    return open;
}

void blenderPrivacyTab(const std::string& blender, const std::string& privacy) {
    ImGui::SeparatorText("Blender");
    ImGui::TextWrapped("%s", blender.c_str());
    ImGui::SeparatorText("Privacy");
    ImGui::TextWrapped("%s", privacy.c_str());
}

} // namespace

void aboutWindow(SceneUi& ui) {
    // The content does not change while the program runs: put together once (the licence texts are long).
    static const AboutIdentity identity = aboutIdentity();
    static const std::vector<ThirdPartyNotice> notices = thirdPartyNotices();
    static const std::string licence = woolShaderLicenceText(), others = otherLicencesStatement();
    static const std::vector<knit::ModelCredit> credits = exampleModelCredits();
    static const std::string teddy = teddyStatement(), models = exampleModelsStatement();
    static const std::string blender = blenderStatement(), privacy = privacyStatement();

    // Centred on the screen each time it opens (it can be moved while open); its size is kept between sessions.
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(uiPx(640), uiPx(520)), ImGuiCond_FirstUseEver);
    // A tab asked for by --do about:<n> is selected once.
    const int request = std::exchange(ui.panels.aboutTabRequest, -1);
    if (ImGui::Begin("About Wool Shader", &ui.panels.aboutOpen) && ImGui::BeginTabBar("##aboutTabs")) {
        if (beginTab("Wool Shader", request, 0)) {
            identityTab(identity);
            ImGui::EndTabItem();
        }
        if (beginTab("Licence", request, 1)) {
            licenceTab(licence, others);
            ImGui::EndTabItem();
        }
        if (beginTab("Third-party software", request, 2)) {
            thirdPartyTab(notices);
            ImGui::EndTabItem();
        }
        if (beginTab("Example models", request, 3)) {
            exampleModelsTab(credits, teddy, models);
            ImGui::EndTabItem();
        }
        if (beginTab("Blender and privacy", request, 4)) {
            blenderPrivacyTab(blender, privacy);
            ImGui::EndTabItem();
        }
        ImGui::EndTabBar();
    }
    ImGui::End();
}

} // namespace viewer
