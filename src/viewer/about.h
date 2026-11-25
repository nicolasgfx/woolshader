#pragma once
// What the About window says (about_window.cpp draws it): who made Wool Shader and how to reach them,
// its licence, the third-party software in the program with each licence text word for word, the
// example models' credits, and the statements on Blender and on privacy. No ImGui: the tests read all
// of it. The licence texts are compiled in at build time (about_texts.h).
#include <string>
#include <vector>

#include "settings/scene_registry.h"

namespace viewer {

// The program and its author.
struct AboutIdentity {
    std::string product;        // "Wool Shader"
    std::string summary;        // what it does, in one sentence
    std::string copyright;      // "© 2026 Nicolas Menzel"
    std::string email;          // the author's contact address
    std::string portfolioUrl;   // the author's rendering portfolio
    std::string repositoryUrl;  // the source code
    std::string build;          // the commit the program is built from
};
AboutIdentity aboutIdentity();

// One library, or other work, that the program contains or ships next to it.
struct ThirdPartyNotice {
    std::string name;       // "Dear ImGui"
    std::string usedFor;    // what the program uses it for
    std::string licence;    // "MIT"
    std::string version;    // the pinned version
    std::string url;        // where its source is
    std::string fetchName;  // its FetchContent name in cmake/dependencies.cmake; "" when it is not fetched that way
    std::string text;       // its copyright notes and licence text, word for word
};
// In the order the window lists them: the interface, the files, the import, the GPU.
std::vector<ThirdPartyNotice> thirdPartyNotices();

// The project's LICENSE (MIT), and what is under another licence: the Blender script in blender/ and the SLIM files in third_party/slim.
std::string woolShaderLicenceText();
std::string otherLicencesStatement();

// The credits of the example models (the presets of knit::sceneRegistry()), in the registry's order.
std::vector<knit::ModelCredit> exampleModelCredits();
// What applies to the teddy bears: their model, its licence, what an export contains of it.
std::string teddyStatement();
// What applies to the example models in general and to models the user imports.
std::string exampleModelsStatement();
// How the program relates to Blender: a separate program it starts only for the .blend export.
std::string blenderStatement();
// What the program does with data: no network, nothing collected.
std::string privacyStatement();

} // namespace viewer
