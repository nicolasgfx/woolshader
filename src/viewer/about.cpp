// What the About window says (about.h). The words are here, the licence texts come from the build
// (about_texts.h), and the example models' credits from the scene registry, so a credit is written once.
#include "viewer/about.h"

#include <algorithm>

#include "viewer/about_texts.h"

namespace viewer {

namespace {

// An embedded licence file as text (about_texts.h stores its UTF-8 bytes unsigned).
std::string text(const unsigned char* bytes) { return std::string(reinterpret_cast<const char*>(bytes)); }

} // namespace

AboutIdentity aboutIdentity() {
    AboutIdentity a;
    a.product = "Wool Shader";
    a.summary = "Knitted yarn and wool fibre geometry from 3D models, for rendering.";
    a.copyright = "\xC2\xA9 2026 Nicolas Menzel";  // the copyright sign in UTF-8
    a.email = "drnicolasmenzel@gmail.com";
    a.portfolioUrl = "https://github.com/nicolasgfx/Portfolio";
    a.repositoryUrl = "https://github.com/nicolasgfx/woolshader";
    a.build = kBuildIdentity;
    return a;
}

std::vector<ThirdPartyNotice> thirdPartyNotices() {
    std::vector<ThirdPartyNotice> v;
    v.push_back({"Dear ImGui", "the interface", "MIT", "v1.91.8", "https://github.com/ocornut/imgui", "imgui",
                 text(kImguiLicenceText) +
                     "\nDear ImGui contains stb_truetype, stb_rect_pack and stb_textedit by Sean Barrett (MIT or public domain, "
                     "the text under stb below).\n"});
    v.push_back({"GLFW", "the window and its input", "zlib", "3.3.9", "https://github.com/glfw/glfw", "glfw", text(kGlfwLicenceText)});
    v.push_back({"Roboto Medium", "the interface's font, next to the program (from Dear ImGui's misc/fonts)", "Apache-2.0",
                 "the copy in Dear ImGui v1.91.8", "https://fonts.google.com/specimen/Roboto", "",
                 std::string("Roboto by Christian Robertson (Google), licensed under the Apache License 2.0, as Dear ImGui's "
                             "docs/FONTS.md states.\n\n") +
                     text(kOpenVdbLicenceText)});
    v.push_back({"cgltf", "reading GLB files", "MIT", "v1.14 (52c23814)", "https://github.com/jkuhlmann/cgltf", "cgltf", text(kCgltfLicenceText)});
    v.push_back({"stb", "decoding images (stb_image)", "MIT or public domain", "2c980bb5", "https://github.com/nothings/stb", "stb", text(kStbLicenceText)});
    v.push_back({"OpenVDB", "the import's voxel remesh", "Apache-2.0", "v13.0.0", "https://github.com/AcademySoftwareFoundation/openvdb",
                 "openvdb", text(kOpenVdbLicenceText)});
    v.push_back({"oneTBB", "OpenVDB's threads", "Apache-2.0", "v2022.3.0", "https://github.com/uxlfoundation/oneTBB", "tbb", text(kTbbLicenceText)});
    // The copyright line repeats the files' own notices word for word, as the MPL-2.0 asks.
    v.push_back({"SLIM solver", "the import's flatten", "MPL-2.0", "eight files from intern/slim of Blender 5.1.1",
                 "https://github.com/blender/blender/tree/v5.1.1/intern/slim/intern", "",
                 std::string("Copyright 2013 Alec Jacobson, 2016 Michael Rabinovich, 2023 Blender Authors (the files' own headers). "
                             "Derived from libigl; the method is Rabinovich, Poranne, Panozzo, Sorkine-Hornung 2017, Scalable "
                             "Locally Injective Mappings.\n\n"
                             "Wool Shader's copy, with its changes, is in third_party/slim/src of Wool Shader's source, "
                             "https://github.com/nicolasgfx/woolshader. These files and every change to them stay under the "
                             "Mozilla Public License 2.0; the rest of Wool Shader is MIT. The changes, from "
                             "third_party/slim/CHANGES.md:\n\n") +
                     text(kSlimChangesText) + "\n" + text(kMplLicenceText)});
    v.push_back({"Eigen", "SLIM and the import's flatten, compiled under the namespace name SlimEigen", "MPL-2.0",
                 "8a1083e9bf41", "https://gitlab.com/libeigen/eigen", "slim_eigen",
                 text(kEigenLicenceNote) + "\n\n" + text(kMplLicenceText)});
    // The stitch-mesh experiment's libraries (the Mark step's builds, through jobs and stitchmesh).
    v.push_back({"RIM (robust_hex_dominant_meshing)", "the stitch-mesh experiment's field-aligned remesher (Method B)", "MIT", "83dd574d",
                 "https://github.com/gaoxifeng/robust_hex_dominant_meshing", "rim",
                 std::string("Changed: patched at build time by third_party/rim/rim.cmake of Wool Shader's source "
                             "(third_party/rim/CHANGES.md).\n\n") +
                     text(kRimLicenceText) +
                     "\nRIM's dedge, subdivide and bvh units are part of Instant Meshes (Wenzel Jakob, Daniele Panozzo, Marco Tarini, "
                     "Olga Sorkine-Hornung), under its BSD licence:\n\n" +
                     text(kInstantMeshesLicenceText) + "\nRIM includes Tomas Moller's triangle-triangle test:\n\n" +
                     text(kMollerTriTriLicenceText)});
    v.push_back({"pcg32", "random numbers for RIM", "Apache-2.0", "0ef13e68", "https://github.com/wjakob/pcg32", "pcg32",
                 std::string("The notice at the top of pcg32.h (Wenzel Jakob; the PCG generator by Melissa O'Neill), then the Apache License "
                             "2.0.\n\n") +
                     text(kPcg32Notice).substr(0, text(kPcg32Notice).find("#ifndef")) + "\n" + text(kOpenVdbLicenceText)});
    v.push_back({"HiGHS", "the stitch-mesh experiment's integer programs (Method B's labels and directions)", "MIT", "364c83a5",
                 "https://github.com/ERGO-Code/HiGHS", "highs",
                 text(kHighsLicenceText) + "\nHiGHS compiles its bundled filereaderlp (MIT):\n\n" + text(kFilereaderlpLicenceText) +
                     "\nand its bundled pdqsort (Orson Peters, zlib licence).\n"});
    v.push_back({"NVIDIA CUDA runtime", "computing on the GPU; linked into the program", "NVIDIA CUDA Toolkit EULA",
                 std::string("CUDA Toolkit ") + kCudaVersion, "https://developer.nvidia.com/cuda-toolkit", "",
                 "The CUDA runtime library (cudart_static) is Copyright NVIDIA Corporation. It is distributed under the NVIDIA "
                 "CUDA Toolkit End User License Agreement: https://docs.nvidia.com/cuda/eula/\n\n"
                 "NVIDIA and CUDA are trademarks or registered trademarks of NVIDIA Corporation.\n"});
    return v;
}

std::string woolShaderLicenceText() { return text(kWoolShaderLicenceText); }

std::string otherLicencesStatement() {
    return "The MIT License above covers Wool Shader, except two parts of its source under their own licences:\n\n"
           "- blender/ppt_package_to_blend.py, the script of the Blender file export, is under GPL-3.0-or-later "
           "(blender/LICENSE), because it uses Blender's Python API, which Blender publishes under the GNU General Public "
           "License. It runs in Blender, a separate program, and is not compiled into Wool Shader (Blender and privacy).\n"
           "- third_party/slim/src, eight files of the SLIM solver from Blender's source, is under the Mozilla Public "
           "License 2.0 (Third-party software).\n\n"
           "The libraries in the program keep their own licences (Third-party software). The example models are not covered "
           "by the MIT License (Example models).";
}

std::vector<knit::ModelCredit> exampleModelCredits() {
    std::vector<knit::ModelCredit> v;
    for (const knit::SceneEntry& e : knit::sceneRegistry()) {
        if (!e.isPreset()) continue;
        // Presets of one model (the wool swatch and the blue swatch) credit it once.
        const bool listed = std::any_of(v.begin(), v.end(), [&](const knit::ModelCredit& c) { return c.title == e.credit.title; });
        if (!listed) v.push_back(e.credit);
    }
    return v;
}

std::string teddyStatement() {
    const knit::SceneEntry* teddy = knit::findScene("teddy");
    if (!teddy || !teddy->isPreset()) return "";
    const knit::ModelCredit& c = teddy->credit;
    return "The teddy bears are made from \"" + c.title + "\" by " + c.author + " (" + c.source + "), licensed under " + c.licence +
           " (" + c.licenceUrl +
           ").\n\n"
           "Wool Shader changed the model: it knits it as sewn pieces of yarn, adds fibres and gives it its own materials. The "
           "model's file is not exported as it is, but an export contains parts of it: the parts that are not knitted (the "
           "eyes, the noses and the bow) as the original meshes with their texture, and behind each knitted part its surface "
           "moved behind the yarn (the Filling). Under Fur, the whole mesh and its texture are exported under the fur.\n\n"
           "Each export carries CREDITS.txt with this credit. Keep it with the export when you share it: the licence requires "
           "the credit, a link to the licence and a note that the model was changed. Neither " +
           c.author + " nor Sketchfab endorses Wool Shader.";
}

std::string exampleModelsStatement() {
    return "Each example model keeps its own licence, named with its credit above. Yarn and fibres generated from a model are "
           "work made from it, so its licence and credit apply to them as well.\n\n"
           "The wool swatch is the author's own scene, under CC BY-NC 4.0: it may be shared with the credit, but not used "
           "commercially. Its needle's wood texture is \"Wood005\" from ambientCG (https://ambientcg.com/view?id=Wood005), "
           "CC0 1.0.\n\n"
           "Wool Shader adds no credit for a model you import. Its licence is yours to follow, and it applies to what you "
           "export from it.";
}

std::string blenderStatement() {
    return "Blender is a separate program and is not part of Wool Shader. Wool Shader uses it for one thing: writing a Blender "
           "file, the Export step's Blender format, when Blender 5.1 or later is installed. Everything else, the import "
           "included, works without it. Blender is free software under the GNU General Public License.\n\n"
           "For that export Wool Shader starts Blender with one script, blender/ppt_package_to_blend.py, the exported package "
           "folder and the output path. The script uses Blender's Python API, so it is released under GPL-3.0-or-later "
           "(blender/LICENSE); it is the only part of Wool Shader under the GPL. Blender itself is not compiled into Wool "
           "Shader or shipped with it.\n\n"
           "Blender is a registered trademark of the Blender Foundation. Wool Shader is not made, endorsed or supported by the "
           "Blender Foundation.";
}

std::string privacyStatement() {
    return "Wool Shader collects no data and makes no network connections. Models, settings and exports stay on this computer. "
           "The only programs it starts run on this computer too: Blender, for a Blender file export, and the photon path "
           "tracer for the Render preview. It talks to them through files, a local pipe and their standard input.";
}

} // namespace viewer
