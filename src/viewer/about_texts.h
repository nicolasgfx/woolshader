#pragma once
// The texts the About window shows word for word, compiled into the program at build time
// (src/viewer/CMakeLists.txt): the licence files of this project and of the libraries it is built
// with, read from the checkout and from the dependencies' fetched sources (cmake/embed_texts.cmake),
// and the commit it is built from with the CUDA Toolkit version (cmake/build_identity.cmake). Each is
// a 0-terminated UTF-8 text with "\n" line ends. about.h puts them together with their names and
// credits.

namespace viewer {

extern const unsigned char kWoolShaderLicenceText[];  // LICENSE: MIT
extern const unsigned char kImguiLicenceText[];       // Dear ImGui: MIT
extern const unsigned char kGlfwLicenceText[];        // GLFW: zlib
extern const unsigned char kCgltfLicenceText[];       // cgltf: MIT
extern const unsigned char kStbLicenceText[];         // stb: MIT or public domain
extern const unsigned char kOpenVdbLicenceText[];     // OpenVDB: Apache-2.0
extern const unsigned char kTbbLicenceText[];         // oneTBB: Apache-2.0
extern const unsigned char kEigenLicenceNote[];       // Eigen's COPYING.README: which of its licences apply where
extern const unsigned char kMplLicenceText[];         // Eigen's COPYING.MPL2: MPL-2.0, also the licence of the SLIM files
extern const unsigned char kRimLicenceText[];         // robust_hex_dominant_meshing (RIM): MIT
extern const unsigned char kHighsLicenceText[];       // HiGHS: MIT
extern const unsigned char kFilereaderlpLicenceText[]; // HiGHS's bundled filereaderlp: MIT
extern const unsigned char kInstantMeshesLicenceText[]; // Instant Meshes, compiled as part of RIM: BSD (third_party/rim)
extern const unsigned char kMollerTriTriLicenceText[];  // Moller's triangle-triangle test, included by RIM: MIT (third_party/rim)
extern const unsigned char kPcg32Notice[];            // pcg32.h, whose header comment carries its Apache-2.0 notice
extern const unsigned char kSlimChangesText[];        // third_party/slim/CHANGES.md: every change made to the SLIM files
extern const char kBuildIdentity[];          // git describe --always --dirty at build time, or "unknown"
extern const char kCudaVersion[];            // the CUDA Toolkit the runtime library comes from, e.g. "13.4.48"

} // namespace viewer
