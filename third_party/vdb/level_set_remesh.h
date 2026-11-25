#pragma once
// A closed triangle surface rebuilt as even quads through a narrow-band level set (Museth 2013, VDB): the surface is converted
// to a signed distance grid and its zero level is meshed again. This is the only file that sees OpenVDB; the interface takes
// plain arrays, so no OpenVDB or TBB header reaches its callers.
#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace knit {

struct LevelSetInput {
    std::vector<std::array<float, 3>> positions;
    std::vector<std::array<uint32_t, 3>> triangles;
    float voxelSize = 1;
};

// OpenVDB's arrays as they come: the points, then the quads and the triangles of the zero level.
struct LevelSetOutput {
    std::vector<std::array<float, 3>> positions;
    std::vector<std::array<uint32_t, 4>> quads;
    std::vector<std::array<uint32_t, 3>> triangles;
};

// The input through a level set of voxel `voxelSize` (OpenVDB's default half-band of 3 voxels) and back to a mesh at isovalue 0
// with adaptivity 0. False with OpenVDB's message if it throws.
[[nodiscard]] bool remeshThroughLevelSet(const LevelSetInput& in, LevelSetOutput& out, std::string* err);

} // namespace knit
