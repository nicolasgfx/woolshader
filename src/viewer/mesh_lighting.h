// How the viewer's meshes are lit, Flat or Studio.
//
// Studio lighting: the viewer's shaded meshes with depth cues the fixed-function headlight lacks.
//
//   1. normal / depth pass: view-space normal and distance of every pixel (an off-screen buffer);
//   2. screen-space ambient occlusion (SSAO): how much of each pixel's hemisphere nearby surfaces
//      cover, at two radii: the yarn diameter (the grooves between neighbouring yarns) and four
//      times that (an ear against the head, an arm against the body);
//   3. a 4 x 4 depth-aware blur of the occlusion;
//   4. per-pixel lighting from a camera-relative rig: a key light from the upper left (about 40
//      degrees off the view axis, so tubes show a lit side and a shaded side), a weak fill from
//      the right, and a sky / ground hemisphere (world +z up) that the occlusion darkens.
//
// GL 2.1 with GLSL 1.20 shaders; the entry points above 1.1 come from gl_functions.h.
//
// Flat lighting is the fixed-function headlight, also the fallback without shaders. Both draw the same
// MeshLooks from buffer objects uploaded once per mesh (renderer.cpp).
#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "viewer/gl_functions.h"
#include "model/vec.h"

namespace viewer {

// A mesh on the GPU (Renderer owns them).
struct MeshBuffers {
    GLuint vertices = 0, indices = 0;
    // Byte offsets in `vertices` (-1: absent).
    long long normals = -1, colours = -1, plyUv = -1, albedoUv = -1;
    GLsizei indexCount = 0;
    long long bytes = 0;  // both buffers' sizes (Renderer::residentBytes)
};

enum class MeshTexture { None, Albedo, Ply };
enum class MeshFill { Solid, Wireframe };

// A GLSL 1.20 shader stage compiled, and a program linked from two (0 with the compiler's log in *err): the
// Studio lighting's passes and the Fur step's discs (fur_discs.h).
GLuint compileShader(GLenum type, const char* source, const char* name, std::string* err);
GLuint linkProgram(const char* vertex, const char* fragment, const char* name, std::string* err);

// How one mesh is drawn.
struct MeshLook {
    const MeshBuffers* buffers = nullptr;
    float color[3] = {0.7f, 0.7f, 0.7f};  // sRGB, as the fixed-function display takes it
    bool vertexColour = false;              // the mesh's vertex colours instead of color
    MeshTexture texture = MeshTexture::None;
    GLuint textureId = 0;                   // the albedo texture or the ply texture
    // Yarn tubes: the normals point out of the tube, used as they are (buildTubeMesh's triangles
    // wind clockwise seen from outside, so the face orientation would turn them inwards). Other
    // meshes: the normal is turned towards the viewer (either winding; open surfaces).
    bool outwardNormals = false;
    float alpha = 1;  // below 1: drawn see-through after the opaque meshes (drawMeshesSeeThrough)
};

// Sets the look's arrays up from its buffers; unbindMesh leaves no buffer bound for ImGui's client
// arrays.
void bindMesh(const MeshBuffers& b, const MeshLook& look);
void unbindMesh();

// With the matrices currently loaded: the fixed-function headlight.
void drawMeshesFlat(const std::vector<MeshLook>& looks, MeshFill fill);
// The same, blended by each look's alpha without writing depth (after the opaque meshes; a few pieces around the
// Mark step's isolated one, unsorted: where they overlap, the order shows).
void drawMeshesSeeThrough(const std::vector<MeshLook>& looks);

// The frame: the viewport (framebuffer pixels) and the camera's projection as the occlusion pass
// reconstructs view positions from it.
struct StudioView {
    int x = 0, y = 0, width = 1, height = 1;
    double tanHalfFov = 0.3;        // tan(vertical field of view / 2)
    double farDistance = 1000;      // [scene mm] the far plane (stored distances are fractions of it)
    knit::Vec3 skyInView{0, 1, 0};  // world +z in view coordinates (the hemisphere light's sky)
    double occlusionRadius = 1;     // [scene mm] the near radius (the far one is 4 times it)
    double occlusionStrength = 1;   // 0: no occlusion, 1: full
    int clipX = 0, clipWidth = 1;   // the framebuffer columns the lit pass writes (one half while comparing)
};

class StudioLighting {
public:
    StudioLighting() = default;
    StudioLighting(const StudioLighting&) = delete;
    StudioLighting& operator=(const StudioLighting&) = delete;
    ~StudioLighting();

    // Compiles the shaders (needs the GL entry points); false states why and leaves Flat lighting
    // the only choice.
    [[nodiscard]] bool init(std::string* err);
    bool ready() const { return ready_; }

    // Draws the meshes with the matrices currently loaded (GL_PROJECTION, GL_MODELVIEW) into the
    // default framebuffer's viewport, depth-tested, writing only the columns clipX .. clipX + clipWidth (the
    // scissor test is left on); it restores the fixed-function state (no
    // program, buffer, framebuffer bound; texture unit 0 active). False: nothing drawn (the
    // off-screen buffers failed at this size); the caller draws Flat.
    [[nodiscard]] bool draw(const StudioView& view, const std::vector<MeshLook>& looks);

private:
    struct Program {
        GLuint id = 0;
        std::map<std::string, GLint> uniforms;
        GLint at(const std::string& name);
    };
    struct Target {  // an off-screen colour texture (with depth for the normal / depth pass)
        GLuint framebuffer = 0, texture = 0, depth = 0;
    };

    [[nodiscard]] bool resizeTargets(int width, int height);
    void drawNormalDepth(const std::vector<MeshLook>& looks);
    void computeOcclusion(const StudioView& view);
    void blurOcclusion(const StudioView& view);
    void drawLit(const StudioView& view, const std::vector<MeshLook>& looks);
    void releaseTargets();

    bool ready_ = false;
    Program normalDepth_, occlusion_, blur_, lit_;
    Target normalDepthTarget_, occlusionTarget_, blurTarget_;
    int width_ = 0, height_ = 0;
    double farDistance_ = 1000;  // the frame's StudioView::farDistance (the normal / depth pass)
    GLuint rotations_ = 0;
    std::vector<float> kernel_;  // 16 x 3
};

// The sample pattern of the screen-space ambient occlusion, GL-free so its properties are tested
// (tests_ui.cpp).
//
// n sample offsets in the unit hemisphere around +z (the surface normal at the shaded point).
// Sample i is scaled by 0.1 + 0.9 (i / n)^2, so most samples lie close to the point and near
// occluders (the groove between two touching yarns) count most; the common SSAO kernel, e.g.
// J. de Vries, LearnOpenGL "SSAO". Every sample keeps z >= 0.15 |p|: offsets in the tangent plane
// would test the shaded surface against itself. Deterministic for a seed.
std::vector<knit::Vec3> ssaoKernel(int n, uint32_t seed);

// The 4 x 4 per-pixel kernel rotations: unit vectors in the xy plane. The 4 x 4 blur after the
// occlusion pass averages over one tile, so the rotations trade banding for no visible noise.
std::vector<knit::Vec3> ssaoRotations(uint32_t seed);

} // namespace viewer
