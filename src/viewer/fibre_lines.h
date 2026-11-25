// Fibre lines in the viewer. The shown result's DisplayLines live in GL buffer objects (uploaded
// once per result; the colours rebuilt only when the fibre colours change) instead of client
// arrays every frame, and are drawn chunk by chunk (knit::FibreChunk):
//   1. a chunk outside the view frustum is skipped;
//   2. with a segment budget (while the camera moves) every visible chunk draws the same share
//      of its fibres, a prefix of its random order: an even thinning.
// The model's far side is not culled: a test on the fibres' root normals lost hair hanging past
// the silhouette (shiba, teddy) and skipped almost nothing on wool.
#pragma once

#include <memory>
#include <vector>

#include "viewer/camera.h"
#include "viewer/gl_functions.h"
#include "jobs/built_scene.h"

namespace viewer {

// Linear colour -> the sRGB value the fixed-function display shows (object colours are linear
// albedo, as the renderer takes them).
float toSrgb(float linear);

// Where a fibre population's display colour comes from: a fixed colour (fur: the fur colour; wool in
// one colour: yellow along the yarn, blue stray, the two populations told apart), its part's colour,
// or its own (the texture colour at its stitch or root).
enum class FibreColourFrom { Fixed, Part, Own };

// What the fibre colours depend on besides the lines (the buffers are recoloured when it changes).
struct FibreColours {
    FibreColourFrom fur = FibreColourFrom::Fixed, oriented = FibreColourFrom::Fixed, stray = FibreColourFrom::Fixed;
    Vec3 furColour;
};
bool sameColours(const FibreColours& a, const FibreColours& b);
// The settings' colours, as the export writes the fibres.
FibreColours fibreColoursOf(const knit::SceneSettings& s);
// The same colours with the two wool populations told apart: yellow along the yarn, blue stray,
// whatever colour they export in. Hair keeps its own.
FibreColours populationColours(FibreColours colours);
// The display colour of vertex v of fibre lines (sRGB, as the meshes are shown).
Vec3 fibreVertexColour(const knit::DisplayLines& lines, size_t v, const FibreColours& colours);

// The camera as the frustum test takes it.
struct FibreCamera {
    Mat4 clip{};           // projection x view: scene mm -> clip space
    bool cull = true;      // false: every chunk drawn (a model transform the matrices do not hold)
};
FibreCamera fibreCamera(const Mat4& projection, const Mat4& view);

// Whether any part of the chunk's bounding box can be inside the view frustum.
bool chunkInFrustum(const knit::FibreChunk& chunk, const Mat4& clip);

class FibreLines {
public:
    FibreLines() = default;
    FibreLines(const FibreLines&) = delete;
    FibreLines& operator=(const FibreLines&) = delete;
    ~FibreLines();

    // Draws `lines` (held by `owner`, kept alive while its buffers exist) as GL_LINES with the
    // matrices currently loaded, depth-tested. segmentBudget > 0: at most about that many segments
    // (every visible chunk thinned alike); 0: all visible ones. Returns the segments drawn.
    long long draw(const std::shared_ptr<const knit::BuiltScene>& owner, const knit::DisplayLines& lines,
                   const FibreColours& colours, const FibreCamera& camera, long long segmentBudget);

    // Frees the buffers unless they hold the lines of `shown`.
    void releaseUnless(const knit::BuiltScene* shown);

private:
    void upload(const knit::DisplayLines& lines);
    void recolour(const knit::DisplayLines& lines, const FibreColours& colours);
    void release();

    std::shared_ptr<const knit::BuiltScene> owner_;
    const knit::DisplayLines* lines_ = nullptr;
    GLuint positions_ = 0, colours_ = 0;  // buffer objects: float xyz, RGBA8 per vertex
    bool coloured_ = false;
    FibreColours colourKey_;              // the colours rgba_ holds
    std::vector<const knit::FibreChunk*> visible_;  // the frame's chunks in view (kept: no allocation per frame)
};

} // namespace viewer
