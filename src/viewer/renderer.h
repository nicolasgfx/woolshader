#pragma once
// One frame of the 3D view: the meshes under the chosen lighting, then the fibres. Everything on the
// GPU is here, under one lifetime rule: it belongs to a scene that may be on screen (MeshOwners), and
// a scene's meshes and textures stay uploaded while that scene is listed.
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "viewer/camera.h"
#include "viewer/fibre_lines.h"
#include "viewer/fur_discs.h"
#include "viewer/gl_functions.h"
#include "viewer/mesh_lighting.h"
#include "viewer/view_content.h"

namespace viewer {

// The viewport in framebuffer pixels (GL: origin bottom left).
struct ViewportPixels {
    int x = 0, y = 0, width = 1, height = 1;
};
// The interface's viewportRect (in the interface's units, origin top left) in pixels.
ViewportPixels viewportPixels(const SceneUi& ui, int framebufferWidth, int framebufferHeight, float displayWidth, float displayHeight);

struct FrameTarget {
    int width = 1, height = 1;  // the framebuffer
    ViewportPixels viewport;
    bool cameraMoving = false;  // the fibres thin
};

class Renderer {
public:
    Renderer() = default;
    Renderer(const Renderer&) = delete;
    Renderer& operator=(const Renderer&) = delete;
    ~Renderer();

    // Loads the GL entry points and compiles the shaders (needs the context). False: Studio cannot
    // run, *why says so; Flat still works.
    [[nodiscard]] bool init(std::string* why);
    // The frame: the shown content over the whole viewport, or the source in its left half and the shown
    // content in its right half. Both halves take the whole viewport's camera and projection.
    void draw(const SceneUi& ui, const ViewportFrame& frame, const Camera& camera, const FrameTarget& target);
    void releaseUnused(const SceneUi& ui);
    // The Render preview's newest frame into its texture (ui.renderPreview), when it is new.
    void uploadRenderPreview(SceneUi& ui);
    long long fibreSegmentsDrawn() const { return fibreSegments_; }
    const std::string& discNote() const { return furDiscs_.diagnosticNote(); }  // why discs show as marks only ("" none)
    long long residentBytes() const;          // the uploaded mesh buffers and albedo textures (--bench)
    double uploadMs() const { return uploadMs_; }  // the time spent uploading them, all frames (--bench)

private:
    // The framebuffer columns one pass may write (the scissor): the whole viewport, or one half.
    struct ViewColumns {
        int x = 0, width = 1;
    };
    void drawContent(const SceneUi& ui, const ViewContent& c, const Camera& camera, const FrameTarget& target, ViewColumns columns);
    const MeshBuffers* buffersOf(const knit::DisplayMesh& d, const MeshOwners& owners);
    GLuint albedoOf(const ViewContent& c, const knit::DisplayMesh& d);
    std::vector<MeshLook> meshLooks(const SceneUi& ui, const ViewContent& c, const MeshOwners& owners);
    void drawFocusOutline(const SceneUi& ui, const ViewContent& c, const MeshOwners& owners);
    void keepPlyTexture(const std::shared_ptr<const knit::TexturePreview>& preview);
    void release(const MeshOwners& owners);

    MeshOwners kept_;  // whatever is uploaded belongs to one of these, which it keeps alive
    std::map<const knit::DisplayMesh*, MeshBuffers> buffers_;
    std::map<const knit::DisplayTexture*, GLuint> albedo_;
    std::map<GLuint, long long> albedoBytes_;  // each albedo texture's size with its mipmaps
    GLuint plyTexture_ = 0;
    GLuint previewTexture_ = 0;  // the Render preview's frame
    std::shared_ptr<const knit::TexturePreview> plyPreview_;
    StudioLighting studio_;
    FibreLines fibres_;
    FurDiscs furDiscs_;
    long long fibreSegments_ = 0;
    double uploadMs_ = 0;
};

// The back buffer (the whole window, interface included) as a PNG.
[[nodiscard]] bool writeBackBufferPng(const std::string& path, int width, int height);

} // namespace viewer
