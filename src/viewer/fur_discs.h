#pragma once
// The Fur step's disc preview on the GPU (jobs/fur_preview.h builds the data). Each loop disc is one triangle in its
// own plane, circumscribed about the disc's circle; the fragment stage discards outside the circle (x^2 + y^2 > 1 in
// the disc's local coordinates) before anything is written, so the depth is the disc's own: no impostor, correct
// against the meshes and each other, also for a disc larger than the screen or crossing the near plane (clipped,
// not culled). Both sides are lit alike (|n . l|): an edge-on or back view keeps every disc. The loops' direction
// marks and the fine fibres are lines. A fibre is drawn when its budget key lies below the key limit (1: all): an
// even thinning while the camera moves, as the fibre lines thin (fibre_lines.h).
//
// Without the shader (an old driver) the discs cannot be cut to circles, so they are not drawn as solid triangles:
// only the marks and the fine fibres show, and diagnosticNote() says why.
#include <cstdint>
#include <memory>
#include <string>

#include "jobs/fur_preview.h"
#include "viewer/camera.h"
#include "viewer/gl_functions.h"

namespace viewer {

class FurDiscs {
public:
    FurDiscs() = default;
    FurDiscs(const FurDiscs&) = delete;
    FurDiscs& operator=(const FurDiscs&) = delete;
    ~FurDiscs();

    // Compiles the disc shader (needs the context, Renderer::init). False: discs fall back to their marks.
    bool init(std::string* why);
    // Draws the preview with the camera's matrices loaded (uploading it first when it is new); keyLimit in (0, 1].
    // Returns the milliseconds spent uploading this frame (0: already on the GPU).
    double draw(const std::shared_ptr<const knit::FurPreview>& preview, float keyLimit);
    // Why the discs are not drawn as circles ("" while the shader works).
    const std::string& diagnosticNote() const { return note_; }
    long long residentBytes() const { return bytes_; }

private:
    void upload(const knit::FurPreview& p);
    void release();

    GLuint program_ = 0, lineProgram_ = 0;
    GLint keyLimit_ = -1, lineKeyLimit_ = -1;
    GLuint buffers_[10] = {};  // disc pos, normal, local, colour; mark pos, key, colour; needle pos, key, colour
    std::shared_ptr<const knit::FurPreview> shown_;  // the uploaded preview, kept alive with its buffers
    long long bytes_ = 0;
    std::string note_;
};

} // namespace viewer
