// The disc preview on the GPU (fur_discs.h has the method).
#include "viewer/fur_discs.h"

#include <chrono>

#include "viewer/mesh_lighting.h"

namespace viewer {

namespace {

// A vertex whose fibre's budget key reaches the limit is moved outside the clip volume: it is thinned away.
const char* kDiscVertex = R"(
#version 120
uniform float keyLimit;
varying vec3 local;
varying vec3 normal;
varying vec4 colour;
void main() {
    local = gl_MultiTexCoord0.xyz;
    normal = gl_NormalMatrix * gl_Normal;
    colour = gl_Color;
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    if (gl_MultiTexCoord0.z >= keyLimit) gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
}
)";

// Outside the disc's circle nothing is written (no colour, no depth); inside, the Studio key light on both sides alike.
const char* kDiscFragment = R"(
#version 120
varying vec3 local;
varying vec3 normal;
varying vec4 colour;
void main() {
    if (dot(local.xy, local.xy) > 1.0) discard;
    float light = abs(dot(normalize(normal), vec3(-0.4215, 0.4917, 0.7660)));
    gl_FragColor = vec4(colour.rgb * (0.45 + 0.55 * light), 1.0);
}
)";

const char* kLineVertex = R"(
#version 120
uniform float keyLimit;
varying vec4 colour;
void main() {
    colour = gl_Color;
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
    if (gl_MultiTexCoord0.x >= keyLimit) gl_Position = vec4(0.0, 0.0, 2.0, 1.0);
}
)";

const char* kLineFragment = R"(
#version 120
varying vec4 colour;
void main() { gl_FragColor = colour; }
)";

enum Buffer { kDiscPos, kDiscNormal, kDiscLocal, kDiscColour, kMarkPos, kMarkKey, kMarkColour, kNeedlePos, kNeedleKey, kNeedleColour, kBuffers };

// Lines from three buffers (positions, keys as texture coordinate s, colours), thinned by the bound program.
void drawLines(const GLuint* buffers, Buffer pos, Buffer key, Buffer colour, long long count, bool keyed) {
    if (count <= 0) return;
    gl::BindBuffer(gl::kArrayBuffer, buffers[pos]);
    glVertexPointer(3, GL_FLOAT, 0, nullptr);
    if (keyed) {
        gl::BindBuffer(gl::kArrayBuffer, buffers[key]);
        glTexCoordPointer(1, GL_FLOAT, 0, nullptr);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
    }
    gl::BindBuffer(gl::kArrayBuffer, buffers[colour]);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, nullptr);
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(count * 2));
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
}

} // namespace

FurDiscs::~FurDiscs() {
    release();
    if (program_) gl::DeleteProgram(program_);
    if (lineProgram_) gl::DeleteProgram(lineProgram_);
}

bool FurDiscs::init(std::string* why) {
    std::string err;
    program_ = linkProgram(kDiscVertex, kDiscFragment, "disc shader", &err);
    lineProgram_ = program_ ? linkProgram(kLineVertex, kLineFragment, "disc line shader", &err) : 0;
    if (!program_ || !lineProgram_) {
        note_ = "The disc shader is not available (" + err + "): the discs show as their direction marks only.";
        if (why) *why = note_;
        return false;
    }
    keyLimit_ = gl::GetUniformLocation(program_, "keyLimit");
    lineKeyLimit_ = gl::GetUniformLocation(lineProgram_, "keyLimit");
    return true;
}

void FurDiscs::release() {
    if (buffers_[0]) gl::DeleteBuffers(kBuffers, buffers_);
    for (GLuint& b : buffers_) b = 0;
    bytes_ = 0;
    shown_.reset();
}

void FurDiscs::upload(const knit::FurPreview& p) {
    if (!buffers_[0]) gl::GenBuffers(kBuffers, buffers_);
    const auto put = [&](Buffer b, const void* data, size_t bytes) {
        gl::BindBuffer(gl::kArrayBuffer, buffers_[b]);
        gl::BufferData(gl::kArrayBuffer, static_cast<std::ptrdiff_t>(bytes), bytes ? data : nullptr, gl::kStaticDraw);
        bytes_ += static_cast<long long>(bytes);
    };
    bytes_ = 0;
    put(kDiscPos, p.discPos.data(), p.discPos.size() * sizeof(float));
    put(kDiscNormal, p.discNormal.data(), p.discNormal.size() * sizeof(float));
    put(kDiscLocal, p.discLocal.data(), p.discLocal.size() * sizeof(float));
    put(kDiscColour, p.discColour.data(), p.discColour.size() * sizeof(uint32_t));
    put(kMarkPos, p.markPos.data(), p.markPos.size() * sizeof(float));
    put(kMarkKey, p.markKey.data(), p.markKey.size() * sizeof(float));
    put(kMarkColour, p.markColour.data(), p.markColour.size() * sizeof(uint32_t));
    put(kNeedlePos, p.needlePos.data(), p.needlePos.size() * sizeof(float));
    put(kNeedleKey, p.needleKey.data(), p.needleKey.size() * sizeof(float));
    put(kNeedleColour, p.needleColour.data(), p.needleColour.size() * sizeof(uint32_t));
    gl::BindBuffer(gl::kArrayBuffer, 0);
}

double FurDiscs::draw(const std::shared_ptr<const knit::FurPreview>& preview, float keyLimit) {
    if (!preview) return 0;
    double uploadMs = 0;
    if (preview != shown_) {
        const auto t0 = std::chrono::steady_clock::now();
        upload(*preview);
        glFinish();  // the upload's real cost, for the latency notes (once per new preview)
        uploadMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        shown_ = preview;
    }
    const knit::FurPreview& p = *preview;
    glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glDisable(GL_TEXTURE_2D);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    if (program_ && p.loops > 0) {
        gl::UseProgram(program_);
        gl::Uniform1f(keyLimit_, keyLimit);
        gl::BindBuffer(gl::kArrayBuffer, buffers_[kDiscPos]);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        gl::BindBuffer(gl::kArrayBuffer, buffers_[kDiscNormal]);
        glNormalPointer(GL_FLOAT, 0, nullptr);
        glEnableClientState(GL_NORMAL_ARRAY);
        gl::BindBuffer(gl::kArrayBuffer, buffers_[kDiscLocal]);
        glTexCoordPointer(3, GL_FLOAT, 0, nullptr);
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        gl::BindBuffer(gl::kArrayBuffer, buffers_[kDiscColour]);
        glColorPointer(4, GL_UNSIGNED_BYTE, 0, nullptr);
        glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(p.loops * 3));
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    }
    if (lineProgram_) {
        gl::UseProgram(lineProgram_);
        gl::Uniform1f(lineKeyLimit_, keyLimit);
    }
    drawLines(buffers_, kMarkPos, kMarkKey, kMarkColour, p.loops, lineProgram_ != 0);
    drawLines(buffers_, kNeedlePos, kNeedleKey, kNeedleColour, p.fibres, lineProgram_ != 0);
    gl::UseProgram(0);
    gl::BindBuffer(gl::kArrayBuffer, 0);
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glPopAttrib();
    return uploadMs;
}

} // namespace viewer
