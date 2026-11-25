// Flat lighting, and Studio: normal / depth pass, SSAO, blur, per-pixel lit meshes (mesh_lighting.h).
#include "viewer/mesh_lighting.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

namespace viewer {

namespace {

constexpr int kKernelSize = 16;

// --- shaders (GLSL 1.20: the fixed-function attributes and matrices of a GL 2.1 context)

// The meshes' vertex stage, shared by the normal / depth and the lit pass: view-space position
// and normal, the vertex colour and the texture coordinates.
const char* kMeshVertex = R"(
#version 120
varying vec3 viewPosition;
varying vec3 viewNormal;
varying vec4 colour;
varying vec2 texCoord;
void main() {
    vec4 p = gl_ModelViewMatrix * gl_Vertex;
    viewPosition = p.xyz;
    viewNormal = gl_NormalMatrix * gl_Normal;
    colour = gl_Color;
    texCoord = gl_MultiTexCoord0.xy;
    gl_Position = gl_ProjectionMatrix * p;
}
)";

// Normal (xyz) and distance along the view axis (w; 0 = no surface) of the nearest surface,
// stored in [0, 1] (n / 2 + 1/2, distance / far plane): a GL 2.1 context may clamp the colour
// written even to a float target. The normal as MeshLook::outwardNormals states.
const char* kNormalDepthFragment = R"(
#version 120
uniform float farDistance;
uniform int outwardNormals;
varying vec3 viewPosition;
varying vec3 viewNormal;
void main() {
    vec3 n = normalize(viewNormal);
    if (outwardNormals == 0 && dot(n, viewPosition) > 0.0) n = -n;
    gl_FragColor = vec4(n * 0.5 + 0.5, -viewPosition.z / farDistance);
}
)";

// A full-screen quad (vertices in clip space already).
const char* kQuadVertex = R"(
#version 120
varying vec2 uv;
void main() {
    uv = gl_Vertex.xy * 0.5 + 0.5;
    gl_Position = gl_Vertex;
}
)";

// Screen-space ambient occlusion (Crytek 2007; normal-oriented hemisphere as in LearnOpenGL
// "SSAO"): each kernel sample around the pixel's view position is projected to the screen; it is
// occluded when the surface stored there lies in front of it. The range check fades occluders
// far in front (a foreground object is not a groove). Two radii: the yarn's grooves (near) and
// the contact between parts (far, counted at half weight).
const char* kOcclusionFragment = R"(
#version 120
uniform sampler2D normalDepth;
uniform sampler2D rotations;
uniform vec3 kernel[16];
uniform vec2 tanHalf;
uniform float radiusNear;
uniform float radiusFar;
uniform float farDistance;
varying vec2 uv;

vec3 viewPositionAt(vec2 t, float depth) { return vec3((t * 2.0 - 1.0) * tanHalf * depth, -depth); }

float occlusion(vec3 p, mat3 tbn, float radius) {
    float occluded = 0.0;
    for (int i = 0; i < 16; ++i) {
        vec3 s = p + tbn * (kernel[i] * radius);
        vec2 t = s.xy / (-s.z * tanHalf) * 0.5 + 0.5;
        if (t.x <= 0.0 || t.y <= 0.0 || t.x >= 1.0 || t.y >= 1.0) continue;
        float d = texture2D(normalDepth, t).w * farDistance;
        if (d <= 0.0) continue;
        if (-s.z - d > 0.02 * radius) occluded += smoothstep(0.0, 1.0, radius / abs(-p.z - d));
    }
    return 1.0 - occluded / 16.0;
}

void main() {
    vec4 nd = texture2D(normalDepth, uv);
    if (nd.w <= 0.0) {
        gl_FragColor = vec4(1.0);
        return;
    }
    vec3 n = normalize(nd.xyz * 2.0 - 1.0);
    vec3 p = viewPositionAt(uv, nd.w * farDistance);
    vec3 r = texture2D(rotations, gl_FragCoord.xy / 4.0).xyz * 2.0 - 1.0;
    vec3 t = r - n * dot(r, n);
    if (dot(t, t) < 1e-6) t = vec3(-r.y, r.x, 0.0);
    t = normalize(t);
    mat3 tbn = mat3(t, cross(n, t), n);
    float near = occlusion(p, tbn, radiusNear);
    float far = occlusion(p, tbn, radiusFar);
    gl_FragColor = vec4(vec3(near * mix(1.0, far, 0.5)), 1.0);
}
)";

// The mean over the 4 x 4 pixels of one rotation tile, of the pixels at about the same depth
// (no occlusion bleeding across silhouettes).
const char* kBlurFragment = R"(
#version 120
uniform sampler2D occlusion;
uniform sampler2D normalDepth;
uniform vec2 texel;
uniform float depthRange;
uniform float farDistance;
varying vec2 uv;
void main() {
    float d0 = texture2D(normalDepth, uv).w * farDistance;
    float sum = 0.0, weight = 0.0;
    for (int y = -2; y < 2; ++y)
        for (int x = -2; x < 2; ++x) {
            vec2 t = uv + vec2(float(x), float(y)) * texel;
            float d = texture2D(normalDepth, t).w * farDistance;
            if (d <= 0.0 || abs(d - d0) > depthRange) continue;
            sum += texture2D(occlusion, t).r;
            weight += 1.0;
        }
    gl_FragColor = vec4(vec3(weight > 0.0 ? sum / weight : 1.0), 1.0);
}
)";

// Per-pixel lighting in linear colour (the colours and textures are sRGB): a sky / ground
// hemisphere darkened by the occlusion, a wrapped-diffuse key light (wrap 0.2: the terminator
// stays soft on thin tubes) that the occlusion darkens at half strength, a weak fill and a broad,
// faint Blinn highlight (wool is not glossy).
const char* kLitFragment = R"(
#version 120
uniform sampler2D surfaceTexture;
uniform int useTexture;
uniform sampler2D occlusion;
uniform vec2 viewportOrigin;
uniform vec2 viewportSize;
uniform float occlusionStrength;
uniform vec3 sky;
uniform vec3 keyDirection;
uniform vec3 fillDirection;
uniform int outwardNormals;
varying vec3 viewPosition;
varying vec3 viewNormal;
varying vec4 colour;
varying vec2 texCoord;

const vec3 skyColour = vec3(0.30, 0.32, 0.35);
const vec3 groundColour = vec3(0.10, 0.09, 0.08);
const vec3 keyColour = vec3(0.82, 0.79, 0.74);
const vec3 fillColour = vec3(0.10, 0.11, 0.12);

void main() {
    vec3 base = colour.rgb;
    if (useTexture != 0) base *= texture2D(surfaceTexture, texCoord).rgb;
    vec3 albedo = pow(base, vec3(2.2));
    vec3 n = normalize(viewNormal);
    if (outwardNormals == 0 && dot(n, viewPosition) > 0.0) n = -n;
    float ao = mix(1.0, texture2D(occlusion, (gl_FragCoord.xy - viewportOrigin) / viewportSize).r, occlusionStrength);
    vec3 ambient = mix(groundColour, skyColour, 0.5 + 0.5 * dot(n, sky)) * ao;
    float nl = dot(n, keyDirection);
    vec3 key = keyColour * clamp((nl + 0.2) / 1.2, 0.0, 1.0) * mix(1.0, ao, 0.5);
    vec3 fill = fillColour * max(dot(n, fillDirection), 0.0);
    vec3 h = normalize(keyDirection + normalize(-viewPosition));
    vec3 highlight = nl > 0.0 ? keyColour * 0.04 * pow(max(dot(n, h), 0.0), 16.0) : vec3(0.0);
    vec3 c = albedo * (ambient + key + fill) + highlight * ao;
    gl_FragColor = vec4(pow(c, vec3(1.0 / 2.2)), 1.0);
}
)";

// The rig in view coordinates (x right, y up, z towards the viewer): the key from the upper left,
// 40 degrees off the view axis; the fill from the lower right.
constexpr float kKeyDirection[3] = {-0.4215f, 0.4917f, 0.7660f};
constexpr float kFillDirection[3] = {0.6247f, -0.2082f, 0.7526f};

} // namespace

GLuint compileShader(GLenum type, const char* source, const char* name, std::string* err) {
    const GLuint s = gl::CreateShader(type);
    gl::ShaderSource(s, 1, &source, nullptr);
    gl::CompileShader(s);
    GLint ok = 0;
    gl::GetShaderiv(s, gl::kCompileStatus, &ok);
    if (ok) return s;
    char log[2048] = {};
    gl::GetShaderInfoLog(s, sizeof log, nullptr, log);
    *err = std::string(name) + ": " + log;
    gl::DeleteShader(s);
    return 0;
}

GLuint linkProgram(const char* vertex, const char* fragment, const char* name, std::string* err) {
    const GLuint vs = compileShader(gl::kVertexShader, vertex, name, err);
    if (!vs) return 0;
    const GLuint fs = compileShader(gl::kFragmentShader, fragment, name, err);
    if (!fs) {
        gl::DeleteShader(vs);
        return 0;
    }
    const GLuint p = gl::CreateProgram();
    gl::AttachShader(p, vs);
    gl::AttachShader(p, fs);
    gl::LinkProgram(p);
    gl::DeleteShader(vs);  // flagged; freed with the program
    gl::DeleteShader(fs);
    GLint ok = 0;
    gl::GetProgramiv(p, gl::kLinkStatus, &ok);
    if (ok) return p;
    char log[2048] = {};
    gl::GetProgramInfoLog(p, sizeof log, nullptr, log);
    *err = std::string(name) + ": " + log;
    gl::DeleteProgram(p);
    return 0;
}

namespace {

void drawFullScreenQuad() {
    glBegin(GL_QUADS);
    glVertex2f(-1, -1);
    glVertex2f(1, -1);
    glVertex2f(1, 1);
    glVertex2f(-1, 1);
    glEnd();
}

} // namespace

// --- the mesh arrays, and Flat

void bindMesh(const MeshBuffers& b, const MeshLook& look) {
    gl::BindBuffer(gl::kArrayBuffer, b.vertices);
    gl::BindBuffer(gl::kElementArrayBuffer, b.indices);
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, nullptr);
    if (b.normals >= 0) {
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, reinterpret_cast<const void*>(b.normals));
    }
    if (look.vertexColour && b.colours >= 0) {
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(3, GL_FLOAT, 0, reinterpret_cast<const void*>(b.colours));
    } else {
        glColor4f(look.color[0], look.color[1], look.color[2], look.alpha);
    }
    const long long uv = look.texture == MeshTexture::Albedo ? b.albedoUv : look.texture == MeshTexture::Ply ? b.plyUv : -1;
    if (uv >= 0) {
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glTexCoordPointer(2, GL_FLOAT, 0, reinterpret_cast<const void*>(uv));
    }
}

void unbindMesh() {
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    gl::BindBuffer(gl::kArrayBuffer, 0);
    gl::BindBuffer(gl::kElementArrayBuffer, 0);
}

// Flat lighting: the fixed-function headlight (a directional light fixed to the camera at
// (0.3, 0.5, 1) in eye space, Gouraud shading). Two-sided lighting (the back face's normal flipped)
// except on yarn tubes, whose outward normals hold as they are (MeshLook::outwardNormals).
void drawMeshesFlat(const std::vector<MeshLook>& looks, MeshFill fill) {
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    const float lp[4] = {0.3f, 0.5f, 1.0f, 0.0f}, amb[4] = {0.28f, 0.28f, 0.28f, 1}, dif[4] = {0.78f, 0.78f, 0.78f, 1};
    glLightfv(GL_LIGHT0, GL_POSITION, lp);
    glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
    glPopMatrix();
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glPolygonMode(GL_FRONT_AND_BACK, fill == MeshFill::Wireframe ? GL_LINE : GL_FILL);
    for (const MeshLook& look : looks) {
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, look.outwardNormals ? GL_FALSE : GL_TRUE);
        if (look.texture != MeshTexture::None) {
            glEnable(GL_TEXTURE_2D);
            glBindTexture(GL_TEXTURE_2D, look.textureId);
            glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        }
        bindMesh(*look.buffers, look);
        glDrawElements(GL_TRIANGLES, look.buffers->indexCount, GL_UNSIGNED_INT, nullptr);
        unbindMesh();
        glDisable(GL_TEXTURE_2D);
    }
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void drawMeshesSeeThrough(const std::vector<MeshLook>& looks) {
    if (looks.empty()) return;
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    drawMeshesFlat(looks, MeshFill::Solid);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
}

// --- Studio

GLint StudioLighting::Program::at(const std::string& name) {
    if (uniforms.count(name)) return uniforms.at(name);
    const GLint loc = gl::GetUniformLocation(id, name.c_str());
    uniforms[name] = loc;
    return loc;
}

StudioLighting::~StudioLighting() {
    if (!ready_) return;
    releaseTargets();
    for (Program* p : {&normalDepth_, &occlusion_, &blur_, &lit_}) gl::DeleteProgram(p->id);
    glDeleteTextures(1, &rotations_);
}

bool StudioLighting::init(std::string* err) {
    normalDepth_.id = linkProgram(kMeshVertex, kNormalDepthFragment, "normal/depth shader", err);
    occlusion_.id = normalDepth_.id ? linkProgram(kQuadVertex, kOcclusionFragment, "occlusion shader", err) : 0;
    blur_.id = occlusion_.id ? linkProgram(kQuadVertex, kBlurFragment, "blur shader", err) : 0;
    lit_.id = blur_.id ? linkProgram(kMeshVertex, kLitFragment, "lighting shader", err) : 0;
    if (!lit_.id) return false;
    for (const knit::Vec3& k : ssaoKernel(kKernelSize, 7))
        kernel_.insert(kernel_.end(), {float(k.x), float(k.y), float(k.z)});
    unsigned char rgb[16 * 3];
    const std::vector<knit::Vec3> rot = ssaoRotations(11);
    for (int i = 0; i < 16; ++i) {
        rgb[3 * i] = static_cast<unsigned char>((rot[i].x * 0.5 + 0.5) * 255 + 0.5);
        rgb[3 * i + 1] = static_cast<unsigned char>((rot[i].y * 0.5 + 0.5) * 255 + 0.5);
        rgb[3 * i + 2] = 128;
    }
    glGenTextures(1, &rotations_);
    glBindTexture(GL_TEXTURE_2D, rotations_);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 4, 4, 0, GL_RGB, GL_UNSIGNED_BYTE, rgb);
    glBindTexture(GL_TEXTURE_2D, 0);
    // The off-screen buffers must be renderable here (a float colour target needs
    // ARB_texture_float): checked once at a small size.
    if (!resizeTargets(64, 64)) {
        *err = "the off-screen buffers are not supported (framebuffer incomplete)";
        return false;
    }
    ready_ = true;
    return true;
}

void StudioLighting::releaseTargets() {
    for (Target* t : {&normalDepthTarget_, &occlusionTarget_, &blurTarget_}) {
        if (t->framebuffer) gl::DeleteFramebuffers(1, &t->framebuffer);
        if (t->texture) glDeleteTextures(1, &t->texture);
        if (t->depth) gl::DeleteRenderbuffers(1, &t->depth);
        *t = Target{};
    }
    width_ = height_ = 0;
}

bool StudioLighting::resizeTargets(int width, int height) {
    if (width == width_ && height == height_) return true;
    releaseTargets();
    bool complete = true;
    for (Target* t : {&normalDepthTarget_, &occlusionTarget_, &blurTarget_}) {
        const bool normalDepth = t == &normalDepthTarget_;
        glGenTextures(1, &t->texture);
        glBindTexture(GL_TEXTURE_2D, t->texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, normalDepth ? GL_NEAREST : GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, normalDepth ? GL_NEAREST : GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, gl::kClampToEdge);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, gl::kClampToEdge);
        // The distance needs float precision (a millimetre groove at a metre's distance).
        if (normalDepth) glTexImage2D(GL_TEXTURE_2D, 0, gl::kRgba32f, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);
        else glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        glBindTexture(GL_TEXTURE_2D, 0);
        gl::GenFramebuffers(1, &t->framebuffer);
        gl::BindFramebuffer(gl::kFramebuffer, t->framebuffer);
        gl::FramebufferTexture2D(gl::kFramebuffer, gl::kColorAttachment0, GL_TEXTURE_2D, t->texture, 0);
        if (normalDepth) {
            gl::GenRenderbuffers(1, &t->depth);
            gl::BindRenderbuffer(gl::kRenderbuffer, t->depth);
            gl::RenderbufferStorage(gl::kRenderbuffer, gl::kDepthComponent24, width, height);
            gl::BindRenderbuffer(gl::kRenderbuffer, 0);
            gl::FramebufferRenderbuffer(gl::kFramebuffer, gl::kDepthAttachment, gl::kRenderbuffer, t->depth);
        }
        complete = complete && gl::CheckFramebufferStatus(gl::kFramebuffer) == gl::kFramebufferComplete;
    }
    gl::BindFramebuffer(gl::kFramebuffer, 0);
    width_ = width;
    height_ = height;
    return complete;
}

void StudioLighting::drawNormalDepth(const std::vector<MeshLook>& looks) {
    gl::BindFramebuffer(gl::kFramebuffer, normalDepthTarget_.framebuffer);
    glViewport(0, 0, width_, height_);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    gl::UseProgram(normalDepth_.id);
    gl::Uniform1f(normalDepth_.at("farDistance"), float(farDistance_));
    for (const MeshLook& look : looks) {
        if (look.buffers->normals < 0) continue;
        gl::Uniform1i(normalDepth_.at("outwardNormals"), look.outwardNormals ? 1 : 0);
        bindMesh(*look.buffers, look);
        glDrawElements(GL_TRIANGLES, look.buffers->indexCount, GL_UNSIGNED_INT, nullptr);
        unbindMesh();
    }
}

void StudioLighting::computeOcclusion(const StudioView& view) {
    gl::BindFramebuffer(gl::kFramebuffer, occlusionTarget_.framebuffer);
    glViewport(0, 0, width_, height_);
    glDisable(GL_DEPTH_TEST);
    gl::UseProgram(occlusion_.id);
    gl::ActiveTexture(gl::kTexture0);
    glBindTexture(GL_TEXTURE_2D, normalDepthTarget_.texture);
    gl::ActiveTexture(gl::kTexture0 + 1);
    glBindTexture(GL_TEXTURE_2D, rotations_);
    gl::Uniform1i(occlusion_.at("normalDepth"), 0);
    gl::Uniform1i(occlusion_.at("rotations"), 1);
    gl::Uniform3fv(occlusion_.at("kernel"), kKernelSize, kernel_.data());
    const double aspect = double(view.width) / view.height;
    gl::Uniform2f(occlusion_.at("tanHalf"), float(view.tanHalfFov * aspect), float(view.tanHalfFov));
    gl::Uniform1f(occlusion_.at("radiusNear"), float(view.occlusionRadius));
    gl::Uniform1f(occlusion_.at("radiusFar"), float(4 * view.occlusionRadius));
    gl::Uniform1f(occlusion_.at("farDistance"), float(view.farDistance));
    drawFullScreenQuad();
    glBindTexture(GL_TEXTURE_2D, 0);
    gl::ActiveTexture(gl::kTexture0);
}

void StudioLighting::blurOcclusion(const StudioView& view) {
    gl::BindFramebuffer(gl::kFramebuffer, blurTarget_.framebuffer);
    glViewport(0, 0, width_, height_);
    gl::UseProgram(blur_.id);
    gl::ActiveTexture(gl::kTexture0);
    glBindTexture(GL_TEXTURE_2D, occlusionTarget_.texture);
    gl::ActiveTexture(gl::kTexture0 + 1);
    glBindTexture(GL_TEXTURE_2D, normalDepthTarget_.texture);
    gl::Uniform1i(blur_.at("occlusion"), 0);
    gl::Uniform1i(blur_.at("normalDepth"), 1);
    gl::Uniform2f(blur_.at("texel"), 1.0f / width_, 1.0f / height_);
    gl::Uniform1f(blur_.at("depthRange"), float(view.occlusionRadius));
    gl::Uniform1f(blur_.at("farDistance"), float(view.farDistance));
    drawFullScreenQuad();
    glBindTexture(GL_TEXTURE_2D, 0);
    gl::ActiveTexture(gl::kTexture0);
}

void StudioLighting::drawLit(const StudioView& view, const std::vector<MeshLook>& looks) {
    gl::BindFramebuffer(gl::kFramebuffer, 0);
    glViewport(view.x, view.y, view.width, view.height);
    glEnable(GL_SCISSOR_TEST);  // only the columns of this pass (one half while comparing)
    glScissor(view.clipX, view.y, view.clipWidth, view.height);
    glEnable(GL_DEPTH_TEST);
    gl::UseProgram(lit_.id);
    gl::ActiveTexture(gl::kTexture0 + 1);
    glBindTexture(GL_TEXTURE_2D, blurTarget_.texture);
    gl::ActiveTexture(gl::kTexture0);
    gl::Uniform1i(lit_.at("surfaceTexture"), 0);
    gl::Uniform1i(lit_.at("occlusion"), 1);
    gl::Uniform2f(lit_.at("viewportOrigin"), float(view.x), float(view.y));
    gl::Uniform2f(lit_.at("viewportSize"), float(view.width), float(view.height));
    gl::Uniform1f(lit_.at("occlusionStrength"), float(view.occlusionStrength));
    gl::Uniform3f(lit_.at("sky"), float(view.skyInView.x), float(view.skyInView.y), float(view.skyInView.z));
    gl::Uniform3f(lit_.at("keyDirection"), kKeyDirection[0], kKeyDirection[1], kKeyDirection[2]);
    gl::Uniform3f(lit_.at("fillDirection"), kFillDirection[0], kFillDirection[1], kFillDirection[2]);
    for (const MeshLook& look : looks) {
        const MeshBuffers* b = look.buffers;
        if (b->normals < 0) continue;
        const bool textured = look.textureId && (look.texture == MeshTexture::Albedo ? b->albedoUv >= 0 : b->plyUv >= 0);
        gl::Uniform1i(lit_.at("useTexture"), textured ? 1 : 0);
        gl::Uniform1i(lit_.at("outwardNormals"), look.outwardNormals ? 1 : 0);
        glBindTexture(GL_TEXTURE_2D, textured ? look.textureId : 0);
        bindMesh(*b, look);
        glDrawElements(GL_TRIANGLES, b->indexCount, GL_UNSIGNED_INT, nullptr);
        unbindMesh();
    }
    glBindTexture(GL_TEXTURE_2D, 0);
    gl::ActiveTexture(gl::kTexture0 + 1);
    glBindTexture(GL_TEXTURE_2D, 0);
    gl::ActiveTexture(gl::kTexture0);
}

bool StudioLighting::draw(const StudioView& view, const std::vector<MeshLook>& looks) {
    if (!ready_) return false;
    if (!resizeTargets(view.width, view.height)) {
        gl::BindFramebuffer(gl::kFramebuffer, 0);
        return false;
    }
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    farDistance_ = view.farDistance;
    // The off-screen passes cover the whole viewport, so the occlusion next to a comparison's divider
    // sees its real neighbours; drawLit clips.
    glDisable(GL_SCISSOR_TEST);
    drawNormalDepth(looks);
    computeOcclusion(view);
    blurOcclusion(view);
    drawLit(view, looks);
    gl::UseProgram(0);
    return true;
}

// --- the sample pattern

namespace {

// A uniform number in [0, 1) from the generator's 24 high bits (the same on every platform,
// unlike std::uniform_real_distribution).
double uniform01(std::mt19937& rng) { return double(rng() >> 8) / 16777216.0; }

} // namespace

std::vector<knit::Vec3> ssaoKernel(int n, uint32_t seed) {
    std::mt19937 rng(seed);
    std::vector<knit::Vec3> kernel;
    while (static_cast<int>(kernel.size()) < n) {
        const knit::Vec3 p{2 * uniform01(rng) - 1, 2 * uniform01(rng) - 1, uniform01(rng)};
        const double len = knit::norm(p);
        if (len > 1 || len < 1e-3 || p.z < 0.15 * len) continue;  // uniform directions in the cone z >= 0.15 |p|
        const double t = double(kernel.size()) / n;
        const double scale = (0.1 + 0.9 * t * t) * (0.25 + 0.75 * uniform01(rng));
        kernel.push_back(p * (scale / len));
    }
    return kernel;
}

std::vector<knit::Vec3> ssaoRotations(uint32_t seed) {
    std::mt19937 rng(seed);
    std::vector<knit::Vec3> rotations;
    for (int i = 0; i < 16; ++i) {
        const double a = 2 * knit::kPi * uniform01(rng);
        rotations.push_back({std::cos(a), std::sin(a), 0});
    }
    return rotations;
}

} // namespace viewer
