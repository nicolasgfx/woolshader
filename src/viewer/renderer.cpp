// One frame of the 3D view, and what it keeps on the GPU (renderer.h).
#include "viewer/renderer.h"
#include "viewer/overlay_lines.h"

#include <GLFW/glfw3.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "io/image_io.h"

namespace viewer {

namespace {

// While the camera moves the fibres are thinned to this many segments.
constexpr long long kMovingSegments = 300000;
constexpr float kPlainNodeGrey = 0.38f;  // an open file's unticked component
constexpr float kFocusOutline[3] = {0.26f, 0.52f, 0.90f};  // an open file's focused component: the interface's accent

// Keeps the first missing name and goes on: buffers still load on a driver without framebuffers.
GLFWglproc find(const char* name, std::string* missing) {
    const GLFWglproc p = glfwGetProcAddress(name);
    if (!p && missing->empty()) *missing = name;
    return p;
}

// Every entry point of gl_functions.h from the current context; false names the first one missing.
bool loadGlFunctions(std::string* missing) {
    using namespace gl;
    CreateShader = reinterpret_cast<decltype(CreateShader)>(find("glCreateShader", missing));
    ShaderSource = reinterpret_cast<decltype(ShaderSource)>(find("glShaderSource", missing));
    CompileShader = reinterpret_cast<decltype(CompileShader)>(find("glCompileShader", missing));
    GetShaderiv = reinterpret_cast<decltype(GetShaderiv)>(find("glGetShaderiv", missing));
    GetShaderInfoLog = reinterpret_cast<decltype(GetShaderInfoLog)>(find("glGetShaderInfoLog", missing));
    DeleteShader = reinterpret_cast<decltype(DeleteShader)>(find("glDeleteShader", missing));
    CreateProgram = reinterpret_cast<decltype(CreateProgram)>(find("glCreateProgram", missing));
    AttachShader = reinterpret_cast<decltype(AttachShader)>(find("glAttachShader", missing));
    LinkProgram = reinterpret_cast<decltype(LinkProgram)>(find("glLinkProgram", missing));
    GetProgramiv = reinterpret_cast<decltype(GetProgramiv)>(find("glGetProgramiv", missing));
    GetProgramInfoLog = reinterpret_cast<decltype(GetProgramInfoLog)>(find("glGetProgramInfoLog", missing));
    UseProgram = reinterpret_cast<decltype(UseProgram)>(find("glUseProgram", missing));
    DeleteProgram = reinterpret_cast<decltype(DeleteProgram)>(find("glDeleteProgram", missing));
    GetUniformLocation = reinterpret_cast<decltype(GetUniformLocation)>(find("glGetUniformLocation", missing));
    Uniform1i = reinterpret_cast<decltype(Uniform1i)>(find("glUniform1i", missing));
    Uniform1f = reinterpret_cast<decltype(Uniform1f)>(find("glUniform1f", missing));
    Uniform2f = reinterpret_cast<decltype(Uniform2f)>(find("glUniform2f", missing));
    Uniform3f = reinterpret_cast<decltype(Uniform3f)>(find("glUniform3f", missing));
    Uniform3fv = reinterpret_cast<decltype(Uniform3fv)>(find("glUniform3fv", missing));
    GenBuffers = reinterpret_cast<decltype(GenBuffers)>(find("glGenBuffers", missing));
    BindBuffer = reinterpret_cast<decltype(BindBuffer)>(find("glBindBuffer", missing));
    BufferData = reinterpret_cast<decltype(BufferData)>(find("glBufferData", missing));
    BufferSubData = reinterpret_cast<decltype(BufferSubData)>(find("glBufferSubData", missing));
    DeleteBuffers = reinterpret_cast<decltype(DeleteBuffers)>(find("glDeleteBuffers", missing));
    GenFramebuffers = reinterpret_cast<decltype(GenFramebuffers)>(find("glGenFramebuffers", missing));
    BindFramebuffer = reinterpret_cast<decltype(BindFramebuffer)>(find("glBindFramebuffer", missing));
    FramebufferTexture2D = reinterpret_cast<decltype(FramebufferTexture2D)>(find("glFramebufferTexture2D", missing));
    CheckFramebufferStatus = reinterpret_cast<decltype(CheckFramebufferStatus)>(find("glCheckFramebufferStatus", missing));
    DeleteFramebuffers = reinterpret_cast<decltype(DeleteFramebuffers)>(find("glDeleteFramebuffers", missing));
    GenRenderbuffers = reinterpret_cast<decltype(GenRenderbuffers)>(find("glGenRenderbuffers", missing));
    BindRenderbuffer = reinterpret_cast<decltype(BindRenderbuffer)>(find("glBindRenderbuffer", missing));
    RenderbufferStorage = reinterpret_cast<decltype(RenderbufferStorage)>(find("glRenderbufferStorage", missing));
    FramebufferRenderbuffer = reinterpret_cast<decltype(FramebufferRenderbuffer)>(find("glFramebufferRenderbuffer", missing));
    DeleteRenderbuffers = reinterpret_cast<decltype(DeleteRenderbuffers)>(find("glDeleteRenderbuffers", missing));
    ActiveTexture = reinterpret_cast<decltype(ActiveTexture)>(find("glActiveTexture", missing));
    return missing->empty();
}

// The buffer objects of a mesh: one vertex buffer holding its arrays one after the other.
MeshBuffers uploadMesh(const knit::DisplayMesh& d) {
    MeshBuffers b;
    const size_t n = d.pos.size() / 3;
    long long size = 0;
    auto place = [&](const std::vector<float>& v, size_t perVertex) {
        if (v.size() != n * perVertex || v.empty()) return -1LL;
        const long long at = size;
        size += static_cast<long long>(v.size() * sizeof(float));
        return at;
    };
    place(d.pos, 3);
    b.normals = place(d.nrm, 3);
    b.colours = place(d.vcol, 3);
    b.plyUv = place(d.uv, 2);
    b.albedoUv = place(d.albedoUv, 2);
    gl::GenBuffers(1, &b.vertices);
    gl::BindBuffer(gl::kArrayBuffer, b.vertices);
    gl::BufferData(gl::kArrayBuffer, size, nullptr, gl::kStaticDraw);
    gl::BufferSubData(gl::kArrayBuffer, 0, static_cast<long long>(d.pos.size() * sizeof(float)), d.pos.data());
    if (b.normals >= 0) gl::BufferSubData(gl::kArrayBuffer, b.normals, static_cast<long long>(d.nrm.size() * sizeof(float)), d.nrm.data());
    if (b.colours >= 0) gl::BufferSubData(gl::kArrayBuffer, b.colours, static_cast<long long>(d.vcol.size() * sizeof(float)), d.vcol.data());
    if (b.plyUv >= 0) gl::BufferSubData(gl::kArrayBuffer, b.plyUv, static_cast<long long>(d.uv.size() * sizeof(float)), d.uv.data());
    if (b.albedoUv >= 0)
        gl::BufferSubData(gl::kArrayBuffer, b.albedoUv, static_cast<long long>(d.albedoUv.size() * sizeof(float)), d.albedoUv.data());
    gl::BindBuffer(gl::kArrayBuffer, 0);
    gl::GenBuffers(1, &b.indices);
    gl::BindBuffer(gl::kElementArrayBuffer, b.indices);
    gl::BufferData(gl::kElementArrayBuffer, static_cast<long long>(d.idx.size() * sizeof(uint32_t)), d.idx.data(), gl::kStaticDraw);
    gl::BindBuffer(gl::kElementArrayBuffer, 0);
    b.indexCount = static_cast<GLsizei>(d.idx.size());
    b.bytes = size + static_cast<long long>(d.idx.size() * sizeof(uint32_t));
    return b;
}

// A repeating, mipmapped texture from tightly packed rows.
GLuint uploadTexture(int width, int height, GLenum format, const void* pixels) {
    GLuint id = 0;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, 0x8191 /* GL_GENERATE_MIPMAP (GL 1.4) */, GL_TRUE);
    glTexImage2D(GL_TEXTURE_2D, 0, static_cast<GLint>(format), width, height, 0, format, GL_UNSIGNED_BYTE, pixels);
    glBindTexture(GL_TEXTURE_2D, 0);
    return id;
}

bool holdsMesh(const std::vector<knit::DisplayMesh>& meshes, const knit::DisplayMesh& d) {
    return !meshes.empty() && &d >= meshes.data() && &d < meshes.data() + meshes.size();
}

bool listsScene(const MeshOwners& o, const knit::BuiltScene* s) {
    for (const std::shared_ptr<const knit::BuiltScene>& p : o.scenes)
        if (p.get() == s) return true;
    return false;
}

bool listsTexture(const MeshOwners& o, const knit::DisplayTexture* t) {
    for (const std::shared_ptr<const knit::BuiltScene>& s : o.scenes)
        for (const std::shared_ptr<const knit::DisplayTexture>& p : s->albedo)
            if (p.get() == t) return true;
    return false;
}

// Keeps the scene that holds the mesh alive (the mesh address stays this mesh's); false: none does.
bool keepOwnerOf(MeshOwners& kept, const knit::DisplayMesh& d, const MeshOwners& owners) {
    for (const std::shared_ptr<const knit::BuiltScene>& s : owners.scenes)
        if (holdsMesh(s->meshes, d)) {
            if (!listsScene(kept, s.get())) kept.scenes.push_back(s);
            return true;
        }
    if (!owners.preview || !holdsMesh(owners.preview->meshes, d)) return false;
    kept.preview = owners.preview;
    return true;
}

// The camera's matrices, and the import preview's file -> library transform: scale, the up axis
// (Y-up: (x, y, z) -> (x, -z, y)), centring.
void loadCamera(const Mat4& P, const Mat4& V, const ViewContent& c) {
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixd(P.m);
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixd(V.m);
    if (!c.preview) return;
    glTranslated(c.preview->offset.x, c.preview->offset.y, c.preview->offset.z);
    glScaled(c.preview->scale, c.preview->scale, c.preview->scale);
    if (c.preview->yUp) {
        const double m[16] = {1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1};
        glMultMatrixd(m);
    }
}

// The ambient occlusion's near radius [scene mm]: the shown yarn's diameter (the grooves between
// neighbouring yarns), else 1 % of the content's size (imported meshes, hair bases).
double occlusionRadius(const ViewContent& c, const knit::SceneSettings& s) {
    bool yarn = false;
    for (const knit::DisplayMesh* d : c.meshes) yarn |= d->kind == knit::DisplayMesh::Kind::Yarn;
    if (!yarn) return std::max(0.01 * knit::norm(c.bounds.hi - c.bounds.lo), 0.1);
    double D = 0;
    if (c.meshOwner)
        for (const knit::PartRow& p : c.meshOwner->parts) D = std::max(D, p.knitted ? p.yarnD : 0.0);
    return D > 0 ? D : s.diameter;
}

} // namespace

ViewportPixels viewportPixels(const SceneUi& ui, int framebufferWidth, int framebufferHeight, float displayWidth, float displayHeight) {
    const ViewportRect r = viewportRect(ui, displayWidth, displayHeight);
    const float sx = displayWidth > 0 ? framebufferWidth / displayWidth : 1, sy = displayHeight > 0 ? framebufferHeight / displayHeight : 1;
    ViewportPixels v;
    v.x = static_cast<int>(std::lround(r.x * sx));
    v.width = std::max(1, static_cast<int>(std::lround(r.w * sx)));
    v.height = std::max(1, static_cast<int>(std::lround(r.h * sy)));
    v.y = std::max(0, static_cast<int>(std::lround((displayHeight - r.y - r.h) * sy)));
    return v;
}

Renderer::~Renderer() {
    release(MeshOwners{});
    if (plyTexture_) glDeleteTextures(1, &plyTexture_);
}

bool Renderer::init(std::string* why) {
    std::string missing;
    if (!loadGlFunctions(&missing)) {
        *why = "OpenGL entry point " + missing + " is missing";
        return false;
    }
    std::string discWhy;
    furDiscs_.init(&discWhy);  // without it the discs show as their marks; discNote() says why
    return studio_.init(why);
}

// The mesh's buffers, uploaded on first use; null when no listed owner holds the mesh.
const MeshBuffers* Renderer::buffersOf(const knit::DisplayMesh& d, const MeshOwners& owners) {
    if (buffers_.count(&d)) return &buffers_.at(&d);
    if (!keepOwnerOf(kept_, d, owners)) return nullptr;
    const double start = glfwGetTime();
    MeshBuffers& b = buffers_[&d] = uploadMesh(d);
    uploadMs_ += 1000 * (glfwGetTime() - start);
    return &b;
}

GLuint Renderer::albedoOf(const ViewContent& c, const knit::DisplayMesh& d) {
    if (d.albedo < 0 || !c.meshOwner || d.albedo >= static_cast<int>(c.meshOwner->albedo.size())) return 0;
    const knit::DisplayTexture* t = c.meshOwner->albedo[d.albedo].get();
    if (albedo_.count(t)) return albedo_.at(t);
    // Rows top first; the meshes' albedo uv are in that image space (scene_builder.h).
    const double start = glfwGetTime();
    const GLuint id = uploadTexture(t->width, t->height, GL_RGBA, t->rgba.data());
    uploadMs_ += 1000 * (glfwGetTime() - start);
    albedoBytes_[id] = static_cast<long long>(t->width) * t->height * 4 * 4 / 3;  // the mipmaps add a third
    return albedo_[t] = id;
}

// The ply-map preview as a repeating texture (rows flipped so that t = v, ply_profile.h).
void Renderer::keepPlyTexture(const std::shared_ptr<const knit::TexturePreview>& preview) {
    if (preview == plyPreview_) return;
    plyPreview_ = preview;
    if (plyTexture_) glDeleteTextures(1, &plyTexture_);
    plyTexture_ = 0;
    if (!preview || preview->rgb.empty()) return;
    std::vector<uint8_t> flipped(preview->rgb.size());
    const size_t row = size_t(preview->width) * 3;
    for (int y = 0; y < preview->height; ++y) std::memcpy(&flipped[size_t(y) * row], &preview->rgb[size_t(preview->height - 1 - y) * row], row);
    plyTexture_ = uploadTexture(preview->width, preview->height, GL_RGB, flipped.data());
}

// How each shown mesh is drawn. Yarn: the CURRENT settings' colour (material edits show at once),
// the part's object colour or the slot colour; per-vertex stitch colours only while the texture
// (colourwork) source is chosen. A base-colour texture, else the opted-in ply texture on yarn.
std::vector<MeshLook> Renderer::meshLooks(const SceneUi& ui, const ViewContent& c, const MeshOwners& owners) {
    const knit::SceneSettings& s = ui.s;
    std::vector<MeshLook> looks;
    for (size_t i = 0; i < c.meshes.size(); ++i) {
        const knit::DisplayMesh& d = *c.meshes[i];
        if (d.idx.empty()) continue;
        MeshLook look;
        look.alpha = i < c.meshAlpha.size() ? c.meshAlpha[i] : 1.0f;
        look.buffers = buffersOf(d, owners);
        if (!look.buffers) continue;
        if (previewNodeIsPlain(ui, d)) {
            for (int k = 0; k < 3; ++k) look.color[k] = kPlainNodeGrey;
            looks.push_back(look);
            continue;
        }
        const bool yarn = d.kind == knit::DisplayMesh::Kind::Yarn;
        if (yarn && d.hasPartColor && s.yarnColourSource != "single") {
            for (int k = 0; k < 3; ++k) look.color[k] = toSrgb(d.partColor[k]);
        } else if (yarn) {
            const Vec3& col = s.yarnWhite.color;
            look.color[0] = float(col.x), look.color[1] = float(col.y), look.color[2] = float(col.z);
        } else {
            for (int k = 0; k < 3; ++k) look.color[k] = d.color[k];
        }
        look.outwardNormals = yarn;
        look.vertexColour = !d.vcol.empty() && d.vcol.size() == d.pos.size() && (!yarn || s.yarnColourSource == "texture");
        const GLuint albedo = albedoOf(c, d);
        if (albedo && d.albedoUv.size() * 3 == d.pos.size() * 2 && (!yarn || ui.view.yarnMaps)) {
            look.texture = MeshTexture::Albedo;
            look.textureId = albedo;
        } else if (yarn && ui.view.plyTexture && plyTexture_ && d.uv.size() * 3 == d.pos.size() * 2) {
            look.texture = MeshTexture::Ply;
            look.textureId = plyTexture_;
        }
        looks.push_back(look);
    }
    return looks;
}

// The Import step's focused component, drawn again as accent lines over the shaded meshes (the depth
// buffer of the pass before hides what is behind; the offset keeps the lines in front of their own faces).
void Renderer::drawFocusOutline(const SceneUi& ui, const ViewContent& c, const MeshOwners& owners) {
    std::vector<MeshLook> focused;
    for (const knit::DisplayMesh* d : c.meshes) {
        if (d->idx.empty() || !previewNodeIsFocused(ui, *d)) continue;
        MeshLook look;
        look.buffers = buffersOf(*d, owners);
        if (!look.buffers) continue;
        for (int k = 0; k < 3; ++k) look.color[k] = kFocusOutline[k];
        focused.push_back(look);
    }
    if (focused.empty()) return;
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_POLYGON_OFFSET_LINE);
    glPolygonOffset(-1.0f, -1.0f);
    drawMeshesFlat(focused, MeshFill::Wireframe);
    glDisable(GL_POLYGON_OFFSET_LINE);
    glDepthFunc(GL_LESS);
}

// One pass: this content with the frame's camera, written only into its columns.
void Renderer::drawContent(const SceneUi& ui, const ViewContent& c, const Camera& camera, const FrameTarget& target, ViewColumns columns) {
    const ViewportPixels& v = target.viewport;
    glViewport(v.x, v.y, v.width, v.height);
    glEnable(GL_SCISSOR_TEST);
    glScissor(columns.x, v.y, columns.width, v.height);
    const Mat4 P = camera.proj(double(v.width) / v.height), V = camera.view();
    loadCamera(P, V, c);  // the source keeps its own preview transform
    glEnable(GL_DEPTH_TEST);
    const MeshOwners owners = meshOwners(ui);
    const std::vector<MeshLook> all = meshLooks(ui, c, owners);
    std::vector<MeshLook> looks, seeThrough;  // the opaque meshes first, then the see-through ones over them
    for (const MeshLook& l : all) {
        if (l.alpha < 1) seeThrough.push_back(l);
        else looks.push_back(l);
    }
    bool drawn = false;
    if (ui.view.lighting == Lighting::Studio && !ui.view.wire && studio_.ready()) {
        StudioView sv;
        sv.x = v.x, sv.y = v.y, sv.width = v.width, sv.height = v.height;
        sv.clipX = columns.x, sv.clipWidth = columns.width;
        sv.tanHalfFov = std::tan(camera.fovDeg * knit::kPi / 360.0);
        sv.farDistance = camera.farPlane();
        sv.skyInView = {V.m[8], V.m[9], V.m[10]};  // the view matrix (column-major) applied to world +z
        sv.occlusionRadius = occlusionRadius(c, ui.s);
        sv.occlusionStrength = ui.view.occlusion;
        drawn = studio_.draw(sv, looks);
    }
    if (!drawn) drawMeshesFlat(looks, ui.view.wire ? MeshFill::Wireframe : MeshFill::Solid);
    drawMeshesSeeThrough(seeThrough);
    drawFocusOutline(ui, c, owners);
    if (c.overlay) {
        loadCamera(P, V, c);  // the overlay is in the preview's file coordinates, as the meshes
        drawOverlay(*c.overlay);
    }
    if (c.furDiscs) {
        // The disc preview (fur_discs.h), a quarter of the fibres while the camera moves.
        loadCamera(P, V, c);
        uploadMs_ += furDiscs_.draw(c.furDiscs, target.cameraMoving ? 0.25f : 1.0f);
    }
    if (c.fibres && !c.fibres->pos.empty()) {
        // Unlit lines (fibre_lines.h), coloured from the CURRENT materials.
        FibreCamera fibreCam = fibreCamera(P, V);
        fibreCam.cull = !c.preview;
        const FibreColours exported = fibreColoursOf(ui.s);
        const FibreColours colours = ui.view.fibrePopulations ? populationColours(exported) : exported;
        fibreSegments_ += fibres_.draw(ui.results.fibres, *c.fibres, colours, fibreCam, target.cameraMoving ? kMovingSegments : 0);
    }
}

void Renderer::draw(const SceneUi& ui, const ViewportFrame& frame, const Camera& camera, const FrameTarget& target) {
    const knit::BuiltScene* shown = ui.results.partial ? ui.results.partial.get() : ui.results.knit.get();
    if (shown) keepPlyTexture(shown->texture);
    glDisable(GL_SCISSOR_TEST);
    glViewport(0, 0, target.width, target.height);
    glClearColor(0.11f, 0.12f, 0.14f, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // once: the halves never overlap
    fibreSegments_ = 0;
    const ViewportPixels& v = target.viewport;
    if (frame.source) {
        const int half = v.width / 2;
        drawContent(ui, *frame.source, camera, target, ViewColumns{v.x, half});
        drawContent(ui, frame.shown, camera, target, ViewColumns{v.x + half, v.width - half});
    } else {
        drawContent(ui, frame.shown, camera, target, ViewColumns{v.x, v.width});
    }
    glDisable(GL_SCISSOR_TEST);
    glViewport(0, 0, target.width, target.height);
}

long long Renderer::residentBytes() const {
    long long bytes = 0;
    for (const std::pair<const GLuint, long long>& texture : albedoBytes_) bytes += texture.second;
    for (const std::pair<const knit::DisplayMesh* const, MeshBuffers>& entry : buffers_) bytes += entry.second.bytes;
    bytes += furDiscs_.residentBytes();  // the disc preview (fur_discs.h)
    return bytes;
}

void Renderer::uploadRenderPreview(SceneUi& ui) {
    RenderPreviewView& v = ui.renderPreview;
    const std::shared_ptr<const knit::PreviewFrame>& f = v.state.frame;
    if (!f || f->rgba.empty() || v.state.frameNumber == v.uploadedFrame) return;
    if (!previewTexture_) glGenTextures(1, &previewTexture_);
    glBindTexture(GL_TEXTURE_2D, previewTexture_);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, 0x812F /* GL_CLAMP_TO_EDGE */);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, 0x812F);
    if (f->width == v.textureW && f->height == v.textureH)
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, f->width, f->height, GL_RGBA, GL_UNSIGNED_BYTE, f->rgba.data());
    else
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, f->width, f->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, f->rgba.data());
    glBindTexture(GL_TEXTURE_2D, 0);
    v.texture = previewTexture_;
    v.textureW = f->width;
    v.textureH = f->height;
    v.uploadedFrame = v.state.frameNumber;
}

void Renderer::releaseUnused(const SceneUi& ui) {
    release(meshOwners(ui));
    fibres_.releaseUnless(ui.results.fibres.get());
}

void Renderer::release(const MeshOwners& owners) {
    auto releaseMeshes = [&](const std::vector<knit::DisplayMesh>& meshes) {
        for (const knit::DisplayMesh& d : meshes) {
            if (!buffers_.count(&d)) continue;
            const MeshBuffers& b = buffers_.at(&d);
            gl::DeleteBuffers(1, &b.vertices);
            gl::DeleteBuffers(1, &b.indices);
            buffers_.erase(&d);
        }
    };
    bool sceneLeft = false;
    for (size_t i = 0; i < kept_.scenes.size();) {
        if (listsScene(owners, kept_.scenes[i].get())) {
            ++i;
            continue;
        }
        releaseMeshes(kept_.scenes[i]->meshes);
        kept_.scenes.erase(kept_.scenes.begin() + static_cast<long long>(i));
        sceneLeft = true;
    }
    if (kept_.preview && kept_.preview != owners.preview) {
        releaseMeshes(kept_.preview->meshes);
        kept_.preview.reset();
    }
    if (!sceneLeft) return;
    // A texture goes when no kept scene lists it: the coverage overlay shares the hair's.
    std::erase_if(albedo_, [&](const std::pair<const knit::DisplayTexture* const, GLuint>& entry) {
        if (listsTexture(kept_, entry.first)) return false;
        glDeleteTextures(1, &entry.second);
        albedoBytes_.erase(entry.second);
        return true;
    });
}

bool writeBackBufferPng(const std::string& path, int width, int height) {
    std::vector<unsigned char> pixels(size_t(width) * height * 3), flipped(pixels.size());
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadBuffer(GL_BACK);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    for (int y = 0; y < height; ++y) std::memcpy(&flipped[size_t(y) * width * 3], &pixels[size_t(height - 1 - y) * width * 3], size_t(width) * 3);
    return knit::writePngRgb8(path, width, height, flipped);
}

} // namespace viewer
