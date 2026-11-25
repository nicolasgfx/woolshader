// The OpenGL entry points above version 1.1 that the viewer's Studio lighting uses: shaders (2.0),
// buffer objects (1.5), framebuffer objects (3.0 / ARB_framebuffer_object) and texture units
// (1.3). Windows' opengl32 exports only 1.1, so they are loaded at run time (glfwGetProcAddress)
// once the context exists (Renderer::init); a driver without one of them keeps the viewer on Flat
// lighting.
#pragma once

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif
#include <GL/gl.h>

#include <cstddef>

#ifndef APIENTRY
#define APIENTRY
#endif

namespace viewer::gl {

// Enum values from the OpenGL registry (glext.h names in the comments).
constexpr GLenum kFragmentShader = 0x8B30;        // GL_FRAGMENT_SHADER
constexpr GLenum kVertexShader = 0x8B31;          // GL_VERTEX_SHADER
constexpr GLenum kCompileStatus = 0x8B81;         // GL_COMPILE_STATUS
constexpr GLenum kLinkStatus = 0x8B82;            // GL_LINK_STATUS
constexpr GLenum kArrayBuffer = 0x8892;           // GL_ARRAY_BUFFER
constexpr GLenum kElementArrayBuffer = 0x8893;    // GL_ELEMENT_ARRAY_BUFFER
constexpr GLenum kStaticDraw = 0x88E4;            // GL_STATIC_DRAW
constexpr GLenum kFramebuffer = 0x8D40;           // GL_FRAMEBUFFER
constexpr GLenum kRenderbuffer = 0x8D41;          // GL_RENDERBUFFER
constexpr GLenum kColorAttachment0 = 0x8CE0;      // GL_COLOR_ATTACHMENT0
constexpr GLenum kDepthAttachment = 0x8D00;       // GL_DEPTH_ATTACHMENT
constexpr GLenum kFramebufferComplete = 0x8CD5;   // GL_FRAMEBUFFER_COMPLETE
constexpr GLenum kDepthComponent24 = 0x81A6;      // GL_DEPTH_COMPONENT24
constexpr GLenum kRgba32f = 0x8814;               // GL_RGBA32F
constexpr GLenum kTexture0 = 0x84C0;              // GL_TEXTURE0
constexpr GLenum kClampToEdge = 0x812F;           // GL_CLAMP_TO_EDGE

inline GLuint(APIENTRY* CreateShader)(GLenum type) = nullptr;
inline void(APIENTRY* ShaderSource)(GLuint shader, GLsizei count, const char* const* source, const GLint* length) = nullptr;
inline void(APIENTRY* CompileShader)(GLuint shader) = nullptr;
inline void(APIENTRY* GetShaderiv)(GLuint shader, GLenum name, GLint* value) = nullptr;
inline void(APIENTRY* GetShaderInfoLog)(GLuint shader, GLsizei size, GLsizei* length, char* log) = nullptr;
inline void(APIENTRY* DeleteShader)(GLuint shader) = nullptr;
inline GLuint(APIENTRY* CreateProgram)() = nullptr;
inline void(APIENTRY* AttachShader)(GLuint program, GLuint shader) = nullptr;
inline void(APIENTRY* LinkProgram)(GLuint program) = nullptr;
inline void(APIENTRY* GetProgramiv)(GLuint program, GLenum name, GLint* value) = nullptr;
inline void(APIENTRY* GetProgramInfoLog)(GLuint program, GLsizei size, GLsizei* length, char* log) = nullptr;
inline void(APIENTRY* UseProgram)(GLuint program) = nullptr;
inline void(APIENTRY* DeleteProgram)(GLuint program) = nullptr;
inline GLint(APIENTRY* GetUniformLocation)(GLuint program, const char* name) = nullptr;
inline void(APIENTRY* Uniform1i)(GLint location, GLint v) = nullptr;
inline void(APIENTRY* Uniform1f)(GLint location, GLfloat v) = nullptr;
inline void(APIENTRY* Uniform2f)(GLint location, GLfloat x, GLfloat y) = nullptr;
inline void(APIENTRY* Uniform3f)(GLint location, GLfloat x, GLfloat y, GLfloat z) = nullptr;
inline void(APIENTRY* Uniform3fv)(GLint location, GLsizei count, const GLfloat* v) = nullptr;
inline void(APIENTRY* GenBuffers)(GLsizei n, GLuint* ids) = nullptr;
inline void(APIENTRY* BindBuffer)(GLenum target, GLuint id) = nullptr;
inline void(APIENTRY* BufferData)(GLenum target, std::ptrdiff_t size, const void* data, GLenum usage) = nullptr;
inline void(APIENTRY* BufferSubData)(GLenum target, std::ptrdiff_t offset, std::ptrdiff_t size, const void* data) = nullptr;
inline void(APIENTRY* DeleteBuffers)(GLsizei n, const GLuint* ids) = nullptr;
inline void(APIENTRY* GenFramebuffers)(GLsizei n, GLuint* ids) = nullptr;
inline void(APIENTRY* BindFramebuffer)(GLenum target, GLuint id) = nullptr;
inline void(APIENTRY* FramebufferTexture2D)(GLenum target, GLenum attachment, GLenum texTarget, GLuint texture, GLint level) = nullptr;
inline GLenum(APIENTRY* CheckFramebufferStatus)(GLenum target) = nullptr;
inline void(APIENTRY* DeleteFramebuffers)(GLsizei n, const GLuint* ids) = nullptr;
inline void(APIENTRY* GenRenderbuffers)(GLsizei n, GLuint* ids) = nullptr;
inline void(APIENTRY* BindRenderbuffer)(GLenum target, GLuint id) = nullptr;
inline void(APIENTRY* RenderbufferStorage)(GLenum target, GLenum format, GLsizei width, GLsizei height) = nullptr;
inline void(APIENTRY* FramebufferRenderbuffer)(GLenum target, GLenum attachment, GLenum rbTarget, GLuint renderbuffer) = nullptr;
inline void(APIENTRY* DeleteRenderbuffers)(GLsizei n, const GLuint* ids) = nullptr;
inline void(APIENTRY* ActiveTexture)(GLenum unit) = nullptr;

} // namespace viewer::gl
