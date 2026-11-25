// The Mark step's overlay (overlay_lines.h).
#include "viewer/overlay_lines.h"

#include "viewer/gl_functions.h"

namespace viewer {

namespace {

void drawArrays(GLenum mode, const std::vector<float>& pos, const std::vector<uint8_t>& rgba) {
    if (pos.empty()) return;
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, pos.data());
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, rgba.data());
    glDrawArrays(mode, 0, static_cast<GLsizei>(pos.size() / 3));
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
}

} // namespace

void drawOverlay(const OverlayGeometry& g) {
    if (g.linePos.empty() && g.pointPos.empty()) return;
    glDisable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    glLineWidth(2.5f);
    glPointSize(8.0f);
    // behind the surface: faint
    glDepthFunc(GL_GREATER);
    drawArrays(GL_LINES, g.linePos, g.lineRgbaHidden);
    drawArrays(GL_POINTS, g.pointPos, g.pointRgbaHidden);
    // in front: bright, pulled a little towards the camera so a path on the surface wins its depth test
    glDepthFunc(GL_LEQUAL);
    glDepthRange(0.0, 0.9995);
    drawArrays(GL_LINES, g.linePos, g.lineRgba);
    drawArrays(GL_POINTS, g.pointPos, g.pointRgba);
    glDepthRange(0.0, 1.0);
    glDepthFunc(GL_LESS);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glLineWidth(1.0f);
    glPointSize(1.0f);
}

} // namespace viewer
