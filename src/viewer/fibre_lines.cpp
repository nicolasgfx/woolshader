// The fibre lines in GL buffers, culled and thinned per chunk (fibre_lines.h).
#include "viewer/fibre_lines.h"

#include <algorithm>
#include <cmath>

#include "jobs/previews.h"
#include "model/colour.h"

namespace viewer {

namespace {

// Vertex v's colour from a per-vertex linear array, as the meshes are shown.
Vec3 shownColour(const std::vector<float>& linear, size_t v) { return {toSrgb(linear[3 * v]), toSrgb(linear[3 * v + 1]), toSrgb(linear[3 * v + 2])}; }

// A source as the fibres take it: the texture their own colour, "object" their part's.
FibreColourFrom colourFromSource(const std::string& source) {
    return source == "texture" ? FibreColourFrom::Own : source == "object" ? FibreColourFrom::Part : FibreColourFrom::Fixed;
}

uint8_t toByte(double c) { return static_cast<uint8_t>(std::lround(std::clamp(c, 0.0, 1.0) * 255)); }

} // namespace

FibreColours fibreColoursOf(const knit::SceneSettings& s) {
    FibreColours colours;
    colours.fur = colourFromSource(s.pileColourSource);
    colours.furColour = s.pileMaterial.color;
    const FibreColourFrom yarn = colourFromSource(s.yarnColourSource);
    colours.oriented = s.fibreInheritsYarn ? yarn : FibreColourFrom::Fixed;
    colours.stray = s.fibreRandomInheritsYarn ? yarn : FibreColourFrom::Fixed;
    return colours;
}

FibreColours populationColours(FibreColours colours) {
    colours.oriented = colours.stray = FibreColourFrom::Fixed;
    return colours;
}

Vec3 fibreVertexColour(const knit::DisplayLines& lines, size_t v, const FibreColours& colours) {
    const uint8_t population = lines.population[v];
    const FibreColourFrom from = population == 2 ? colours.fur : population == 1 ? colours.stray : colours.oriented;
    if (from == FibreColourFrom::Own && lines.ownColor.size() == lines.partColor.size()) return shownColour(lines.ownColor, v);
    if (from != FibreColourFrom::Fixed) return shownColour(lines.partColor, v);
    return population == 2 ? colours.furColour : population == 1 ? knit::kStrayFibreColour : knit::kOrientedFibreColour;
}

bool sameColours(const FibreColours& a, const FibreColours& b) {
    return a.fur == b.fur && a.oriented == b.oriented && a.stray == b.stray && a.furColour.x == b.furColour.x && a.furColour.y == b.furColour.y &&
           a.furColour.z == b.furColour.z;
}

float toSrgb(float c) { return float(knit::srgbEncode(std::clamp(c, 0.0f, 1.0f))); }

FibreCamera fibreCamera(const Mat4& projection, const Mat4& view) {
    FibreCamera camera;
    camera.clip = multiply(projection, view);
    return camera;
}

// The plane test of a box against clip space: false when the whole box is outside one plane. The
// planes come from the rows of the clip matrix (Gribb and Hartmann, "Fast extraction of viewing
// frustum planes from the world-view-projection matrix", 2001): row3 +- row0/1/2.
bool chunkInFrustum(const knit::FibreChunk& chunk, const Mat4& clip) {
    const double* m = clip.m;
    for (int axis = 0; axis < 3; ++axis)
        for (int sign = -1; sign <= 1; sign += 2) {
            const double a = m[3] + sign * m[axis], b = m[7] + sign * m[4 + axis], c = m[11] + sign * m[8 + axis],
                         d = m[15] + sign * m[12 + axis];
            // The box corner farthest along the plane normal.
            const double x = a >= 0 ? chunk.hi.x : chunk.lo.x, y = b >= 0 ? chunk.hi.y : chunk.lo.y, z = c >= 0 ? chunk.hi.z : chunk.lo.z;
            if (a * x + b * y + c * z + d < 0) return false;
        }
    return true;
}

FibreLines::~FibreLines() { release(); }

void FibreLines::release() {
    if (positions_) gl::DeleteBuffers(1, &positions_);
    if (colours_) gl::DeleteBuffers(1, &colours_);
    positions_ = colours_ = 0;
    owner_.reset();
    lines_ = nullptr;
    coloured_ = false;
}

void FibreLines::releaseUnless(const knit::BuiltScene* shown) {
    if (owner_ && owner_.get() != shown) release();
}

void FibreLines::upload(const knit::DisplayLines& lines) {
    gl::GenBuffers(1, &positions_);
    gl::BindBuffer(gl::kArrayBuffer, positions_);
    gl::BufferData(gl::kArrayBuffer, static_cast<std::ptrdiff_t>(lines.pos.size() * sizeof(float)), lines.pos.data(), gl::kStaticDraw);
    gl::GenBuffers(1, &colours_);
    gl::BindBuffer(gl::kArrayBuffer, 0);
}

void FibreLines::recolour(const knit::DisplayLines& lines, const FibreColours& colours) {
    const size_t n = lines.pos.size() / 3;
    const bool havePopulation = lines.population.size() == n && lines.partColor.size() == n * 3;
    std::vector<uint8_t> rgba(n * 4);
    for (size_t v = 0; v < n; ++v) {
        const Vec3 k = havePopulation ? fibreVertexColour(lines, v, colours) : knit::kOrientedFibreColour;
        rgba[4 * v] = toByte(k.x);
        rgba[4 * v + 1] = toByte(k.y);
        rgba[4 * v + 2] = toByte(k.z);
        rgba[4 * v + 3] = 255;
    }
    gl::BindBuffer(gl::kArrayBuffer, colours_);
    gl::BufferData(gl::kArrayBuffer, static_cast<std::ptrdiff_t>(rgba.size()), rgba.data(), gl::kStaticDraw);
    gl::BindBuffer(gl::kArrayBuffer, 0);
    coloured_ = true;
    colourKey_ = colours;
}

long long FibreLines::draw(const std::shared_ptr<const knit::BuiltScene>& owner, const knit::DisplayLines& lines,
                           const FibreColours& colours, const FibreCamera& camera, long long segmentBudget) {
    if (&lines != lines_) {
        release();
        owner_ = owner;
        lines_ = &lines;
        upload(lines);
    }
    if (!coloured_ || !sameColours(colourKey_, colours)) recolour(lines, colours);
    // 1. The chunks in the view frustum.
    visible_.clear();
    long long visibleSegments = 0;
    for (const knit::FibreChunk& c : lines.chunks) {
        if (camera.cull && !chunkInFrustum(c, camera.clip)) continue;
        visible_.push_back(&c);
        visibleSegments += c.vertexCount / 2;
    }
    // 2. The same share of every visible chunk (a prefix: the chunk's fibres are in random order).
    const double share = segmentBudget > 0 && visibleSegments > segmentBudget ? double(segmentBudget) / double(visibleSegments) : 1.0;
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    glLineWidth(1.0f);
    gl::BindBuffer(gl::kArrayBuffer, positions_);
    glVertexPointer(3, GL_FLOAT, 0, nullptr);
    gl::BindBuffer(gl::kArrayBuffer, colours_);
    glColorPointer(4, GL_UNSIGNED_BYTE, 0, nullptr);
    gl::BindBuffer(gl::kArrayBuffer, 0);
    // 3. One draw call per run of consecutive chunks.
    long long drawn = 0;
    GLint runFirst = 0;
    GLsizei runCount = 0;
    for (const knit::FibreChunk* c : visible_) {
        const GLsizei count = share < 1 ? static_cast<GLsizei>(c->vertexCount * share) & ~1 : static_cast<GLsizei>(c->vertexCount);
        if (count == 0) continue;
        if (runCount > 0 && runFirst + runCount != static_cast<GLint>(c->firstVertex)) {
            glDrawArrays(GL_LINES, runFirst, runCount);
            runCount = 0;
        }
        if (runCount == 0) runFirst = static_cast<GLint>(c->firstVertex);
        runCount += count;
        drawn += count / 2;
    }
    if (runCount > 0) glDrawArrays(GL_LINES, runFirst, runCount);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
    return drawn;
}

} // namespace viewer
