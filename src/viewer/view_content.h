#pragma once
// What the viewport shows: the window area it is drawn in, the meshes and fibres of the view on show with
// their caption, which of the toolbar's stages can be shown and which Display rows apply, the comparison
// with the source, the scenes whose buffers may stay on the GPU, the counts overlay, and picking (a click's
// ray taken back through the import preview's transform). No ImGui and no GL: the renderer, main.cpp and
// the interface read it; view_content.cpp holds the code (const SceneUi& only).
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "viewer/workflow.h"

namespace viewer {

// The viewport: the window area the model is drawn in, between the panel, the workflow strip +
// view toolbar and the status bar + the viewport's own status line, in the units drawSceneUi lays out
// in (origin top left, y down). The renderer, the mouse handling and the viewport labels use it.
struct ViewportRect {
    float x = 0, y = 0, w = 1, h = 1;
    bool contains(float px, float py) const { return px >= x && px < x + w && py >= y && py < y + h; }
};
ViewportRect viewportRect(const SceneUi& ui, float width, float height);

// The import preview's file -> library transform: scale, the up axis, centring.
struct PreviewTransform {
    double scale = 1;
    bool yUp = true;
    knit::Vec3 offset;
};
// A point in the imported file's coordinates placed in the library as the viewport draws it (the inverse of
// knit::toFileCoordinates with the same scale, offset and up axis).
knit::Vec3 toLibrary(const PreviewTransform& p, const knit::Vec3& v);

// What the viewport shows and how.
struct ViewContent {
    std::vector<const knit::DisplayMesh*> meshes;
    std::vector<float> meshAlpha;                     // per mesh, 1 opaque (empty: all opaque); the Mark step's see-through pieces
    const knit::BuiltScene* meshOwner = nullptr;      // the scene the meshes' albedo indices refer to
    const OverlayGeometry* overlay = nullptr;         // lines and points over the meshes (the Mark step's markings)
    const knit::DisplayLines* fibres = nullptr;
    std::shared_ptr<const knit::FurPreview> furDiscs;  // the Fur step's disc preview (fur_discs.h)
    std::optional<PreviewTransform> preview;          // the import preview only
    Bounds bounds{{-100, -100, 0}, {100, 100, 200}};  // library mm
    std::string caption;
    // The view the content belongs to: each has its own camera ("model", "reference-wool",
    // "reference-fur"), framed automatically the first time only.
    std::string viewKey = "model";
};
// The view on show: the open import file, the Mark step's pieces, the Fibres step's reference sample, else the chosen stage.
ViewContent viewContent(const SceneUi& ui);
// The model view of one stage, whatever stage is chosen: Mesh is the source (the import preview or the
// result's base meshes, never fibres), Knit the yarn, Fibres the yarn with its fibres or the hair.
ViewContent stageContent(const SceneUi& ui, ViewMode mode);
// The wool fibres the Fibres view draws: only those of the knit on show, never over a partial build
// (its fibres belong to the previous knit).
const knit::FibreResult* shownFibres(const SceneUi& ui);

// --- the toolbar's stages and Display rows

// The stage the segment marks as shown: Mesh while an import file is open (its preview is a source),
// Fibres for a Knit chosen before the Knit step was skipped (the hair), else the chosen one. Neither
// changes the chosen stage.
ViewMode markedViewMode(const SceneUi& ui);

// Whether the Source / Knit / Fibres button can show its stage, and why not.
struct ViewModeState {
    bool available = false;
    std::string reason;  // "" when available
};
ViewModeState viewModeState(const SceneUi& ui, ViewMode mode);

// The rows of the Display popover that apply to what the viewport shows. A row that does not apply is
// hidden and keeps its value.
struct DisplayRows {
    bool occlusion = false;         // Studio chosen
    bool plyTexture = false;        // the knit on show has yarn
    std::string plyTextureWhyNot;   // "" usable
    bool yarnMaps = false;          // the knit on show carries baked cavity maps on its tubes
    bool unknittedParts = false;    // the knit on show has unknitted parts (the scene's, not the filtered view's)
    bool fibrePopulations = false;  // wool fibres on the model
    bool compare = false;           // Knit or Fibres on the model view
    std::string compareWhyNot;      // "" when the source and the stage can both be drawn
};
DisplayRows displayRows(const SceneUi& ui);

// What the viewport draws this frame. One camera and the whole viewport's projection; while comparing,
// only the columns each half may write differ. `bounds` covers both (Fit, the first framing).
struct ViewportFrame {
    ViewContent shown;                  // the whole viewport, or its right half while comparing
    std::optional<ViewContent> source;  // the left half: only while the comparison can be drawn
    Bounds bounds;
};
ViewportFrame viewportFrame(const SceneUi& ui);

// The viewport's status line, left side: the object's name and its size.
struct ObjectLine {
    std::string name;  // the open or imported file's name, else the scene's label
    std::string size;  // "140 x 326 x 84 mm"; "" while there is nothing to measure
};
ObjectLine objectLine(const SceneUi& ui);

// Where a toolbar popover goes: the whole rectangle inside the viewport (interface units). It starts just
// inside the viewport's top, right-aligned with its button where it fits; its content scrolls past maxHeight.
struct PopoverRect {
    float x = 0, y = 0, width = 1, maxHeight = 1;
};
PopoverRect popoverRect(const ViewportRect& viewport, float anchorRight, float preferredWidth, float margin);

// --- the import preview, the GPU's scenes, the counts

// An open file's component left out of the knit (nodeKeptAsMesh): the viewport draws it grey.
bool previewNodeIsPlain(const SceneUi& ui, const knit::DisplayMesh& d);
// An open file's focused component (ImportFlow::focusedNode): the viewport outlines it.
bool previewNodeIsFocused(const SceneUi& ui, const knit::DisplayMesh& d);
std::optional<int> nodeUnderPixel(const SceneUi& ui, const Camera& camera, const ViewPixel& pixel);  // index into ImportedMesh::nodes
// The Mark step's pick: the nearest visible piece under the pixel, the point as a file triangle with barycentric
// weights, and its position (file coordinates).
struct SurfaceHit {
    int piece = -1;
    knit::SurfacePoint point;
    knit::Vec3 position;
};
std::optional<SurfaceHit> surfaceUnderPixel(const SceneUi& ui, const Camera& camera, const ViewPixel& pixel);
// The scenes whose meshes and textures may be on screen. Parked results are not listed: they stay in
// memory, not on the GPU.
struct MeshOwners {
    std::vector<std::shared_ptr<const knit::BuiltScene>> scenes;
    std::shared_ptr<const knit::FibrePreview> preview;
};
MeshOwners meshOwners(const SceneUi& ui);
// The canvas overlay's counts: the triangles shown (by kind), the fibres / hairs generated (with
// their segments) and how many are shown; the reference patch's own counts are labelled as such.
std::string viewCounts(const SceneUi& ui, const ViewContent& c);

} // namespace viewer
