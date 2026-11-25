#pragma once
// The viewer's camera, with no GLFW and no GL: the matrices, the turntable moves, what a frame of mouse
// does to them, and one camera per view. main.cpp reads the mouse, renderer.cpp reads the matrices.
#include <map>
#include <optional>
#include <set>
#include <string>

#include "model/vec.h"

namespace viewer {

using knit::Vec3;

struct Mat4 {
    double m[16];  // column-major
    static Mat4 identity();
};
Mat4 multiply(const Mat4& a, const Mat4& b);

struct Bounds {
    Vec3 lo, hi;
};

// Camera::orbit's Shift-drag: the drag moves along one axis only.
enum class OrbitLock { Free, Yaw, Pitch };

struct Camera {
    Vec3 target{0, 0, 0};
    Vec3 forward{0, 0, -1};  // eye -> target
    Vec3 up{0, 1, 0};
    double dist = 80;
    double fovDeg = 35;

    Vec3 eye() const { return target - forward * dist; }
    Vec3 right() const;
    Mat4 view() const;
    Mat4 proj(double aspect) const;
    double farPlane() const { return dist * 10 + 100; }  // proj's far clipping distance
    Vec3 heading() const;
    void orbit(double dxPix, double dyPix, OrbitLock lock = OrbitLock::Free);
    void level();
    void pan(double dxPix, double dyPix, double viewportH);
    void zoom(double wheel);
    void look(const Vec3& fwd, const Vec3& upHint);
};

// The view presets. ThreeQuarter: from the front and the Side preset's side at once, looking a little down.
enum class CameraView { Front, Side, Top, Back, ThreeQuarter };
// "front", "side", "top", "back", "threequarter"
std::optional<CameraView> cameraViewNamed(const std::string& name);

// A pixel inside the viewport, top-left origin, with the viewport's size.
struct ViewPixel {
    double x = 0, y = 0, width = 1, height = 1;
};

struct Ray {
    Vec3 origin, direction;
};

Ray viewRay(const Camera& camera, const ViewPixel& pixel);  // from the eye through the pixel

void frameBox(Camera& cam, const Bounds& box);
void lookFrom(Camera& cam, CameraView view);
// The preset whose orientation the camera has (forward and up), none after an orbit; pan and zoom keep it.
std::optional<CameraView> presetOf(const Camera& cam);

// What the interface asks of the camera; applied after the frame, when the content is current.
struct CameraRequests {
    std::optional<CameraView> view;
    bool frame = false;    // stays set until there is content to frame
    bool upright = false;  // level the camera (pitch 0, +z up), keeping its heading
    std::optional<Bounds> focus;
};

// One frame of the mouse, in framebuffer pixels.
struct MouseState {
    double x = 0, y = 0, wheel = 0;
    bool down[3] = {false, false, false};  // left, right, middle (GLFW's order)
    bool shift = false;
    bool overViewport = false;             // inside the viewport and not over the interface
};

// Turns the mouse into camera moves. A drag acts when it starts over the viewport: left orbits
// (Shift: its first 4 pixels choose one axis), right and middle pan; the wheel zooms.
class CameraDrag {
public:
    bool update(Camera& cam, const MouseState& mouse, double viewportHeight);  // true: the camera moved
    // A left press over the viewport released within 4 pixels, where it was released, once. An orbit is no click.
    std::optional<knit::Vec2> takeClick();

private:
    bool held_[3] = {false, false, false};  // down, and it started over the viewport
    double lastX_ = 0, lastY_ = 0;
    OrbitLock lock_ = OrbitLock::Free;      // the Shift-drag's axis, once chosen
    double lockDx_ = 0, lockDy_ = 0;        // the Shift-drag's motion before the axis is chosen
    double pressX_ = 0, pressY_ = 0;        // where the left button went down
    bool dragged_ = false;                  // it has since moved too far to be a click
    std::optional<knit::Vec2> click_;
};

// One camera per view (ViewContent::viewKey): framed automatically the first time it shows content;
// afterwards only Frame and the view buttons move it.
class CameraRig {
public:
    Camera& camera() { return camera_; }
    const Camera& camera() const { return camera_; }
    void show(const std::string& viewKey, bool hasContent, const Bounds& content, CameraView startView);
    void apply(CameraRequests& requests, const Bounds& content, bool hasContent);
    void noteMoved(double now) { movedAt_ = now; }
    bool moving(double now) const;  // the fibres thin while it moves, and shortly after

private:
    Camera camera_;
    std::map<std::string, Camera> parked_;
    std::set<std::string> framed_;
    std::string viewKey_ = "model";
    double movedAt_ = -1;  // [s] the last orbit, pan or zoom
};

} // namespace viewer
