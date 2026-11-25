// The camera (camera.h).
#include "viewer/camera.h"

#include <algorithm>
#include <cmath>

namespace viewer {

using knit::cross;
using knit::dot;
using knit::normalized;

Mat4 Mat4::identity() {
    Mat4 r{};
    r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1;
    return r;
}

Mat4 multiply(const Mat4& a, const Mat4& b) {
    Mat4 r{};
    for (int col = 0; col < 4; ++col)
        for (int row = 0; row < 4; ++row)
            for (int k = 0; k < 4; ++k) r.m[4 * col + row] += a.m[4 * k + row] * b.m[4 * col + k];
    return r;
}

namespace {
Vec3 rotate(const Vec3& v, const Vec3& axis, double ang) {
    Vec3 k = normalized(axis);
    double c = std::cos(ang), s = std::sin(ang);
    return v * c + cross(k, v) * s + k * (dot(k, v) * (1 - c));
}
} // namespace

Vec3 Camera::right() const { return normalized(cross(forward, up)); }

Mat4 Camera::view() const {
    Vec3 f = normalized(forward), r = right(), u = cross(r, f);
    Vec3 e = eye();
    Mat4 v = Mat4::identity();
    v.m[0] = r.x; v.m[4] = r.y; v.m[8] = r.z;
    v.m[1] = u.x; v.m[5] = u.y; v.m[9] = u.z;
    v.m[2] = -f.x; v.m[6] = -f.y; v.m[10] = -f.z;
    v.m[12] = -dot(r, e); v.m[13] = -dot(u, e); v.m[14] = dot(f, e);
    return v;
}

Mat4 Camera::proj(double aspect) const {
    Mat4 p{};
    double n = std::max(0.01, dist * 0.01), f = farPlane();
    double t = std::tan(fovDeg * knit::kPi / 360.0);
    p.m[0] = 1 / (t * aspect); p.m[5] = 1 / t;
    p.m[10] = -(f + n) / (f - n); p.m[11] = -1; p.m[14] = -2 * f * n / (f - n);
    return p;
}

// Turntable: a horizontal drag turns around the world vertical (+z, the library's up), a
// vertical drag tilts between straight down and straight up. The camera's up is rebuilt from
// +z every step, so the view never rolls and the object stays upright.
void Camera::orbit(double dx, double dy, OrbitLock lock) {
    const double k = 0.006;
    const Vec3 z{0, 0, 1};
    if (lock == OrbitLock::Pitch) dx = 0;
    if (lock == OrbitLock::Yaw) dy = 0;
    const double pitch = std::asin(std::clamp(normalized(forward).z, -1.0, 1.0)) - dy * k;
    const double newPitch = std::clamp(pitch, -knit::kPi / 2, knit::kPi / 2);
    const Vec3 h = rotate(heading(), z, -dx * k);
    forward = normalized(h * std::cos(newPitch) + z * std::sin(newPitch));
    up = normalized(h * -std::sin(newPitch) + z * std::cos(newPitch));
}

// The horizontal direction the camera faces; looking straight down or up, the one its up
// vector shows (the Top view looks down with +y at the top of the screen).
Vec3 Camera::heading() const {
    Vec3 f = normalized(forward), u = normalized(up);
    Vec3 h{f.x, f.y, 0};
    if (knit::norm(h) < 1e-9) h = f.z < 0 ? Vec3{u.x, u.y, 0} : Vec3{-u.x, -u.y, 0};
    if (knit::norm(h) < 1e-9) h = {0, 1, 0};
    return normalized(h);
}

// Level the view: pitch 0, facing the same way, with +z up.
void Camera::level() {
    forward = heading();
    up = {0, 0, 1};
}

void Camera::pan(double dx, double dy, double viewportH) {
    double scale = 2 * dist * std::tan(fovDeg * knit::kPi / 360.0) / std::max(1.0, viewportH);
    Vec3 r = right(), u = cross(r, normalized(forward));
    target = target - r * (dx * scale) + u * (dy * scale);
}

void Camera::zoom(double wheel) { dist = std::clamp(dist * std::pow(0.88, wheel), 0.5, 2000.0); }

void Camera::look(const Vec3& fwd, const Vec3& upHint) {
    forward = normalized(fwd);
    Vec3 r = normalized(cross(forward, upHint));
    up = normalized(cross(r, forward));
}

Ray viewRay(const Camera& camera, const ViewPixel& pixel) {
    const Vec3 f = normalized(camera.forward), r = camera.right(), u = cross(r, f);  // as Camera::view
    const double t = std::tan(camera.fovDeg * knit::kPi / 360.0), aspect = pixel.width / pixel.height;
    const double right = (2 * pixel.x / pixel.width - 1) * t * aspect, up = (1 - 2 * pixel.y / pixel.height) * t;
    return {camera.eye(), f + r * right + u * up};
}

std::optional<CameraView> cameraViewNamed(const std::string& name) {
    if (name == "front") return CameraView::Front;
    if (name == "side") return CameraView::Side;
    if (name == "top") return CameraView::Top;
    if (name == "back") return CameraView::Back;
    if (name == "threequarter") return CameraView::ThreeQuarter;
    return std::nullopt;
}

void frameBox(Camera& cam, const Bounds& box) {
    cam.target = (box.lo + box.hi) * 0.5;
    cam.dist = std::max(1.6 * knit::norm(box.hi - box.lo), 5.0);
}

void lookFrom(Camera& cam, CameraView view) {
    switch (view) {
    case CameraView::Side: cam.look({-1, 0, 0}, {0, 0, 1}); break;
    case CameraView::Top: cam.look({0, 0, -1}, {0, 1, 0}); break;
    case CameraView::Back: cam.look({0, -1, 0}, {0, 0, 1}); break;
    case CameraView::Front: cam.look({0, 1, 0}, {0, 0, 1}); break;  // glTF -Z front / Z-up front view -> library +y view
    case CameraView::ThreeQuarter: {
        // 30 degrees from Front towards the side the Side preset shows, 15 degrees down.
        const double turn = 30 * knit::kPi / 180, tilt = 15 * knit::kPi / 180;
        cam.look({-std::sin(turn) * std::cos(tilt), std::cos(turn) * std::cos(tilt), -std::sin(tilt)}, {0, 0, 1});
        break;
    }
    }
}

std::optional<CameraView> presetOf(const Camera& cam) {
    for (const CameraView view : {CameraView::Front, CameraView::Side, CameraView::Top, CameraView::Back, CameraView::ThreeQuarter}) {
        Camera preset;
        lookFrom(preset, view);
        if (knit::norm(normalized(cam.forward) - preset.forward) < 1e-6 && knit::norm(normalized(cam.up) - preset.up) < 1e-6) return view;
    }
    return std::nullopt;
}

std::optional<knit::Vec2> CameraDrag::takeClick() {
    const std::optional<knit::Vec2> click = click_;
    click_.reset();
    return click;
}

bool CameraDrag::update(Camera& cam, const MouseState& mouse, double viewportHeight) {
    bool moved = false;
    for (int b = 0; b < 3; ++b) {
        const bool pressed = mouse.down[b] && !held_[b];
        const bool on = mouse.down[b] && (held_[b] || mouse.overViewport);
        if (pressed && on && b == 0) {
            lock_ = OrbitLock::Free;
            lockDx_ = lockDy_ = 0;
            pressX_ = mouse.x;
            pressY_ = mouse.y;
            dragged_ = false;
        }
        if (b == 0 && held_[0]) {
            dragged_ = dragged_ || std::abs(mouse.x - pressX_) + std::abs(mouse.y - pressY_) > 4;
            if (!mouse.down[0] && !dragged_) click_ = knit::Vec2{mouse.x, mouse.y};
        }
        if (on && held_[b]) {
            double dx = mouse.x - lastX_, dy = mouse.y - lastY_;
            if (b == 0 && mouse.shift && lock_ == OrbitLock::Free) {
                lockDx_ += dx;
                lockDy_ += dy;
                if (std::abs(lockDx_) + std::abs(lockDy_) >= 4) {
                    lock_ = std::abs(lockDx_) >= std::abs(lockDy_) ? OrbitLock::Yaw : OrbitLock::Pitch;
                    dx = lockDx_;
                    dy = lockDy_;
                } else {
                    dx = dy = 0;
                }
            }
            if (b == 0) cam.orbit(dx, dy, mouse.shift ? lock_ : OrbitLock::Free);
            else cam.pan(dx, dy, viewportHeight);
            moved = moved || dx != 0 || dy != 0;
        }
        held_[b] = on;
    }
    lastX_ = mouse.x;
    lastY_ = mouse.y;
    if (mouse.overViewport && mouse.wheel != 0) {
        cam.zoom(mouse.wheel);
        moved = true;
    }
    return moved;
}

void CameraRig::show(const std::string& viewKey, bool hasContent, const Bounds& content, CameraView startView) {
    if (viewKey != viewKey_) {  // another view: its own camera
        parked_[viewKey_] = camera_;
        if (parked_.count(viewKey)) camera_ = parked_.at(viewKey);
        viewKey_ = viewKey;
    }
    if (framed_.count(viewKey) || !hasContent) return;
    lookFrom(camera_, startView);
    if (viewKey == "reference-wool") camera_.look({0, 0.6, -0.8}, {0, 0, 1});  // the flat patch from above, obliquely
    frameBox(camera_, content);
    framed_.insert(viewKey);
}

void CameraRig::apply(CameraRequests& requests, const Bounds& content, bool hasContent) {
    if (requests.view) {
        lookFrom(camera_, *requests.view);
        frameBox(camera_, content);
        requests.view.reset();
    }
    if (requests.upright) {
        camera_.level();
        requests.upright = false;
    }
    if (requests.frame && hasContent) {
        frameBox(camera_, content);
        requests.frame = false;
    }
    if (requests.focus) {
        frameBox(camera_, *requests.focus);
        requests.focus.reset();
    }
}

// While the camera moves (and this long after), the fibres are thinned.
constexpr double kMovingHoldSeconds = 0.15;

bool CameraRig::moving(double now) const { return movedAt_ >= 0 && now - movedAt_ < kMovingHoldSeconds; }

} // namespace viewer
