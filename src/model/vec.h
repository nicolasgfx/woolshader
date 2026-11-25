#pragma once
// Minimal 2D/3D vector types. Units are millimetres throughout the library; coordinates.h says where a point is.
#include <algorithm>
#include <cmath>
#include <type_traits>

namespace knit {

struct Vec3 {
    double x = 0, y = 0, z = 0;
    constexpr Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
    constexpr Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vec3 operator-() const { return {-x, -y, -z}; }
    constexpr Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    constexpr Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
    constexpr bool operator==(const Vec3&) const = default;
};
inline constexpr Vec3 operator*(double s, const Vec3& v) { return v * s; }
inline constexpr double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
inline Vec3 normalized(const Vec3& a) {
    double n = norm(a);
    return n > 0 ? a / n : Vec3{};
}
inline Vec3 lerp(const Vec3& a, const Vec3& b, double t) { return a + (b - a) * t; }
inline Vec3 componentMin(const Vec3& a, const Vec3& b) { return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)}; }
inline Vec3 componentMax(const Vec3& a, const Vec3& b) { return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)}; }

// Digests hash Vec3's raw bytes (curvesDigest, meshDigest, the knit cache): its layout is frozen.
static_assert(sizeof(Vec3) == 3 * sizeof(double) && std::is_trivially_copyable_v<Vec3>);

struct Vec2 {
    double x = 0, y = 0;
    constexpr Vec2() = default;
    constexpr Vec2(double x_, double y_) : x(x_), y(y_) {}
    constexpr Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }
    constexpr Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }
    constexpr Vec2 operator*(double s) const { return {x * s, y * s}; }
};
inline constexpr double dot(const Vec2& a, const Vec2& b) { return a.x * b.x + a.y * b.y; }
inline constexpr double cross(const Vec2& a, const Vec2& b) { return a.x * b.y - a.y * b.x; }
inline double norm(const Vec2& a) { return std::sqrt(dot(a, a)); }

constexpr double kPi = 3.14159265358979323846;

} // namespace knit
