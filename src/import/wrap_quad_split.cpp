// How a quad is split into two triangles (quadSplit). In the flatten library, which the import
// component links, so the flatten and the stages above it share one definition.
#include "import/wrap_stages.h"

#include <algorithm>
#include <cmath>

namespace knit {

namespace {

constexpr double kDegenerateTriangle = 1e-12;  // a triangle's |cross| below this times the squared longer diagonal is degenerate
constexpr double kEqualDiagonals = 1e-12;      // diagonals whose relative difference is below this count as equal

// The quad's normal by Newell's method: the sum of the cross products of consecutive corners, written per component so that it
// does not depend on where the origin lies.
Vec3 newellNormal(const std::array<Vec3, 4>& c) {
    Vec3 n;
    for (int i = 0; i < 4; ++i) {
        const Vec3& a = c[i];
        const Vec3& b = c[(i + 1) % 4];
        n.x += (a.y - b.y) * (a.z + b.z);
        n.y += (a.z - b.z) * (a.x + b.x);
        n.z += (a.x - b.x) * (a.y + b.y);
    }
    return n;
}

// A split is valid when both its triangles have area and face the same way as the quad's normal.
bool splitIsValid(const std::array<Vec3, 4>& c, const std::array<std::array<int, 3>, 2>& split, const Vec3& normal, double longerDiagonal) {
    for (const std::array<int, 3>& t : split) {
        const Vec3 doubleArea = cross(c[t[1]] - c[t[0]], c[t[2]] - c[t[0]]);
        if (!(norm(doubleArea) > kDegenerateTriangle * longerDiagonal * longerDiagonal)) return false;
        if (!(dot(doubleArea, normal) > 0)) return false;
    }
    return true;
}

} // namespace

std::array<std::array<int, 3>, 2> quadSplit(const std::array<Vec3, 4>& corners) {
    const std::array<std::array<int, 3>, 2> alongFirst = {{{0, 1, 2}, {0, 2, 3}}};   // diagonal 0-2
    const std::array<std::array<int, 3>, 2> alongSecond = {{{0, 1, 3}, {1, 2, 3}}};  // diagonal 1-3
    const double first = norm(corners[2] - corners[0]);
    const double second = norm(corners[3] - corners[1]);
    const double longer = std::max(first, second);
    const bool equal = std::fabs(first - second) <= kEqualDiagonals * longer;
    const bool firstIsShorter = equal || first < second;
    const std::array<std::array<int, 3>, 2>& shorter = firstIsShorter ? alongFirst : alongSecond;
    const std::array<std::array<int, 3>, 2>& other = firstIsShorter ? alongSecond : alongFirst;
    const Vec3 normal = newellNormal(corners);
    if (splitIsValid(corners, shorter, normal, longer)) return shorter;
    if (splitIsValid(corners, other, normal, longer)) return other;
    return alongFirst;
}

} // namespace knit
