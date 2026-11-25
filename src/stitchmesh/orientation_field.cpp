#include "stitchmesh/orientation_field.h"

#include <array>
#include <chrono>
#include <cmath>

#include "rim_remesh.h"

namespace knit {

bool computeOrientationField(const TriMesh& m, const OrientationFieldRequest& r, OrientationField& out, std::string* err) {
    out = OrientationField{};
    const auto start = std::chrono::steady_clock::now();
    std::vector<std::array<double, 3>> pos;
    std::vector<std::array<int, 3>> tri(m.tri.begin(), m.tri.end());
    for (const Vec3& p : m.pos) pos.push_back({p.x, p.y, p.z});
    knit_rim::Options o;
    o.spacing = r.spacing;
    o.rosy = r.rosy;
    o.extrinsic = r.extrinsic;
    o.sweepsPerLevel = r.sweepsPerLevel;
    o.seed = r.seed;
    o.fieldsOnly = true;
    for (const Vec3& g : r.guide) o.orientation.push_back({g.x, g.y, g.z});
    knit_rim::Result res;
    if (!knit_rim::remesh(pos, tri, o, res, err)) return false;
    out.log = res.log;
    out.levelSizes = res.levelSizes;
    out.finestVertices = res.inputVertices;
    if (res.inputVertices != (int) m.pos.size()) {
        if (err) *err = "orientation field: RIM subdivided the input (" + std::to_string(res.inputVertices) + " vertices instead of " +
                        std::to_string(m.pos.size()) + "), so the field is not per input vertex";
        return false;
    }
    for (int i = 0; i < (int) m.pos.size(); ++i) {
        const Vec3 n = m.vnormal.empty() ? Vec3{} : m.vnormal[i];
        Vec3 q{res.fieldQ[i][0], res.fieldQ[i][1], res.fieldQ[i][2]};
        q = q - n * dot(n, q);
        out.q.push_back(normalized(q));
    }
    for (const auto& s : res.orientationSingularPoints) out.singularities.push_back({s[0], s[1], s[2]});
    out.ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    return true;
}

std::vector<Vec3> ringAxisGuide(const TriMesh& m, const Vec3& axis) {
    const Vec3 a = normalized(axis);
    std::vector<Vec3> out;
    for (size_t i = 0; i < m.pos.size(); ++i) {
        const Vec3 radial = m.pos[i] - a * dot(a, m.pos[i]);
        Vec3 around = cross(a, radial);
        if (!m.vnormal.empty()) around = around - m.vnormal[i] * dot(m.vnormal[i], around);
        out.push_back(normalized(around));
    }
    return out;
}

double directionDeviation(const Vec3& a, const Vec3& b, const Vec3& normal, double symmetryDegrees) {
    const Vec3 n = normalized(normal);
    const Vec3 ta = normalized(a - n * dot(n, a)), tb = normalized(b - n * dot(n, b));
    double deg = std::atan2(dot(cross(ta, tb), n), dot(ta, tb)) * 180.0 / kPi;   // signed, (-180, 180]
    deg = std::fmod(std::fabs(deg), symmetryDegrees);
    return std::min(deg, symmetryDegrees - deg);
}

} // namespace knit
