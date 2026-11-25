// The static settle of a fibre on a mesh and the stuck test (fur.h has the method and its sources).
#include "fibres/fur.h"

#include <algorithm>
#include <cmath>

#include "fibres/strands_file.h"
#include "model/mesh_query.h"

namespace knit {

namespace {

constexpr double kGravity = 9810.0;  // [mm/s^2]

// The minimal rotation taking unit a to unit b, applied to v (v unchanged when a == b; a
// half turn about any axis normal to a when b == -a).
Vec3 minimalRotation(const Vec3& a, const Vec3& b, const Vec3& v) {
    const Vec3 axis = cross(a, b);
    const double s = norm(axis), c = dot(a, b);
    if (s < 1e-12) {
        if (c > 0) return v;
        Vec3 k = normalized(cross(a, std::fabs(a.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}));
        return rotateAbout(v, k, kPi);
    }
    return rotateAbout(v, axis / s, std::atan2(s, c));
}

// Moves p out of the body, to `clearance` from its closest point along that face's normal (the local
// pseudo-normal side decides inside, MeshQuery::closest). Whether it moved.
bool pushOutOfBody(const MeshQuery& body, Vec3& p, double clearance, double reach) {
    const MeshQuery::Hit hit = body.closest(p, reach);
    if (hit.tri < 0) return false;
    if (!hit.inside && hit.dist >= clearance) return false;
    p = hit.p + body.faceNormals()[static_cast<size_t>(hit.tri)] * clearance;
    return true;
}

bool insideBody(const MeshQuery& body, const Vec3& p, double reach) {
    const MeshQuery::Hit hit = body.closest(p, reach);
    return hit.tri >= 0 && hit.inside;
}

// One strand's static settle (see fur.h): x[0], x[1] stay, the rest settles. radii: per vertex; offset: the
// body's clearance beyond each radius. Returns the largest last-step move [mm].
double settleStrand(std::vector<Vec3>& x, const std::vector<double>& radii, const StrandSettleParams& s, const MeshQuery& body,
                    double offset) {
    const int nv = static_cast<int>(x.size());
    if (nv < 3) return 0;
    // Rest state: segment lengths and directions (the groom), and the shape stiffness per vertex.
    std::vector<double> L(static_cast<size_t>(nv), 0), k(static_cast<size_t>(nv), 0);
    std::vector<Vec3> R(static_cast<size_t>(nv));
    for (int i = 1; i < nv; ++i) {
        L[i] = norm(x[i] - x[i - 1]);
        R[i] = (x[i] - x[i - 1]) / L[i];
    }
    for (int i = 2; i < nv; ++i) k[i] = s.stiffness + (s.tipStiffness - s.stiffness) * double(i - 2) / std::max(1, nv - 3);
    const double dt = s.timeStep;
    const Vec3 g{0, 0, -kGravity * s.gravity};
    std::vector<Vec3> v(static_cast<size_t>(nv)), p(static_cast<size_t>(nv)), d(static_cast<size_t>(nv) + 1);
    double residual = 0;
    for (int step = 0; step < s.steps; ++step) {
        // Prediction (the two root vertices are fixed).
        p[0] = x[0];
        p[1] = x[1];
        for (int i = 2; i < nv; ++i) {
            v[i] = v[i] * s.velocityKeep + g * dt;
            p[i] = x[i] + v[i] * dt;
        }
        // Local shape: toward the rest direction carried by the previous segment's rotation.
        for (int i = 2; i < nv; ++i) {
            const Vec3 prevCur = normalized(p[i - 1] - p[i - 2]);
            const Vec3 goal = p[i - 1] + minimalRotation(R[i - 1], prevCur, R[i]) * L[i];
            p[i] += (goal - p[i]) * k[i];
        }
        // Follow-The-Leader with the body: from the root outward, each vertex at its exact length
        // from its (final) parent, pushed out of the body and put back at that length (a few
        // rounds: the push and the length agree once the vertex rests against the surface).
        // d = the total correction of the vertex.
        followTheLeader(p, L, radii, body, offset, &d);
        d[nv] = Vec3{};
        // DFTL velocity (Mueller et al. 2012): the next vertex's correction is damped out.
        residual = 0;
        for (int i = 2; i < nv; ++i) {
            v[i] = (p[i] - x[i]) / dt - d[i + 1] * (s.damping / dt);
            residual = std::max(residual, norm(p[i] - x[i]));
            x[i] = p[i];
        }
    }
    return residual;
}

std::vector<double> hairRadii(const Hair& h) {
    const int nv = static_cast<int>(h.points.size());
    std::vector<double> radii(static_cast<size_t>(nv));
    for (int i = 0; i < nv; ++i) radii[i] = taperedRadius(h.rootRadius, h.tipRadius, i, nv - 1);
    return radii;
}

} // namespace

Vec3 rotateAbout(const Vec3& v, const Vec3& k, double a) {
    const double c = std::cos(a), s = std::sin(a);
    return v * c + cross(k, v) * s + k * (dot(k, v) * (1 - c));
}

void followTheLeader(std::vector<Vec3>& x, const std::vector<double>& L, const std::vector<double>& radii, const MeshQuery& body,
                     double offset, std::vector<Vec3>* d) {
    const int nv = static_cast<int>(x.size());
    for (int i = 2; i < nv; ++i) {
        const Vec3 before = x[i];
        x[i] = x[i - 1] + normalized(x[i] - x[i - 1]) * L[i];
        const double reach = radii[i] + offset + L[i];
        for (int round = 0; round < 3 && pushOutOfBody(body, x[i], radii[i] + offset, reach); ++round)
            x[i] = x[i - 1] + normalized(x[i] - x[i - 1]) * L[i];
        if (d) (*d)[i] = x[i] - before;
    }
}

double visibleLength(const Hair& h) {
    double L = 0;
    for (size_t i = 2; i < h.points.size(); ++i) L += norm(h.points[i] - h.points[i - 1]);
    return L;
}

bool hairInsideBody(const Hair& h, const MeshQuery& body) {
    const int nv = static_cast<int>(h.points.size());
    const std::vector<double> radii = hairRadii(h);
    for (int i = 2; i < nv; ++i) {
        const double seg = norm(h.points[i] - h.points[i - 1]);
        if (insideBody(body, h.points[i], radii[i] + seg)) return true;
        // The first visible segment starts on the surface: a crossing at its root does not count.
        if (body.segmentCrosses(h.points[i - 1], h.points[i], i == 2 ? kRootCrossingSkip : 0.0)) return true;
    }
    return false;
}

bool settleHair(Hair& h, const StrandSettleParams& settle, double gravityFullLength, double collisionOffset, const MeshQuery& body,
                double gravityLength) {
    if (h.points.size() < 3) return true;
    StrandSettleParams s = settle;
    if (gravityFullLength > 0) s.gravity *= std::min(1.0, std::pow(gravityLength / gravityFullLength, 3.0));
    h.residual = settleStrand(h.points, hairRadii(h), s, body, collisionOffset);
    // Final check: no vertex or segment point on the body's inner side (local pseudo-normal of the closest feature).
    return !hairInsideBody(h, body);
}

} // namespace knit
