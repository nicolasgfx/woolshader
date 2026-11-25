#pragma once
// GPU backend: the surface model (rings, ringAt, evalSurface, yarnCoordinates), the
// SurfaceIndex queries and the mesh queries, ported from knitting/tube/yarn_surface.cpp / model/mesh_query.cpp / knitting/layout/chart.cpp
// in single precision (math.cuh: precision policy).
#include "gpu/geometry/math.cuh"

namespace knit::gpu {

struct DeviceRing {
    double s, u;
    f3 c, T, n1, n2;
    float a, b, phi;
};
__device__ __forceinline__ DeviceRing loadRing(const RingRecord& r) {
    DeviceRing o;
    o.s = r.s;
    o.u = r.u;
    o.c = fromF(r.c);
    o.T = fromF(r.T);
    o.n1 = fromF(r.n1);
    o.n2 = fromF(r.n2);
    o.a = r.a;
    o.b = r.b;
    o.phi = r.phi;
    return o;
}

__device__ __forceinline__ float wrapAxis(float x) {
    while (x > 0.5f * kPiF) x -= kPiF;
    while (x <= -0.5f * kPiF) x += kPiF;
    return x;
}
__device__ __forceinline__ void ringAxes(const DeviceRing& r, f3& ea, f3& eb) {
    float s, c;
    sincosf(r.phi, &s, &c);
    ea = r.n1 * c + r.n2 * s;
    eb = r.n2 * c - r.n1 * s;
}

// One yarn's surface on the device: its curve and its rings.
struct YarnSurfaceView {
    CurveView curve;
    const RingRecord* rings;
    int nRings;
    double ds, length;
    float radius, ext;
};
__device__ __forceinline__ YarnSurfaceView devYarnSurf(const CurveSetView& cs, const SurfaceView& sm, int y) {
    const YarnSurfaceRecord& ys = sm.yarns[y];
    return {devCurve(cs, y), sm.rings + ys.ringOff, ys.nRings, ys.ds, ys.length, ys.radius, ys.ext};
}

// yarn_surface.cpp ringAt: interpolated frame and axes, centre and tangent on the spline at the
// arc-length parameter (paramAtArcLength).
__device__ inline DeviceRing ringAt(const YarnSurfaceView& ys, double s) {
    const int N = ys.nRings - 1;
    s = clampd(s, 0.0, ys.length);
    int k = static_cast<int>(floor(s / ys.ds));
    k = k < 0 ? 0 : (k > N - 1 ? N - 1 : k);
    const DeviceRing A = loadRing(ys.rings[k]), B = loadRing(ys.rings[k + 1]);
    const float f = static_cast<float>(clampd((s - A.s) / ys.ds, 0.0, 1.0));
    DeviceRing r;
    r.s = s;
    r.u = paramAtArcLength(ys.curve, s);
    r.c = curveEval(ys.curve, r.u);
    r.T = normalized(curveD1(ys.curve, r.u));
    const f3 n = A.n1 * (1 - f) + B.n1 * f;
    r.n1 = normalized(n - r.T * dot(n, r.T));
    r.n2 = cross(r.T, r.n1);
    r.a = A.a * (1 - f) + B.a * f;
    r.b = A.b * (1 - f) + B.b * f;
    r.phi = A.phi + f * wrapAxis(B.phi - A.phi);
    return r;
}

struct DeviceSurfacePoint {
    f3 p, n;
    DeviceRing ring;
};
// yarn_surface.cpp evalSurface.
__device__ inline DeviceSurfacePoint evalSurface(const YarnSurfaceView& ys, double s, float theta) {
    DeviceSurfacePoint sp;
    sp.ring = ringAt(ys, s);
    f3 ea, eb;
    ringAxes(sp.ring, ea, eb);
    float sn, c;
    sincosf(theta - sp.ring.phi, &sn, &c);
    sp.p = sp.ring.c + ea * (sp.ring.a * c) + eb * (sp.ring.b * sn);
    sp.n = normalized(ea * (c / sp.ring.a) + eb * (sn / sp.ring.b));
    return sp;
}

struct YarnCoordinates {
    f3 X, dS, dTheta, dLambda;
};
// yarn_surface.cpp yarnCoordinates (u interpolated linearly between the model rings' u, in double).
__device__ inline YarnCoordinates yarnCoordinates(const YarnSurfaceView& ys, double s, float theta, float lambda) {
    const int N = ys.nRings - 1;
    s = clampd(s, 0.0, ys.length);
    int k = static_cast<int>(floor(s / ys.ds));
    k = k < 0 ? 0 : (k > N - 1 ? N - 1 : k);
    const DeviceRing A = loadRing(ys.rings[k]), B = loadRing(ys.rings[k + 1]);
    const double fd = clampd((s - A.s) / ys.ds, 0.0, 1.0);
    const float f = static_cast<float>(fd);
    const float dus = static_cast<float>((B.u - A.u) / ys.ds);
    const double u = A.u + fd * (B.u - A.u);
    const f3 d1 = curveD1(ys.curve, u), d2 = curveD2(ys.curve, u);
    const f3 c = curveEval(ys.curve, u);
    const float sp = norm(d1);
    const f3 T = d1 / sp, dT = (d2 - T * dot(d2, T)) * (dus / sp), dc = d1 * dus;
    const float invDs = static_cast<float>(1.0 / ys.ds);
    const f3 m = A.n1 * (1 - f) + B.n1 * f, dm = (B.n1 - A.n1) * invDs;
    const float mT = dot(m, T);
    const f3 w = m - T * mT, dw = dm - dT * mT - T * (dot(dm, T) + dot(m, dT));
    const float wl = norm(w);
    const f3 n1 = w / wl, dn1 = (dw - n1 * dot(n1, dw)) / wl;
    const f3 n2 = cross(T, n1), dn2 = cross(dT, n1) + cross(T, dn1);
    const float a = A.a * (1 - f) + B.a * f, b = A.b * (1 - f) + B.b * f;
    const float da = (B.a - A.a) * invDs, db = (B.b - A.b) * invDs;
    const float dphi = wrapAxis(B.phi - A.phi) * invDs, phi = A.phi + f * wrapAxis(B.phi - A.phi);
    float sph, cp;
    sincosf(phi, &sph, &cp);
    const f3 ea = n1 * cp + n2 * sph, eb = n2 * cp - n1 * sph;
    const f3 dea = dn1 * cp + dn2 * sph + eb * dphi, deb = dn2 * cp - dn1 * sph - ea * dphi;
    float si, co;
    sincosf(theta - phi, &si, &co);
    YarnCoordinates out;
    out.dLambda = ea * (a * co) + eb * (b * si);
    out.X = c + out.dLambda * lambda;
    out.dTheta = (eb * (b * co) - ea * (a * si)) * lambda;
    out.dS = dc + (dea * (a * co) + ea * (da * co + a * si * dphi) + deb * (b * si) + eb * (db * si - b * co * dphi)) * lambda;
    return out;
}

// --- SurfaceIndex
struct Exclude {
    int yarn = -1;
    double s0 = 0, s1 = -1;
    __device__ bool covers(int y, double s) const { return y == yarn && s >= s0 && s <= s1; }
};

struct SurfaceIndexQuery {
    SurfaceIndexView idx;
    SurfaceView sm;
};

__device__ __forceinline__ bool pointInEnclosure(const SurfaceIndexQuery& I, f3 q, float inflate, const Exclude& ex) {
    bool hit = false;
    for (GridCursor cursor(I.idx.grid, q - f3{inflate, inflate, inflate}, q + f3{inflate, inflate, inflate}); cursor.next();) {
        const int i = cursor.item();
        const CapsuleRecord& c = I.idx.caps[i];
        float t;
        const float d = pointSegment(q, fromF(c.a), fromF(c.b), t);
        if (ex.covers(c.yarn, c.s0 + t * (c.s1 - c.s0))) continue;
        if (d < c.R + inflate) { hit = true; break; }
    }
    return hit;
}

__device__ __forceinline__ bool segmentHitsEnclosure(const SurfaceIndexQuery& I, f3 a, f3 b, float inflate, const Exclude& ex) {
    bool hit = false;
    const f3 g{inflate, inflate, inflate};
    for (GridCursor cursor(I.idx.grid, vmin(a, b) - g, vmax(a, b) + g); cursor.next();) {
        const int i = cursor.item();
        const CapsuleRecord& c = I.idx.caps[i];
        float s, t;
        const float d = segmentDistance(a, b, fromF(c.a), fromF(c.b), s, t);
        if (ex.covers(c.yarn, c.s0 + t * (c.s1 - c.s0))) continue;
        if (d < c.R + inflate) { hit = true; break; }
    }
    return hit;
}

struct Hit {
    bool inside = false;
    int yarn = -1;
    double s = 0;
    float depth = 0;
    f3 normal{0, 0, 0};  // the strand's outward normal in the ring plane at q (the ellipse's gradient)
};

// q inside capsule c's strand at its nearest ring: the hit's yarn, ring arc length, depth and normal.
__device__ __forceinline__ bool ringInside(const SurfaceIndexQuery& I, const CapsuleRecord& c, f3 q, float inflate, Hit& hit) {
    const YarnSurfaceRecord& ys = I.sm.yarns[c.yarn];
    const RingRecord* rings = I.sm.rings + ys.ringOff;
    float t;
    pointSegment(q, fromF(c.a), fromF(c.b), t);
    const double sEst = c.s0 + t * (c.s1 - c.s0);
    const int N = ys.nRings - 1;
    int k0 = static_cast<int>(llround(sEst / ys.ds));
    k0 = k0 < 0 ? 0 : (k0 > N ? N : k0);
    int best = -1;
    float bd = 3.0e38f;
    const int ka = k0 - 6 < 0 ? 0 : k0 - 6, kb = k0 + 6 > N ? N : k0 + 6;
    for (int k = ka; k <= kb; ++k) {
        const float d = fabsf(dot(q - fromF(rings[k].c), fromF(rings[k].T)));
        if (d < bd) { bd = d; best = k; }
    }
    if (best < 0 || bd > ys.ds) return false;
    const DeviceRing R = loadRing(rings[best]);
    f3 v = q - R.c;
    v = v - R.T * dot(v, R.T);
    f3 ea, eb;
    ringAxes(R, ea, eb);
    const float A = R.a + inflate, B = R.b + inflate;
    const float x = dot(v, ea) / A, y = dot(v, eb) / B;
    const float qq = x * x + y * y;
    if (qq >= 1.0f) return false;
    // The gradient of x^2 + y^2 in the ring plane; on the axis itself any direction is out: ea.
    const f3 grad = ea * (x / A) + eb * (y / B);
    hit.inside = true;
    hit.yarn = c.yarn;
    hit.s = R.s;
    hit.depth = (1 - sqrtf(qq)) * fminf_(A, B);
    hit.normal = norm(grad) > 1e-12f ? normalized(grad) : ea;
    return true;
}

__device__ __forceinline__ Hit exactInside(const SurfaceIndexQuery& I, f3 q, float inflate, const Exclude& only, const Exclude& ex) {
    Hit h;
    for (GridCursor cursor(I.idx.grid, q - f3{inflate, inflate, inflate}, q + f3{inflate, inflate, inflate}); cursor.next();) {
        const int i = cursor.item();
        const CapsuleRecord& c = I.idx.caps[i];
        if (only.yarn >= 0 && c.yarn != only.yarn) continue;
        float t;
        if (pointSegment(q, fromF(c.a), fromF(c.b), t) >= c.R + inflate) continue;
        Hit ring;
        if (!ringInside(I, c, q, inflate, ring)) continue;
        if (only.yarn >= 0 && !only.covers(c.yarn, ring.s)) continue;
        if (ex.covers(c.yarn, ring.s)) continue;
        if (!h.inside || ring.depth > h.depth) h = ring;
    }
    return h;
}

__device__ __forceinline__ float nearestAxisDistance(const SurfaceIndexQuery& I, f3 q, float radius, int yarn) {
    float best = 3.0e38f;  // CPU: 1e300 (none)
    for (GridCursor cursor(I.idx.grid, q - f3{radius, radius, radius}, q + f3{radius, radius, radius}); cursor.next();) {
        const int i = cursor.item();
        const CapsuleRecord& c = I.idx.caps[i];
        if (c.yarn != yarn) continue;
        float t;
        const float d = pointSegment(q, fromF(c.a), fromF(c.b), t);
        if (d < radius) best = fminf_(best, d);
    }
    return best;
}

// --- meshes
struct MeshHit {
    int tri = -1;
    float bary[3] = {1, 0, 0};
    float dist = 3.0e38f;
    f3 p{0, 0, 0};
};
// MeshQuery::closest (expanding search, the strict minimum; prunes of evaluated and of
// farther boxes, as mesh.cpp). maxDist: the CPU's 1e300 is 3e38.
__device__ inline MeshHit meshClosest(const MeshView& m, f3 q) {
    MeshHit best;
    const f3 lo = fromF(m.lo), hi = fromF(m.hi);
    const f3 d{fmaxf_(fmaxf_(lo.x - q.x, 0.0f), q.x - hi.x), fmaxf_(fmaxf_(lo.y - q.y, 0.0f), q.y - hi.y),
               fmaxf_(fmaxf_(lo.z - q.z, 0.0f), q.z - hi.z)};
    const float rMax = norm(d) + norm(hi - lo);
    const float margin = 1e-6f * (1 + fmaxf_(fmaxf_(fabsf(q.x), fabsf(q.y)), fabsf(q.z)) + norm(hi - lo));  // float: CPU 1e-9
    bool havePrev = false;
    f3 pLo{0, 0, 0}, pHi{0, 0, 0};
    for (float r = fmaxf_(m.grid.cell, norm(d));; r *= 2) {
        const float rr = fminf_(r, rMax);
        const f3 qLo = q - f3{rr, rr, rr}, qHi = q + f3{rr, rr, rr};
        for (GridCursor cursor(m.grid, qLo, qHi); cursor.next();) {
            const int f = cursor.item();
            const f3 a = fromF(m.grid.boxLo + 3 * f), b = fromF(m.grid.boxHi + 3 * f);
            if (havePrev && !(a.x > pHi.x || b.x < pLo.x || a.y > pHi.y || b.y < pLo.y || a.z > pHi.z || b.z < pLo.z)) continue;
            if (best.tri >= 0) {
                const float dx = fmaxf_(fmaxf_(a.x - q.x, 0.0f), q.x - b.x), dy = fmaxf_(fmaxf_(a.y - q.y, 0.0f), q.y - b.y),
                            dz = fmaxf_(fmaxf_(a.z - q.z, 0.0f), q.z - b.z);
                const float lim = best.dist + margin;
                if (dx * dx + dy * dy + dz * dz > lim * lim) continue;
            }
            const int* t = m.tri + 3 * f;
            float bary[3];
            const f3 p = closestOnTriangle(q, fromF(m.pos + 3 * t[0]), fromF(m.pos + 3 * t[1]), fromF(m.pos + 3 * t[2]), bary);
            const float dist = norm(q - p);
            if (dist < best.dist) {
                best.tri = f;
                best.dist = dist;
                best.p = p;
                best.bary[0] = bary[0]; best.bary[1] = bary[1]; best.bary[2] = bary[2];
            }
        }
        if (best.dist <= rr || rr >= rMax) break;
        havePrev = true;
        pLo = qLo;
        pHi = qHi;
    }
    return best;
}

// One triangle's closest point to q.
__device__ __forceinline__ MeshHit triangleHit(const MeshView& m, f3 q, int f) {
    MeshHit h;
    const int* t = m.tri + 3 * f;
    h.tri = f;
    h.p = closestOnTriangle(q, fromF(m.pos + 3 * t[0]), fromF(m.pos + 3 * t[1]), fromF(m.pos + 3 * t[2]), h.bary);
    h.dist = norm(q - h.p);
    return h;
}
// How far v lies outside [lo, hi]; 0 inside.
__device__ __forceinline__ float outsideInterval(float v, float lo, float hi) { return fmaxf_(fmaxf_(lo - v, 0.0f), v - hi); }

// meshClosest's result, found from a triangle near the answer. The fabric-back test asks along a fibre, and the
// previous sample's closest triangle bounds the distance: U = |q - startTri's closest point|. Only the grid cells
// that meet the ball of radius U + margin are visited, where meshClosest scans a cube that doubles from the distance
// to the mesh's box: about (4 H / cell)^2 cells for a point H above a flat fabric.
// The result is meshClosest's, bit for bit:
//  - every triangle whose box lies within best + margin of q is evaluated (the cells meeting the ball list them all);
//  - if a second triangle comes within margin of the best, the tie is left to meshClosest, whose visiting order
//    decides it;
//  - otherwise the best triangle T* is unique by more than margin, and meshClosest evaluates and keeps it too: it can
//    leave T* out only by stopping on a cube that excludes T*'s box, which needs its own best below d* + margin, and
//    a triangle that close is evaluated here and would have been a tie. One triangle and one q give one MeshHit.
__device__ inline MeshHit meshClosestFrom(const MeshView& m, f3 q, int startTri) {
    const f3 lo = fromF(m.lo), hi = fromF(m.hi);
    const float margin = 1e-6f * (1 + fmaxf_(fmaxf_(fabsf(q.x), fabsf(q.y)), fabsf(q.z)) + norm(hi - lo));  // as meshClosest
    MeshHit best = triangleHit(m, q, startTri);
    if (!(best.dist < 3.0e38f)) return meshClosest(m, q);  // a degenerate start triangle (NaN), never meshClosest's answer
    float second = 3.0e38f;  // the nearest other triangle
    const GridView& g = m.grid;
    // The ball, two margins wider: a closest point on a cell's boundary (cellOf rounds) still meets its cell.
    const float R = best.dist + 2 * margin;
    const float inv = 1.0f / g.cell;
    const long long nj = g.hi[1] - g.lo[1] + 1, nk = g.hi[2] - g.lo[2] + 1;
    long long k0 = static_cast<long long>(floorf((q.z - R) * inv)), k1 = static_cast<long long>(floorf((q.z + R) * inv));
    k0 = k0 > g.lo[2] ? k0 : g.lo[2];
    k1 = k1 < g.hi[2] ? k1 : g.hi[2];
    for (long long k = k0; k <= k1; ++k) {
        const float dz = outsideInterval(q.z, k * g.cell, (k + 1) * g.cell);
        if (dz > R) continue;
        const float ry = sqrtf(fmaxf_(R * R - dz * dz, 0.0f));
        long long j0 = static_cast<long long>(floorf((q.y - ry) * inv)), j1 = static_cast<long long>(floorf((q.y + ry) * inv));
        j0 = j0 > g.lo[1] ? j0 : g.lo[1];
        j1 = j1 < g.hi[1] ? j1 : g.hi[1];
        for (long long j = j0; j <= j1; ++j) {
            const float dy = outsideInterval(q.y, j * g.cell, (j + 1) * g.cell);
            if (dz * dz + dy * dy > R * R) continue;
            const float rx = sqrtf(fmaxf_(R * R - dz * dz - dy * dy, 0.0f));
            long long i0 = static_cast<long long>(floorf((q.x - rx) * inv)), i1 = static_cast<long long>(floorf((q.x + rx) * inv));
            i0 = i0 > g.lo[0] ? i0 : g.lo[0];
            i1 = i1 < g.hi[0] ? i1 : g.hi[0];
            for (long long i = i0; i <= i1; ++i) {
                const long long c = ((i - g.lo[0]) * nj + (j - g.lo[1])) * nk + (k - g.lo[2]);
                for (unsigned at = g.cellStart[c]; at < g.cellStart[c + 1]; ++at) {
                    const int f = static_cast<int>(g.entries[at] & ((1u << 29) - 1));
                    if (f == best.tri) continue;  // the start triangle, or the best met again in another cell
                    const f3 a = fromF(g.boxLo + 3 * f), b = fromF(g.boxHi + 3 * f);
                    const float bx = outsideInterval(q.x, a.x, b.x), by = outsideInterval(q.y, a.y, b.y), bz = outsideInterval(q.z, a.z, b.z);
                    const float lim = best.dist + margin;
                    if (bx * bx + by * by + bz * bz > lim * lim) continue;
                    const MeshHit h = triangleHit(m, q, f);
                    if (h.dist < best.dist) {
                        second = best.dist;
                        best = h;
                    } else if (h.dist < second) {
                        second = h.dist;
                    }
                }
            }
        }
    }
    if (second <= best.dist + margin) return meshClosest(m, q);
    return best;
}

// MeshQuery::inside: closest feature's angle-weighted pseudo-normal (face, edge or vertex).
__device__ inline bool meshInside(const MeshView& m, f3 q) {
    if (q.x < m.lo[0] || q.y < m.lo[1] || q.z < m.lo[2] || q.x > m.hi[0] || q.y > m.hi[1] || q.z > m.hi[2]) return false;
    const MeshHit h = meshClosest(m, q);
    if (h.tri < 0) return false;
    const int* t = m.tri + 3 * h.tri;
    int zeros = 0;
    for (int k = 0; k < 3; ++k) zeros += h.bary[k] == 0.0f;
    f3 n = fromF(m.faceN + 3 * h.tri);
    if (zeros == 2) {
        for (int k = 0; k < 3; ++k)
            if (h.bary[k] == 1.0f) n = fromF(m.pvertN + 3 * t[k]);
    } else if (zeros == 1) {
        for (int k = 0; k < 3; ++k)
            if (h.bary[k] == 0.0f) n = fromF(m.edgeN + 9 * h.tri + 3 * k);
    }
    return h.dist > 0 && dot(q - h.p, n) < 0;
}

// Chart::closestPoint / MeshSurface::project normal at a closest point (h.tri >= 0): the barycentric
// interpolation of the area-weighted vertex normals (MeshSurface: the face normal where it
// vanishes).
__device__ inline f3 interpolatedNormal(const MeshView& m, const MeshHit& h) {
    const int* t = m.tri + 3 * h.tri;
    const f3 n = fromF(m.vnormal + 3 * t[0]) * h.bary[0] + fromF(m.vnormal + 3 * t[1]) * h.bary[1] + fromF(m.vnormal + 3 * t[2]) * h.bary[2];
    if (norm(n) < kTinyLen) return fromF(m.faceN + 3 * h.tri);
    return normalized(n);
}
// The normal at the closest point of q; zero without a mesh.
__device__ inline f3 meshOutwardNormal(const MeshView& m, f3 q) {
    const MeshHit h = meshClosest(m, q);
    if (h.tri < 0) return f3{0, 0, 0};
    return interpolatedNormal(m, h);
}
// The signed height of q above the mesh, along the interpolated normal at the closest point, for samples
// along a polyline: the search starts at `tri`, the previous sample's closest triangle (-1: none yet), which
// is updated. meshClosestFrom finds the same closest point as meshClosest.
__device__ inline float meshHeightFrom(const MeshView& m, f3 q, int& tri) {
    const MeshHit h = tri >= 0 ? meshClosestFrom(m, q, tri) : meshClosest(m, q);
    if (h.tri < 0) return 0;
    tri = h.tri;
    return dot(q - h.p, interpolatedNormal(m, h));
}

// --- mesh continuation
// GPU backend: MeshSurface::continueFrom (src/knitting/layout/base_surface.cpp) on the device, single precision:
// steepest descent of |q - x| over the triangle adjacency from the previous accepted foot, with
// the walk budget |q - q_prev| + lambda w, the normal-continuity limit, the open-boundary and
// the depth-band rules. Results as knit::ProjResult (0 Owned, 1 Budget, 2 NormalTurn,
// 3 OpenBoundary, 4 DepthBand); the state is advanced only on Owned.
//
// Float constants (CPU in brackets): barycentric zero 1e-12 [1e-12] (closestOnTriangle returns
// exact zeros on edges / vertices, so the test is the same); required improvement of a step
// 1e-6 (1 + d) mm [1e-12 (1 + d)]: below float resolution of d the CPU's margin is 0 ulp, the
// float walk needs a few ulp to stop on coplanar fans.

struct ProjectionState {
    int tri;
    float b[3];
    f3 foot, n, query;
};

struct DeviceContinuationParams {
    float lambdaW, minNormalCos, band;
    int maxSteps;
};

// MeshSurface::normalAt: interpolated vertex normals, the face normal where they vanish.
__device__ __forceinline__ f3 meshNormalAt(const MeshView& m, int tri, const float b[3]) {
    const int* t = m.tri + 3 * tri;
    f3 n = fromF(m.vnormal + 3 * t[0]) * b[0] + fromF(m.vnormal + 3 * t[1]) * b[1] + fromF(m.vnormal + 3 * t[2]) * b[2];
    if (norm(n) < kTinyLen) n = fromF(m.faceN + 3 * tri);
    return normalized(n);
}

__device__ __forceinline__ f3 triClosest(const MeshView& m, int t, f3 q, float b[3]) {
    const int* v = m.tri + 3 * t;
    return closestOnTriangle(q, fromF(m.pos + 3 * v[0]), fromF(m.pos + 3 * v[1]), fromF(m.pos + 3 * v[2]), b);
}

__device__ inline int continueFrom(const MeshView& m, const AdjacencyView& adj, const ProjectionState& prev, f3 q, ProjectionState& next,
                            const DeviceContinuationParams& cp) {
    next = prev;
    if (prev.tri < 0) return 1;  // Budget
    const float budget = norm(q - prev.query) + cp.lambdaW;
    int t = prev.tri;
    float b[3];
    f3 x = triClosest(m, t, q, b);
    float d = norm(q - x), walked = norm(x - prev.foot);
    if (dot(meshNormalAt(m, t, b), prev.n) < cp.minNormalCos) return 2;
    if (walked > budget) return 1;
    const float zero = 1e-12f;
    int steps = 0;
    for (; steps < cp.maxSteps; ++steps) {
        int nz = 0, corner = -1;
        for (int k = 0; k < 3; ++k)
            if (b[k] <= zero) ++nz;
            else corner = corner < 0 ? k : corner;
        int best = -1;
        float bestD = d - 1e-6f * (1 + d), bb[3] = {0, 0, 0};
        f3 bx{0, 0, 0};
        auto tryTri = [&](int c) {
            float bc[3];
            const f3 xc = triClosest(m, c, q, bc);
            const float dc = norm(q - xc);
            if (dc < bestD) {
                best = c;
                bestD = dc;
                bx = xc;
                bb[0] = bc[0]; bb[1] = bc[1]; bb[2] = bc[2];
            }
        };
        if (nz == 1) {
            for (int k = 0; k < 3; ++k)
                if (b[k] <= zero) {
                    const int nb = adj.nbr[3 * t + (k + 1) % 3];
                    if (nb >= 0) tryTri(nb);
                }
        } else if (nz == 2) {
            const int v = m.tri[3 * t + corner];
            for (int at = adj.vtStart[v]; at < adj.vtStart[v + 1]; ++at)
                if (adj.vtList[at] != t) tryTri(adj.vtList[at]);
        }
        if (best < 0) break;
        walked += norm(bx - x);
        t = best;
        x = bx;
        d = bestD;
        b[0] = bb[0]; b[1] = bb[1]; b[2] = bb[2];
        if (dot(meshNormalAt(m, t, b), prev.n) < cp.minNormalCos) return 2;
        if (walked > budget) return 1;
    }
    if (steps >= cp.maxSteps) return 1;
    for (int k = 0; k < 3; ++k)
        if (b[k] <= zero && adj.nbr[3 * t + (k + 1) % 3] < 0) return 3;  // OpenBoundary
    const f3 n = meshNormalAt(m, t, b);
    if (fabsf(dot(q - x, n)) > cp.band) return 4;  // DepthBand
    next.tri = t;
    next.b[0] = b[0]; next.b[1] = b[1]; next.b[2] = b[2];
    const int* v = m.tri + 3 * t;
    next.foot = fromF(m.pos + 3 * v[0]) * b[0] + fromF(m.pos + 3 * v[1]) * b[1] + fromF(m.pos + 3 * v[2]) * b[2];
    next.n = n;
    next.query = q;
    return 0;  // Owned
}

// --- the same surface in double precision: the solver's band
enum BandMode : int { kBandNone = 0, kBandMesh, kBandBind, kBandCone, kBandPlane };

// The band's surface in double precision (the CPU's arithmetic on the CPU's mesh): a smooth,
// deterministic band value. The float device ports (surface.cuh, surface.cuh) carry
// ~1e-5 mm of noise at 100 mm coordinates, larger than the merit changes near convergence: a band
// row stayed violated by ~1e-5, its multiplier grew every outer iteration and the solver stalled
// (linear-accuracy experiment, analysis/linacc). Triangle candidates come from the resident float
// grid (box queries widened by a guard); every distance is double.
struct DoubleMesh {
    const double* pos = nullptr;     // 3 per vertex (TriMesh::pos)
    const double* vnormal = nullptr; // 3 per vertex (TriMesh::vnormal)
    const double* faceN = nullptr;   // 3 per triangle (TriMesh::faceNormal)
    const int* tri = nullptr;
    GridView grid;                   // triangle boxes (float), the resident mesh's
    AdjacencyView adj{};              // MeshSurface's neighbours
    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
};
struct DoubleProjectionState {  // knit::ProjState in double
    int tri;
    double b[3];
    d3 foot, n, query;
};
struct DoubleContinuationParams {
    double lambdaW, minNormalCos, band;
    int maxSteps;
};
struct BandParams {
    int mode = kBandNone;
    bool faceFallback = false;  // MeshSurface::normalAt (face normal where the vertex normals vanish); Chart: none
    DoubleMesh mesh;
    DoubleContinuationParams cont{};
    double R0 = 0, R1 = 0, H = 0;
};

// mesh.cpp closestOnTriangle (Ericson 5.1.5), double.
__device__ inline d3 closestOnTriangle(d3 p, d3 a, d3 b, d3 c, double bary[3]) {
    const d3 ab = b - a, ac = c - a, ap = p - a;
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) { bary[0] = 1; bary[1] = bary[2] = 0; return a; }
    const d3 bp = p - b;
    const double d3_ = dot(ab, bp), d4 = dot(ac, bp);
    if (d3_ >= 0 && d4 <= d3_) { bary[1] = 1; bary[0] = bary[2] = 0; return b; }
    const double vc = d1 * d4 - d3_ * d2;
    if (vc <= 0 && d1 >= 0 && d3_ <= 0) {
        const double v = d1 / (d1 - d3_);
        bary[0] = 1 - v; bary[1] = v; bary[2] = 0;
        return a + ab * v;
    }
    const d3 cp = p - c;
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) { bary[2] = 1; bary[0] = bary[1] = 0; return c; }
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) {
        const double w = d2 / (d2 - d6);
        bary[0] = 1 - w; bary[1] = 0; bary[2] = w;
        return a + ac * w;
    }
    const double va = d3_ * d6 - d5 * d4;
    if (va <= 0 && (d4 - d3_) >= 0 && (d5 - d6) >= 0) {
        const double w = (d4 - d3_) / ((d4 - d3_) + (d5 - d6));
        bary[0] = 0; bary[1] = 1 - w; bary[2] = w;
        return b + (c - b) * w;
    }
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom, w = vc * denom;
    bary[0] = 1 - v - w; bary[1] = v; bary[2] = w;
    return a + ab * v + ac * w;
}
__device__ __forceinline__ d3 dvec(const double* p, int i) { return {p[3 * i], p[3 * i + 1], p[3 * i + 2]}; }
__device__ __forceinline__ d3 triClosest(const DoubleMesh& m, int t, d3 q, double b[3]) {
    const int* v = m.tri + 3 * t;
    return closestOnTriangle(q, dvec(m.pos, v[0]), dvec(m.pos, v[1]), dvec(m.pos, v[2]), b);
}
// MeshSurface::normalAt (faceFallback) / Chart::at's normal (normalized interpolation).
__device__ inline d3 meshNormalAt(const DoubleMesh& m, int t, const double b[3], bool faceFallback) {
    const int* v = m.tri + 3 * t;
    d3 n = dvec(m.vnormal, v[0]) * b[0] + dvec(m.vnormal, v[1]) * b[1] + dvec(m.vnormal, v[2]) * b[2];
    const double l = norm(n);
    if (faceFallback && l < 1e-12) return dvec(m.faceN, t);
    return l > 0 ? n / l : n;
}

// Chart::closestPoint / MeshQuery-like global closest point: expanding box queries from the grid cell
// size, every triangle whose (float) box meets the query box (widened by a guard) evaluated in double;
// stops once the best distance lies inside the queried radius.
__device__ inline int meshClosest(const DoubleMesh& m, d3 q, double b[3], d3& foot) {
    int best = -1;
    double bestD = 1e300;
    const double ext = sqrt((m.hi[0] - m.lo[0]) * (m.hi[0] - m.lo[0]) + (m.hi[1] - m.lo[1]) * (m.hi[1] - m.lo[1]) +
                            (m.hi[2] - m.lo[2]) * (m.hi[2] - m.lo[2]));
    const double dOut = sqrt(fmax(fmax(m.lo[0] - q.x, 0.0), q.x - m.hi[0]) * fmax(fmax(m.lo[0] - q.x, 0.0), q.x - m.hi[0]) +
                             fmax(fmax(m.lo[1] - q.y, 0.0), q.y - m.hi[1]) * fmax(fmax(m.lo[1] - q.y, 0.0), q.y - m.hi[1]) +
                             fmax(fmax(m.lo[2] - q.z, 0.0), q.z - m.hi[2]) * fmax(fmax(m.lo[2] - q.z, 0.0), q.z - m.hi[2]));
    const double rMax = dOut + ext + 1;
    for (double r = fmax(double(m.grid.cell), dOut);; r *= 2) {
        const double rr = fmin(r, rMax);
        const double guard = 1e-4 * (1 + rr);  // the float boxes' rounding
        const f3 lo = tof(q - d3{rr + guard, rr + guard, rr + guard}), hi = tof(q + d3{rr + guard, rr + guard, rr + guard});
        for (GridCursor cursor(m.grid, lo, hi); cursor.next();) {
            const int t = cursor.item();
            double bc[3];
            const d3 x = triClosest(m, t, q, bc);
            const double d = norm(q - x);
            if (d < bestD || (d == bestD && t < best)) {
                bestD = d;
                best = t;
                foot = x;
                b[0] = bc[0]; b[1] = bc[1]; b[2] = bc[2];
            }
        }
        if ((best >= 0 && bestD <= rr) || rr >= rMax) break;
    }
    return best;
}

// MeshSurface::continueFrom (base_surface.cpp), double, the CPU's constants. 0 Owned, 1 Budget,
// 2 NormalTurn, 3 OpenBoundary, 4 DepthBand; next is written on Owned only.
__device__ inline int continueFrom(const DoubleMesh& m, const DoubleProjectionState& prev, d3 q, DoubleProjectionState& next, const DoubleContinuationParams& cp) {
    if (prev.tri < 0) return 1;
    const double budget = norm(q - prev.query) + cp.lambdaW;
    int t = prev.tri;
    double b[3];
    d3 x = triClosest(m, t, q, b);
    double d = norm(q - x), walked = norm(x - prev.foot);
    if (dot(meshNormalAt(m, t, b, true), prev.n) < cp.minNormalCos) return 2;
    if (walked > budget) return 1;
    const double zero = 1e-12;
    int steps = 0;
    for (; steps < cp.maxSteps; ++steps) {
        int nz = 0, corner = -1;
        for (int k = 0; k < 3; ++k)
            if (b[k] <= zero) ++nz;
            else corner = corner < 0 ? k : corner;
        int best = -1;
        double bestD = d - 1e-12 * (1 + d), bb[3] = {0, 0, 0};
        d3 bx{0, 0, 0};
        auto tryTri = [&](int c) {
            double bc[3];
            const d3 xc = triClosest(m, c, q, bc);
            const double dc = norm(q - xc);
            if (dc < bestD) {
                best = c;
                bestD = dc;
                bx = xc;
                bb[0] = bc[0]; bb[1] = bc[1]; bb[2] = bc[2];
            }
        };
        if (nz == 1) {
            for (int k = 0; k < 3; ++k)
                if (b[k] <= zero) {
                    const int nb = m.adj.nbr[3 * t + (k + 1) % 3];
                    if (nb >= 0) tryTri(nb);
                }
        } else if (nz == 2) {
            const int v = m.tri[3 * t + corner];
            for (int at = m.adj.vtStart[v]; at < m.adj.vtStart[v + 1]; ++at)
                if (m.adj.vtList[at] != t) tryTri(m.adj.vtList[at]);
        }
        if (best < 0) break;
        walked += norm(bx - x);
        t = best;
        x = bx;
        d = bestD;
        b[0] = bb[0]; b[1] = bb[1]; b[2] = bb[2];
        if (dot(meshNormalAt(m, t, b, true), prev.n) < cp.minNormalCos) return 2;
        if (walked > budget) return 1;
    }
    if (steps >= cp.maxSteps) return 1;
    for (int k = 0; k < 3; ++k)
        if (b[k] <= zero && m.adj.nbr[3 * t + (k + 1) % 3] < 0) return 3;
    const d3 n = meshNormalAt(m, t, b, true);
    if (fabs(dot(q - x, n)) > cp.band) return 4;
    next.tri = t;
    next.b[0] = b[0]; next.b[1] = b[1]; next.b[2] = b[2];
    const int* v = m.tri + 3 * t;
    next.foot = dvec(m.pos, v[0]) * b[0] + dvec(m.pos, v[1]) * b[1] + dvec(m.pos, v[2]) * b[2];
    next.n = n;
    next.query = q;
    return 0;
}

// The band's surface at x: foot and unit normal; false when the binding's continuation from the
// refresh state is not Owned (or no projection was found).
__device__ inline bool bandSurface(const BandParams& bp, const DoubleProjectionState* state, int g, d3 xd, d3& foot, d3& nrm) {
    if (bp.mode == kBandBind) {
        DoubleProjectionState next;
        if (continueFrom(bp.mesh, state[g], xd, next, bp.cont) != 0) return false;
        foot = next.foot;
        nrm = next.n;
        return true;
    }
    if (bp.mode == kBandMesh) {
        double b[3];
        const int t = meshClosest(bp.mesh, xd, b, foot);
        if (t < 0) return false;
        nrm = meshNormalAt(bp.mesh, t, b, bp.faceFallback);
        return true;
    }
    if (bp.mode == kBandCone) {
        const d3 r{xd.x, xd.y, 0};
        const double rho = norm(r);
        const d3 radial = rho > 0 ? r / rho : d3{1, 0, 0};
        const d3 A{bp.R0, 0, 0}, d{bp.R1 - bp.R0, 0, bp.H}, Q{rho, 0, xd.z};
        const double a = dot(Q - A, d) / dot(d, d);
        const double slope = (bp.R1 - bp.R0) / bp.H;
        const d3 nn = radial - d3{0, 0, slope};
        nrm = nn / norm(nn);
        foot = radial * (bp.R0 + (bp.R1 - bp.R0) * a) + d3{0, 0, a * bp.H};
        return true;
    }
    if (bp.mode == kBandPlane) {
        foot = {xd.x, xd.y, 0};
        nrm = {0, 0, 1};
        return true;
    }
    return false;
}

} // namespace knit::gpu
