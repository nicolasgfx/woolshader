#pragma once
// The cross-section fit of the contact flattening (yarn_surface.h buildSurface). The CPU works in double
// and the GPU's surface_build.cu in float, so every type and function exists twice with the same
// expressions: plain twins, resolved by their argument types.
//
// Every ring collects requirements (h, psi): its support along the direction psi (in the ring's
// (n1, n2) frame) may be at most h.
//   * contacts: a neighbouring centreline at distance d gives h = d / 2 (both sides take their
//     half, so the supports meet along the line between the centres);
//   * fold: the principal normal gives h = rho / (1 + kFoldMargin), so the swept cross-section
//     stays locally injective (yarn_surface.h checkSurface's fold criterion).
// A set keeps up to kMaxReq directions (directions closer than kMergeAngle, as axis lines, merge
// into the smaller h; a full set replaces its largest h). All h are below r: a larger one binds
// nothing. fit then gives the ellipse in closed form:
//   1. primary p = the smallest h: b = max(h_p, lower), the a-axis perpendicular to psi_p;
//   2. every other requirement at theta = psi_i - psi_p bounds
//        a <= sqrt((h_i^2 - b^2 cos^2 theta) / sin^2 theta)        (support^2 = a^2 sin^2 + b^2 cos^2);
//   3. a = max(min(r, bounds), lower) (a bound below `lower` is a clamp).
// Geometric approximation, not mechanics (yarn_surface.h).
#include <cmath>

#if defined(__CUDACC__)
#define KNIT_CFIT_HD __host__ __device__ __forceinline__
#else
#define KNIT_CFIT_HD inline
#endif

namespace knit::cfit {

constexpr int kMaxReq = 8;
constexpr double kMergeAngle = 0.17453292519943295;  // 10 degrees
constexpr double kFoldMargin = 0.02;                 // fold requirement h = rho / (1 + kFoldMargin)
constexpr double kPiD = 3.14159265358979323846;

struct Requirement {
    double h, psi;
};
struct RequirementF {
    float h, psi;
};

struct RequirementSet {
    int n = 0;
    Requirement e[kMaxReq];
};
struct RequirementSetF {
    int n = 0;
    RequirementF e[kMaxReq];
};

struct Section {
    double a, b, phi;
    bool flattened = false, clamped = false;
};
struct SectionF {
    float a, b, phi;
    bool flattened = false, clamped = false;
};

// |psi1 - psi2| between two axis lines, in [0, pi/2].
KNIT_CFIT_HD double axisDistance(double psi1, double psi2) {
    const double pi = kPiD;
    double x = psi1 - psi2;
    x -= pi * floor(x / pi);  // [0, pi)
    return x > pi / 2 ? pi - x : x;
}
KNIT_CFIT_HD float axisDistance(float psi1, float psi2) {
    const float pi = static_cast<float>(kPiD);
    float x = psi1 - psi2;
    x -= pi * floorf(x / pi);  // [0, pi)
    return x > pi / 2 ? pi - x : x;
}

// Keeps the requirement in the set. A full set gives up its largest h for a smaller new one, or
// else the new one: `dropped` receives the h given up, and stays as it was otherwise.
KNIT_CFIT_HD void insert(RequirementSet& S, double h, double psi, double& dropped) {
    for (int i = 0; i < S.n; ++i)
        if (axisDistance(S.e[i].psi, psi) < kMergeAngle) {
            if (h < S.e[i].h) S.e[i] = {h, psi};
            return;
        }
    if (S.n < kMaxReq) {
        S.e[S.n] = {h, psi};
        ++S.n;
        return;
    }
    int w = 0;
    for (int i = 1; i < S.n; ++i)
        if (S.e[i].h > S.e[w].h) w = i;
    if (h < S.e[w].h) {
        dropped = S.e[w].h;
        S.e[w] = {h, psi};
    } else {
        dropped = h;
    }
}
// The device does not count what a full set gives up.
KNIT_CFIT_HD void insert(RequirementSetF& S, float h, float psi) {
    for (int i = 0; i < S.n; ++i)
        if (axisDistance(S.e[i].psi, psi) < static_cast<float>(kMergeAngle)) {
            if (h < S.e[i].h) S.e[i] = {h, psi};
            return;
        }
    if (S.n < kMaxReq) {
        S.e[S.n] = {h, psi};
        ++S.n;
        return;
    }
    int w = 0;
    for (int i = 1; i < S.n; ++i)
        if (S.e[i].h > S.e[w].h) w = i;
    if (h < S.e[w].h) S.e[w] = {h, psi};
}

// The smoothing max-filter: a neighbouring ring's requirement at weight w (1 at the ring itself)
// relaxes towards r.
KNIT_CFIT_HD double smoothed(double h, double w, double r) { return r - (r - h) * w; }
KNIT_CFIT_HD float smoothed(float h, float w, float r) { return r - (r - h) * w; }

KNIT_CFIT_HD Section fit(const RequirementSet& S, double r, double lower) {
    Section o{r, r, 0.0};
    if (S.n == 0) return o;
    int p = 0;
    for (int i = 1; i < S.n; ++i)
        if (S.e[i].h < S.e[p].h) p = i;
    const double hp = S.e[p].h;
    if (!(hp < r)) return o;
    o.flattened = true;
    o.clamped = hp < lower;
    o.b = hp < lower ? lower : hp;  // std::max(hp, lower)
    o.phi = S.e[p].psi + 0.5 * kPiD;
    double a = r;
    for (int i = 0; i < S.n; ++i) {
        if (i == p) continue;
        const double c = cos(S.e[i].psi - S.e[p].psi), s2 = 1 - c * c;
        if (!(s2 > 1e-6)) continue;  // the same axis: b <= h_p <= h_i
        const double num = S.e[i].h * S.e[i].h - o.b * o.b * c * c;
        const double bound = num > 0 ? sqrt(num / s2) : 0.0;
        if (bound < a) a = bound;
    }
    if (a < lower) {
        a = lower;
        o.clamped = true;
    }
    o.a = a;
    return o;
}
KNIT_CFIT_HD SectionF fit(const RequirementSetF& S, float r, float lower) {
    SectionF o{r, r, 0.0f};
    if (S.n == 0) return o;
    int p = 0;
    for (int i = 1; i < S.n; ++i)
        if (S.e[i].h < S.e[p].h) p = i;
    const float hp = S.e[p].h;
    if (!(hp < r)) return o;
    o.flattened = true;
    o.clamped = hp < lower;
    o.b = hp < lower ? lower : hp;  // std::max(hp, lower)
    o.phi = S.e[p].psi + 0.5f * static_cast<float>(kPiD);
    float a = r;
    for (int i = 0; i < S.n; ++i) {
        if (i == p) continue;
        const float c = cosf(S.e[i].psi - S.e[p].psi), s2 = 1 - c * c;
        if (!(s2 > static_cast<float>(1e-6))) continue;  // the same axis: b <= h_p <= h_i
        const float num = S.e[i].h * S.e[i].h - o.b * o.b * c * c;
        const float bound = num > 0 ? sqrtf(num / s2) : 0.0f;
        if (bound < a) a = bound;
    }
    if (a < lower) {
        a = lower;
        o.clamped = true;
    }
    o.a = a;
    return o;
}

} // namespace knit::cfit
