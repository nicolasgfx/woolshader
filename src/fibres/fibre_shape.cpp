// One fibre as pure mathematics: draws, curves and the arc-length fit.
#include "fibres/fibre_shape.h"

#include <algorithm>
#include <cmath>

#include "knitting/tube/yarn_surface.h"
#include "model/mesh.h"
#include "model/hash.h"

namespace knit {

namespace {

constexpr double kDeg = kPi / 180;
constexpr int kSolveChords = 16;  // chords of the table while solving for the span
constexpr int kStraySideChannel = 30;  // hash channel: which way along the yarn a spread stray starts

} // namespace

double normalCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

double normalQuantile(double p) {
    static const double a[] = {-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
                               1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00};
    static const double b[] = {-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
                               6.680131188771972e+01, -1.328068155288572e+01};
    static const double c[] = {-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
                               -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00};
    static const double d[] = {7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
                               3.754408661907416e+00};
    p = std::clamp(p, 1e-300, 1 - 1e-16);
    double x;
    if (p < 0.02425) {
        const double q = std::sqrt(-2 * std::log(p));
        x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else if (p > 1 - 0.02425) {
        const double q = std::sqrt(-2 * std::log1p(-p));
        x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
    } else {
        const double q = p - 0.5, r = q * q;
        x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
            (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
    }
    for (int it = 0; it < 2; ++it) {
        const double e = normalCdf(x) - p, pdf = std::exp(-0.5 * x * x) / std::sqrt(2 * kPi);
        if (pdf > 0) x -= e / pdf;
    }
    return x;
}

double fibreTargetLength(const FibreParams& fp, uint64_t seed, int yarn, int i) {
    if (fp.lengthMean <= 0) return fp.length;
    const double u = hash01(seed, yarn, i, 6);
    if (fp.lengthSpread > 0) {
        // ln L = mu + sigma z, z restricted to the interval the bounds give.
        const double mu = std::log(fp.lengthMean), sg = fp.lengthSpread;
        const double Fa = normalCdf((std::log(fp.lengthMin) - mu) / sg), Fb = normalCdf((std::log(fp.lengthMax) - mu) / sg);
        const double L = std::exp(mu + sg * normalQuantile(Fa + u * (Fb - Fa)));
        return std::clamp(L, fp.lengthMin, fp.lengthMax);  // rounding only: the quantile is inside
    }
    return std::clamp(-fp.lengthMean * std::log1p(-u), fp.lengthMin, fp.lengthMax);
}

bool hasRandomOverrides(const FibreParams& fp) {
    return fp.randomLengthMean > 0 || fp.randomLengthMin > 0 || fp.randomLengthMax > 0 || fp.randomLengthSpread >= 0 ||
           fp.randomThickness > 0 || fp.randomThicknessSigma >= 0;
}

FibreParams populationParams(const FibreParams& fp, FibrePopulation p) {
    FibreParams v = fp;
    if (p != FibrePopulation::Random) return v;
    if (fp.randomLengthMean > 0) v.lengthMean = fp.randomLengthMean;
    if (fp.randomLengthMin > 0) v.lengthMin = fp.randomLengthMin;
    if (fp.randomLengthMax > 0) v.lengthMax = fp.randomLengthMax;
    if (fp.randomLengthSpread >= 0) v.lengthSpread = fp.randomLengthSpread;
    if (fp.randomThickness > 0) v.thickness = fp.randomThickness;
    if (fp.randomThicknessSigma >= 0) v.thicknessSigma = fp.randomThicknessSigma;
    return v;
}

FibrePopulation drawPopulation(const FibreParams& fp, uint64_t seed, int yarn, int i, const IslandWeight& island) {
    const double share = double(island.fuzzy) * kFuzzyStrays + island.matted;
    return hash01(seed, yarn, i, 17) < fp.randomFraction + (1 - fp.randomFraction) * share ? FibrePopulation::Random : FibrePopulation::Directed;
}

FibreDraw drawFibre(const FibreParams& fp, uint64_t seed, int yarn, int i, const IslandWeight& island) {
    FibreDraw d;
    d.sgn = hash01(seed, yarn, i, 7) < 0.5 ? 1 : -1;
    d.psi = fp.handedness * (fp.helixDeg + fp.helixSigmaDeg * hashNormal(seed, yarn, i, 8, 10, 2.5)) * kDeg;
    d.lift = (fp.liftMinDeg + (fp.liftMaxDeg - fp.liftMinDeg) * hash01(seed, yarn, i, 11)) * kDeg;
    d.bowTheta = (hash01(seed, yarn, i, 12) < 0.5 ? -1 : 1) * (0.5 + 0.5 * hash01(seed, yarn, i, 13)) * fp.bowTheta;
    d.bowRadial = (2 * hash01(seed, yarn, i, 14) - 1) * fp.bowRadial;
    d.population = drawPopulation(fp, seed, yarn, i, island);
    if (d.population == FibrePopulation::Directed && fp.crimpAmplitude > 0) {
        d.crimpPhase = 2 * kPi * hash01(seed, yarn, i, 28);
        d.crimpPlane = (2 * hash01(seed, yarn, i, 29) - 1) * 0.25 * kPi;
    }
    // A random fibre draws its length and thickness factor from its own parameters, on the same channels.
    const FibreParams lp = populationParams(fp, d.population);
    d.length = fibreTargetLength(lp, seed, yarn, i) * (1 + kFuzzyLength * island.fuzzy) * (1 - kMattedLength * island.matted);
    d.radiusFactor = std::exp(lp.thicknessSigma * hashNormal(seed, yarn, i, 15, 16, 2.0));
    d.key = hash01(seed, yarn, i, 9);
    // The start direction's azimuth around the surface normal, 0 = along the yarn (randomCurve).
    // At spread 180 every direction (the original formula, so its fibres stay bit-identical);
    // below, within +-spread of one of the two directions along the yarn.
    if (fp.randomSpreadDeg >= 180) {
        d.azimuth = 2 * kPi * hash01(seed, yarn, i, 18);
    } else {
        const double side = hash01(seed, yarn, i, kStraySideChannel) < 0.5 ? 0.0 : kPi;
        d.azimuth = side + fp.randomSpreadDeg * kDeg * (2 * hash01(seed, yarn, i, 18) - 1);
    }
    d.elev = (fp.randomElevMinDeg + (fp.randomElevMaxDeg - fp.randomElevMinDeg) * hash01(seed, yarn, i, 19)) * kDeg * (1 - kMattedFlat * island.matted);
    d.bendAngle = 2 * kPi * hash01(seed, yarn, i, 20);
    d.bendAmp = (hash01(seed, yarn, i, 21) < 0.5 ? -1 : 1) * (0.5 + 0.5 * hash01(seed, yarn, i, 22)) * fp.randomBend;
    if (d.population == FibrePopulation::Random && fp.randomCurlRadius > 0) {
        // Each fibre its own turn length and radius; the helix never takes more than 90 % of the
        // length (the axis must keep some), and a short fibre gets a bend, not a knot.
        const double turnLength = fp.randomCurlLength * (0.7 + 0.8 * hash01(seed, yarn, i, 23));
        d.curlTurns = std::min(d.length / turnLength, fp.maxSegments / double(kCurlChordsPerTurn));
        d.curlRadius = std::min({fp.randomCurlRadius * (0.5 + 0.5 * hash01(seed, yarn, i, 24)),
                                 0.9 * d.length / (2 * kPi * d.curlTurns), 0.15 * d.length});
        d.curlPhase = 2 * kPi * hash01(seed, yarn, i, 25);
        d.curlDrift = 2 * kPi * hash01(seed, yarn, i, 26);
        d.curlSense = hash01(seed, yarn, i, 27) < 0.5 ? -1 : 1;
    }
    return d;
}

YarnCoords DirectedCurve::coords(double u) const {
    const double slopeTheta = sgn * (S * p / r), slopeLambda = S * q / r;
    const double bow = std::sin(kPi * u), dBow = kPi * std::cos(kPi * u);
    YarnCoords y;
    y.s = sE + sgn * S * u;
    y.theta = thetaE + slopeTheta * u + bowTheta * bow;
    y.lambda = lambda0 + slopeLambda * u + bowLambda * bow;
    y.ds = sgn * S;
    y.dTheta = slopeTheta + bowTheta * dBow;
    y.dLambda = slopeLambda + bowLambda * dBow;
    if (hug > 0) {
        const double level = std::exp(-slopeLambda / hug * u);
        y.lambda = lambda0 + hug * (1 - level) + bowLambda * bow;
        y.dLambda = slopeLambda * level + bowLambda * dBow;
    }
    if (rampWinding && S > 0) {
        // The integral of (1 - e^-au)^2 = 1 - 2 e^-au + e^-2au from 0 to u.
        const double a = S / kWindRampLength, e = std::exp(-a * u);
        const double wound = u - 2 * (1 - e) / a + (1 - e * e) / (2 * a);
        y.theta = thetaE + slopeTheta * wound + bowTheta * bow;
        y.dTheta = slopeTheta * (1 - e) * (1 - e) + bowTheta * dBow;
    }
    if (crimpAmplitude > 0 && S > 0) {
        const double cycles = S / crimpWavelength, ramp = std::clamp(0.5 / cycles, 0.05, 0.5);
        const double t = std::min(1.0, u / ramp), f = t * t * (3 - 2 * t), dF = u < ramp ? 6 * t * (1 - t) / ramp : 0;
        const double w = 2 * kPi * cycles, wave = std::sin(w * u + crimpPhase), dWave = w * std::cos(w * u + crimpPhase);
        const double offset = crimpAmplitude * f * wave / r, dOffset = crimpAmplitude * (dF * wave + f * dWave) / r;
        y.theta += offset * std::cos(crimpPlane);
        y.lambda += offset * std::sin(crimpPlane);
        y.dTheta += dOffset * std::cos(crimpPlane);
        y.dLambda += dOffset * std::sin(crimpPlane);
    }
    return y;
}

Vec3 DirectedCurve::X(double u) const {
    const YarnCoords y = coords(u);
    return yarnPosition(*cs, *ys, y.s, y.theta, y.lambda);
}

Vec3 DirectedCurve::dX(double u) const {
    const YarnCoords y = coords(u);
    const YarnCoordPoint point = yarnCoordinates(*cs, *ys, y.s, y.theta, y.lambda);
    return point.dS * y.ds + point.dTheta * y.dTheta + point.dLambda * y.dLambda;
}

RandomCurve randomCurve(const FibreDraw& d, const Vec3& A, const Vec3& n, const Vec3& T) {
    Vec3 t1 = T - n * dot(T, n);
    t1 = norm(t1) > 1e-12 ? normalized(t1) : normalized(cross(n, std::fabs(n.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}));
    const Vec3 t2 = cross(n, t1);
    RandomCurve c;
    c.A = A;
    c.d0 = normalized((t1 * std::cos(d.azimuth) + t2 * std::sin(d.azimuth)) * std::cos(d.elev) + n * std::sin(d.elev));
    // Basis perpendicular to d0: e1 towards n (d0 is never parallel to n below 90 deg elevation).
    Vec3 e1 = n - c.d0 * dot(n, c.d0);
    e1 = norm(e1) > 1e-9 ? normalized(e1) : t1;
    const Vec3 e2 = cross(c.d0, e1);
    c.w = e1 * std::cos(d.bendAngle) + e2 * std::sin(d.bendAngle);
    c.a = d.bendAmp;
    c.e1 = e1;
    c.e2 = e2;
    c.rho = d.curlRadius;
    if (c.rho > 0) {
        c.turns = d.curlTurns;
        c.phase = d.curlPhase;
        c.drift = d.curlDrift;
        c.sense = d.curlSense;
        c.rampEnd = std::clamp(0.5 / c.turns, 0.05, 0.5);
    }
    return c;
}

Vec3 RandomCurve::X(double u) const {
    const Vec3 axis = A + (d0 * u + w * (a * u * u)) * S;
    return rho > 0 ? axis + curl(u) : axis;
}

Vec3 RandomCurve::dX(double u) const {
    const Vec3 axis = (d0 + w * (2 * a * u)) * S;
    return rho > 0 ? axis + dCurl(u) : axis;
}

Vec3 RandomCurve::curl(double u) const {
    const double t = 2 * kPi * turns * u, r = std::min(1.0, u / rampEnd);
    const double radius = rho * r * r * (3 - 2 * r) * (1 + 0.3 * std::cos(0.37 * t + drift));
    const double th = phase + sense * t + 0.7 * std::sin(0.6 * t + drift);
    return (e1 * std::cos(th) + e2 * std::sin(th)) * radius;
}

Vec3 RandomCurve::dCurl(double u) const {
    const double dt = 2 * kPi * turns, t = dt * u, r = std::min(1.0, u / rampEnd);
    const double ramp = r * r * (3 - 2 * r), dRamp = u < rampEnd ? 6 * r * (1 - r) / rampEnd : 0;
    const double wave = 1 + 0.3 * std::cos(0.37 * t + drift), dWave = -0.3 * 0.37 * dt * std::sin(0.37 * t + drift);
    const double radius = rho * ramp * wave, dRadius = rho * (dRamp * wave + ramp * dWave);
    const double th = phase + sense * t + 0.7 * std::sin(0.6 * t + drift), dTh = sense * dt + 0.7 * 0.6 * dt * std::cos(0.6 * t + drift);
    const Vec3 radial = e1 * std::cos(th) + e2 * std::sin(th), around = e2 * std::cos(th) - e1 * std::sin(th);
    return radial * dRadius + around * (radius * dTh);
}

double SmoothCurve::setSpan(double S, int chords) {
    if (random) {
        stray.S = S;
    } else {
        directed.S = S;
        const double cap = 0.9 * S * directed.q / (kPi * directed.r);  // keeps d lambda / du > 0
        directed.bowLambda = directed.hug > 0 ? 0 : std::copysign(std::min(std::fabs(bowRadial / directed.r), cap), bowRadial);
    }
    const int step = kTableChords / chords;
    arc[0] = 0;
    for (int k = 0; k <= kTableChords; k += step) {
        table[k] = point(double(k) / kTableChords);
        if (k) arc[k] = arc[k - step] + norm(table[k] - table[k - step]);
    }
    return arc[kTableChords];
}

bool SmoothCurve::solveLength(double target) {
    if (maxSpan <= 0) return false;
    double S0 = std::min(lastSpan > 0 ? lastSpan : target, maxSpan), l0 = setSpan(S0, kSolveChords) - target;
    if (l0 < 0 && S0 >= maxSpan) return false;
    double S1 = std::min(maxSpan, S0 * (target / std::max(1e-9, l0 + target))), l1 = setSpan(S1, kSolveChords) - target;
    for (int it = 0; it < 12 && std::fabs(l1) > 1e-6 * target; ++it) {
        double S2 = l1 != l0 ? S1 - l1 * (S1 - S0) / (l1 - l0) : S1 * target / (l1 + target);
        S2 = std::clamp(S2, 1e-6, maxSpan);
        S0 = S1;
        l0 = l1;
        S1 = S2;
        l1 = setSpan(S1, kSolveChords) - target;
        if (l1 < 0 && S1 >= maxSpan) return false;
    }
    if (std::fabs(l1) > 1e-4 * target) return false;
    lastSpan = S1;
    setSpan(S1, kTableChords);
    return true;
}

std::vector<Vec3> SmoothCurve::equalArcVertices(int nv) const {
    std::vector<Vec3> vertices;
    for (int k = 1; k <= nv; ++k) {
        const double target = length() * k / nv;
        int j = 1;
        while (j < kTableChords && arc[j] < target) ++j;
        const double t = (target - arc[j - 1]) / std::max(1e-15, arc[j] - arc[j - 1]);
        const double u = (j - 1 + std::clamp(t, 0.0, 1.0)) / kTableChords;
        vertices.push_back(k == nv ? table[kTableChords] : point(u));
    }
    return vertices;
}

double SmoothCurve::deviationFrom(const std::vector<Vec3>& polyline) const {
    double deviation = 0;
    for (int k = 0; k <= kTableChords; ++k) {
        double nearest = 1e300;
        for (size_t j = 1; j < polyline.size(); ++j)
            nearest = std::min(nearest, closestSegmentPoints(table[k], table[k], polyline[j - 1], polyline[j]).distance);
        deviation = std::max(deviation, nearest);
    }
    return deviation;
}

} // namespace knit
