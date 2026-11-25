// The .pptstrands writer and reader, and the renderer's rules every written strand keeps (strands_file.h has the layout).
#include "fibres/strands_file.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>

#include "model/vec.h"

namespace knit {

namespace {

constexpr char kMagic[8] = {'P', 'P', 'T', 'S', 'T', 'R', 'N', 'D'};

struct Header {
    char magic[8];
    uint32_t version;
    uint32_t basis;
    uint64_t vertexCount;
    uint64_t strandCount;
    uint64_t reserved;
};
static_assert(sizeof(Header) == 40, "the .pptstrands header is 40 bytes");

bool fail(std::string* err, const std::string& msg) {
    if (err) *err = msg;
    return false;
}

struct P3 { double x, y, z; };
P3 at(const StrandSetF& s, size_t v) { return {s.xyz[3 * v], s.xyz[3 * v + 1], s.xyz[3 * v + 2]}; }
double dist(const P3& a, const P3& b) { return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z)); }

// Circumradius of the triangle (a, b, c); +inf for collinear points.
double circumradius(const P3& a, const P3& b, const P3& c) {
    const double ab = dist(a, b), bc = dist(b, c), ca = dist(c, a);
    const P3 u{b.x - a.x, b.y - a.y, b.z - a.z}, w{c.x - a.x, c.y - a.y, c.z - a.z};
    const P3 x{u.y * w.z - u.z * w.y, u.z * w.x - u.x * w.z, u.x * w.y - u.y * w.x};
    const double area2 = std::sqrt(x.x * x.x + x.y * x.y + x.z * x.z);
    return area2 > 0 ? ab * bc * ca / (2 * area2) : 1e300;
}

std::string at(const char* what, size_t k) { return std::string(what) + " " + std::to_string(k); }

} // namespace

void StrandSetF::add(const std::vector<float>& pts, const std::vector<float>& radii) {
    if (first.empty()) first.push_back(0);
    xyz.insert(xyz.end(), pts.begin(), pts.end());
    radius.insert(radius.end(), radii.begin(), radii.end());
    first.push_back(static_cast<uint32_t>(radius.size()));
}

double taperedRadius(double r0, double r1, int k, int segments) {
    const double n = double(segments);
    const double tau = std::max(r1 / r0, std::pow(kTaperStep, -n));
    return r0 * std::pow(tau, double(k) / n);
}

const char* strandFaultName(StrandFault f) {
    switch (f) {
    case StrandFault::None: return "none";
    case StrandFault::Structure: return "structure";
    case StrandFault::NonFinite: return "non_finite";
    case StrandFault::RadiusFloor: return "radius_floor";
    case StrandFault::SegmentShort: return "segment_short";
    case StrandFault::SegmentLong: return "segment_long";
    case StrandFault::RadiusStep: return "radius_step";
    case StrandFault::Curvature: return "curvature";
    }
    return "unknown";
}

StrandCheck checkStrands(const StrandSetF& s, const StrandRules& rules) {
    StrandCheck c;
    auto bad = [&](StrandFault fault, const std::string& msg) {
        if (c.ok) {
            c.error = msg;
            c.fault = fault;
        }
        c.ok = false;
    };
    const size_t V = s.vertices();
    if (s.xyz.size() != 3 * V) bad(StrandFault::Structure, "xyz holds " + std::to_string(s.xyz.size()) + " floats for " + std::to_string(V) + " radii");
    if (s.first.size() < 2) bad(StrandFault::Structure, "no strands");
    if (!c.ok) return c;
    if (V > (1ull << 31)) bad(StrandFault::Structure, "more than 2^31 vertices");
    if (s.first.front() != 0 || s.first.back() != V) bad(StrandFault::Structure, "strand table does not run from 0 to V");
    for (size_t v = 0; v < V && c.ok; ++v) {
        const P3 p = at(s, v);
        const double r = s.radius[v];
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z)) bad(StrandFault::NonFinite, at("non-finite position at vertex", v));
        if (!std::isfinite(r) || !(r > 0)) bad(StrandFault::NonFinite, at("non-finite or non-positive radius at vertex", v));
        if (r < rules.minRadius) bad(StrandFault::RadiusFloor, at("radius below the minimum at vertex", v));
        c.minRadius = std::min(c.minRadius, r);
        c.maxRadius = std::max(c.maxRadius, r);
        const double m = std::max({std::fabs(p.x), std::fabs(p.y), std::fabs(p.z)});
        c.maxCoord = std::max(c.maxCoord, m);
        const double prec = r > 0 ? rules.precisionK * m / r : 1e300;
        c.worstPrecision = std::max(c.worstPrecision, prec);
        c.precisionWarnings += prec > 1;
    }
    for (size_t k = 0; k + 1 < s.first.size() && c.ok; ++k) {
        const uint32_t a = s.first[k], b = s.first[k + 1];
        if (b < a) { bad(StrandFault::Structure, at("strand table decreases at strand", k)); break; }
        if (b - a < 2) { bad(StrandFault::Structure, at("fewer than two vertices in strand", k)); break; }
        for (uint32_t v = a; v + 1 < b && c.ok; ++v) {
            const double ra = s.radius[v], rb = s.radius[v + 1], rl = std::min(ra, rb), rh = std::max(ra, rb);
            const double len = dist(at(s, v), at(s, v + 1));
            c.minSegOverR = std::min(c.minSegOverR, len / rl);
            c.maxSegOverR = std::max(c.maxSegOverR, len / rl);
            c.maxRadiusRatio = std::max(c.maxRadiusRatio, rh / rl);
            if (len < rules.minSegOverR * rl) bad(StrandFault::SegmentShort, at("segment shorter than minSegOverR * r in strand", k));
            if (len < rules.minSegOverRMax * rh) bad(StrandFault::SegmentShort, at("segment shorter than minSegOverRMax * the larger r in strand", k));
            if (len > rules.maxSegOverR * rl) bad(StrandFault::SegmentLong, at("segment longer than maxSegOverR * r in strand", k));
            if (rh > (1 + rules.maxRadiusStep) * rl) bad(StrandFault::RadiusStep, at("radius step above maxRadiusStep in strand", k));
            if (v + 2 < b) {
                const double rt = std::min(rl, double(s.radius[v + 2]));
                const double R = circumradius(at(s, v), at(s, v + 1), at(s, v + 2));
                c.minCurvOverR = std::min(c.minCurvOverR, R / rt);
                if (R < rt) bad(StrandFault::Curvature, at("discrete curvature radius below the fibre radius in strand", k));
            }
        }
    }
    return c;
}

bool writeStrandsFile(const std::string& path, const StrandSetF& s, const StrandRules& rules, std::string* err,
                      StrandCheck* check) {
    StrandCheck c = checkStrands(s, rules);
    if (check) *check = c;
    if (!c.ok) return fail(err, "refusing to write " + path + ": " + c.error);
    if (s.strands() == 0) return fail(err, "refusing to write an empty strand set to " + path);
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(path).parent_path(), ec);
    std::ofstream f(path, std::ios::binary);
    if (!f) return fail(err, "cannot open " + path);
    Header h{};
    std::memcpy(h.magic, kMagic, 8);
    h.version = 1;
    h.basis = 0;
    h.vertexCount = s.vertices();
    h.strandCount = s.strands();
    h.reserved = 0;
    f.write(reinterpret_cast<const char*>(&h), sizeof h);
    f.write(reinterpret_cast<const char*>(s.xyz.data()), static_cast<std::streamsize>(s.xyz.size() * sizeof(float)));
    f.write(reinterpret_cast<const char*>(s.radius.data()), static_cast<std::streamsize>(s.radius.size() * sizeof(float)));
    f.write(reinterpret_cast<const char*>(s.first.data()), static_cast<std::streamsize>(s.first.size() * sizeof(uint32_t)));
    if (!f) return fail(err, "short write to " + path);
    return true;
}

bool readStrandsFile(const std::string& path, StrandSetF& out, std::string* err) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return fail(err, "cannot open " + path);
    Header h{};
    if (!f.read(reinterpret_cast<char*>(&h), sizeof h)) return fail(err, "truncated header");
    if (std::memcmp(h.magic, kMagic, 8) != 0) return fail(err, "bad magic");
    if (h.version != 1) return fail(err, "version " + std::to_string(h.version));
    if (h.basis != 0) return fail(err, "basis " + std::to_string(h.basis));
    if (h.vertexCount > (1ull << 31) || h.strandCount >= h.vertexCount + 1 || h.reserved != 0)
        return fail(err, "implausible counts");
    std::error_code ec;
    const uint64_t bytes = std::filesystem::file_size(path, ec);
    const uint64_t expect = 40 + 16 * h.vertexCount + 4 * (h.strandCount + 1);
    if (ec || bytes != expect)
        return fail(err, std::to_string(bytes) + " bytes, counts declare " + std::to_string(expect));
    StrandSetF s;
    s.xyz.resize(3 * h.vertexCount);
    s.radius.resize(h.vertexCount);
    s.first.resize(h.strandCount + 1);
    f.read(reinterpret_cast<char*>(s.xyz.data()), static_cast<std::streamsize>(s.xyz.size() * sizeof(float)));
    f.read(reinterpret_cast<char*>(s.radius.data()), static_cast<std::streamsize>(s.radius.size() * sizeof(float)));
    f.read(reinterpret_cast<char*>(s.first.data()), static_cast<std::streamsize>(s.first.size() * sizeof(uint32_t)));
    if (!f) return fail(err, "truncated body");
    char extra;
    if (f.read(&extra, 1)) return fail(err, "trailing bytes");
    if (s.first.size() < 2 || s.first.front() != 0 || s.first.back() != h.vertexCount)
        return fail(err, "strand table does not run from 0 to V");
    for (size_t k = 0; k + 1 < s.first.size(); ++k)
        if (s.first[k + 1] < s.first[k] || s.first[k + 1] - s.first[k] < 2) return fail(err, at("bad strand", k));
    for (size_t v = 0; v < s.vertices(); ++v) {
        if (!std::isfinite(s.xyz[3 * v]) || !std::isfinite(s.xyz[3 * v + 1]) || !std::isfinite(s.xyz[3 * v + 2]))
            return fail(err, at("non-finite vertex", v));
        if (!std::isfinite(s.radius[v]) || !(s.radius[v] >= 1e-6f)) return fail(err, at("bad radius at vertex", v));
    }
    out = std::move(s);
    return true;
}

StrandSetF keepStrands(const StrandSetF& s, const std::vector<float>& keys, float q) {
    StrandSetF out;
    std::vector<float> pts, radii;
    for (size_t i = 0; i < s.strands() && i < keys.size(); ++i) {
        if (!(keys[i] < q)) continue;
        const uint32_t a = s.first[i], b = s.first[i + 1];
        pts.assign(s.xyz.begin() + 3 * a, s.xyz.begin() + 3 * b);
        radii.assign(s.radius.begin() + a, s.radius.begin() + b);
        out.add(pts, radii);
    }
    return out;
}

} // namespace knit
