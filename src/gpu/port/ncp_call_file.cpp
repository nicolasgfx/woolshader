// The relax call file (ncp_call_file.h): a flat little-endian record written field by field in the
// order below, each vector as its length followed by its elements, then the FNV-1a 64 digest of all
// preceding bytes. The two parameter structs are plain data written whole, guarded by their sizes.
#include "gpu/port/ncp_call_file.h"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <type_traits>

#include "knitting/layout/base_surface.h"

namespace knit::gpu {

namespace {

constexpr char kMagic[8] = {'N', 'C', 'P', 'C', 'A', 'L', 'L', '1'};
constexpr uint64_t kFnvOffset = 1469598103934665603ull, kFnvPrime = 1099511628211ull;
static_assert(std::is_trivially_copyable_v<KnitProblemParams> && std::is_trivially_copyable_v<NcpSolveParams>, "written as plain bytes");

// Appends to a binary file and digests what it wrote.
class CallWriter {
public:
    explicit CallWriter(const std::string& path) : f_(path, std::ios::binary) {}
    bool ok() const { return static_cast<bool>(f_); }
    void bytes(const void* p, size_t n) {
        const unsigned char* c = static_cast<const unsigned char*>(p);
        for (size_t i = 0; i < n; ++i) digest_ = (digest_ ^ c[i]) * kFnvPrime;
        f_.write(static_cast<const char*>(p), static_cast<std::streamsize>(n));
    }
    void count(size_t n) {
        const uint64_t v = n;
        bytes(&v, sizeof v);
    }
    void integer(long long v) { bytes(&v, sizeof v); }
    void number(double v) { bytes(&v, sizeof v); }
    void doubles(const std::vector<double>& v) {
        count(v.size());
        bytes(v.data(), v.size() * sizeof(double));
    }
    void points(const std::vector<Vec3>& v) {
        count(v.size());
        for (const Vec3& p : v) {
            number(p.x);
            number(p.y);
            number(p.z);
        }
    }
    void flags(const std::vector<char>& v) {
        count(v.size());
        bytes(v.data(), v.size());
    }
    void finish() {
        const uint64_t d = digest_;
        f_.write(reinterpret_cast<const char*>(&d), sizeof d);
        f_.flush();
    }

private:
    std::ofstream f_;
    uint64_t digest_ = kFnvOffset;
};

// Reads what CallWriter wrote; any short read or absurd length marks the file bad.
class CallReader {
public:
    explicit CallReader(const std::string& path) : f_(path, std::ios::binary) { good_ = static_cast<bool>(f_); }
    bool good() const { return good_; }
    void bytes(void* p, size_t n) {
        if (!good_) return;
        f_.read(static_cast<char*>(p), static_cast<std::streamsize>(n));
        if (static_cast<size_t>(f_.gcount()) != n) {
            good_ = false;
            return;
        }
        const unsigned char* c = static_cast<const unsigned char*>(p);
        for (size_t i = 0; i < n; ++i) digest_ = (digest_ ^ c[i]) * kFnvPrime;
    }
    size_t count() {
        uint64_t v = 0;
        bytes(&v, sizeof v);
        if (v > (uint64_t(1) << 32)) good_ = false;  // no field of a relax call is this long
        return good_ ? static_cast<size_t>(v) : 0;
    }
    long long integer() {
        long long v = 0;
        bytes(&v, sizeof v);
        return v;
    }
    double number() {
        double v = 0;
        bytes(&v, sizeof v);
        return v;
    }
    void doubles(std::vector<double>& v) {
        v.resize(count());
        bytes(v.data(), v.size() * sizeof(double));
    }
    void points(std::vector<Vec3>& v) {
        v.resize(count());
        for (Vec3& p : v) {
            p.x = number();
            p.y = number();
            p.z = number();
        }
    }
    void flags(std::vector<char>& v) {
        v.resize(count());
        bytes(v.data(), v.size());
    }
    // The stored digest equals the digest of everything read before it.
    bool digestMatches() {
        if (!good_) return false;
        const uint64_t expected = digest_;
        uint64_t stored = 0;
        f_.read(reinterpret_cast<char*>(&stored), sizeof stored);
        return static_cast<size_t>(f_.gcount()) == sizeof stored && stored == expected;
    }

private:
    std::ifstream f_;
    bool good_ = false;
    uint64_t digest_ = kFnvOffset;
};

void writeControls(CallWriter& w, const Controls& c) {
    w.count(c.size());
    for (const std::vector<Vec3>& y : c) w.points(y);
}

void readControls(CallReader& r, Controls& c) {
    c.resize(r.count());
    for (std::vector<Vec3>& y : c) r.points(y);
}

}  // namespace

bool writeNcpCall(const std::string& path, const RelaxRequest& request, const KnitProblemParams& problem, const NcpSolveParams& solve, std::string* err) {
    const MeshSurface* mesh = dynamic_cast<const MeshSurface*>(request.surface);
    if (request.surface && !mesh) {
        if (err) *err = "relax call capture: only a mesh band (MeshSurface) or no band can be written";
        return false;
    }
    if (request.bind) {
        if (err) *err = "relax call capture: a surface binding cannot be written (version 1)";
        return false;
    }
    CallWriter w(path);
    if (!w.ok()) {
        if (err) *err = "relax call capture: cannot write " + path;
        return false;
    }
    w.bytes(kMagic, sizeof kMagic);
    w.count(request.cs.yarns.size());
    for (const YarnPath& y : request.cs.yarns) {
        const YarnPath::Tables& t = y.tables();
        w.points(y.curve.ctrl);
        w.doubles(t.tu);
        w.doubles(t.ts);
        w.doubles(t.preLo);
        w.doubles(t.preHi);
        w.number(t.length);
        w.number(t.lengthErr);
    }
    writeControls(w, request.P);
    writeControls(w, request.Prest);
    w.count(request.freeCtrl.size());
    for (const std::vector<char>& y : request.freeCtrl) w.flags(y);
    w.integer(mesh ? 1 : 0);
    if (mesh) {
        const TriMesh& m = mesh->mesh();
        w.points(m.pos);
        w.count(m.tri.size());
        for (const std::array<int, 3>& t : m.tri)
            for (int k = 0; k < 3; ++k) w.integer(t[k]);
        w.points(m.vnormal);
    }
    const RelaxMoveParams& mp = request.mp;
    for (double v : {mp.thr, mp.kappaMax, mp.window, mp.lengthStiffness, mp.bandHalf, mp.maxDisplacement}) w.number(v);
    for (long long v : {static_cast<long long>(mp.iterations), static_cast<long long>(mp.refresh), static_cast<long long>(mp.samplesPerSpan),
                        static_cast<long long>(mp.call.round), static_cast<long long>(mp.call.stage)})
        w.integer(v);
    const std::string phase = mp.call.phase ? mp.call.phase : "";
    w.count(phase.size());
    w.bytes(phase.data(), phase.size());
    w.count(sizeof problem);
    w.bytes(&problem, sizeof problem);
    w.count(sizeof solve);
    w.bytes(&solve, sizeof solve);
    w.finish();
    if (!w.ok()) {
        if (err) *err = "relax call capture: writing " + path + " failed";
        return false;
    }
    return true;
}

bool readNcpCall(const std::string& path, NcpCallRecord& out, std::string* err) {
    CallReader r(path);
    if (!r.good()) {
        if (err) *err = "relax call file: cannot open " + path;
        return false;
    }
    char magic[8] = {};
    r.bytes(magic, sizeof magic);
    if (std::memcmp(magic, kMagic, sizeof kMagic) != 0) {
        if (err) *err = "relax call file: " + path + " is not a version 1 relax call";
        return false;
    }
    out = NcpCallRecord{};
    out.curves.yarns.resize(r.count());
    for (YarnPath& y : out.curves.yarns) {
        YarnPath::Tables t;
        r.points(y.curve.ctrl);
        r.doubles(t.tu);
        r.doubles(t.ts);
        r.doubles(t.preLo);
        r.doubles(t.preHi);
        t.length = r.number();
        t.lengthErr = r.number();
        y.setTables(std::move(t));
    }
    out.curves.kind = CurveKind::Surface;
    readControls(r, out.P);
    readControls(r, out.Prest);
    out.freeCtrl.resize(r.count());
    for (std::vector<char>& y : out.freeCtrl) r.flags(y);
    out.hasMesh = r.integer() == 1;
    if (out.hasMesh) {
        r.points(out.mesh.pos);
        out.mesh.tri.resize(r.count());
        for (std::array<int, 3>& t : out.mesh.tri)
            for (int k = 0; k < 3; ++k) t[k] = static_cast<int>(r.integer());
        r.points(out.mesh.vnormal);
    }
    RelaxMoveParams& mp = out.move;
    mp.thr = r.number();
    mp.kappaMax = r.number();
    mp.window = r.number();
    mp.lengthStiffness = r.number();
    mp.bandHalf = r.number();
    mp.maxDisplacement = r.number();
    mp.iterations = static_cast<int>(r.integer());
    mp.refresh = static_cast<int>(r.integer());
    mp.samplesPerSpan = static_cast<int>(r.integer());
    mp.call.round = static_cast<int>(r.integer());
    mp.call.stage = static_cast<int>(r.integer());
    out.phase.resize(r.count());
    r.bytes(out.phase.data(), out.phase.size());
    mp.call.phase = nullptr;  // the caller points it at out.phase once the record stays put
    const bool problemSize = r.count() == sizeof out.problem;
    r.bytes(&out.problem, sizeof out.problem);
    const bool solveSize = r.count() == sizeof out.solve;
    r.bytes(&out.solve, sizeof out.solve);
    if (!problemSize || !solveSize) {
        if (err) *err = "relax call file: " + path + " was written by a build with other parameter records";
        return false;
    }
    if (!r.digestMatches()) {
        if (err) *err = "relax call file: " + path + " is truncated or corrupted (digest)";
        return false;
    }
    return true;
}

} // namespace knit::gpu
