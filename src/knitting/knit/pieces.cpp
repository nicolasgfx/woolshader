// The scene pieces manifest, the gauge and the scene <-> library mapping (pieces.h).
#include "knitting/knit/pieces.h"

#include "knitting/layout/chart.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

namespace knit {

bool readObj(const std::string& path, TriMesh& m, std::string* err) {
    std::ifstream in(path);
    if (!in) {
        if (err) *err = "cannot open " + path;
        return false;
    }
    m = TriMesh{};
    std::vector<Vec2> vt;
    std::string line;
    int lineNo = 0;
    auto fail = [&](const std::string& what) {
        if (err) *err = path + ":" + std::to_string(lineNo) + ": " + what;
        return false;
    };
    while (std::getline(in, line)) {
        ++lineNo;
        if (line.size() < 2 || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "v") {
            Vec3 p;
            if (!(ss >> p.x >> p.y >> p.z)) return fail("bad vertex");
            m.pos.push_back(p);
        } else if (tag == "vt") {
            Vec2 t;
            if (!(ss >> t.x >> t.y)) return fail("bad texture coordinate");
            vt.push_back(t);
        } else if (tag == "f") {
            std::vector<int> vi, ti;
            std::string tok;
            while (ss >> tok) {
                int v = 0, t = 0;
                const char* s = tok.c_str();
                char* e = nullptr;
                v = static_cast<int>(std::strtol(s, &e, 10));
                if (*e == '/' && e[1] != '/') t = static_cast<int>(std::strtol(e + 1, &e, 10));
                if (v < 0) v = static_cast<int>(m.pos.size()) + v + 1;
                if (t < 0) t = static_cast<int>(vt.size()) + t + 1;
                if (v < 1 || v > static_cast<int>(m.pos.size())) return fail("vertex index out of range");
                if (t > static_cast<int>(vt.size())) return fail("texture index out of range");
                vi.push_back(v - 1);
                ti.push_back(t - 1);
            }
            if (vi.size() < 3) return fail("face with fewer than 3 vertices");
            for (size_t k = 1; k + 1 < vi.size(); ++k) {
                m.tri.push_back({vi[0], vi[k], vi[k + 1]});
                if (ti[0] >= 0 && ti[k] >= 0 && ti[k + 1] >= 0) m.uv.push_back({vt[ti[0]], vt[ti[k]], vt[ti[k + 1]]});
            }
        }
    }
    if (!m.uv.empty() && m.uv.size() != m.tri.size()) return fail("texture coordinates on some faces only");
    m.computeNormals();
    return true;
}

namespace {

bool parseVec(const std::string& s, Vec3& v) {
    const char* p = s.c_str();
    char* e = nullptr;
    double c[3];
    for (int k = 0; k < 3; ++k) {
        c[k] = std::strtod(p, &e);
        if (e == p || (k < 2 && *e != ',')) return false;
        p = e + 1;
    }
    v = {c[0], c[1], c[2]};
    return *e == 0;
}

std::string dirOf(const std::string& path) {
    size_t k = path.find_last_of("/\\");
    return k == std::string::npos ? std::string(".") : path.substr(0, k);
}

} // namespace

std::string piecePath(const Manifest& m, const std::string& file) { return m.dir + "/" + file; }

bool readManifest(const std::string& path, Manifest& m, std::string* err) {
    std::ifstream in(path);
    if (!in) {
        if (err) *err = "cannot open " + path;
        return false;
    }
    m = Manifest{};
    m.dir = dirOf(path);
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string kind, file, tok;
        ss >> kind >> file;
        std::map<std::string, std::string> kv;
        while (ss >> tok) {
            size_t e = tok.find('=');
            if (e != std::string::npos) kv[tok.substr(0, e)] = tok.substr(e + 1);
        }
        auto bad = [&](const char* what) {
            if (err) *err = path + ":" + std::to_string(lineNo) + ": " + what;
            return false;
        };
        if (kind == "piece") {
            PieceEntry p;
            p.file = file;
            p.object = kv["object"];
            p.material = kv["material"];
            p.island = std::atoi(kv["island"].c_str());
            p.rings = std::atoi(kv["rings"].c_str());
            p.segments = std::atoi(kv["segments"].c_str());
            p.area = std::atof(kv["area"].c_str());
            if (!parseVec(kv["centroid"], p.centroid) || !parseVec(kv["poleA"], p.poleA) ||
                !parseVec(kv["poleB"], p.poleB))
                return bad("bad vector field");
            m.pieces.push_back(p);
        } else if (kind == "extras") {
            ExtrasEntry x;
            x.file = file;
            x.object = kv["object"];
            x.islands = std::atoi(kv["islands"].c_str());
            m.extras.push_back(x);
        } else {
            return bad("unknown record");
        }
    }
    return true;
}

bool computeGauge(const Manifest& m, const GaugeRequest& request, Gauge& g, std::string* err) {
    const std::string& file = request.file;
    int rings = 0;
    for (const PieceEntry& pe : m.pieces)
        if (pe.file == file) rings = pe.rings;
    if (!rings) {
        if (err) *err = "gauge piece " + file + " is not in the manifest";
        return false;
    }
    TriMesh mesh;
    if (!readObj(piecePath(m, file), mesh, err)) return false;
    Chart ch(mesh, rings);
    g.file = file;
    g.stitches = request.stitches;
    g.w = request.w;
    const Chart::RingAt longest = ch.longestRing();
    g.aAt = longest.a;
    g.ringLength = longest.length;
    g.sigma = request.stitches * request.w / g.ringLength;
    return true;
}

bool loadPiece(const Manifest& m, int index, double sigma, LoadedPiece& out, std::string* err) {
    out = LoadedPiece{};
    out.entry = m.pieces[index];
    if (!readObj(piecePath(m, out.entry.file), out.mesh, err)) return false;
    out.map.sigma = sigma;
    out.map.centre = out.mesh.centroid();
    out.mesh.transform(sigma, out.map.centre);
    return true;
}

} // namespace knit
