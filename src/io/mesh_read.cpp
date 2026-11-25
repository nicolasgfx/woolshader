// Foreign meshes, OBJ with its MTL and PLY, read within the import limits (mesh_read.h has the bounds).
#include "io/mesh_read.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <map>
#include <sstream>
#include <string_view>

#include "io/files.h"
#include "model/cancel.h"

namespace knit {

namespace fs = std::filesystem;

namespace {

const char kDecimateHint[] = "Decimate the mesh first (with a modelling program's decimate tool); detail below a stitch (about 5 mm) does not show in the knit.";

std::string lower(std::string s) {
    for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

std::string mib(uint64_t bytes) {
    char b[32];
    std::snprintf(b, sizeof b, "%.1f MB", double(bytes) / 1048576.0);
    return b;
}

bool finite(const Vec3& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }

// Line cursor over a buffer, with a cancellation checkpoint and progress every MiB.
struct Lines {
    const std::string& s;
    size_t pos = 0;
    int line = 0;
    size_t lastTick = 0;
    const ImportProgress* progress;
    const char* stage;
    Lines(const std::string& buf, const ImportProgress* p, const char* st) : s(buf), progress(p), stage(st) {}
    bool next(std::string_view& out) {
        if (pos >= s.size()) return false;
        size_t e = s.find('\n', pos);
        if (e == std::string::npos) e = s.size();
        size_t end = e;
        if (end > pos && s[end - 1] == '\r') --end;
        out = std::string_view(s).substr(pos, end - pos);
        pos = e + 1;
        ++line;
        if (pos - lastTick >= (1u << 20)) {
            lastTick = pos;
            cancelCheckpoint();
            if (progress && *progress) (*progress)(stage, double(pos) / double(std::max<size_t>(1, s.size())));
        }
        return true;
    }
};

std::vector<std::string_view> split(std::string_view l) {
    std::vector<std::string_view> t;
    size_t i = 0;
    while (i < l.size()) {
        while (i < l.size() && (l[i] == ' ' || l[i] == '\t')) ++i;
        if (i >= l.size()) break;
        size_t j = i;
        while (j < l.size() && l[j] != ' ' && l[j] != '\t') ++j;
        t.push_back(l.substr(i, j - i));
        i = j;
    }
    return t;
}

bool parseDouble(std::string_view t, double& v) {
    std::string s(t);
    char* end = nullptr;
    v = std::strtod(s.c_str(), &end);
    return end && *end == 0 && !s.empty();
}

bool parseLong(std::string_view t, long long& v) {
    std::string s(t);
    char* end = nullptr;
    v = std::strtoll(s.c_str(), &end, 10);
    return end && *end == 0 && !s.empty();
}

std::string trim(std::string_view v) {
    size_t a = 0, b = v.size();
    while (a < b && std::isspace(static_cast<unsigned char>(v[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(v[b - 1]))) --b;
    return std::string(v.substr(a, b - a));
}

std::string resolve(const fs::path& dir, const std::string& ref) {
    std::string r = ref;
    for (char& c : r)
        if (c == '\\') c = '/';
    fs::path p(r);
    if (p.is_relative()) p = dir / p;
    std::error_code ec;
    const fs::path c = fs::weakly_canonical(p, ec);
    return (ec ? p.lexically_normal() : c).string();
}

std::vector<std::string> ignoredLines(const std::map<std::string, int>& ignored, const char* where) {
    std::vector<std::string> out;
    for (const auto& [k, n] : ignored) out.push_back(std::string(where) + " " + k + " (" + countText(n) + "x)");
    return out;
}

// A referenced file: a dependency of the import when it exists, listed as missing when not.
bool noteDependency(ImportedMesh& m, const std::string& path) {
    std::error_code ec;
    const bool found = fs::exists(path, ec);
    (found ? m.dependencies : m.missingDependencies).push_back(path);
    return found;
}

// --- MTL

// map_Kd [options] file: options of the MTL spec with their argument counts; the file name is
// the rest of the line (it may contain spaces).
std::string mapFileName(std::string_view rest) {
    static const std::map<std::string, int> opt = {{"-blendu", 1}, {"-blendv", 1}, {"-boost", 1}, {"-mm", 2}, {"-o", 3},      {"-s", 3},
                                                   {"-t", 3},      {"-texres", 1}, {"-clamp", 1}, {"-bm", 1}, {"-imfchan", 1}, {"-type", 1}};
    std::string s = trim(rest);
    for (;;) {
        if (s.empty() || s[0] != '-') return s;
        const size_t sp = s.find_first_of(" \t");
        if (sp == std::string::npos) return "";
        const std::string o = s.substr(0, sp);
        s = trim(std::string_view(s).substr(sp));
        const auto it = opt.find(o);
        const int n = it == opt.end() ? 1 : it->second;
        for (int k = 0; k < n; ++k) {
            // -o / -s / -t take 1 to 3 numbers: stop at the first non-number.
            const size_t e = s.find_first_of(" \t");
            const std::string tok = s.substr(0, e);
            double d;
            if ((o == "-o" || o == "-s" || o == "-t") && k > 0 && !parseDouble(tok, d)) break;
            s = e == std::string::npos ? std::string() : trim(std::string_view(s).substr(e));
        }
    }
}

// The index of a named material, added if the file has not named it yet.
int materialIndex(ImportedMesh& m, std::map<std::string, int>& index, const std::string& name) {
    const auto [it, added] = index.emplace(name, static_cast<int>(m.materials.size()));
    if (added) {
        m.materials.push_back({});
        m.materials.back().name = name;
    }
    return it->second;
}

void readMtl(const std::string& path, ImportedMesh& m, std::map<std::string, int>& matIndex, std::map<std::string, int>& ignored) {
    std::string buf, e;
    if (!readWholeFile(path, buf, &e)) return;
    const fs::path dir = fs::path(path).parent_path();
    Lines L(buf, nullptr, "");
    std::string_view l;
    ImportMaterial* cur = nullptr;
    while (L.next(l)) {
        const std::vector<std::string_view> t = split(l);
        if (t.empty() || t[0][0] == '#') continue;
        const std::string k(t[0]);
        if (k == "newmtl" && t.size() >= 2) {
            cur = &m.materials[materialIndex(m, matIndex, trim(l.substr(l.find(t[1]))))];
        } else if (k == "Kd" && cur && t.size() >= 4) {
            double r, g, b;
            if (parseDouble(t[1], r) && parseDouble(t[2], g) && parseDouble(t[3], b)) {
                cur->kd = {std::clamp(r, 0.0, 1.0), std::clamp(g, 0.0, 1.0), std::clamp(b, 0.0, 1.0)};
                cur->hasKd = true;
            }
        } else if (k == "map_Kd" && cur) {
            cur->textureRef = mapFileName(l.substr(t[0].size()));
            if (!cur->textureRef.empty()) {
                const std::string p = resolve(dir, cur->textureRef);
                if (noteDependency(m, p)) cur->texture = p;
                else cur->textureMissing = true;
            }
        } else {
            ++ignored[k];
        }
    }
}

// --- OBJ

// What reading one OBJ file has gathered, and where it is.
struct ObjFile {
    const std::string& path;
    const ImportLimits& lim;
    ImportedMesh& m;
    std::string* err;
    const Lines& lines;
    std::vector<Vec2> vts;
    std::map<std::string, int> matIndex, ignored, mtlIgnored;
    int curMat = -1;
    std::string curNode = "(unnamed)";  // the last o or g; it becomes a node with its first face
    std::map<std::string, int> nodeIndex;
    bool anyFaceUv = false;
    long long vnCount = 0;
    std::vector<int> faceVertices, faceUvs;

    std::string where() const { return path + ", line " + std::to_string(lines.line); }
    bool refuse(const std::string& why) const { return fail(err, where() + ": " + why); }
};

bool readObjVertex(ObjFile& o, const std::vector<std::string_view>& t) {
    if (t.size() < 4) return o.refuse("a vertex needs 3 coordinates");
    Vec3 p;
    if (!parseDouble(t[1], p.x) || !parseDouble(t[2], p.y) || !parseDouble(t[3], p.z)) return o.refuse("unreadable vertex coordinate");
    if (!finite(p)) return o.refuse("a vertex coordinate is not a finite number (NaN or infinity)");
    ImportedMesh& m = o.m;
    if (m.pos.size() >= o.lim.maxVertices)
        return fail(o.err, "more than " + countText(o.lim.maxVertices) + " vertices (" + o.where() + "): the import limit. " + kDecimateHint);
    m.pos.push_back(p);
    if (t.size() >= 7) {  // the common "v x y z r g b" extension
        Vec3 c;
        if (parseDouble(t[4], c.x) && parseDouble(t[5], c.y) && parseDouble(t[6], c.z)) {
            if (m.vcolor.size() + 1 < m.pos.size()) m.vcolor.resize(m.pos.size() - 1, Vec3{0.8, 0.8, 0.8});
            m.vcolor.push_back({std::clamp(c.x, 0.0, 1.0), std::clamp(c.y, 0.0, 1.0), std::clamp(c.z, 0.0, 1.0)});
        }
    }
    return true;
}

bool readObjTexCoord(ObjFile& o, const std::vector<std::string_view>& t) {
    if (t.size() < 2) return o.refuse("a texture coordinate needs at least u");
    Vec2 uv;
    if (!parseDouble(t[1], uv.x) || (t.size() >= 3 && !parseDouble(t[2], uv.y))) return o.refuse("unreadable texture coordinate");
    if (!std::isfinite(uv.x) || !std::isfinite(uv.y)) return o.refuse("a texture coordinate is not finite");
    o.vts.push_back(uv);
    return true;
}

// The node of the current name, made at its first face: o and g are one rule and a name without faces is no node, as in
// three.js's OBJLoader. Faces under one name stay one node, wherever they stand in the file.
int nodeOfFace(ObjFile& o) {
    int& node = o.nodeIndex.try_emplace(o.curNode, -1).first->second;
    if (node < 0) {
        node = static_cast<int>(o.m.nodes.size());
        o.m.nodes.push_back({o.curNode});
    }
    return node;
}

// One corner of a face, "v", "v/vt" or "v/vt/vn", with negative indices counted from the end.
bool readObjCorner(ObjFile& o, std::string_view tok) {
    const size_t s1 = tok.find('/');
    long long vi = 0, ti = 0;
    if (!parseLong(tok.substr(0, s1), vi) || vi == 0) return o.refuse("unreadable vertex index '" + std::string(tok) + "'");
    vi = vi < 0 ? static_cast<long long>(o.m.pos.size()) + vi : vi - 1;
    if (vi < 0 || vi >= static_cast<long long>(o.m.pos.size()))
        return o.refuse("vertex index " + std::string(tok.substr(0, s1)) + " refers to no vertex (" + countText(o.m.pos.size()) + " defined so far)");
    o.faceVertices.push_back(static_cast<int>(vi));
    if (s1 == std::string_view::npos) return true;
    const size_t s2 = tok.find('/', s1 + 1);
    const std::string_view tt = tok.substr(s1 + 1, s2 == std::string_view::npos ? std::string_view::npos : s2 - s1 - 1);
    if (tt.empty()) return true;
    if (!parseLong(tt, ti) || ti == 0) return o.refuse("unreadable texture index '" + std::string(tok) + "'");
    ti = ti < 0 ? static_cast<long long>(o.vts.size()) + ti : ti - 1;
    if (ti < 0 || ti >= static_cast<long long>(o.vts.size())) return o.refuse("texture index " + std::string(tt) + " refers to no texture coordinate");
    o.faceUvs.push_back(static_cast<int>(ti));
    return true;
}

bool readObjFace(ObjFile& o, const std::vector<std::string_view>& t) {
    if (t.size() < 4) return o.refuse("a face needs at least 3 corners");
    const int n = static_cast<int>(t.size()) - 1;
    if (n > o.lim.maxPolygon)
        return o.refuse("a face with " + std::to_string(n) + " corners (the limit is " + std::to_string(o.lim.maxPolygon) + "); triangulate the mesh before importing");
    o.faceVertices.clear();
    o.faceUvs.clear();
    for (int c = 1; c <= n; ++c)
        if (!readObjCorner(o, t[c])) return false;
    const std::vector<int>& fv = o.faceVertices;
    const std::vector<int>& fvt = o.faceUvs;
    const bool faceUv = static_cast<int>(fvt.size()) == n;
    if (!fvt.empty() && !faceUv) return o.refuse("texture indices on some corners of a face only");
    ImportedMesh& m = o.m;
    ++m.faces;
    if (n > 3) ++m.polygons;
    m.maxPolygon = std::max(m.maxPolygon, n);
    if (m.tri.size() + static_cast<size_t>(n - 2) > o.lim.maxTriangles)
        return fail(o.err, "more than " + countText(o.lim.maxTriangles) + " triangles (" + o.where() + "): the import limit. " + kDecimateHint);
    if (faceUv) o.anyFaceUv = true;
    else ++m.facesWithoutUv;
    if (n == 4) m.quadFirstTriangle.push_back(static_cast<int>(m.tri.size()));
    const int node = nodeOfFace(o);
    for (int c = 1; c + 1 < n; ++c) {  // fan
        m.tri.push_back({fv[0], fv[c], fv[c + 1]});
        m.triMaterial.push_back(o.curMat);
        m.nodeOfTriangle.push_back(node);
        if (faceUv) m.uv.push_back({o.vts[fvt[0]], o.vts[fvt[c]], o.vts[fvt[c + 1]]});
        else m.uv.push_back({Vec2{}, Vec2{}, Vec2{}});
    }
    return true;
}

bool readObjFile(const std::string& path, const std::string& buf, const ImportLimits& lim, ImportedMesh& m, std::string* err,
                 const ImportProgress& progress) {
    m.format = "OBJ";
    const fs::path dir = fs::path(path).parent_path();
    Lines L(buf, &progress, "reading");
    ObjFile o{path, lim, m, err, L};
    std::string_view l;
    while (L.next(l)) {
        const std::vector<std::string_view> t = split(l);
        if (t.empty() || t[0][0] == '#') continue;
        const std::string_view k = t[0];
        if (k == "v") {
            if (!readObjVertex(o, t)) return false;
        } else if (k == "vt") {
            if (!readObjTexCoord(o, t)) return false;
        } else if (k == "vn") {
            ++o.vnCount;
        } else if (k == "f") {
            if (!readObjFace(o, t)) return false;
        } else if (k == "usemtl") {  // perhaps before (or without) its definition: Kd unknown
            o.curMat = materialIndex(m, o.matIndex, t.size() >= 2 ? trim(l.substr(l.find(t[1]))) : "");
        } else if (k == "mtllib") {
            const std::string p = resolve(dir, trim(l.substr(k.size())));
            if (noteDependency(m, p)) readMtl(p, m, o.matIndex, o.mtlIgnored);
        } else if (k == "o" || k == "g") {
            std::string name = t.size() >= 2 ? trim(l.substr(l.find(t[1]))) : "(unnamed)";
            std::replace(name.begin(), name.end(), '|', '_');  // | separates names in the settings (import.plain_nodes)
            o.curNode = name;
        } else if (k != "s") {  // smoothing: no effect on the knit
            ++o.ignored[std::string(k)];
        }
    }
    if (o.vnCount) m.ignored.push_back("vn normals (" + countText(o.vnCount) + "; recomputed from the geometry)");
    for (const std::string& s : ignoredLines(o.ignored, "OBJ statement")) m.ignored.push_back(s);
    for (const std::string& s : ignoredLines(o.mtlIgnored, "MTL property")) m.ignored.push_back(s + "; only the base colour (Kd, map_Kd) is used");
    if (!o.anyFaceUv) m.uv.clear();
    if (m.vcolor.size() != m.pos.size()) m.vcolor.clear();
    return true;
}

// --- PLY

struct PlyProp {
    std::string name, type, countType;  // countType non-empty: a list
    bool list() const { return !countType.empty(); }
};
struct PlyElement {
    std::string name;
    uint64_t count = 0;
    std::vector<PlyProp> props;
};

int typeSize(const std::string& t) {
    if (t == "char" || t == "uchar" || t == "int8" || t == "uint8") return 1;
    if (t == "short" || t == "ushort" || t == "int16" || t == "uint16") return 2;
    if (t == "int" || t == "uint" || t == "int32" || t == "uint32" || t == "float" || t == "float32") return 4;
    if (t == "double" || t == "float64") return 8;
    return 0;
}

bool isFloatType(const std::string& t) { return t == "float" || t == "float32" || t == "double" || t == "float64"; }

struct BinReader {
    const std::string& s;
    size_t pos;
    bool bigEndian;
    bool ok = true;
    double read(const std::string& t) {
        const int n = typeSize(t);
        if (pos + n > s.size()) {
            ok = false;
            return 0;
        }
        unsigned char b[8];
        std::memcpy(b, s.data() + pos, n);
        pos += n;
        if (bigEndian) std::reverse(b, b + n);
        if (t == "char" || t == "int8") return double(static_cast<int8_t>(b[0]));
        if (t == "uchar" || t == "uint8") return double(b[0]);
        if (t == "short" || t == "int16") { int16_t v; std::memcpy(&v, b, 2); return v; }
        if (t == "ushort" || t == "uint16") { uint16_t v; std::memcpy(&v, b, 2); return v; }
        if (t == "int" || t == "int32") { int32_t v; std::memcpy(&v, b, 4); return v; }
        if (t == "uint" || t == "uint32") { uint32_t v; std::memcpy(&v, b, 4); return v; }
        if (t == "float" || t == "float32") { float v; std::memcpy(&v, b, 4); return v; }
        double v;
        std::memcpy(&v, b, 8);
        return v;
    }
};

// The data section's values one by one, from ascii or binary: everything after this is one path.
class PlyValues {
public:
    PlyValues(const std::string& buf, size_t start, bool ascii, bool bigEndian) : ascii_(ascii), binary_{buf, start, bigEndian} {
        if (ascii_) text_.str(buf.substr(start));
    }
    double next(const std::string& type) {
        double v = 0;
        if (ascii_) text_ >> v;
        else v = binary_.read(type);
        return v;
    }
    bool ok() const { return ascii_ ? static_cast<bool>(text_) : binary_.ok; }
    const char* failure() const { return ascii_ ? "unreadable" : "truncated"; }

private:
    bool ascii_;
    std::istringstream text_;
    BinReader binary_;
};

struct PlyHeader {
    std::vector<PlyElement> elements;
    const PlyElement* vertex = nullptr;
    const PlyElement* face = nullptr;
    std::string textureFile;
    size_t dataStart = 0;
    bool ascii = false, bigEndian = false;
};

bool parsePlyHeader(const std::string& path, const std::string& buf, ImportedMesh& m, PlyHeader& h, std::string* err) {
    const size_t endHeader = buf.find("end_header");
    if (buf.compare(0, 3, "ply") != 0 || endHeader == std::string::npos)
        return fail(err, path + ": not a PLY file (no 'ply' magic or no 'end_header')");
    h.dataStart = buf.find('\n', endHeader);
    if (h.dataStart == std::string::npos) return fail(err, path + ": the header has no end");
    ++h.dataStart;
    std::istringstream hs(buf.substr(0, endHeader));
    std::string line, format;
    while (std::getline(hs, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::istringstream ls(line);
        std::string k;
        ls >> k;
        if (k == "format") {
            std::string ver;
            ls >> format >> ver;
        } else if (k == "comment" || k == "obj_info") {
            std::string w;
            ls >> w;
            if (k == "comment" && (w == "TextureFile" || w == "texturefile")) std::getline(ls >> std::ws, h.textureFile);
        } else if (k == "element") {
            PlyElement e;
            std::string cnt;
            ls >> e.name >> cnt;
            long long c = 0;
            if (!parseLong(cnt, c) || c < 0) return fail(err, path + ": bad element count '" + cnt + "'");
            e.count = static_cast<uint64_t>(c);
            h.elements.push_back(e);
        } else if (k == "property") {
            if (h.elements.empty()) return fail(err, path + ": a property before any element");
            PlyProp p;
            std::string a;
            ls >> a;
            if (a == "list") {
                ls >> p.countType >> p.type >> p.name;
                if (!typeSize(p.countType) || !typeSize(p.type)) return fail(err, path + ": unknown list type in '" + line + "'");
            } else {
                p.type = a;
                ls >> p.name;
                if (!typeSize(p.type)) return fail(err, path + ": unknown property type in '" + line + "'");
            }
            h.elements.back().props.push_back(p);
        }
    }
    if (format == "ascii") h.ascii = true, m.format = "PLY ascii";
    else if (format == "binary_little_endian") m.format = "PLY binary little endian";
    else if (format == "binary_big_endian") h.bigEndian = true, m.format = "PLY binary big endian";
    else return fail(err, path + ": unknown PLY format '" + format + "'");
    for (const PlyElement& e : h.elements) {
        if (e.name == "vertex") h.vertex = &e;
        if (e.name == "face") h.face = &e;
    }
    if (!h.vertex) return fail(err, path + ": no vertex element");
    if (!h.face || h.face->count == 0)
        return fail(err, path + ": no faces (a point cloud). Knitting needs a surface: reconstruct one first (e.g. Poisson reconstruction in MeshLab).");
    return true;
}

// The header's counts against the limits, before anything is allocated.
bool checkPlyLimits(const std::string& path, const std::string& buf, const PlyHeader& h, const ImportLimits& lim, std::string* err) {
    if (h.vertex->count > lim.maxVertices)
        return fail(err, path + " declares " + countText(h.vertex->count) + " vertices; the import limit is " + countText(lim.maxVertices) +
                             ". Decimate the mesh first (detail below a stitch, about 5 mm, does not show in the knit).");
    if (h.face->count > lim.maxTriangles)
        return fail(err, path + " declares " + countText(h.face->count) + " faces; the import limit is " + countText(lim.maxTriangles) +
                             " triangles. Decimate the mesh first.");
    if (h.ascii) return true;
    // Binary: the smallest payload the header implies must fit the file (truncation).
    long double need = 0;
    for (const PlyElement& e : h.elements) {
        uint64_t row = 0;
        for (const PlyProp& p : e.props) row += p.list() ? typeSize(p.countType) : typeSize(p.type);
        need += static_cast<long double>(row) * e.count;
    }
    if (need > static_cast<long double>(buf.size() - h.dataStart))
        return fail(err, path + " is truncated: the header declares at least " + mib(static_cast<uint64_t>(need)) + " of data, the file holds " +
                             mib(buf.size() - h.dataStart) + " after the header");
    return true;
}

// Which property of the vertex and face elements holds what (-1: none).
struct PlyColumns {
    int x = -1, y = -1, z = -1, u = -1, v = -1, r = -1, g = -1, b = -1;
    int indices = -1, texcoord = -1;
    bool vertexUv() const { return u >= 0 && v >= 0; }
    bool vertexColour() const { return r >= 0 && g >= 0 && b >= 0; }
};

bool findPlyColumns(const std::string& path, const PlyHeader& h, PlyColumns& c, std::map<std::string, int>& ignored, std::string* err) {
    for (size_t k = 0; k < h.vertex->props.size(); ++k) {
        const std::string n = h.vertex->props[k].name;
        const int i = static_cast<int>(k);
        if (n == "x") c.x = i;
        else if (n == "y") c.y = i;
        else if (n == "z") c.z = i;
        else if (n == "s" || n == "u" || n == "texture_u" || n == "texture_s") c.u = i;
        else if (n == "t" || n == "v" || n == "texture_v" || n == "texture_t") c.v = i;
        else if (n == "red" || n == "r" || n == "diffuse_red") c.r = i;
        else if (n == "green" || n == "g" || n == "diffuse_green") c.g = i;
        else if (n == "blue" || n == "b" || n == "diffuse_blue") c.b = i;
        else if (n == "nx" || n == "ny" || n == "nz") ++ignored["vertex normals (recomputed)"];
        else ++ignored["vertex property " + n];
    }
    if (c.x < 0 || c.y < 0 || c.z < 0) return fail(err, path + ": the vertex element has no x / y / z");
    for (size_t k = 0; k < h.face->props.size(); ++k) {
        const std::string n = h.face->props[k].name;
        if ((n == "vertex_indices" || n == "vertex_index") && h.face->props[k].list()) c.indices = static_cast<int>(k);
        else if (n == "texcoord" && h.face->props[k].list()) c.texcoord = static_cast<int>(k);
        else ++ignored["face property " + n];
    }
    if (c.indices < 0) return fail(err, path + ": the face element has no vertex_indices list");
    for (const PlyElement& e : h.elements)
        if (e.name != "vertex" && e.name != "face" && e.count) ++ignored["element " + e.name];
    return true;
}

// One element row: its scalar properties, and the list values the mesh needs.
struct PlyRow {
    std::vector<double> scalars;
    std::vector<int> indices;
    std::vector<double> texcoord;
};

bool readPlyRow(const std::string& path, const PlyHeader& h, const PlyColumns& c, const PlyElement& e, uint64_t r, PlyValues& data, PlyRow& row, std::string* err) {
    row.scalars.assign(e.props.size(), 0);
    row.indices.clear();
    row.texcoord.clear();
    for (size_t k = 0; k < e.props.size(); ++k) {
        const PlyProp& p = e.props[k];
        if (!p.list()) {
            row.scalars[k] = data.next(p.type);
            if (!data.ok()) return fail(err, path + ": " + data.failure() + " data in element " + e.name + " " + countText(r));
            continue;
        }
        const double cnt = data.next(p.countType);
        if (!data.ok() || cnt < 0 || cnt > 4096) return fail(err, path + ": " + data.failure() + " list in element " + e.name + " " + countText(r));
        const int n = static_cast<int>(cnt);
        for (int i = 0; i < n; ++i) {
            const double v = data.next(p.type);
            if (!data.ok()) return fail(err, path + ": truncated list data in element " + e.name);
            if (&e == h.face && static_cast<int>(k) == c.indices) row.indices.push_back(static_cast<int>(v));
            else if (&e == h.face && static_cast<int>(k) == c.texcoord) row.texcoord.push_back(v);
        }
    }
    return true;
}

bool addPlyVertex(const std::string& path, const PlyHeader& h, const PlyColumns& c, const PlyRow& row, uint64_t r, std::vector<Vec2>& vertexUv, ImportedMesh& m, std::string* err) {
    const Vec3 p{row.scalars[c.x], row.scalars[c.y], row.scalars[c.z]};
    if (!finite(p)) return fail(err, path + ": vertex " + countText(r) + " has a coordinate that is not finite (NaN or infinity)");
    m.pos.push_back(p);
    if (c.vertexUv()) vertexUv.push_back({row.scalars[c.u], row.scalars[c.v]});
    if (c.vertexColour()) {
        const double scale = !isFloatType(h.vertex->props[c.r].type) ? 1.0 / 255.0 : 1.0;
        m.vcolor.push_back({std::clamp(row.scalars[c.r] * scale, 0.0, 1.0), std::clamp(row.scalars[c.g] * scale, 0.0, 1.0),
                            std::clamp(row.scalars[c.b] * scale, 0.0, 1.0)});
    }
    return true;
}

bool addPlyFace(const std::string& path, const PlyHeader& h, const PlyColumns& c, const ImportLimits& lim, const PlyRow& row, uint64_t r,
                const std::vector<Vec2>& vertexUv, bool& anyFaceUv, ImportedMesh& m, std::string* err) {
    const std::vector<int>& idx = row.indices;
    const std::vector<double>& tc = row.texcoord;
    const int n = static_cast<int>(idx.size());
    if (n < 3) return fail(err, path + ": face " + countText(r) + " has " + std::to_string(n) + " corners");
    if (n > lim.maxPolygon)
        return fail(err, path + ": face " + countText(r) + " has " + std::to_string(n) + " corners (the limit is " + std::to_string(lim.maxPolygon) + "); triangulate the mesh first");
    for (int v : idx)
        if (v < 0 || v >= static_cast<int>(m.pos.size()))
            return fail(err, path + ": face " + countText(r) + " refers to vertex " + std::to_string(v) + ", which does not exist (" + countText(m.pos.size()) + " vertices)");
    ++m.faces;
    if (n > 3) ++m.polygons;
    m.maxPolygon = std::max(m.maxPolygon, n);
    if (m.tri.size() + static_cast<size_t>(n - 2) > lim.maxTriangles)
        return fail(err, path + ": more than " + countText(lim.maxTriangles) + " triangles after triangulation; decimate the mesh first");
    const bool faceUv = static_cast<int>(tc.size()) == 2 * n;
    if (faceUv) anyFaceUv = true;
    if (n == 4) m.quadFirstTriangle.push_back(static_cast<int>(m.tri.size()));
    for (int k = 1; k + 1 < n; ++k) {
        m.tri.push_back({idx[0], idx[k], idx[k + 1]});
        m.triMaterial.push_back(h.textureFile.empty() ? -1 : 0);
        if (faceUv) m.uv.push_back({Vec2{tc[0], tc[1]}, Vec2{tc[2 * k], tc[2 * k + 1]}, Vec2{tc[2 * k + 2], tc[2 * k + 3]}});
        else if (c.vertexUv()) m.uv.push_back({vertexUv[idx[0]], vertexUv[idx[k]], vertexUv[idx[k + 1]]});
        else m.uv.push_back({Vec2{}, Vec2{}, Vec2{}});
    }
    if (!faceUv && !c.vertexUv()) ++m.facesWithoutUv;
    return true;
}

bool readPlyElements(const std::string& path, const std::string& buf, const PlyHeader& h, const PlyColumns& c, const ImportLimits& lim, ImportedMesh& m,
                     std::string* err, const ImportProgress& progress) {
    m.pos.reserve(h.vertex->count);
    std::vector<Vec2> vertexUv;
    if (c.vertexUv()) vertexUv.reserve(h.vertex->count);
    if (c.vertexColour()) m.vcolor.reserve(h.vertex->count);
    PlyValues data(buf, h.dataStart, h.ascii, h.bigEndian);
    uint64_t done = 0;
    const uint64_t total = std::max<uint64_t>(1, h.vertex->count + h.face->count);
    PlyRow row;
    bool anyFaceUv = false;
    for (const PlyElement& e : h.elements) {
        for (uint64_t r = 0; r < e.count; ++r) {
            if ((++done & 0x3fff) == 0) {
                cancelCheckpoint();
                if (progress) progress("reading", double(done) / double(total));
            }
            if (!readPlyRow(path, h, c, e, r, data, row, err)) return false;
            if (&e == h.vertex && !addPlyVertex(path, h, c, row, r, vertexUv, m, err)) return false;
            if (&e == h.face && !addPlyFace(path, h, c, lim, row, r, vertexUv, anyFaceUv, m, err)) return false;
        }
        if (&e == h.vertex && c.vertexUv()) anyFaceUv = true;
    }
    if (!anyFaceUv) {
        m.uv.clear();
        m.facesWithoutUv = 0;
    }
    return true;
}

void attachPlyTexture(const std::string& path, const std::string& textureFile, ImportedMesh& m) {
    ImportMaterial mat;
    mat.name = "texture";
    mat.textureRef = textureFile;
    const std::string p = resolve(fs::path(path).parent_path(), textureFile);
    if (noteDependency(m, p)) mat.texture = p;
    else mat.textureMissing = true;
    m.materials.push_back(mat);
}

bool readPlyFile(const std::string& path, const std::string& buf, const ImportLimits& lim, ImportedMesh& m, std::string* err,
                 const ImportProgress& progress) {
    PlyHeader h;
    if (!parsePlyHeader(path, buf, m, h, err) || !checkPlyLimits(path, buf, h, lim, err)) return false;
    PlyColumns columns;
    std::map<std::string, int> ignored;
    if (!findPlyColumns(path, h, columns, ignored, err) || !readPlyElements(path, buf, h, columns, lim, m, err, progress)) return false;
    if (!h.textureFile.empty()) attachPlyTexture(path, h.textureFile, m);
    for (const std::string& s : ignoredLines(ignored, "PLY")) m.ignored.push_back(s);
    return true;
}

} // namespace

std::string countText(uint64_t n) {
    std::string s = std::to_string(n), out;
    for (size_t i = 0; i < s.size(); ++i) {
        if (i && (s.size() - i) % 3 == 0) out += ',';
        out += s[i];
    }
    return out;
}

bool readMeshFile(const std::string& path, const ImportLimits& lim, ImportedMesh& out, std::string* err, const ImportProgress& progress) {
    out = ImportedMesh{};
    out.path = path;
    std::error_code ec;
    const fs::file_status st = fs::status(path, ec);
    if (ec || !fs::exists(st)) return fail(err, path + " does not exist");
    if (!fs::is_regular_file(st)) return fail(err, path + " is not a file");
    out.fileBytes = fs::file_size(path, ec);
    if (ec) return fail(err, "cannot read the size of " + path);
    const std::string ext = lower(fs::path(path).extension().string());
    if (ext != ".obj" && ext != ".ply")
        return fail(err, fs::path(path).filename().string() + ": only OBJ and PLY files can be imported" +
                             (ext == ".glb" || ext == ".gltf" || ext == ".fbx" || ext == ".stl"
                                  ? " (export the model as OBJ from your modelling program)"
                                  : ""));
    // The file size first: nothing is read from a file beyond the limit.
    if (out.fileBytes > lim.maxFileBytes)
        return fail(err, fs::path(path).filename().string() + " is " + mib(out.fileBytes) + "; the import limit is " + mib(lim.maxFileBytes) +
                             " of mesh data. Knitting cost follows the surface area, not the triangle count, so a decimated copy (with a "
                             "modelling program's decimate tool, to under 2 million triangles) knits the same.");
    if (out.fileBytes == 0) return fail(err, path + " is empty");
    if (progress) progress("reading", 0);
    std::string buf;
    if (!readWholeFile(path, buf, err)) return false;
    cancelCheckpoint();
    const bool ok = ext == ".obj" ? readObjFile(path, buf, lim, out, err, progress) : readPlyFile(path, buf, lim, out, err, progress);
    if (!ok) return false;
    if (out.tri.empty()) return fail(err, path + ": no faces. Knitting needs a surface (a point cloud or a set of lines cannot be knitted).");
    return true;
}

std::vector<ColourImage> loadImportTextures(ImportedMesh& m, const ImportLimits& lim, std::vector<std::string>* warnings) {
    std::vector<ColourImage> out;
    std::map<std::string, size_t> seen;
    uint64_t bytes = 0;
    for (ImportMaterial& mat : m.materials) {
        if (mat.texture.empty()) continue;
        if (seen.count(mat.texture)) continue;
        cancelCheckpoint();
        const std::string file = fs::path(mat.texture).filename().string();
        int w = 0, h = 0;
        std::string why;
        if (!readImageSize(mat.texture, w, h, &why)) {
            if (warnings) warnings->push_back(file + ": not a readable image (" + why + "); material " + mat.name + " uses its Kd colour");
            mat.texture.clear();
            continue;
        }
        if (w > lim.maxTextureSide || h > lim.maxTextureSide || bytes + uint64_t(w) * h * 4 > lim.maxTextureBytes) {
            if (warnings) warnings->push_back(file + " is " + std::to_string(w) + " x " + std::to_string(h) + " pixels, beyond the texture limit (" +
                                              std::to_string(lim.maxTextureSide) + " per side, 512 MB in all); material " + mat.name + " uses its Kd colour");
            mat.texture.clear();
            continue;
        }
        ColourImage t;
        if (!decodeImageFile(mat.texture, t, &why)) {
            if (warnings) warnings->push_back(file + ": cannot decode (" + why + ")");
            mat.texture.clear();
            continue;
        }
        bytes += t.rgba.size();
        seen[mat.texture] = out.size();
        out.push_back(std::move(t));
    }
    return out;
}

} // namespace knit
