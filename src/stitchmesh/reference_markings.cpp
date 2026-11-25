#include "stitchmesh/reference_markings.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "model/hash.h"
#include "model/mesh_query.h"

namespace knit::refmark {

namespace {

struct TypeName {
    MarkingType type;
    const char* name;
};
const TypeName kTypes[] = {{MarkingType::StartPoint, "start_point"},       {MarkingType::FinishBoundary, "finish_boundary"},
                           {MarkingType::GrowthDirection, "growth_direction"}, {MarkingType::RowStartPath, "row_start_path"},
                           {MarkingType::Cut, "cut"},                        {MarkingType::AttachmentRegion, "attachment_region"}};

bool typeFrom(const std::string& s, MarkingType& t) {
    for (const TypeName& n : kTypes)
        if (s == n.name) {
            t = n.type;
            return true;
        }
    return false;
}

bool provenanceFrom(const std::string& s, Provenance& p) {
    if (s == "manual") p = Provenance::Manual;
    else if (s == "automatic_proposal") p = Provenance::AutomaticProposal;
    else if (s == "user_edited") p = Provenance::UserEdited;
    else return false;
    return true;
}

// "key=value" tokens of a line after its first word.
std::map<std::string, std::string> tokens(std::istringstream& in) {
    std::map<std::string, std::string> out;
    std::string tok;
    while (in >> tok) {
        const size_t eq = tok.find('=');
        if (eq != std::string::npos) out[tok.substr(0, eq)] = tok.substr(eq + 1);
    }
    return out;
}

bool fail(std::string* err, const std::string& what, int line) {
    if (err) *err = "markings line " + std::to_string(line) + ": " + what;
    return false;
}

// Barycentric coordinates of q in triangle (a, b, c) projected onto its plane, and the distance of q from that plane.
void baryOf(const Vec3& q, const Vec3& a, const Vec3& b, const Vec3& c, double out[3], double* planeDistance) {
    const Vec3 e0 = b - a, e1 = c - a, d = q - a;
    const double d00 = dot(e0, e0), d01 = dot(e0, e1), d11 = dot(e1, e1), d20 = dot(d, e0), d21 = dot(d, e1);
    const double den = d00 * d11 - d01 * d01;
    if (!(std::fabs(den) > 0)) {
        out[0] = 1;
        out[1] = out[2] = 0;
        if (planeDistance) *planeDistance = 1e300;
        return;
    }
    const double v = (d11 * d20 - d01 * d21) / den, w = (d00 * d21 - d01 * d20) / den;
    out[0] = 1 - v - w;
    out[1] = v;
    out[2] = w;
    if (planeDistance) {
        const Vec3 n = cross(e0, e1);
        *planeDistance = std::fabs(dot(d, n)) / norm(n);
    }
}

} // namespace

const char* toString(MarkingType t) {
    for (const TypeName& n : kTypes)
        if (n.type == t) return n.name;
    return "?";
}

const char* toString(Provenance p) {
    switch (p) {
    case Provenance::Manual: return "manual";
    case Provenance::AutomaticProposal: return "automatic_proposal";
    case Provenance::UserEdited: return "user_edited";
    }
    return "?";
}

std::string Marking::setting(const std::string& key, const std::string& fallback) const {
    for (const auto& [k, v] : settings)
        if (k == key) return v;
    return fallback;
}

const Marking* MarkingSet::find(MarkingType t) const {
    for (const Marking& m : markings)
        if (m.type == t) return &m;
    return nullptr;
}

const Marking* MarkingSet::byId(const std::string& id) const {
    for (const Marking& m : markings)
        if (m.id == id) return &m;
    return nullptr;
}

// The format, one statement per line:
//   stockinette-markings 1
//   component name=<label> node=<node> material=<material> piece=<k> source_digest=<hex> import_identity=<hex>
//   reference digest=<hex> vertices=<n> triangles=<n>
//   frame <free text to the end of the line>
//   marking id=<id> type=<type> provenance=<provenance> closed=0|1
//   setting <key>=<value>
//   point <triangle> <b0> <b1> <b2>
//   triangles <t> <t> ...
//   note <free text to the end of the line>
//   end
// Labels and settings carry no spaces (the writer replaces them with '_').
bool writeMarkings(const std::string& path, const MarkingSet& s, std::string* err) {
    std::ofstream f(path, std::ios::binary);
    if (!f) {
        if (err) *err = "cannot write " + path;
        return false;
    }
    auto clean = [](std::string t) {
        for (char& c : t)
            if (c == ' ' || c == '\t' || c == '\n' || c == '\r') c = '_';
        return t;
    };
    char b[256];
    f << "stockinette-markings " << s.version << "\n";
    f << "component name=" << clean(s.component.name) << " node=" << clean(s.component.node) << " material=" << clean(s.component.material)
      << " piece=" << s.component.piece << " source_digest=" << s.component.sourceDigest << " import_identity=" << s.component.importIdentity << "\n";
    f << "reference digest=" << s.reference.digest << " vertices=" << s.reference.vertices << " triangles=" << s.reference.triangles << "\n";
    f << "frame " << s.reference.frame << "\n";
    for (const Marking& m : s.markings) {
        f << "marking id=" << clean(m.id) << " type=" << toString(m.type) << " provenance=" << toString(m.provenance) << " closed=" << (m.closed ? 1 : 0) << "\n";
        for (const auto& [k, v] : m.settings) f << "setting " << clean(k) << "=" << clean(v) << "\n";
        for (const SurfacePoint& p : m.points) {
            std::snprintf(b, sizeof b, "point %d %.17g %.17g %.17g\n", p.tri, p.bary[0], p.bary[1], p.bary[2]);
            f << b;
        }
        if (!m.triangles.empty()) {
            f << "triangles";
            for (int t : m.triangles) f << " " << t;
            f << "\n";
        }
        if (!m.note.empty()) f << "note " << m.note << "\n";
        f << "end\n";
    }
    return static_cast<bool>(f);
}

bool readMarkings(const std::string& path, MarkingSet& s, std::string* err) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        if (err) *err = "cannot read " + path;
        return false;
    }
    s = MarkingSet{};
    std::string line;
    int no = 0;
    Marking* open = nullptr;
    bool header = false;
    while (std::getline(f, line)) {
        ++no;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::istringstream in(line);
        std::string word;
        in >> word;
        if (!header) {
            if (word != "stockinette-markings" || !(in >> s.version)) return fail(err, "the first statement must be 'stockinette-markings <version>'", no);
            if (s.version != 1) return fail(err, "version " + std::to_string(s.version) + " is not supported (this reader knows version 1)", no);
            header = true;
            continue;
        }
        if (word == "component") {
            const auto t = tokens(in);
            const auto get = [&](const char* k) { const auto it = t.find(k); return it == t.end() ? std::string() : it->second; };
            s.component.name = get("name");
            s.component.node = get("node");
            s.component.material = get("material");
            s.component.sourceDigest = get("source_digest");
            s.component.importIdentity = get("import_identity");
            s.component.piece = get("piece").empty() ? -1 : std::stoi(get("piece"));
        } else if (word == "reference") {
            const auto t = tokens(in);
            const auto it = t.find("digest");
            s.reference.digest = it == t.end() ? "" : it->second;
            s.reference.vertices = t.count("vertices") ? std::stoi(t.at("vertices")) : 0;
            s.reference.triangles = t.count("triangles") ? std::stoi(t.at("triangles")) : 0;
        } else if (word == "frame") {
            std::getline(in >> std::ws, s.reference.frame);
        } else if (word == "marking") {
            if (open) return fail(err, "a marking starts before the previous one ended", no);
            const auto t = tokens(in);
            Marking m;
            m.id = t.count("id") ? t.at("id") : "";
            if (!t.count("type") || !typeFrom(t.at("type"), m.type)) return fail(err, "unknown marking type", no);
            if (!t.count("provenance") || !provenanceFrom(t.at("provenance"), m.provenance))
                return fail(err, "provenance must be manual, automatic_proposal or user_edited", no);
            m.closed = t.count("closed") && t.at("closed") == "1";
            s.markings.push_back(m);
            open = &s.markings.back();
        } else if (word == "setting" || word == "point" || word == "triangles" || word == "note") {
            if (!open) return fail(err, "'" + word + "' outside a marking", no);
            if (word == "setting") {
                std::string kv;
                in >> kv;
                const size_t eq = kv.find('=');
                if (eq == std::string::npos) return fail(err, "a setting needs key=value", no);
                open->settings.push_back({kv.substr(0, eq), kv.substr(eq + 1)});
            } else if (word == "point") {
                SurfacePoint p;
                if (!(in >> p.tri >> p.bary[0] >> p.bary[1] >> p.bary[2])) return fail(err, "a point needs a triangle and three barycentric coordinates", no);
                open->points.push_back(p);
            } else if (word == "triangles") {
                int t;
                while (in >> t) open->triangles.push_back(t);
            } else {
                std::getline(in >> std::ws, open->note);
            }
        } else if (word == "end") {
            if (!open) return fail(err, "'end' without a marking", no);
            open = nullptr;
        } else {
            return fail(err, "unknown statement '" + word + "'", no);
        }
    }
    if (!header) return fail(err, "empty file", no);
    if (open) return fail(err, "the last marking has no 'end'", no);
    return true;
}

std::string referenceDigest(const TriMesh& m) {
    uint64_t h = kFnvBasis;
    char b[96];
    for (const Vec3& p : m.pos) {
        const int n = std::snprintf(b, sizeof b, "%.6f %.6f %.6f;", p.x, p.y, p.z);
        h = fnv1a64Bytes(b, static_cast<size_t>(n), h);
    }
    for (const auto& t : m.tri) {
        const int n = std::snprintf(b, sizeof b, "%d %d %d;", t[0], t[1], t[2]);
        h = fnv1a64Bytes(b, static_cast<size_t>(n), h);
    }
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(h));
    return b;
}

Vec3 pointPosition(const TriMesh& m, const SurfacePoint& p) {
    const auto& t = m.tri[p.tri];
    return p.bary[0] * m.pos[t[0]] + p.bary[1] * m.pos[t[1]] + p.bary[2] * m.pos[t[2]];
}

int commonTriangle(const TriMesh& m, const SurfacePoint& a, const SurfacePoint& b) {
    if (a.tri == b.tri) return a.tri;
    const Vec3 pa = pointPosition(m, a), pb = pointPosition(m, b);
    const double scale = norm(m.pos[m.tri[a.tri][1]] - m.pos[m.tri[a.tri][0]]) + 1e-12;
    // triangles sharing a vertex with a's triangle: the ones a point on a's edges or corners can also lie in
    std::set<int> corners(m.tri[a.tri].begin(), m.tri[a.tri].end());
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t) {
        bool near = false;
        for (int v : m.tri[t]) near = near || corners.count(v) > 0;
        if (!near) continue;
        bool both = true;
        for (const Vec3& q : {pa, pb}) {
            double w[3], d = 0;
            baryOf(q, m.pos[m.tri[t][0]], m.pos[m.tri[t][1]], m.pos[m.tri[t][2]], w, &d);
            if (d > 1e-7 * scale || w[0] < -1e-7 || w[1] < -1e-7 || w[2] < -1e-7) both = false;
        }
        if (both) return t;
    }
    return -1;
}

namespace {

// q inside triangle t (on its plane, within tolerance), with its barycentric coordinates
bool inTriangle(const TriMesh& m, int t, const Vec3& q, double scale, double out[3]) {
    double d = 0;
    baryOf(q, m.pos[m.tri[t][0]], m.pos[m.tri[t][1]], m.pos[m.tri[t][2]], out, &d);
    return d <= 1e-7 * scale && out[0] >= -1e-9 && out[1] >= -1e-9 && out[2] >= -1e-9;
}

} // namespace

bool densifyPath(const TriMesh& m, std::vector<SurfacePoint>& path, int* inserted) {
    if (inserted) *inserted = 0;
    std::map<std::pair<int, int>, std::vector<int>> edgeTris;
    std::vector<std::vector<int>> trisOfVertex(m.pos.size());
    for (int t = 0; t < static_cast<int>(m.tri.size()); ++t)
        for (int k = 0; k < 3; ++k) {
            edgeTris[{std::min(m.tri[t][k], m.tri[t][(k + 1) % 3]), std::max(m.tri[t][k], m.tri[t][(k + 1) % 3])}].push_back(t);
            trisOfVertex[m.tri[t][k]].push_back(t);
        }
    std::vector<SurfacePoint> out;
    for (size_t i = 0; i + 1 < path.size(); ++i) {
        out.push_back(path[i]);
        if (commonTriangle(m, path[i], path[i + 1]) >= 0) continue;
        const Vec3 q = pointPosition(m, path[i + 1]);
        Vec3 c = pointPosition(m, path[i]);
        const double scale = norm(q - c) + 1e-12;
        // the triangle the segment enters: one containing c and a point a little way towards q
        int tri = -1;
        {
            const Vec3 ahead = c + (q - c) * 1e-6;
            std::set<int> around;
            for (int v : m.tri[path[i].tri])
                for (int t : trisOfVertex[v]) around.insert(t);
            double w[3];
            for (int t : around)
                if (inTriangle(m, t, c, scale, w) && inTriangle(m, t, ahead, scale, w)) {
                    tri = t;
                    break;
                }
        }
        if (tri < 0) return false;
        bool reached = false;
        for (int guard = 0; guard < 100000 && !reached; ++guard) {
            double wq[3];
            if (inTriangle(m, tri, q, scale, wq)) {
                reached = true;
                break;
            }
            // leave tri along c -> q: the barycentric coordinate that reaches 0 first (beyond the current point)
            const auto& t = m.tri[tri];
            double bc[3], bq[3];
            baryOf(c, m.pos[t[0]], m.pos[t[1]], m.pos[t[2]], bc, nullptr);
            baryOf(q, m.pos[t[0]], m.pos[t[1]], m.pos[t[2]], bq, nullptr);
            double best = 1e300;
            int edge = -1;
            for (int k = 0; k < 3; ++k) {
                const double d = bc[k] - bq[k];
                if (d <= 1e-15) continue;   // this coordinate does not decrease
                const double s = std::max(0.0, bc[k]) / d;
                if (s < best) {
                    best = s;
                    edge = k;
                }
            }
            if (edge < 0 || best >= 1) return false;
            SurfacePoint p;
            p.tri = tri;
            for (int k = 0; k < 3; ++k) p.bary[k] = bc[k] + best * (bq[k] - bc[k]);
            p.bary[edge] = 0;
            for (int k = 0; k < 3; ++k) p.bary[k] = std::max(0.0, p.bary[k]);
            const double sum = p.bary[0] + p.bary[1] + p.bary[2];
            for (int k = 0; k < 3; ++k) p.bary[k] /= sum;
            out.push_back(p);
            if (inserted) ++*inserted;
            c = pointPosition(m, p);
            const int a = t[(edge + 1) % 3], b = t[(edge + 2) % 3];
            const std::vector<int>& two = edgeTris[{std::min(a, b), std::max(a, b)}];
            if (two.size() != 2) return false;
            tri = two[0] == tri ? two[1] : two[0];
        }
        if (!reached) return false;
    }
    out.push_back(path.back());
    path = out;
    return true;
}

std::string MarkingCheck::text() const {
    std::string t = std::string("markings ") + toString(status) + ": " + std::to_string(points) + " surface points; provenance: " + std::to_string(proposals) +
                    " automatic proposals, " + std::to_string(edited) + " user-edited, " + std::to_string(manual) + " manual\n";
    for (const std::string& r : refusals) t += "  REFUSED: " + r + "\n";
    for (const std::string& o : open) t += "  open: " + o + "\n";
    return t;
}

MarkingCheck validateMarkings(const MarkingSet& s, const TriMesh& reference, const ComponentIdentity& expected) {
    MarkingCheck c;
    auto refuse = [&](const std::string& why) { c.refusals.push_back(why); };
    if (s.version != 1) refuse("version " + std::to_string(s.version) + " is not version 1");
    if (!(s.component == expected))
        refuse("the markings name component '" + s.component.node + "/" + s.component.material + "/" + std::to_string(s.component.piece) + "' of import " +
               s.component.importIdentity + " (source " + s.component.sourceDigest.substr(0, 12) + "), not the expected '" + expected.node + "/" + expected.material +
               "/" + std::to_string(expected.piece) + "' of import " + expected.importIdentity + " (source " + expected.sourceDigest.substr(0, 12) + ")");
    const std::string digest = referenceDigest(reference);
    if (s.reference.digest != digest || s.reference.vertices != static_cast<int>(reference.pos.size()) || s.reference.triangles != static_cast<int>(reference.tri.size()))
        refuse("the markings are bound to reference " + s.reference.digest + " (" + std::to_string(s.reference.vertices) + " vertices, " +
               std::to_string(s.reference.triangles) + " triangles), not to this one, " + digest + " (" + std::to_string(reference.pos.size()) + ", " +
               std::to_string(reference.tri.size()) + "): rebind them (checked) or mark them again");
    const int nT = static_cast<int>(reference.tri.size());
    std::set<std::string> ids;
    std::map<MarkingType, int> count;
    for (const Marking& m : s.markings) {
        const std::string who = "marking '" + m.id + "' (" + toString(m.type) + ")";
        if (m.id.empty() || !ids.insert(m.id).second) refuse(who + ": a missing or repeated id");
        count[m.type] += 1;
        (m.provenance == Provenance::AutomaticProposal ? c.proposals : m.provenance == Provenance::UserEdited ? c.edited : c.manual) += 1;
        bool pointsValid = true;
        for (const SurfacePoint& p : m.points) {
            ++c.points;
            const double sum = p.bary[0] + p.bary[1] + p.bary[2];
            const bool finite = std::isfinite(p.bary[0]) && std::isfinite(p.bary[1]) && std::isfinite(p.bary[2]);
            if (p.tri < 0 || p.tri >= nT || !finite || std::fabs(sum - 1) > 1e-9 || p.bary[0] < -1e-9 || p.bary[1] < -1e-9 || p.bary[2] < -1e-9) pointsValid = false;
        }
        if (!pointsValid) refuse(who + ": a point outside the reference (triangle id out of range, or barycentric coordinates not in the triangle)");
        for (int t : m.triangles)
            if (t < 0 || t >= nT) {
                refuse(who + ": a region triangle id out of range");
                break;
            }
        switch (m.type) {
        case MarkingType::StartPoint:
        case MarkingType::FinishBoundary:
            if (m.points.size() != 1) refuse(who + ": needs exactly one point");
            break;
        case MarkingType::RowStartPath:
        case MarkingType::Cut:
            if (m.points.size() < 2) refuse(who + ": a path needs at least two points");
            else if (pointsValid)
                for (size_t i = 0; i + 1 < m.points.size(); ++i)
                    if (commonTriangle(reference, m.points[i], m.points[i + 1]) < 0) {
                        refuse(who + ": points " + std::to_string(i) + " and " + std::to_string(i + 1) +
                               " share no triangle (a path's segments must lie on the surface; densify it across the triangles it passes)");
                        break;
                    }
            if (m.type == MarkingType::RowStartPath && m.closed) refuse(who + ": a row-start path runs from the start to the finish and cannot be closed");
            if (m.type == MarkingType::Cut) c.open.push_back(who + ": a physical cut is declared, and the mesher does not consume cuts yet (unsupported construction)");
            break;
        case MarkingType::GrowthDirection: {
            const Marking* from = s.byId(m.setting("from"));
            const Marking* to = s.byId(m.setting("to"));
            if (!from || !to || from->type != MarkingType::StartPoint || to->type != MarkingType::FinishBoundary)
                refuse(who + ": needs from=<start marking> and to=<finish marking>");
            break;
        }
        case MarkingType::AttachmentRegion:
            if (m.triangles.empty() || m.setting("other").empty()) refuse(who + ": needs triangles and other=<component>");
            if (m.setting("join", "unresolved") != "unresolved")
                refuse(who + ": an attachment region is a join candidate; join=" + m.setting("join") + " is not something a marking can establish");
            break;
        }
    }
    if (count[MarkingType::StartPoint] != 1) refuse("exactly one start_point is needed, found " + std::to_string(count[MarkingType::StartPoint]));
    if (count[MarkingType::FinishBoundary] != 1) refuse("exactly one finish_boundary is needed, found " + std::to_string(count[MarkingType::FinishBoundary]));
    if (count[MarkingType::RowStartPath] > 1) refuse("at most one row_start_path, found " + std::to_string(count[MarkingType::RowStartPath]));
    if (count[MarkingType::GrowthDirection] > 1) refuse("at most one growth_direction");
    if (c.proposals > 0) c.open.push_back(std::to_string(c.proposals) + " marking(s) are automatic proposals, not reviewed by a person");
    c.status = !c.refusals.empty() ? Status::Fail : !c.open.empty() ? Status::Unresolved : Status::Pass;
    return c;
}

bool rebindMarkings(const MarkingSet& in, const TriMesh& oldRef, const TriMesh& newRef, double tolerance, MarkingSet& out, RebindReport& report) {
    report = RebindReport{};
    out = in;
    if (in.reference.digest != referenceDigest(oldRef)) {
        report.detail = "the markings are not bound to the old reference given";
        return false;
    }
    const MeshQuery q(newRef);
    for (Marking& m : out.markings) {
        if (!m.triangles.empty()) {
            ++report.refused;
            report.detail += "region '" + m.id + "' cannot be rebound by closest points; ";
            continue;
        }
        for (SurfacePoint& p : m.points) {
            ++report.points;
            const Vec3 x = pointPosition(oldRef, p);
            const MeshQuery::Hit h = q.closest(x);
            report.maxDistance = std::max(report.maxDistance, h.dist);
            if (h.tri < 0 || h.dist > tolerance) {
                ++report.refused;
                continue;
            }
            const TrianglePoint tp = closestOnTriangle(x, newRef.pos[newRef.tri[h.tri][0]], newRef.pos[newRef.tri[h.tri][1]], newRef.pos[newRef.tri[h.tri][2]]);
            p.tri = h.tri;
            for (int k = 0; k < 3; ++k) p.bary[k] = tp.bary[k];
            const double sum = p.bary[0] + p.bary[1] + p.bary[2];
            for (int k = 0; k < 3; ++k) p.bary[k] /= sum;
        }
        if (m.type == MarkingType::RowStartPath || m.type == MarkingType::Cut) {
            int added = 0;
            if (!densifyPath(newRef, m.points, &added)) {
                ++report.refused;
                report.detail += "path '" + m.id + "' could not be traced on the new reference; ";
            }
            report.inserted += added;
        }
        m.note += (m.note.empty() ? "" : " | ") + std::string("rebound from reference ") + in.reference.digest;
    }
    out.reference.digest = referenceDigest(newRef);
    out.reference.vertices = static_cast<int>(newRef.pos.size());
    out.reference.triangles = static_cast<int>(newRef.tri.size());
    char b[160];
    std::snprintf(b, sizeof b, "%d points, max distance %.3g mm (tolerance %.3g), %d refused, %d path points inserted by densifyPath", report.points, report.maxDistance,
                  tolerance, report.refused, report.inserted);
    report.detail += b;
    return report.refused == 0;
}

} // namespace knit::refmark
