// The "knit-markings" text format, v5 (v1 to v4 read and migrated) (marking_file.h).
#include "io/marking_file.h"

#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <sstream>

#include "io/files.h"

namespace knit {

namespace {

const char kHeader[] = "knit-markings v5";
const char kHeaderV4[] = "knit-markings v4";   // read and migrated: v5 only adds the markings' knit_component (v4 means 0)
const char kHeaderV3[] = "knit-markings v3";   // read and migrated: v4 only adds the proposals' setup and basis
const char kHeaderV2[] = "knit-markings v2";   // read and migrated: v3 only adds the set's preparation (v2 means facets)
const char kHeaderV1[] = "knit-markings v1";   // read and migrated: v2 only adds per-marking settings

// Free text stays on one line.
std::string oneLine(std::string s) {
    for (char& c : s)
        if (c == '\n' || c == '\r') c = ' ';
    return s;
}

std::string trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r");
    if (a == std::string::npos) return "";
    const size_t b = s.find_last_not_of(" \t\r");
    return s.substr(a, b - a + 1);
}

bool parseInt(const std::string& text, int& out) {
    char* end = nullptr;
    const long v = std::strtol(text.c_str(), &end, 10);
    if (end == text.c_str() || *end != '\0') return false;
    out = static_cast<int>(v);
    return true;
}

std::string lineError(int line, const std::string& what) { return "markings file, line " + std::to_string(line) + ": " + what; }

} // namespace

std::string writeMarkings(const MarkingSet& s) {
    std::ostringstream o;
    char b[160];
    o << kHeader << "\n";
    o << "source = " << oneLine(s.source) << "\n";
    std::snprintf(b, sizeof b, "digest = %016" PRIx64 "\n", s.digest);
    o << b;
    o << "next_id = " << s.nextId << "\n";
    o << "preparation = " << oneLine(s.preparation) << "\n";
    for (const SurfaceMarking& m : s.markings) {
        o << "marking " << m.id << "\n";
        o << "meaning = " << toString(m.meaning) << "\n";
        o << "shape = " << toString(m.shape) << "\n";
        o << "origin = " << toString(m.origin) << "\n";
        o << "proposal = " << toString(m.proposal) << "\n";
        o << "proposal_id = " << m.proposalId << "\n";
        o << "group = " << m.alternativeGroup << "\n";
        if (m.setup > 0) o << "setup = " << m.setup << "\n";
        if (m.knitComponent > 0) o << "knit_component = " << m.knitComponent << "\n";
        if (m.basis != 0) {
            std::snprintf(b, sizeof b, "basis = %016" PRIx64 "\n", m.basis);
            o << b;
        }
        if (!m.note.empty()) o << "note = " << oneLine(m.note) << "\n";
        if (!m.reason.empty()) o << "reason = " << oneLine(m.reason) << "\n";
        if (!m.unsuitable.empty()) o << "unsuitable = " << oneLine(m.unsuitable) << "\n";
        for (const auto& [k, v] : m.settings) o << "setting " << oneLine(k) << " = " << oneLine(v) << "\n";
        for (const SurfacePoint& p : m.points) {
            std::snprintf(b, sizeof b, "point %d %.17g %.17g %.17g\n", p.tri, p.bary[0], p.bary[1], p.bary[2]);
            o << b;
        }
        o << "end\n";
    }
    return o.str();
}

bool readMarkings(const std::string& text, MarkingSet& out, std::string* err) {
    out = MarkingSet{};
    std::istringstream in(text);
    std::string line;
    int n = 0;
    if (!std::getline(in, line)) return fail(err, "not a markings file (it is empty)");
    if (trim(line) == kHeader) out.version = 5;
    else if (trim(line) == kHeaderV4) out.version = 4;
    else if (trim(line) == kHeaderV3) out.version = 3;
    else if (trim(line) == kHeaderV2) out.version = 2;
    else if (trim(line) == kHeaderV1) out.version = 1;
    else if (trim(line).rfind("knit-markings v", 0) == 0) return fail(err, "markings file version \"" + trim(line) + "\" is newer than this program reads (v1 to v5)");
    else return fail(err, "not a markings file (the first line must be \"" + std::string(kHeader) + "\")");
    ++n;
    SurfaceMarking* m = nullptr;
    while (std::getline(in, line)) {
        ++n;
        const std::string t = trim(line);
        if (t.empty() || t[0] == '#') continue;
        if (t.rfind("marking ", 0) == 0) {
            if (m) return fail(err, lineError(n, "a marking starts before the previous one ends"));
            SurfaceMarking fresh;
            if (!parseInt(trim(t.substr(8)), fresh.id) || fresh.id <= 0) return fail(err, lineError(n, "a marking needs a positive id"));
            for (const SurfaceMarking& other : out.markings)
                if (other.id == fresh.id) return fail(err, lineError(n, "marking id " + std::to_string(fresh.id) + " is used twice"));
            out.markings.push_back(fresh);
            m = &out.markings.back();
            continue;
        }
        if (t == "end") {
            if (!m) return fail(err, lineError(n, "\"end\" without a marking"));
            m = nullptr;
            continue;
        }
        if (t.rfind("point ", 0) == 0) {
            if (!m) return fail(err, lineError(n, "a point outside a marking"));
            std::istringstream ps(t.substr(6));
            SurfacePoint p;
            std::string rest;
            if (!(ps >> p.tri >> p.bary[0] >> p.bary[1] >> p.bary[2]) || (ps >> rest) || p.tri < 0)
                return fail(err, lineError(n, "a point is a triangle index and three weights"));
            m->points.push_back(p);
            continue;
        }
        if (t.rfind("setting ", 0) == 0) {
            if (!m) return fail(err, lineError(n, "a setting outside a marking"));
            if (out.version < 2) return fail(err, lineError(n, "settings need a v2 file"));
            const size_t e = t.find('=');
            const std::string key = e == std::string::npos ? "" : trim(t.substr(8, e - 8));
            if (key.empty() || key.find(' ') != std::string::npos) return fail(err, lineError(n, "a setting is \"setting <key> = <value>\""));
            m->settings.push_back({key, trim(t.substr(e + 1))});
            continue;
        }
        const size_t eq = t.find('=');
        if (eq == std::string::npos) return fail(err, lineError(n, "expected \"key = value\""));
        const std::string key = trim(t.substr(0, eq)), value = trim(t.substr(eq + 1));
        bool ok = true;
        if (!m) {
            if (key == "source") out.source = value;
            else if (key == "digest") {
                char* end = nullptr;
                out.digest = std::strtoull(value.c_str(), &end, 16);
                ok = end != value.c_str() && *end == '\0';
            }
            else if (key == "next_id") ok = parseInt(value, out.nextId);
            else if (key == "preparation") {
                if (out.version < 3) return fail(err, lineError(n, "the preparation needs a v3 file"));
                if (value != "facets" && value != "loop_limit") return fail(err, lineError(n, "preparation \"" + value + "\" is not facets or loop_limit"));
                out.preparation = value;
            }
            else return fail(err, lineError(n, "unknown key \"" + key + "\""));
        } else {
            if (key == "meaning") ok = parseMeaning(value, m->meaning);
            else if (key == "shape") ok = parseShape(value, m->shape);
            else if (key == "origin") ok = parseOrigin(value, m->origin);
            else if (key == "proposal") ok = parseProposalState(value, m->proposal);
            else if (key == "proposal_id") ok = parseInt(value, m->proposalId);
            else if (key == "group") ok = parseInt(value, m->alternativeGroup);
            else if (key == "setup" || key == "basis") {
                if (out.version < 4) return fail(err, lineError(n, "\"" + key + "\" needs a v4 file"));
                if (key == "setup") ok = parseInt(value, m->setup) && m->setup >= 0;
                else {
                    char* end = nullptr;
                    m->basis = std::strtoull(value.c_str(), &end, 16);
                    ok = end != value.c_str() && *end == '\0';
                }
            }
            else if (key == "knit_component") {
                if (out.version < 5) return fail(err, lineError(n, "\"knit_component\" needs a v5 file"));
                ok = parseInt(value, m->knitComponent) && m->knitComponent >= 0;
            }
            else if (key == "note") m->note = value;
            else if (key == "reason") m->reason = value;
            else if (key == "unsuitable") m->unsuitable = value;
            else return fail(err, lineError(n, "unknown key \"" + key + "\" in a marking"));
        }
        if (!ok) return fail(err, lineError(n, "\"" + value + "\" is not a value of " + key));
    }
    if (m) return fail(err, "markings file: the last marking has no \"end\"");
    for (const SurfaceMarking& x : out.markings)
        if (x.id >= out.nextId) out.nextId = x.id + 1;  // never reuse an id, even when next_id was written by hand
    return true;
}

bool saveMarkingFile(const std::string& path, const MarkingSet& s, std::string* err) { return replaceFile(path, writeMarkings(s), err); }

bool loadMarkingFile(const std::string& path, MarkingSet& out, std::string* err) {
    std::string text;
    if (!readWholeFile(path, text, err)) return false;
    return readMarkings(text, out, err);
}

} // namespace knit
