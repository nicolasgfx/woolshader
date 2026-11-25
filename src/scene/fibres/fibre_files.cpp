// Saved fibres: the adapter between the fibre files and the shared result (fibre_files.h).
#include "scene/fibres/fibre_files.h"
#include "scene/fibres/part_fibre_job.h"
#include "scene/stage_params.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <utility>
#include <vector>

#include "fibres/fibre_checks.h"
#include "io/files.h"
#include "scene/scene.h"
#include "settings/settings_mapping.h"
#include "fibres/fibre_budget.h"
#include "scene/fibres/scene_fibres.h"

namespace knit {

// --- the file format, fibre-cache v6

namespace {

void putFibre(ByteWriter& w, const FibreSet& fs, const Fibre& f) {
    const std::span<const Vec3f> pts = fs.pointsOf(f);
    w.i64(f.yarnId); w.i64(f.parentYarn); w.i64(static_cast<int>(f.population)); w.f64(f.s); w.f64(f.theta); w.vec3(toVec3(f.anchor)); w.vec3(toVec3(f.normal));
    w.i64(int64_t(pts.size()));
    for (const Vec3f& p : pts) w.vec3(toVec3(p));
    w.i64(f.truncated); w.i64(1); w.i64(0); w.vec3(toVec3(pts[0]));  // a buried root, not collinear, at pts[0]
    w.f64(f.targetLength); w.f64(f.radius0); w.f64(f.radius1);
    w.i64(1); w.f64(f.emergeL); w.f64(f.attach); w.f64(f.attachArc); w.f64(f.key);
    const float fl[8] = {f.psi, f.lift, f.smoothLength, f.emergeAngle, f.helixActual, f.maxDev, f.minTurn, f.maxTurn};
    w.bytes(fl, sizeof fl);
    w.i64(f.candidate);
}
void getFibre(ByteReader& q, FibreSet& fs) {
    Fibre f;
    f.yarnId = int(q.i64()); f.parentYarn = int(q.i64()); f.population = static_cast<FibrePopulation>(q.i64()); f.s = q.f64(); f.theta = float(q.f64());
    f.anchor = toVec3f(q.vec3()); f.normal = toVec3f(q.vec3());
    std::vector<Vec3f> pts(q.count(24));
    for (Vec3f& p : pts) p = toVec3f(q.vec3());
    f.truncated = q.i64() != 0;
    q.i64(); q.i64(); q.vec3();  // retired: has root, collinear root, root
    f.targetLength = float(q.f64()); f.radius0 = float(q.f64()); f.radius1 = float(q.f64());
    q.i64();  // retired: polyline
    f.emergeL = float(q.f64()); f.attach = float(q.f64()); f.attachArc = float(q.f64()); f.key = q.f64();
    float fl[8];
    q.bytes(fl, sizeof fl);
    f.psi = fl[0]; f.lift = fl[1]; f.smoothLength = fl[2]; f.emergeAngle = fl[3]; f.helixActual = fl[4]; f.maxDev = fl[5];
    f.minTurn = fl[6]; f.maxTurn = fl[7];
    f.candidate = int(q.i64());
    fs.addFibre(f, std::span<const Vec3f>(pts));
}
void putReport(ByteWriter& w, const FibreReport& r) {
    for (int FibreReport::*const counter : fibreReportCounters()) w.i64(r.*counter);
    w.f64(r.totalLength); w.f64(r.ms);
    w.i64(int64_t(r.pop.size()));
    for (const FibreReport& p : r.pop) putReport(w, p);
}
FibreReport getReport(ByteReader& q) {
    FibreReport r;
    for (int FibreReport::*const counter : fibreReportCounters()) r.*counter = int(q.i64());
    r.totalLength = q.f64(); r.ms = q.f64();
    r.pop.resize(q.count(8));
    for (FibreReport& p : r.pop) p = getReport(q);
    return r;
}


const char kFibreMagic[] = "fibre-cache v6\n";

// What a file holds beyond the fibres and is not read back: the candidate source of a budget run
// and the results of the fibre checks. A file written here has them empty.
void writeEmptyCandidateSource(ByteWriter& w) {
    w.i64(0); w.i64(0);    // candidates, covered
    w.i64(0); w.i64(0); w.i64(0); w.i64(0);  // has source, placement, seed, keyed
    w.ints({}); w.ints({});  // per strand: yarn id, count
}
void skipCandidateSource(ByteReader& q) {
    for (int k = 0; k < 6; ++k) q.i64();
    q.ints();
    q.ints();
}
void writeEmptyCheckResults(ByteWriter& w) {
    w.i64(0);                                     // the gap rays did not run
    for (int k = 0; k < 5; ++k) w.i64(0);         // open fraction
    w.f64(0); w.f64(0);
    w.doubles({});
    for (int k = 0; k < 8; ++k) w.i64(0);         // attachment
    w.f64(0); w.ints({}); w.f64(0); w.f64(0);
    w.i64(0);                                     // no reference check
}
void skipCheckResults(ByteReader& q) {
    q.i64();
    for (int k = 0; k < 5; ++k) q.i64();
    q.f64(); q.f64();
    q.doubles();
    for (int k = 0; k < 8; ++k) q.i64();
    q.f64(); q.ints(); q.f64(); q.f64();
    const int64_t references = q.i64();
    if (references < 0 || references > 2) q.refuse();
    for (int k = 0; q.ok() && k < references; ++k) {
        for (int v = 0; v < 5; ++v) q.i64();
        q.f64(); q.f64(); q.ints(); q.ints(); q.f64(); q.f64();
    }
}

} // namespace

bool writeFibreCache(const std::string& path, const FibreFileIdentity& identity, const FibreSet& fs, std::string* err) {
    ByteWriter w;
    w.bytes(kFibreMagic, sizeof kFibreMagic - 1);
    w.i64(0);  // the attachment check did not run
    w.string64(identity.key);
    w.f64(double(identity.keyBound));
    writeEmptyCandidateSource(w);
    w.i64(int64_t(fs.fibres.size()));
    for (const Fibre& f : fs.fibres) putFibre(w, fs, f);
    putReport(w, fs.rep);
    writeEmptyCheckResults(w);
    w.i64(static_cast<int64_t>(identity.fibreDigest));
    w.i64(0);  // check digest
    w.f64(0);  // generation time
    w.bytes(kFibreMagic, sizeof kFibreMagic - 1);
    return replaceFile(path, w.written(), err);
}

std::optional<FibreCacheEntry> readFibreCache(const std::string& path, std::string* err) {
    std::string b;
    if (!readWholeFile(path, b, err)) return std::nullopt;
    const size_t m = sizeof kFibreMagic - 1;
    FibreCacheEntry e;
    if (b.size() < 2 * m || b.compare(0, m, kFibreMagic) != 0) {
        (void)fail(err, path + ": not a saved fibre file of this version (fibre-cache v6; older ones are no longer read: grow the fibres again)");
        return std::nullopt;
    }
    ByteReader q(b, m);
    q.i64();  // whether the attachment check ran
    e.identity.key = q.string64();
    e.identity.keyBound = float(q.f64());
    skipCandidateSource(q);
    const size_t fibres = q.count(64);
    for (size_t k = 0; k < fibres; ++k) getFibre(q, e.fs);
    e.fs.rep = getReport(q);
    skipCheckResults(q);
    e.identity.fibreDigest = static_cast<uint64_t>(q.i64());
    q.i64();  // check digest
    q.f64();  // generation time
    if (!q.ok() || q.at() + m != b.size() || b.compare(q.at(), m, kFibreMagic) != 0) {
        if (err) *err = path + ": truncated or corrupt saved fibre file";
        return std::nullopt;
    }
    return e;
}

namespace {
// Two key values agree: equal text, or numbers equal to rounding. Derived values such as the
// serialized margin differ in the last bits between the former route and this one; a changed setting
// differs by far more.
bool sameKeyValue(const std::string& a, const std::string& b) {
    if (a == b) return true;
    char* endA = nullptr;
    char* endB = nullptr;
    const double x = std::strtod(a.c_str(), &endA), y = std::strtod(b.c_str(), &endB);
    return *endA == 0 && *endB == 0 && std::fabs(x - y) <= 1e-12 * std::max(std::fabs(x), std::fabs(y));
}
} // namespace

std::string staleFibreKeyField(const std::string& have, const std::string& want) {
    std::istringstream wanted(want);
    std::string line;
    while (std::getline(wanted, line)) {
        const size_t eq = line.find(" = ");
        if (eq == std::string::npos) continue;
        const std::string name = line.substr(0, eq), label = name + " = ";
        // Searched anywhere, not per line: the former writer's fixed buffer cut the line before
        // "fibre params". A file older than the placement line grew at the legacy placement.
        const size_t at = have.find(label);
        if (at == std::string::npos) {
            if (name == "fibre placement" && line.substr(eq + 3) == "legacy") continue;
            return name;
        }
        const size_t from = at + label.size();
        std::istringstream stored(have.substr(from, have.find('\n', from) - from)), current(line.substr(eq + 3));
        std::string storedValue, currentValue;
        // `want` may stop early: the former export tube line continues with values that no longer exist.
        while (current >> currentValue)
            if (!(stored >> storedValue) || !sameKeyValue(storedValue, currentValue)) return name;
    }
    return "";
}

// --- a scene's fibres

namespace {

std::string fibreFilePath(const std::string& dir, const std::string& partName) { return dir + "/" + partName + ".fibrecache"; }

// A wool result of the fibre step: what the files can hold.
bool checkWoolOnly(const SceneSettings& s, std::string* err) {
    if (knitWorkflow(s)) return true;
    if (err) *err = "saved fibres hold wool fibres only (no fur)";
    return false;
}

} // namespace

bool saveSceneFibres(const std::string& dir, const Scene& t, const SceneSettings& s, const SceneFibres& f, std::string* err) {
    if (!checkWoolOnly(s, err)) return false;
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);
    FibreParams fpBase;
    if (!fibreParamsFromSettings(s, fpBase, err)) return false;
    for (const PartFibres& pf : f.parts) {
        const std::optional<std::string> key = partFibreKey(t.parts[pf.part], pf.part, s, fpBase, err);
        if (!key) return false;
        const FibreFileIdentity identity{*key, pf.coveredBelow, fibreSnapshotChunkDigest(pf.fs, ExecPolicy{true})};
        if (!writeFibreCache(fibreFilePath(dir, pf.name), identity, pf.fs, err)) return false;
    }
    return true;
}

bool loadSceneFibres(const std::string& dir, const Scene& t, const SceneSettings& s, SceneFibres& out, std::string* err) {
    if (!checkWoolOnly(s, err)) return false;
    FibreParams fpBase;
    if (!fibreParamsFromSettings(s, fpBase, err)) return false;
    out = SceneFibres{};
    long long fibreCount = 0;
    for (size_t i = 0; i < t.parts.size(); ++i) {
        const ScenePart& r = t.parts[i];
        if (!r.hasStitches()) continue;
        const std::string path = fibreFilePath(dir, r.name);
        auto refuse = [&](const std::string& why) {
            if (err) *err = path + ": saved fibres refused (" + why + ")";
            return false;
        };
        std::optional<FibreCacheEntry> saved = readFibreCache(path, err);
        if (!saved) return false;
        const FibreFileIdentity& identity = saved->identity;
        const std::optional<std::string> want = partFibreKey(r, static_cast<int>(i), s, fpBase, err);
        if (!want) return false;
        if (const std::string field = staleFibreKeyField(identity.key, *want); !field.empty()) return refuse("stale: " + field);
        if (fibreSnapshotChunkDigest(saved->fs, ExecPolicy{true}) != identity.fibreDigest) return refuse("the fibres differ from their stored digest");
        fibreCount += static_cast<long long>(saved->fs.fibres.size());
        PartFibres pf{static_cast<int>(i), r.name, std::move(saved->fs), {}};
        pf.coveredBelow = identity.keyBound;
        out.budgetMode = out.budgetMode || identity.keyBound < 2.0f;
        out.parts.push_back(std::move(pf));
    }
    const std::string budgets = applyObjectBudgets(t, s.segmentBudget, out.parts);
    for (const PartFibres& pf : out.parts)
        if (pf.keyBelow > pf.coveredBelow) {
            char b[240];
            std::snprintf(b, sizeof b, "the saved fibres of %s cover key < %.6f: fibres.segment_budget %lld needs more of them, so they must be grown",
                          pf.name.c_str(), double(pf.coveredBelow), s.segmentBudget);
            if (err) *err = b;
            return false;
        }
    char b[200];
    std::snprintf(b, sizeof b, "fibres: %lld from %zu saved parts (%s; digests verified)\n", fibreCount, out.parts.size(), dir.c_str());
    out.log = b + budgets;
    return true;
}

} // namespace knit
