// Publisher: a generated folder is staged beside its target and appears complete or not at all (publish.h).
#include "io/publish.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <vector>

#include "io/files.h"
#include "model/cancel.h"

namespace knit {

namespace fs = std::filesystem;

namespace {

// 16 lower-case hex digits: the only tag recovery recognises.
std::string randomTag() {
    static std::mt19937_64 rng(std::random_device{}() ^
                               static_cast<uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count()));
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(rng()));
    return b;
}

bool isTag16(const std::string& s) {
    if (s.size() != 16) return false;
    for (char c : s)
        if (!((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'))) return false;
    return true;
}

// Absolute, with the PARENT resolved; the last component is kept as named (a link or junction
// as the target must be seen as one, never followed to its destination).
fs::path canonical(const std::string& p) {
    std::error_code ec;
    fs::path a = fs::absolute(fs::path(p), ec).lexically_normal();
    if (!a.has_filename()) a = a.parent_path();
    const fs::path parent = fs::weakly_canonical(a.parent_path(), ec);
    return ec ? a : parent / a.filename();
}

// Relative paths of all regular files under dir (sorted), the marker and the owner token excluded.
std::vector<std::string> listFiles(const fs::path& dir) {
    std::vector<std::string> out;
    std::error_code ec;
    for (fs::recursive_directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec))
        if (it->is_regular_file(ec)) {
            const std::string rel = fs::relative(it->path(), dir, ec).generic_string();
            if (rel != kPublishMarker && rel != kPublishOwner) out.push_back(rel);
        }
    std::sort(out.begin(), out.end());
    return out;
}

// The directory's own type, links not followed (a symlink or junction is never a directory here).
fs::file_type ownType(const fs::path& p) {
    std::error_code ec;
    return fs::symlink_status(p, ec).type();
}

// The owner token names `canonTarget`.
bool ownedBy(const fs::path& dir, const fs::path& canonTarget) {
    std::ifstream f(dir / kPublishOwner, std::ios::binary);
    std::string magic, tag, target;
    if (!std::getline(f, magic) || magic != "knit-publish-owner 1") return false;
    if (!(f >> tag) || tag != "target") return false;
    std::getline(f >> std::ws, target);
    return fs::path(target) == fs::path(canonTarget.generic_string());
}

// The marker of a publication: its key and every listed file with its size and checksum.
struct MarkerCheck {
    bool present = false, verified = false;
    std::string key, why;
    std::vector<std::string> listed;              // the marker's files, sorted
    std::vector<std::string> changed, untracked;  // listed files that differ / files not listed
};

MarkerCheck checkMarker(const fs::path& target) {
    MarkerCheck m;
    std::ifstream f(target / kPublishMarker);
    if (!f) {
        m.why = "no completion marker";
        return m;
    }
    m.present = true;
    std::string line;
    if (!std::getline(f, line) || line.rfind("key ", 0) != 0) {
        m.why = "bad marker";
        return m;
    }
    m.key = line.substr(4);
    bool bad = false;
    while (std::getline(f, line)) {
        std::istringstream ls(line);
        std::string tag, hash, rel;
        unsigned long long size = 0;
        ls >> tag >> size >> hash;
        std::getline(ls >> std::ws, rel);
        if (tag != "file") {
            m.why = "bad marker line";
            bad = true;
            break;
        }
        m.listed.push_back(rel);
        const fs::path p = target / rel;
        std::error_code ec;
        if (!fs::exists(p, ec) || fs::file_size(p, ec) != size) {
            m.changed.push_back(rel);
            continue;
        }
        char b[40];
        std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(fnv1aFile(p.string())));
        if (hash != b) m.changed.push_back(rel);
    }
    std::sort(m.listed.begin(), m.listed.end());
    for (const std::string& rel : listFiles(target))
        if (!std::binary_search(m.listed.begin(), m.listed.end(), rel)) m.untracked.push_back(rel);
    if (!bad && m.changed.empty() && m.untracked.empty()) m.verified = true;
    else if (!bad && !m.changed.empty()) m.why = m.changed.front() + ": missing or changed";
    else if (!bad) m.why = "unlisted files";
    return m;
}

std::string firstFew(const std::vector<std::string>& v) {
    std::string s;
    for (size_t i = 0; i < v.size() && i < 5; ++i) s += (i ? ", " : "") + v[i];
    if (v.size() > 5) s += ", ... (" + std::to_string(v.size()) + " in all)";
    return s;
}

bool emptyDirectory(const fs::path& p) {
    std::error_code ec;
    return fs::is_empty(p, ec) && !ec;
}

// May the existing target be replaced? (Missing: yes.)
bool replaceable(const fs::path& target, PublishPolicy policy, std::string* why) {
    auto no = [&](const std::string& w) {
        if (why) *why = w;
        return false;
    };
    const fs::file_type t = ownType(target);
    if (t == fs::file_type::not_found) return true;
    if (t != fs::file_type::directory)
        return no(target.string() + " exists and is not a plain folder (a file, link or junction): choose another folder");
    if (policy == PublishPolicy::IntoFolder) return true;  // written into: its other files stay
    if (emptyDirectory(target)) return true;
    const MarkerCheck m = checkMarker(target);
    if (policy == PublishPolicy::CacheEntry) {
        if (m.present) return true;  // the cache's own entry (possibly damaged): regenerated
        return no(target.string() + " is not a cache entry (no completion marker)");
    }
    if (!m.present)
        return no(target.string() + " already exists and is not an export of this tool: choose another folder name, or move it away");
    if (!m.untracked.empty())
        return no(target.string() + " contains files that are not part of the previous export (" + firstFew(m.untracked) +
                  "): move them out of the folder, or choose another folder name");
    if (!m.changed.empty())
        return no(target.string() + " contains export files changed since the previous export (" + firstFew(m.changed) +
                  "): move them out of the folder, or choose another folder name");
    return true;
}

// Entries of `folder` that can only have been made by a Publisher of this target: the exact
// name pattern (`prefix` + 16 hex), a plain folder (no link or junction), and this target's
// owner token.
std::vector<fs::path> ownedEntries(const fs::path& folder, const std::string& prefix, const fs::path& target, std::string* log) {
    std::vector<fs::path> out;
    std::error_code ec;
    for (fs::directory_iterator it(folder, ec), end; !ec && it != end; it.increment(ec)) {
        const std::string name = it->path().filename().string();
        if (name.rfind(prefix, 0) != 0) continue;
        if (!isTag16(name.substr(prefix.size())) || ownType(it->path()) != fs::file_type::directory || !ownedBy(it->path(), target)) {
            if (log) *log += "ignored foreign " + name + "; ";
            continue;
        }
        out.push_back(it->path());
    }
    return out;
}

std::vector<fs::path> ownedSiblings(const fs::path& target, const std::string& infix, std::string* log) {
    const fs::path parent = target.has_parent_path() ? target.parent_path() : fs::path(".");
    return ownedEntries(parent, target.filename().string() + infix, target, log);
}

// An IntoFolder commit: every staged file moves to its place in the target, the marker last. A
// file of that name goes to <staging>/.replaced first. A folder or a link of that name is left
// alone and fails the commit.
bool moveIntoFolder(const fs::path& staging, const fs::path& target, std::string* err) {
    std::vector<std::string> files = listFiles(staging);
    files.push_back(kPublishMarker);
    for (const std::string& rel : files) {
        const fs::path from = staging / rel, to = target / rel, saved = staging / kPublishReplaced / rel;
        std::error_code ec;
        const fs::file_type t = ownType(to);
        if (t != fs::file_type::not_found && t != fs::file_type::regular) {
            if (err) *err = to.string() + " is in the way (a folder or a link, not a file): move it away, or choose another folder";
            return false;
        }
        if (t == fs::file_type::regular) {
            fs::create_directories(saved.parent_path(), ec);
            fs::rename(to, saved, ec);
            if (ec) {
                if (err) *err = "cannot replace " + to.string() + ": " + ec.message() + " (is it open in another program?)";
                return false;
            }
        }
        fs::create_directories(to.parent_path(), ec);
        fs::rename(from, to, ec);
        if (ec) {
            if (err) *err = "cannot write " + to.string() + ": " + ec.message();
            return false;
        }
    }
    return true;
}

// Undoes an unfinished IntoFolder commit, whether it failed or a crash stopped it: the staged
// marker lists the export's files, and one that left the staging directory is in the target. Its
// presence means the commit did not finish, since the marker moves last. The caller removes the
// staging directory.
void rollBackIntoFolder(const fs::path& staging, const fs::path& target) {
    const MarkerCheck m = checkMarker(staging);
    if (!m.present) return;
    std::vector<std::string> files = m.listed;
    files.push_back(kPublishMarker);
    for (const std::string& rel : files) {
        const fs::path moved = target / rel, saved = staging / kPublishReplaced / rel;
        std::error_code ec;
        if (!fs::exists(staging / rel, ec)) fs::remove(moved, ec);
        if (fs::exists(saved, ec)) fs::rename(saved, moved, ec);
        // The folders the commit made (a folder that holds anything else is not removed).
        for (fs::path d = moved.parent_path(); d != target && fs::remove(d, ec); d = d.parent_path()) {}
    }
}

} // namespace

Publisher::~Publisher() { abort(); }

bool Publisher::begin(const std::string& target, std::string* err, PublishPolicy policy) {
    abort();
    const fs::path canon = canonical(target);
    std::string why;
    if (!replaceable(canon, policy, &why)) {
        if (err) *err = why;
        return false;
    }
    target_ = canon.string();
    policy_ = policy;
    created_ = false;
    std::error_code ec;
    if (policy == PublishPolicy::IntoFolder) {
        created_ = fs::create_directories(canon, ec);
        if (ec) {
            if (err) *err = "cannot create " + target_ + ": " + ec.message();
            return false;
        }
        staging_ = (canon / (".staging-" + randomTag())).string();
    } else {
        staging_ = target_ + ".staging-" + randomTag();
    }
    if (!fs::create_directory(staging_, ec) || ec) {
        if (err) *err = "cannot create " + staging_ + ": " + (ec ? ec.message() : std::string("exists"));
        if (created_) fs::remove(target_, ec);
        return false;
    }
    open_ = true;
    if (!markPublishOwned(staging_, target_, err)) {
        abort();
        return false;
    }
    return true;
}

void Publisher::abort() {
    if (!open_) return;
    std::error_code ec;
    fs::remove_all(staging_, ec);
    if (created_) fs::remove(target_, ec);  // not recursive: only while it is empty
    open_ = false;
}

// The completion marker, written last: the key, then the size and checksum of every file.
bool Publisher::writeMarker(const std::string& key) const {
    std::ostringstream m;
    m << "key " << key << "\n";
    for (const std::string& rel : listFiles(staging_)) {
        std::error_code ec;
        const std::uintmax_t size = fs::file_size(fs::path(staging_) / rel, ec);
        char b[40];
        std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(fnv1aFile((fs::path(staging_) / rel).string())));
        m << "file " << size << " " << b << " " << rel << "\n";
    }
    std::ofstream f(fs::path(staging_) / kPublishMarker, std::ios::binary);
    f << m.str();
    f.close();
    return static_cast<bool>(f);
}

bool Publisher::commitIntoFolder(std::string* err) {
    if (!moveIntoFolder(staging_, target_, err)) {
        rollBackIntoFolder(staging_, target_);
        abort();
        return false;
    }
    open_ = false;
    std::error_code ec;
    fs::remove_all(staging_, ec);  // the replaced files; a leftover goes at the next recovery
    return true;
}

bool Publisher::swapIntoTarget(std::string* err) {
    // The target may have changed since begin(): checked again right before the swap.
    std::string why;
    if (!replaceable(target_, policy_, &why)) {
        if (err) *err = why;
        abort();
        return false;
    }
    std::error_code ec;
    fs::path old;
    if (ownType(target_) == fs::file_type::directory) {
        // The previous output carries the owner token BEFORE it is moved aside, so a crash in
        // between leaves a sibling that recovery recognises as ours.
        if (!markPublishOwned(target_, target_, err)) {
            abort();
            return false;
        }
        old = target_ + ".old-" + randomTag();
        fs::rename(target_, old, ec);
        if (ec) {
            fs::remove(fs::path(target_) / kPublishOwner, ec);
            if (err) *err = "cannot move the previous " + target_ + " aside: " + ec.message();
            abort();
            return false;
        }
    }
    fs::rename(staging_, target_, ec);
    if (ec) {
        // Put the previous output back; the staging directory is discarded.
        std::error_code ec2;
        if (!old.empty()) {
            fs::rename(old, target_, ec2);
            fs::remove(fs::path(target_) / kPublishOwner, ec2);
        }
        if (err) *err = "cannot publish " + target_ + ": " + ec.message();
        abort();
        return false;
    }
    open_ = false;
    fs::remove(fs::path(target_) / kPublishOwner, ec);
    if (!old.empty()) fs::remove_all(old, ec);
    return true;
}

bool Publisher::commit(const std::string& key, std::string* err) {
    if (!open_) {
        if (err) *err = "commit without begin";
        return false;
    }
    if (!writeMarker(key)) {
        if (err) *err = "cannot write the completion marker";
        abort();
        return false;
    }
    return policy_ == PublishPolicy::IntoFolder ? commitIntoFolder(err) : swapIntoTarget(err);
}

std::string recoverPublished(const std::string& target) {
    std::string log;
    std::error_code ec;
    const fs::path t = canonical(target);
    if (ownType(t) == fs::file_type::not_found) {
        // The newest owned backup that is a complete publication is restored.
        std::vector<fs::path> olds = ownedSiblings(t, ".old-", &log);
        std::sort(olds.begin(), olds.end(), [](const fs::path& a, const fs::path& b) {
            std::error_code e1, e2;
            return fs::last_write_time(a, e1) > fs::last_write_time(b, e2);
        });
        for (const fs::path& o : olds) {
            const MarkerCheck m = checkMarker(o);
            if (!m.verified) continue;
            fs::rename(o, t, ec);
            if (!ec) {
                fs::remove(t / kPublishOwner, ec);
                log += "restored " + o.filename().string() + "; ";
            }
            break;
        }
    }
    for (const fs::path& s : ownedSiblings(t, ".staging-", nullptr)) {
        fs::remove_all(s, ec);
        log += "removed " + s.filename().string() + "; ";
    }
    // An IntoFolder staging directory: a commit a crash stopped is undone first.
    for (const fs::path& s : ownedEntries(t, ".staging-", t, nullptr)) {
        rollBackIntoFolder(s, t);
        fs::remove_all(s, ec);
        log += "removed " + s.filename().string() + "; ";
    }
    if (ownType(t) == fs::file_type::directory)
        for (const fs::path& o : ownedSiblings(t, ".old-", nullptr)) {
            fs::remove_all(o, ec);
            log += "removed " + o.filename().string() + "; ";
        }
    return log;
}

bool isPublished(const std::string& target, const std::string& key, std::string* why) {
    const MarkerCheck m = checkMarker(fs::path(target));
    if (!m.present) {
        if (why) *why = m.why;
        return false;
    }
    if (m.key != key) {
        if (why) *why = "different key";
        return false;
    }
    if (!m.verified) {
        if (why) *why = m.why.find("unlisted") != std::string::npos ? "unlisted files" : m.why.find("changed") != std::string::npos ? m.why + " (checksum differs)" : m.why;
        return false;
    }
    return true;
}

bool markPublishOwned(const std::string& dir, const std::string& target, std::string* err) {
    std::ofstream f(fs::path(dir) / kPublishOwner, std::ios::binary);
    f << "knit-publish-owner 1\ntarget " << canonical(target).generic_string() << "\n";
    f.close();
    if (!f) {
        if (err) *err = "cannot write the owner token in " + dir;
        return false;
    }
    return true;
}

} // namespace knit
