#pragma once
// S6: filesystem publication of generated directories (exports, texture map
// sets, caches). Writers never open their final paths:
//   1. begin(): checks that the target may be replaced (PublishPolicy) and creates a unique
//      staging directory next to it, <target>.staging-<16 hex>, holding an owner token
//      (.knit-publish-owner) that names the target;
//   2. the writer fills it (and validates what it wrote);
//   3. commit(key): the completion marker (key + size and FNV-1a 64 of every file) is
//      written LAST; the replacement rule is checked again; an existing target gets the
//      owner token and is renamed to <target>.old-<16 hex>, the staging directory is renamed
//      to the target, the old one deleted. abort() (cancel, failure) deletes the staging
//      directory; the previous published output stays byte-identical.
// The export target lives in a USER folder (an imported file's parent), so nothing else is
// ever touched: recovery acts only on siblings with the exact name pattern AND this target's
// owner token (a user's "export.old-x" folder is left alone); a target with files that are
// not part of the previous publication, or a folder that is no publication, is never
// replaced; a link or junction as the target is refused (never followed).
// recover() at startup: a missing target with an owned, complete .old-* sibling is
// restored; owned leftover .staging-* siblings are deleted. A CACHE HIT is a target whose
// marker exists, carries the requested key and whose checksums verify; staging and
// unmarked directories never are.
// A folder the user chose as the export's destination (PublishPolicy::IntoFolder) is never
// replaced or moved: the staging directory is made INSIDE it, .staging-<16 hex>, and commit
// moves each staged file to its place, the marker last. A file of the same name goes to
// <staging>/.replaced first, so a failed move puts everything back, and recover() does the
// same after a crash between the moves. Any other file in the folder stays as it is.
#include <cstdint>
#include <string>

namespace knit {

// What an existing target may be replaced by a publication:
//   Strict      (the default): a missing target, an empty folder, or a complete earlier
//               publication without untracked or changed files; anything else is refused;
//   CacheEntry  (the tool's own caches): also a damaged earlier publication (marker present);
//   IntoFolder  (the export job): any plain folder, or a missing one (created), written into.
enum class PublishPolicy { Strict, CacheEntry, IntoFolder };

class Publisher {
public:
    ~Publisher();  // aborts an uncommitted staging directory
    bool begin(const std::string& target, std::string* err, PublishPolicy policy = PublishPolicy::Strict);
    const std::string& staging() const { return staging_; }
    bool commit(const std::string& key, std::string* err);
    void abort();

private:
    bool writeMarker(const std::string& key) const;
    bool commitIntoFolder(std::string* err);
    bool swapIntoTarget(std::string* err);

    std::string target_, staging_;
    PublishPolicy policy_ = PublishPolicy::Strict;
    bool open_ = false;
    bool created_ = false;  // IntoFolder: begin() made the folder, abort() removes it again
};

inline constexpr const char* kPublishMarker = ".publish-complete";

// Holds the files an IntoFolder commit replaced, inside its staging directory, until it is done.
inline constexpr const char* kPublishReplaced = ".replaced";

// Startup recovery of one target (see above). Returns a short log line ("" if nothing done).
std::string recoverPublished(const std::string& target);
// True iff target is a complete publication of `key` (marker present, key equal, every
// listed file present with its size and checksum, no unlisted file).
bool isPublished(const std::string& target, const std::string& key, std::string* why = nullptr);

// The owner token of a publication directory (`dir` is the target itself, or one of its
// staging / backup siblings): names the target it belongs to. Publisher writes it; exposed
// for the crash-recovery tests.
inline constexpr const char* kPublishOwner = ".knit-publish-owner";
bool markPublishOwned(const std::string& dir, const std::string& target, std::string* err);

} // namespace knit
