// Publication safety (knit_tests_publish): the export target is a folder inside a USER folder
// (the imported file's parent), so publication must never delete or replace anything it did
// not create.
//   * recovery touches only siblings named <target>.staging-<16 hex> / <target>.old-<16 hex>
//     that carry this target's owner token (.knit-publish-owner);
//   * an existing target is replaced only when it is an empty directory or a complete
//     publication without untracked files;
//   * the exported files' name is one plain path component.
// An export into a folder the user chose (IntoFolder) replaces only the files it writes:
//   * every other file stays, and nothing appears next to the folder;
//   * a failed commit, or recovery after a crash during one, restores what it replaced.
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

#include "harness.h"
#include "export/export.h"
#include "io/publish.h"

using namespace knit;
namespace fs = std::filesystem;

namespace {

std::string freshDir(const std::string& name) {
    const fs::path p = fs::temp_directory_path() / "knit_publish_tests" / name;
    std::error_code ec;
    fs::remove_all(p, ec);
    fs::create_directories(p, ec);
    return (p / "export").string();  // the target inside a "user" parent folder
}

void writeText(const std::string& path, const std::string& text) {
    fs::create_directories(fs::path(path).parent_path());
    std::ofstream f(path, std::ios::binary);
    f << text;
}

std::string readText(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

bool publishText(const std::string& t, const std::string& key, const std::string& text, std::string* err) {
    Publisher p;
    if (!p.begin(t, err)) return false;
    writeText(p.staging() + "/scene.json", text);
    return p.commit(key, err);
}

std::string parentOf(const std::string& t) { return fs::path(t).parent_path().string(); }

// An export of two files into the folder `t`.
bool exportInto(const std::string& t, const std::string& text, std::string* err) {
    Publisher p;
    if (!p.begin(t, err, PublishPolicy::IntoFolder)) return false;
    writeText(p.staging() + "/scene.json", text);
    writeText(p.staging() + "/meshes/part.ply", text);
    return p.commit("k", err);
}

// Nothing of the publication left: no staging directory in the folder or beside it.
void checkNoLeftovers(const std::string& t) {
    for (const std::string& dir : {t, parentOf(t)})
        for (const auto& e : fs::directory_iterator(dir)) CHECK(e.path().filename().string().find(".staging-") == std::string::npos);
}

} // namespace

// Directories next to the target that merely LOOK like publication leftovers belong to the
// user: recovery neither restores nor deletes them.
TEST(publish_recovery_keeps_foreign_siblings) {
    const std::string t = freshDir("foreign_siblings");
    std::string err;
    CHECK(publishText(t, "k1", "one", &err));
    const std::string names[] = {t + ".old-x", t + ".staging-abc", t + ".old-0123456789abcdef", t + ".staging-fedcba9876543210",
                                 t + ".staging-user notes"};
    for (const std::string& n : names) writeText(n + "/keep.txt", "user data");
    const std::string log = recoverPublished(t);
    for (const std::string& n : names) CHECK(readText(n + "/keep.txt") == "user data");
    CHECK(readText(t + "/scene.json") == "one" && isPublished(t, "k1"));
    // The same with the target missing: an unowned ".old-" directory is not "restored".
    fs::remove_all(t);
    recoverPublished(t);
    CHECK(!fs::exists(t));
    for (const std::string& n : names) CHECK(readText(n + "/keep.txt") == "user data");
}

// A crash between moving the previous output aside and renaming the staging directory: the
// OWNED leftovers are recovered (the previous output restored, the staging removed).
TEST(publish_recovery_owned_leftovers) {
    const std::string t = freshDir("owned_leftovers");
    std::string err;
    CHECK(publishText(t, "k1", "one", &err));
    const std::string old = t + ".old-0123456789abcdef", staging = t + ".staging-fedcba9876543210";
    CHECK(markPublishOwned(t, t, &err));
    fs::rename(t, old);
    fs::create_directories(staging);
    writeText(staging + "/scene.json", "half");
    CHECK(markPublishOwned(staging, t, &err));
    const std::string log = recoverPublished(t);
    CHECK(!log.empty());
    CHECK(readText(t + "/scene.json") == "one" && isPublished(t, "k1"));
    CHECK(!fs::exists(old) && !fs::exists(staging));
}

// An owner token naming ANOTHER target does not make a sibling ours.
TEST(publish_recovery_token_of_other_target) {
    const std::string t = freshDir("other_token");
    std::string err;
    CHECK(publishText(t, "k1", "one", &err));
    const std::string staging = t + ".staging-fedcba9876543210";
    fs::create_directories(staging);
    writeText(staging + "/keep.txt", "x");
    CHECK(markPublishOwned(staging, parentOf(t) + "/other", &err));
    recoverPublished(t);
    CHECK(readText(staging + "/keep.txt") == "x");
}

// A publication into which the user put files is not replaced (the files would be lost).
TEST(publish_refuses_untracked_files) {
    const std::string t = freshDir("untracked");
    std::string err;
    CHECK(publishText(t, "k1", "one", &err));
    writeText(t + "/my notes.txt", "mine");
    err.clear();
    CHECK(!publishText(t, "k2", "two", &err));
    CHECK(err.find("my notes.txt") != std::string::npos);
    CHECK(readText(t + "/my notes.txt") == "mine" && readText(t + "/scene.json") == "one");
    for (const auto& e : fs::directory_iterator(parentOf(t))) CHECK(e.path().filename().string().find(".staging-") == std::string::npos);
}

// A folder that is not a publication (no marker) is never replaced; an empty one is.
TEST(publish_refuses_foreign_target) {
    const std::string t = freshDir("foreign_target");
    writeText(t + "/thesis.docx", "important");
    std::string err;
    CHECK(!publishText(t, "k1", "one", &err));
    CHECK(!err.empty() && readText(t + "/thesis.docx") == "important" && !fs::exists(t + "/scene.json"));
    const std::string e = freshDir("empty_target");
    fs::create_directories(e);
    CHECK(publishText(e, "k1", "one", &err));
    CHECK(readText(e + "/scene.json") == "one");
}

// Files next to the target in the parent folder are never touched by a publication.
TEST(publish_keeps_parent_files) {
    const std::string t = freshDir("parent_files");
    writeText(parentOf(t) + "/model.obj", "v 0 0 0");
    writeText(parentOf(t) + "/export_old/keep.txt", "k");
    std::string err;
    CHECK(publishText(t, "k1", "one", &err));
    CHECK(publishText(t, "k2", "two", &err));
    CHECK(readText(parentOf(t) + "/model.obj") == "v 0 0 0" && readText(parentOf(t) + "/export_old/keep.txt") == "k");
    CHECK(readText(t + "/scene.json") == "two");
}

TEST(publish_export_name_rules) {
    std::string why;
    for (const char* ok : {"bunny_knit", "Teddy knit 2", "a.b"}) CHECK(validateExportName(ok, &why));
    for (const char* bad : {"", ".", "..", "a/b", "a\\b", "c:", "CON", "nul.txt", "LPT1", "name.", "name ", "a*b", "a?b", "a<b", "a|b", "a\"b"})
        CHECK(!validateExportName(bad, &why) && !why.empty());
    CHECK(!validateExportName(std::string(121, 'a'), &why));
}

// The chosen folder holds the user's files: the export lands directly in it, and they stay,
// nested ones too. Files of the same name are replaced; an earlier export's others stay as well.
TEST(publish_into_folder_keeps_other_files) {
    const std::string t = freshDir("into_folder");
    writeText(t + "/thesis.docx", "important");
    writeText(t + "/meshes/my mesh.ply", "mine");
    writeText(t + "/scene.json", "theirs");
    std::string err;
    CHECK(exportInto(t, "one", &err));
    CHECK(readText(t + "/thesis.docx") == "important" && readText(t + "/meshes/my mesh.ply") == "mine");
    CHECK(readText(t + "/scene.json") == "one" && readText(t + "/meshes/part.ply") == "one");
    CHECK(fs::exists(t + "/" + kPublishMarker) && !fs::exists(t + "/" + kPublishOwner));
    CHECK(exportInto(t, "two", &err));
    CHECK(readText(t + "/thesis.docx") == "important" && readText(t + "/scene.json") == "two");
    checkNoLeftovers(t);
}

// A folder that does not exist is created; a failed export takes it away again, and leaves an
// existing one with what it held.
TEST(publish_into_folder_creates_and_aborts) {
    const std::string t = freshDir("into_new");
    std::string err;
    {
        Publisher p;
        CHECK(p.begin(t, &err, PublishPolicy::IntoFolder) && fs::is_directory(t));
        writeText(p.staging() + "/scene.json", "half");
    }
    CHECK(!fs::exists(t));
    CHECK(exportInto(t, "one", &err));
    {
        Publisher p;
        CHECK(p.begin(t, &err, PublishPolicy::IntoFolder));
        writeText(p.staging() + "/scene.json", "half");
    }
    CHECK(readText(t + "/scene.json") == "one");
    checkNoLeftovers(t);
}

// A commit that fails half way puts everything back. "report.json" sorts after "meshes/...",
// which has moved by then: a folder of that name, the user's, stops the commit and stays.
TEST(publish_into_folder_failure_restores) {
    const std::string t = freshDir("into_failure");
    writeText(t + "/keep.txt", "k");
    writeText(t + "/meshes/old.ply", "old");
    writeText(t + "/report.json/inside.txt", "a folder");
    std::string err;
    Publisher p;
    CHECK(p.begin(t, &err, PublishPolicy::IntoFolder));
    writeText(p.staging() + "/meshes/old.ply", "new");
    writeText(p.staging() + "/meshes/added.ply", "new");
    writeText(p.staging() + "/nested/added.ply", "new");
    writeText(p.staging() + "/report.json", "new");
    CHECK(!p.commit("k", &err));
    std::printf("    %s\n", err.c_str());
    CHECK(err.find("report.json") != std::string::npos);
    CHECK(readText(t + "/keep.txt") == "k" && readText(t + "/meshes/old.ply") == "old" && readText(t + "/report.json/inside.txt") == "a folder");
    CHECK(!fs::exists(t + "/meshes/added.ply") && !fs::exists(t + "/nested") && !fs::exists(t + "/" + kPublishMarker));
    checkNoLeftovers(t);
}

// A file another program holds open cannot be replaced: the export fails and says which one, and
// the earlier files stay.
TEST(publish_into_folder_open_file) {
#ifdef _WIN32
    const std::string t = freshDir("into_open");
    std::string err;
    CHECK(exportInto(t, "one", &err));
    {
        std::ifstream held(t + "/scene.json", std::ios::binary);  // no delete sharing
        CHECK(!exportInto(t, "two", &err));
    }
    std::printf("    %s\n", err.c_str());
    CHECK(err.find("scene.json") != std::string::npos);
    CHECK(readText(t + "/scene.json") == "one" && readText(t + "/meshes/part.ply") == "one");
    checkNoLeftovers(t);
#else
    SKIP("Windows file sharing");
#endif
}

// A crash during the moves: the staged marker says the commit did not finish, so recovery puts
// the replaced files back and removes what had arrived, and leaves the user's alone.
TEST(publish_into_folder_recovery) {
    const std::string t = freshDir("into_crash");
    std::string err;
    CHECK(exportInto(t, "one", &err));
    writeText(t + "/keep.txt", "k");
    // The state after the first move: meshes/part.ply replaced, the rest still staged.
    const std::string s = t + "/.staging-fedcba9876543210";
    writeText(s + "/scene.json", "two");
    CHECK(markPublishOwned(s, t, &err));
    writeText(s + "/" + kPublishMarker, "key k\nfile 3 0000000000000000 meshes/part.ply\nfile 3 0000000000000000 scene.json\n");
    writeText(s + "/" + kPublishReplaced + "/meshes/part.ply", readText(t + "/meshes/part.ply"));
    writeText(t + "/meshes/part.ply", "two");
    CHECK(!recoverPublished(t).empty());
    CHECK(readText(t + "/meshes/part.ply") == "one" && readText(t + "/scene.json") == "one" && readText(t + "/keep.txt") == "k");
    CHECK(!fs::exists(s));
    // The same name without this target's token is the user's folder.
    writeText(s + "/keep.txt", "k");
    recoverPublished(t);
    CHECK(readText(s + "/keep.txt") == "k");
}

// A file, or a link, is no folder to write into.
TEST(publish_into_folder_refuses_a_file) {
    const std::string t = freshDir("into_file");
    writeText(t, "a file");
    std::string err;
    CHECK(!exportInto(t, "one", &err));
    CHECK(!err.empty() && readText(t) == "a file");
}

// A junction (or symlink) as the target is refused: replacing it would act on its destination.
TEST(publish_refuses_junction_target) {
#ifdef _WIN32
    const std::string t = freshDir("junction");
    const std::string real = parentOf(t) + "/real";
    writeText(real + "/keep.txt", "k");
    const std::string cmd = "cmd /c mklink /J \"" + fs::path(t).make_preferred().string() + "\" \"" + fs::path(real).make_preferred().string() + "\" >nul 2>&1";
    if (std::system(cmd.c_str()) != 0 || !fs::exists(t)) SKIP("cannot create a junction here");
    std::string err;
    CHECK(!publishText(t, "k1", "one", &err));
    CHECK(readText(real + "/keep.txt") == "k" && !fs::exists(real + "/scene.json"));
    fs::remove(t);  // the junction only
#else
    SKIP("Windows junctions");
#endif
}

int main(int argc, char** argv) {
    int failed = 0, run = 0, skipped = 0;
    const std::string only = argc > 1 ? argv[1] : "";
    if (only == "--list") {
        for (const auto& c : knit_test::registry()) std::printf("%s\n", c.name);
        return 0;
    }
    for (const auto& c : knit_test::registry()) {
        if (!only.empty() && only != c.name) continue;
        ++run;
        try {
            c.fn();
            std::printf("[PASS] %s\n", c.name);
        } catch (const knit_test::Skipped& s) {
            ++skipped;
            std::printf("[SKIP] %s: %s\n", c.name, s.why);
        } catch (const std::exception& e) {
            ++failed;
            std::printf("[FAIL] %s: %s\n", c.name, e.what());
        }
    }
    if (run == 0) {
        std::printf("no test named %s\n", only.c_str());
        return 1;
    }
    std::printf("%d run, %d failed, %d skipped\n", run, failed, skipped);
    if (failed) return 1;
    return skipped == run ? 77 : 0;  // 77: every case run was skipped (CTest SKIP_RETURN_CODE)
}
