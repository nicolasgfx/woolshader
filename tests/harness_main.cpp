// The main of the harness-based test programs (knit_tests_gpu, knit_tests_export, knit_tests_import,
// knit_tests_architecture): every source registers its cases with harness.h.
// Run one case: <program> <name>; all: <program>; list: <program> --list. A case
// that needs a CUDA device SKIPs (never passes) without one.
#include <cstdio>
#include <cstring>
#include <exception>

#include "harness.h"

int main(int argc, char** argv) {
    using knit_test::registry;
    if (argc > 1 && std::strcmp(argv[1], "--list") == 0) {
        for (const knit_test::Case& c : registry()) std::printf("%s\n", c.name);
        return 0;
    }
    int failed = 0, ran = 0, skipped = 0;
    for (const knit_test::Case& c : registry()) {
        if (argc > 1 && std::strcmp(argv[1], c.name) != 0) continue;
        ++ran;
        try {
            c.fn();
            std::printf("[ OK ] %s\n", c.name);
        } catch (const knit_test::Skipped& s) {
            ++skipped;
            std::printf("[SKIP] %s: %s\n", c.name, s.why);
        } catch (const std::exception& e) {
            ++failed;
            std::printf("[FAIL] %s: %s\n", c.name, e.what());
        }
        std::fflush(stdout);
    }
    if (ran == 0) { std::printf("no test named %s\n", argc > 1 ? argv[1] : "?"); return 2; }
    if (skipped) std::printf("%d skipped\n", skipped);
    if (failed) return 1;
    return skipped == ran ? 77 : 0;  // 77: every case run was skipped (CTest SKIP_RETURN_CODE)
}
