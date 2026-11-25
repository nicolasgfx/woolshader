#pragma once
// Test registry and checks of knit_tests_s6 (same conventions as tests.cpp, which keeps
// its own copy: S6 tests live in a separate executable so the S5 / fibre-perf branches,
// which edit tests.cpp, merge without conflicts). Run one case: knit_tests_s6 <name>.
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "model/env.h"

namespace knit_test {

struct Case { const char* name; void (*fn)(); };
inline std::vector<Case>& registry() { static std::vector<Case> r; return r; }
struct Reg { Reg(const char* n, void (*f)()) { registry().push_back({n, f}); } };
// Thrown by SKIP(): reported as skipped, never as passed.
struct Skipped { const char* why; };

// The folder of teddy knit caches (*.knitcache) that some GPU cases load: the environment variable
// KNIT_TEST_KNITCACHE; empty when it is unset, and those cases SKIP with that reason.
inline std::string knitCacheFixtureDir() {
    const auto dir = knit::envVar("KNIT_TEST_KNITCACHE");
    return dir ? *dir : std::string();
}

} // namespace knit_test

#define TEST(name) static void name(); static knit_test::Reg reg_##name(#name, name); static void name()
#define CHECK(cond) do { if (!(cond)) { char b_[512]; std::snprintf(b_, sizeof b_, "%s:%d: CHECK(%s)", __FILE__, __LINE__, #cond); throw std::runtime_error(b_); } } while (0)
#define CHECK_NEAR(a, b, tol) do { double a_ = (a), b_ = (b); if (!(std::fabs(a_ - b_) <= (tol))) { char m_[512]; std::snprintf(m_, sizeof m_, "%s:%d: |%s - %s| = |%.17g - %.17g| > %g", __FILE__, __LINE__, #a, #b, a_, b_, double(tol)); throw std::runtime_error(m_); } } while (0)
#define SKIP(why) throw knit_test::Skipped{why}

namespace knit_test {
// A case that needs a model from scenes/ (only the wool swatch's model ships with the source) passes the path that
// findProjectFile found; an empty path means this checkout has no such model, so the case skips.
inline void requireScene(const std::string& resolvedPath) {
    if (resolvedPath.empty()) SKIP("the model is not in this checkout (only the wool swatch's model ships with the source)");
}
} // namespace knit_test
