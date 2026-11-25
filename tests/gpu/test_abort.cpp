// GPU port tests: the cooperative abort of long kernels (runtime.h abortFlagDevice, waitForKernels),
// on the curvature kernel. A curve with coincident controls has
// zero-speed spans that never become tight, so their subdivision runs to the full depth: one
// launch takes long enough to cancel inside it. Checked: the cancel stops the launch promptly and
// throws (nothing of it is returned), the calling thread's flag is clear afterwards, and the same
// thread's next launch completes with the result of a run that was never cancelled (a flag left
// raised would end the subdivision early and change the bounds).
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <functional>
#include <memory>
#include <thread>
#include <vector>

#include "gpu_common.h"
#include "model/cancel.h"
#include "knitting/checks/checks.h"
#include "model/env.h"
#include "gpu/validation/checks.h"
#include "gpu/device/runtime.h"

using namespace knit;

namespace {

CurveSet degenerateLine() {
    CurveSet cs;
    cs.kind = CurveKind::Surface;
    YarnPath y;
    y.id = 0;
    for (int i = 0; i <= 6; ++i) y.curve.ctrl.push_back({double(i), 0.1 * (i % 2), 0});
    for (int k = 0; k < 4; ++k) y.curve.ctrl.push_back({7, 0, 0});  // a span collapsed to a point
    for (int i = 8; i <= 12; ++i) y.curve.ctrl.push_back({double(i), 0.1 * (i % 2), 0});
    y.finalize();
    cs.yarns.push_back(std::move(y));
    return cs;
}

bool sameResults(const std::vector<CurvatureResult>& a, const std::vector<CurvatureResult>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i].kappaHi != b[i].kappaHi || a[i].speedLo != b[i].speedLo || a[i].leaves != b[i].leaves ||
            a[i].spanKappaHi != b[i].spanKappaHi || a[i].regular != b[i].regular || a[i].noFold != b[i].noFold)
            return false;
    return true;
}

double msSince(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
}

} // namespace

TEST(gpu_abort_cancel_then_rerun) {
    REQUIRE_GPU();
    const CurveSet cs = degenerateLine();
    const Allowances al = Allowances::forExtent(cs.extent());
    const double radius = 0.5;
    const gpu::UploadedCurves d = gpu::uploadCurves(cs);
    // (1) The reference: a fresh thread, whose flag is never raised.
    std::vector<CurvatureResult> ref;
    double msRef = 0;
    std::thread([&] {
        const auto t0 = std::chrono::steady_clock::now();
        ref = gpu::checkCurvatureOnDevice(d, radius, al);
        msRef = msSince(t0);
    }).join();
    CHECK(!ref.empty());
    std::printf("    degenerate span: one curvature launch %.0f ms, %d leaves\n", msRef, ref[0].leaves);
    CHECK(msRef > 400);  // long enough to cancel 100 ms in
    // (2) Cancelled 100 ms into the launch on this thread: it throws, promptly.
    CancelToken token;
    bool threw = false;
    double msStop = 0;
    {
        ThreadCancelScope scope(&token);
        std::thread canceller([&] {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            token.cancel();
        });
        const auto t0 = std::chrono::steady_clock::now();
        try {
            gpu::checkCurvatureOnDevice(d, radius, al);
        } catch (const Cancelled&) {
            threw = true;
        }
        msStop = msSince(t0) - 100;
        canceller.join();
    }
    std::printf("    cancelled 100 ms in: stopped %.0f ms after the request\n", msStop);
    CHECK(threw);
    CHECK(msStop < 1000);
    CHECK(gpu::abortFlagValue() == 0);  // the flag was cleared after the abort
    // (3) The same thread again, no cancel: completes, identical to the reference.
    const auto t1 = std::chrono::steady_clock::now();
    const std::vector<CurvatureResult> again = gpu::checkCurvatureOnDevice(d, radius, al);
    std::printf("    rerun on the same thread: %.0f ms, identical to the reference: %s\n", msSince(t1),
                sameResults(ref, again) ? "yes" : "NO");
    CHECK(sameResults(ref, again));
    CHECK(gpu::abortFlagValue() == 0);
}

// A cancelled job leaves the device's FIFO admission queue (runtime.h Reservation): behind
// another waiter its ticket is skipped when the queue reaches it, at the front at once; it holds
// nothing, its thread's nesting depth is restored (the thread's next request is a real one),
// the requests behind it are served, and nothing stays reserved. Synchronized on the queue's
// state (queuedRequests), not on sleeps.
TEST(gpu_reservation_cancel_leaves_queue) {
    REQUIRE_GPU();
    const size_t all = gpu::Reservation::exclusiveBytes();
    CHECK(gpu::memStats().reserved == 0 && gpu::queuedRequests() == 0);
    auto waitFor = [](const char* what, const std::function<bool()>& cond) {
        const auto t0 = std::chrono::steady_clock::now();
        while (!cond()) {
            if (msSince(t0) > 10000) {
                std::printf("    timeout waiting for %s\n", what);
                return false;
            }
            std::this_thread::yield();
        }
        return true;
    };
    // (1) Behind another waiter: A holds the device, C waits, B waits behind C, D behind B.
    {
        auto a = std::make_unique<gpu::Reservation>(all);
        std::atomic<bool> cDone{false}, dDone{false}, bCancelled{false}, bNextReal{false};
        std::thread c([&] {
            gpu::Reservation r(all);
            cDone = true;
        });
        CHECK(waitFor("C queued", [] { return gpu::queuedRequests() == 1; }));
        CancelToken tb;
        std::thread b([&] {
            ThreadCancelScope scope(&tb);
            try {
                gpu::Reservation r(all);
            } catch (const Cancelled&) {
                bCancelled = true;
            }
            if (gpu::leaseDepth() == 0) {
                ThreadCancelScope none(nullptr);
                gpu::Reservation r2(1);  // waits its turn like any request
                bNextReal = !r2.nested();
            }
        });
        CHECK(waitFor("B queued", [] { return gpu::queuedRequests() == 2; }));
        std::thread d([&] {
            gpu::Reservation r(1);
            dDone = true;
        });
        CHECK(waitFor("D queued", [] { return gpu::queuedRequests() == 3; }));
        tb.cancel();
        CHECK(waitFor("B cancelled", [&] { return bCancelled.load(); }));
        // B's ticket waits behind C to be skipped; B's next request queues behind D (C, D, B's next).
        CHECK(waitFor("B's next request queued", [] { return gpu::queuedRequests() == 3; }));
        CHECK(gpu::abandonedTickets() == 1);
        CHECK(!cDone && !dDone);  // A still holds the device
        a.reset();
        c.join();
        d.join();
        b.join();
        CHECK(cDone && dDone && bCancelled && bNextReal);
        CHECK(gpu::abandonedTickets() == 0 && gpu::queuedRequests() == 0 && gpu::memStats().reserved == 0);
    }
    // (2) At the front: B is next but there is no room; cancelled, the queue moves on at once.
    {
        auto a = std::make_unique<gpu::Reservation>(all);
        CancelToken tb;
        std::atomic<bool> bCancelled{false};
        std::thread b([&] {
            ThreadCancelScope scope(&tb);
            try {
                gpu::Reservation r(all);
            } catch (const Cancelled&) {
                bCancelled = true;
            }
            CHECK(gpu::leaseDepth() == 0);
        });
        CHECK(waitFor("B queued", [] { return gpu::queuedRequests() == 1; }));
        tb.cancel();
        b.join();
        CHECK(bCancelled && gpu::abandonedTickets() == 0 && gpu::queuedRequests() == 0);
        a.reset();
    }
    // (3) A fresh request afterwards is granted at once, and nothing stays reserved.
    {
        gpu::Reservation r(all);
        CHECK(!r.nested() && gpu::memStats().reserved == all);
    }
    CHECK(gpu::memStats().reserved == 0 && gpu::leaseDepth() == 0);
}

// gpu::AllocationCap: a capped allocation that would raise the tracked total
// above the cap fails (capHit, nothing counted); another thread's allocation counts towards the
// total; a failed cudaMalloc gives its reserved bytes back; caps nest and end with their scope.
TEST(gpu_allocation_cap_counts_all_threads) {
    REQUIRE_GPU();
    const size_t mb = size_t(1) << 20, base = gpu::memStats().current;
    std::string err;
    {
        gpu::AllocationCap cap(base + 8 * mb);
        CHECK(!gpu::AllocationCap::capHit());
        void* a = gpu::devAlloc(4 * mb, &err);
        CHECK(a != nullptr);
        CHECK(gpu::devAlloc(6 * mb, &err) == nullptr);  // 4 + 6 > 8
        CHECK(gpu::AllocationCap::capHit() && err.find("device memory cap") != std::string::npos);
        CHECK(gpu::memStats().current == base + 4 * mb);
        void* c = nullptr;
        std::thread([&] {  // no cap on this thread: allocated, and counted
            std::string e;
            c = gpu::devAlloc(3 * mb, &e);
        }).join();
        CHECK(c != nullptr);
        CHECK(gpu::devAlloc(2 * mb, &err) == nullptr);  // 4 + 3 + 2 > 8: the other thread's bytes count
        gpu::devFree(c);
        void* d = gpu::devAlloc(2 * mb, &err);  // 4 + 2 <= 8
        CHECK(d != nullptr);
        {
            gpu::AllocationCap inner(base + 7 * mb);  // the innermost cap applies
            CHECK(gpu::devAlloc(2 * mb, &err) == nullptr);
        }
        void* f = gpu::devAlloc(2 * mb, &err);  // the outer cap again: 4 + 2 + 2 <= 8
        CHECK(f != nullptr);
        gpu::devFree(a);
        gpu::devFree(d);
        gpu::devFree(f);
    }
    CHECK(gpu::devAlloc(size_t(1) << 50, &err) == nullptr);  // cudaMalloc fails: its bytes are given back
    CHECK(gpu::memStats().current == base);
    void* g = gpu::devAlloc(16 * mb, &err);  // no cap after the scope
    CHECK(g != nullptr);
    gpu::devFree(g);
    CHECK(gpu::memStats().current == base);
}

// gpu::solveCapBytes: without KNIT_GPU_BUDGET_MB, the tracked memory plus the
// device's free memory now, minus a reserve: never below what is tracked, never above tracked + total.
TEST(gpu_solve_cap_follows_free_memory) {
    REQUIRE_GPU();
    if (envVar("KNIT_GPU_BUDGET_MB").has_value()) {
        std::printf("    KNIT_GPU_BUDGET_MB is set: the cap is the explicit budget (%zu MB)\n", gpu::solveCapBytes() >> 20);
        CHECK(gpu::solveCapBytes() == gpu::deviceBudget());
        return;
    }
    const size_t cur = gpu::memStats().current, cap = gpu::solveCapBytes(), total = gpu::deviceProperties().totalBytes;
    std::printf("    tracked %.0f MB, cap %.0f MB, device %.0f MB\n", cur / 1048576.0, cap / 1048576.0, total / 1048576.0);
    CHECK(cap >= cur);
    CHECK(cap <= cur + total);
    void* p = gpu::devAlloc(size_t(64) << 20, nullptr);  // tracked memory moves with the cap
    CHECK(p != nullptr);
    const size_t cap2 = gpu::solveCapBytes();
    gpu::devFree(p);
    std::printf("    with 64 MB more tracked: cap %.0f MB\n", cap2 / 1048576.0);
    CHECK(cap2 + (size_t(16) << 20) >= cap && cap2 <= cap + (size_t(16) << 20));  // free fell by ~64 MB, tracked rose by 64 MB
}
