// Residency (gpu/device.h): the priority lease of a wool GPU job over the live preview's memory.
// A FakeTenant stands for the preview; setAvailabilityForTests stands for the device's readings, so
// "the tenant gave its memory back" can be observed by the lease exactly as it would be on the card.
// The reservation cases need a device (the allocator); the lease cases do not.
#include <atomic>
#include <cstdio>
#include <chrono>
#include <memory>
#include <mutex>
#include <thread>

#include "gpu_common.h"
#include "gpu/device/runtime.h"
#include "model/cancel.h"

using namespace knit;
using namespace knit::gpu;

namespace {

// The device as the lease reads it: `available` bytes, changed by the fake tenant's release.
std::atomic<size_t> gAvailable{0};

MemoryAvailability fakeReading() {
    MemoryAvailability a;
    a.budgetKnown = true;
    a.availableBytes = gAvailable.load();
    return a;
}

constexpr size_t MiB = size_t(1) << 20;

// The preview: holds `held` bytes; release() takes `releaseMs`, then makes them available.
class FakeTenant : public ResidencyTenant {
public:
    explicit FakeTenant(size_t held, int releaseMs = 0) : held_(held), releaseMs_(releaseMs) {}
    size_t reclaimableBytes() override { return held_.load(); }
    size_t pendingBytes() override { return held_.load() > 0 ? pending.load() : 0; }
    size_t release(int gracefulMs) override {
        ++releases;
        releasing = true;
        lastGracefulMs = gracefulMs;
        std::this_thread::sleep_for(std::chrono::milliseconds(releaseMs_));
        const size_t freed = held_.exchange(0);
        gAvailable += freed;
        return freed;
    }
    void mayRebuild() override {
        ++rebuildSignals;
        rebuildSawPriority = priorityHeld();
    }
    std::atomic<int> releases{0}, rebuildSignals{0}, lastGracefulMs{0};
    std::atomic<bool> rebuildSawPriority{true};
    std::atomic<bool> releasing{false};   // a release() has started
    std::atomic<size_t> pending{0};       // a "load in flight" (pendingBytes)

private:
    std::atomic<size_t> held_;
    int releaseMs_;
};

// Installs the fake readings for one case and restores the real ones.
struct FakeDevice {
    explicit FakeDevice(size_t available) {
        gAvailable = available;
        setAvailabilityForTests(&fakeReading);
    }
    ~FakeDevice() { setAvailabilityForTests(nullptr); }
};

} // namespace

// The job does not fit while the preview holds memory: the lease releases the tenant, re-reads,
// and is granted. The graceful timeout reaches the tenant.
TEST(lease_reclaims_then_grants) {
    FakeDevice device(1000 * MiB);
    FakeTenant preview(6000 * MiB);
    TenantRegistration reg(&preview);
    {
        PriorityLease lease(4000 * MiB);
        CHECK(preview.releases == 1);
        CHECK(preview.lastGracefulMs == kGracefulReleaseMs);
        CHECK(lease.freedBytes() == 6000 * MiB);
        CHECK(priorityHeld());
    }
    CHECK(!priorityHeld());
}

// The job fits: the tenant is not released (the preview stays warm), but it still may not rebuild
// while the lease lives.
TEST(lease_fits_without_release) {
    FakeDevice device(8000 * MiB);
    FakeTenant preview(3000 * MiB);
    TenantRegistration reg(&preview);
    {
        PriorityLease lease(4000 * MiB);
        CHECK(preview.releases == 0);
        CHECK(priorityHeld());
    }
    CHECK(preview.rebuildSignals == 1);
}

// Even after the release the job does not fit: the lease reports DeviceBudgetExceeded with the
// numbers, after at most the recovery wait, and leaves no priority behind.
TEST(lease_short_reports_budget) {
    FakeDevice device(500 * MiB);
    FakeTenant preview(1000 * MiB);
    TenantRegistration reg(&preview);
    const auto t0 = std::chrono::steady_clock::now();
    bool reported = false;
    try {
        PriorityLease lease(8000 * MiB);
    } catch (const DeviceBudgetExceeded& e) {
        reported = std::string(e.what()).find("8000 MiB") != std::string::npos &&
                   std::string(e.what()).find("gave back 1000 MiB") != std::string::npos;
    }
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    CHECK(reported);
    CHECK(s < kRecoveryWaitMs / 1000.0 + 1.0);
    CHECK(!priorityHeld());
    CHECK(preview.rebuildSignals == 1);  // the failed lease ended: the preview may come back
}

// The tenant hears "may rebuild" only when the LAST lease ends, and no lease is held then.
TEST(lease_blocks_rebuild) {
    FakeDevice device(100000 * MiB);
    FakeTenant preview(0);
    TenantRegistration reg(&preview);
    std::atomic<bool> innerDone{false};
    {
        PriorityLease outer(MiB);
        std::thread other([&] {
            PriorityLease second(MiB);
            innerDone = true;
        });
        other.join();
        CHECK(innerDone);
        CHECK(preview.rebuildSignals == 0);  // the outer lease still lives
        {
            PriorityLease nested(MiB);  // same thread: a no-op
            CHECK(nested.nested());
        }
        CHECK(preview.rebuildSignals == 0);
    }
    CHECK(preview.rebuildSignals == 1);
    CHECK(!preview.rebuildSawPriority);
}

// Two jobs need memory at once: the tenant is released once; the second lease waits for that
// reclaim and is granted from it.
TEST(lease_concurrent_coalesce) {
    FakeDevice device(1000 * MiB);
    FakeTenant preview(6000 * MiB, 300);
    TenantRegistration reg(&preview);
    std::atomic<int> granted{0};
    std::thread a([&] {
        PriorityLease l(3000 * MiB);
        ++granted;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    });
    while (!preview.releasing) std::this_thread::sleep_for(std::chrono::milliseconds(5));
    std::thread b([&] {
        PriorityLease l(3000 * MiB);
        ++granted;
    });
    a.join();
    b.join();
    CHECK(granted == 2);
    CHECK(preview.releases == 1);
}

// A job cancelled while its lease waits for another lease's reclaim leaves with Cancelled and no
// priority; the reclaim itself finishes.
TEST(lease_cancel_during_reclaim) {
    FakeDevice device(0);
    FakeTenant preview(6000 * MiB, 800);
    TenantRegistration reg(&preview);
    std::thread first([&] { PriorityLease l(1000 * MiB); });
    while (!preview.releasing) std::this_thread::sleep_for(std::chrono::milliseconds(5));
    CancelToken token;
    bool cancelled = false;
    std::thread second([&] {
        ThreadCancelScope scope(&token);
        try {
            PriorityLease l(1000 * MiB);
        } catch (const Cancelled&) {
            cancelled = true;
        }
    });
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    token.cancel();
    second.join();
    first.join();
    CHECK(cancelled);
    CHECK(preview.releases == 1);
    CHECK(!priorityHeld());
}

// Unregistering waits for a release in flight: the tenant never disappears under a call.
TEST(tenant_unregister_waits_inflight) {
    FakeDevice device(0);
    FakeTenant preview(2000 * MiB, 400);
    auto reg = std::make_unique<TenantRegistration>(&preview);
    std::thread job([&] { PriorityLease l(1000 * MiB); });
    while (!preview.releasing) std::this_thread::sleep_for(std::chrono::milliseconds(5));
    const auto t0 = std::chrono::steady_clock::now();
    reg.reset();
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    job.join();
    std::printf("unregister waited %.0f ms\n", ms);
    CHECK(ms > 300);  // it waited for (most of) the 400 ms release that was running
    CHECK(preview.releases == 1);
}

// A job without an estimate (kJobBudget) releases the tenants when the budget does not fit, and runs with
// what there is afterwards: it never reports, even when the card stays short (other processes).
TEST(lease_job_budget_never_reports) {
    REQUIRE_GPU();  // jobBytes() is the allocator's budget
    FakeDevice device(10 * MiB);
    FakeTenant preview(300 * MiB);
    TenantRegistration reg(&preview);
    const auto t0 = std::chrono::steady_clock::now();
    {
        PriorityLease lease(kJobBudget);  // no throw
        CHECK(preview.releases == 1);
    }
    FakeTenant idle(0);
    TenantRegistration reg2(&idle);
    {
        PriorityLease lease(kJobBudget);  // nothing to release: granted at once
        CHECK(idle.releases == 0);
    }
    const double s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    CHECK(s < kRecoveryWaitMs / 1000.0 + 1.0);
}

// Without a tenant a lease reads nothing and does nothing but hold priority (today's behaviour).
TEST(lease_without_tenant_reads_nothing) {
    FakeDevice device(0);  // would refuse any job if it were read
    PriorityLease lease(kJobBudget);
    CHECK(priorityHeld());
}

// An allocator reservation in a thread without a lease takes one while a tenant exists: an
// oversized request is admitted only after the tenant gave its memory back.
TEST(reservation_oversize_asks_tenant) {
    REQUIRE_GPU();
    FakeDevice device(100 * MiB);
    FakeTenant preview(Reservation::exclusiveBytes() + (size_t(1) << 30));
    TenantRegistration reg(&preview);
    {
        Reservation r(Reservation::exclusiveBytes());
        CHECK(preview.releases == 1);
        CHECK(priorityHeld());
    }
    CHECK(!priorityHeld());
    CHECK(memStats().reserved == 0);
}

// Without a tenant the reservation path is unchanged: no lease, no reading.
TEST(reservation_unchanged_without_tenant) {
    REQUIRE_GPU();
    FakeDevice device(0);
    {
        Reservation r(Reservation::exclusiveBytes());
        CHECK(!priorityHeld());
    }
    CHECK(memStats().reserved == 0);
}

// A load in flight keeps allocating after the lease is granted: the lease counts the tenant's pending
// growth, and releases it although the job alone would fit now.
TEST(lease_counts_pending_growth) {
    FakeDevice device(5000 * MiB);
    FakeTenant preview(300 * MiB);
    preview.pending = 4000 * MiB;
    TenantRegistration reg(&preview);
    {
        PriorityLease lease(2000 * MiB);
        CHECK(preview.releases == 1);
    }
    FakeTenant quiet(300 * MiB);
    TenantRegistration reg2(&quiet);
    preview.pending = 0;
    {
        PriorityLease lease(2000 * MiB);
        CHECK(quiet.releases == 0);
    }
}
