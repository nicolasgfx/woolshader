// Residency (gpu/device.h): the priority lease of a wool GPU job over the tenants that hold device
// memory outside this process's allocator (the live preview's renderer), and the live budget reading
// that decides whether a job fits.
//
// State: one LeaseState (mutex + condition variable) that is NOT the allocator's lock (runtime.cu):
// tenant calls block for up to a second or more, and the allocator must never wait on them.
// - leases: living outermost leases. While > 0, priorityHeld() is true.
// - reclaiming: one reclaim at a time; a second lease that needs memory waits for the running one
//   (reclaimRound counts the finished ones) instead of releasing the tenants again.
// - callsInFlight: calls into tenants running outside the lock; TenantRegistration's destructor waits
//   for them, so a tenant never disappears under a call.
#include "gpu/device.h"
#include "gpu/device/runtime.h"
#include "model/cancel.h"

#include <cuda_runtime.h>

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <dxgi1_4.h>
#pragma comment(lib, "dxgi.lib")
#endif

namespace knit::gpu {
namespace {

struct LeaseState {
    std::mutex m;
    std::condition_variable cv;
    int leases = 0;
    bool reclaiming = false;
    long long reclaimRound = 0;
    size_t lastFreed = 0;
    std::vector<ResidencyTenant*> tenants;
    int callsInFlight = 0;
};

LeaseState& state() {
    static LeaseState s;
    return s;
}

thread_local int tlsPriorityDepth = 0;

MemoryAvailability (*gTestReading)() = nullptr;

void throwIfCancelled() {
    const CancelToken* token = threadCancelToken();
    if (token && token->cancelled()) throw Cancelled{};
}

std::string mib(size_t bytes) {
    char b[32];
    std::snprintf(b, sizeof b, "%.0f MiB", bytes / 1048576.0);
    return b;
}

#ifdef _WIN32
// The DXGI adapter of the CUDA device (matched by LUID), found once. nullptr: no match.
IDXGIAdapter3* adapterOfTheDevice() {
    static IDXGIAdapter3* adapter = nullptr;
    static std::once_flag once;
    std::call_once(once, [] {
        const DeviceProperties& p = deviceProperties();
        if (!p.ok || !p.luidValid) return;
        IDXGIFactory1* factory = nullptr;
        if (FAILED(CreateDXGIFactory1(__uuidof(IDXGIFactory1), reinterpret_cast<void**>(&factory)))) return;
        IDXGIAdapter1* a = nullptr;
        for (UINT i = 0; factory->EnumAdapters1(i, &a) != DXGI_ERROR_NOT_FOUND; ++i) {
            DXGI_ADAPTER_DESC1 d;
            a->GetDesc1(&d);
            if (std::memcmp(&d.AdapterLuid, p.luid, sizeof(LUID)) == 0)
                a->QueryInterface(__uuidof(IDXGIAdapter3), reinterpret_cast<void**>(&adapter));
            a->Release();
            if (adapter) break;
        }
        factory->Release();
    });
    return adapter;
}
// NVML (the driver's nvml.dll, opened at run time): the device's used memory across ALL processes.
// On WDDM neither cudaMemGetInfo nor the DXGI budget reflects another process's memory -- cudaMemGetInfo
// answers this process's budget minus its own usage (measured: 11 071 MiB "free" while the preview's
// renderer held 7.6 GiB of the 12 GiB card). NVML does.
struct Nvml {
    void* device = nullptr;
    int (*getMemoryInfo)(void*, unsigned long long*) = nullptr;   // nvmlMemory_t {total, free, used}
};

const Nvml& nvml() {
    static Nvml n;
    static std::once_flag once;
    std::call_once(once, [] {
        HMODULE lib = LoadLibraryA("nvml.dll");
        if (!lib) return;
        auto init = reinterpret_cast<int (*)()>(GetProcAddress(lib, "nvmlInit_v2"));
        auto byBus = reinterpret_cast<int (*)(const char*, void**)>(GetProcAddress(lib, "nvmlDeviceGetHandleByPciBusId_v2"));
        auto info = reinterpret_cast<int (*)(void*, unsigned long long*)>(GetProcAddress(lib, "nvmlDeviceGetMemoryInfo"));
        if (!init || !byBus || !info || init() != 0) return;
        char bus[64] = {};
        if (cudaDeviceGetPCIBusId(bus, sizeof bus, 0) != cudaSuccess) {
            cudaGetLastError();  // no bus id: no NVML reading
            return;
        }
        void* device = nullptr;
        if (byBus(bus, &device) != 0) return;
        n.device = device;
        n.getMemoryInfo = info;
    });
    return n;
}

// The card's physically free memory, or false when NVML cannot answer.
bool physicallyFree(size_t& freeBytes) {
    const Nvml& n = nvml();
    if (!n.getMemoryInfo) return false;
    unsigned long long m[3] = {0, 0, 0};
    if (n.getMemoryInfo(n.device, m) != 0) return false;
    freeBytes = static_cast<size_t>(m[1]);
    return true;
}
#endif

// Releases every tenant once, or waits for the reclaim another lease is running. Returns the bytes
// that reclaim freed.
size_t reclaim(int gracefulMs) {
    LeaseState& s = state();
    std::unique_lock<std::mutex> lk(s.m);
    if (s.reclaiming) {
        const long long round = s.reclaimRound;
        while (s.reclaimRound == round) {
            s.cv.wait_for(lk, std::chrono::milliseconds(20));
            if (s.reclaimRound == round) {
                lk.unlock();
                throwIfCancelled();
                lk.lock();
            }
        }
        return s.lastFreed;
    }
    s.reclaiming = true;
    const std::vector<ResidencyTenant*> tenants = s.tenants;
    ++s.callsInFlight;
    lk.unlock();
    size_t freed = 0;
    for (ResidencyTenant* t : tenants) {
        try {
            if (t->reclaimableBytes() > 0) freed += t->release(gracefulMs);
        } catch (...) {
            // A tenant that fails to release frees nothing; the re-measure decides.
        }
    }
    lk.lock();
    s.reclaiming = false;
    --s.callsInFlight;
    ++s.reclaimRound;
    s.lastFreed = freed;
    lk.unlock();
    s.cv.notify_all();
    return freed;
}

// What the tenants are about to take on their own (loads in flight).
size_t pendingGrowth() {
    LeaseState& s = state();
    std::vector<ResidencyTenant*> tenants;
    {
        std::lock_guard<std::mutex> lk(s.m);
        tenants = s.tenants;
        ++s.callsInFlight;
    }
    size_t pending = 0;
    for (ResidencyTenant* t : tenants) pending += std::min(t->pendingBytes(), kUnknownGrowth);
    {
        std::lock_guard<std::mutex> lk(s.m);
        --s.callsInFlight;
    }
    s.cv.notify_all();
    return pending;
}

// True when `bytes` fit beside the tenants' pending growth.
bool fits(size_t bytes, size_t pending, const MemoryAvailability& a) {
    return pending <= a.availableBytes && bytes <= a.availableBytes - pending;
}

// The end of an outermost lease: the last one lets the tenants rebuild.
void endLease() {
    LeaseState& s = state();
    std::vector<ResidencyTenant*> tenants;
    {
        std::lock_guard<std::mutex> lk(s.m);
        if (--s.leases > 0) return;
        tenants = s.tenants;
        ++s.callsInFlight;
    }
    for (ResidencyTenant* t : tenants) t->mayRebuild();
    {
        std::lock_guard<std::mutex> lk(s.m);
        --s.callsInFlight;
    }
    s.cv.notify_all();
}

} // namespace

MemoryAvailability memoryAvailability() {
    if (gTestReading) return gTestReading();
    MemoryAvailability a;
    a.headroomBytes = kHeadroomBytes;
    if (!deviceProperties().ok) return a;
    size_t freeB = 0, totalB = 0;
    if (cudaMemGetInfo(&freeB, &totalB) != cudaSuccess) {
        cudaGetLastError();  // clears the sticky error; nothing is known to be available
        return a;
    }
#ifdef _WIN32
    size_t nvmlFree = 0;
    if (physicallyFree(nvmlFree)) freeB = std::min(freeB, nvmlFree);
#endif
    a.physicalFreeBytes = freeB;
    size_t limit = freeB;
#ifdef _WIN32
    if (IDXGIAdapter3* adapter = adapterOfTheDevice()) {
        DXGI_QUERY_VIDEO_MEMORY_INFO info = {};
        if (SUCCEEDED(adapter->QueryVideoMemoryInfo(0, DXGI_MEMORY_SEGMENT_GROUP_LOCAL, &info))) {
            a.budgetKnown = true;
            a.budgetBytes = static_cast<size_t>(info.Budget);
            a.usageBytes = static_cast<size_t>(info.CurrentUsage);
            const size_t left = a.budgetBytes > a.usageBytes ? a.budgetBytes - a.usageBytes : 0;
            limit = std::min(limit, left);
        }
    }
#endif
    a.availableBytes = limit > a.headroomBytes ? limit - a.headroomBytes : 0;
    return a;
}

TenantRegistration::TenantRegistration(ResidencyTenant* tenant) : tenant_(tenant) {
    LeaseState& s = state();
    std::lock_guard<std::mutex> lk(s.m);
    s.tenants.push_back(tenant);
}

TenantRegistration::~TenantRegistration() {
    LeaseState& s = state();
    std::unique_lock<std::mutex> lk(s.m);
    s.tenants.erase(std::remove(s.tenants.begin(), s.tenants.end(), tenant_), s.tenants.end());
    s.cv.wait(lk, [&] { return s.callsInFlight == 0; });
}

PriorityLease::PriorityLease(size_t additionalBytes, int gracefulMs) {
    if (tlsPriorityDepth++ > 0) {
        nested_ = true;
        return;
    }
    LeaseState& s = state();
    bool anyTenant = false;
    {
        std::lock_guard<std::mutex> lk(s.m);
        ++s.leases;
        anyTenant = !s.tenants.empty();
    }
    if (!anyTenant) return;  // nothing outside the allocator to ask: behaviour as without leases
    const bool withoutEstimate = additionalBytes == kJobBudget;
    if (withoutEstimate) additionalBytes = jobBytes();
    try {
        MemoryAvailability a = memoryAvailability();
        if (fits(additionalBytes, pendingGrowth(), a)) return;
        freed_ = reclaim(gracefulMs);
        if (withoutEstimate && freed_ == 0) return;  // nothing came back: run with what there is
        // The freed memory shows in the readings after the holder's process let go of it: read again
        // until the job fits, or report.
        const auto until = std::chrono::steady_clock::now() + std::chrono::milliseconds(kRecoveryWaitMs);
        for (;;) {
            a = memoryAvailability();
            if (fits(additionalBytes, pendingGrowth(), a)) return;
            if (std::chrono::steady_clock::now() >= until) break;
            throwIfCancelled();
            std::this_thread::sleep_for(std::chrono::milliseconds(50));
        }
        if (withoutEstimate) return;  // the tenants gave back what they held: the job's own caps decide the rest
        throw DeviceBudgetExceeded("the GPU job needs " + mib(additionalBytes) + " but " + mib(a.availableBytes) +
                                   " are available after the live preview gave back " + mib(freed_) + " (budget " +
                                   mib(a.budgetBytes) + ", in use " + mib(a.usageBytes) + ", physically free " +
                                   mib(a.physicalFreeBytes) + ", headroom " + mib(a.headroomBytes) + ")");
    } catch (...) {
        endLease();
        --tlsPriorityDepth;
        throw;
    }
}

PriorityLease::~PriorityLease() {
    --tlsPriorityDepth;
    if (!nested_) endLease();
}

bool priorityHeld() {
    LeaseState& s = state();
    std::lock_guard<std::mutex> lk(s.m);
    return s.leases > 0;
}

bool tenantsRegistered() {
    LeaseState& s = state();
    std::lock_guard<std::mutex> lk(s.m);
    return !s.tenants.empty();
}

bool threadHoldsPriority() { return tlsPriorityDepth > 0; }

size_t jobBytes() { return deviceBudget(); }

void setAvailabilityForTests(MemoryAvailability (*reading)()) { gTestReading = reading; }

} // namespace knit::gpu
