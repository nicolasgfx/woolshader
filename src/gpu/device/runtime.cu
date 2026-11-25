// GPU backend: the device service (runtime.h).
#include "gpu/device/runtime.h"

#include <algorithm>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <map>
#include <set>
#include <mutex>
#include <thread>
#include <unordered_map>

#include <cuda_runtime.h>

#include "model/cancel.h"

namespace knit::gpu {

namespace {

struct Service {
    std::mutex m;
    std::condition_variable cv;
    MemStats stats;
    size_t budget = 0;
    unsigned long long ticketNext = 0, ticketServing = 0;  // FIFO admission
    std::set<unsigned long long> abandoned;                 // tickets whose waiter was cancelled (skipped)
    std::unordered_map<void*, size_t> sizes;               // live allocations
    std::map<uint64_t, std::shared_ptr<Resident>> resident;
};
Service& svc() {
    static Service* s = new Service;  // never destroyed: allocations may outlive static teardown order
    return *s;
}

} // namespace

namespace {

// The calling thread's abort flag in host-mapped memory (UVA: the host pointer is valid on the
// device). Released at thread exit; errors then are ignored (the context may be gone).
struct AbortFlag {
    int* host = nullptr;
    AbortFlag() {
        // Handled: without the mapped flag a cancel waits for the kernel's end instead of aborting it.
        if (cudaHostAlloc(reinterpret_cast<void**>(&host), sizeof(int), cudaHostAllocMapped) != cudaSuccess) host = nullptr;
        if (host) *host = 0;
    }
    ~AbortFlag() {
        if (host) cudaFreeHost(host);  // cleanup in a destructor: never throws
    }
};
AbortFlag& abortFlag() {
    thread_local AbortFlag f;
    return f;
}

} // namespace

const volatile int* abortFlagDevice() {
    static const int never = 0;
    AbortFlag& f = abortFlag();
    if (!f.host) return &never;  // mapping failed: kernels run to completion (cancel waits)
    // Called as a launch argument: the cleared flag is written before the launch, so the kernel
    // never sees an abort left over from an earlier wait (waitForKernels clears it as well).
    *f.host = 0;
    int* d = nullptr;
    // Handled: a flag the device cannot read never aborts (as without the mapped flag).
    if (cudaHostGetDevicePointer(reinterpret_cast<void**>(&d), f.host, 0) != cudaSuccess) return &never;
    return d;
}

int abortFlagValue() {
    AbortFlag& f = abortFlag();
    return f.host ? *f.host : 0;
}

namespace {
thread_local long long tChecks = 0, tFailAt = 0;  // this thread's checks so far; the one to fail (0: none)
}

void failCheckForTests(long long n) { tFailAt = n > 0 ? tChecks + n : 0; }
long long checksForTests() { return tChecks; }

void countCheck(const char* what) {
    ++tChecks;
    if (tFailAt == 0 || tChecks != tFailAt) return;
    tFailAt = 0;
    throw DeviceError(std::string(what) + ": injected failure");
}

void check(int cudaStatus, const char* what) {
    countCheck(what);
    if (cudaStatus == cudaSuccess) return;
    throw DeviceError(std::string(what) + ": " + cudaGetErrorString(static_cast<cudaError_t>(cudaStatus)));
}

void checkLaunch(const char* what) { check(cudaGetLastError(), what); }

void waitForKernels(const char* what) { waitForStream(what, nullptr); }

void waitForStream(const char* what, CUstream_st* stream) {
    const knit::CancelToken* token = knit::threadCancelToken();
    AbortFlag& f = abortFlag();
    bool aborted = false;
    for (;;) {
        const cudaError_t e = cudaStreamQuery(stream);
        if (e == cudaSuccess) break;
        if (e != cudaErrorNotReady) {
            if (f.host) *f.host = 0;
            check(e, what);
        }
        if (!aborted && token && token->cancelled()) {
            if (f.host) *f.host = 1;
            aborted = true;
        }
        std::this_thread::sleep_for(std::chrono::microseconds(200));
    }
    if (f.host) *f.host = 0;
    if (aborted || (token && token->cancelled())) throw knit::Cancelled{};
}

bool reportFailure(const DeviceError& e, std::string* err) {
    if (err) *err = e.what();
    return false;
}

double elapsedMs(Clock::time_point since) { return std::chrono::duration<double, std::milli>(Clock::now() - since).count(); }

const DeviceProperties& deviceProperties() {
    static DeviceProperties info;
    static std::once_flag once;
    std::call_once(once, [] {
        int n = 0;
        cudaError_t e = cudaGetDeviceCount(&n);
        if (e != cudaSuccess || n == 0) {
            info.error = e != cudaSuccess ? std::string("no CUDA device: ") + cudaGetErrorString(e) : "no CUDA device";
            return;
        }
        cudaDeviceProp p;
        e = cudaGetDeviceProperties(&p, 0);
        if (e != cudaSuccess) {
            info.error = std::string("cudaGetDeviceProperties: ") + cudaGetErrorString(e);
            return;
        }
        e = cudaFree(nullptr);  // context creation (cold start) here, outside any timing
        if (e != cudaSuccess) {
            info.error = std::string("CUDA context: ") + cudaGetErrorString(e);
            return;
        }
        info.name = p.name;
        info.sm = p.major * 10 + p.minor;
        info.multiProcessors = p.multiProcessorCount;
        info.totalBytes = p.totalGlobalMem;
        std::memcpy(info.luid, p.luid, sizeof info.luid);  // meaningful on Windows only (luidDeviceNodeMask 0 elsewhere)
        info.luidValid = p.luidDeviceNodeMask != 0;
        // Informational (the device summary): a failed query leaves 0.
        if (cudaDriverGetVersion(&info.driver) != cudaSuccess) info.driver = 0;
        if (cudaRuntimeGetVersion(&info.runtime) != cudaSuccess) info.runtime = 0;
        size_t freeB = 0, totalB = 0;
        e = cudaMemGetInfo(&freeB, &totalB);  // the budget's base: without it the device is not usable
        if (e != cudaSuccess) {
            info.error = std::string("cudaMemGetInfo: ") + cudaGetErrorString(e);
            return;
        }
        size_t budget = freeB / 4 * 3;
        if (const char* env = std::getenv("KNIT_GPU_BUDGET_MB")) {
            const long long mb = std::atoll(env);
            if (mb > 0) budget = static_cast<size_t>(mb) << 20;
        }
        std::lock_guard<std::mutex> lk(svc().m);
        svc().budget = budget;
        info.ok = true;
    });
    return info;
}

size_t deviceBudget() {
    deviceProperties();
    std::lock_guard<std::mutex> lk(svc().m);
    return svc().budget;
}

MemStats memStats() {
    std::lock_guard<std::mutex> lk(svc().m);
    return svc().stats;
}

void resetPeak() {
    std::lock_guard<std::mutex> lk(svc().m);
    svc().stats.peak = svc().stats.current;
    svc().stats.peakReserved = svc().stats.reserved;
}

size_t solveCapBytes() {
    if (const char* env = std::getenv("KNIT_GPU_BUDGET_MB"))
        if (std::atoll(env) > 0) return deviceBudget();
    deviceProperties();
    size_t freeB = 0, totalB = 0;
    if (cudaMemGetInfo(&freeB, &totalB) != cudaSuccess) {  // handled: the cap falls back below
        cudaGetLastError();  // clears the sticky error
        return deviceBudget();
    }
    const size_t reserve = std::max(size_t(768) << 20, totalB / 100 * 6);
    const size_t cur = memStats().current;
    return freeB > reserve ? cur + (freeB - reserve) : cur;
}

namespace {
thread_local size_t tlCap = 0;  // AllocationCap: 0 = none
thread_local bool tlCapHit = false;
} // namespace

AllocationCap::AllocationCap(size_t bytes) : prev_(tlCap) {
    tlCap = bytes;
    tlCapHit = false;
}
AllocationCap::~AllocationCap() { tlCap = prev_; }
bool AllocationCap::capHit() { return tlCapHit; }
size_t AllocationCap::limit() { return tlCap; }

void* devAlloc(size_t bytes, std::string* err) {
    // The bytes are counted before cudaMalloc (and given back if it fails), in the same critical
    // section as the cap check: a capped check sees every tracked allocation in flight on other threads.
    {
        std::lock_guard<std::mutex> lk(svc().m);
        MemStats& s = svc().stats;
        if (tlCap > 0 && s.current + bytes > tlCap) {
            tlCapHit = true;
            if (err) {
                char b[240];
                std::snprintf(b, sizeof b, "device memory cap: %.1f MB more would raise the tracked total from %.1f MB to %.1f MB, above the cap of %.0f MB",
                              bytes / 1048576.0, s.current / 1048576.0, (s.current + bytes) / 1048576.0, tlCap / 1048576.0);
                *err = b;
            }
            return nullptr;
        }
        s.current += bytes;
        if (s.current > s.peak) s.peak = s.current;
    }
    void* p = nullptr;
    const cudaError_t e = cudaMalloc(&p, bytes);
    std::lock_guard<std::mutex> lk(svc().m);
    if (e != cudaSuccess) {
        svc().stats.current -= bytes;
        cudaGetLastError();  // clear the sticky-free error state of a failed allocation
        if (err) {
            char b[200];
            std::snprintf(b, sizeof b, "cudaMalloc(%zu bytes): %s", bytes, cudaGetErrorString(e));
            *err = b;
        }
        return nullptr;
    }
    ++svc().stats.allocations;
    svc().sizes[p] = bytes;
    return p;
}

void devFree(void* p) {
    if (!p) return;
    cudaFree(p);  // cleanup (destructors call this): never throws
    std::lock_guard<std::mutex> lk(svc().m);
    auto it = svc().sizes.find(p);
    if (it != svc().sizes.end()) {
        svc().stats.current -= it->second;
        svc().sizes.erase(it);
    }
}

void DeviceBuffer::alloc(size_t bytes) {
    reset();
    std::string why;
    p_ = devAlloc(bytes ? bytes : 1, &why);
    if (!p_) {
        if (AllocationCap::capHit()) throw DeviceBudgetExceeded(why);
        throw DeviceError(why);
    }
    n_ = bytes;
}

bool DeviceBuffer::reserve(size_t bytes) {
    if (p_ && n_ >= bytes) return false;
    alloc(std::max<size_t>(bytes + bytes / 4, 64));
    return true;
}

void DeviceBuffer::upload(const void* host, size_t bytes) {
    alloc(bytes);
    if (bytes == 0) return;
    check(cudaMemcpy(p_, host, bytes, cudaMemcpyHostToDevice), "upload");
    std::lock_guard<std::mutex> lk(svc().m);
    svc().stats.uploadedBytes += static_cast<long long>(bytes);
}

void DeviceBuffer::upload(const std::vector<int>& v) { upload(v.data(), v.size() * sizeof(int)); }
void DeviceBuffer::upload(const std::vector<unsigned>& v) { upload(v.data(), v.size() * sizeof(unsigned)); }
void DeviceBuffer::upload(const std::vector<long long>& v) { upload(v.data(), v.size() * sizeof(long long)); }
void DeviceBuffer::upload(const std::vector<float>& v) { upload(v.data(), v.size() * sizeof(float)); }
void DeviceBuffer::upload(const std::vector<double>& v) { upload(v.data(), v.size() * sizeof(double)); }
void DeviceBuffer::upload(const std::vector<unsigned char>& v) { upload(v.data(), v.size()); }

void DeviceBuffer::download(void* host, size_t bytes, size_t offset) const {
    if (bytes == 0) return;
    check(cudaMemcpy(host, static_cast<const char*>(p_) + offset, bytes, cudaMemcpyDeviceToHost), "download");
    std::lock_guard<std::mutex> lk(svc().m);
    svc().stats.downloadedBytes += static_cast<long long>(bytes);
}
void DeviceBuffer::download(std::vector<int>& v, size_t count) const {
    v.resize(count);
    download(v.data(), count * sizeof(int));
}
void DeviceBuffer::download(std::vector<unsigned>& v, size_t count) const {
    v.resize(count);
    download(v.data(), count * sizeof(unsigned));
}
void DeviceBuffer::download(std::vector<float>& v, size_t count) const {
    v.resize(count);
    download(v.data(), count * sizeof(float));
}
void DeviceBuffer::download(std::vector<double>& v, size_t count) const {
    v.resize(count);
    download(v.data(), count * sizeof(double));
}
void DeviceBuffer::download(std::vector<unsigned char>& v, size_t count) const {
    v.resize(count);
    download(v.data(), count);
}

namespace {
thread_local int tlsLeaseDepth = 0;  // reservations held by this thread

// Advances the served ticket past cancelled waiters (s.m held).
void skipAbandoned(Service& s) {
    for (auto it = s.abandoned.find(s.ticketServing); it != s.abandoned.end(); it = s.abandoned.find(s.ticketServing)) {
        s.abandoned.erase(it);
        ++s.ticketServing;
    }
}
} // namespace

int leaseDepth() { return tlsLeaseDepth; }

size_t abandonedTickets() {
    std::lock_guard<std::mutex> lk(svc().m);
    return svc().abandoned.size();
}

size_t queuedRequests() {
    std::lock_guard<std::mutex> lk(svc().m);
    return static_cast<size_t>(svc().ticketNext - svc().ticketServing) - svc().abandoned.size();
}

size_t Reservation::exclusiveBytes() { return std::max<size_t>(deviceBudget(), 1); }

Reservation::Reservation(size_t bytes) : bytes_(bytes) {
    // Before the ticket and outside the allocator's lock: a lease may wait for the tenants.
    // A request that runs alone (exclusiveBytes) has no size of its own: kJobBudget.
    if (tlsLeaseDepth == 0 && !threadHoldsPriority() && tenantsRegistered())
        implicitLease_ = new PriorityLease(bytes >= exclusiveBytes() ? kJobBudget : bytes);
    if (tlsLeaseDepth++ > 0) {
        nested_ = true;
        return;
    }
    const auto t0 = Clock::now();
    Service& s = svc();
    deviceProperties();
    std::unique_lock<std::mutex> lk(s.m);
    const unsigned long long ticket = s.ticketNext++;
    bool waited = false;
    // FIFO: wait for the turn, then for room (a request above the budget waits for an empty device).
    // A cancelled job (cancel.h, the thread's token) leaves the queue: its ticket is skipped, it
    // holds nothing, the thread's lease depth is restored, the others are woken; then Cancelled.
    const knit::CancelToken* token = knit::threadCancelToken();
    auto ready = [&] {
        const bool turn = s.ticketServing == ticket;
        const size_t avail = s.budget > s.stats.reserved + s.stats.resident ? s.budget - s.stats.reserved - s.stats.resident : 0;
        const bool room = bytes <= avail || s.stats.reserved == 0;
        if (!(turn && room)) waited = true;
        return turn && room;
    };
    while (!ready()) {
        if (token && token->cancelled()) {
            s.abandoned.insert(ticket);
            skipAbandoned(s);
            lk.unlock();
            --tlsLeaseDepth;
            s.cv.notify_all();
            delete implicitLease_;
            implicitLease_ = nullptr;
            throw knit::Cancelled{};
        }
        if (token) s.cv.wait_for(lk, std::chrono::milliseconds(20));
        else s.cv.wait(lk);
    }
    ++s.ticketServing;
    skipAbandoned(s);
    s.stats.reserved += bytes;
    if (s.stats.reserved > s.stats.peakReserved) s.stats.peakReserved = s.stats.reserved;
    if (waited) ++s.stats.waits;
    held_ = true;
    lk.unlock();
    s.cv.notify_all();
    waitedMs_ = elapsedMs(t0);
}

Reservation::~Reservation() {
    --tlsLeaseDepth;
    if (held_) {
        Service& s = svc();
        {
            std::lock_guard<std::mutex> lk(s.m);
            s.stats.reserved -= bytes_;
        }
        s.cv.notify_all();
    }
    delete implicitLease_;
}

std::shared_ptr<Resident> residentGet(uint64_t key) {
    std::lock_guard<std::mutex> lk(svc().m);
    auto it = svc().resident.find(key);
    return it == svc().resident.end() ? nullptr : it->second;
}

void residentPut(uint64_t key, std::shared_ptr<Resident> r) {
    std::lock_guard<std::mutex> lk(svc().m);
    std::shared_ptr<Resident>& slot = svc().resident[key];
    if (slot) svc().stats.resident -= slot->bytes;
    slot = std::move(r);
    svc().stats.resident += slot->bytes;
}

void releaseResident() {
    std::map<uint64_t, std::shared_ptr<Resident>> drop;
    {
        std::lock_guard<std::mutex> lk(svc().m);
        drop.swap(svc().resident);
        svc().stats.resident = 0;
    }
    drop.clear();  // frees outside the lock (devFree takes it)
    svc().cv.notify_all();
}

LaunchPacer::LaunchPacer(long long start, long long minN, long long maxN)
    : n_(std::clamp(start, std::max(1ll, minN), std::max(1ll, maxN))), minN_(std::max(1ll, minN)), maxN_(std::max(1ll, maxN)), smallest_(n_) {
    target_ = 150;
    if (const char* e = std::getenv("KNIT_GPU_LAUNCH_MS")) {
        const double v = std::atof(e);
        if (v > 0) target_ = v;
    }
}

void LaunchPacer::measured(long long n, double ms) {
    slowest_ = std::max(slowest_, ms);
    if (n <= 0) return;
    const double perItem = std::max(ms, 1e-3) / static_cast<double>(n);
    const double want = target_ / perItem;
    const double cap = static_cast<double>(n_) * 4;
    n_ = std::clamp(static_cast<long long>(std::min(want, cap)), minN_, maxN_);
    smallest_ = std::min(smallest_, n_);
}

bool available(std::string* why) {
    const DeviceProperties& d = deviceProperties();
    if (!d.ok && why) *why = d.error;
    return d.ok;
}

std::string deviceSummary() {
    const DeviceProperties& d = deviceProperties();
    if (!d.ok) return "no usable CUDA device: " + d.error;
    char b[256];
    std::snprintf(b, sizeof b, "%s (sm_%d, %d SMs, %.1f GB), CUDA runtime %d.%d, driver %d.%d, float32 geometry, budget %.0f MB",
                  d.name.c_str(), d.sm, d.multiProcessors, d.totalBytes / 1073741824.0, d.runtime / 1000, d.runtime % 1000 / 10,
                  d.driver / 1000, d.driver % 1000 / 10, deviceBudget() / 1048576.0);
    return b;
}

} // namespace knit::gpu
