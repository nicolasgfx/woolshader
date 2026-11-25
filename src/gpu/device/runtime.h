#pragma once
// GPU backend: the device service (one per process). Every device allocation goes through
// it (current / peak bytes); concurrent jobs (the viewer's workers)
// reserve their working set first, and a reservation that does not fit the budget waits
// until earlier jobs release theirs (FIFO), so N jobs never hold N uncontrolled uploads.
// Resident data (meshes shared by several pieces) is kept across jobs, counted in the budget.
// No knit types: included by the CUDA translation units and the host adapters.
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "gpu/device.h"

struct CUstream_st;  // cudaStream_t, without the CUDA headers

namespace knit::gpu {

// Every CUDA call's status goes through check(), unless
// it is handled explicitly or is cleanup that must not throw (architecture_gpu_cuda_calls_checked).
void check(int cudaStatus, const char* what);  // throws DeviceError("what: reason") unless cudaSuccess
void checkLaunch(const char* what);            // the last launch's error, if any
// Waits for the device. When the calling thread's job is cancelled it raises the abort flag, waits for
// the kernels to stop and throws Cancelled: nothing they wrote may be used.
void waitForKernels(const char* what);
// The same for the work of one stream (a graph launched on it).
void waitForStream(const char* what, CUstream_st* stream);
bool reportFailure(const DeviceError& e, std::string* err);  // writes *err, returns false

using Clock = std::chrono::steady_clock;
double elapsedMs(Clock::time_point since);

struct DeviceProperties {
    bool ok = false;
    std::string name, error;
    int sm = 0, multiProcessors = 0, driver = 0, runtime = 0;
    size_t totalBytes = 0;
    // The adapter's LUID (Windows), so DXGI can read this process's budget on the same device.
    unsigned char luid[8] = {};
    bool luidValid = false;
};
// Initializes the CUDA context once (thread-safe); later calls return the same result.
const DeviceProperties& deviceProperties();

// Budget [bytes] for reservations: KNIT_GPU_BUDGET_MB, else 75 % of the free memory at
// initialization.
size_t deviceBudget();

struct MemStats {
    size_t current = 0, peak = 0, reserved = 0, peakReserved = 0, resident = 0;
    long long allocations = 0, waits = 0;
    long long uploadedBytes = 0, downloadedBytes = 0;  // through DeviceBuffer
};
MemStats memStats();
void resetPeak();

// Cooperative abort of long kernels (cancellation). A kernel whose per-thread work can run for
// seconds (adaptive subdivision on degenerate spans) polls the calling thread's flag, which lives
// in host-mapped memory, and returns early when it is raised (waitForKernels).
// The flag is host-mapped pinned memory (cudaHostAllocMapped, not write-combined): host stores
// before a launch are visible to it, and a store made while it runs becomes visible to its
// volatile polls (a late view only delays the abort). abortFlagDevice also clears the flag.
const volatile int* abortFlagDevice();
int abortFlagValue();  // tests: the calling thread's flag as the host sees it (0: none pending)

// The memory cap of one solve (AllocationCap): KNIT_GPU_BUDGET_MB when set,
// else what the process already tracks plus the device's free memory now, minus a reserve for
// untracked memory (the larger of 768 MB and 6 % of the device: CUDA context growth, library
// workspaces, the display). Free memory is read at the call, so other processes' use counts.
size_t solveCapBytes();

// Tracked device allocation (throws nothing; nullptr with `err` on failure).
void* devAlloc(size_t bytes, std::string* err);
void devFree(void* p);

// A cap on the tracked device memory of the process, checked by this thread's allocations while
// the cap lives: an allocation of this thread that would raise the tracked total (all threads,
// allocations in flight included) above `bytes` fails with a message instead of overcommitting
// the device. On Windows (WDDM) an overcommitted GPU pages to system memory and keeps running
// 10-30x slower.
// Scope: it limits this thread's allocations only. Other threads' allocations are counted but not
// refused, and untracked memory (CUDA context, library workspaces, other processes) is not counted;
// the pipelines rely on the knitting slot and the exclusive lease for the rest. Caps nest; the
// innermost applies. Constructing one clears the thread's capHit flag.
class AllocationCap {
public:
    explicit AllocationCap(size_t bytes);
    ~AllocationCap();
    AllocationCap(const AllocationCap&) = delete;
    AllocationCap& operator=(const AllocationCap&) = delete;
    // The calling thread's last refused allocation (since its innermost cap was made) hit the cap.
    static bool capHit();
    // The calling thread's innermost cap [bytes]; 0: none (a solve sizes its capacity to fit it).
    static size_t limit();

private:
    size_t prev_;
};

// A reservation of `bytes` of the budget, released on destruction. A request larger than
// the whole budget is granted alone (it waits until nothing else is reserved). Reservations
// nest per thread: one taken while the same thread holds another is a no-op (the outer one
// covers the call), so an entry point calling another cannot deadlock on itself.
// A waiting request whose thread's job is cancelled (cancel.h) leaves the FIFO queue and throws
// knit::Cancelled: it holds nothing, its ticket is skipped for the requests behind it, and the
// thread's nesting depth is as before the request.
class Reservation {
public:
    Reservation() = default;
    explicit Reservation(size_t bytes);
    ~Reservation();
    Reservation(const Reservation&) = delete;
    Reservation& operator=(const Reservation&) = delete;
    double waitedMs() const { return waitedMs_; }
    bool nested() const { return nested_; }

    // The whole budget: parts without a working-set estimate (validation, surface model,
    // gate (d), tube mesh, plies, the tiled lattice and boundary checks) run alone on the
    // device, so N jobs stay within the budget (plus resident data) at the cost of waiting.
    static size_t exclusiveBytes();

private:
    size_t bytes_ = 0;
    bool held_ = false, nested_ = false;
    double waitedMs_ = 0;
    // The outermost reservation of a thread that holds no PriorityLease takes one while a tenant is
    // registered (gpu/device.h): a request, oversized or not, is then admitted only after the tenants
    // gave back what it needs. Without tenants nothing changes.
    PriorityLease* implicitLease_ = nullptr;
};

// True while the calling thread holds a PriorityLease (residency.cpp).
bool threadHoldsPriority();

// What every checker does first: counts the check, and throws the failure a test injected there.
void countCheck(const char* what);
// Tests: the calling thread's n-th next check throws DeviceError("<what>: injected failure"); 0: off.
void failCheckForTests(long long n);
long long checksForTests();  // tests: this thread's checks so far

// Tests: the calling thread's reservation nesting depth; cancelled tickets not yet skipped.
int leaseDepth();
size_t abandonedTickets();
size_t queuedRequests();  // waiting (not cancelled) requests

// Launch pacing. On a GPU that drives a display, a kernel launch that runs longer than the
// driver's watchdog (Windows TDR, 2 s by default) resets the device; the process sees
// "unspecified launch failure". Batched kernels therefore size each batch from the measured
// time of the previous one: target KNIT_GPU_LAUNCH_MS (default 150 ms), growth at most x4 per
// batch, shrinking at once. Results never depend on the batch partition (outputs are per item,
// compacted in item order or sorted by key).
class LaunchPacer {
public:
    LaunchPacer(long long start, long long minN, long long maxN);
    long long next() const { return n_; }
    void measured(long long n, double ms);
    long long smallest() const { return smallest_; }
    double slowestMs() const { return slowest_; }

private:
    long long n_, minN_, maxN_, smallest_;
    double target_, slowest_ = 0;
};

// A device buffer owned by RAII (tracked allocation).
class DeviceBuffer {
public:
    DeviceBuffer() = default;
    ~DeviceBuffer() { reset(); }
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    DeviceBuffer(DeviceBuffer&& o) noexcept : p_(o.p_), n_(o.n_) { o.p_ = nullptr; o.n_ = 0; }
    DeviceBuffer& operator=(DeviceBuffer&& o) noexcept {
        if (this != &o) { reset(); p_ = o.p_; n_ = o.n_; o.p_ = nullptr; o.n_ = 0; }
        return *this;
    }
    void alloc(size_t bytes);                     // zero bytes are allowed
    bool reserve(size_t bytes);                   // keeps at least bytes (grows by 25 %); true when it reallocated
    // Allocate and copy from the host (synchronous).
    void upload(const void* host, size_t bytes);
    void upload(const std::vector<int>& v);
    void upload(const std::vector<unsigned>& v);
    void upload(const std::vector<long long>& v);
    void upload(const std::vector<float>& v);
    void upload(const std::vector<double>& v);
    void upload(const std::vector<unsigned char>& v);
    // Copy to the host (synchronous).
    void download(void* host, size_t bytes, size_t offset = 0) const;
    void download(std::vector<int>& v, size_t count) const;
    void download(std::vector<unsigned>& v, size_t count) const;
    void download(std::vector<float>& v, size_t count) const;
    void download(std::vector<double>& v, size_t count) const;
    void download(std::vector<unsigned char>& v, size_t count) const;
    void reset() {
        if (p_) devFree(p_);
        p_ = nullptr;
        n_ = 0;
    }
    // The contents as an array of plain numbers. A buffer of records is cast once, where its view is filled.
    void* data() const { return p_; }
    int* ints() const { return static_cast<int*>(p_); }
    unsigned* uints() const { return static_cast<unsigned*>(p_); }
    long long* longs() const { return static_cast<long long*>(p_); }
    unsigned long long* ulongs() const { return static_cast<unsigned long long*>(p_); }
    float* floats() const { return static_cast<float*>(p_); }
    double* doubles() const { return static_cast<double*>(p_); }
    unsigned char* chars() const { return static_cast<unsigned char*>(p_); }
    size_t bytes() const { return n_; }

private:
    void* p_ = nullptr;
    size_t n_ = 0;
};

// Resident data shared across jobs, keyed by a 64-bit content digest: built once by `make`,
// kept until releaseResident(). Counted in memStats().resident.
struct Resident {
    virtual ~Resident() = default;
    size_t bytes = 0;
};
std::shared_ptr<Resident> residentGet(uint64_t key);
void residentPut(uint64_t key, std::shared_ptr<Resident> r);

} // namespace knit::gpu
