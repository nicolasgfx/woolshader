#pragma once
// The compute device as the rest of the program sees it: is it usable, what is it, and the failures a
// caller can meet. Every entry point of the backend (gpu/knitting.h, gpu/surfaces.h,
// gpu/fibre_generator.h) reports a missing device instead of silently running the CPU code. CUDA-free
// and free of knit types: device/ includes it; the way into this folder is ARCHITECTURE.md.
#include <stdexcept>
#include <string>

namespace knit::gpu {

// A device failure is not control flow: the code below the public headers throws, every buffer is RAII,
// and each public entry point turns the exception into false + err. Cancelled (model/cancel.h) passes
// through.
class DeviceError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};
// The allocation would overcommit the device (AllocationCap, solveCapBytes). The controls are
// unchanged. The tiled pipeline leaves such a part unknitted with this message; elsewhere it fails the run.
class DeviceBudgetExceeded : public DeviceError {
public:
    using DeviceError::DeviceError;
};

// A usable CUDA device; the reason otherwise.
[[nodiscard]] bool available(std::string* why = nullptr);
// "<device> (sm_89, 46 SMs, 12.0 GB), CUDA runtime 13.4, driver 13.3, float32" or the reason.
std::string deviceSummary();

// Frees resident data (meshes shared by several pieces) that no running job holds. The job queue calls it when a job
// starts on another scene.
void releaseResident();

// --- Residency: device memory held outside this process's allocator (the live render preview) ---
//
// The live preview's renderer runs in its own process and can hold most of the card. A wool GPU job has
// priority over it: before the job runs, the job takes a PriorityLease, which asks every registered
// ResidencyTenant to give its memory back when the job's bytes do not fit, and which keeps every tenant
// from loading, rebuilding or resuming until the last lease ends.

// How much device memory THIS process may still use. Windows (WDDM) decides residency against each
// process's budget, not against the card's free memory: an allocation past the budget is demoted to
// system memory although the card has room. The usable figure is therefore the smaller of the budget
// left (DXGI Budget - CurrentUsage, read live: budgets change with the other processes) and the
// physically free memory of the card (all processes), minus a headroom.
struct MemoryAvailability {
    bool budgetKnown = false;         // DXGI answered (Windows, the CUDA device's adapter)
    size_t budgetBytes = 0;           // DXGI Budget of this process
    size_t usageBytes = 0;            // DXGI CurrentUsage of this process
    size_t physicalFreeBytes = 0;     // the card's free memory across all processes (NVML; on WDDM
                                      // cudaMemGetInfo answers this process's budget instead)
    size_t headroomBytes = 0;         // kept back (kHeadroomBytes)
    size_t availableBytes = 0;        // min(budget - usage, physical free) - headroom, 0 at worst
};
constexpr size_t kHeadroomBytes = size_t(512) << 20;
MemoryAvailability memoryAvailability();

// A holder of device memory outside the allocator: the live preview (its renderer process). The
// methods are called from the thread that takes a lease, never under an allocator lock.
class ResidencyTenant {
public:
    virtual ~ResidencyTenant() = default;
    // Bytes release() can really free now (0 when nothing is held).
    virtual size_t reclaimableBytes() = 0;
    // Bytes the tenant is about to take on its own: a scene load in flight keeps allocating after a lease
    // was granted, so the lease counts it against the availability (kUnknownGrowth: not forecast).
    virtual size_t pendingBytes() { return 0; }
    // Gives the memory back and returns when it is confirmed free (the holder said so, or its process
    // was terminated and has exited). `gracefulMs`: how long the holder may take before it is
    // terminated. Returns the bytes freed.
    virtual size_t release(int gracefulMs) = 0;
    // The last lease ended: the tenant may load and rebuild again.
    virtual void mayRebuild() = 0;
};

// Registers a tenant for its lifetime. The destructor waits for calls into the tenant in flight.
class TenantRegistration {
public:
    explicit TenantRegistration(ResidencyTenant* tenant);
    ~TenantRegistration();
    TenantRegistration(const TenantRegistration&) = delete;
    TenantRegistration& operator=(const TenantRegistration&) = delete;

private:
    ResidencyTenant* tenant_;
};

// pendingBytes() of a load whose size is not known: the lease releases the tenant.
constexpr size_t kUnknownGrowth = size_t(1) << 40;

// How long a tenant may take to give its memory back gracefully before it is terminated (F3: the
// renderer releases in ~0.1 s between samples, but a scene load cannot be interrupted).
constexpr int kGracefulReleaseMs = 1000;
// How long a lease waits, after the tenants returned, for the freed memory to show in the readings.
constexpr int kRecoveryWaitMs = 3000;

// A wool GPU job's claim, held for the whole job. While any lease lives, priorityHeld() is true and the
// tenants must not load, rebuild or resume. Construction: if `additionalBytes` do not fit in
// memoryAvailability(), every tenant is released (one reclaim at a time; a second lease waits for the
// running one), the availability is read again until it fits or kRecoveryWaitMs pass, and if it still
// does not fit the lease throws DeviceBudgetExceeded with the numbers (the job reports it, it does not
// wait) -- except for kJobBudget (see there). A waiting lease whose job is cancelled throws knit::Cancelled. Leases nest per thread: an inner
// one is a no-op.
class PriorityLease {
public:
    explicit PriorityLease(size_t additionalBytes, int gracefulMs = kGracefulReleaseMs);
    ~PriorityLease();
    PriorityLease(const PriorityLease&) = delete;
    PriorityLease& operator=(const PriorityLease&) = delete;
    size_t freedBytes() const { return freed_; }   // what the tenants gave back for this lease
    bool nested() const { return nested_; }

private:
    bool nested_ = false;
    size_t freed_ = 0;
};

// True while any PriorityLease lives (in any thread).
bool priorityHeld();
// True while a tenant is registered.
bool tenantsRegistered();
// The allocator's budget (deviceBudget): what a job without an estimate of its own is measured against.
size_t jobBytes();
// As `additionalBytes`, for a job without an estimate (every SceneBuilder job, and an allocator request
// that runs alone on the device): the tenants are released when jobBytes() does not fit, and then the job
// runs with what there is, as it would without a preview -- such a lease never reports
// DeviceBudgetExceeded (jobBytes is a scheduling budget taken at start, not memory the job must have).
constexpr size_t kJobBudget = ~size_t(0);

// Tests: replaces the readings behind memoryAvailability() (nullptr: the real ones).
void setAvailabilityForTests(MemoryAvailability (*reading)());

} // namespace knit::gpu
