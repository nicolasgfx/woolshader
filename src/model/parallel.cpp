// The process-wide executor with its CPU budget of H slots, and parallelFor (parallel.h).
#include "model/parallel.h"

#include <algorithm>
#include <condition_variable>
#include <exception>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

#include "model/cancel.h"

namespace knit {

namespace {

// One parallelFor call: its chunks, how many started / finished / are running.
struct Range {
    const std::function<void(int64_t, int64_t, int64_t)>* fn = nullptr;
    int64_t n = 0, chunk = 1, chunks = 0;
    int64_t next = 0;      // next chunk to start
    int64_t finished = 0;  // finished (or cancelled) chunks
    int running = 0;
    int maxConcurrency = 0;
    std::exception_ptr error;
    int64_t errorChunk = std::numeric_limits<int64_t>::max();
    bool startable() const { return next < chunks && (maxConcurrency <= 0 || running < maxConcurrency); }
    bool done() const { return finished == chunks; }
};

// This thread holds one of the H slots. Read and written only by its own thread, while it holds the executor's m_.
thread_local bool tHoldsSlot = false;

class Executor {
public:
    explicit Executor(int h) : slotsFree_(h) {
        for (int i = 0; i < h; ++i) workers_.emplace_back([this] { workerLoop(); });
    }
    ~Executor() {
        {
            std::lock_guard<std::mutex> lk(m_);
            stop_ = true;
        }
        cv_.notify_all();
        for (std::thread& t : workers_) t.join();
    }

    void run(Range& r) {
        std::unique_lock<std::mutex> lk(m_);
        if (stop_) throw std::logic_error("parallelFor after executor shutdown");
        open_.push_back(&r);
        cv_.notify_all();
        const bool hadSlot = tHoldsSlot;
        bool holding = hadSlot;
        while (!r.done()) {
            Range* pick = r.startable() ? &r : firstStartable();
            if (pick && (holding || slotsFree_ > 0)) {
                if (!holding) {
                    --slotsFree_;
                    holding = tHoldsSlot = true;
                }
                runChunk(*pick, lk);
                continue;
            }
            // Nothing this thread can start: give the slot back while waiting.
            if (holding) {
                ++slotsFree_;
                holding = tHoldsSlot = false;
                cv_.notify_all();
            }
            cv_.wait(lk);
        }
        // Return holding exactly what the caller held.
        if (hadSlot && !holding) {
            cv_.wait(lk, [&] { return slotsFree_ > 0; });
            --slotsFree_;
        } else if (!hadSlot && holding) {
            ++slotsFree_;
            cv_.notify_all();
        }
        tHoldsSlot = hadSlot;
        lk.unlock();
        if (r.error) std::rethrow_exception(r.error);
    }

private:
    std::mutex m_;
    std::condition_variable cv_;
    std::vector<Range*> open_;  // ranges that may still have chunks to start, in submission order
    int slotsFree_;
    bool stop_ = false;
    std::vector<std::thread> workers_;

    Range* firstStartable() {
        for (Range* r : open_)
            if (r->startable()) return r;
        return nullptr;
    }

    // Runs one chunk of r (the lock is held on entry and on return, released while running).
    void runChunk(Range& r, std::unique_lock<std::mutex>& lk) {
        const int64_t k = r.next++;
        ++r.running;
        if (r.next == r.chunks) open_.erase(std::find(open_.begin(), open_.end(), &r));
        lk.unlock();
        std::exception_ptr err;
        try {
            const int64_t b = k * r.chunk, e = std::min(r.n, b + r.chunk);
            (*r.fn)(b, e, k);
        } catch (...) {
            err = std::current_exception();
        }
        lk.lock();
        --r.running;
        ++r.finished;
        if (err && k < r.errorChunk) {
            r.error = err;
            r.errorChunk = k;
        }
        if (err && r.next < r.chunks) {  // cancel the chunks that have not started
            r.finished += r.chunks - r.next;
            r.next = r.chunks;
            open_.erase(std::find(open_.begin(), open_.end(), &r));
        }
        cv_.notify_all();
    }

    void workerLoop() {
        std::unique_lock<std::mutex> lk(m_);
        for (;;) {
            cv_.wait(lk, [&] { return stop_ || (slotsFree_ > 0 && firstStartable()); });
            if (stop_) return;
            --slotsFree_;
            tHoldsSlot = true;
            while (Range* r = firstStartable()) runChunk(*r, lk);
            tHoldsSlot = false;
            ++slotsFree_;
            cv_.notify_all();
        }
    }
};

// The process-wide executor. gExecMutex guards gBudget and gExec. setThreadBudget sets gBudget only while gExec does
// not exist (a later call is ignored); the first parallel call sets it to the hardware threads if nobody did, then
// creates gExec, which lives until static destruction: its destructor joins the workers.
std::mutex gExecMutex;
int gBudget = 0;
std::unique_ptr<Executor> gExec;

Executor& executor() {
    std::lock_guard<std::mutex> lk(gExecMutex);
    if (!gExec) {
        if (gBudget <= 0) gBudget = std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
        gExec = std::make_unique<Executor>(gBudget);
    }
    return *gExec;
}

} // namespace

void setThreadBudget(int h) {
    std::lock_guard<std::mutex> lk(gExecMutex);
    if (!gExec && h > 0) gBudget = h;
}

int threadBudget() {
    std::lock_guard<std::mutex> lk(gExecMutex);
    return gBudget > 0 ? gBudget : std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
}

void parallelFor(int64_t n, int64_t chunk, const std::function<void(int64_t, int64_t, int64_t)>& fn, const ExecPolicy& policy) {
    if (n <= 0) return;
    chunk = std::max<int64_t>(1, chunk);
    const int64_t chunks = (n + chunk - 1) / chunk;
    if (!policy.parallel || chunks == 1) {
        for (int64_t k = 0; k < chunks; ++k) fn(k * chunk, std::min(n, (k + 1) * chunk), k);
        return;
    }
    // The caller's cancellation token (cancel.h) holds inside every chunk of this range,
    // on whichever thread runs it; always set (also to none), since a waiting caller runs
    // other ranges' chunks too.
    const CancelToken* token = threadCancelToken();
    const std::function<void(int64_t, int64_t, int64_t)> scoped = [&fn, token](int64_t b, int64_t e, int64_t k) {
        ThreadCancelScope scope(token);
        fn(b, e, k);
    };
    Range r;
    r.fn = &scoped;
    r.n = n;
    r.chunk = chunk;
    r.chunks = chunks;
    r.maxConcurrency = policy.maxConcurrency;
    executor().run(r);
}

} // namespace knit
