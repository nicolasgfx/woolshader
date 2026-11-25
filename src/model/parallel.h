#pragma once
// One process-wide executor with a CPU budget of H slots.
//
// Every thread that executes work holds a slot, so at most H chunks run at once:
// executor workers and external callers (a thread calling parallelFor from outside the
// executor first acquires a slot). parallelFor is SYNCHRONOUS: it returns only after every
// chunk of its range has finished, so whatever the chunks capture outlives them. While its
// range is incomplete the calling thread runs queued chunks itself (its own range first,
// then others), and gives its slot back only while it has nothing to run, so nested calls
// cannot deadlock and a held slot is never idle. The first exception (lowest chunk index)
// cancels the range's chunks that have not started; running ones finish; it is rethrown on
// the caller, and the executor stays usable. The caller's cancellation token (cancel.h)
// holds inside every chunk of its range.
//
// Parallel execution is always opt-in (ExecPolicy::parallel, default false = serial on
// the calling thread). Passing a parallel policy to an API declares that every callback the
// call receives or stores may be called concurrently (no shared mutable state).
#include <cstdint>
#include <functional>

namespace knit {

struct ExecPolicy {
    bool parallel = false;   // false: serial, on the calling thread, in chunk order
    int maxConcurrency = 0;  // > 0: at most this many chunks of the range run at once
    int64_t chunk = 0;       // items per chunk (0: the API's default)
};

// The budget H (default: hardware threads). Set before the first parallel call; a later
// call is ignored once the executor has started.
void setThreadBudget(int h);
int threadBudget();

// Runs fn(begin, end, chunkIndex) for the chunks [k chunk, min((k + 1) chunk, n)) of [0, n),
// k = 0 .. ceil(n / chunk) - 1, and returns when all have finished. Serial policy: in chunk
// order on the calling thread. n <= 0: returns at once.
void parallelFor(int64_t n, int64_t chunk, const std::function<void(int64_t begin, int64_t end, int64_t chunkIndex)>& fn,
                 const ExecPolicy& policy);

} // namespace knit
