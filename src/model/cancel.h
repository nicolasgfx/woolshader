#pragma once
// Cooperative cancellation of long computations. A worker thread installs its token
// (ThreadCancelScope); long loops call cancelCheckpoint(), which throws Cancelled once the token is
// cancelled. Without a token (the command line, the tests) a checkpoint does nothing. The token is
// read atomically, so checkpoints are safe in parallel loops. A cancel takes effect at the next
// checkpoint, so the longest interval between two checkpoints of a thread is the worst wait.
#include <atomic>

namespace knit {

struct Cancelled {};

class CancelToken {
public:
    void cancel() { flag_.store(true, std::memory_order_relaxed); }
    bool cancelled() const { return flag_.load(std::memory_order_relaxed); }

private:
    std::atomic<bool> flag_{false};
};

// The token of the current thread (nullptr: none).
inline const CancelToken*& threadCancelToken() {
    thread_local const CancelToken* token = nullptr;
    return token;
}

class ThreadCancelScope {
public:
    explicit ThreadCancelScope(const CancelToken* token) : previous_(threadCancelToken()) { threadCancelToken() = token; }
    ~ThreadCancelScope() { threadCancelToken() = previous_; }

private:
    const CancelToken* previous_;
};

inline void cancelCheckpoint() {
    const CancelToken* token = threadCancelToken();
    if (token && token->cancelled()) throw Cancelled{};
}

} // namespace knit
