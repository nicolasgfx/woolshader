#pragma once
// A SERIAL stand-in for the Threading Building Blocks subset that RIM (third_party/rim) uses. Every construct runs on the
// calling thread in index order, so the remesher's output is deterministic and independent of the
// machine; the locks are no-ops (nothing runs concurrently).
#include <algorithm>
#include <cstddef>
#include <functional>
#include <queue>
#include <vector>

namespace tbb {

template <typename T>
class blocked_range {
public:
    using const_iterator = T;
    blocked_range(T begin, T end, std::size_t grain = 1) : b_(begin), e_(end), g_(grain) {}
    T begin() const { return b_; }
    T end() const { return e_; }
    std::size_t size() const { return static_cast<std::size_t>(e_ - b_); }
    bool empty() const { return !(b_ < e_); }
    std::size_t grainsize() const { return g_; }

private:
    T b_, e_;
    std::size_t g_;
};

template <typename Range, typename Body>
void parallel_for(const Range& r, const Body& body) {
    body(r);
}
template <typename Index, typename Function>
void parallel_for(Index first, Index last, const Function& f) {
    for (Index i = first; i < last; ++i) f(i);
}

template <typename Range, typename Value, typename RealBody, typename Reduction>
Value parallel_reduce(const Range& r, const Value& identity, const RealBody& body, const Reduction&) {
    return body(r, identity);
}
template <typename Range, typename Value, typename RealBody, typename Reduction>
Value parallel_deterministic_reduce(const Range& r, const Value& identity, const RealBody& body, const Reduction&) {
    return body(r, identity);
}

template <typename It>
void parallel_sort(It first, It last) {
    std::sort(first, last);
}
template <typename It, typename Compare>
void parallel_sort(It first, It last, const Compare& comp) {
    std::sort(first, last, comp);
}

class spin_mutex {
public:
    class scoped_lock {
    public:
        scoped_lock() = default;
        explicit scoped_lock(spin_mutex&) {}
        void acquire(spin_mutex&) {}
        void release() {}
    };
    void lock() {}
    void unlock() {}
};

template <typename T>
class concurrent_vector : public std::vector<T> {
public:
    using std::vector<T>::vector;
};

template <typename T, typename Compare = std::less<T>>
class concurrent_priority_queue {
public:
    void push(const T& v) { q_.push(v); }
    bool try_pop(T& v) {
        if (q_.empty()) return false;
        v = q_.top();
        q_.pop();
        return true;
    }
    bool empty() const { return q_.empty(); }
    std::size_t size() const { return q_.size(); }
    void clear() { q_ = {}; }

private:
    std::priority_queue<T, std::vector<T>, Compare> q_;
};

// The classic task API (bvh.cpp's BVH build; the remesher does not build a BVH: it extracts without
// smoothing). Serial: spawn runs the task at once; a task returning itself (recycled) runs again.
// Tasks are not freed (the path is unused at run time).
struct task_alloc_proxy {};

class task {
public:
    virtual ~task() = default;
    virtual task* execute() = 0;
    static task_alloc_proxy allocate_root() { return {}; }
    task_alloc_proxy allocate_child() { return {}; }
    task_alloc_proxy allocate_continuation() { return {}; }
    void set_ref_count(int) {}
    void recycle_as_child_of(task&) {}
    static void spawn(task& t) { run(t); }
    static void spawn_root_and_wait(task& t) { run(t); }

private:
    static void run(task& t) {
        task* cur = &t;
        while (cur->execute() == cur) {
        }
    }
};

class empty_task : public task {
public:
    task* execute() override { return nullptr; }
};

class task_scheduler_init {
public:
    static const int automatic = -1;
    explicit task_scheduler_init(int = automatic) {}
    static int default_num_threads() { return 1; }
};

} // namespace tbb

inline void* operator new(std::size_t n, const tbb::task_alloc_proxy&) { return ::operator new(n); }
inline void operator delete(void* p, const tbb::task_alloc_proxy&) { ::operator delete(p); }
