#pragma once
// Stand-in for a small header that the MPL-2.0 SLIM files include; this one is our own (MIT). A serial loop: the results do
// not depend on it, and the import runs its charts in parallel one level up.
#include <cstdint>

namespace slim_support {

class IndexRange {
public:
    explicit IndexRange(int64_t size) : size_(size) {}

    class Iterator {
    public:
        explicit Iterator(int64_t at) : at_(at) {}
        int64_t operator*() const { return at_; }
        Iterator& operator++() {
            ++at_;
            return *this;
        }
        bool operator!=(const Iterator& other) const { return at_ != other.at_; }

    private:
        int64_t at_;
    };

    Iterator begin() const { return Iterator(0); }
    Iterator end() const { return Iterator(size_); }

private:
    int64_t size_;
};

namespace threading {

template <typename Function>
void parallel_for(IndexRange range, int64_t, const Function& function) {
    function(range);
}

} // namespace threading
} // namespace slim_support
