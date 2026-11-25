#pragma once
// A serial stand-in for pss::parallel_stable_sort (wjakob/pss, which needs the TBB task API):
// std::stable_sort, the same order.
#include <algorithm>

namespace pss {

template <class RandomAccessIterator>
void parallel_stable_sort(RandomAccessIterator first, RandomAccessIterator last) {
    std::stable_sort(first, last);
}
template <class RandomAccessIterator, class Compare>
void parallel_stable_sort(RandomAccessIterator first, RandomAccessIterator last, Compare comp) {
    std::stable_sort(first, last, comp);
}

} // namespace pss
