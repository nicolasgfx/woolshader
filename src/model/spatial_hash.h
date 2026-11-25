#pragma once
// Uniform 3D grid hash of axis-aligned boxes: candidate enumeration only.
//
// Every item is registered in all cells its box overlaps. query(lo, hi) visits each item whose box
// overlaps [lo, hi] in at least one shared cell, exactly once; the caller still runs its own exact
// test, so the hash never decides a result: a broad phase using it returns the same pairs as an
// all-pairs loop with the same box filter (tested: hash_matches_bruteforce).
//
// Each item is visited once per query without any mutable state, so const queries may run
// concurrently: an item is reported only in its FIRST shared cell in query order (i, then j, then k),
// which is the cell (max(i0_item, i0_query), max(j0_item, j0_query), max(k0_item, k0_query)). So in
// a cell of the item's box, the axis i passes iff i is the item's lowest or the query's lowest i.
// Each bucket entry stores three flags "this is the item's lowest cell in i / j / k", each query cell
// has the same three flags for the query, and the item is visited iff every axis has one of the two
// (flags | cell == 7; tested: hash_canonical_cell_order).
//
// Supported range: |coordinate / cell| < 2^20 on every axis (21-bit cell keys; outside it distant
// cells alias), fewer than 2^29 items. A query's cell range is clipped to the range of occupied
// cells, and finalize() may replace the map by a dense array over that range, for speed. Neither
// changes which cells are visited, their order, or the order within a bucket.
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

#include "model/vec.h"

namespace knit {

class SpatialHash3 {
public:
    explicit SpatialHash3(double cell = 1.0) : cell_(cell > 0 ? cell : 1.0) {}

    double cell() const { return cell_; }

    // Items get consecutive ids 0, 1, ... in insertion order. Not after finalize().
    int insert(const Vec3& lo, const Vec3& hi);

    // Optional, after the last insert: a dense lookup over the occupied cell range if it has at most
    // kDenseCells cells (else the map stays). Same visits.
    void finalize();

    // Visits every item whose box shares a cell with [lo, hi] and whose box overlaps [lo, hi]
    // (componentwise, closed), once, in its first shared cell. Thread-safe for concurrent const use.
    void query(const Vec3& lo, const Vec3& hi, const std::function<void(int)>& visit) const;

    const Vec3& lo(int id) const { return lo_[id]; }
    const Vec3& hi(int id) const { return hi_[id]; }

private:
    using Bucket = std::vector<std::uint32_t>;  // item id in the low 29 bits, its flags above

    std::int64_t cellOf(double v) const;
    size_t denseIndex(std::int64_t i, std::int64_t j, std::int64_t k) const;
    const Bucket* bucketAt(std::int64_t i, std::int64_t j, std::int64_t k) const;  // nullptr: empty

    double cell_;
    std::unordered_map<std::int64_t, Bucket> cells_;
    std::vector<Vec3> lo_, hi_;
    // Occupied cell range (empty while lo > hi) and the optional dense lookup over it.
    std::int64_t occLo_[3] = {INT64_MAX, INT64_MAX, INT64_MAX}, occHi_[3] = {INT64_MIN, INT64_MIN, INT64_MIN};
    std::vector<std::int32_t> dense_;  // cell -> bucket index, -1 empty
    std::vector<Bucket> denseBuckets_;
};

} // namespace knit
