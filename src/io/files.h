#pragma once
// What every file format stands on: a record as bytes in memory, a whole file in or out, a file's
// hash and a project file's path. A format turns its record into bytes with ByteWriter, back with
// ByteReader, and only these functions touch the disk.
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "model/vec.h"

namespace knit {

// A little-endian record being built in memory. Widths are named because formats differ: .kmesh
// writes 32-bit indices and string lengths, the knit formats 64-bit integers everywhere.
class ByteWriter {
public:
    void bytes(const void* p, size_t n);
    void u8(uint8_t v);
    void i32(int32_t v);
    void u32(uint32_t v);
    void i64(int64_t v);
    void u64(uint64_t v);
    void f64(double v);                          // bit for bit
    void vec2(const Vec2& v);
    void vec3(const Vec3& v);
    void string32(const std::string& s);         // u32 length (.kmesh)
    void string64(const std::string& s);         // i64 length (the knit formats)
    void doubles(const std::vector<double>& v);  // i64 count, then the values
    void ints(const std::vector<int>& v);        // i64 count, each int as an i64
    const std::string& written() const { return out_; }

private:
    std::string out_;
};

// Reads what ByteWriter wrote from `in`, which must outlive it. Past the end every read gives zero
// and ok() stays false, so a reader asks once, after the last field, instead of after each one.
class ByteReader {
public:
    ByteReader(const std::string& in, size_t start) : in_(in), at_(start) {}
    bool ok() const { return ok_; }
    size_t at() const { return at_; }
    void refuse() { ok_ = false; }
    void bytes(void* p, size_t n);
    uint8_t u8();
    int32_t i32();
    uint32_t u32();
    int64_t i64();
    uint64_t u64();
    double f64();
    Vec2 vec2();
    Vec3 vec3();
    std::string string32();
    std::string string64();
    std::vector<double> doubles();
    std::vector<int> ints();
    // An i64 count of items the remaining bytes can hold; otherwise 0, and the reader refuses. A
    // corrupt file therefore never allocates beyond its own size.
    size_t count(size_t bytesPerItem);

private:
    const std::string& in_;
    size_t at_;
    bool ok_ = true;
};

// Sets *err when it is given; returns false, so a refusal is one line: return fail(err, "...").
[[nodiscard]] bool fail(std::string* err, const std::string& message);

[[nodiscard]] bool readWholeFile(const std::string& path, std::string& out, std::string* err);
// Written beside the target and renamed over it, so a reader never sees a partial file.
[[nodiscard]] bool replaceFile(const std::string& path, const std::string& bytes, std::string* err);

// FNV-1a 64 of a file's bytes, standard basis (asset content hashes).
uint64_t fnv1aFile(const std::string& path, bool* ok = nullptr);

// An existing file for a repo-relative path, from the working directory and then the source tree;
// "" if there is none. Scenes, presets and tools are all found this way.
std::string findProjectFile(const std::string& relativePath);
// The source tree configured at build time (KNIT_SOURCE_DIR), "" in a build without it. The define
// reaches this function alone.
std::string sourceTreeDir();
// A file's FNV-1a digest as 16 hex digits (fnv1aFile); "" when it cannot be read.
std::string fileHashHex(const std::string& path);

} // namespace knit
