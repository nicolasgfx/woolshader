// ByteWriter, ByteReader, whole-file read and atomic replace, file hashes and project paths (files.h).
#include "io/files.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>

#include "model/hash.h"

namespace knit {

void ByteWriter::bytes(const void* p, size_t n) { out_.append(static_cast<const char*>(p), n); }
void ByteWriter::u8(uint8_t v) { bytes(&v, sizeof v); }
void ByteWriter::i32(int32_t v) { bytes(&v, sizeof v); }
void ByteWriter::u32(uint32_t v) { bytes(&v, sizeof v); }
void ByteWriter::i64(int64_t v) { bytes(&v, sizeof v); }
void ByteWriter::u64(uint64_t v) { bytes(&v, sizeof v); }
void ByteWriter::f64(double v) { bytes(&v, sizeof v); }

void ByteWriter::vec2(const Vec2& v) {
    f64(v.x);
    f64(v.y);
}

void ByteWriter::vec3(const Vec3& v) {
    f64(v.x);
    f64(v.y);
    f64(v.z);
}

void ByteWriter::string32(const std::string& s) {
    u32(static_cast<uint32_t>(s.size()));
    bytes(s.data(), s.size());
}

void ByteWriter::string64(const std::string& s) {
    i64(static_cast<int64_t>(s.size()));
    bytes(s.data(), s.size());
}

void ByteWriter::doubles(const std::vector<double>& v) {
    i64(static_cast<int64_t>(v.size()));
    bytes(v.data(), v.size() * sizeof(double));
}

void ByteWriter::ints(const std::vector<int>& v) {
    i64(static_cast<int64_t>(v.size()));
    for (int x : v) i64(x);
}

void ByteReader::bytes(void* p, size_t n) {
    if (!ok_ || at_ + n > in_.size()) {
        ok_ = false;
        std::memset(p, 0, n);
        return;
    }
    std::memcpy(p, in_.data() + at_, n);
    at_ += n;
}

uint8_t ByteReader::u8() {
    uint8_t v = 0;
    bytes(&v, sizeof v);
    return v;
}

int32_t ByteReader::i32() {
    int32_t v = 0;
    bytes(&v, sizeof v);
    return v;
}

uint32_t ByteReader::u32() {
    uint32_t v = 0;
    bytes(&v, sizeof v);
    return v;
}

int64_t ByteReader::i64() {
    int64_t v = 0;
    bytes(&v, sizeof v);
    return v;
}

uint64_t ByteReader::u64() {
    uint64_t v = 0;
    bytes(&v, sizeof v);
    return v;
}

double ByteReader::f64() {
    double v = 0;
    bytes(&v, sizeof v);
    return v;
}

Vec2 ByteReader::vec2() {
    Vec2 v;
    v.x = f64();
    v.y = f64();
    return v;
}

Vec3 ByteReader::vec3() {
    Vec3 v;
    v.x = f64();
    v.y = f64();
    v.z = f64();
    return v;
}

std::string ByteReader::string32() {
    const uint32_t n = u32();
    if (!ok_ || n > in_.size() - at_) {
        ok_ = false;
        return {};
    }
    std::string s(n, '\0');
    bytes(s.data(), n);
    return s;
}

std::string ByteReader::string64() {
    std::string s(count(1), '\0');
    bytes(s.data(), s.size());
    return s;
}

std::vector<double> ByteReader::doubles() {
    std::vector<double> v(count(sizeof(double)));
    bytes(v.data(), v.size() * sizeof(double));
    return v;
}

std::vector<int> ByteReader::ints() {
    std::vector<int> v(count(sizeof(int64_t)));
    for (int& x : v) x = static_cast<int>(i64());
    return v;
}

size_t ByteReader::count(size_t bytesPerItem) {
    const int64_t n = i64();
    const size_t left = in_.size() - std::min(at_, in_.size());
    if (n < 0 || static_cast<uint64_t>(n) > left / std::max<size_t>(bytesPerItem, 1)) {
        ok_ = false;
        return 0;
    }
    return static_cast<size_t>(n);
}

bool fail(std::string* err, const std::string& message) {
    if (err) *err = message;
    return false;
}

bool readWholeFile(const std::string& path, std::string& out, std::string* err) {
    std::ifstream f(path, std::ios::binary);
    if (!f) return fail(err, "cannot open " + path);
    f.seekg(0, std::ios::end);
    const std::streamoff size = f.tellg();
    f.seekg(0);
    out.resize(static_cast<size_t>(size));
    if (size > 0) f.read(&out[0], size);
    if (!f) return fail(err, "cannot read " + path);
    return true;
}

bool replaceFile(const std::string& path, const std::string& bytes, std::string* err) {
    const std::string tmp = path + ".tmp";
    {
        std::ofstream f(tmp, std::ios::binary);
        f.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
        if (!f) return fail(err, "cannot write " + tmp);
    }
    std::error_code ec;
    std::filesystem::rename(tmp, path, ec);
    if (ec) return fail(err, "cannot rename " + tmp + ": " + ec.message());
    return true;
}

uint64_t fnv1aFile(const std::string& path, bool* ok) {
    std::ifstream f(path, std::ios::binary);
    if (ok) *ok = static_cast<bool>(f);
    uint64_t h = kFnvStandardBasis;
    char buf[1 << 16];
    while (f) {
        f.read(buf, sizeof buf);
        h = fnv1a(buf, static_cast<size_t>(f.gcount()), h);
    }
    return h;
}

std::string findProjectFile(const std::string& relativePath) {
    std::error_code ec;
    if (std::filesystem::exists(relativePath, ec)) return relativePath;
    const std::string tree = sourceTreeDir();
    if (tree.empty()) return "";
    const std::string inSourceTree = tree + "/" + relativePath;
    if (std::filesystem::exists(inSourceTree, ec)) return inSourceTree;
    return "";
}

std::string sourceTreeDir() {
#ifdef KNIT_SOURCE_DIR
    return KNIT_SOURCE_DIR;
#else
    return "";
#endif
}

std::string fileHashHex(const std::string& path) {
    bool ok = false;
    const uint64_t h = fnv1aFile(path, &ok);
    if (!ok) return "";
    char b[24];
    std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(h));
    return b;
}

} // namespace knit
