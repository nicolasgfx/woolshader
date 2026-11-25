#pragma once
// Stand-in for nanogui's serializer (RIM's save / load of its state; never called here).
#include <string>
namespace nanogui {
class Serializer {
public:
    template <typename T> void set(const std::string&, const T&) {}
    template <typename T> bool get(const std::string&, T&) { return false; }
    void push(const std::string&) {}
    void pop() {}
};
} // namespace nanogui
