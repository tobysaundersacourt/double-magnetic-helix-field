#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>

namespace utils {

inline void clear(const std::string& filename) {
    std::ofstream file(filename, std::ios::trunc);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for clearing.\n";
    }
}

template <typename T>
void append(const std::string& filename, const T& data) {
    if (data.empty()) return;
    std::ofstream file(filename, std::ios::app);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for appending.\n";
        return;
    }
    // Preserve enough precision that post-processing does not acquire artificial
    // staircase/ripple structure from truncated coordinates or field values.
    file << std::setprecision(17);
    auto it = data.begin();
    const auto end = data.end();
    while (it != end) {
        file << *it;
        auto next_it = it;
        ++next_it;
        file << (next_it != end ? "," : "\n");
        it = next_it;
    }
}

template<typename T>
T acos(T value, T pi) {
    if (value <= -1.0) return pi;
    if (value >= 1.0) return 0;
    return std::acos(value);
}

} // namespace utils
#endif
