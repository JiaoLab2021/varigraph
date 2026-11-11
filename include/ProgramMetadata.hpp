#pragma once

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>

namespace program {

inline constexpr const char* name = "varigraph";
inline constexpr const char* description = "Genotype diploid and polyploid samples using a pangenome graph.";
inline constexpr const char* version = "v1.1-r1";
inline constexpr const char* build_date = "2026-11-11";
inline constexpr const char* author = "Zezhen Du";
inline constexpr const char* email = "dzz0539@gmail.com or dzz0539@163.com";

inline constexpr const char* version_note = "Unifies CPU/GPU commands, help and logging, and improves task scheduling and error handling.";

inline std::string cmdline(int argc, char** argv) {
    std::ostringstream oss;
    std::copy(argv, argv + argc, std::ostream_iterator<const char*>(oss, " "));
    std::string s = oss.str();
    if (!s.empty() && s.back() == ' ')
        s.pop_back();
    return s;
}

}
