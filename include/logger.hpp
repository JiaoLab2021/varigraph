#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <string>

inline bool DEBUG_ENABLED = false;

inline std::ostream& diagnostic_stream() { return std::cerr; }
inline std::ostream& output_stream() { return std::cout; }

namespace logger_detail {

inline const auto start_wall = std::chrono::steady_clock::now();
inline const std::clock_t start_cpu = std::clock();

inline std::string format_function(const char* function) {
    std::string text(function);
    text.resize(10, '.');
    return text;
}

inline std::string format_number(double number, std::size_t width) {
    std::string text = std::to_string(number);
    const auto dot = text.find('.');
    width = std::max(width, dot == std::string::npos ? text.size() : dot);
    text.resize(width, '.');
    return text;
}

inline std::ostream& write_prefix(char level, const char* function) {
    const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - start_wall).count();
    const double cpu = static_cast<double>(std::clock() - start_cpu) / CLOCKS_PER_SEC;
    const std::string prefix = std::string("[") + level + "::" + format_function(function) + "::"
        + format_number(wall, 6) + "*" + format_number(wall > 0.0 ? cpu / wall : 0.0, 4) + "] ";
    return diagnostic_stream().write(prefix.data(), static_cast<std::streamsize>(prefix.size()));
}

}

inline void set_debug(bool enabled) { DEBUG_ENABLED = enabled; }

#define log_stream() logger_detail::write_prefix('M', __func__)
#define stage_stream() (diagnostic_stream().put('\n'), logger_detail::write_prefix('M', __func__))
#define warning_stream() logger_detail::write_prefix('W', __func__)
#define error_stream() logger_detail::write_prefix('E', __func__)
#define debug_block() if (!DEBUG_ENABLED) {} else
#define debug_stream() debug_block() logger_detail::write_prefix('D', __func__)
