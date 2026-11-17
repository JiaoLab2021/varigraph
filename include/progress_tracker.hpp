#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <string>

#include "logger.hpp"

class ProgressTracker {
public:
    explicit ProgressTracker(std::size_t total)
        : ProgressTracker(total, 0) {}

    static ProgressTracker Every(std::size_t step) {
        return ProgressTracker(0, step ? step : 1);
    }

    void hit(std::size_t count = 1) {
        if (count == 0 || (step_ == 0 && total_ == 0)) return;
        const auto done = done_.fetch_add(count, std::memory_order_relaxed) + count;
        if (step_ != 0) {
            if (done < next_count_.load(std::memory_order_relaxed)) return;
            std::lock_guard<std::mutex> lock(print_mutex_);
            const auto current = done_.load(std::memory_order_relaxed);
            if (finished_ || current < next_count_.load(std::memory_order_relaxed)) return;
            log_stream() << "  # " << current << " done\n";
            next_count_.store((current / step_ + 1) * step_, std::memory_order_relaxed);
            return;
        }
        if (done > total_ || percentage(done) < next_percent_.load(std::memory_order_relaxed)) return;
        std::lock_guard<std::mutex> lock(print_mutex_);
        const auto current = done_.load(std::memory_order_relaxed);
        if (current > total_) return;
        const auto percent = percentage(current);
        if (percent < next_percent_.load(std::memory_order_relaxed)) return;
        log_stream() << "  # " << std::setw(3) << percent << "% (" << std::setw(width_) << current << '/' << total_ << ") done\n";
        next_percent_.store((percent / 10 + 1) * 10, std::memory_order_relaxed);
    }

    void finish() {
        if (step_ == 0) return;
        std::lock_guard<std::mutex> lock(print_mutex_);
        if (finished_) return;
        log_stream() << "  # all " << done_.load(std::memory_order_relaxed) << " done\n";
        finished_ = true;
    }

private:
    ProgressTracker(std::size_t total, std::size_t step)
        : total_(total), step_(step), width_(number_width(total)), next_count_(step) {}

    static int number_width(std::size_t number) {
        std::ostringstream text;
        text.imbue(diagnostic_stream().getloc());
        text << number;
        return static_cast<int>(text.str().size());
    }

    unsigned percentage(std::size_t done) const {
        return done == total_ ? 100 : std::min(99u, static_cast<unsigned>(100.0 * done / total_));
    }

    const std::size_t total_;
    const std::size_t step_;
    const int width_;
    std::atomic<std::size_t> done_{0};
    std::atomic<unsigned> next_percent_{10};
    std::atomic<std::size_t> next_count_;
    std::mutex print_mutex_;
    bool finished_ = false;
};
