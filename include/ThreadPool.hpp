#pragma once

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <functional>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

// A bounded worker pool; the owner submits work and drains it before destruction.
class ThreadPool {
public:
    explicit ThreadPool(std::size_t threads, std::size_t capacity = 0)
        : thread_count_(std::max<std::size_t>(threads, 1)),
          capacity_(capacity ? capacity : 2 * thread_count_) {
        workers_.reserve(thread_count_);
        try {
            for (std::size_t i = 0; i < thread_count_; ++i) {
                workers_.emplace_back([this] { work(); });
            }
        } catch (...) {
            stop();
            throw;
        }
    }

    ~ThreadPool() { stop(); }
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    template<class F, class... Args>
    auto submit(F&& function, Args&&... args) -> std::future<std::invoke_result_t<F, Args...>> {
        using Result = std::invoke_result_t<F, Args...>;
        auto task = std::make_shared<std::packaged_task<Result()>>(
            std::bind(std::forward<F>(function), std::forward<Args>(args)...));
        auto future = task->get_future();
        {
            std::unique_lock<std::mutex> lock(mutex_);
            space_.wait(lock, [this] { return stopped_ || tasks_.size() < capacity_; });
            if (stopped_) throw std::runtime_error("submit on stopped ThreadPool");
            tasks_.emplace([task] { (*task)(); });
        }
        ready_.notify_one();
        return future;
    }

    // Borrowed loop state stays alive until every submitted worker has finished.
    template<class F>
    void parallel_for(std::size_t count, F&& function) {
        if (count == 0) return;
        const std::size_t workers = std::min(count, thread_count_);
        std::atomic<std::size_t> next{0};
        std::vector<std::future<void>> futures;
        futures.reserve(workers);
        std::exception_ptr error;
        try {
            for (std::size_t i = 0; i < workers; ++i) {
                futures.emplace_back(submit([&next, count, &function] {
                    for (;;) {
                        const auto index = next.fetch_add(1, std::memory_order_relaxed);
                        if (index >= count) return;
                        std::invoke(function, index);
                    }
                }));
            }
        } catch (...) {
            error = std::current_exception();
        }
        for (auto& future : futures) {
            try {
                future.get();
            } catch (...) {
                if (!error) error = std::current_exception();
            }
        }
        if (error) std::rethrow_exception(error);
    }

    // Called by the owner; queued work is completed before the workers are joined.
    void stop() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopped_ = true;
        }
        ready_.notify_all();
        space_.notify_all();
        for (auto& worker : workers_) {
            if (worker.joinable()) worker.join();
        }
    }

    std::size_t thread_count() const { return thread_count_; }
    std::size_t queue_capacity() const { return capacity_; }

private:
    void work() {
        for (;;) {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                ready_.wait(lock, [this] { return stopped_ || !tasks_.empty(); });
                if (tasks_.empty()) return;
                task = std::move(tasks_.front());
                tasks_.pop();
            }
            space_.notify_one();
            task();  // packaged_task delivers exceptions through its future.
        }
    }

    const std::size_t thread_count_;
    const std::size_t capacity_;
    std::mutex mutex_;
    std::condition_variable ready_;
    std::condition_variable space_;
    std::queue<std::function<void()>> tasks_;
    bool stopped_ = false;
    std::vector<std::thread> workers_;
};
