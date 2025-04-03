#pragma once
#include <condition_variable>
#include <mutex>

class semaphore {
  std::mutex mutex_;
  std::condition_variable condition_;
  unsigned long count_ = 0; // Initialized as locked.

public:
  semaphore(unsigned long count = 0) : count_(count) {}
  void release(unsigned long count = 1) {
    std::lock_guard<decltype(mutex_)> lock(mutex_);
    for (unsigned long i = 0; i < count; ++i) {
      ++count_;
      condition_.notify_one();
    }
  }

  void acquire() {
    std::unique_lock<decltype(mutex_)> lock(mutex_);
    while (!count_) // Handle spurious wake-ups.
      condition_.wait(lock);
    --count_;
  }

  bool try_acquire() {
    std::lock_guard<decltype(mutex_)> lock(mutex_);
    if (count_) {
      --count_;
      return true;
    }
    return false;
  }
};