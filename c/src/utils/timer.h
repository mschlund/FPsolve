#pragma once

#include <chrono>

class Timer {
  typedef std::chrono::high_resolution_clock clock;

  public:
    typedef std::chrono::microseconds Microseconds;
    typedef std::chrono::milliseconds Milliseconds;

    Timer() : start_(), stop_() {}

    void Start() {
      start_ = clock::now();
      stop_ = start_;
    }

    void Stop() {
      stop_ = clock::now();
    }

    Microseconds GetMicroseconds() const {
      return std::chrono::duration_cast<Microseconds>(stop_ - start_);
    }

    Milliseconds GetMilliseconds() const {
      return std::chrono::duration_cast<Milliseconds>(stop_ - start_);
    }

  private:
    clock::time_point start_;
    clock::time_point stop_;
};
