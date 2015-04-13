#ifndef TPMC_TEST_TIMER_HH
#define TPMC_TEST_TIMER_HH

#include <chrono>

namespace tpmc_test {
  class Timer {
    typedef std::chrono::high_resolution_clock clock;
    typedef std::chrono::duration<double, std::chrono::seconds::period> duration;
  public:
    Timer() : init_(clock::now()), last_(init_) {
    }

    duration lap() {
      clock::time_point now = clock::now();
      duration result = std::chrono::duration_cast<duration>(now-last_);
      last_ = now;
      return result;
    }

    duration total() const {
      return std::chrono::duration_cast<duration>(clock::now()-init_);
    }
  private:
    clock::time_point init_;
    clock::time_point last_;
  };
}

#endif // TPMC_TEST_TIMER_HH
