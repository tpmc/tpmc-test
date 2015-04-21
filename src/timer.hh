#ifndef TPMC_TEST_TIMER_HH
#define TPMC_TEST_TIMER_HH

#include <chrono>

namespace tpmc_test {
  class Timer {
    typedef std::chrono::high_resolution_clock clock;
    typedef std::chrono::duration<double, std::chrono::seconds::period> duration;
  public:

    Timer();
    duration lap();
    duration total() const;
  private:
    clock::time_point init_;
    clock::time_point last_;
  };
}

#endif // TPMC_TEST_TIMER_HH
