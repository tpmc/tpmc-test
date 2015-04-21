#include "timer.hh"

namespace tpmc_test
{
  Timer::Timer()
      : init_(clock::now())
      , last_(init_)
  {
  }
  Timer::duration Timer::lap()
  {
    clock::time_point now = clock::now();
    duration result = std::chrono::duration_cast<duration>(now - last_);
    last_ = now;
    return result;
  }

  Timer::duration Timer::total() const
  {
    return std::chrono::duration_cast<duration>(clock::now() - init_);
  }
}
