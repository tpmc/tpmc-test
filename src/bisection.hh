#ifndef TPMC_TEST_BISECTION_HH
#define TPMC_TEST_BISECTION_HH

#include <functional>
#include <limits>

namespace tpmc_test
{
  template <class Domain, class Range>
  class Bisection
  {
  public:
    explicit Bisection(const Domain& tolerance,
                       unsigned int maxIterations = std::numeric_limits<unsigned int>::max())
        : tolerance_(tolerance)
        , maxIterations_(maxIterations)
    {
    }

    Domain apply(std::function<Range(Domain)> f, Domain low, Domain high) const
    {
      Range lowValue = f(low);
      Range highValue = f(high);
      Domain middle;
      for (unsigned int iteration = 0; iteration < maxIterations_; ++iteration) {
        middle = 0.5 * (high + low);
        Range middleValue = f(middle);
        if (std::abs(middleValue) < tolerance_  || std::abs(high-low) < tolerance_)
          break;
        if (middleValue * lowValue > 0) {
          low = middle;
          lowValue = middleValue;
        } else {
          high = middle;
          highValue = middleValue;
        }
      }
      return middle;
    }

  private:
    Domain tolerance_;
    unsigned int maxIterations_;
  };
}

#endif // TPMC_TEST_BISECTION_HH
