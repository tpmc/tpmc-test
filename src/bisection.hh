#ifndef TPMC_TEST_BISECTION_HH
#define TPMC_TEST_BISECTION_HH

#include <functional>
#include <limits>
#include <iostream>
#include "exceptions.hh"

namespace tpmc_test
{
  template <class Domain, class Range>
  class Bisection
  {
  public:
    explicit Bisection(const Domain& domainTolerance, const Range& rangeTolerance,
                       unsigned int maxIterations = std::numeric_limits<unsigned int>::max())
        : domainTolerance_(domainTolerance)
        , rangeTolerance_(rangeTolerance)
        , maxIterations_(maxIterations)
    {
    }

    Domain apply(std::function<Range(Domain)> f, Domain low, Domain high) const
    {
      Range lowValue = f(low);
      Range highValue = f(high);
      if (lowValue * highValue >= 0.0) {
        throw IllegalArgumentException(
            "precondition not met: values of low and high have the same sign");
      }
      Domain middle;
      for (unsigned int iteration = 0; iteration < maxIterations_; ++iteration) {
        middle = 0.5 * (high + low);
        Range middleValue = f(middle);
        if (std::abs(middleValue) < rangeTolerance_ || std::abs(high - low) < domainTolerance_)
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
    Domain domainTolerance_;
    Range rangeTolerance_;
    unsigned int maxIterations_;
  };
}

#endif // TPMC_TEST_BISECTION_HH
