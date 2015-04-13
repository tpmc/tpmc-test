#ifndef TPMC_TEST_GEOMETRY_HH
#define TPMC_TEST_GEOMETRY_HH

#include <exception>
#include <tpmc/fieldtraits.hh>

namespace tpmc_test
{
  // compute the d-1 dimension measure of a face
  template <class field_type, class InputIterator>
  field_type area(InputIterator begin, InputIterator end)
  {
    typedef typename std::iterator_traits<InputIterator>::value_type domain_type;
    std::size_t corners = std::distance(begin, end);
    switch (corners) {
    case 3:
      // use the cross product formula
      return (begin[2] - begin[0]).cross(begin[1] - begin[0]).norm() * 0.5;
    case 4:
      // for now, just used the mean of both triangulations
      field_type first = area<field_type>(begin, begin + 3) + area<field_type>(begin + 1, end);
      std::vector<domain_type> reordered{ { begin[1], begin[0], begin[3], begin[2] } };
      field_type second = area<field_type>(reordered.begin(), reordered.begin() + 3)
                          + area<field_type>(reordered.begin() + 1, reordered.end());
      return 0.5*(first+second);
    }
    throw std::exception();
  }
}

#endif // TPMC_TEST_GEOMETRY_HH
