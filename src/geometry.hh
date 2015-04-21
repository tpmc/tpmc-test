#ifndef TPMC_TEST_GEOMETRY_HH
#define TPMC_TEST_GEOMETRY_HH

#include <Eigen/Dense>
#include <tpmc/fieldtraits.hh>
#include "exceptions.hh"

namespace tpmc_test
{
  double quadrilateralArea(const Eigen::Matrix<double, 3, 1>& v0,
                           const Eigen::Matrix<double, 3, 1>& v1,
                           const Eigen::Matrix<double, 3, 1>& v2,
                           const Eigen::Matrix<double, 3, 1>& v3);
  double triangleArea(const Eigen::Matrix<double, 3, 1>& v0, const Eigen::Matrix<double, 3, 1>& v1,
                      const Eigen::Matrix<double, 3, 1>& v2);

  // compute the 2 dimension measure of a face
  template <class field_type, class InputIterator>
  field_type area(InputIterator begin, InputIterator end)
  {
    std::size_t corners = std::distance(begin, end);
    switch (corners) {
    case 3:
      return triangleArea(begin[0], begin[1], begin[2]);
    case 4:
      return quadrilateralArea(begin[0], begin[1], begin[2], begin[3]);
    }
    throw NotImplementedException("area not implemented for given face");
  }
}

#endif // TPMC_TEST_GEOMETRY_HH
