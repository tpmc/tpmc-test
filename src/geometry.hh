#ifndef TPMC_TEST_GEOMETRY_HH
#define TPMC_TEST_GEOMETRY_HH

#include <exception>
#include <Eigen/Dense>
#include <tpmc/fieldtraits.hh>
#include "quadrature.hh"

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
      /*
      // for now, just used the mean of both triangulations
      field_type first = area<field_type>(begin, begin + 3) + area<field_type>(begin + 1, end);
      std::vector<domain_type> reordered{ { begin[1], begin[0], begin[3], begin[2] } };
      field_type second = area<field_type>(reordered.begin(), reordered.begin() + 3)
                          + area<field_type>(reordered.begin() + 1, reordered.end());
      return 0.5*(first+second);
      */
      domain_type a = begin[0]-begin[1]-begin[2]+begin[3];
      domain_type b = begin[1]-begin[0];
      domain_type c = begin[2]-begin[0];
      domain_type d = begin[0];
      field_type ca = a.squaredNorm();
      field_type cb = a.dot(b);
      field_type cc = b.squaredNorm();
      field_type cd = a.dot(c);
      field_type ce = b.dot(c);
      field_type cf = c.squaredNorm();
      auto ol = [ca,cb,cc](const Eigen::Matrix<field_type,2,1>& x) { return ca*x(1)*x(1)+2*cb*x(1)+cc; };
      auto ur = [ca,cd,cf](const Eigen::Matrix<field_type,2,1>& x) { return ca*x(0)*x(0)+2*cd*x(0)+cf; };
      auto diag = [ca,cd,cb,ce](const Eigen::Matrix<field_type,2,1>& x) { return ca*x(0)*x(1)+cd*x(1)+cb*x(0)+ce; };
      auto det = [ol,ur,diag](const Eigen::Matrix<field_type,2,1>& x) { 
        field_type d = diag(x);
        return std::sqrt(ol(x)*ur(x)-d*d);
      };
      const tpmc_test::QuadratureRule<double,2>& rule = tpmc_test::QuadratureRules<double,2>::get(2);
      typedef tpmc_test::QuadratureRule<double,2>::point_type::domain_type domain_type;
      double sum = 0;
      for (auto qp: rule) {
        sum += det(qp.position())*qp.weight();
      }
      return sum;
    }
    throw std::exception();
  }
}

#endif // TPMC_TEST_GEOMETRY_HH
