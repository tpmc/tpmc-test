#include "geometry.hh"
#include "quadrature.hh"

namespace tpmc_test
{
  double triangleArea(const Eigen::Matrix<double, 3, 1>& v0, const Eigen::Matrix<double, 3, 1>& v1,
                      const Eigen::Matrix<double, 3, 1>& v2)
  {
    // return area using crossproduct formula
    return (v2 - v0).cross(v1 - v0).norm() * 0.5;
  }

  double quadrilateralArea(const Eigen::Matrix<double, 3, 1>& v0,
                           const Eigen::Matrix<double, 3, 1>& v1,
                           const Eigen::Matrix<double, 3, 1>& v2,
                           const Eigen::Matrix<double, 3, 1>& v3)
  {
    typedef double field_type;
    typedef Eigen::Matrix<field_type, 3, 1> domain_type;
    domain_type a = v0 - v1 - v2 + v3;
    domain_type b = v1 - v0;
    domain_type c = v2 - v0;
    domain_type d = v0;
    field_type ca = a.squaredNorm();
    field_type cb = a.dot(b);
    field_type cc = b.squaredNorm();
    field_type cd = a.dot(c);
    field_type ce = b.dot(c);
    field_type cf = c.squaredNorm();
    auto ol = [ca, cb, cc](const Eigen::Matrix<field_type, 2, 1>& x) {
      return ca * x(1) * x(1) + 2 * cb * x(1) + cc;
    };
    auto ur = [ca, cd, cf](const Eigen::Matrix<field_type, 2, 1>& x) {
      return ca * x(0) * x(0) + 2 * cd * x(0) + cf;
    };
    auto diag = [ca, cd, cb, ce](const Eigen::Matrix<field_type, 2, 1>& x) {
      return ca * x(0) * x(1) + cd * x(1) + cb * x(0) + ce;
    };
    auto det = [ol, ur, diag](const Eigen::Matrix<field_type, 2, 1>& x) {
      field_type d = diag(x);
      return std::sqrt(ol(x) * ur(x) - d * d);
    };
    const tpmc_test::QuadratureRule<field_type, 2>& rule
        = tpmc_test::QuadratureRules<field_type, 2>::get(2);
    field_type sum = 0;
    for (auto qp : rule) {
      sum += det(qp.position()) * qp.weight();
    }
    return sum;
  }
}
