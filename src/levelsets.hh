#ifndef TPMC_TEST_LEVELSETS_HH
#define TPMC_TEST_LEVELSETS_HH

#include <Eigen/Dense>
#include <algorithm>
#include <functional>
#include <map>
#include <exception>
#include <boost/property_tree/ptree.hpp>
#include <tpmc/fieldtraits.hh>

namespace tpmc_test
{
  namespace
  {
    template <class domain_type, class T>
    domain_type parseVector(const std::vector<T>& vec)
    {
      domain_type result;
      for (int i = 0; i < vec.size(); ++i)
        result[i] = vec[i];
      return result;
    }

    template <typename T>
    std::vector<T> as_vector(boost::property_tree::ptree const& pt,
                             boost::property_tree::ptree::key_type const& key)
    {
      std::vector<T> r;
      for (auto& item : pt.get_child(key))
        r.push_back(item.second.get_value<T>());
      return r;
    }

    template <typename domain_type>
    domain_type as_domain_type(boost::property_tree::ptree const& pt,
                               boost::property_tree::ptree::key_type const& key)
    {
      typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;
      std::vector<field_type> v = as_vector<field_type>(pt, key);
      return parseVector<domain_type>(v);
    }

    template <class domain_type>
    struct RotationMatrix;

    template <class T>
    struct RotationMatrix<Eigen::Matrix<T, 3, 1> >
    {
      typedef Eigen::Matrix<T, 3, 1> domain_type;
      typedef Eigen::Matrix<T, 3, 3> matrix_type;

      static matrix_type construct(const domain_type& angles)
      {
        T cx = std::cos(angles(0)), cy = std::cos(angles(1)), cz = std::cos(angles(2));
        T sx = std::sin(angles(0)), sy = std::sin(angles(1)), sz = std::sin(angles(2));
        matrix_type matrix;
        matrix << cy* cz, cz* sx* sy - cx* sz, cx* cz* sy + sx* sz, cy* sz, cx* cz + sx* sy* sz,
            -cz* sx + cx* sy* sz, -sy, cy* sx, cx* cy;
        return matrix;
      }
    };
  }
  // tangle cube
  template <class domain_type>
  std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
  tangleCubeLevelSet(double psquare, double pconst, const domain_type& center)
  {
    typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;
    return [center, psquare, pconst](const domain_type& p) {
      field_type x = p(0) - center(0);
      field_type xx = x * x;
      field_type y = p(1) - center(1);
      field_type yy = y * y;
      field_type z = p(2) - center(2);
      field_type zz = z * z;
      return xx * xx + yy * yy + zz * zz - psquare * (xx + yy + zz) + pconst;
    };
  }
  // tangle cube
  template <class domain_type>
  std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
  tangleCubeLevelSet(const boost::property_tree::ptree& tree)
  {
    double psquare = tree.get<double>("squareCoefficient", 5);
    double pconst = tree.get<double>("constCoefficient", 11.8);
    domain_type center = as_domain_type<domain_type>(tree, "center");
    return tangleCubeLevelSet(psquare, pconst, center);
  }

  // sphere
  template <class domain_type>
  std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
  sphereLevelSet(double radius, const domain_type& center)
  {
    return [center, radius](const domain_type& p) { return (p - center).norm() - radius; };
  }

  // sphere
  template <class domain_type>
  std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
  sphereLevelSet(const boost::property_tree::ptree& tree)
  {
    double radius = tree.get<double>("radius", 1);
    domain_type center = as_domain_type<domain_type>(tree, "center");
    return sphereLevelSet(radius, center);
  }

  // torus
  template <class domain_type>
  std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
  torusLevelSet(double torusRadius, double ringRadius, const domain_type& center,
                const domain_type& normal)
  {
    typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;
    return [torusRadius, ringRadius, center, normal](const domain_type& p) {
      field_type pdotn = normal.dot(p - center);
      domain_type projection = p - center - pdotn * normal;
      field_type first = projection.norm() - torusRadius;
      return std::sqrt(first * first + pdotn * pdotn) - ringRadius;
    };
  }

  // torus
  template <class domain_type>
  std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
  torusLevelSet(const boost::property_tree::ptree& tree)
  {
    double torusRadius = tree.get<double>("torusRadius");
    double ringRadius = tree.get<double>("ringRadius");
    domain_type center = as_domain_type<domain_type>(tree, "center");
    domain_type normal = as_domain_type<domain_type>(tree, "normal");
    return torusLevelSet(torusRadius, ringRadius, center, normal);
  }

  template <class domain_type>
  std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
  constructLevelSet(const boost::property_tree::ptree& tree)
  {
    std::string type = tree.get<std::string>("type");
    if (type == "sphere")
      return sphereLevelSet<domain_type>(tree);
    if (type == "torus")
      return torusLevelSet<domain_type>(tree);
    if (type == "tangleCube")
      return tangleCubeLevelSet<domain_type>(tree);
    throw std::exception();
  }

  template <class domain_type, class OutIterator>
  void constructLevelSets(const boost::property_tree::ptree& tree, OutIterator out)
  {
    for (auto& t : tree)
      *out++ = constructLevelSet<domain_type>(t.second);
  }

  template <class domain_type, class C>
  std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
  constructCombinedLevelSet(C levelSets, const boost::property_tree::ptree& tree)
  {
    std::vector<int> indices = as_vector<int>(tree, "indices");
    typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;
    typedef std::function<field_type(field_type, field_type)> Reduction;
    // note: use lambda here due to failing overload resolution in std::min
    std::map<std::string, Reduction> reductions{
      { "min", [](field_type x, field_type y) { return std::min(x, y); } },
      { "max", [](field_type x, field_type y) { return std::max(x, y); } },
      { "multiplies", [](field_type x, field_type y) { return x * y; } }
    };
    Reduction reduction = reductions[tree.get<std::string>("reduction")];
    field_type init = tree.get<field_type>("init");
    domain_type translation = as_domain_type<domain_type>(tree, "translation");
    typedef RotationMatrix<domain_type> RM;
    typedef typename RM::matrix_type matrix_type;
    domain_type rot = as_domain_type<domain_type>(tree, "rotation");
    matrix_type rotation = RM::construct(rot);
    return [levelSets, reduction, init, translation, rotation](const domain_type& x) {
      domain_type tx = rotation * x - translation;
      std::vector<field_type> values;
      for (auto& f : levelSets)
        values.push_back(f(tx));
      return std::accumulate(values.begin(), values.end(), init, reduction);
    };
  }
}

#endif // TPMC_TEST_LEVELSETS_HH
