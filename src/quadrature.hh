#ifndef TPMC_TEST_QUADRATURE_HH
#define TPMC_TEST_QUADRATURE_HH

#include <Eigen/Dense>
#include <vector>
#include <map>

namespace tpmc_test
{
  template <class T, int dim>
  class QuadratureRules;

  template <class T, int dim>
  class QuadratureRule;

  template <class T, int dim>
  class QuadraturePoint;

  namespace
  {
    template <class T, int dim>
    struct QuadratureRuleFactory;

    template <class T>
    struct QuadratureRuleFactory<T, 1>
    {
      template <class C>
      static void fill(C& container)
      {
        typedef QuadraturePoint<T, 1> point_type;
        typedef typename point_type::domain_type domain_type;
        {
          QuadratureRule<T, 1> rule;
          rule.points_.emplace_back(domain_type::Constant(0.5), 1.0);
          container[1] = rule;
        }
        {
          QuadratureRule<T, 1> rule;
          rule.points_.emplace_back(domain_type::Constant(0.5 - 1.0 / (2.0 * std::sqrt(3.0))), 0.5);
          rule.points_.emplace_back(domain_type::Constant(0.5 + 1.0 / (2.0 * std::sqrt(3.0))), 0.5);
          container[2] = rule;
        }
      }
    };

    template <class T, int dim>
    struct QuadratureRuleFactory
    {
      template <class C>
      static void fill(C& container)
      {
        typedef QuadraturePoint<T, dim> point_type;
        typedef typename point_type::domain_type domain_type;
        std::vector<unsigned int> orders{ { 1, 2 } };
        for (unsigned int order : orders) {
          QuadratureRule<T, dim> rule;
          const QuadratureRule<T, dim - 1>& low_rule = QuadratureRules<T, dim - 1>::get(order);
          const QuadratureRule<T, 1>& one_rule = QuadratureRules<T, 1>::get(order);
          for (auto low_point : low_rule) {
            domain_type d;
            d.segment(0, dim - 1) = low_point.position();
            for (auto one_point : one_rule) {
              d(dim - 1) = one_point.position()(0);
              rule.points_.emplace_back(d, low_point.weight() * one_point.weight());
            }
          }
          container[order] = rule;
        }
      }
    };
  }

  template <class T, int dim>
  class QuadraturePoint
  {
  public:
    typedef Eigen::Matrix<T, dim, 1> domain_type;
    typedef typename domain_type::Scalar field_type;

    QuadraturePoint(const domain_type& position, field_type weight)
        : position_(position)
        , weight_(weight)
    {
    }

    const domain_type& position() const { return position_; }

    field_type weight() const { return weight_; }

  private:
    domain_type position_;
    field_type weight_;
  };

  template <class T, int dim>
  class QuadratureRule
  {
  public:
    typedef QuadraturePoint<T, dim> point_type;
    typedef std::vector<point_type> container_type;
    typedef typename container_type::const_iterator const_iterator;

    const_iterator begin() const { return points_.begin(); }
    const_iterator end() const { return points_.end(); }

  private:
    template <class U, int d>
    friend class QuadratureRuleFactory;

    std::vector<QuadraturePoint<T, dim> > points_;
  };

  template <class T, int dim>
  class QuadratureRules
  {
    typedef std::map<unsigned int, QuadratureRule<T, dim> > container_type;

  public:
    static const QuadratureRule<T, dim>& get(unsigned int order)
    {
      return instance().getImpl(order);
    }

  private:
    QuadratureRules() { QuadratureRuleFactory<T, dim>::fill(rules_); }

    const QuadratureRule<T, dim>& getImpl(unsigned int order) const
    {
      typename container_type::const_iterator it = rules_.lower_bound(order);
      if (it != rules_.end()) {
        return it->second;
      } else {
        throw std::exception();
      }
    }

    static const QuadratureRules& instance()
    {
      static QuadratureRules instance;
      return instance;
    }

    container_type rules_;
  };
}

#endif // TPMC_TEST_QUADRATURE_HH
