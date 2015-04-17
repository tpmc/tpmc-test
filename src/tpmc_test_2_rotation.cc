#include <iostream>
#include <iomanip>
#include <numeric>
#include <boost/range/algorithm/transform.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <Eigen/Dense>
#include <tpmc/marchingcubes.hh>
#include <tpmc/fieldtraits.hh>
#include "grid.hh"
#include "geometry.hh"
#include "levelsets.hh"
#include "bisection.hh"

namespace tpmc
{
  template <class T, int dim>
  struct FieldTraits<Eigen::Matrix<T, dim, 1> >
  {
    typedef T field_type;
  };
}

template <int dim, class domain_type>
class OptimizationFunction
{
public:
  typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;

  OptimizationFunction(const tpmc_test::Grid<dim>& grid, int angleDegree)
      : grid_(grid)
  {
    field_type angle = angleDegree / 180.0 * M_PI;
    firstCenter_ << 0.5, 0.5, 0.375;
    firstNormal_ << std::cos(angle), std::sin(angle), 0;
    secondCenter_ << 0.5, 0.5, 0.625;
    secondNormal_ << -std::sin(angle), std::cos(angle), 0;
  }

  field_type operator()(field_type x) const
  {
    field_type ringRadius = 0.5 * (0.25 - x);

    auto firstTorus = tpmc_test::torusLevelSet(0.25, ringRadius, firstCenter_, firstNormal_);
    // second torus
    domain_type secondCenter, secondNormal;
    auto secondTorus = tpmc_test::torusLevelSet(0.25, ringRadius, secondCenter_, secondNormal_);
    // combined level sets
    auto tori = [firstTorus, secondTorus](const domain_type& x) {
      return std::min(firstTorus(x), secondTorus(x));
    };
    std::vector<std::vector<std::pair<unsigned int, unsigned int> > > localGroups(
        grid_.vertexCounts().prod());
    typedef std::pair<unsigned int, short> localGroupType;
    typedef std::map<localGroupType,int> rankType;
    typedef boost::associative_property_map<rankType> rankMap;
    typedef std::map<localGroupType,localGroupType> parentType;
    typedef boost::associative_property_map<parentType> parentMap;
    rankType rank;
    rankMap prank(rank);
    parentType parent;
    parentMap pparent(parent);
    boost::disjoint_sets<rankMap, parentMap> dsets(prank, pparent);
    // loop over all elements
    for (auto e : grid_.elements()) {
      unsigned int linear = tpmc_test::GridIndex<dim>(grid_.elementCounts(), e.index()).linear();
      // compute corner values
      std::vector<field_type> cornerValues;
      boost::transform(e.corners(), std::back_inserter(cornerValues), tori);
      // fetch key
      std::size_t key = mc33_.getKey(std::begin(cornerValues), std::end(cornerValues));
      // compute vertices
      std::vector<domain_type> vertices;
      mc33_.getVertices(std::begin(cornerValues), std::end(cornerValues), key,
                        std::back_inserter(vertices));
      std::vector<int> vertexGroups;
      mc33_.getVertexGroups(tpmc::makeGeometryType(dim, e.cornerCount()), key,
                            std::back_inserter(vertexGroups));
      for (unsigned int i = 0; i< e.cornerCount(); ++i) {
        localGroups[grid_.cornerIndex(e, i)].push_back(std::make_pair(linear , vertexGroups[i]));
        dsets.make_set(std::make_pair(linear, vertexGroups[i]));
      }
    }
    for (int i = 0; i<localGroups.size(); ++i) {
      for (int j = 1; j<localGroups[i].size(); ++j) {
        dsets.union_set(localGroups[i][0], localGroups[i][j]);
      }
    }
    unsigned int count = 0;
    for (auto& x: parent) {
      count += x.first == x.second;
    }
    return count;
  }

private:
  tpmc_test::Grid<dim> grid_;
  tpmc::MarchingCubes<field_type, dim, domain_type> mc33_;
  domain_type firstNormal_;
  domain_type firstCenter_;
  domain_type secondNormal_;
  domain_type secondCenter_;
};

int main()
{
  // seed random generator
  srand(time(0));

  // define general grid properties
  const int dim = 3;
  const unsigned int numberOfElements = 32;
  typedef tpmc_test::Grid<dim>::domain_type domain_type;
  typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;
  domain_type low = domain_type::Zero();
  domain_type high = domain_type::Ones();

  // construct Grid
  typedef tpmc_test::Grid<dim>::dimension_type dimension_type;
  dimension_type elements = dimension_type::Constant(numberOfElements);
  tpmc_test::Grid<dim> grid(elements, low, high);

  // construct marching cubes class
  field_type ie = 1.0/numberOfElements;

  for (int angleDegree = 0; angleDegree < 180; ++angleDegree) {
    OptimizationFunction<dim, domain_type> opti(grid, angleDegree);
    auto lambda = [opti](field_type x) { return opti(x)-2.5; };
    tpmc_test::Bisection<field_type, field_type> bisection(1e-5);
    field_type v = bisection.apply(lambda, 0.5*ie, 2.5*ie);
    std::cout << angleDegree << " " << std::setprecision(15) << v/ie << "\n";
  }
}
