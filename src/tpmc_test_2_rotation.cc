#include <iostream>
#include <iomanip>
#include <numeric>
#include <boost/range/algorithm/transform.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <Eigen/Dense>
#include <tpmc/marchingcubes.hh>
#include <tpmc/fieldtraits.hh>
#include "grid.hh"
#include "levelsets.hh"
#include "bisection.hh"
#include "io.hh"

namespace tpmc
{
  template <class T, int dim>
  struct FieldTraits<Eigen::Matrix<T, dim, 1> >
  {
    typedef T field_type;
  };
}

// simple functor which computes the number of connected components in a
// tpmc reconstruction of two entangled tori
template <int dim, class domain_type>
class ConnectedComponentsFunctor
{
public:
  typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;

  ConnectedComponentsFunctor(const tpmc_test::Grid<dim>& grid, int angleDegree,
                             tpmc::AlgorithmType algorithmType)
      : grid_(grid)
      , mc33_(algorithmType)
  {
    // precompute normal and centers
    field_type angleRad = angleDegree / 180.0 * M_PI;
    firstCenter_ << 0.5, 0.5, 0.375;
    firstNormal_ << std::cos(angleRad), std::sin(angleRad), 0;
    secondCenter_ << 0.5, 0.5, 0.625;
    secondNormal_ << -std::sin(angleRad), std::cos(angleRad), 0;
  }

  // return the number of connected components for two towi levelsets
  // with a distance of x
  unsigned int operator()(field_type x) const
  {
    field_type ringRadius = 0.5 * (0.25 - x);

    // first torus
    auto firstTorus = tpmc_test::torusLevelSet(0.25, ringRadius, firstCenter_, firstNormal_);
    // second torus
    domain_type secondCenter, secondNormal;
    auto secondTorus = tpmc_test::torusLevelSet(0.25, ringRadius, secondCenter_, secondNormal_);
    // combined level sets
    auto tori = [firstTorus, secondTorus](const domain_type& x) {
      return std::min(firstTorus(x), secondTorus(x));
    };
    // for each grid vertex, store a list of element local groups (as (elementIndex,localGroupIndex)
    // pairs) the vertex belongs to
    std::vector<std::vector<std::pair<unsigned int, unsigned int> > > localGroups(
        grid_.vertexCounts().prod());
    // use a union find datastructure for global connected components
    typedef std::pair<unsigned int, short> localGroupType;
    typedef std::map<localGroupType, int> rankType;
    typedef boost::associative_property_map<rankType> rankMap;
    typedef std::map<localGroupType, localGroupType> parentType;
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
      // fetch vertex groups
      std::vector<int> vertexGroups;
      mc33_.getVertexGroups(tpmc::makeGeometryType(dim, e.cornerCount()), key,
                            std::back_inserter(vertexGroups));
      // insert local groups into set and initialize new set
      for (unsigned int i = 0; i < e.cornerCount(); ++i) {
        localGroups[grid_.cornerIndex(e, i)].push_back(std::make_pair(linear, vertexGroups[i]));
        dsets.make_set(std::make_pair(linear, vertexGroups[i]));
      }
    }
    // merge all local groups belonging to a same vertex
    for (std::size_t i = 0; i < localGroups.size(); ++i) {
      for (std::size_t j = 1; j < localGroups[i].size(); ++j) {
        dsets.union_set(localGroups[i][0], localGroups[i][j]);
      }
    }
    // the global number of connected components is now the number of disjoints sets, computable by
    // the number of representatives
    unsigned int count = 0;
    for (auto& x : parent) {
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

int main(int argc, char** argv)
{
  // parse command line parameters
  std::string inifile = tpmc_test::parseCMDLineParameters(argc,argv);

  // read configuration
  std::cout << "reading configuration " << inifile << std::endl;
  auto param = tpmc_test::parseIniFile(inifile);
  const unsigned int angleStepSize = param["angleStepSize"].to_uint();
  const unsigned int numberOfElements = param["numberOfElements"].to_uint();
  const std::string referenceFile = param["referenceFile"];
  const double bisectionTolerance = param["bisectionTolerance"].to_double();
  const double fuzzyTolerance = param["fuzzyTolerance"].to_double();
  const std::string outputFilename = param["outputFilename"];

  // define general grid properties
  const int dim = 3;
  typedef tpmc_test::Grid<dim>::domain_type domain_type;
  typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;
  domain_type low = domain_type::Zero();
  domain_type high = domain_type::Ones();

  // construct Grid
  typedef tpmc_test::Grid<dim>::dimension_type dimension_type;
  dimension_type elements = dimension_type::Constant(numberOfElements);
  tpmc_test::Grid<dim> grid(elements, low, high);
  // grid width
  const field_type h = 1.0 / numberOfElements;
  // construct bisection class
  tpmc_test::Bisection<field_type, field_type> bisection(bisectionTolerance, bisectionTolerance, 30);

  // read reference file
  std::vector<std::string> referenceValues;
  bool checkReferenceSolution = (referenceFile != "");
  if (checkReferenceSolution)
    referenceValues = tpmc_test::readFile( tpmc_test::pathInfo(inifile).first + referenceFile );
  auto reference = referenceValues.begin();

  std::ofstream output(outputFilename);
  bool success = true;
  // loop through all angles
  for (int angleDegree = 0; angleDegree <= 180; angleDegree+=angleStepSize) {
    // create functor for full tpmc
    ConnectedComponentsFunctor<dim, domain_type> ccFullTPMC(grid, angleDegree,
                                                            tpmc::AlgorithmType::fullTPMC);
    auto lambdaFullTPMC = [ccFullTPMC](field_type x) { return ccFullTPMC(x) == 3 ? -1.0 : 1.0; };
    // create functor for simple tpmc
    ConnectedComponentsFunctor<dim, domain_type> ccSimpleTPMC(grid, angleDegree,
                                                              tpmc::AlgorithmType::simpleMC);
    auto lambdaSimpleTPMC =
        [ccSimpleTPMC](field_type x) { return ccSimpleTPMC(x) == 3 ? -1.0 : 1.0; };
    // compute topology change via bisection
    field_type vFullTPMC = bisection.apply(lambdaFullTPMC, 0.8 * h, 2.5 * h);
    field_type vSimpleTPMC = bisection.apply(lambdaSimpleTPMC, 0.8 * h, 2.5 * h);
    // give visual feed back
    std::cout << "." << std::flush;
    // check result against reference file
    if (checkReferenceSolution)
    {
      std::vector<tpmc_test::ini_value> refValues = tpmc_test::ini_value(*reference).to_vector();
      if (refValues[0].to_int() != angleDegree)
        throw tpmc_test::TpmcTestException("data in reference file seems to be for a different test");
      success &= std::abs(1.0 - refValues[1].to_double()/(vFullTPMC / h))   < fuzzyTolerance;
      success &= std::abs(1.0 - refValues[2].to_double()/(vSimpleTPMC / h)) < fuzzyTolerance;
      reference++;
    }

    // output to log file
    output << angleDegree << " " << std::setprecision(15) << vFullTPMC / h << " "
           << vSimpleTPMC / h << std::endl;
  }
  std::cout << std::endl;

  if (!success)
  {
    std::cerr << "Failed to reproduce reference solution" << std::endl;
    return -1;
  }
}
