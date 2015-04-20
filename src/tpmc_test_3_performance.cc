#include <iostream>
#include <iomanip>
#include <numeric>
#include <random>
#include <boost/range/algorithm/transform.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <Eigen/Dense>
#include <tpmc/marchingcubes.hh>
#include <tpmc/fieldtraits.hh>
#include "timer.hh"
#include "grid.hh"
#include "geometry.hh"

namespace tpmc
{
  template <class T, int dim>
  struct FieldTraits<Eigen::Matrix<T, dim, 1> >
  {
    typedef T field_type;
  };
}

// simple functor which computes the time needed for calculating the keys for all elements
template <int dim, class domain_type>
class KeyTimeFunctor
{
public:
  typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;

  KeyTimeFunctor(const tpmc_test::Grid<dim>& grid, tpmc::AlgorithmType algorithmType)
      : grid_(grid)
      , mc33_(algorithmType)
  {
  }

  double operator()(const std::vector<field_type>& data) const
  {
    double time = 0;
    // loop through all elements
    for (const auto& element : grid_.elements()) {
      // extract corner values of current elements
      std::vector<field_type> cornerValues;
      for (unsigned int i = 0; i < element.cornerCount(); ++i) {
        cornerValues.push_back(data[grid_.cornerIndex(element, i)]);
      }
      // calculate key
      tpmc_test::Timer timer;
      std::size_t key = mc33_.getKey(cornerValues.begin(), cornerValues.end());
      // add current time to total count
      time += timer.total().count();
    }
    return time;
  }

private:
  tpmc_test::Grid<dim> grid_;
  tpmc::MarchingCubes<field_type, dim, domain_type> mc33_;
};

int main()
{
  // define general grid properties
  const int dim = 3;
  const std::vector<unsigned int> numbersOfElements{ { 16, 32, 64, 128, 256 } };
  const unsigned int numberOfRandomRuns = 30;
  const unsigned int numberOfRunsPerDataset = 30;
  typedef tpmc_test::Grid<dim>::domain_type domain_type;
  typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;
  domain_type low = domain_type::Zero();
  domain_type high = domain_type::Ones();

  const unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<field_type> distribution(-1.0, 1.0);

  std::map<unsigned int, std::vector<field_type> > timeRatios;

  for (auto numberOfElements : numbersOfElements) {
    // construct Grid
    typedef tpmc_test::Grid<dim>::dimension_type dimension_type;
    dimension_type elements = dimension_type::Constant(numberOfElements);
    tpmc_test::Grid<dim> grid(elements, low, high);
    // construct functors
    KeyTimeFunctor<dim, domain_type> keyTimeFullTPMC(grid, tpmc::AlgorithmType::fullTPMC);
    KeyTimeFunctor<dim, domain_type> keyTimeSimpleTPMC(grid, tpmc::AlgorithmType::simpleMC);

    std::vector<field_type> data(grid.vertexCounts().prod());
    tpmc_test::Timer timer;
    for (unsigned int i = 0; i < numberOfRandomRuns; ++i) {
      // generate random data
      std::generate(data.begin(), data.end(), [&]() { return distribution(generator); });
      // calculate time for fulltpmc and simpletpmc
      field_type timeFullTPMC = 0.0;
      field_type timeSimpleTPMC = 0.0;
      for (unsigned int j = 0; j < numberOfRunsPerDataset; ++j) {
        timeFullTPMC += keyTimeFullTPMC(data);
        timeSimpleTPMC += keyTimeSimpleTPMC(data);
      }
      timeFullTPMC /= numberOfRunsPerDataset;
      timeSimpleTPMC /= numberOfRunsPerDataset;
      // add time ratio to output
      timeRatios[numberOfElements].push_back(timeFullTPMC / timeSimpleTPMC);
    }
    std::cout << "time for " << numberOfRunsPerDataset << " runs of " << numberOfRandomRuns
              << " random datasets with " << numberOfElements
              << " elements: " << timer.total().count() << "s\n";
  }
  // output statistics
  std::cout << "N min max mean std\n";
  for (auto x : timeRatios) {
    field_type min
        = std::accumulate(x.second.begin(), x.second.end(), std::numeric_limits<field_type>::max(),
                          [](field_type a, field_type b) { return std::min(a, b); });
    field_type max
        = std::accumulate(x.second.begin(), x.second.end(), std::numeric_limits<field_type>::min(),
                          [](field_type a, field_type b) { return std::max(a, b); });
    field_type mean = std::accumulate(x.second.begin(), x.second.end(), 0.0) / x.second.size();
    field_type st = 0.0;
    for (auto v : x.second) {
      st += (v - mean) * (v - mean);
    }
    st = std::sqrt(st / (x.second.size() - 1));
    std::cout << x.first << " " << min << " " << max << " " << mean << " " << st << "\n";
  }
}
