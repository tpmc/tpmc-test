#include <map>
#include <iostream>
#include <Eigen/Dense>
#include <tpmc/fieldtraits.hh>
#include <tpmc/marchingcubes.hh>
#include "grid.hh"
#include "timer.hh"
#include "geometry.hh"

namespace tpmc
{
  template <class T, int dim>
  struct FieldTraits<Eigen::Matrix<T, dim, 1> >
  {
    typedef T field_type;
  };
}

struct TimeResult
{
  // the time for computing the key for each grid element
  double key;
  // the time for computing the surface area of the interface
  double iteration;
  // total time, including grid iteration
  double total;
  // area of the interface
  double area;

  double elements;

  TimeResult& operator+=(const TimeResult& other)
  {
    key += other.key;
    iteration += other.iteration;
    total += other.total;
    area += other.area;
    elements += other.elements;
    return *this;
  }

  TimeResult& operator/=(double v)
  {
    key /= v;
    iteration /= v;
    total /= v;
    area /= v;
    elements /= v;
    return *this;
  }
};

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

  TimeResult operator()(const std::vector<field_type>& data) const
  {
    TimeResult result{ 0.0, 0.0, 0.0, 0.0, 0 };
    tpmc_test::Timer timer;
    tpmc_test::Timer iterationTimer;
    // loop through all elements
    for (const auto& element : grid_.elements()) {
      // extract corner values of current elements
      std::vector<field_type> cornerValues;
      for (unsigned int i = 0; i < element.cornerCount(); ++i) {
        cornerValues.push_back(data[grid_.cornerIndex(element, i)]);
      }
      tpmc_test::Timer keyTimer;
      // calculate key
      std::size_t key = mc33_.getKey(cornerValues.begin(), cornerValues.end());
      result.key += keyTimer.total().count();

      result.iteration += iterationTimer.lap().count();

      // compute surface area
      std::vector<domain_type> vertices;
      mc33_.getVertices(std::begin(cornerValues), std::end(cornerValues), key,
                        std::back_inserter(vertices));
      std::vector<std::vector<int> > faces;
      mc33_.getElements(tpmc::makeGeometryType(dim, element.cornerCount()), key,
                        tpmc::ReconstructionType::Interface, std::back_inserter(faces));
      result.elements += faces.size();
      auto toGlobal =
          [&](int c) { return c < 0 ? element.corner(-c - 1) : element.global(vertices[c]); };
      for (const auto& face : faces) {
        // transform face to global coordinates
        std::vector<domain_type> global;
        std::transform(face.begin(), face.end(), std::back_inserter(global), toGlobal);
        // add the area of global
        result.area += tpmc_test::area<field_type>(global.begin(), global.end());
      }

      iterationTimer.lap();
    }
    result.total = timer.total().count();
    return result;
  }

  const tpmc::MarchingCubes<field_type, dim, domain_type>& getMC() const { return mc33_; }

private:
  tpmc_test::Grid<dim> grid_;
  tpmc::MarchingCubes<field_type, dim, domain_type> mc33_;
};

template <class Map>
void printStatistics(std::ostream& output, const Map& m)
{
  typedef typename decltype(m.begin()->second)::value_type field_type;
  output << "N min max mean std\n";
  for (auto x : m) {
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
    output << x.first << " " << min << " " << max << " " << mean << " " << st << "\n";
  }
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Error: no output file provided. call:\n" << argv[0] << " <outputfilename>\n";
    return -1;
  }
  const std::string outputfilename = argv[1];
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

  std::map<unsigned int, std::vector<field_type> > timeRatiosKey;
  std::map<unsigned int, std::vector<field_type> > timeRatiosIteration;
  std::map<unsigned int, std::vector<field_type> > timeRatiosIntegration;
  std::map<unsigned int, std::vector<field_type> > elementCountRatios;
  std::map<unsigned int, unsigned long> keyGenerations;
  std::map<unsigned int, unsigned long> faceTests;
  std::map<unsigned int, unsigned long> centerTests;

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
      TimeResult acctimeResultFullTPMC{ 0.0, 0.0, 0.0 , 0.0, 0.0};
      TimeResult acctimeResultSimpleTPMC{ 0.0, 0.0, 0.0 , 0.0, 0.0};
      for (unsigned int j = 0; j < numberOfRunsPerDataset; ++j) {
        acctimeResultFullTPMC += keyTimeFullTPMC(data);
        acctimeResultSimpleTPMC += keyTimeSimpleTPMC(data);
      }
      acctimeResultFullTPMC /= numberOfRunsPerDataset;
      acctimeResultSimpleTPMC /= numberOfRunsPerDataset;
      // add time ratio to output
      timeRatiosKey[numberOfElements].push_back(acctimeResultFullTPMC.key
                                                / acctimeResultSimpleTPMC.key);
      timeRatiosIteration[numberOfElements].push_back(acctimeResultFullTPMC.iteration
                                                      / acctimeResultSimpleTPMC.iteration);
      timeRatiosIntegration[numberOfElements].push_back(acctimeResultFullTPMC.total
                                                        / acctimeResultSimpleTPMC.total);
      elementCountRatios[numberOfElements].push_back(static_cast<double>(acctimeResultFullTPMC.elements)
                                                     / acctimeResultSimpleTPMC.elements);
    }
    keyGenerations[numberOfElements] = keyTimeFullTPMC.getMC().profKeyGenerations();
    faceTests[numberOfElements] = keyTimeFullTPMC.getMC().profFaceTests();
    centerTests[numberOfElements] = keyTimeFullTPMC.getMC().profCenterTests();

    std::cout << "time for " << numberOfRunsPerDataset << " runs of " << numberOfRandomRuns
              << " random datasets with " << numberOfElements
              << " elements: " << timer.total().count() << "s\n";
  }
  // output statistics
  std::ofstream output(outputfilename);
  output << "time ratios key\n";
  printStatistics(output, timeRatiosKey);
  output << "\n";
  output << "time ratios iteration\n";
  printStatistics(output, timeRatiosIteration);
  output << "\n";
  output << "time ratios total\n";
  printStatistics(output, timeRatiosIntegration);
  output << "\n";
  output << "element count ratios\n";
  printStatistics(output, elementCountRatios);
  output << "\n";
  output << "N keyGen faceTests relFaceTests centerTests relCenterTests\n";
  for (auto numberOfElements : numbersOfElements) {
    output << numberOfElements << " " << keyGenerations[numberOfElements] << " "
           << faceTests[numberOfElements] << " "
           << static_cast<double>(faceTests[numberOfElements]) / keyGenerations[numberOfElements]
           << " " << centerTests[numberOfElements] << " "
           << static_cast<double>(centerTests[numberOfElements]) / keyGenerations[numberOfElements]
           << "\n";
  }
}
