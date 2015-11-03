#include <iostream>
#include <numeric>
#include <boost/range/algorithm/transform.hpp>
#include <Eigen/Dense>
#include <tpmc/marchingcubes.hh>
#include <tpmc/fieldtraits.hh>
#include "grid.hh"
#include "geometry.hh"
#include "levelsets.hh"
#include "timer.hh"
#include "io.hh"
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

namespace tpmc
{
  template <class T, int dim>
  struct FieldTraits<Eigen::Matrix<T, dim, 1> >
  {
    typedef T field_type;
  };
}

int main(int argc, char** argv)
{
  auto path_info = tpmc_test::pathInfo(argv[0]);
  std::string inifile =
    TOSTRING(TPMC_TEST_DIR) + path_info.second + ".ini";
  if (argc >= 2) {
    inifile = argv[1];
  }

  // read configuration
  std::cout << "reading configuration " << inifile << std::endl;
  auto param = tpmc_test::parseIniFile(inifile);
  const unsigned int numberOfLevels = param["numberOfLevels"].to_uint();
  const unsigned int numberOfRandomShifts = param["numberOfRandomShifts"].to_uint();
  const std::string referenceFile = param["referenceFile"];
  const double fuzzyTolerance = param["fuzzyTolerance"].to_double();
  const std::string outputFilename = param["outputFilename"];

  // seed random generator
  srand(time(0));
  tpmc_test::Timer timer;

  // define general grid properties
  const int dim = 3;
  std::vector<unsigned int> numbersOfElements = {{16}};
  for (unsigned int l = 1; l<numberOfLevels; l++)
      numbersOfElements.push_back(numbersOfElements.back() * 2);
  typedef tpmc_test::Grid<dim>::domain_type domain_type;
  domain_type low = domain_type::Zero();
  domain_type high = domain_type::Ones();

  // construct level sets
  typedef tpmc::FieldTraits<domain_type>::field_type field_type;
  // first torus
  domain_type firstCenter, firstNormal;
  firstCenter << 0.5, 0.5, 0.375;
  firstNormal << 1, 0, 0;
  auto firstTorus = tpmc_test::torusLevelSet(0.25, 0.075, firstCenter, firstNormal);
  // second torus
  domain_type secondCenter, secondNormal;
  secondCenter << 0.5, 0.5, 0.625;
  secondNormal << 0, 1, 0;
  auto secondTorus = tpmc_test::torusLevelSet(0.25, 0.075, secondCenter, secondNormal);
  // combined level sets
  auto tori = [firstTorus, secondTorus](const domain_type& x) {
    return std::min(firstTorus(x), secondTorus(x));
  };
  // analytic surface
  field_type referenceSurface = 8.0 * M_PI * M_PI * 0.25 * 0.075;

  // construct marching cubes class
  tpmc::MarchingCubes<field_type, dim, domain_type> mc33;

  // storage for relative errors
  std::map<unsigned int, std::vector<field_type> > relativeErrors;

  std::cout << "initialization: " << timer.lap().count() << "s\n";

  // loop over all numberOfElements
  for (auto numberOfElements : numbersOfElements) {
    // construct grid
    typedef typename tpmc_test::Grid<dim>::dimension_type dimension_type;
    dimension_type elements = dimension_type::Constant(numberOfElements);
    tpmc_test::Grid<dim> grid(elements, low, high);

    // loop over all shifts
    for (unsigned int i = 0; i < numberOfRandomShifts; ++i) {
      // generate random shift
      domain_type shift = (1.0 / numberOfElements) * domain_type::Random();
      // construct shifted level set
      auto shiftedLevelSet = [tori, shift](const domain_type& x) { return tori(x + shift); };

      field_type area = 0.0;

      // loop over all elements
      for (auto e : grid.elements()) {
        // compute corner values
        std::vector<field_type> cornerValues;
        boost::transform(e.corners(), std::back_inserter(cornerValues), shiftedLevelSet);
        // fetch key
        std::size_t key = mc33.getKey(std::begin(cornerValues), std::end(cornerValues));
        // compute vertices
        std::vector<domain_type> vertices;
        mc33.getVertices(std::begin(cornerValues), std::end(cornerValues), key,
                         std::back_inserter(vertices));
        // retrieve faces
        std::vector<std::vector<int> > faces;
        mc33.getElements(tpmc::makeGeometryType(dim, e.cornerCount()), key,
                         tpmc::ReconstructionType::Interface, std::back_inserter(faces));
        // construct global transformation
        auto toGlobal = [&](int c) { return c < 0 ? e.corner(-c - 1) : e.global(vertices[c]); };
        // loop over all faces
        for (const auto& face : faces) {
          // transform face to global coordinates
          std::vector<domain_type> global;
          std::transform(face.begin(), face.end(), std::back_inserter(global), toGlobal);
          // add the area of global
          area += tpmc_test::area<field_type>(global.begin(), global.end());
        }
      }
      // insert relative error to the global set
      relativeErrors[numberOfElements].push_back(std::abs(area - referenceSurface)
                                                 / referenceSurface);
    }

    std::cout << "tests for " << numberOfElements << " elements: " << timer.lap().count() << "s\n";
  }


  // read reference file
  std::vector<std::string> referenceValues;
  bool checkReferenceSolution = (referenceFile != "");
  if (checkReferenceSolution)
    referenceValues = tpmc_test::readFile( tpmc_test::pathInfo(inifile).first + referenceFile );
  auto reference = referenceValues.begin();

  // output statistics
  bool success = true;
  std::ofstream output(outputFilename);
  output << "N min max mean std\n";
  reference++; // skip title line
  for (auto x : relativeErrors) {
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

    // check result against reference file
    if (checkReferenceSolution)
    {
      std::vector<tpmc_test::ini_value> refValues = tpmc_test::ini_value(*reference).to_vector();
      if (refValues[0].to_uint() != x.first)
        throw tpmc_test::TpmcTestException("data in reference file seems to be for a different test");
      success &= std::abs(1.0 - refValues[1].to_double()/min)  < fuzzyTolerance;
      success &= std::abs(1.0 - refValues[2].to_double()/max)  < fuzzyTolerance;
      success &= std::abs(1.0 - refValues[3].to_double()/mean) < fuzzyTolerance;
      reference++;
    }

    output << x.first << " " << min << " " << max << " " << mean << " " << st << "\n";
  }
  output.close();

  std::cout << "statistics: " << timer.lap().count() << "s\n";
  std::cout << "total: " << timer.total().count() << "s\n";

  if (!success)
  {
    std::cerr << "Failed to reproduce reference solution" << std::endl;
    return -1;
  }
}
