#include <iostream>
#include <numeric>
#include <boost/range/algorithm/transform.hpp>
#include <Eigen/Dense>
#include <tpmc/marchingcubes.hh>
#include <tpmc/fieldtraits.hh>
#include "grid.hh"
#include "geometry.hh"
#include "levelsets.hh"

namespace tpmc
{
  template <class T, int dim>
  struct FieldTraits<Eigen::Matrix<T, dim, 1> >
  {
    typedef T field_type;
  };
}

int main()
{
  // seed random generator
  srand(time(0));

  // define general grid properties
  const int dim = 3;
  const std::vector<unsigned int> numbersOfElements{ { 16, 32, 64, 128, 256 } };
  const unsigned int numberOfRandomShifts = 30;
  typedef tpmc_test::Grid<dim>::domain_type domain_type;
  domain_type low = domain_type::Zero();
  domain_type high = domain_type::Ones();

  // construct level sets
  typedef tpmc::FieldTraits<domain_type>::field_type field_type;
  // first torus
  domain_type firstCenter, firstNormal;
  firstCenter << 0.375, 0.5, 0.5;
  firstNormal << 0, 1, 0;
  auto firstTorus = tpmc_test::torusLevelSet(0.25, 0.075, firstCenter, firstNormal);
  // second torus
  domain_type secondCenter, secondNormal;
  secondCenter << 0.625, 0.5, 0.5;
  secondNormal << 0, 0, 1;
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

  // loop over all numberOfElements
  for (auto numberOfElements : numbersOfElements) {
    // construct grid
    std::array<std::size_t, dim> elements;
    elements.fill(numberOfElements);
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
  }
  // output statistics
  std::cout << "N min max mean std\n";
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
    std::cout << x.first << " " << min << " " << max << " " << mean << " " << st << "\n";
  }
}
