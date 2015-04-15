#include <iostream>
#include <boost/range/algorithm/transform.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <Eigen/Dense>
#include <tpmc/marchingcubes.hh>
#include <tpmc/thresholdfunctor.hh>
#include <tpmc/fieldtraits.hh>
#include "grid.hh"
#include "geometry.hh"
#include "hexahedron.hh"
#include "vtkwriter.hh"
#include "timer.hh"
#include "levelsets.hh"
#include "quadrature.hh"

namespace tpmc
{
  template <class T, int dim>
  struct FieldTraits<Eigen::Matrix<T, dim, 1> >
  {
    typedef T field_type;
  };
}

void test() {
  const int dim = 3;
  const tpmc_test::QuadratureRule<double,3>& rule = tpmc_test::QuadratureRules<double,3>::get(2);
  typedef tpmc_test::QuadratureRule<double,3>::point_type::domain_type domain_type;
  auto func = [](domain_type x) { return x(0)*x(1)*x(2); };
  double sum = 0;
  for (auto qp: rule) {
    sum += func(qp.position())*qp.weight();
  }
  std::cout << sum << "\n";
  std::vector<domain_type> coords(4);
  coords[0] << 0,0,0;
  coords[1] << 1,0,0;
  coords[2] << 0,1,0;
  coords[3] << 1,2,0;
  std::cout << "single area: " << tpmc_test::area<double>(coords.begin(), coords.end()) << "\n";
}

int main(int argc, char** argv)
{
  test();

  tpmc_test::Timer timer;
  const int dim = 3;

  // parse config file
  const std::string configFile = argv[1];
  boost::property_tree::ptree tree;
  boost::property_tree::read_json(configFile, tree);

  // construct Grid
  const double gridSize = tree.get<double>("grid.size");
  std::array<std::size_t, dim> elements;
  elements.fill(tree.get<int>("grid.elements"));
  typedef tpmc_test::Grid<dim>::domain_type domain_type;
  tpmc_test::Grid<dim> grid(elements, domain_type::Constant(-gridSize),
                            domain_type::Constant(gridSize));

  // construct level set
  typedef tpmc::FieldTraits<domain_type>::field_type field_type;
  typedef std::function<field_type(domain_type)> LevelSet;
  std::vector<LevelSet> levelSets;
  tpmc_test::constructLevelSets<domain_type>(tree.get_child("levelSets"),
                                             std::back_inserter(levelSets));
  boost::optional<boost::property_tree::ptree &> combinedTree = tree.get_child_optional("combinedLevelSet");
  auto levelSet = combinedTree
                      ? tpmc_test::constructCombinedLevelSet<domain_type>(levelSets, *combinedTree)
                      : levelSets[0];

  // construct marching cubes object
  tpmc::MarchingCubes<field_type, dim, domain_type> mc33;

  std::vector<std::vector<domain_type> > interface;

  field_type area = 0.0;

  // loop over all elements
  for (auto e : grid.elements()) {
    // compute corner values
    std::vector<field_type> cornerValues;
    boost::transform(e.corners(), std::back_inserter(cornerValues), levelSet);
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
    // transform faces to global coordinates and add them to global list
    auto toGlobal = [&](int c) { return c < 0 ? e.corner(-c - 1) : e.global(vertices[c]); };
    for (const auto& face : faces) {
      std::vector<domain_type> global;
      std::transform(face.begin(), face.end(), std::back_inserter(global), toGlobal);
      area += tpmc_test::area<field_type>(global.begin(), global.end());
      interface.push_back(global);
    }
  }
  std::cout << "interface reconstruction: " << timer.lap().count() << "s\n";
  tpmc_test::VTKWriter<dim>().write("tpmc_test.vtk", interface.begin(), interface.end());
  std::cout << "vtk output: " << timer.lap().count() << "s\n";
  std::cout << "total: " << timer.total().count() << "s\n";
  std::cout << "area: " << std::abs(area-4*M_PI) << "\n";
}
