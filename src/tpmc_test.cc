#include <iostream>
#include <boost/range/algorithm/transform.hpp>
#include <Eigen/Dense>
#include <tpmc/marchingcubes.hh>
#include <tpmc/thresholdfunctor.hh>
#include <tpmc/fieldtraits.hh>
#include "grid.hh"
#include "geometry.hh"
#include "hexahedron.hh"
#include "vtkwriter.hh"
#include "timer.hh"

int ipow(int base, int exp)
{
  int result = 1;
  while (exp) {
    if (exp & 1)
      result *= base;
    exp >>= 1;
    base *= base;
  }

  return result;
}

namespace tpmc
{
  template <class T, int dim>
  struct FieldTraits<Eigen::Matrix<T, dim, 1> >
  {
    typedef T field_type;
  };
}

// tangle cube
template <class domain_type>
std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
tangleCubeLevelSet(const domain_type& center)
{
  typedef typename tpmc::FieldTraits<domain_type>::field_type field_type;
  return [=](const domain_type& p) {
    field_type x = p(0) - center(0);
    field_type xx = x * x;
    field_type y = p(1) - center(1);
    field_type yy = y * y;
    field_type z = p(2) - center(2);
    field_type zz = z * z;
    return xx * xx - 5 * xx + yy * yy - 5 * yy + zz * zz - 5 * zz + 11.8;
  };
}

template <class domain_type>
std::function<typename tpmc::FieldTraits<domain_type>::field_type(domain_type)>
sphereLevelSet(const domain_type& center)
{
  return [=](const domain_type& p) {
    return (p-center).norm()-1.0;
  };
}

int main(int argc, char** argv)
{
  tpmc_test::Timer timer;
  const int dim = 3;
  const int elementCount = std::stoi(argv[1]);
  const double a = std::stod(argv[2]);
  std::array<std::size_t, dim> elements;
  elements.fill(elementCount);

  // construct Grid
  typedef tpmc_test::Grid<dim>::domain_type domain_type;
  tpmc_test::Grid<dim> grid(elements, domain_type::Constant(-a), domain_type::Constant(a));
  typedef domain_type::Scalar field_type;

  // construct level set
  const domain_type levelSetCenter = domain_type::Constant(0);
  //auto levelSet = tangleCubeLevelSet<domain_type>(levelSetCenter);
  auto levelSet = sphereLevelSet<domain_type>(levelSetCenter);

  // construct marching cubes object
  tpmc::MarchingCubes<field_type, dim, domain_type, tpmc::ThresholdFunctor<field_type> > mc33;

  std::vector<std::vector<domain_type> > interface;

  field_type area = 0.0;

  // loop over all elements
  for (auto e : grid.elements()) {
    // compute corner values
    std::vector<field_type> cornerValues;
    boost::transform(e.corners(), std::back_inserter(cornerValues), levelSet);
    // fetch key
    std::size_t key = mc33.getKey(cornerValues.begin(), cornerValues.end());
    // compute vertices
    std::vector<domain_type> vertices;
    mc33.getVertices(cornerValues.begin(), cornerValues.end(), key, std::back_inserter(vertices));
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
  //std::cout << "interface reconstruction: " << timer.lap().count() << "s\n";
  //tpmc_test::VTKWriter<dim>().write("tpmc_test.vtk", interface.begin(), interface.end());
  //std::cout << "vtk output: " << timer.lap().count() << "s\n";
  //std::cout << "total: " << timer.total().count() << "s\n";
  std::cout << elementCount << " " << std::abs(area-4.0*M_PI) << "\n";
}
