#ifndef TPMC_TEST_VTKWRITER_HH
#define TPMC_TEST_VTKWRITER_HH

#include <fstream>
#include <iomanip>
#include <boost/range/iterator_range.hpp>
#include "grid.hh"

namespace tpmc_test
{
  template <int dim>
  class VTKWriter
  {
  public:
    template <class InputIterator>
    void write(const std::string& filename, InputIterator begin, InputIterator end)
    {
      // count coordinates
      std::size_t count = 0;
      for (const auto& element : boost::make_iterator_range(begin, end)) {
        count += element.size();
      }
      // write header
      std::ofstream stream(filename);
      stream << "# vtk DataFile Version 2.0\n"
             << "tpmc\n"
             << "\n"
             << "ASCII\n"
             << "DATASET UNSTRUCTURED_GRID\n"
             << "POINTS " << count << " float\n";
      // write coordinates
      for (const auto& element : boost::make_iterator_range(begin, end)) {
        for (const auto& x : element) {
          for (int d = 0; d < dim; ++d) {
            stream << std::setprecision(20) << (d > 0 ? " " : "") << x[d];
          }
          for (int d = dim; d < 3; ++d)
            stream << " " << 0;
          stream << "\n";
        }
      }
      // write elements
      std::size_t elementCount = std::distance(begin, end);
      stream << "CELLS " << elementCount << " " << count + elementCount << "\n";
      int offset = 0;
      for (const auto& element : boost::make_iterator_range(begin, end)) {
        stream << element.size();
        for (int j = 0; j < element.size(); ++j) {
          stream << " " << offset + j;
        }
        stream << "\n";
        offset += element.size();
      }
      stream << "CELL_TYPES " << elementCount << "\n";
      for (const auto& element : boost::make_iterator_range(begin, end)) {
        int key = 0;
        switch (element.size()) {
        case 2:
          key = 3;
          break; // line
        case 3:
          key = 5;
          break; // triangle
        case 4:
          key = 8;
          break; // square
        }
        stream << key << "\n";
      }
      // write cell index data
    }
  };
}

#endif // TPMC_TEST_VTKWRITER_HH
