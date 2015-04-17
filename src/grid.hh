#ifndef TPMC_TEST_GRID_HH
#define TPMC_TEST_GRID_HH

#include <array>
#include <Eigen/Dense>
#include <boost/iterator/iterator_facade.hpp>
#include "hexahedron.hh"

namespace tpmc_test
{
  template <int dim>
  class Grid;

  template <int dim>
  class GridIndex
  {
  public:
    typedef Eigen::Matrix<unsigned int, dim,1> index_type;

    explicit GridIndex(const index_type& counts)
        : counts_(counts), current_(index_type::Constant(0))
    {
    }

    GridIndex(const index_type& counts, const index_type& current)
        : counts_(counts)
        , current_(current)
    {
    }

    const index_type& current() const { return current_; }

    const index_type& counts() const { return counts_; }

    GridIndex& operator++()
    {
      ++current_(0);
      for (int i = 0; i < dim - 1 && current_(i) == counts_(i); ++i) {
        current_(i) = 0;
        ++current_(i + 1);
      }
    }

    GridIndex operator++(int)
    {
      GridIndex old(*this);
      ++(*this);
      return old;
    }

    GridIndex& operator--()
    {
      int i = 0;
      for (int i = 0; i < dim && current_(i) == 0; ++i) {
        current_(i) = counts_(i) - 1;
      }
      if (i < dim) {
        --current_(i);
      }
    }

    GridIndex operator--(int)
    {
      GridIndex old(*this);
      --(*this);
      return old;
    }

    void advance(unsigned int offset)
    {
      current_(0) += offset;
      for (int i = 0; i < dim - 1 && current_(i) >= counts_(i); ++i) {
        while (current_(i) >= counts_(i)) {
          current_(i) -= counts_(i);
          ++current_(i + 1);
        }
      }
    }

    unsigned int linear() const
    {
      unsigned int result(current_(dim - 1));
      for (int i = dim - 1; i > 0; --i) {
        result *= counts_(i);
        result += current_(i - 1);
      }
      return result;
    }

  private:
    index_type counts_;
    index_type current_;
  };

  template <int dim>
  bool operator==(const GridIndex<dim>& a, const GridIndex<dim>& b)
  {
    return a.current() == b.current() && a.counts() == b.counts();
  }

  template <int dim>
  class Grid
  {
  public:
    typedef double field_type;
    typedef Eigen::Matrix<field_type, dim, 1> domain_type;
    typedef Eigen::Matrix<unsigned int, dim, 1> dimension_type;
    typedef Hexahedron<dim> element_type;

    class ElementIterator : public boost::iterator_facade<ElementIterator, element_type,
                                                          boost::random_access_traversal_tag,
                                                          element_type, std::ptrdiff_t>
    {
      typedef boost::iterator_facade<ElementIterator, element_type,
                                     boost::random_access_traversal_tag, element_type,
                                     std::ptrdiff_t> BaseT;

    public:
      typedef typename BaseT::value_type value_type;
      typedef typename BaseT::difference_type difference_type;
      typedef typename BaseT::reference reference;
      typedef std::random_access_iterator_tag iterator_category;

      ElementIterator(const Grid* grid, const typename GridIndex<dim>::index_type& currentIndex)
          : grid_(grid)
          , currentIndex_(grid->elementCounts(), currentIndex)
      {
      }

    private:
      void increment() { ++currentIndex_; }

      void decrement() { --currentIndex_; }

      void advance(difference_type offset) { currentIndex_.advance(offset); }

      difference_type distance_to(const ElementIterator& other) const
      {
        unsigned int lthis = currentIndex_.linear();
        unsigned int lother = other.currentIndex_.linear();
        return lthis >= lother ? lthis - lother : lother - lthis;
      }

      bool equal(const ElementIterator& other) const
      {
        return grid_ == other.grid_ && currentIndex_ == other.currentIndex_;
      }

      reference dereference() const { return grid_->element(currentIndex_.current()); }

      friend class boost::iterator_core_access;
      const Grid* grid_;
      GridIndex<dim> currentIndex_;
    };

    Grid(const dimension_type& elementCounts, const domain_type& lowerLeft = domain_type::Zero(),
         const domain_type& upperRight = domain_type::Ones())
        : elementCounts_(elementCounts)
        , vertexCounts_(elementCounts_ + dimension_type::Ones())
        , lowerLeft_(lowerLeft)
        , upperRight_(upperRight)
    {
    }

    const dimension_type& elementCounts() const { return elementCounts_; }
    const dimension_type& vertexCounts() const { return vertexCounts_; }
    const domain_type& lowerLeft() const { return lowerLeft_; }
    const domain_type& upperRight() const { return upperRight_; }

    element_type element(const typename GridIndex<dim>::index_type& index) const
    {
      domain_type s
          = (upperRight_ - lowerLeft_).cwiseQuotient(elementCounts_.template cast<field_type>());
      domain_type ll = lowerLeft_ + s.cwiseProduct(index.template cast<field_type>());
      return element_type(ll, s, index);
    }

    unsigned int cornerIndex(const element_type& element, unsigned int corner) const {
      dimension_type vertexIndex(element.index());
      for (int i = 0; i<dim; ++i) {
        if (corner & (1<<i)) {
          ++vertexIndex(i);
        }
      }
      return GridIndex<dim>(vertexCounts_, vertexIndex).linear();
    }

    boost::iterator_range<ElementIterator> elements() const
    {
      typedef typename GridIndex<dim>::index_type index_type;
      index_type beginIndex = index_type::Constant(0);
      index_type endIndex = beginIndex;
      endIndex(dim - 1)= elementCounts_(dim - 1);
      return boost::iterator_range<ElementIterator>(ElementIterator(this, beginIndex),
                                                    ElementIterator(this, endIndex));
    }

  private:
    dimension_type elementCounts_;
    dimension_type vertexCounts_;
    domain_type lowerLeft_;
    domain_type upperRight_;
  };
}

#endif
