#ifndef TPMC_TEST_GRID_HH
#define TPMC_TEST_GRID_HH

#include <array>
#include <Eigen/Dense>
#include <boost/iterator/iterator_facade.hpp>
#include "hexahedron.hh"

namespace tpmc_test
{
  template <int dim>
  class Grid
  {
  public:
    typedef double field_type;
    typedef Eigen::Matrix<field_type, dim, 1> domain_type;
    typedef Hexahedron<dim> element_type;
    typedef std::array<std::size_t, dim> element_index_type;

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

      ElementIterator(const Grid* grid, const element_index_type& currentIndex)
          : grid_(grid)
          , currentIndex_(currentIndex)
      {
      }

    private:
      void increment()
      {
        ++currentIndex_[0];
        for (int i = 0; i < dim - 1 && currentIndex_[i] == grid_->elementCounts()[i]; ++i) {
          currentIndex_[i] = 0;
          ++currentIndex_[i + 1];
        }
      }

      void decrement()
      {
        int i = 0;
        for (int i = 0; i < dim && currentIndex_[i] == 0; ++i) {
          currentIndex_[i] = grid_->elementCounts()[i] - 1;
        }
        if (i < dim) {
          --currentIndex_[i];
        }
      }

      void advance(difference_type offset)
      {
        currentIndex_[0] += offset;
        for (int i = 0; i < dim - 1 && currentIndex_[i] >= grid_->elementCounts()[i]; ++i) {
          while (currentIndex_[i] >= grid_->elementCounts()[i]) {
            currentIndex_[i] -= grid_->elementCounts()[i];
            ++currentIndex_[i + 1];
          }
        }
      }

      difference_type distance_to(const ElementIterator& other) const
      {
        std::size_t lthis = linear(currentIndex_);
        std::size_t lother = linear(other.currentIndex_);
        return lthis >= lother ? lthis - lother : lother - lthis;
      }

      bool equal(const ElementIterator& other) const
      {
        return grid_ == other.grid_ && currentIndex_ == other.currentIndex_;
      }

      reference dereference() const { return grid_->element(currentIndex_); }

      friend class boost::iterator_core_access;
      const Grid* grid_;
      std::array<std::size_t, dim> currentIndex_;
    };

    Grid(const element_index_type& elementCounts,
         const domain_type& lowerLeft = domain_type::Zero(),
         const domain_type& upperRight = domain_type::Ones())
        : elementCounts_(elementCounts)
        , lowerLeft_(lowerLeft)
        , upperRight_(upperRight)
    {
    }

    const element_index_type& elementCounts() const { return elementCounts_; }
    const domain_type& lowerLeft() const { return lowerLeft_; }
    const domain_type& upperRight() const { return upperRight_; }

    std::size_t linear(const element_index_type& index) const
    {
      std::size_t result(index[dim - 1]);
      for (int i = dim - 1; i > 0; --i) {
        result *= elementCounts_[i];
        result += index[i - 1];
      }
      return result;
    }

    element_type element(const element_index_type& index) const
    {
      domain_type s = upperRight_ - lowerLeft_;
      for (int i = 0; i < dim; ++i) {
        s(i) /= elementCounts_[i];
      }
      domain_type ll = lowerLeft_;
      for (int i = 0; i < dim; ++i) {
        ll(i) += s(i) * index[i];
      }
      return element_type(ll, s);
    }

    boost::iterator_range<ElementIterator> elements() const
    {
      element_index_type beginIndex;
      std::fill(beginIndex.begin(), beginIndex.end(), 0);
      element_index_type endIndex = beginIndex;
      endIndex[dim - 1] = elementCounts_[dim - 1];
      return boost::iterator_range<ElementIterator>(ElementIterator(this, beginIndex),
                                                    ElementIterator(this, endIndex));
    }

  private:
    element_index_type elementCounts_;
    domain_type lowerLeft_;
    domain_type upperRight_;
  };
}

#endif
