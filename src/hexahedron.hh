#ifndef TPMC_TEST_HEXAHEDRON_HH
#define TPMC_TEST_HEXAHEDRON_HH

#include <Eigen/Dense>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>

namespace tpmc_test
{
  template <int dim>
  class Hexahedron
  {
  public:
    typedef double field_type;
    typedef Eigen::Matrix<field_type, dim, 1> domain_type;

    class CornerIterator : public boost::iterator_facade<CornerIterator, domain_type,
                                                         boost::random_access_traversal_tag,
                                                         domain_type, std::ptrdiff_t>
    {
      typedef boost::iterator_facade<CornerIterator, domain_type,
                                     boost::random_access_traversal_tag, domain_type,
                                     std::ptrdiff_t> BaseT;

    public:
      typedef typename BaseT::value_type value_type;
      typedef typename BaseT::difference_type difference_type;
      typedef typename BaseT::reference reference;
      typedef std::random_access_iterator_tag iterator_category;

      CornerIterator(const Hexahedron* hexahedron, std::size_t currentIndex)
          : hexahedron_(hexahedron)
          , currentIndex_(currentIndex)
      {
      }

    private:
      void increment() { ++currentIndex_; }

      void decrement() { --currentIndex_; }

      void advance(difference_type offset) { currentIndex_ += offset; }

      difference_type distance_to(const CornerIterator& other) const
      {
        return other.currentIndex_ >= currentIndex_ ? other.currentIndex_ - currentIndex_
                                                    : currentIndex_ - other.currentIndex_;
      }

      bool equal(const CornerIterator& other) const
      {
        return hexahedron_ == other.hexahedron_ && currentIndex_ == other.currentIndex_;
      }

      reference dereference() const { return hexahedron_->corner(currentIndex_); }

      friend class boost::iterator_core_access;
      const Hexahedron* hexahedron_;
      std::size_t currentIndex_;
    };

    Hexahedron(const domain_type& lowerLeft, const domain_type& size)
        : lowerLeft_(lowerLeft)
        , size_(size)
    {
    }

    std::size_t cornerCount() const { return 1 << dim; }
    domain_type corner(std::size_t index) const
    {
      domain_type result = lowerLeft_;
      for (int i = 0; i < dim; ++i) {
        if (index & (1 << i)) {
          result[i] += size_[i];
        }
      }
      return result;
    }

    boost::iterator_range<CornerIterator> corners() const
    {
      return boost::iterator_range<CornerIterator>(CornerIterator(this, 0),
                                                   CornerIterator(this, cornerCount()));
    }

    domain_type global(const domain_type& l) const {
      return l.cwiseProduct(size_)+lowerLeft_;
    }

    const domain_type& lowerLeft() const { return lowerLeft_; }

  private:
    domain_type lowerLeft_;
    domain_type size_;
  };
}

#endif // TPMC_TEST_HEXAHEDRON_HH

