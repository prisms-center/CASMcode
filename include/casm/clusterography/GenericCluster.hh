#ifndef CASM_GenericCluster
#define CASM_GenericCluster

#include <vector>

#include "casm/container/Permutation.hh"
#include "casm/container/algorithm.hh"
#include "casm/misc/Comparisons.hh"

namespace CASM {

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /* -- GenericCluster Declarations ------------------------------------- */

  /*
  template<typename Derived>
  struct traits<GenericCluster<Derived> > {
    typedef typename traits<Derived>::MostDerived MostDerived;
    typedef typename traits<Derived>::Element Element;
    typedef typename traits<Derived>::InvariantsType InvariantsType;
  };
  */

  /// \brief A CRTP base class for a cluster of anything
  ///
  /// GenericCluster provides a standard interface for clusters of various type. In particular, it
  /// is commonly necessary to iterate over elements in the cluster, access elements in the cluster,
  /// sort/permute elements in the cluster into a canonical representation, and compare clusters.
  ///
  /// To implement this, there must be a traits class for the particular derived type of cluster
  /// with the following members:
  /// - typename traits<MostDerived>::Element
  ///
  /// The derived cluster type must implement public methods:
  /// - std::vector<Element>& MostDerived::elements();
  /// - const std::vector<Element>& MostDerived::elements() const;
  ///
  /// Optionally, the derived cluster type may specialize protected methods:
  /// - MostDerived& sort_impl();
  /// - bool is_sorted_impl() const;
  /// - Permutation sort_permutation_impl() const;
  /// - bool compare_impl(MostDerived const &B) const;
  ///
  /// To implement the CRTP pattern:
  /// - _Base must inherit from CRTPBase<MostDerived>
  ///
  /// \ingroup Clusterography
  ///
  template<typename Base>
  class GenericCluster : public Comparisons<Base> {

  public:

    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    typedef typename traits<MostDerived>::Element Element;
    typedef Index size_type;

    typedef typename std::vector<Element>::value_type value_type;
    typedef typename std::vector<Element>::iterator iterator;
    typedef typename std::vector<Element>::const_iterator const_iterator;

    /// \brief Iterator to first element in the cluster
    iterator begin() {
      return derived().elements().begin();
    }

    /// \brief Iterator to first element in the cluster
    const_iterator begin() const {
      return derived().elements().begin();
    }

    /// \brief Iterator to the past-the-last element in the cluster
    iterator end() {
      return derived().elements().end();
    }

    /// \brief Iterator to the past-the-last element in the cluster
    const_iterator end() const {
      return derived().elements().end();
    }

    /// \brief Iterator to first element in the cluster
    const_iterator cbegin() const {
      return derived().elements().cbegin();
    }

    /// \brief Iterator to the past-the-last element in the cluster
    const_iterator cend() const {
      return derived().elements().cend();
    }

    /// \brief Number of elements in the cluster
    size_type size() const {
      return derived().elements().size();
    }

    /// \brief Access an element in the cluster by index
    value_type &operator[](size_type index) {
      return derived().element(index);
    }

    /// \brief Access an element in the cluster by index
    value_type const &operator[](size_type index) const {
      return derived().element(index);
    }

    /// \brief Access an element in the cluster by index
    value_type &element(size_type index) {
      return derived().elements()[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    value_type const &element(size_type index) const {
      return derived().elements()[index];
    }

    MostDerived &sort() {
      return derived().sort_impl();
    }

    MostDerived sorted() const {
      MostDerived tmp {derived()};
      return tmp.sort();
    }

    bool is_sorted() const {
      return derived().is_sorted_impl();
    }

    Permutation sort_permutation() const {
      return derived().sort_permutation_impl();
    }

    bool operator<(MostDerived const &B) const {
      return derived().compare_impl(B);
    }

  protected:

    /// \brief Construct an empty GenericCluster
    GenericCluster() {}

    MostDerived &sort_impl() {
      std::sort(begin(), end());
      return derived();
    }

    bool is_sorted_impl() const {
      return std::is_sorted(begin(), end());
    }

    Permutation sort_permutation_impl() const {
      if(size() == 0)
        return Permutation(0);
      std::vector<Index> ind = sequence(Index(0), Index(size() - 1));

      std::sort(ind.begin(), ind.end(),
      [&](const Index & a, const Index & b) {
        return (this->derived().element(a) < this->derived().element(b));
      });

      return Permutation(std::move(ind));
    }

    bool compare_impl(MostDerived const &B) const {
      if(size() != B.size()) {
        return size() < B.size();
      }
      return lexicographical_compare(begin(), end(), B.begin(), B.end());
    }

  };
}

#endif
