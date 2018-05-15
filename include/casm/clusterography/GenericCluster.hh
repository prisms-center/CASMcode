#ifndef CASM_GenericCluster
#define CASM_GenericCluster

#include <vector>

#include "casm/container/algorithm.hh"
#include "casm/clusterography/ClusterDecl.hh"
#include "casm/symmetry/SymCompare.hh"
#include "casm/symmetry/PermuteIterator.hh"
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
  /// - Needs a traits<MostDerived>::Element type
  /// - Needs a traits<MostDerived>::InvariantsType type
  /// - Requires implementations:
  ///   - bool operator<(Element A, Element B);
  ///   - std::vector<Element>& MostDerived::elements_impl();
  ///   - const std::vector<Element>& MostDerived::elements_impl() const;
  /// - Optional implementations:
  ///   - MostDerived& MostDerived::sort_impl();
  ///   - bool MostDerived::is_sorted_impl() const;
  /// - _Base must inherit from CRTP_Base<MostDerived>
  ///
  /// \ingroup Clusterography
  ///
  template<typename _Base>
  class GenericCluster : public SymComparable<Comparisons<_Base>> {

  public:

    typedef SymComparable<Comparisons<_Base>> Base;
    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    typedef typename traits<MostDerived>::Element Element;
    typedef typename traits<MostDerived>::InvariantsType InvariantsType;
    typedef typename traits<MostDerived>::size_type size_type;

    typedef typename std::vector<Element>::value_type value_type;
    typedef typename std::vector<Element>::iterator iterator;
    typedef typename std::vector<Element>::const_iterator const_iterator;

    /// \brief Iterator to first UnitCellCoord in the cluster
    iterator begin() {
      this->reset_invariants();
      return derived().elements().begin();
    }

    /// \brief Iterator to first UnitCellCoord in the cluster
    const_iterator begin() const {
      return derived().elements().begin();
    }

    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    iterator end() {
      this->reset_invariants();
      return derived().elements().end();
    }

    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    const_iterator end() const {
      return derived().elements().end();
    }

    /// \brief Iterator to first UnitCellCoord in the cluster
    const_iterator cbegin() const {
      return derived().elements().cbegin();
    }

    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    const_iterator cend() const {
      return derived().elements().cend();
    }

    /// \brief Number of UnitCellCoords in the cluster
    size_type size() const {
      return derived().elements().size();
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    value_type &operator[](size_type index) {
      this->reset_invariants();
      return derived().element(index);
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type &operator[](size_type index) const {
      return derived().element(index);
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    value_type &element(size_type index) {
      this->reset_invariants();
      return derived().elements()[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type &element(size_type index) const {
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

    bool operator<(const GenericCluster &B) const {
      if(size() != B.size()) {
        return size() < B.size();
      }
      return lexicographical_compare(begin(), end(), B.begin(), B.end());
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

  };


  /*
  template<typename Derived>
  struct traits<ElementWiseSymCluster<Derived> > {
    typedef typename traits<Derived>::MostDerived MostDerived;
    typedef typename traits<Derived>::Element Element;
    typedef typename traits<Derived>::InvariantsType InvariantsType;
  };
  */

  /// \brief CRTP-Base cluster class to apply_sym on an element-by-element basis
  ///
  /// - Requires:
  ///   - Element& Element::apply_sym(const SymOp& el);
  ///
  /// \ingroup Clusterography
  ///
  template<typename Base>
  class ElementWiseSymApply : public Base {

  public:

    typedef typename Base::MostDerived MostDerived;
    using Base::derived;

    /// \brief ElementWiseSymCluster applies symmetry element-by-element
    MostDerived &apply_sym(const SymOp &op) {
      for(auto &e : derived()) {
        e.apply_sym(op);
      }
      return derived();
    }

    /// \brief ElementWiseSymCluster applies symmetry element-by-element
    MostDerived &apply_sym(const PermuteIterator &it) {
      return apply_sym(it.sym_op());
    }

  };

}

#endif
