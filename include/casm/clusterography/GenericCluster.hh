#ifndef CASM_GenericCluster
#define CASM_GenericCluster

#include <vector>
#include "casm/misc/CASM_TMP.hh"

namespace CASM {

  template <typename ClusterType> class ClusterInvariants;

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /* -- GenericCluster Declarations ------------------------------------- */

  /// \brief A CRTP base class for a cluster of anything
  ///
  /// - Needs a CASM_TMP::traits<Derived>::Element type
  /// - Needs a CASM_TMP::traits<Derived>::InvariantsType type
  /// - Needs Derived::apply_sym_impl to be implemented
  ///
  /// \ingroup Clusterography
  ///
  template<typename Derived>
  class GenericCluster {

  public:

    typedef unsigned int size_type;
    typedef typename CASM_TMP::traits<Derived>::Element Element;
    typedef typename CASM_TMP::traits<Derived>::InvariantsType InvariantsType;
    typedef typename std::vector<Element>::value_type value_type;
    typedef typename std::vector<Element>::iterator iterator;
    typedef typename std::vector<Element>::const_iterator const_iterator;

    /// \brief Construct an empty GenericCluster
    GenericCluster() {}

    /// \brief Construct a GenericCluster with a range of Element
    template<typename InputIterator>
    GenericCluster(InputIterator _begin,
                   InputIterator _end) :
      m_element(_begin, _end) {}


    /// \brief Iterator to first UnitCellCoord in the cluster
    iterator begin() {
      m_invariants.reset();
      return m_element.begin();
    }

    /// \brief Iterator to first UnitCellCoord in the cluster
    const_iterator begin() const {
      return m_element.begin();
    }

    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    iterator end() {
      m_invariants.reset();
      return m_element.end();
    }

    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    const_iterator end() const {
      return m_element.end();
    }

    /// \brief Iterator to first UnitCellCoord in the cluster
    const_iterator cbegin() const {
      return m_element.cbegin();
    }

    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    const_iterator cend() const {
      return m_element.cend();
    }

    /// \brief Number of UnitCellCoords in the cluster
    size_type size() const {
      return m_element.size();
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    value_type &operator[](size_type index) {
      m_invariants.reset();
      return m_element[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type &operator[](size_type index) const {
      return m_element[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    value_type &element(size_type index) {
      m_invariants.reset();
      return m_element[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type &element(size_type index) const {
      return m_element[index];
    }

    /// \brief Access vector of elements
    std::vector<Element> &elements() {
      m_invariants.reset();
      return m_element;
    }

    /// \brief const Access vector of elements
    const std::vector<Element> &elements() const {
      return m_element;
    }

    const InvariantsType &invariants() const {
      if(!m_invariants) {
        m_invariants = notstd::make_cloneable<InvariantsType>(derived());
      }
      return *m_invariants;
    }

  protected:

    Derived &derived() {
      return *static_cast<Derived *>(this);
    }

    const Derived &derived() const {
      return *static_cast<const Derived *>(this);
    }

  private:

    std::vector<Element> m_element;

    mutable notstd::cloneable_ptr<InvariantsType> m_invariants;

  };

  /// \brief CRTP-Base cluster class to apply_sym on an element-by-element basis
  ///
  /// - Needs a CASM_TMP::traits<Derived>::Element type
  /// - Needs a CASM_TMP::traits<Derived>::InvariantsType type
  ///
  /// \ingroup Clusterography
  ///
  template<typename Derived>
  class ElementWiseSymCluster : public GenericCluster<Derived> {

  public:

    /// \brief Construct an empty ElementWiseSymCluster
    ElementWiseSymCluster() {}

    /// \brief Construct a GenericCluster with a range of Element
    template<typename InputIterator>
    ElementWiseSymCluster(InputIterator _begin,
                          InputIterator _end) :
      GenericCluster<Derived>(_begin, _end) {}

  private:

    /// \brief ElementWiseSymCluster applies symmetry element-by-element
    Derived &apply_sym_impl(const SymOp &op) {
      for(auto &e : this->derived()) {
        e.apply_sym(op);
      }
      return this->derived();
    }

  };


  /* -- ClusterSymCompare Declaration ------------------------------------- */

  template<typename Derived> class ClusterSymCompare;

  namespace CASM_TMP {

    /// \brief Traits class for any ClusterSymCompare derived class
    ///
    /// \ingroup IntegralCluster
    ///
    template<typename Derived>
    struct traits<ClusterSymCompare<Derived> > {
      typedef typename traits<Derived>::MostDerived MostDerived;
      typedef typename traits<Derived>::Element Element;
      typedef typename traits<Derived>::InvariantsType InvariantsType;
    };
  }

  /// \brief CRTP Base class for Cluster comparisons
  ///
  /// Implements:
  /// - 'invariants_compare_impl' using 'compare'
  /// - 'inter_orbit_compare_impl' using 'cluster_inter_orbit_compare'
  /// - 'apply_sym_impl' (does nothing)
  ///
  /// Does not implement:
  /// - 'prepare_impl'
  /// - 'intra_orbit_compare_impl'
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare
  ///   - ClusterSymCompare
  ///     - IntegralClusterSymCompare (implements 'intra_orbit_compare_impl')
  ///       - LocalSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - PrimPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///       - ScelPeriodicSymCompare<IntegralCluster> (implements 'prepare_impl')
  ///
  /// \ingroup Clusterography
  ///
  template<typename Derived>
  class ClusterSymCompare : public SymCompare<ClusterSymCompare<Derived> > {

  public:

    /// Element refers to Cluster, not element of Cluster
    typedef typename CASM_TMP::traits<Derived>::MostDerived MostDerived;
    typedef typename CASM_TMP::traits<Derived>::Element Element;
    typedef Element ClusterType;
    typedef typename CASM_TMP::traits<Derived>::InvariantsType InvariantsType;

    /// \brief Return tolerance
    double tol() const {
      return m_tol;
    }

  protected:

    friend class SymCompare<ClusterSymCompare<Derived> >;

    /// \brief Constructor
    ///
    /// \param tol Tolerance for inter_orbit_compare of site-to-site distances
    ///
    ClusterSymCompare(double tol):
      SymCompare<ClusterSymCompare<Derived> >(),
      m_tol(tol) {

    }

    // intra_orbit_compare_impl : implement in Derived

    /// \brief Orders 'prepared' elements in the same orbit
    ///
    /// - Returns 'true' to indicate A < B
    /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
    /// - Assumes elements are 'prepared' before being compared
    /// Implementation:
    /// - First compares by number of sites in cluster
    /// - Then compare all displacements, from longest to shortest
    bool invariants_compare_impl(const InvariantsType &A, const InvariantsType &B) const {
      return compare(A, B, tol());
    }

    /// \brief Orders orbit prototypes, breaking invariants_compare ties
    ///
    /// - Returns 'true' to indicate A < B
    /// - Equivalence is indicated by \code !compare(A,B) && !compare(B,A) \endcode
    /// - Assumes elements are in canonical form
    ///
    /// Implementation:
    /// - First, check invariants_compare
    /// - Break ties with, std::lexicographical_compare of elements in A and B
    bool inter_orbit_compare_impl(const ClusterType &A, const ClusterType &B) const {

      // first compare invariants
      if(this->invariants_compare(A.invariants(), B.invariants())) {
        return true;
      }
      if(this->invariants_compare(B.invariants(), A.invariants())) {
        return false;
      }

      // next lexicographical_compare of Element in A and B
      return this->intra_orbit_compare(A, B);
    }

    /// \brief Apply symmetry to this
    ///
    /// - Affects no change
    void apply_sym_impl(const SymOp &op) {
      return;
    }

  private:

    double m_tol;

  };

}

#endif