#ifndef CASM_GenericCluster
#define CASM_GenericCluster

#include <vector>

namespace CASM {

  template <typename ClusterType> class ClusterInvariants;

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /* -- GenericCluster Declarations ------------------------------------- */

  template<typename Derived> class GenericCluster;

  /// \brief Traits class for GenericCluster
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename Derived>
  struct traits<GenericCluster<Derived> > {
    typedef typename traits<Derived>::MostDerived MostDerived;
    typedef typename traits<Derived>::Element Element;
    typedef typename traits<Derived>::InvariantsType InvariantsType;
  };

  /// \brief A CRTP base class for a cluster of anything
  ///
  /// - Needs a traits<Derived>::Element type
  /// - Needs a traits<Derived>::InvariantsType type
  /// - Needs Derived::apply_sym_impl to be implemented
  ///
  /// \ingroup Clusterography
  ///
  template<typename Derived>
  class GenericCluster : public SymComparable<GenericCluster<Derived> > {

  public:

    typedef unsigned int size_type;
    typedef typename traits<Derived>::MostDerived MostDerived;
    typedef typename traits<Derived>::Element Element;
    typedef typename traits<Derived>::InvariantsType InvariantsType;
    typedef typename std::vector<Element>::value_type value_type;
    typedef typename std::vector<Element>::iterator iterator;
    typedef typename std::vector<Element>::const_iterator const_iterator;

    /// \brief Iterator to first UnitCellCoord in the cluster
    iterator begin() {
      this->_reset_invariants();
      return m_element.begin();
    }

    /// \brief Iterator to first UnitCellCoord in the cluster
    const_iterator begin() const {
      return m_element.begin();
    }

    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    iterator end() {
      this->_reset_invariants();
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
      this->_reset_invariants();
      return m_element[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type &operator[](size_type index) const {
      return m_element[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    value_type &element(size_type index) {
      this->_reset_invariants();
      return m_element[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type &element(size_type index) const {
      return m_element[index];
    }

    /// \brief Access vector of elements
    std::vector<Element> &elements() {
      this->_reset_invariants();
      return m_element;
    }

    /// \brief const Access vector of elements
    const std::vector<Element> &elements() const {
      return m_element;
    }

  protected:

    /// \brief Construct an empty GenericCluster
    GenericCluster() {}

    /// \brief Construct a GenericCluster with a range of Element
    template<typename InputIterator>
    GenericCluster(InputIterator _begin,
                   InputIterator _end) :
      m_element(_begin, _end) {}

  private:

    std::vector<Element> m_element;

  };

  template<typename Derived> class ElementWiseSymCluster;

  /// \brief Traits class for any ClusterSymCompare derived class
  ///
  /// \ingroup IntegralCluster
  ///
  template<typename Derived>
  struct traits<ElementWiseSymCluster<Derived> > {
    typedef typename traits<Derived>::MostDerived MostDerived;
    typedef typename traits<Derived>::Element Element;
    typedef typename traits<Derived>::InvariantsType InvariantsType;
  };

  /// \brief CRTP-Base cluster class to apply_sym on an element-by-element basis
  ///
  /// - Needs a traits<Derived>::Element type
  /// - Needs a traits<Derived>::InvariantsType type
  ///
  /// \ingroup Clusterography
  ///
  template<typename Derived>
  class ElementWiseSymCluster : public GenericCluster<ElementWiseSymCluster<Derived> > {

  public:

    typedef typename traits<Derived>::MostDerived MostDerived;

    /// \brief ElementWiseSymCluster applies symmetry element-by-element
    MostDerived &apply_sym(const SymOp &op) {
      this->derived().apply_sym_impl(op);
      return this->derived();
    }

  protected:

    /// \brief Construct an empty ElementWiseSymCluster
    ElementWiseSymCluster() :
      GenericCluster<ElementWiseSymCluster<Derived> >() {}

    /// \brief Construct a GenericCluster with a range of Element
    template<typename InputIterator>
    ElementWiseSymCluster(InputIterator _begin,
                          InputIterator _end) :
      GenericCluster<ElementWiseSymCluster<Derived> >(_begin, _end) {}


    /// \brief ElementWiseSymCluster applies symmetry element-by-element
    MostDerived &apply_sym_impl(const SymOp &op) {
      for(auto &e : this->derived()) {
        e.apply_sym(op);
      }
      return this->derived();
    }

  };


  /* -- ClusterSymCompare Declaration ------------------------------------- */

  template<typename Derived> class ClusterSymCompare;

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

  /// \brief CRTP Base class for Cluster comparisons
  ///
  /// Implements:
  /// - 'invariants_compare_impl' using 'compare'
  /// - 'apply_sym_impl' (does nothing)
  ///
  /// Does not implement:
  /// - 'prepare_impl'
  /// - 'compare_impl'
  ///
  /// The ClusterSymCompare hierarchy:
  /// - SymCompare
  ///   - ClusterSymCompare
  ///     - IntegralClusterSymCompare (implements 'compare_impl')
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
    typedef typename traits<Derived>::MostDerived MostDerived;
    typedef typename traits<Derived>::Element Element;
    typedef Element ClusterType;
    typedef typename traits<Derived>::InvariantsType InvariantsType;

    /// \brief Return tolerance
    double tol() const {
      return m_tol;
    }

  protected:

    friend class SymCompare<ClusterSymCompare<Derived> >;

    /// \brief Constructor
    ///
    /// \param tol Tolerance for invariants_compare of site-to-site distances
    ///
    ClusterSymCompare(double tol):
      SymCompare<ClusterSymCompare<Derived> >(),
      m_tol(tol) {

    }

    // compare_impl : implement in Derived

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
