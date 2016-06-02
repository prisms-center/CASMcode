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
  /// - Needs a CASM_TMP::traits<DerivedCluster>::Element type
  /// - Needs a valid ClusterInvariants<DerivedCluster> to help comparisons
  /// - Needs DerivedCluster::apply_sym to be implemented
  /// \ingroup Clusterography
  ///
  template<typename DerivedCluster>
  class GenericCluster {

  public:

    typedef unsigned int size_type;
    typedef typename CASM_TMP::traits<DerivedCluster>::Element Element;
    typedef typename std::vector<Element>::value_type value_type;
    typedef typename std::vector<Element>::iterator iterator;
    typedef typename std::vector<Element>::const_iterator const_iterator;
    typedef ClusterInvariants<DerivedCluster> InvariantsType;

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

    DerivedCluster &derived() {
      return *static_cast<DerivedCluster *>(this);
    }

    const DerivedCluster &derived() const {
      return *static_cast<const DerivedCluster *>(this);
    }

  private:

    std::vector<Element> m_element;

    mutable notstd::cloneable_ptr<InvariantsType> m_invariants;

  };

  /// \brief CRTP-Base cluster class to apply_sym on an element-by-element basis
  ///
  /// - Needs a CASM_TMP::traits<DerivedCluster>::Element type
  /// - Needs a valid ClusterInvariants<DerivedCluster> to help comparisons
  ///
  template<typename DerivedCluster>
  class ElementWiseSymCluster : public GenericCluster<DerivedCluster> {

  public:

    /// \brief Construct an empty ElementWiseSymCluster
    ElementWiseSymCluster() {}

    /// \brief Construct a GenericCluster with a range of Element
    template<typename InputIterator>
    ElementWiseSymCluster(InputIterator _begin,
                          InputIterator _end) :
      GenericCluster<DerivedCluster>(_begin, _end) {}

    /// \brief ElementWiseSymCluster applies symmetry element-by-element
    DerivedCluster &apply_sym(const SymOp &op) {
      for(auto &e : this->derived()) {
        e.apply_sym(op);
      }
      return this->derived();
    }

  };

  /// \brief Default intra orbit comparison of clusters
  ///
  /// Implements:
  /// \code
  /// std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare);
  /// \endcode
  ///
  /// \ingroup Clusterography
  ///
  template<typename ClusterType, typename Compare = std::less<typename CASM_TMP::traits<ClusterType>::Element> >
  bool cluster_intra_orbit_compare(const ClusterType &A,
                                   const ClusterType &B,
                                   Compare compare = Compare()) {
    return std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare);
  }

  /// \brief Default comparison of orbit prototypes
  ///
  /// \param A, B clusters to check for A < B
  /// \param tol Tolerance for comparison of element-to-element distances
  /// \param compare_element Compare concept functor for lexicographical comparison of cluster elements
  /// \param compare_invariants Compare concept functor for cluster invariants
  ///
  /// Compares:
  /// - cluster size
  /// - element-to-element distances (max to min)
  /// - lexicographical_compare of elements, via compare
  ///
  /// Implements:
  /// \code
  /// // first compare cluster size
  /// if(A.size() != B.size()) {
  ///   return A.size() < B.size();
  /// }
  ///
  /// // next compare invariants
  /// if( compare_invariants(ClusterInvariants<Element>(A), ClusterInvariants<Element>(B), tol) ) {
  ///   return true;
  /// }
  /// if( compare_invariants(ClusterInvariants<Element>(B), ClusterInvariants<Element>(A), tol) ) {
  ///   return false;
  /// }
  ///
  /// // next lexicographical_compare of UnitCellCoord in A and B
  /// return std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare_element);
  /// \endcode
  ///
  /// \ingroup Clusterography
  ///
  template<typename ClusterType,
           typename CompareElement = std::less<typename CASM_TMP::traits<ClusterType>::Element>,
           typename CompareInvariants = FloatCompare>
  bool cluster_inter_orbit_compare(const ClusterType &A,
                                   const ClusterType &B,
                                   double tol,
                                   CompareElement compare_element = CompareElement(),
                                   CompareInvariants compare_invariants = FloatCompare(TOL)) {

    // first compare cluster size
    if(A.size() != B.size()) {
      return A.size() < B.size();
    }

    // next compare invariants
    if(compare_invariants(A.invariants(), B.invariants())) {
      return true;
    }
    if(compare_invariants(B.invariants(), A.invariants())) {
      return false;
    }

    // next lexicographical_compare of Element in A and B
    return std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare_element);
  }

}

#endif