#ifndef CASM_GenericCluster
#define CASM_GenericCluster

#include <vector>

namespace CASM {

  template <typename ClusterType> ClusterInvariants;

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /* -- GenericCluster Declarations ------------------------------------- */

  /// \brief A cluster of anything
  ///
  /// - Needs a valid ClusterInvariants<GenericCluster<Element> > to help comparisons
  ///
  /// \ingroup Clusterography
  ///
  template<typename Element>
  class GenericCluster {

  public:

    typedef unsigned int size_type;
    typedef std::vector<Element>::value_type value_type;
    typedef std::vector<Element>::iterator iterator;
    typedef std::vector<Element>::const_iterator const_iterator;
    typedef ClusterInvariants<GenericCluster<Element> > InvariantsType;

    /// \brief Construct an empty GenericCluster
    GenericCluster() {}

    /// \brief Construct a GenericCluster with a range of Element
    template<typename InputIterator>
    GenericCluster(InputIterator _begin,
                   InputIterator _end) :
      m_element(_uccoord_begin, _uccoord_end) {}

    /// \brief Default move constructor
    GenericCluster(GenericCluster &&other) = default;

    /// \brief Default move assignment
    GenericCluster &operator=(GenericCluster && other) = default;


    /// \brief Return a reference to the primitive Structure
    const Structure &prim() const {
      return *m_prim_ptr;
    }

    /// \brief Iterator to first UnitCellCoord in the cluster
    iterator begin() {
      m_invarients.reset();
      return m_element.begin();
    }

    /// \brief Iterator to first UnitCellCoord in the cluster
    const_iterator begin() const {
      return m_element.begin();
    }

    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    iterator end() {
      m_invarients.reset();
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
      m_invarients.reset();
      return m_element[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type &operator[](size_type index) const {
      return m_element[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    value_type &element(size_type index) {
      m_invarients.reset();
      return m_element[index];
    }

    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type &element(size_type index) const {
      return m_element[index];
    }

    /// \brief Access vector of elements
    std::vector<Element> &elements() {
      m_invarients.reset();
      return m_element;
    }

    /// \brief const Access vector of elements
    const std::vector<Element> &elements() const {
      return m_element;
    }

    const InvariantsType &invariants() const {
      if(!m_invariants) {
        m_invariants = make_cloneable<InvariantsType>(*this);
      }
      return *m_invariants;
    }

    GenericCluster<CoordType> &apply_sym(const SymOp &op) {
      for(auto it = begin(); it != end(); ++it) {
        apply(op, *it);
      }
    }


  private:

    std::vector<Element> m_element;

    notstd::cloneable_ptr<InvariantsType> m_invariants;

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
  template<typename Element, typename Compare = std::less<Element> >
  bool cluster_intra_orbit_compare(const GenericCluster<Element> &A,
                                   const GenericCluster<Element> &B,
                                   Compare compare = Compare()) {
    std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare);
  }

  /// \brief Default comparison of orbit prototypes
  ///
  /// \param A, B clusters to check for A < B
  /// \param tol Tolerance for comparison of element-to-element distances
  /// \param compare Compare concept functor for lexicographical comparison of cluster elements
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
  /// if( compare(ClusterInvariants<Element>(A), ClusterInvariants<Element>(B), tol) ) {
  ///   return true;
  /// }
  /// if( compare(ClusterInvariants<Element>(B), ClusterInvariants<Element>(A), tol) ) {
  ///   return false;
  /// }
  ///
  /// // next lexicographical_compare of UnitCellCoord in A and B
  /// return std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare);
  /// \endcode
  ///
  /// \ingroup Clusterography
  ///
  template<typename Element, typename Compare = std::less<Element> >
  bool cluster_inter_orbit_compare(const GenericCluster<Element> &A,
                                   const GenericCluster<Element> &B,
                                   double tol,
                                   Compare compare = Compare()) {

    // first compare cluster size
    if(A.size() != B.size()) {
      return A.size() < B.size();
    }

    // next compare invariants
    if(compare(A.invariants(), B.invariants(), tol)) {
      return true;
    }
    if(compare(B.invariants(), A.invariants(), tol)) {
      return false;
    }

    // next lexicographical_compare of Element in A and B
    return std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare);
  }

}

#endif