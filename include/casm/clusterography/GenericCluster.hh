#ifndef CASM_GenericCluster
#define CASM_GenericCluster

#include <vector>

#include "casm/crystallography/UnitCellCoord.hh"

namespace casm {
  
  /** \defgroup Clusterography
      
      \brief Functions and classes related to clusters
  */
  
  /* -- GenericCluster Declarations ------------------------------------- */
  
  class Structure;
  class SymOp;
  
  /// \brief A cluster of UnitCellCoord
  ///
  /// \ingroup Clusterography
  ///
  template<typename CoordType>
  class GenericCluster {
    
    public:
    
    typedef unsigned int size_type;
    typedef std::vector<CoordType>::value_type value_type;
    typedef std::vector<CoordType>::iterator iterator;
    typedef std::vector<CoordType>::const_iterator const_iterator;
    
    typedef Structure PrimType;
    
    /// \brief Construct an empty UnitCellCoordCluster
    explicit GenericCluster(const PrimType& _prim) :
      m_prim_ptr(&_prim) {}
    
    /// \brief Construct a UnitCellCoordCluster with a range of UnitCellCoord
    template<typename InputIterator>
    GenericCluster(const PrimType& _prim,
                   InputIterator _uccoord_begin,
                   InputIterator _uccoord_end) :
      m_prim_ptr(&_prim),
      m_site(_uccoord_begin, _uccoord_end) {}
    
    /// \brief Return a reference to the primitive Structure
    const Structure& prim() const {
      return *m_prim_ptr;
    }
    
    /// \brief Iterator to first UnitCellCoord in the cluster
    iterator begin() {
      return m_site.begin();
    }
    
    /// \brief Iterator to first UnitCellCoord in the cluster
    const_iterator begin() const {
      return m_site.begin();
    }
    
    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    iterator end() {
      return m_site.end();
    }
    
    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    const_iterator end() const {
      return m_site.end();
    }
    
    /// \brief Iterator to first UnitCellCoord in the cluster
    const_iterator cbegin() const {
      return m_site.cbegin();
    }
    
    /// \brief Iterator to the past-the-last UnitCellCoord in the cluster
    const_iterator cend() const {
      return m_site.cend();
    }
    
    /// \brief Number of UnitCellCoords in the cluster
    size_type size() const {
      return m_site.size();
    }
    
    /// \brief Access a UnitCellCoord in the cluster by index
    value_type& operator[](size_type index) {
      return m_site[index];
    }
    
    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type& operator[](size_type index) const {
      return m_site[index];
    }
    
    /// \brief Access a UnitCellCoord in the cluster by index
    value_type& site(size_type index) {
      return m_site[index];
    }
    
    /// \brief Access a UnitCellCoord in the cluster by index
    const value_type& site(size_type index) const {
      return m_site[index];
    }
    
    /// \brief Access vector of sites
    std::vector<CoordType>& sites() {
      return m_site;
    }
    
    /// \brief const Access vector of sites
    const std::vector<CoordType>& sites() const {
      return m_site;
    }
    
    
    /// \brief Translate the cluster by a UnitCell translation
    GenericCluster<CoordType>& operator+=(UnitCell trans) {
      for(auto it = begin(); it != end(); ++it) {
        *it += trans;
      }
      return *this;
    }
    
    /// \brief Translate the UnitCellCoordCluster by a UnitCell translation
    GenericCluster<CoordType>& operator-=(UnitCell trans) {
      for(auto it = begin(); it != end(); ++it) {
        *it -= trans;
      }
      return *this;
    }
    
    
    GenericCluster<CoordType>& apply_sym(const SymOp& op) {
      for(auto it=begin(); it!=end(); ++it) {
        apply(op, *it);
      }
    }
    
    
    private:
    
    const PrimType* m_prim_ptr;
    
    std::vector<CoordType> m_site;
    
  };
  
  /// \brief Translate a cluster 
  template<typename CoordType>
  GenericCluster<CoordType> operator+(GenericCluster<CoordType> cluster, UnitCell trans);
  
  /// \brief Translate a cluster 
  ///
  /// \relates GenericCluster
  ///
  template<typename CoordType>
  GenericCluster<CoordType> operator-(GenericCluster<CoordType> cluster, UnitCell trans);
  
  
  /// \brief Default intra orbit comparison of clusters
  ///
  /// Implements:
  /// \code
  /// std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare);
  /// \endcode
  template<typename CoordType, typename Compare = std::less<CoordType> >
  bool cluster_intra_orbit_compare(const GenericCluster<CoordType>& A, 
                                   const GenericCluster<CoordType>& B, 
                                   Compare compare = Compare()) {
    std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare);
  }
  
  /// \brief Default comparison of orbit prototypes
  ///
  /// \param A, B clusters to check for A < B
  /// \param compare Compare concept functor for lexicographical comparison of cluster sites
  /// \param tol Tolerance for comparison of site-to-site distances
  ///
  /// Compares:
  /// - cluster size
  /// - site-to-site distances (max to min)
  /// - lexicographical_compare of sites, via compare
  ///
  /// Implements:
  /// \code
  /// // first compare cluster size
  /// if(A.size() != B.size()) {
  ///   return A.size() < B.size();
  /// }
  /// 
  /// // next compare invariants
  /// if( compare(ClusterInvariants(A), ClusterInvariants(B), tol) ) {
  ///   return true;
  /// }
  /// if( compare(ClusterInvariants(B), ClusterInvariants(A), tol) ) {
  ///   return false;
  /// }
  /// 
  /// // next lexicographical_compare of UnitCellCoord in A and B
  /// return std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare);
  /// \endcode
  template<typename CoordType, typename Compare = std::less<CoordType> >
  bool cluster_inter_orbit_compare(const GenericCluster<CoordType>& A, 
                                   const GenericCluster<CoordType>& B, 
                                   Compare compare = Compare(),
                                   double tol) {
    
    // first compare cluster size
    if(A.size() != B.size()) {
      return A.size() < B.size();
    }
    
    // next compare invariants
    if( compare(ClusterInvariants(A), ClusterInvariants(B), tol) ) {
      return true;
    }
    if( compare(ClusterInvariants(B), ClusterInvariants(A), tol) ) {
      return false;
    }
    
    // next lexicographical_compare of UnitCellCoord in A and B
    return std::lexicographical_compare(A.begin(), A.end(), B.begin(), B.end(), compare);
  }
  
  
  
  /* -- GenericCluster Definitions ------------------------------------- */
  
  /// \brief Translate a cluster 
  ///
  /// \relates GenericCluster
  ///
  template<typename CoordType>
  GenericCluster<CoordType> operator+(GenericCluster<CoordType> cluster, UnitCell trans) {
    return cluster += trans;
  }
  
  /// \brief Translate a cluster 
  ///
  /// \relates GenericCluster
  ///
  template<typename CoordType>
  GenericCluster<CoordType> operator-(GenericCluster<CoordType> cluster, UnitCell trans) {
    return cluster -= trans;
  }
  
}

#endif