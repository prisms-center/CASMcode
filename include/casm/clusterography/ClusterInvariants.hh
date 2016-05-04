#ifndef CASM_ClusterInvariants
#define CASM_ClusterInvariants

#include <vector>

namespace casm {
  
  /** \defgroup Clusterography
      
      \brief Functions and classes related to clusters
  */
  
  /* -- ClusterInvariants Declarations ------------------------------------- */
  
  /// \brief Stores cluster invariants: number of sites and site distances 
  ///
  /// Expects:
  /// - \code ClusterType::size() \endcode 
  /// - \code ClusterType::coordinate(size_type index) \endcode
  ///
  /// \ingroup Clusterography
  ///
  class ClusterInvariants {
    
    public:
    
    /// \brief Construct and calculate cluster invariants
    template<typename ClusterType>
    explicit ClusterInvariants(const ClusterType &cluster) {
      
      // save size of cluster
      m_size = cluster.size();
        
      // calculate distances between points
      std::vector<double> disp;
      for(int i=0; i<m_size; i++) {
        for(int j=i+1; j<m_size; j++) {
          m_disp.push_back((cluster.coordinate(i) - cluster.coordinate(j)).norm());
        }
      }
      std::sort(m_disp.begin(), m_disp.end());
    }
    
    /// \brief Number of UnitCellCoords in the cluster
    int size() const {
      return m_size;
    }
    
    /// \brief const Access displacements between each UnitCellCoord in the cluster, sorted in ascending order
    const std::vector<double>& displacement() const {
      return m_disp;
    }
    
    
    
    private:
    
    /// \brief Number of UnitCellCoords in cluster
    int m_size;
    
    /// \brief Displacement between each pair of UnitCellCoords, sorted in ascending order
    std::vector<double> m_disp;
    
  };
  
  /// \brief Check if ClusterInvariants are equal
  bool almost_equal(const ClusterInvariants &A, const ClusterInvariants &B, double tol);
  
  /// \brief Compare ClusterInvariants
  bool compare(const ClusterInvariants &A, const ClusterInvariants &B, double tol);
  
  /// \brief Print ClusterInvariants
  std::ostream& operator<<(std::ostream &sout, const ClusterInvariants &invariants);
  
}

#endif