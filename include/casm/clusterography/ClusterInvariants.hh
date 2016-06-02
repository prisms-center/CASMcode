#ifndef CASM_ClusterInvariants
#define CASM_ClusterInvariants

#include <vector>
#include <iostream>

namespace CASM {

  /** \defgroup Clusterography

      \brief Functions and classes related to clusters
  */

  /* -- ClusterInvariants Declarations ------------------------------------- */

  /// \brief Stores cluster invariants: number of sites and site distances
  ///
  /// Default expects (as for CoordCluster):
  /// - \code ClusterType::size() \endcode
  /// - \code ClusterType::coordinate(size_type index) \endcode
  ///
  /// \ingroup Clusterography
  ///
  template<typename ClusterType>
  class ClusterInvariants {

  public:

    /// \brief Construct and calculate cluster invariants
    explicit ClusterInvariants(const ClusterType &cluster) {

      // save size of cluster
      m_size = cluster.size();

      // calculate distances between points
      std::vector<double> disp;
      for(int i = 0; i < m_size; i++) {
        for(int j = i + 1; j < m_size; j++) {
          m_disp.push_back((cluster.coordinate(i) - cluster.coordinate(j)).const_cart().norm());
        }
      }
      std::sort(m_disp.begin(), m_disp.end());
    }

    /// \brief Number of elements in the cluster
    int size() const {
      return m_size;
    }

    /// \brief const Access displacements between coordinates in the cluster, sorted in ascending order
    const std::vector<double> &displacement() const {
      return m_disp;
    }

    std::unique_ptr<ClusterInvariants<ClusterType> > clone() const {
      return notstd::make_unique<ClusterInvariants<ClusterType> >(*this);
    }

  private:

    /// \brief Number of UnitCellCoords in cluster
    int m_size;

    /// \brief Displacement between each pair of UnitCellCoords, sorted in ascending order
    std::vector<double> m_disp;

  };

  /// \brief Check if ClusterInvariants are equal
  template<typename ClusterType>
  bool almost_equal(const ClusterInvariants<ClusterType> &A, const ClusterInvariants<ClusterType> &B, double tol);

  /// \brief Compare ClusterInvariants
  template<typename ClusterType>
  bool compare(const ClusterInvariants<ClusterType> &A, const ClusterInvariants<ClusterType> &B, double tol);

  /// \brief Print ClusterInvariants
  template<typename ClusterType>
  std::ostream &operator<<(std::ostream &sout, const ClusterInvariants<ClusterType> &invariants);

}

#endif