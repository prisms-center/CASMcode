#ifndef CASM_ClusterInvariants
#define CASM_ClusterInvariants

#include <iostream>
#include <vector>

#include "casm/global/eigen.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

class IntegralCluster;

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
class ClusterInvariants : public notstd::Cloneable {
  CLONEABLE(ClusterInvariants)
 public:
  /// \brief Construct and calculate cluster invariants
  explicit ClusterInvariants(IntegralCluster const &cluster);

  /// \brief Number of elements in the cluster
  int size() const;

  /// \brief const Access displacements between coordinates in the cluster,
  /// sorted in ascending order
  std::vector<double> const &displacement() const;

 private:
  /// \brief Number of UnitCellCoords in cluster
  int m_size;

  /// \brief Displacement between each pair of UnitCellCoords, sorted in
  /// ascending order
  std::vector<double> m_disp;
};

/// \brief Check if ClusterInvariants are equal
bool almost_equal(ClusterInvariants const &A, ClusterInvariants const &B,
                  double tol);

/// \brief Compare ClusterInvariants
bool compare(ClusterInvariants const &A, ClusterInvariants const &B,
             double tol);

/// \brief Stores cluster invariants: number of sites and site distances (using
/// robust_min_dist)
///
/// Default expects (as for CoordCluster):
/// - \code ClusterType::size() \endcode
/// - \code ClusterType::coordinate(size_type index) \endcode
///
/// \ingroup Clusterography
///
class WithinScelClusterInvariants : public notstd::Cloneable {
  CLONEABLE(WithinScelClusterInvariants)
 public:
  /// \brief Construct and calculate cluster invariants, using robust_min_dist
  /// in the supercell lattice
  explicit WithinScelClusterInvariants(IntegralCluster const &cluster,
                                       Eigen::Matrix3l const &transf_mat);

  /// \brief Number of elements in the cluster
  int size() const;

  /// \brief const Access displacements between coordinates in the cluster,
  /// sorted in ascending order, calculated using robust_min_dist in the
  /// supercell lattice
  std::vector<double> const &displacement() const;

 private:
  /// \brief Number of UnitCellCoords in cluster
  int m_size;

  /// \brief Displacement between each pair of UnitCellCoords, sorted in
  /// ascending order
  std::vector<double> m_disp;
};

/// \brief Check if WithinScelClusterInvariants are equal
bool almost_equal(WithinScelClusterInvariants const &A,
                  WithinScelClusterInvariants const &B, double tol);

/// \brief Compare WithinScelClusterInvariants
bool compare(WithinScelClusterInvariants const &A,
             WithinScelClusterInvariants const &B, double tol);

}  // namespace CASM

#endif
