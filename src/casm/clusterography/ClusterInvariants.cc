#include "casm/clusterography/ClusterInvariants.hh"

#include "casm/clusterography/IntegralCluster.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {

/// \brief Construct and calculate cluster invariants
ClusterInvariants::ClusterInvariants(IntegralCluster const &cluster) {
  // save size of cluster
  m_size = cluster.size();

  // calculate distances between points
  std::vector<double> disp;
  for (int i = 0; i < m_size; i++) {
    for (int j = i + 1; j < m_size; j++) {
      m_disp.push_back(
          (cluster.coordinate(i) - cluster.coordinate(j)).const_cart().norm());
    }
  }
  std::sort(m_disp.begin(), m_disp.end());
}

/// \brief Number of elements in the cluster
int ClusterInvariants::size() const { return m_size; }

/// \brief const Access displacements between coordinates in the cluster, sorted
/// in ascending order
const std::vector<double> &ClusterInvariants::displacement() const {
  return m_disp;
}

/// \brief Check if ClusterInvariants are equal
bool almost_equal(ClusterInvariants const &A, ClusterInvariants const &B,
                  double tol) {
  return A.size() == B.size() &&
         std::equal(A.displacement().cbegin(), A.displacement().cend(),
                    B.displacement().cbegin(), [&](double a, double b) {
                      return almost_equal(a, b, tol);
                    });
}

/// \brief Compare ClusterInvariants
///
/// \returns True if A < B, to specified tolerance
///
/// - First compares by number of sites in cluster
/// - Then compare all displacements, from longest to shortest
///
bool compare(ClusterInvariants const &A, ClusterInvariants const &B,
             double tol) {
  // first sort by number of sites in cluster
  if (A.size() < B.size()) {
    return true;
  }
  if (A.size() > B.size()) {
    return false;
  }

  // all displacements
  for (int i = A.displacement().size() - 1; i >= 0; i--) {
    if (almost_equal(A.displacement()[i], B.displacement()[i], tol)) {
      continue;
    }
    if (A.displacement()[i] < B.displacement()[i]) {
      return true;
    }
    if (A.displacement()[i] > B.displacement()[i]) {
      return false;
    }
  }
  return false;
}

/// \brief Construct and calculate cluster invariants, using robust_min_dist in
/// the supercell lattice
WithinScelClusterInvariants::WithinScelClusterInvariants(
    IntegralCluster const &cluster, Eigen::Matrix3l const &transf_mat) {
  /// \brief Supercell lattice to use for periodic min dist
  xtal::Lattice scel_lat =
      make_superlattice(cluster.prim().lattice(), transf_mat);

  // save size of cluster
  m_size = cluster.size();

  // calculate distances between points
  std::vector<xtal::Coordinate> coords;
  for (int i = 0; i < m_size; ++i) {
    coords.push_back(cluster.coordinate(i));
    coords.back().set_lattice(scel_lat, CART);
  }
  std::vector<double> disp;
  for (int i = 0; i < m_size; i++) {
    for (int j = i + 1; j < m_size; j++) {
      m_disp.push_back(coords[i].robust_min_dist(coords[j]));
    }
  }
  std::sort(m_disp.begin(), m_disp.end());
}

/// \brief Number of elements in the cluster
int WithinScelClusterInvariants::size() const { return m_size; }

/// \brief const Access displacements between coordinates in the cluster, sorted
/// in ascending order
std::vector<double> const &WithinScelClusterInvariants::displacement() const {
  return m_disp;
}

/// \brief Check if ClusterInvariants are equal
bool almost_equal(WithinScelClusterInvariants const &A,
                  WithinScelClusterInvariants const &B, double tol) {
  return A.size() == B.size() &&
         std::equal(A.displacement().cbegin(), A.displacement().cend(),
                    B.displacement().cbegin(), [&](double a, double b) {
                      return almost_equal(a, b, tol);
                    });
}

/// \brief Compare WithinScelClusterInvariants
///
/// \returns True if A < B, to specified tolerance
///
/// - First compares by number of sites in cluster
/// - Then compare all displacements, from longest to shortest
///
bool compare(WithinScelClusterInvariants const &A,
             WithinScelClusterInvariants const &B, double tol) {
  // first sort by number of sites in cluster
  if (A.size() < B.size()) {
    return true;
  }
  if (A.size() > B.size()) {
    return false;
  }

  // all displacements
  for (int i = A.displacement().size() - 1; i >= 0; i--) {
    if (almost_equal(A.displacement()[i], B.displacement()[i], tol)) {
      continue;
    }
    if (A.displacement()[i] < B.displacement()[i]) {
      return true;
    }
    if (A.displacement()[i] > B.displacement()[i]) {
      return false;
    }
  }
  return false;
}

}  // namespace CASM
