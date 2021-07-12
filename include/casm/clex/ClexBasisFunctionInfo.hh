#ifndef CASM_ClexBasisFunctionInfo
#define CASM_ClexBasisFunctionInfo

#include "casm/clusterography/IntegralCluster.hh"

namespace CASM {

class ClexBasis;

/// Clex basis function info
///
/// Stores info for each cluster expansion basis function
struct ClexBasisFunctionInfo {
  ClexBasisFunctionInfo(Index _linear_orbit_index, Index _linear_function_index,
                        IntegralCluster _prototype, Index _multiplicity,
                        std::vector<Index> _cluster_invariant_group_indices)
      : linear_orbit_index(_linear_orbit_index),
        linear_function_index(_linear_function_index),
        prototype(_prototype),
        multiplicity(_multiplicity),
        cluster_invariant_group_indices(_cluster_invariant_group_indices) {}

  /// \brief Linear index of the cluster orbit
  Index linear_orbit_index;

  /// \brief Linear index of the cluster basis function
  Index linear_function_index;

  /// \brief Prototype cluster
  IntegralCluster prototype;

  /// \brief Number of equivalent clusters
  Index multiplicity;

  /// \brief Indices of prim factor group operations that map this cluster onto
  /// itself according to the symmetry used to generate clusters
  std::vector<Index> cluster_invariant_group_indices;
};

/// Make ClexBasisFunctionInfo
template <typename ClusterOrbitIterator>
std::vector<ClexBasisFunctionInfo> make_clex_basis_function_info(
    ClexBasis const &clex_basis, ClusterOrbitIterator begin,
    ClusterOrbitIterator end);

}  // namespace CASM

#endif
