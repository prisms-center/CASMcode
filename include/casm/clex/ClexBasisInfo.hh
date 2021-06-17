#ifndef CASM_ClexBasisInfo
#define CASM_ClexBasisInfo

#include <map>
#include <set>

#include "casm/clusterography/IntegralCluster_impl.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/global/definitions.hh"

namespace CASM {

struct ClexBasisSpecs;
class PrimNeighborList;
class Structure;

/// Clex basis function info
///
/// Stored for each cluster expansion basis function in ClexBasisInfo
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

/// Useful ClexBasis info
///
/// - Construct with `make_clex_basis_info`
struct ClexBasisInfo {
  /// \brief Information on basis functions
  std::vector<ClexBasisFunctionInfo> basis_function_info;

  /// \brief Maximum possible number of sites involved in evaluation
  Index neighborhood_size;

  /// \brief Number of point correlations
  Index n_point_corr;

  /// \brief Point correlation neighborhoods
  ///
  /// A map of UnitCellCoord (translated as necessary to the canonical
  /// unit cell by the SymCompare method), to the set of UnitCellCoord that it
  /// is in clusters with. For cluster function evaluation, this gives the
  /// neighborhood of sites whose DoF values are needed to evaluate clusters
  /// involving a particular site.
  std::map<xtal::UnitCellCoord, std::set<xtal::UnitCellCoord>>
      site_dependency_neighborhoods;
};

/// Make ClexBasisInfo
ClexBasisInfo make_clex_basis_info(
    std::shared_ptr<Structure const> const &shared_prim,
    ClexBasisSpecs const &basis_set_specs,
    PrimNeighborList &prim_neighbor_list);

}  // namespace CASM

#endif
