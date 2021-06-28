#ifndef CASM_NeighborhoodInfo
#define CASM_NeighborhoodInfo

#include <map>
#include <set>
#include <vector>

#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/global/definitions.hh"

namespace CASM {

class PrimNeighborList;

// Useful neighborhood info for a set of cluster orbits
struct NeighborhoodInfo {
  /// Make NeighborhoodInfo
  template <typename ClusterOrbitIterator>
  NeighborhoodInfo(PrimNeighborList &prim_neighbor_list,
                   ClusterOrbitIterator begin, ClusterOrbitIterator end);

  /// \brief Size of neighbor list required to include all sites in
  /// unitcellcoord_neighborhood.
  Index neighbor_list_size;

  /// \brief UnitCellCoord for neighborhood of origin unit cell
  std::set<xtal::UnitCellCoord> unitcellcoord_neighborhood;

  /// \brief UnitCell for neighborhood of origin unit cell
  std::set<xtal::UnitCell> unitcell_neighborhood;

  /// \brief Number of distinct point correlations per unit cell
  Index n_point_corr;

  /// \brief UnitCellCoord of distinct point correlations in origin unit cell.
  /// Size should be equal to n_point_corr.
  std::vector<xtal::UnitCellCoord> point_corr_unitcellcoord;

  /// \brief Asymmetric unit indices of distinct point correlations. Size
  /// should be equal to n_point_corr.
  std::vector<Index> point_corr_asymmetric_unit_indices;

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

class ClusterSpecs;
class IntegralCluster;

std::unique_ptr<NeighborhoodInfo> make_neighborhood_info(
    ClusterSpecs const &cluster_specs, PrimNeighborList &prim_neighbor_list);

std::unique_ptr<NeighborhoodInfo> make_neighborhood_info(
    ClusterSpecs const &cluster_specs,
    std::vector<IntegralCluster> const &generating_elements,
    PrimNeighborList &prim_neighbor_list);

}  // namespace CASM

#endif
