#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/SupercellClusterOrbits.hh"
#include "casm/container/Permutation.hh"
#include "casm/symmetry/PermuteIterator.hh"
#include "casm/symmetry/SupercellSymInfo.hh"

namespace CASM {

  // Note:
  //     On their own `make_cluster_site_indices` and `make_cluster` only require
  //     `xtal::UnitCellCoordIndexConverter`, and not `SupercellSymInfo`. However, in this context
  //     they are always going to be used with `SupercellSymInfo` present to provide the permute
  //     group.

  /// Return site indices of cluster sites
  std::set<Index> make_cluster_site_indices(IntegralCluster const &cluster,
                                            SupercellSymInfo const &sym_info) {
    std::set<Index> cluster_site_indices;
    for(xtal::UnitCellCoord const &uccoord : cluster) {
      cluster_site_indices.insert(sym_info.unitcellcoord_index_converter()(uccoord));
    }
    return cluster_site_indices;
  }

  /// Return cluster from cluster site indices
  IntegralCluster make_cluster(std::set<Index> const &cluster_site_indices,
                               Structure const &prim,
                               SupercellSymInfo const &sym_info) {
    IntegralCluster cluster {prim};
    for(Index site_index : cluster_site_indices) {
      cluster.elements().push_back(sym_info.unitcellcoord_index_converter()(site_index));
    }
    return cluster;
  }

  /// Rather than permute values, here we want to permute cluster site indices
  ///
  /// Since, when permuting values, value[i] = value[permute[i]], or,
  ///    site_index_before = permute[site_index_after],
  /// to permute site indices use inverse_permute = permute.inverse, and then
  ///    site_index_after = inverse_permute[site_index_before]
  ///
  std::set<Index> permute_cluster_site_indices(Permutation const &inverse_permute,
                                               std::set<Index> const &cluster_site_indices) {
    std::set<Index> cluster_site_indices_after;
    for(Index site_index : cluster_site_indices) {
      cluster_site_indices_after.insert(inverse_permute[site_index]);
    }
    return cluster_site_indices_after;
  }

  /// Sort std::set<Index> by size, then value (lexicographical compare)
  bool ClusterSiteIndicesCompare::operator()(std::set<Index> const &LHS, std::set<Index> const &RHS) const {
    if(LHS.size() < RHS.size()) {
      return true;
    }
    if(LHS.size() > RHS.size()) {
      return false;
    }
    return LHS < RHS;
  }

}
