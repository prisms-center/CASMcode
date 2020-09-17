#ifndef CASM_SupercellClusterOrbits
#define CASM_SupercellClusterOrbits

#include <set>
#include <vector>

#include "casm/global/definitions.hh"

namespace CASM {

  class IntegralCluster;
  class Structure;
  class SupercellSymInfo;
  class Permutation;
  class PermuteIterator;


  /// Return site indices of cluster sites
  std::set<Index> make_cluster_site_indices(IntegralCluster const &cluster,
                                            SupercellSymInfo const &sym_info);

  /// Return cluster from cluster site indices
  IntegralCluster make_cluster(std::set<Index> const &cluster_site_indices,
                               Structure const &prim,
                               SupercellSymInfo const &sym_info);

  /// Return true if the permutation does not mix cluster sites and other sites
  bool cluster_site_indices_are_invariant(PermuteIterator const &permute_it,
                                          std::set<Index> const &cluster_site_indices);

  /// Make inverse permutations
  template<typename PermuteIteratorIt>
  std::vector<Permutation> make_inverse_permutations(
    PermuteIteratorIt begin,
    PermuteIteratorIt end);

  /// Rather than permute values, here we want to permute cluster site indices
  std::set<Index> permute_cluster_site_indices(Permutation const &inverse_permute,
                                               std::set<Index> const &cluster_site_indices);

  /// Apply inverse permutations and return the maximum equivalent cluster_site_indices
  template<typename InversePermutationIterator>
  std::set<Index> make_canonical_cluster_site_indices(
    InversePermutationIterator inverse_permutations_begin,
    InversePermutationIterator inverse_permutations_end,
    std::set<Index> const &cluster_site_indices);

  /// Sort std::set<Index> by size, then value (lexicographical compare)
  struct ClusterSiteIndicesCompare {
    bool operator()(std::set<Index> const &LHS, std::set<Index> const &RHS) const;
  };

  /// Return "within_scel" orbit generators, as sets of cluster site indices
  template<typename InversePermutationIterator, typename ElementIterator>
  std::set<std::set<Index>, ClusterSiteIndicesCompare> make_orbit_generators_under_periodic_boundary_conditions(
    SupercellSymInfo const &sym_info,
    InversePermutationIterator inverse_permutations_begin,
    InversePermutationIterator inverse_permutations_end,
    ElementIterator element_begin,
    ElementIterator element_end);

}

#endif
