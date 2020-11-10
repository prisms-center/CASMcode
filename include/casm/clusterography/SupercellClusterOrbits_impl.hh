#ifndef CASM_SupercellClusterOrbits_impl
#define CASM_SupercellClusterOrbits_impl

#include "casm/clusterography/SupercellClusterOrbits.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

  /// Make inverse permutations
  ///
  /// Inverse permutations allow fast permutation of cluster site indices. For the details see
  /// `permute_cluster_site_indices`.
  template<typename PermuteIteratorIt>
  std::vector<Permutation> make_inverse_permutations(
    PermuteIteratorIt begin,
    PermuteIteratorIt end) {
    std::vector<Permutation> inverse_permutations;
    for(auto permute_it = begin; permute_it != end; ++permute_it) {
      inverse_permutations.push_back(permute_it->combined_permute().inverse());
    }
    return inverse_permutations;
  }

  /// Apply inverse permutations and return the maximum equivalent cluster_site_indices
  template<typename InversePermutationIterator>
  std::set<Index> make_canonical_cluster_site_indices(
    InversePermutationIterator inverse_permutations_begin,
    InversePermutationIterator inverse_permutations_end,
    std::set<Index> const &cluster_site_indices) {
    auto max = cluster_site_indices;
    auto ipermute_it = inverse_permutations_begin;
    for(; ipermute_it != inverse_permutations_end; ++ipermute_it) {
      auto test = permute_cluster_site_indices(*ipermute_it, cluster_site_indices);
      if(test > max) {
        max = test;
      }
    }
    return max;
  }

  /// Return orbit generators, as sets of cluster site indices, under periodic boundary conditions
  ///
  /// Given a range of elements (IntegralCluster), convert to sets of cluster site indices,
  /// and apply permute group operations (using inverse permutations to permute site indices
  /// instead of site values) to identify unique cluster site indices under periodic boundary
  /// conditions
  ///
  /// Note:
  /// - A supercell/configuration may have lower symmetry than the prim, resulting in orbit
  ///   breaking. The input elements must have been generated using the supercell/configuration
  ///   symmetry in order to generate all the unique orbit generating elements under periodic
  ///   boundary conditions. Using "prim_periodic" orbit prototypes is allowed, but then the output
  ///   of this method will not reflect any symmetry breaking do to the supercell/configuration.
  ///
  /// \param sym_info SupercellSymInfo for the supercell in which cluster orbits are generated
  /// \param inverse_permutations_begin, inverse_permutations_end The inverse permutations of the
  ///        permute group that will be used to generate cluster orbits in the supercell (i.e
  ///        the result of `make_inverse_permutations` for the range [sym_info.permute_begin(),
  ///        sym_info.permute_end()) or a subgroup (i.e. `Configuration::factor_group()`).
  /// \param element_begin, element_end Iterators that provide orbit elements (IntegralCluster).
  ///
  template<typename InversePermutationIterator, typename ElementIterator>
  std::set<std::set<Index>, ClusterSiteIndicesCompare> make_orbit_generators_under_periodic_boundary_conditions(
    SupercellSymInfo const &sym_info,
    InversePermutationIterator inverse_permutations_begin,
    InversePermutationIterator inverse_permutations_end,
    ElementIterator element_begin,
    ElementIterator element_end) {

    // Using set::set<Index> for cluster site indices:
    // - removes duplicate sites / clusters that may occur due to periodic boundary conditions
    // - site indices transform via inverse permutations
    // - sorts cluster site indices for easy comparison to find the canonical cluster site indices
    std::set<std::set<Index>, ClusterSiteIndicesCompare> orbit_generators;
    for(auto it = element_begin; it != element_end; ++it) {
      auto element = make_cluster_site_indices(*it, sym_info);
      auto canonical_element = make_canonical_cluster_site_indices(
                                 inverse_permutations_begin,
                                 inverse_permutations_end,
                                 element);
      orbit_generators.insert(canonical_element);
    }

    return orbit_generators;
  }

}

#endif
