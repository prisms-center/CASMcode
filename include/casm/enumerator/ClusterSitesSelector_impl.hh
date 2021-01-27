#ifndef CASM_enumerator_ClusterSitesSelector_impl
#define CASM_enumerator_ClusterSitesSelector_impl

#include "casm/clusterography/SupercellClusterOrbits_impl.hh"
#include "casm/enumerator/ClusterSitesSelector.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/symmetry/SubOrbits_impl.hh"

namespace CASM {

/// Create ConfigEnumInput with unique cluster sites selected, starting from
/// "prim_periodic" orbits
///
template <typename Inserter>
Inserter select_cluster_sites(
    ConfigEnumInput const &reference_config_enum_input,
    std::vector<PrimPeriodicIntegralClusterOrbit> const &orbits,
    Inserter result) {
  if (orbits.size() == 0) {
    return result;
  }

  Configuration const &ref_config = reference_config_enum_input.configuration();
  std::shared_ptr<Structure const> shared_prim =
      ref_config.supercell().shared_prim();
  SupercellSymInfo const &ref_supercell_sym_info =
      ref_config.supercell().sym_info();
  Configuration prim_ref_config = ref_config.primitive();
  SupercellSymInfo const &prim_ref_supercell_sym_info =
      prim_ref_config.supercell().sym_info();
  PrimPeriodicSymCompare<IntegralCluster> prim_periodic_sym_compare =
      orbits.begin()->sym_compare();

  // Construct orbit generating elements for orbits with
  // prim_ref_supercell_sym_info.factor_group() symmetry
  // - shared_prim->factor_group() -> prim_ref_supercell_sym_info.factor_group()
  // symmetry breaking
  std::vector<IntegralCluster> prim_ref_supercell_generators;
  make_suborbit_generators(
      shared_prim->factor_group().begin(), shared_prim->factor_group().end(),
      prim_ref_supercell_sym_info.factor_group().begin(),
      prim_ref_supercell_sym_info.factor_group().end(),
      prototype_iterator(orbits.begin()), prototype_iterator(orbits.end()),
      prim_periodic_sym_compare,
      std::back_inserter(prim_ref_supercell_generators));

  // Construct orbit generating elements for orbits with
  // prim_ref_config.factor_group() symmetry
  // - prim_ref_supercell_sym_info.factor_group() ->
  // prim_ref_config.factor_group() symmetry breaking
  ScelPeriodicSymCompare<IntegralCluster> prim_ref_supercell_sym_compare{
      shared_prim, prim_ref_supercell_sym_info.transformation_matrix_to_super(),
      shared_prim->lattice().tol()};
  std::vector<PermuteIterator> prim_ref_config_fg =
      prim_ref_config.factor_group();
  std::vector<IntegralCluster> prim_ref_config_orbit_generators;
  make_suborbit_generators(
      prim_ref_supercell_sym_info.permute_begin(),
      prim_ref_supercell_sym_info.permute_end(), prim_ref_config_fg.begin(),
      prim_ref_config_fg.end(), prim_ref_supercell_generators.begin(),
      prim_ref_supercell_generators.end(), prim_ref_supercell_sym_compare,
      std::back_inserter(prim_ref_config_orbit_generators));

  // Construct orbit generating elements for orbits with
  // ref_config.factor_group() symmetry under periodic boundary conditions
  // prim_ref_config.factor_group() -> ref_config.factor_group() w/ periodic
  // boundar conditions
  std::vector<PermuteIterator> ref_config_fg = ref_config.factor_group();
  std::vector<Permutation> ref_config_inverse_permutations =
      make_inverse_permutations(ref_config_fg.begin(), ref_config_fg.end());
  auto ref_config_orbit_generators_pbc =
      make_orbit_generators_under_periodic_boundary_conditions(
          ref_supercell_sym_info, ref_config_inverse_permutations.begin(),
          ref_config_inverse_permutations.end(),
          prim_ref_config_orbit_generators.begin(),
          prim_ref_config_orbit_generators.end());

  // Construct ConfigEnumInput with ref_config & cluster site indices
  for (auto const &cluster_site_indices : ref_config_orbit_generators_pbc) {
    *result++ = ConfigEnumInput{ref_config, cluster_site_indices};
  }
  return result;
}

}  // namespace CASM

#endif
