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
    Configuration const &reference_config,
    std::vector<PrimPeriodicIntegralClusterOrbit> const &orbits,
    Inserter result) {
  if (orbits.size() == 0) {
    return result;
  }

  Configuration const &ref_config = reference_config;
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

namespace {
bool _is_group(std::vector<PermuteIterator> const &permute_group) {
  std::set<PermuteIterator> check;
  for (auto const &it_a : permute_group) {
    for (auto const &it_b : permute_group) {
      check.insert(it_a * it_b);
    }
  }
  return check.size() == permute_group.size();
}

/// \brief Return all configurations that are equivalent to reference_config
///     (as infinite crystals), but inequivalent under permute_group
///
/// Usage: use this to generate the background configurations for distinct
/// combinations of background + phenomenal cluster, where permute_group is
/// the invariant subgroup of the phenomenal cluster consistent with the
/// reference_configuration's supercell.
///
/// Note: this does not consider distinct phenomenal clusters in the supercell.
///
/// \param reference_config A reference configuration
/// \param permute_group A permute group in the supercell of `reference_config`.
///
/// \returns All configurations that are equivalent to reference_config
///     (generated using `make_all_super_configurations` so they are
///     equivalent as infinite crystals under the prim factor group), but
///     inequivalent under the action of permute_group.
std::set<Configuration> _make_inequivalent_reference_configurations(
    Configuration const &reference_config,
    std::vector<PermuteIterator> const &permute_group) {
  std::shared_ptr<Supercell const> shared_supercell =
      reference_config.shared_supercell();
  if (shared_supercell == nullptr) {
    shared_supercell = std::make_shared<Supercell>(
        reference_config.supercell().shared_prim(),
        reference_config.supercell().sym_info().supercell_lattice());
  }
  std::vector<Configuration> all_configs;
  make_all_super_configurations(reference_config, shared_supercell,
                                std::back_inserter(all_configs));
  std::set<Configuration> inequivalent_configs;
  auto begin = permute_group.begin();
  auto end = permute_group.end();
  for (auto const &config : all_configs) {
    inequivalent_configs.insert(config.canonical_form(begin, end));
  }
  return inequivalent_configs;
}

template <typename Inserter>
Inserter _select_local_cluster_sites_fixed_background(
    Configuration const &reference_config,
    std::vector<LocalIntegralClusterOrbit> const &orbits,
    std::vector<PermuteIterator> const &generating_permute_group_in_supercell,
    Inserter result) {
  if (orbits.size() == 0) {
    return result;
  }

  Configuration const &ref_config = reference_config;
  SupercellSymInfo const &ref_supercell_sym_info =
      ref_config.supercell().sym_info();

  // Make generating permute group consistent w/ configuration
  Configuration prim_ref_config = ref_config.primitive();
  std::vector<PermuteIterator> prim_ref_config_fg =
      prim_ref_config.factor_group();
  std::set<PermuteIterator> prim_ref_config_fg_set{prim_ref_config_fg.begin(),
                                                   prim_ref_config_fg.end()};
  std::vector<PermuteIterator> generating_permute_group_in_config =
      make_allowed_permute(generating_permute_group_in_supercell.begin(),
                           generating_permute_group_in_supercell.end(),
                           prim_ref_config_fg_set);

  // Convert generating permute group in config to SymGroup
  SymGroup generating_group_in_config =
      make_sym_group(generating_permute_group_in_config.begin(),
                     generating_permute_group_in_config.end(),
                     ref_supercell_sym_info.prim_lattice());

  // Get orbit generators after symmetry breaking due to:
  // - local generating group in prim -> local generating group in config
  std::vector<IntegralCluster> orbit_generators_in_config;
  make_suborbit_generators(
      orbits[0].generating_group().begin(), orbits[0].generating_group().end(),
      generating_group_in_config.begin(), generating_group_in_config.end(),
      prototype_iterator(orbits.begin()), prototype_iterator(orbits.end()),
      orbits[0].sym_compare(), std::back_inserter(orbit_generators_in_config));

  // Construct orbit generating elements for orbits
  // under periodic boundary conditions
  std::vector<Permutation> inverse_permutations =
      make_inverse_permutations(generating_permute_group_in_config.begin(),
                                generating_permute_group_in_config.end());
  //
  auto orbit_generators_pbc =
      make_orbit_generators_under_periodic_boundary_conditions(
          ref_supercell_sym_info, inverse_permutations.begin(),
          inverse_permutations.end(), orbit_generators_in_config.begin(),
          orbit_generators_in_config.end());

  // Construct ConfigEnumInput with ref_config & cluster site indices
  for (auto const &cluster_site_indices : orbit_generators_pbc) {
    *result++ = ConfigEnumInput{ref_config, cluster_site_indices};
  }

  return result;
}

}  // namespace

/// Create ConfigEnumInput with unique cluster sites selected, starting from
/// local orbits
template <typename Inserter>
Inserter select_local_cluster_sites(
    Configuration const &reference_config,
    std::vector<LocalIntegralClusterOrbit> const &orbits,
    bool include_all_inequivalent_reference_configs, Inserter result) {
  if (orbits.size() == 0) {
    return result;
  }

  // Get objects to be used below
  SymGroup const &generating_group = orbits[0].generating_group();
  SymGroup const &prim_fg =
      reference_config.supercell().shared_prim()->factor_group();
  SupercellSymInfo const &ref_supercell_sym_info =
      reference_config.supercell().sym_info();

  // Make generating permute group consistent w/ supercell
  std::vector<PermuteIterator> generating_permute_group_in_supercell =
      make_local_permute_group(generating_group, prim_fg,
                               ref_supercell_sym_info);

  if (include_all_inequivalent_reference_configs) {
    std::set<Configuration> backgrounds =
        _make_inequivalent_reference_configurations(
            reference_config, generating_permute_group_in_supercell);
    for (auto const &background : backgrounds) {
      result = _select_local_cluster_sites_fixed_background(
          background, orbits, generating_permute_group_in_supercell, result);
    }
  } else {
    result = _select_local_cluster_sites_fixed_background(
        reference_config, orbits, generating_permute_group_in_supercell,
        result);
  }

  return result;
}

}  // namespace CASM

#endif
