#ifndef CASM_enumerator_ClusterSitesSelector_impl
#define CASM_enumerator_ClusterSitesSelector_impl

#include "casm/clex/Supercell.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/enumerator/ClusterSitesSelector.hh"

namespace CASM {

  /// Create ConfigEnumInput with unique cluster sites selected, not implemented in general
  ///
  /// Note:
  /// - ClusterSitesSelector is only implemented for "within_scel" and "prim_periodic" orbits. For other
  ///   orbit types this will throw.
  /// - TODO: Implement for "local" orbits
  template<typename OrbitType, typename Inserter>
  Inserter select_cluster_sites(
    ConfigEnumInput const &reference_config_enum_input,
    std::vector<OrbitType> const &orbits,
    Inserter result) {
    std::stringstream msg;
    msg << "Error constructing ConfigEnumInput: Unsupported cluster_specs type. Please use one of: "
        << "\"periodic_max_length\" or \"within_scel_max_length\".";
    throw libcasm_runtime_error(msg.str());
  }

  /// Create ConfigEnumInput with unique cluster sites selected, starting from "within_scel" orbits
  ///
  /// For each orbit:
  /// - Copy reference_config_enum_input (keeping any sites already selected)
  /// - Set all orbit prototype sites on the copy
  /// - Insert the copy into the container pointed to by the `result` insert iterator
  template<typename Inserter>
  Inserter select_cluster_sites(
    ConfigEnumInput const &reference_config_enum_input,
    std::vector<WithinScelIntegralClusterOrbit> const &orbits,
    Inserter result) {
    for(auto const &orbit : orbits) {
      ConfigEnumInput tmp {reference_config_enum_input};
      tmp.select_sites(orbit.prototype());
      *result++ = tmp;
    }
    return result;
  }

  /// Create ConfigEnumInput with unique cluster sites selected, starting from "prim_periodic" orbits
  ///
  /// Method:
  /// - Converts "prim_periodic" orbits to "within_scel" using
  ///   `make_within_scel_orbits_from_prim_periodic` with the reference_config_enum_input supercell.
  /// - Then calls `select_cluster_sites` with "within_scel" orbits
  template<typename Inserter>
  Inserter select_cluster_sites(
    ConfigEnumInput const &reference_config_enum_input,
    std::vector<PrimPeriodicIntegralClusterOrbit> const &orbits,
    Inserter result) {

    // If we have PrimPeriodicIntegralClusterOrbit, then we need to translate all clusters
    // throughout the supercell and apply symmetry to find the orbits in the supercell.
    std::vector<WithinScelIntegralClusterOrbit> within_scel_orbits;
    Supercell const &reference_supercell = reference_config_enum_input.configuration().supercell();

    within_scel_orbits = make_within_scel_orbits_from_prim_periodic(
                           reference_supercell.shared_prim(),
                           reference_supercell.sym_info().transformation_matrix_to_super(),
                           make_invariant_group(reference_config_enum_input),
                           orbits);

    // Now use these to construct ConfigEnumInput with appropriate sites set
    return select_cluster_sites(reference_config_enum_input, within_scel_orbits, result);
  }


  template<typename Inserter>
  ClusterSitesSelector<Inserter>::ClusterSitesSelector(
    ConfigEnumInput const &_reference_config_enum_input,
    Inserter _inserter):
    reference_config_enum_input(_reference_config_enum_input),
    result(_inserter) {}

  /// Create ConfigEnumInput with unique cluster sites selected, starting orbits
  ///
  /// Delegates to standalone `select_cluster_sites`
  template<typename Inserter>
  template<typename OrbitType>
  void ClusterSitesSelector<Inserter>::operator()(std::vector<OrbitType> const &orbits) const {
    result = select_cluster_sites(reference_config_enum_input, orbits, result);
  }

  template<typename Inserter>
  ClusterSitesSelector<Inserter> make_cluster_sites_selector(
    ConfigEnumInput const &reference_config_enum_input,
    Inserter inserter) {
    return ClusterSitesSelector<Inserter> {reference_config_enum_input, inserter};
  }

}

#endif
