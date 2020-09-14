#include "casm/enumerator/ClusterSitesSelector_impl.hh"

namespace CASM {

  /// Create ConfigEnumInput with unique cluster sites selected
  ///
  /// This will use "cluster_specs" to create orbits, then for each orbit:
  /// - Copy reference_config_enum_input (keeping selected any sites already selected)
  /// - Select all orbit prototype sites on the copy
  /// - Insert the copy into the container that will be returned
  ///
  /// Note:
  /// - This will include a copy of reference_config_enum_input if the null cluster orbit is
  ///   generated, as is typical
  /// - ClusterSitesSelector is only implemented for "within_scel" and "prim_periodic" orbits. For
  ///   other orbit types this will throw.
  /// - TODO: Implement for "local" orbits
  std::vector<ConfigEnumInput> select_cluster_sites(
    ConfigEnumInput const &reference_config_enum_input,
    ClusterSpecs const &cluster_specs) {

    std::vector<ConfigEnumInput> with_cluster_sites;

    auto selector = make_cluster_sites_selector(
                      reference_config_enum_input,
                      std::back_inserter(with_cluster_sites));

    // This handles checking what type of orbits ClusterSpecs generates, making them, and
    //   calling selector(orbits);
    for_all_orbits(cluster_specs, null_log(), selector);

    return with_cluster_sites;
  }

}
