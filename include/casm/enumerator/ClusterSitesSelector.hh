#ifndef CASM_enumerator_ClusterSitesSelector
#define CASM_enumerator_ClusterSitesSelector

#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

  /// Create ConfigEnumInput with unique cluster sites selected, starting from "prim_periodic" orbits
  template<typename Inserter>
  Inserter select_cluster_sites(
    ConfigEnumInput const &reference_config_enum_input,
    std::vector<PrimPeriodicIntegralClusterOrbit> const &orbits,
    Inserter result);

  // /// Functor class to select sites on ConfigEnumInput for any supported orbit type
  // ///
  // /// This is intended for use with the `for_all_orbits` functor. Example:
  // /// \code
  // /// // A vector of reference ConfigEnumInput, before selecting cluster sites
  // /// std::vector<ConfigEnumInput> config_enum_input = ...
  // /// // ClusterSpecs, which can generate "within_scel" or "prim_periodic" orbits
  // /// std::unique_ptr<ClusterSpecs> cluster_specs = ...
  // /// // Will store ConfigEnumInput after selecting cluster sites
  // /// std::vector<ConfigEnumInput> with_cluster_sites;
  // /// for(ConfigEnumInput &input : initial_config_enum_input) {
  // ///   auto selector = make_cluster_sites_selector(input, std::back_inserter(with_cluster_sites));
  // ///   // This handles checking what type of orbits ClusterSpecs generates, making them, and
  // ///   //   calling selector(orbits);
  // ///   for_all_orbits(*cluster_specs, null_log(), selector);
  // /// }
  // /// // If ClusterSpecs generates a vector `orbits`, then expect:
  // /// //     with_cluster_sites.size() == initial_config_enum_input.size() * orbits.size()
  // /// \endcode
  // ///
  // template<typename Inserter>
  // class ClusterSitesSelector {
  // public:
  //   ClusterSitesSelector(ConfigEnumInput const &_reference_config_enum_input, Inserter _inserter);
  //
  //   template<typename OrbitType>
  //   void operator()(std::vector<OrbitType> const &orbits) const;
  //
  // private:
  //   ConfigEnumInput reference_config_enum_input;
  //   mutable Inserter result;
  // };
  //
  // template<typename Inserter>
  // ClusterSitesSelector<Inserter> make_cluster_sites_selector(
  //   ConfigEnumInput const &reference_config_enum_input,
  //   Inserter inserter);

  // /// Create ConfigEnumInput with unique cluster sites selected, starting from "within_scel" orbits
  // template<typename Inserter>
  // template<>
  // void ClusterSitesSelector<Inserter>::operator()<WithinScelIntegralClusterOrbit>(std::vector<WithinScelIntegralClusterOrbit> const &orbits) const;
  //
  // /// Create ConfigEnumInput with unique cluster sites selected, starting from "prim_periodic" orbits
  // template<typename Inserter>
  // template<>
  // void ClusterSitesSelector<Inserter>::operator()<PrimPeriodicIntegralClusterOrbit>(std::vector<PrimPeriodicIntegralClusterOrbit> const &orbits) const;

  // TODO:
}

#endif
