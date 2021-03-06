#ifndef CASM_enumerator_ClusterSitesSelector
#define CASM_enumerator_ClusterSitesSelector

#include "casm/enumerator/ConfigEnumInput.hh"

namespace CASM {

/// Create ConfigEnumInput with unique cluster sites selected, starting from
/// "prim_periodic" orbits
template <typename Inserter>
Inserter select_cluster_sites(
    ConfigEnumInput const &reference_config_enum_input,
    std::vector<PrimPeriodicIntegralClusterOrbit> const &orbits,
    Inserter result);

}  // namespace CASM

#endif
