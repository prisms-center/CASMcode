#ifndef CASM_ClusterSpecs_impl
#define CASM_ClusterSpecs_impl

#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/ClusterSpecs.hh"
#include "casm/crystallography/Site.hh"
#include "casm/symmetry/Orbit.hh"

#include "casm/clusterography/ClusterOrbits_impl.hh"

namespace CASM {

  template<typename FunctorType>
  void for_all_orbits(ClusterSpecs const &cluster_specs,
                      std::vector<IntegralCluster> const &generating_elements,
                      FunctorType const &f) {
    if(cluster_specs.periodicity_type() == CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC) {
      auto orbits = cluster_specs.make_periodic_orbits(generating_elements);
      f(orbits);
      return;
    }
    else if(cluster_specs.periodicity_type() == CLUSTER_PERIODICITY_TYPE::LOCAL) {
      auto orbits = cluster_specs.make_local_orbits(generating_elements);
      f(orbits);
      return;
    }
    else if(cluster_specs.periodicity_type() == CLUSTER_PERIODICITY_TYPE::WITHIN_SCEL) {
      auto orbits = cluster_specs.make_within_scel_orbits(generating_elements);
      f(orbits);
      return;
    }

    throw libcasm_runtime_error("Error: unsupported orbit type");
  }

  template<typename FunctorType>
  void for_all_orbits(ClusterSpecs const &cluster_specs, std::ostream &status, FunctorType const &f) {
    if(cluster_specs.periodicity_type() == CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC) {
      auto orbits = cluster_specs.make_periodic_orbits(status);
      f(orbits);
      return;
    }
    else if(cluster_specs.periodicity_type() == CLUSTER_PERIODICITY_TYPE::LOCAL) {
      auto orbits = cluster_specs.make_local_orbits(status);
      f(orbits);
      return;
    }
    else if(cluster_specs.periodicity_type() == CLUSTER_PERIODICITY_TYPE::WITHIN_SCEL) {
      auto orbits = cluster_specs.make_within_scel_orbits(status);
      f(orbits);
      return;
    }

    throw libcasm_runtime_error("Error: unsupported orbit type");
  }
}

#endif
