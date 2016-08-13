#include "casm/clusterography/ClusterFunctions.hh"

#include <string>

#include "casm/CASM_global_definitions.hh"
#include "casm/clusterography/Orbitree.hh"
//#include "casm/clusterography/HopCluster.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  jsonParser &to_json(const GenericOrbitBranch<SiteCluster> &branch, jsonParser &json) {
    return branch.to_json(json);
  }

  /// Assumes the prototype lattice is already set
  void from_json(GenericOrbitBranch<SiteCluster> &branch, const jsonParser &json) {
    branch.from_json(json);
  }

  jsonParser &to_json(const GenericOrbit<SiteCluster> &orbit, jsonParser &json) {
    return orbit.to_json(json);
  }
  /// Assumes the prototype lattice is already set
  void from_json(GenericOrbit<SiteCluster> &orbit, const jsonParser &json) {
    orbit.from_json(json);
  }

};

