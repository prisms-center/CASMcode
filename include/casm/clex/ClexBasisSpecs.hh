#ifndef CASM_ClexBasisSpecs
#define CASM_ClexBasisSpecs

#include "casm/basis_set/BasisFunctionSpecs.hh"
#include "casm/clusterography/ClusterSpecs.hh"

namespace CASM {

  // --- Specify type of orbits ---

  /// Provides parameters for constructing a cluster expansion basis (ClexBasis)
  struct ClexBasisSpecs {

    ClexBasisSpecs(
      BasisFunctionSpecs const &_basis_function_specs,
      notstd::cloneable_ptr<ClusterSpecs> _cluster_specs):
      basis_function_specs(_basis_function_specs),
      cluster_specs(std::move(_cluster_specs)) {
      if(!cluster_specs) {
        throw libcasm_runtime_error("Error constructing ClexBasisSpecs: Empty \"cluster_specs\".");
      }
    }

    BasisFunctionSpecs basis_function_specs;

    notstd::cloneable_ptr<ClusterSpecs> cluster_specs;

  };

}
#endif
