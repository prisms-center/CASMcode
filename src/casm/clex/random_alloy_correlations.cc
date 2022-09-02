#include "casm/basis_set/DoFTraits.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clex/random_alloy_correlations_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"

namespace CASM {

/// \brief Constructor, using saved orbit prototypes
///
/// \code
/// std::vector<IntegralCluster> prototypes;
/// jsonParser clust_json{primclex.dir().clust(basis_set_name)};
/// read_clust(std::back_inserter(prototypes), clust_json, *shared_prim);
/// \endcode
RandomAlloyCorrCalculator::RandomAlloyCorrCalculator(
    std::shared_ptr<Structure const> const &shared_prim,
    ClexBasisSpecs const &basis_set_specs,
    std::vector<IntegralCluster> const &prototypes)
    : orbits(basis_set_specs.cluster_specs->make_periodic_orbits(prototypes)),
      dof_dict(&DoFType::traits_dict()),
      clex_basis(shared_prim, basis_set_specs, dof_dict) {
  clex_basis.generate(orbits.begin(), orbits.end());
}

/// \brief Constructor, generates orbits from specs
RandomAlloyCorrCalculator::RandomAlloyCorrCalculator(
    std::shared_ptr<Structure const> const &shared_prim,
    ClexBasisSpecs const &basis_set_specs, Log &log)
    : orbits(basis_set_specs.cluster_specs->make_periodic_orbits(log)),
      dof_dict(&DoFType::traits_dict()),
      clex_basis(shared_prim, basis_set_specs, dof_dict) {
  clex_basis.generate(orbits.begin(), orbits.end());
}

/// \brief Calculate random alloy correlations for a given
///     sublattice composition
Eigen::VectorXd RandomAlloyCorrCalculator::operator()(
    std::vector<Eigen::VectorXd> const &sublattice_prob) const {
  return make_random_alloy_correlations(clex_basis, orbits.begin(),
                                        orbits.end(), sublattice_prob);
}

}  // namespace CASM
