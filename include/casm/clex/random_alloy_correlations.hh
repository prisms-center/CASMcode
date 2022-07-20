#ifndef CASM_clex_random_alloy_correlations
#define CASM_clex_random_alloy_correlations

#include "casm/casm_io/Log.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/clusterography/IntegralClusterSymCompareTraits.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

/// \brief Make expected random alloy correlations
template <typename ClusterOrbitIterator>
Eigen::VectorXd make_random_alloy_correlations(
    ClexBasis const &clex_basis, ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    std::vector<Eigen::VectorXd> const &sublattice_prob);

struct RandomAlloyCorrCalculator {
  /// \brief Constructor, using saved orbit prototypes
  RandomAlloyCorrCalculator(std::shared_ptr<Structure const> const &shared_prim,
                            ClexBasisSpecs const &basis_set_specs,
                            std::vector<IntegralCluster> const &prototypes);

  /// \brief Constructor, generates orbits from specs
  RandomAlloyCorrCalculator(std::shared_ptr<Structure const> const &shared_prim,
                            ClexBasisSpecs const &basis_set_specs, Log &log);

  /// \brief Calculate random alloy correlations for a given
  ///     sublattice composition
  Eigen::VectorXd operator()(
      std::vector<Eigen::VectorXd> const &sublattice_prob) const;

  std::vector<PrimPeriodicOrbit<IntegralCluster>> orbits;

  ParsingDictionary<DoFType::Traits> const *dof_dict;

  ClexBasis clex_basis;
};

}  // namespace CASM

#endif
