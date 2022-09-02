#ifndef CASM_clex_random_alloy_correlations_impl
#define CASM_clex_random_alloy_correlations_impl

#include "casm/clex/random_alloy_correlations.hh"
#include "casm/container/Counter.hh"

namespace CASM {

/// \brief Make expected random alloy correlations
///
/// \param clex_basis The ClexBasis contains the cluster
///     expansion basis functions
/// \param begin, end The range of orbits used to generate
///     clex_basis
/// \param sublattice_prob The random alloy site probability
///     for each occupant on each basis site. The value
///     `sublattice_prob[b][i]` is the probability of the i-th
///     occupant on the b-th sublattice.
///
/// \returns A expected correlations for a random alloy with
///     the specified sublattice compositions
///
template <typename ClusterOrbitIterator>
Eigen::VectorXd make_random_alloy_correlations(
    ClexBasis const &clex_basis, ClusterOrbitIterator begin,
    ClusterOrbitIterator end,
    std::vector<Eigen::VectorXd> const &sublattice_prob) {
  auto shared_prim = clex_basis.shared_prim();

  // Currently only prim with only occ DoF are supported
  if (clex_basis.global_bases().size()) {
    throw std::runtime_error(
        "Error in make_random_alloy_correlations: currently only prim with "
        "only "
        "occupation DoF are supported");
  }
  for (auto const &dofset : clex_basis.site_bases()) {
    if (dofset.first != "occ") {
      throw std::runtime_error(
          "Error in make_random_alloy_correlations: currently only prim with "
          "only "
          "occupation DoF are supported");
    }
  }

  // This will hold correlations
  Eigen::VectorXd random_alloy_corr =
      Eigen::VectorXd::Zero(clex_basis.n_functions());

  // collect basis function info
  Index orbit_index = 0;
  Index function_index = 0;

  // Loop over cluster orbits
  for (auto orbit_it = begin; orbit_it != end; ++orbit_it) {
    BasisSet bset_prototype = clex_basis.bset_orbit(orbit_index)[0];

    // construct DoF handles (essentially pointers to the
    //    values evaluated by the basis functions)
    std::vector<DoF::RemoteHandle> handles;

    Index clust_size = orbit_it->prototype().size();
    std::vector<int> s(clust_size, 0);
    std::vector<int> max_occ(clust_size, 0);
    std::vector<int> sublat(clust_size, -1);

    for (int i = 0; i < clust_size; ++i) {
      handles.push_back(DoF::RemoteHandle("occ", "s", i));
      handles[i] = s[i];
      Index b = orbit_it->prototype()[i].sublattice();
      max_occ[i] = shared_prim->basis()[b].occupant_dof().size() - 1;
      sublat[i] = b;
    }
    bset_prototype.register_remotes(handles);

    // Loop over cluster basis functions associated
    //     with one cluster orbit
    for (Index i = 0; i < bset_prototype.size(); ++i) {
      if (clust_size == 0) {
        random_alloy_corr(function_index) = 1.0;
      } else {
        // loop over occupations on the prototype cluster
        Counter<std::vector<int>> counter(std::vector<int>(clust_size, 0),
                                          max_occ,
                                          std::vector<int>(clust_size, 1));
        // calculate expected value = sum_j p_j * phi
        double expected = 0.0;
        while (counter.valid()) {
          double p = 1.0;
          for (int j = 0; j < clust_size; ++j) {
            s[j] = counter[j];
            p *= sublattice_prob[sublat[j]][counter[j]];
          }
          expected += p * bset_prototype[i]->remote_eval();
          ++counter;
        }
        random_alloy_corr(function_index) = expected;
      }
      function_index++;
    }
    orbit_index++;
  }

  return random_alloy_corr;
}

}  // namespace CASM

#endif
