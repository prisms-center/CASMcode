#include "casm/monte_carlo/MonteCorrelations.hh"

#include "casm/clex/Clexulator.hh"
#include "casm/clex/ConfigCorrelations.hh"
#include "casm/monte_carlo/Conversions.hh"
#include "casm/monte_carlo/OccLocation.hh"

namespace CASM {
namespace Monte {

/// \brief Sets change in (extensive) correlations due to occupation
/// changes, restricted to specified correlations
///
/// \param dcorr, Eigen::VectorXd of change in correlations. Will be set to
/// size `clexulator.corr_size()` if necessary.  Only elements corresponding to
/// indices in `correlations_indices` will be modified.
///
void restricted_delta_corr(Eigen::VectorXd &dcorr, OccEvent const &occ_event,
                           Conversions const &convert,
                           ConfigDoF const &configdof,
                           SuperNeighborList const &supercell_neighbor_list,
                           Clexulator const &clexulator,
                           unsigned int const *corr_indices_begin,
                           unsigned int const *corr_indices_end) {
  const OccEvent &e = occ_event;

  // we'll revert changes when we're done so in the end nothing changes
  ConfigDoF &mutable_configdof = const_cast<ConfigDoF &>(configdof);

  if (e.occ_transform.size() == 0) {
    for (auto it = corr_indices_begin; it != corr_indices_end; ++it) {
      *(dcorr.data() + *it) = 0.0;
    }
    return;
  }

  static std::vector<int> curr_occ;
  curr_occ.resize(e.occ_transform.size());

  // first swap
  OccTransform const &t = e.occ_transform[0];
  curr_occ[0] = configdof.occ(t.l);
  int new_occ = convert.occ_index(t.asym, t.to_species);
  restricted_delta_corr(dcorr, t.l, new_occ, configdof, supercell_neighbor_list,
                        clexulator, corr_indices_begin, corr_indices_end);
  mutable_configdof.occ(t.l) = new_occ;

  // subsequent swaps
  static Eigen::VectorXd tmp_dcorr;
  for (Index i = 1; i < e.occ_transform.size(); ++i) {
    OccTransform const &t = e.occ_transform[i];
    curr_occ[i] = configdof.occ(t.l);
    new_occ = convert.occ_index(t.asym, t.to_species);
    restricted_delta_corr(tmp_dcorr, t.l, new_occ, configdof,
                          supercell_neighbor_list, clexulator,
                          corr_indices_begin, corr_indices_end);
    mutable_configdof.occ(t.l) = new_occ;
    for (auto it = corr_indices_begin; it != corr_indices_end; ++it) {
      *(dcorr.data() + *it) += *(tmp_dcorr.data() + *it);
    }
  }

  // revert changes
  for (Index i = 0; i < e.occ_transform.size(); ++i) {
    OccTransform const &t = e.occ_transform[i];
    mutable_configdof.occ(t.l) = curr_occ[i];
  }
}

}  // namespace Monte
}  // namespace CASM
