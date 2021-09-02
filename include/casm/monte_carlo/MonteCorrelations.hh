#ifndef CASM_MonteCorrlations
#define CASM_MonteCorrlations

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

class Clexulator;
class ConfigDoF;
class Conversions;

namespace clexulator {
class SuperNeighborList;
}
using clexulator::SuperNeighborList;

namespace Monte {

struct OccEvent;

/// \brief Sets change in (extensive) correlations due to occupation
/// changes, restricted to specified correlations
void restricted_delta_corr(Eigen::VectorXd &dcorr, OccEvent const &occ_event,
                           Conversions const &convert,
                           ConfigDoF const &configdof,
                           SuperNeighborList const &supercell_neighbor_list,
                           Clexulator const &clexulator,
                           unsigned int const *corr_indices_begin,
                           unsigned int const *corr_indices_end);

}  // namespace Monte
}  // namespace CASM

#endif
