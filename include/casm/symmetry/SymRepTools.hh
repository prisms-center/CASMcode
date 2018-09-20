#ifndef CASM_SymRepTools
#define CASM_SymRepTools

#include "casm/symmetry/SymGroup.hh"
#include "casm/symmetry/SymGroupRep.hh"

namespace CASM {
  namespace SymRepTools {
    using Wedge = std::pair<std::vector<Index>, Eigen::MatrixXd>;

    std::vector<Wedge> irrep_wedges(SymGroup const &_group, SymGroupRepID id);
  }
}
#endif
