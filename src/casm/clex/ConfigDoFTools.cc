#include "casm/clex/ConfigDoF.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  /// Construct zero-valued ConfigDoF
  ConfigDoF make_configdof(Structure const &prim, Index volume) {
    return ConfigDoF(prim.basis().size(),
                     volume,
                     global_dof_info(prim),
                     local_dof_info(prim),
                     prim.occupant_symrep_IDs(),
                     prim.lattice().tol());
  }

}
