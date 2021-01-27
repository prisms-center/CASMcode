#include "casm/clex/ConfigDoFTools.hh"

#include "casm/clex/ConfigDoF_impl.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

/// Construct zero-valued ConfigDoF
ConfigDoF make_configdof(Structure const &prim, Index volume) {
  return ConfigDoF(prim.basis().size(), volume, global_dof_info(prim),
                   local_dof_info(prim), prim.occupant_symrep_IDs(),
                   prim.lattice().tol());
}

/// Construct zero-valued ConfigDoF
ConfigDoF make_configdof(Structure const &prim, Index volume, double tol) {
  return ConfigDoF(prim.basis().size(), volume, global_dof_info(prim),
                   local_dof_info(prim), prim.occupant_symrep_IDs(), tol);
}

/// Construct zero-valued std::unique_ptr<ConfigDoF>
std::unique_ptr<ConfigDoF> make_unique_configdof(Structure const &prim,
                                                 Index volume) {
  return notstd::make_unique<ConfigDoF>(
      prim.basis().size(), volume, global_dof_info(prim), local_dof_info(prim),
      prim.occupant_symrep_IDs(), prim.lattice().tol());
}

/// Construct zero-valued std::unique_ptr<ConfigDoF>
std::unique_ptr<ConfigDoF> make_unique_configdof(Structure const &prim,
                                                 Index volume, double tol) {
  return notstd::make_unique<ConfigDoF>(
      prim.basis().size(), volume, global_dof_info(prim), local_dof_info(prim),
      prim.occupant_symrep_IDs(), tol);
}

/// Construct zero-valued ConfigDoF
ConfigDoF make_configdof(Supercell const &supercell) {
  return make_configdof(supercell.prim(), supercell.volume());
}

/// Construct zero-valued ConfigDoF
ConfigDoF make_configdof(Supercell const &supercell, double tol) {
  return make_configdof(supercell.prim(), supercell.volume(), tol);
}

/// Construct zero-valued std::unique_ptr<ConfigDoF>
std::unique_ptr<ConfigDoF> make_unique_configdof(Supercell const &supercell) {
  return make_unique_configdof(supercell.prim(), supercell.volume());
}

/// Construct zero-valued std::unique_ptr<ConfigDoF>
std::unique_ptr<ConfigDoF> make_unique_configdof(Supercell const &supercell,
                                                 double tol) {
  return make_unique_configdof(supercell.prim(), supercell.volume(), tol);
}

}  // namespace CASM
