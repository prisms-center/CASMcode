#ifndef CASM_config_Configuration_impl
#define CASM_config_Configuration_impl

#include "casm/clexulator/ConfigDoFValuesTools.hh"
#include "casm/configuration/Configuration.hh"

namespace CASM {
namespace config {

Configuration::Configuration(std::shared_ptr<Supercell const> const &_supercell)
    : supercell(_supercell),
      dof_values(clexulator::make_default_config_dof_values(
          _supercell->prim->basicstructure.basis().size(),  // Index N_sublat
          _supercell->superlattice.size(),                  // Index N_volume
          global_dof_info(_supercell->prim->basicstructure),
          local_dof_info(_supercell->prim->basicstructure))) {}

Configuration::Configuration(std::shared_ptr<Supercell const> const &_supercell,
                             clexulator::ConfigDoFValues const &_dof_values)
    : supercell(_supercell), dof_values(_dof_values) {}

}  // namespace config
}  // namespace CASM

#endif
