#ifndef CASM_config_Configuration
#define CASM_config_Configuration

#include "casm/clexulator/ConfigDoFValues.hh"
#include "casm/configuration/Supercell.hh"
#include "casm/configuration/definitions.hh"

namespace CASM {
namespace config {

/// \brief Data structure encapsulating configuration DoF values and all
/// information necessary to apply symmetry operations
struct Configuration {
  Configuration(std::shared_ptr<Supercell const> const &_supercell);

  Configuration(std::shared_ptr<Supercell const> const &_supercell,
                clexulator::ConfigDoFValues const &_dof_values);

  std::shared_ptr<Supercell const> supercell;

  clexulator::ConfigDoFValues dof_values;
};

}  // namespace config
}  // namespace CASM

#endif
