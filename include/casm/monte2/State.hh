#ifndef CASM_monte2_State
#define CASM_monte2_State

#include <map>
#include <string>

#include "casm/clex/Configuration.hh"
#include "casm/global/eigen.hh"
#include "casm/monte2/Definitions.hh"

namespace CASM {
namespace Monte2 {

/// A state of a Monte Carlo calculation
struct State {
  State(Configuration const &_configuration,
        VectorValueMap _conditions = VectorValueMap(),
        VectorValueMap _properties = VectorValueMap())
      : configuration(_configuration),
        conditions(_conditions),
        properties(_properties) {}

  /// Current configuration
  Configuration configuration;

  /// Conditions of the state
  ///
  /// Thermodynamic conditions or calculation constraints, such as temperature,
  /// chemical potential (for grand canonical Monte Carlo), composition (for
  /// canonical Monte Carlo), etc., depending on the type of Monte Carlo
  /// calculation
  VectorValueMap conditions;

  /// Properties of the state
  ///
  /// Properties of the state could be formation_energy, potential_energy,
  /// comp_n, etc., depending on the type of Monte Carlo calculation.
  VectorValueMap properties;
};

}  // namespace Monte2
}  // namespace CASM

#endif
