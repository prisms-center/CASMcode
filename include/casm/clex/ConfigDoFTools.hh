#ifndef CASM_clex_ConfigDoFTools
#define CASM_clex_ConfigDoFTools

#include <memory>

#include "casm/global/definitions.hh"

namespace CASM {

class ConfigDoF;
class Structure;
class Supercell;

/// Construct zero-valued ConfigDoF
ConfigDoF make_configdof(Structure const &prim, Index volume);

/// Construct zero-valued ConfigDoF
ConfigDoF make_configdof(Structure const &prim, Index volume, double tol);

/// Construct zero-valued std::unique_ptr<ConfigDoF>
std::unique_ptr<ConfigDoF> make_unique_configdof(Structure const &prim,
                                                 Index volume);

/// Construct zero-valued std::unique_ptr<ConfigDoF>
std::unique_ptr<ConfigDoF> make_unique_configdof(Structure const &prim,
                                                 Index volume, double tol);

/// Construct zero-valued ConfigDoF
ConfigDoF make_configdof(Supercell const &supercell);

/// Construct zero-valued ConfigDoF
ConfigDoF make_configdof(Supercell const &supercell, double tol);

/// Construct zero-valued std::unique_ptr<ConfigDoF>
std::unique_ptr<ConfigDoF> make_unique_configdof(Supercell const &supercell);

/// Construct zero-valued std::unique_ptr<ConfigDoF>
std::unique_ptr<ConfigDoF> make_unique_configdof(Supercell const &supercell,
                                                 double tol);

}  // namespace CASM

#endif
