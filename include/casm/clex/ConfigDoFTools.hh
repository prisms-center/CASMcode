#ifndef CASM_clex_ConfigDoFTools
#define CASM_clex_ConfigDoFTools

#include "casm/global/definitions.hh"

namespace CASM {

class ConfigDoF;
class Structure;

/// Construct zero-valued ConfigDoF
ConfigDoF make_configdof(Structure const &prim, Index volume);

}  // namespace CASM

#endif
