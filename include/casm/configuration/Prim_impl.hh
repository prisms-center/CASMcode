#ifndef CASM_config_Prim_impl
#define CASM_config_Prim_impl

#include "casm/configuration/Prim.hh"

namespace CASM {
namespace config {

Prim::Prim(BasicStructure const &_basicstructure)
    : basicstructure(_basicstructure), sym_info(_basicstructure) {}

}  // namespace config
}  // namespace CASM

#endif
