#ifndef CASM_core_io_traits_global_enum
#define CASM_core_io_traits_global_enum

#include "casm/casm_io/enum/io_traits.hh"
#include "casm/global/enum.hh"

namespace CASM {

ENUM_TRAITS(COORD_TYPE)
ENUM_TRAITS(PERIODICITY_TYPE)
ENUM_TRAITS(EQUIVALENCE_TYPE)
ENUM_TRAITS(CELL_TYPE)
ENUM_TRAITS(OnError)

}  // namespace CASM

#endif
