#ifndef CASM_core_io_json_global_enum
#define CASM_core_io_json_global_enum

#include "casm/casm_io/enum/json_io.hh"
#include "casm/global/enum.hh"
#include "casm/global/enum/io_traits.hh"

namespace CASM {

ENUM_JSON_IO_DECL(COORD_TYPE)
ENUM_JSON_IO_DECL(PERIODICITY_TYPE)
ENUM_JSON_IO_DECL(EQUIVALENCE_TYPE)
ENUM_JSON_IO_DECL(CELL_TYPE)
ENUM_JSON_IO_DECL(OnError)

}  // namespace CASM

#endif
