#ifndef CASM_core_io_json_global_enum
#define CASM_core_io_json_global_enum

#include "casm/global/enum.hh"
#include "casm/casm_io/json_io/enum.hh"

namespace CASM {

  ENUM_JSON_IO_DECL(COORD_TYPE)
  ENUM_JSON_IO_DECL(PERIODICITY_TYPE)
  ENUM_JSON_IO_DECL(EQUIVALENCE_TYPE)
  ENUM_JSON_IO_DECL(CELL_TYPE)
  ENUM_JSON_IO_DECL(OnError)

}

#endif
