#ifndef CASM_core_io_stream_global_enum
#define CASM_core_io_stream_global_enum

#include "casm/global/enum.hh"
#include "casm/global/enum/io_traits.hh"
#include "casm/casm_io/enum/stream_io.hh"

namespace CASM {

  ENUM_IO_DECL(COORD_TYPE)
  ENUM_IO_DECL(PERIODICITY_TYPE)
  ENUM_IO_DECL(EQUIVALENCE_TYPE)
  ENUM_IO_DECL(CELL_TYPE)
  ENUM_IO_DECL(OnError)

}

#endif
