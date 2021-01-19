#include "casm/global/enum/json_io.hh"

#include "casm/casm_io/json/jsonParser.hh"
#include "casm/global/enum/io_traits.hh"

namespace CASM {

ENUM_JSON_IO_DEF(COORD_TYPE)
ENUM_JSON_IO_DEF(PERIODICITY_TYPE)
ENUM_JSON_IO_DEF(EQUIVALENCE_TYPE)
ENUM_JSON_IO_DEF(CELL_TYPE)
ENUM_JSON_IO_DEF(OnError)

}  // namespace CASM
