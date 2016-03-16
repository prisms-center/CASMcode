#ifndef CASM_jsonIO_global
#define CASM_jsonIO_global

#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/jsonParser.hh"

namespace CASM {

  jsonParser &to_json(const COORD_TYPE &value, jsonParser &json);
  void from_json(COORD_TYPE &value, const jsonParser &json);

  jsonParser &to_json(const PERIODICITY_TYPE &value, jsonParser &json);
  void from_json(PERIODICITY_TYPE &value, const jsonParser &json);

  jsonParser &to_json(const CELL_TYPE &value, jsonParser &json);
  void from_json(CELL_TYPE &value, const jsonParser &json);

}

#endif
