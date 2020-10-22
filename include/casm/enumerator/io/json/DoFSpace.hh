#ifndef CASM_symmetry_io_json_DoFSpace
#define CASM_symmetry_io_json_DoFSpace

#include "casm/global/definitions.hh"


namespace CASM {

  struct DoFSpace;

  jsonParser &to_json(DoFSpace const &dofspace, jsonParser &json, std::string name);

}

#endif
