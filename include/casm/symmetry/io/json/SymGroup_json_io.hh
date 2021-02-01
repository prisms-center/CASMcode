#ifndef CASM_symmetry_io_json_SymGroup
#define CASM_symmetry_io_json_SymGroup

#include "casm/global/definitions.hh"

namespace CASM {

class SymGroup;
class jsonParser;

void write_symop(const SymGroup &grp, Index i, jsonParser &j);

void write_symgroup(const SymGroup &grp, jsonParser &json);

}  // namespace CASM

#endif
