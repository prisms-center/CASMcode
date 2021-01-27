#ifndef CASM_sym_write_prim_symmetry
#define CASM_sym_write_prim_symmetry

#include <string>

namespace CASM {

class APICommandBase;
class PrimClex;
class jsonParser;

/// Describe the default `casm sym` option
std::string write_prim_symmetry_desc();

/// Write/print prim symmetry
void write_prim_symmetry(PrimClex &primclex, jsonParser const &json_options,
                         jsonParser const &cli_options_as_json);

}  // namespace CASM

#endif
