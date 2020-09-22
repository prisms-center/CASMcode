#ifndef CASM_sym_write_prim_symmetry
#define CASM_sym_write_prim_symmetry

namespace CASM {

  class APICommandBase;
  class jsonParser;

  /// Describe the default `casm sym` option
  std::string write_prim_symmetry_desc();

  /// Write/print prim symmetry (ignores JSON input)
  void write_prim_symmetry(APICommandBase const &cmd,
                           jsonParser const &json_options,
                           jsonParser const &cli_options_as_json);

}

#endif
