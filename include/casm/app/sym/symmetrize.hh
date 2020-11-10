#ifndef CASM_sym_symmetrize
#define CASM_sym_symmetrize

#include <string>

namespace CASM {

  class APICommandBase;
  class jsonParser;

  /// Describe the symmetrize method
  std::string symmetrize_desc();

  /// Adjust a structure's lattice and basis to increase factor group symmetry
  void symmetrize(PrimClex &primclex,
                  jsonParser const &json_options,
                  jsonParser const &cli_options_as_json);

}

#endif
