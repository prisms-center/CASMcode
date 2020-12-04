#ifndef CASM_sym_dof_space_analysis
#define CASM_sym_dof_space_analysis

#include <string>

namespace CASM {

  class APICommandBase;
  class PrimClex;
  class jsonParser;

  /// Describe DoF space analysis input
  std::string dof_space_analysis_desc();

  /// Perform DoF space analysis
  void dof_space_analysis(PrimClex &primclex,
                          jsonParser const &json_options,
                          jsonParser const &cli_options_as_json);

}

#endif
