#ifndef CASM_sym_config_space_analysis
#define CASM_sym_config_space_analysis

#include <string>

namespace CASM {

class APICommandBase;
class PrimClex;
class jsonParser;

/// Describe --config-space-analysis options
std::string config_space_analysis_desc();

/// Run config space analysis
void config_space_analysis(PrimClex &primclex, jsonParser const &json_options,
                           jsonParser const &cli_options_as_json);

}  // namespace CASM

#endif
