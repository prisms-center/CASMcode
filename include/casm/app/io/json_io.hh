#ifndef CASM_app_json_io
#define CASM_app_json_io

#include <map>

namespace CASM {

class jsonParser;

/// Check for --settings, then --input and return JSON input. Else return empty
/// JSON object.
template <typename OptionType>
jsonParser make_json_input(const OptionType &opt);

/// Copy from `json_source` to `json_combined`
jsonParser &combine_json_options(
    std::map<std::string, std::string> const &source_to_combined_keys,
    jsonParser const &json_source, jsonParser &json_combined);

}  // namespace CASM

#endif
