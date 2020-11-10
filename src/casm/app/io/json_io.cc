#include "casm/app/io/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {

  /// Copy from `json_source` to `json_combined`
  jsonParser &combine_json_options(
    std::map<std::string, std::string> const &source_to_combined_keys,
    jsonParser const &json_source,
    jsonParser &json_combined) {

    for(auto const &pair : source_to_combined_keys) {
      auto it = json_source.find(pair.first);
      if(it != json_source.end()) {
        json_combined[pair.second] = *it;
      }
    }
    return json_combined;
  }
}
