#ifndef CASM_app_json_io_impl
#define CASM_app_json_io_impl

#include "casm/casm_io/json/jsonParser.hh"

namespace CASM {

/// Check for --settings, then --input and return JSON input. Else return empty
/// JSON object.
template <typename OptionType>
jsonParser make_json_input(const OptionType &opt) {
  if (!opt.settings_path().empty()) {
    return jsonParser{opt.settings_path()};
  } else if (!opt.input_str().empty()) {
    return jsonParser::parse(opt.input_str());
  } else {
    return jsonParser::object();
  }
}

}  // namespace CASM

#endif
