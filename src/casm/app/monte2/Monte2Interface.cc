
#include "casm/app/monte2/Monte2Interface.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

/// Convert `casm monte2` CLI input to JSON
///
/// All are optionally present, if present on command line
/// \code
/// {
///   "desc": <array of string, name of enumeration method to print the
///       description>,
///   "help": <bool, print/return help>, "method": <string, name of enumeration
///       method>,
///   "settings": <string, represents path to settings JSON file>,
///   "input": <string, a JSON string>,
///   "verbosity": <string, to be read by Log::verbosity_level>
/// }
/// \endcode
jsonParser &to_json(const Completer::Monte2Option &enum_opt, jsonParser &json) {
  const auto &vm = enum_opt.vm();

  json.put_obj();
  if (vm.count("desc")) {
    json["desc"] = enum_opt.desc_vec();  // vector<std::string>
  }
  if (vm.count("help")) {
    json["help"] = static_cast<bool>(vm.count("help"));  // bool
  }
  if (vm.count("method")) {
    json["method"] = enum_opt.method();  // str
  }
  if (vm.count("settings")) {
    json["settings"] = enum_opt.settings_path().string();  // str
  }
  if (vm.count("input")) {
    json["input"] = enum_opt.input_str();  // str
  }

  if (vm.count("verbosity") && !vm["verbosity"].defaulted()) {
    json["verbosity"] = enum_opt.verbosity_str();  // str
  }
  return json;
}

/// A standard approach to combine CLI options with user JSON options into one
/// JSON document
///
/// Not all monte2 interfaces need to incorporate the CLI output using this
/// method, but it is the most common way to do so.
jsonParser combine_monte2_json_options(jsonParser const &json_options,
                                       jsonParser const &cli_options_as_json) {
  jsonParser json_combined{json_options};
  if (cli_options_as_json.contains("verbosity")) {
    json_combined["verbosity"] = cli_options_as_json["verbosity"];
  }
  return json_combined;
}

/// Combine --input / --settings JSON with CLI options
ParentInputParser make_monte2_parent_parser(
    Log &log, jsonParser const &json_options,
    jsonParser const &cli_options_as_json) {
  log.indent() << "Input from JSON (--input or --setings):\n"
               << json_options << std::endl
               << std::endl;
  log.indent() << "Input from `casm monte2` options:\n"
               << cli_options_as_json << std::endl
               << std::endl;

  // combine JSON options and CLI options
  jsonParser json_combined =
      combine_monte2_json_options(json_options, cli_options_as_json);

  log.indent() << "Combined Input:\n"
               << json_combined << std::endl
               << std::endl;
  return ParentInputParser{json_combined};
}

}  // namespace CASM
