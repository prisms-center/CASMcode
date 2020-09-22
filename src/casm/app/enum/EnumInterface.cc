
#include "casm/app/enum/EnumInterface.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  /// Convert `casm enum` CLI input to JSON
  ///
  /// All are optionally present, if present on command line
  /// \code
  /// {
  ///   "desc": <array of string, name of enumeration method to print the description>,
  ///   "help": <bool, print/return help>,
  ///   "method": <string, name of enumeration method>,
  ///   "settings": <string, represents path to settings JSON file>,
  ///   "input": <string, a JSON string>,
  ///   "min": <int, min supercell volume to do enumerations>,
  ///   "max": <int, max supercell volume to do enumerations>,
  ///   "filter": <array of string, filter query/select expression, save configs if evaluates true>,
  ///   "all": <bool, use all existing supercells for enumeration>,
  ///   "verbosity": <string, to be read by Log::verbosity_level>,
  ///   "scelnames": <array of string, list of supercell names, context dependent usage>,
  ///   "confignames": <array of string, list of config names, context dependent usage>,
  ///   "dry_run": <bool, print/return method results but do not save results>
  /// }
  /// \endcode
  jsonParser &to_json(const Completer::EnumOption &enum_opt, jsonParser &json) {
    const auto &vm = enum_opt.vm();

    json.put_obj();
    if(vm.count("desc")) {
      json["desc"] = enum_opt.desc_vec(); //vector<std::string>
    }
    if(vm.count("help")) {
      json["help"] = static_cast<bool>(vm.count("help")); //bool
    }
    if(vm.count("method")) {
      json["method"] = enum_opt.method(); //str
    }
    if(vm.count("settings")) {
      json["settings"] = enum_opt.settings_path().string(); //str
    }
    if(vm.count("input")) {
      json["input"] = enum_opt.input_str(); //str
    }

    if(vm.count("min")) {
      json["min"] = enum_opt.min_volume(); //int
    }
    if(vm.count("max")) {
      json["max"] = enum_opt.max_volume(); //int
    }
    if(vm.count("filter")) {
      json["filter"] = enum_opt.filter_strs(); //vector<std::string>
    }
    if(vm.count("all")) {
      json["all"] = enum_opt.all_existing(); //bool
    }
    if(vm.count("verbosity")) {
      json["verbosity"] = enum_opt.verbosity_str(); //str
    }
    if(vm.count("scelnames")) {
      json["scelnames"] = enum_opt.supercell_strs(); //vector<std::string>
    }
    if(vm.count("confignames")) {
      json["confignames"] = enum_opt.config_strs(); //vector<std::string>
    }
    if(vm.count("dry-run")) {
      json["dry_run"] = enum_opt.dry_run(); //bool
    }
    return json;
  }

}
