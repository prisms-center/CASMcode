#include "casm/app/sym/json_io.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  /// Convert `casm sym` CLI input to JSON
  ///
  /// All are optionally present, if present on command line
  /// \code
  /// {
  ///   "desc": <array of string, name of enumeration method to print the description>,
  ///   "help": <bool, print/return help>,
  ///   "print_lattice_point_group": <bool, true to print lattice point group>,
  ///   "print_factor_group": <bool, true to print factor group>,
  ///   "print_crystal_point_group": <bool, true to print crystal point group>,
  ///   "coordinate_mode": <string, coordinate type for printing>,
  ///   "symmetrize": <string, path to POSCAR file to symmetrize>,
  ///   "tol": <number, enforced tolerance for symmetrization>,
  ///   "dof_space_analysis": <bool, true to run dof_space_analysis>,
  ///   "settings": <string, represents path to settings JSON file>,
  ///   "input": <string, a JSON string>,
  ///   "scelnames": <array of string, list of supercell names, context dependent usage>,
  ///   "confignames": <array of string, list of config names, context dependent usage>,
  ///   "selection": <string, configuration selection name/path>,
  ///   "dofs": <array of string, DoF types for dof_space_analysis>,
  ///   "calc_wedge": <bool, true to calculate DoF space wedge>,
  /// }
  /// \endcode
  jsonParser &to_json(const Completer::SymOption &sym_opt, jsonParser &json) {
    const auto &vm = sym_opt.vm();

    json.put_obj();
    if(vm.count("desc")) {
      json["desc"] = static_cast<bool>(vm.count("help")); //bool
    }
    if(vm.count("help")) {
      json["help"] = static_cast<bool>(vm.count("help")); //bool
    }
    if(vm.count("lattice-point-group")) {
      json["print_lattice_point_group"] = static_cast<bool>(vm.count("lattice-point-group")); //bool
    }
    if(vm.count("factor-group")) {
      json["print_factor_group"] = static_cast<bool>(vm.count("factor-group")); //bool
    }
    if(vm.count("crystal-point-group")) {
      json["print_crystal_point_group"] = static_cast<bool>(vm.count("crystal-point-group")); //bool
    }
    if(vm.count("coord")) {
      json["coordinate_mode"] = sym_opt.coordtype_str(); //string
    }
    if(vm.count("symmetrize")) {
      json["symmetrize"] = sym_opt.poscar_path().string(); //string
    }
    if(vm.count("tol")) {
      json["tol"] = sym_opt.tol(); //number
    }
    if(vm.count("dof-space-analysis")) {
      json["dof_space_analysis"] = static_cast<bool>(vm.count("dof-space-analysis")); //bool
    }
    if(vm.count("settings")) {
      json["settings"] = sym_opt.settings_path().string(); //string
    }
    if(vm.count("input")) {
      json["input"] = sym_opt.input_str(); //string
    }
    if(vm.count("scelnames")) {
      json["scelnames"] = sym_opt.supercell_strs(); //vector<std::string>
    }
    if(vm.count("confignames")) {
      json["confignames"] = sym_opt.config_strs(); //vector<std::string>
    }
    if(vm.count("selection")) {
      json["selection"] = sym_opt.selection_path().string(); //std::string
    }
    if(vm.count("dofs")) {
      json["dofs"] = sym_opt.dof_strs(); //vector<std::string>
    }
    if(vm.count("calc-wedge")) {
      json["calc_wedge"] = static_cast<bool>(vm.count("calc-wedge")); //bool
    }
    return json;
  }
}
