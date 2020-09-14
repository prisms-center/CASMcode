
#include "casm/app/enum.hh"
#include "casm/app/enum/EnumInterface_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/ConfigEnumRandomOccupations.hh"
#include "casm/clex/ConfigEnumRandomLocal.hh"
#include "casm/clex/ConfigEnumStrain.hh"
#include "casm/clex/ConfigEnumSiteDoFs.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/SuperConfigEnum.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"

namespace CASM {

  /// \brief Use to construct an InterfaceMap
  std::unique_ptr<InterfaceMap<Completer::EnumOption> > make_enumerator_map() {
    return make_interface_map<Completer::EnumOption>();
  }

  /// \brief Use to construct an EnumeratorMap with standard Enumerators (not plugins)
  std::unique_ptr<EnumeratorMap> make_standard_enumerator_map() {
    std::unique_ptr<EnumeratorMap> emap = make_enumerator_map();

    emap->insert(
      EnumInterface<ScelEnum>(),
      EnumInterface<ConfigEnumAllOccupations>(),
      EnumInterface<ConfigEnumRandomLocal>(),
      EnumInterface<ConfigEnumRandomOccupations>(),
      EnumInterface<ConfigEnumStrain>(),
      EnumInterface<ConfigEnumSiteDoFs>(),
      EnumInterface<SuperConfigEnum>()
    );

    return emap;
  }

  /// Convert `casm enum` CLI input to JSON
  ///
  /// All are optionally present, if present on command line
  /// \code
  /// {
  ///   "desc": <array of string, list of enumeration methods to get help for>,
  ///   "help": <bool, print/return help>,
  ///   "method": <string, name of enumeration method>,
  ///   "min": <int, min supercell volume to do enumerations>,
  ///   "max": <int, max supercell volume to do enumerations>,
  ///   "filter": <array of string, filter query/select expression, save configs if evaluates true>,
  ///   "all": <bool, use all existing supercells for enumeration>,
  ///   "verbosity": <string, to be read by Log::verbosity_level>,
  ///   "settings": <string, represents path to settings JSON file>,
  ///   "input": <string, a JSON string>,
  ///   "scelnames": <array of string, list of supercell names, context dependent usage>,
  ///   "confignames": <array of string, list of config names, context dependent usage>,
  ///   "dry-run": <bool, print/return method results but do not save results>
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
    if(vm.count("settings")) {
      json["settings"] = enum_opt.settings_path().string(); //str
    }
    if(vm.count("input")) {
      json["input"] = enum_opt.input_str(); //str
    }
    if(vm.count("scelnames")) {
      json["scelnames"] = enum_opt.supercell_strs(); //vector<std::string>
    }
    if(vm.count("confignames")) {
      json["confignames"] = enum_opt.config_strs(); //vector<std::string>
    }
    if(vm.count("dry-run")) {
      json["dry-run"] = enum_opt.dry_run(); //bool
    }
    return json;
  }

  /// \brief Standardizes parsing casm enum input options to make ScelEnum JSON input
  jsonParser make_enumerator_scel_enum_input(
    jsonParser kwargs,
    const Completer::EnumOptionBase &enum_opt) {

    if(!kwargs.is_obj()) {
      kwargs = jsonParser::object();
    }

    jsonParser scel_input;
    kwargs.get_if(scel_input, "supercells");

    // check supercell shortcuts
    if(enum_opt.vm().count("min")) {
      scel_input["min"] = enum_opt.min_volume();
    }

    if(enum_opt.vm().count("max")) {
      scel_input["max"] = enum_opt.max_volume();
      scel_input["existing_only"] = false;
    }

    if(enum_opt.all_existing()) {
      scel_input.erase("min");
      scel_input.erase("max");
      scel_input["existing_only"] = true;
    }
    else {
      std::vector<std::string> scelnames;
      kwargs.get_if(scelnames, "scelnames");
      if(enum_opt.vm().count("scelnames")) {
        for(std::string const &scelname : enum_opt.supercell_strs())
          scelnames.push_back(scelname);
      }
      if(!scelnames.empty() && !scel_input.contains("names"))
        scel_input["names"].put_array();

      for(std::string const &scelname : scelnames)
        scel_input["names"].push_back(scelname);
    }
    if(scel_input.size() > 0 && !scel_input.contains("existing_only"))
      scel_input["existing_only"] = true;

    return scel_input;
  }

  /// \brief Standardizes parsing casm enum input options to make an ScelEnumProps
  xtal::ScelEnumProps make_enumerator_scel_enum_props(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOptionBase &enum_opt) {

    return make_scel_enum_props(
             primclex,
             make_enumerator_scel_enum_input(_kwargs, enum_opt));
  }

  /// \brief Standardizes parsing casm enum input options to make a SuperlatticeEnumerator
  ///
  /// See SuperConfigEnum for example documentation
  std::unique_ptr<xtal::SuperlatticeEnumerator> make_enumerator_superlat_enum(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOptionBase &enum_opt) {

    xtal::ScelEnumProps enum_props = make_enumerator_scel_enum_props(
                                       primclex,
                                       _kwargs,
                                       enum_opt);

    Supercell unit_cell(&primclex, enum_props.generating_matrix());

    auto fg = unit_cell.factor_group();
    return notstd::make_unique<xtal::SuperlatticeEnumerator>(
             fg.begin(),
             fg.end(),
             primclex.prim().lattice(),
             enum_props);
  }

  /// \brief Standardizes parsing casm enum input options to make an ScelEnum
  ///
  /// See ConfigEnumAllOccupations for example documentation
  std::unique_ptr<ScelEnum> make_enumerator_scel_enum(
    const PrimClex &primclex,
    const jsonParser &_kwargs,
    const Completer::EnumOptionBase &enum_opt) {

    return notstd::make_unique<ScelEnum>(
             primclex,
             make_enumerator_scel_enum_input(_kwargs, enum_opt));
  }

  std::vector<ConfigEnumInput> make_enumerator_input_configs(
    PrimClex const &primclex,
    jsonParser const &_kwargs,
    Completer::EnumOptionBase const &enum_opt,
    EnumeratorMap const *interface_map) {

    jsonParser json {_kwargs};
    InputParser<std::vector<ConfigEnumInput>> parser {
      json,
      primclex.shared_prim(),
      primclex.db<Supercell>(),
      primclex.db<Configuration>()};

    std::runtime_error error_if_invalid {"Error reading enumerator input from JSON"};
    report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

    return *parser.value;
  }

  /// \brief Standardizes parsing casm enum filter expressions
  ///
  /// See ConfigEnumAllOccupations for example documentation
  std::vector<std::string> make_enumerator_filter_expr(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {

    std::vector<std::string> filter_expr;
    // check shortcuts
    if(enum_opt.vm().count("filter")) {
      filter_expr = enum_opt.filter_strs();
    }
    else if(_kwargs.contains("filter")) {
      filter_expr.push_back(_kwargs["filter"].get<std::string>());
    };

    return filter_expr;
  }

  /// \brief Get dry-run value (default=false)
  bool dry_run(
    const jsonParser &_kwargs,
    const Completer::EnumOption &enum_opt) {
    if(enum_opt.vm().count("dry-run")) {
      return true;
    }
    else {
      return _kwargs.get_if_else<bool>("dry_run", false);
    }
  }
}
