#include "casm/app/casm_functions.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/database/Import.hh"
#include "casm/database/Selection.hh"

#include "casm/completer/Handlers.hh"

namespace CASM {
  namespace Completer {
    UpdateOption::UpdateOption(): OptionHandlerBase("update") {}

    void UpdateOption::initialize() {
      add_help_suboption();

      fs::path _default = "ALL";
      add_configlist_suboption(_default);

      add_configtype_suboption(
        QueryTraits<Configuration>::short_name, config_types_short());

      bool required = false;
      add_settings_suboption(required);
      add_input_suboption(required);

      return;
    }

  }

  void print_names(std::ostream &sout, const UpdaterMap &updaters) {
    sout << "The update type options are:\n\n";

    for(const auto &f : updaters) {
      sout << "  " << f.name() << std::endl;
    }
  }

  // ///////////////////////////////////////
  // 'update' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  /// Update proceeds in two steps.
  ///   1) For each selected configuration for which properties.calc.json exists:
  ///       - read properties.calc.json file
  ///       - map it onto a Configuration of the PrimClex
  ///       - record relaxation data (lattice & basis deformation cost)
  ///       - clear existing properties from database
  ///
  ///   2) Iterate over each update record and do the following:
  ///       - Store all initial configuration -> relaxed configuration mappings,
  ///         with both the initial and relaxed configuration
  ///       - For all relaxed configurations, determine which properties to use:
  ///         - if self mapped (initial config == relaxed config) and calculated,
  ///         use those calculation results;
  ///         - else determine which configuration is the best mapping to the
  ///         relaxed configuration (if any) and use those calculation results
  ///
  int update_command(const CommandArgs &args) {

    Completer::UpdateOption update_opt;
    po::variables_map &vm = update_opt.vm();

    try {
      po::store(po::parse_command_line(args.argc, args.argv, update_opt.desc()), vm);


      /** --help option
       */
      if(vm.count("help")) {
        args.log << std::endl;
        args.log << update_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << update_opt.desc() << std::endl;
        args.log << "DESCRIPTION" << std::endl;
        args.log << "    Updates all configuration properties from training data\n";
        args.log << "\n";

        return 0;
      }

      po::notify(vm);
    }
    catch(po::error &e) {
      args.err_log << update_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return 1;
    }

    catch(std::exception &e) {
      args.err_log << update_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return 1;
    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    jsonParser input;
    if(vm.count("settings")) {
      input = jsonParser {update_opt.settings_path()};
    }
    else if(vm.count("input")) {
      input = jsonParser::parse(update_opt.input_str());
    }

    std::unique_ptr<UpdaterMap> updaters = make_interface_map<Completer::UpdateOption>();
    updaters->insert(DB::UpdateInterface<Configuration>());

    auto it = updaters->find(update_opt.configtype());
    if(it != updaters->end()) {
      return it->run(primclex, input, update_opt);
    }
    else {
      args.err_log << "No match found for --type " << update_opt.configtype() << std::endl;
      print_names(args.log, *updaters);
      return ERR_INVALID_ARG;
    }
  }
}

