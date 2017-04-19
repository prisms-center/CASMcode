#include "casm/app/casm_functions.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"
#include "casm/database/Import.hh"

#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {
    ImportOption::ImportOption(): OptionHandlerBase("import") {}

    const std::vector<fs::path> &ImportOption::pos_vec() const {
      return m_pos_vec;
    }

    const fs::path &ImportOption::batch_path() const {
      return m_batch_path;
    }

    void ImportOption::initialize() {
      add_help_suboption();

      m_desc.add_options()

      ("pos,p",
       po::value<std::vector<fs::path> >(&m_pos_vec)->multitoken()->value_name(ArgHandler::path()),
       "Path(s) to structure(s) being imported (multiple allowed, but no "
       "wild-card matching)")

      ("batch,b",
       po::value<fs::path>(&m_batch_path)->value_name(ArgHandler::path()),
       "Path to batch file, which should list one structure file path per line "
       "(can be used in combination with --pos)")

      ("data,d",
       "Attempt to extract calculation data from the enclosing "
       "directory of the structure files, if it is available")

      ("copy-additional-files",
       "Recursively copy other files from the same directory as the properties.calc.json file.");

      add_configtype_suboption(
        QueryTraits<Configuration>::short_name, config_types_short());
      bool required = false;
      add_settings_suboption(required);
      add_input_suboption(required);


      return;
    }

  }

  void print_names(std::ostream &sout, const ImporterMap &importers) {
    sout << "The import type options are:\n\n";

    for(const auto &f : importers) {
      sout << "  " << f.name() << std::endl;
    }
  }

  // ///////////////////////////////////////
  // 'import' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  /// Import proceeds in two steps.
  ///   1) for each file:
  ///       - read Structure
  ///       - map it onto a Configuration of the PrimClex
  ///       - record relaxation data (lattice & basis deformation cost)
  ///
  ///   2) If data import was requested, iterate over each import record and do
  ///      the following:
  ///       - if multiple imported structures map onto a configuration for which
  ///         there is no calculation data, import calculation data from the
  ///         structure with the lowest mapping cost
  ///       - if one or more imported structuress map onto a configuration for
  ///         which calculation data already exist, do not import any new data
  ///       - if data is imported, the corresponding properties.calc.json file is
  ///         copied into the directory of the mapped configuration. A structure
  ///         file, relaxed_structure.vasp is also written to the directory.
  ///       - relaxed_structure.vasp gives the relaxed structure in a setting and
  ///         orientation that matches the generated POS file
  ///
  int import_command(const CommandArgs &args) {

    /// Set command line options using boost program_options
    Completer::ImportOption import_opt;
    po::variables_map &vm = import_opt.vm();

    try {

      po::store(po::parse_command_line(args.argc, args.argv, import_opt.desc()), import_opt.vm()); // can throw

      /** --help option
       */
      if(vm.count("help")) {
        args.log << std::endl;
        args.log << import_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << import_opt.desc() << std::endl;

        args.log << "DESCRIPTION" << std::endl;
        args.log << "    Import structure specified by --pos. If it doesn't exist make a directory for it and copy data over" << std::endl;
        args.log << "    If a *.json file is specified, it will be interpreted as a 'calc.properties.json' file." << std::endl;
        return 0;
      }

      po::notify(import_opt.vm()); // throws on error, so do after help in case
      // there are any problems

    }
    catch(po::error &e) {
      args.err_log << import_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return 3;
    }
    catch(std::exception &e) {
      args.err_log << import_opt.desc() << std::endl;
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      return 4;

    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);


    std::vector<fs::path> pos_paths;
    auto res = DB::Import<Configuration>::construct_pos_paths(primclex, import_opt, std::back_inserter(pos_paths));
    if(res.second) {
      return res.second;
    }

    jsonParser input;
    if(vm.count("settings")) {
      input = jsonParser {import_opt.settings_path()};
    }
    else if(vm.count("input")) {
      input = jsonParser::parse(import_opt.input_str());
    }

    std::unique_ptr<ImporterMap> importers = make_interface_map<Completer::ImportOption>();
    importers->insert(DB::ImportInterface<Configuration>());

    auto it = importers->find(import_opt.configtype());
    if(it != importers->end()) {
      return it->run(primclex, input, import_opt);
    }
    else {
      args.err_log << "No match found for --type " << import_opt.configtype() << std::endl;
      print_names(args.log, *importers);
      return ERR_INVALID_ARG;
    }
  }

}




