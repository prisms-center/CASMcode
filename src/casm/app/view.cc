#include "casm/app/casm_functions.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/database/Selection.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Configuration.hh"

#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {

    ViewOption::ViewOption(): OptionHandlerBase("view") {}

    void ViewOption::initialize() {
      add_help_suboption();
      add_confignames_suboption();
      add_configlist_nodefault_suboption();

      return;
    }
  }

  int view_command(const CommandArgs &args) {

    fs::path selection;
    po::variables_map vm;
    bool force;
    std::vector<std::string> confignames;



    // Set command line options using boost program_options
    Completer::ViewOption view_opt;

    // allow confignames as positional options
    po::positional_options_description p;
    p.add("confignames", -1);

    try {
      po::store(po::command_line_parser(args.argc - 1, args.argv + 1).options(view_opt.desc()).positional(p).run(), vm);
      //po::store(po::parse_command_line(argc, argv, view_opt.desc()), vm); // can throw

      /** --help option
      */
      if(vm.count("help")) {
        std::cout << "\n";
        std::cout << view_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        std::cout << "\n";
        std::cout << view_opt.desc() << std::endl;
        std::cout << "This allows opening visualization programs directly from \n"
                  "CASM. It iterates over all selected configurations and   \n"
                  "one by one writes a POSCAR and executes                  \n"
                  "   '$VIEW_COMMAND /path/to/POSCAR'                       \n"
                  "where $VIEW_COMMAND is set via 'casm settings --set-view-command'.\n"
                  "A script 'casm.view' is included with can be used to run \n"
                  "a command and then pause 1s, which is useful for opening \n"
                  "POSCARs with VESTA.  An example on Mac might look like:  \n"
                  "  casm settings --set-view-command 'casm.view \"open -a /Applications/VESTA/VESTA.app\"' \n\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      selection = view_opt.selection_path().string();
      confignames = view_opt.config_strs();
    }
    catch(po::error &e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      std::cerr << view_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    DirectoryStructure dir(root);
    ProjectSettings set(root);
    if(set.view_command().empty()) {
      std::cerr << "Error in 'casm view': No command set. Use 'casm settings "
                "--set-view-command' to set the command to open visualization "
                "software. It should take one argument, the path to a POSCAR "
                "to be visualized. For example, to use VESTA on Mac: casm settings --set-view-command 'casm.view \"open -a /Applications/VESTA/VESTA.app\"'.\n";
      return ERR_MISSING_DEPENDS;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    DB::Selection<Configuration> config_select;
    if(!vm.count("config")) {
      config_select = DB::Selection<Configuration>(primclex, "NONE");
    }
    else if(selection == "MASTER") {
      config_select = DB::Selection<Configuration>(primclex);
    }
    else {
      config_select = DB::Selection<Configuration>(primclex, selection);
    }

    // add --confignames (or positional) input
    for(int i = 0; i < confignames.size(); i++) {
      config_select.data()[confignames[i]] = true;
    }

    fs::path tmp_dir = root / ".casm" / "tmp";
    fs::create_directory(tmp_dir);

    // execute the 'casm view' command for each selected configuration
    for(const auto &config : config_select.selected()) {
      // write '.casm/tmp/POSCAR'
      fs::ofstream file;
      fs::path POSCARpath = tmp_dir / "POSCAR";
      file.open(POSCARpath);
      VaspIO::PrintPOSCAR p(config);
      p.sort();
      p.print(file);
      file.close();

      std::cout << config.name() << ":\n";
      Popen popen;
      popen.popen(set.view_command() + " " + POSCARpath.string());
      popen.print(std::cout);

    }

    return 0;
  };

}
