#include "casm/app/casm_functions.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/clex/ConfigSelection.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/crystallography/jsonStruc.hh"

#include "casm/completer/Handlers.hh"

namespace CASM {

  namespace Completer {

    ViewOption::ViewOption(): OptionHandlerBase("view") {}

    void ViewOption::initialize() {
      add_help_suboption();
      add_confignames_suboption();
      add_configlist_nodefault_suboption();
      m_desc.add_options()
      ("relaxed,r", "Attempt to display corresponding relaxed structures to the given configurations.");
      return;
    }
  }

  int view_command(const CommandArgs &args) {

    fs::path selection;
    po::variables_map vm;
    //bool force;
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
        args.log << "\n";
        args.log << view_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log << "\n";
        args.log << view_opt.desc() << std::endl;
        args.log << "This allows opening visualization programs directly from \n"
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
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log << view_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log << "Unhandled Exception reached the top of main: "
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
      args.err_log << "Error in 'casm view': No command set. Use 'casm settings "
                   "--set-view-command' to set the command to open visualization "
                   "software. It should take one argument, the path to a POSCAR "
                   "to be visualized. For example, to use VESTA on Mac: casm settings --set-view-command 'casm.view \"open -a /Applications/VESTA/VESTA.app\"'.\n";
      return ERR_MISSING_DEPENDS;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);

    ConfigSelection<false> config_select;
    if(!vm.count("config")) {
      config_select = ConfigSelection<false>(primclex, "NONE");
    }
    else if(selection == "MASTER") {
      config_select = ConfigSelection<false>(primclex);
    }
    else {
      config_select = ConfigSelection<false>(primclex, selection);
    }

    // add --confignames (or positional) input
    for(int i = 0; i < confignames.size(); i++) {
      config_select.set_selected(confignames[i], true);
    }

    fs::path tmp_dir = root / ".casm" / "tmp";
    fs::create_directory(tmp_dir);

    // execute the 'casm view' command for each selected configuration
    for(auto it = config_select.selected_config_cbegin(); it != config_select.selected_config_cend(); ++it) {
      //for relaxed structure
      if(vm.count("relaxed") && is_calculated(*it)) {
        fs::ofstream file;
        fs::path POSCARpath = tmp_dir / "POSCAR";

        //Make pos_path the path to properties.calc.json constructed from it
        fs::path pos_path = it->calc_properties_path();
        args.log << "Obtaining relaxed structure from:\n";
        args.log << pos_path.string() << std::endl;
        BasicStructure<Site> import_struc;
        jsonParser datajson(pos_path);
        from_json(simple_json(import_struc, "relaxed_"), datajson);

        file.open(POSCARpath);
        VaspIO::PrintPOSCAR p(import_struc);
        p.sort();
        p.print(file);
        file.close();

        args.log << it->name() << " relaxed" << ":\n";
        Popen popen;
        popen.popen(set.view_command() + " " + POSCARpath.string());
        popen.print(std::cout);
      }

      else {
        // write '.casm/tmp/POSCAR'
        fs::ofstream file;
        fs::path POSCARpath = tmp_dir / "POSCAR";
        file.open(POSCARpath);
        VaspIO::PrintPOSCAR p(*it);
        p.sort();
        p.print(file);
        file.close();
        if(vm.count("relaxed")) {
          args.log << "No relaxed structure found." << "\n";
        }
        args.log << it->name() << ":\n";
        Popen popen;
        popen.popen(set.view_command() + " " + POSCARpath.string());
        popen.print(args.log);
      }

    }

    return 0;
  };

}
