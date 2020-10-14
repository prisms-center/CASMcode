#include <boost/filesystem/fstream.hpp>
#include "casm/app/casm_functions.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/system/Popen.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/Configuration_impl.hh"
#include "casm/clex/SimpleStructureTools.hh"

#include "casm/completer/Handlers.hh"
#include "casm/database/DatabaseTypesTraits.hh"

// need to add specializations here
#include "casm/database/ConfigImport.hh"

namespace CASM {

  namespace Completer {

    ViewOption::ViewOption(): OptionHandlerBase("view") {}

    void ViewOption::initialize() {
      add_help_suboption();
      add_confignames_suboption();
      add_configtype_suboption(traits<Configuration>::short_name, DB::config_types_short());
      add_configlist_nodefault_suboption();
      m_desc.add_options()
      ("images,i", po::value<int>(&m_images)->default_value(0), "Number of images between initial and final state when viewing diff_trans_configs.")
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
      po::store(po::command_line_parser(args.argc() - 1, args.argv() + 1).options(view_opt.desc()).positional(p).run(), vm);
      //po::store(po::parse_command_line(argc, argv, view_opt.desc()), vm); // can throw

      /** --help option
      */
      if(vm.count("help")) {
        args.log() << "\n";
        args.log() << view_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log() << "\n";
        args.log() << view_opt.desc() << std::endl;
        args.log() << "This allows opening visualization programs directly from \n"
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
      args.err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log() << view_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log() << "Unhandled Exception reached the top of main: "
                     << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log().error("No casm project found");
      args.err_log() << std::endl;
      return ERR_NO_PROJ;
    }

    ProjectSettings set = open_project_settings(root);
    if(set.view_command().empty()) {
      args.err_log() << "Error in 'casm view': No command set. Use 'casm settings "
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
      if(vm.count("relaxed") && is_calculated(config)) {
        fs::ofstream file;
        fs::path POSCARpath = tmp_dir / "POSCAR";

        //Make pos_path the path to properties.calc.json constructed from it
        fs::path pos_path = calc_properties_path(primclex, config.name());
        args.log() << "Obtaining relaxed structure from:\n";
        args.log() << pos_path.string() << std::endl;
        SimpleStructure import_struc;
        jsonParser datajson(pos_path);
        from_json(import_struc, datajson);

        file.open(POSCARpath);
        VaspIO::PrintPOSCAR p(import_struc);
        p.sort();
        p.print(file);
        file.close();
        args.log() << config.name() << " relaxed" << ":\n";
        Popen popen;
        popen.popen(set.view_command() + " " + POSCARpath.string());
        popen.print(std::cout);
      }
      else {
        // write '.casm/tmp/POSCAR'
        fs::ofstream file;
        fs::path POSCARpath = tmp_dir / "POSCAR";
        file.open(POSCARpath);
        VaspIO::PrintPOSCAR p(make_simple_structure(config));
        p.sort();
        p.print(file);
        file.close();
        if(vm.count("relaxed")) {
          args.log() << "No relaxed structure found." << "\n";
        }
        args.log() << config.name() << ":\n";
        Popen popen;
        popen.popen(set.view_command() + " " + POSCARpath.string());
        popen.print(args.log());
      }
    }

    return 0;
  };

}
