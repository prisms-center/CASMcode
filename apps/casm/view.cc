#include "view.hh"

#include <string>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"

namespace CASM {

  int view_command(int argc, char *argv[]) {

    fs::path selection;
    po::variables_map vm;
    bool force;
    std::vector<std::string> configname;
    
    

    try {

      // Set command line options using boost program_options
      po::options_description desc("'casm view' usage");
      desc.add_options()
      ("help,h", "Print help message")
      
      ("configname",
        po::value<std::vector<std::string> >(&configname)->multitoken(),
        "1 or more configname to view.")
      
      ("config,c", 
        po::value<fs::path>(&selection), 
        "Selected configurations to view.")
        
      ("force,f", "Overrwrite output file");
      
      
      // allow configname as positional options
      po::positional_options_description p;
      p.add("configname", -1);

      try {
        po::store(po::command_line_parser(argc-1, argv+1).options(desc).positional(p).run(), vm);
        //po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "This prepares input files for eci_search. In the specified         \n" <<
                    "configuration list all selected configurations are used for printing\n" <<
                    "an 'energy' file and 'corr.in' file that form the traning data for \n" <<
                    "eci fitting. \n\n";

          return 0;
        }

        po::notify(vm); // throws on error, so do after help in case
        // there are any problems

      }
      catch(po::error &e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << desc << std::endl;
        return ERR_INVALID_ARG;
      }
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }

    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cerr << "Error in 'casm view': No casm project found." << std::endl;
      return ERR_NO_PROJ;
    }

    std::cout << "\n***************************\n" << std::endl;


    DirectoryStructure dir(root);
    ProjectSettings set(root);
    if(set.view_command().empty()) {
      std::cerr << "Error in 'casm view': No command set. Use 'casm settings "
                   "--set-view-command' to set the command to open visualization "
                   "software. It should take one argument, the path to a POSCAR "
                   "to be visualized. For example, to use VESTA on Mac: 'open -a /Applications/VESTA/VESTA.app'.\n";
      return ERR_MISSING_DEPENDS;
    }
    
    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    std::cout << "Reading config list" << std::endl;
    ConfigSelection<false> config_select;
    if(!vm.count("config")) {
      config_select = ConfigSelection<false>::None(primclex);
    }
    else if(selection == "MASTER") {
      config_select = ConfigSelection<false>(primclex);
    }
    else {
      config_select = ConfigSelection<false>(primclex, selection);
    }
    
    std::cout << "Adding specific configurations" << std::endl;
    // add --configname (or positional) input
    for(int i=0; i<configname.size(); i++) {
      config_select.set_selected(configname[i], true);
    }
    
    fs::path tmp_dir = root / ".casm" / "tmp";
    fs::create_directory(tmp_dir);

    std::cout << "Loop over chosen configurations" << std::endl;
    // execute the 'casm view' command for each selected configuration
    for(auto it=config_select.selected_config_cbegin(); it!=config_select.selected_config_cend(); ++it) {
      // write '.casm/tmp/POSCAR'
      fs::ofstream file;
      fs::path POSCARpath = tmp_dir / "POSCAR";
      file.open(POSCARpath);
      it->get_supercell().print(*it, file, FRAC);
      file.close();
      
      std::cout << it->name() << ":\n";
      Popen p;
      p.popen(set.view_command() + " " + POSCARpath.string());
      p.print(std::cout);
      
    }

    return 0;
  };

}
