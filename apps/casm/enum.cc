#include "enum.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/FilteredConfigIterator.hh"

namespace CASM {


  // ///////////////////////////////////////
  // 'enum' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int enum_command(int argc, char *argv[]) {

    //casm enum [—supercell min max] [—config supercell ] [—hopconfigs hop.background]
    //- enumerate supercells and configs and hop local configurations

    int min_vol = 1, max_vol;
    std::vector<std::string> scellname_list, filter_expr;
    //double tol;
    COORD_TYPE coordtype = CASM::CART;
    po::variables_map vm;


    /// Set command line options using boost program_options
    po::options_description desc("'casm enum' usage");
    desc.add_options()
    ("help,h", "Write help documentation")
    ("min", po::value<int>(&min_vol), "Min volume")
    ("max", po::value<int>(&max_vol), "Max volume")
    ("filter", po::value<std::vector<std::string> >(&filter_expr)->multitoken(), "Filter configuration enumeration so that only configurations matching a 'casm query'-type expression are recorded")
    ("scellname,n", po::value<std::vector<std::string> >(&scellname_list)->multitoken(), "Enumerate configs for given supercells")
    ("all,a", "Enumerate configurations for all supercells")
    ("supercells,s", "Enumerate supercells")
    ("configs,c", "Enumerate configurations");

    // currently unused...
    //("tol", po::value<double>(&tol)->default_value(CASM::TOL), "Tolerance used for checking symmetry")
    //("coord", po::value<COORD_TYPE>(&coordtype)->default_value(CASM::CART), "Coord mode: FRAC=0, or (default) CART=1");

    try {
      po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

      /** --help option
       */
      if(vm.count("help")) {
        std::cout << "\n";
        std::cout << desc << std::endl;


        std::cout << "DESCRIPTION" << std::endl;
        std::cout << "    Enumerate supercells and configurations\n";
        std::cout << "    - expects a PRIM file in the project root directory \n";
        std::cout << "    - if --min is given, then --max must be given \n";


        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      if(vm.count("min") && !vm.count("max")) {
        std::cerr << "\n" << desc << "\n" << std::endl;
        std::cerr << "Error in 'casm enum'. If --min is given, --max must also be given." << std::endl;
        return ERR_INVALID_ARG;
      }
      if(vm.count("supercells") + vm.count("configs") != 1) {
        std::cerr << "\n" << desc << "\n" << std::endl;
        std::cerr << "Error in 'casm enum'. Exactly one of either --supercells or --configs must be given." << std::endl;
        return ERR_INVALID_ARG;
      }
      if(vm.count("supercells") && !vm.count("max")) {
        std::cerr << "\n" << desc << "\n" << std::endl;
        std::cerr << "Error in 'casm enum'. If --supercells is given, --max must be given." << std::endl;
        return ERR_INVALID_ARG;
      }
      if(vm.count("configs") && (vm.count("max") + vm.count("all") != 1)) {
        std::cerr << "\n" << desc << "\n" << std::endl;
        std::cerr << "Error in 'casm enum'. If --configs is given, exactly one of either --max or --all must be given." << std::endl;
        return ERR_INVALID_ARG;
      }
    }
    catch(po::error &e) {
      std::cerr << desc << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      std::cerr << desc << std::endl;
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      return ERR_UNKNOWN;

    }

    COORD_MODE C(coordtype);

    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cerr << "Error in 'casm enum': No casm project found." << std::endl;
      return ERR_NO_PROJ;
    }

    std::cout << "\n***************************\n" << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    const DirectoryStructure &dir = primclex.dir();
    const ProjectSettings &set = primclex.settings();
    std::cout << "  DONE." << std::endl << std::endl;

    if(vm.count("supercells")) {
      std::cout << "\n***************************\n" << std::endl;

      std::cout << "Generating supercells from " << min_vol << " to " << max_vol << std::endl << std::endl;
      primclex.generate_supercells(min_vol, max_vol, true);
      std::cout << "\n  DONE." << std::endl << std::endl;

      std::cout << "Write SCEL." << std::endl << std::endl;
      primclex.print_supercells();

    }
    else if(vm.count("configs")) {
      std::cout << "\n***************************\n" << std::endl;
      std::vector<Supercell *> scel_selection;

      //Build the supercell selection based on user input
      if(vm.count("all")) {
        std::cout << "Enumerate all configurations" << std::endl << std::endl;
        for(int j = 0; j < primclex.get_supercell_list().size(); j++)
          scel_selection.push_back(&(primclex.get_supercell(j)));
      }
      else {
        if(vm.count("max")) {
          std::cout << "Enumerate configurations from volume " << min_vol << " to " << max_vol << std::endl << std::endl;
          for(int j = 0; j < primclex.get_supercell_list().size(); j++) {
            if(primclex.get_supercell(j).volume() >= min_vol && primclex.get_supercell(j).volume() <= max_vol)
              scel_selection.push_back(&(primclex.get_supercell(j)));
          }
        }
        if(vm.count("scellname")) {
          Index j;
          std::cout << "Enumerate configurations for named supercells" << std::endl << std::endl;
          for(int i = 0; i < scellname_list.size(); i++) {
            if(!primclex.contains_supercell(scellname_list[i], j)) {
              std::cout << "Error in 'casm enum'. Did not find supercell: " << scellname_list[i] << std::endl;
              return 1;
            }
            scel_selection.push_back(&(primclex.get_supercell(j)));
          }
        }
      }
      //\Finished building supercell selection

      if(scel_selection.empty()) {
        std::cout << "Did not find any supercells. Make sure to 'casm enum --supercells' first!" << std::endl << std::endl;
        return ERR_MISSING_DEPENDS;
      }

      //We have the selection. Now do enumeration
      for(auto it = scel_selection.begin(); it != scel_selection.end(); ++it) {
        std::cout << "  Enumerate configurations for " << (**it).get_name() << " ...  " << std::flush;

        ConfigEnumAllOccupations<Configuration> enumerator(**it);
        Index num_before = (**it).get_config_list().size();
        if(vm.count("filter")) {
          try {
            (**it).add_unique_canon_configs(filter_begin(enumerator.begin(), enumerator.end(), filter_expr, set.config_io()), filter_end(enumerator.end()));
          }
          catch(std::exception &e) {
            std::cerr << "Cannot filter configurations using the expression provided: \n" << e.what() << "\nExiting...\n";
            return ERR_INVALID_ARG;
          }
        }
        else
          (**it).add_unique_canon_configs(enumerator.begin(), enumerator.end());

        std::cout << ((**it).get_config_list().size() - num_before) << " configs." << std::endl;
      }
      std::cout << "  DONE." << std::endl << std::endl;

    }

    std::cout << "Writing config_list..." << std::endl;
    primclex.write_config_list();
    std::cout << "  DONE" << std::endl;

    std::cout << std::endl;

    return 0;
  };

}

