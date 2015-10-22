#include "perturb.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/clex/ConfigEnumStrain.hh"

namespace CASM {


  // ///////////////////////////////////////
  // 'perturb' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int perturb_command(int argc, char *argv[]) {

    double tol = CASM::TOL;
    bool is_trans = false;
    fs::path cspecs_path, abs_cspecs_path;
    fs::path selection;
    COORD_TYPE coordtype;
    po::variables_map vm;
    std::vector<double> init, final;
    double inc_value;
    std::string strain_mode;
    try {

      /// Set command line options using boost program_options
      po::options_description desc("'casm perturb' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      ("cspecs", po::value<fs::path>(&cspecs_path), "Cluster specifications file defining perturbation")
      ("config,c", po::value<fs::path>(&selection),
       "Selected configurations are used reference for generating perturbations. If not specified, or 'MASTER' given, uses master list selection.")
        ("strain,s","Generate strain perturbations")
        ("init",po::value<std::vector<double> > (&init),"Initial vector")
        ("final",po::value<std::vector<double> > (&final),"Final vector")
        ("inc",po::value<double> (&inc_value)->default_value(0.02),"Strain increment")
        ("inc",po::value<std::string> (&strain_mode)->default_value("GL"),"Strain mode name");

        
      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "DESCRIPTION" << std::endl;
          std::cout << "    Generate supercells that are perturbations of a reference\n";
          std::cout << "    configuration.                                           \n";
          std::cout << "    - using the --cspecs option, a bspecs.json type file is  \n";
          std::cout << "      required to determine the extent of the perturbations. \n";
          std::cout << "      Currently only 'orbit_branch_specs' are supported.     \n";
          std::cout << "    - perturbations are generated about selected reference   \n";
          std::cout << "      configurations                                         \n";
          std::cout << std::endl;

          return 0;
        }

        po::notify(vm); // throws on error, so do after help in case
        // there are any problems

      }
      catch(po::error &e) {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << desc << std::endl;
        return 1;
      }
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return 1;

    }

    COORD_MODE C(coordtype);

    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cout << "Error in 'casm perturb': No casm project found." << std::endl;
      return 1;
    }
    fs::current_path(root);


    std::cout << "\n***************************\n" << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    DirectoryStructure dir(root);
    ProjectSettings set(root);

    if(vm.count("strain")){
      Eigen::MatrixXd tproj(2,6);
      tproj << 1.0/sqrt(2.0), -1.0/sqrt(6.0),
        -1.0/sqrt(2.0), -1.0/sqrt(6.0),
        0.0, 2.0/sqrt(6.0),
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0;
      
      Eigen::VectorXd _init_vec(2), _final_vec(2);
      _init_vec << init[0], init[1];
      _final_vec << final[0], final[1];
      
      ConfigSelection<false> config_select;
      if(!vm.count("config") || selection == "MASTER") {
        config_select = ConfigSelection<false>(primclex);
      }
      else {
        config_select = ConfigSelection<false>(primclex, selection);
      }
      
      std::cout << "\n***************************\n" << std::endl;
      
      std::cout << "Generating perturbations about configurations " << std::endl << std::endl;
      
      bool verbose = false;
      bool print = true;
      for(auto it = config_select.selected_config_begin(); it != config_select.selected_config_end(); ++it) {
        ConfigEnumStrain<Configuration> enumerator(it->get_supercell(),*it,tproj, _init_vec, _final_vec, inc_value,strain_mode);
        (it->get_supercell()).add_unique_canon_configs(enumerator.begin(), enumerator.end());        
      }
      return 0;
    }
    // want absolute paths
    abs_cspecs_path = fs::absolute(cspecs_path);



    ConfigSelection<false> config_select;
    if(!vm.count("config") || selection == "MASTER") {
      config_select = ConfigSelection<false>(primclex);
    }
    else {
      config_select = ConfigSelection<false>(primclex, selection);
    }

    std::cout << "\n***************************\n" << std::endl;

    std::cout << "Generating perturbations about configurations " << std::endl << std::endl;

    bool verbose = false;
    bool print = true;
    for(auto it = config_select.selected_config_begin(); it != config_select.selected_config_end(); ++it) {
      std::cout << "  " << it->get_supercell().get_name() << "/" << it->get_id() << std::endl;
      it->get_supercell().enumerate_perturb_configurations(*it, abs_cspecs_path, tol, verbose, print);
    }

    std::cout << std::endl << "  DONE." << std::endl << std::endl;

    std::cout << "Writing config_list..." << std::endl;
    primclex.write_config_list();
    std::cout << "  DONE" << std::endl;

    std::cout << std::endl;

    return 0;
  };

}

