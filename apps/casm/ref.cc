#include "ref.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"
#include "casm/casm_io/json_io/clex.hh"

namespace CASM {
  
  namespace ref_impl {
  
    int initialize_global(fs::path chem_ref_path, const PrimClex& primclex, std::string set_str, double lin_alg_tol) {
    
      try {
        ChemicalReference chem_ref = 
          read_chemical_reference(jsonParser::parse(set_str), primclex.get_prim(), lin_alg_tol);
        
        std::cout << "Initializing the chemical reference to: \n";
        ChemicalReferencePrinter p(std::cout, chem_ref);
        p.print_all();
        write_chemical_reference(chem_ref, chem_ref_path);
        return 0;
      }
      catch(std::exception& e) {
        std::cerr << "Error parsing input for 'casm ref --set ' with: " 
                  << set_str << "\n";
        std::cerr << e.what() << std::endl;
        return ERR_INVALID_ARG;
      }
    }
    
    int update_global(fs::path chem_ref_path, const PrimClex& primclex, std::string set_str, double lin_alg_tol) {
            
      try {
        ChemicalReference chem_ref = 
          read_chemical_reference(chem_ref_path, primclex.get_prim(), lin_alg_tol);
        
        auto input = one_chemical_reference_from_json(primclex.get_prim(), jsonParser::parse(set_str));
        if(input.second.empty()) {
          chem_ref.set_global(input.first);
        }
        else {
          chem_ref.set_global(input.second.begin(), input.second.end(), lin_alg_tol);
        }
        
        std::cout << "Updating the project-wide chemical reference to: \n";
        ChemicalReferencePrinter p(std::cout, chem_ref);
        p.print_global();
        write_chemical_reference(chem_ref, chem_ref_path);
        return 0;
      }
      catch(std::exception& e) {
        std::cerr << "Error parsing input for 'casm ref --set ' with: " 
                  << set_str << "\n";
        std::cerr << e.what() << std::endl;
        return ERR_INVALID_ARG;
      }
    }
    
    int update_config(std::string configname, 
                      fs::path chem_ref_path, 
                      const PrimClex& primclex, 
                      std::string set_str, 
                      double lin_alg_tol) {
    
      if(!fs::exists(chem_ref_path)) {
        std::cerr << "Error using 'casm ref --set --configname': No reference found.\n";
        std::cerr << "  Expected file at: " << chem_ref_path << "\n";
        std::cerr << "Use 'casm ref --set' or 'casm ref --set-auto' to set a project-wide reference first.\n";
        return ERR_MISSING_INPUT_FILE;
      }
      
      ChemicalReference chem_ref = read_chemical_reference(chem_ref_path, primclex.get_prim(), lin_alg_tol);
      
      try {
        const Configuration& config = primclex.configuration(configname);
      }
      catch(...) {
        std::cerr << "Error using 'casm ref --set --configname': \n"
                     "  Could not find configuration with name: " << configname << "\n";
        return ERR_INVALID_ARG;
      }
      
      try {
        auto input = one_chemical_reference_from_json(primclex.get_prim(), jsonParser::parse(set_str));
        if(input.second.empty()) {
          chem_ref.set_config(configname, input.first);
        }
        else {
          chem_ref.set_config(configname, input.second.begin(), input.second.end(), lin_alg_tol);
        }
        
        std::cout << "Updating the " << configname << " specialized reference to: \n";
        ChemicalReferencePrinter p(std::cout, chem_ref);
        p.print_config(configname);
        write_chemical_reference(chem_ref, chem_ref_path);
        return 0;
      }
      catch(std::exception& e) {
        std::cerr << "Error parsing input for 'casm ref --configname --set ' with: " 
                  << set_str << "\n";
        std::cerr << e.what() << std::endl;
        return ERR_INVALID_ARG;
      }
    }
    
    int update_supercell(std::string scelname, 
                        fs::path chem_ref_path, 
                        const PrimClex& primclex, 
                        std::string set_str, 
                        double lin_alg_tol) {
    
      if(!fs::exists(chem_ref_path)) {
        std::cerr << "Error using 'casm ref --set --scelname': No reference found.\n";
        std::cerr << "  Expected file at: " << chem_ref_path << "\n";
        std::cerr << "Use 'casm ref --set' or 'casm ref --set-auto' to set a project-wide reference first.\n";
        return ERR_MISSING_INPUT_FILE;
      }
      
      ChemicalReference chem_ref = read_chemical_reference(chem_ref_path, primclex.get_prim(), lin_alg_tol);
      
      try {
        const Supercell& scel = primclex.get_supercell(scelname);
      }
      catch(...) {
        std::cerr << "Error using 'casm ref --set --scelname': \n"
                     "  Could not find supercell with name: " << scelname << "\n";
        return ERR_INVALID_ARG;
      }
      
      try {
      
        auto input = one_chemical_reference_from_json(primclex.get_prim(), jsonParser::parse(set_str));
        if(input.second.empty()) {
          chem_ref.set_supercell(scelname, input.first);
        }
        else {
          chem_ref.set_supercell(scelname, input.second.begin(), input.second.end(), lin_alg_tol);
        }
        
        std::cout << "Updating the " << scelname << " specialized reference to: \n";
        ChemicalReferencePrinter p(std::cout, chem_ref);
        p.print_supercell(scelname);
        write_chemical_reference(chem_ref, chem_ref_path);
        return 0;
      }
      catch(std::exception& e) {
        std::cerr << "Error parsing input for 'casm ref --scelname --set ' with: " 
                  << set_str << "\n";
        std::cerr << e.what() << std::endl;
        return ERR_INVALID_ARG;
      }
    }
  
  }

  // ///////////////////////////////////////
  // 'ref' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int ref_command(int argc, char *argv[]) {

    po::variables_map vm;
    int choice;
    std::string scelname, configname, set_str;
    int refid, configid;
    double lin_alg_tol = 1e-14;
    
    try {
      po::options_description desc("'casm ref' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      //("composition-space", "Display information on current composition space")
      ("display,d", "Display current reference states")
      ("set-auto", "Automatically set project level reference states using DFT results")
      ("set", po::value<std::string>(&set_str), 
              "Set reference states using user specified compositions and energies "
              "(Default: set project-wide references). \n"
              "See examples below for the form of expected input.") 
      ("scelname", po::value<std::string>(&scelname), 
                   "Use references given via --set for a particular supercell only")
      ("configname", po::value<std::string>(&configname),
                     "Use references given via --set for a particular configuration only");
      
      
      try {
        po::store(po::parse_command_line(argc, argv, desc), vm);

        bool call_help = false;

        //quit out if there are no arguments
        if(!vm.count("help")) {
          if(vm.count("set") + vm.count("display") != 1) {
            std::cout << "Error in 'casm ref'. Please select one of --display \n";
            std::cout << "or --set to use this option." << std::endl;

            call_help = true;
          }

          if(vm.count("set")) {
            if(vm.count("set-auto")) {
              std::cerr << "Error: Please select only one of --set, --set-auto \n";
              
              call_help = true;
            }
            if(vm.count("scelname") + vm.count("configname") > 1) {
              std::cerr << "Error in 'casm ref --set'. Please select only one of --scelname, --configname \n";
              
              call_help = true;
            }
          }
        }

        /** --help option
         */
        if(vm.count("help") || call_help) {
          std::cout << std::endl;
          std::cout << desc << std::endl;

          std::cout << "DESCRIPTION" << std::endl;
          std::cout << "    The chemical reference determines the value of the formation energy  \n"
                       "    and chemical potentials calculated by CASM.                          \n\n"
                       
                       "    Chemical references states are set by specifying a hyperplane in     \n"
                       "    energy/atom - composition (as atom_frac) space. This may be done by  \n"
                       "    specifying the hyperplane explicitly, or by specifying several       \n"
                       "    reference states with energy/atom and composition (as atom_frac) for \n"
                       "    enough states to span the composition space of the allowed occupants \n"
                       "    specified in the prim. For consistency with other CASM projects,     \n"
                       "    additional reference states extending to other compositional         \n"
                       "    dimensions may be included also.                                     \n\n"
                       
                       "    It is also possible to specialize the chemical reference at the      \n"
                       "    supercell or configuration level.                                    \n\n";

          std::cout << "    Examples:\n";
          //std::cout << "      casm ref --composition-space \n";
          //std::cout << "      - Print composition space column matrix of the primitive\n";
          //std::cout << "      - Print null space column matrix\n";
          //std::cout << "\n";
          std::cout << "      casm ref --display \n";
          std::cout << "      - Print chemical reference\n";
          std::cout << "\n";
          std::cout << "      casm ref --set-auto\n";
          std::cout << "      - set all reference states using DFT results for configurations with\n";
          std::cout << "        extreme compositions.\n";
          std::cout << "      - set reference for compositions outside range of this project to 0.0\n";
          std::cout << "\n";
          std::cout << "      casm ref --set '[{\"Zr\":1, \"energy_per_species\":-8.546979385}, {\"Zr\":1, \"O\":1, \"energy_per_species\":-9.090697345}]'";
          std::cout << "      - set Zr and ZrO, with given energy per species, as reference states\n";
          std::cout << "\n";
          std::cout << "      casm ref --scelname SCEL3_3_1_1_0_2_2 --set '[{\"Zr\":1, \"energy_per_species\":-8.546979385}, {\"Zr\":1, \"O\":1, \"energy_per_species\":-9.090697345}]'";
          std::cout << "      - set reference states as specified for configurations in supercell SCEL3_3_1_1_0_2_2\n";
          std::cout << "\n";
          std::cout << "      casm ref --configname SCEL3_3_1_1_0_2_2/2 --set '[{\"Zr\":1, \"energy_per_species\":-8.546979385}, {\"Zr\":1, \"O\":1, \"energy_per_species\":-9.090697345}]'";
          std::cout << "      - set reference states as specified for configuration SCEL3_3_1_1_0_2_2/2\n";
          std::cout << "\n";
          

          if(call_help)
            return ERR_INVALID_ARG;

          return 0;
        }

        po::notify(vm);

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
      std::cout << "Error in 'casm ref': No casm project found." << std::endl;
      return ERR_NO_PROJ;
    }
    fs::current_path(root);

    std::cout << "\n***************************\n" << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    std::string calctype = primclex.settings().calctype();
    std::string ref = primclex.settings().ref();
    fs::path chem_ref_path = primclex.dir().chemical_reference(calctype, ref);
      
    if(vm.count("display")) {
      if(!fs::exists(chem_ref_path)) {
        std::cerr << "Error using 'casm ref --display': No reference found.\n";
        std::cerr << "  Expected file at: " << chem_ref_path << "\n";
        std::cerr << "Use 'casm ref --set' or 'casm ref --set-auto' to set a reference\n";
        return ERR_MISSING_INPUT_FILE;
      }
      
      ChemicalReference chem_ref = read_chemical_reference(chem_ref_path, primclex.get_prim(), lin_alg_tol);
      ChemicalReferencePrinter p(std::cout, chem_ref);
      p.print_all();
      
    }
    else if(vm.count("set-auto")) {
      try {
        std::cout << "  Set reference states automatically.\n\n" << std::endl;
        ChemicalReference chem_ref = auto_chemical_reference(primclex, lin_alg_tol);
        ChemicalReferencePrinter p(std::cout, chem_ref);
        p.print_all();
        write_chemical_reference(chem_ref, chem_ref_path);
        return 0;
      }
      catch(std::exception& e) {
        std::cerr << "Error setting reference states automatically.\n\n";
        std::cerr << e.what() << std::endl;
        return ERR_UNKNOWN;
      }
    }
    else if(vm.count("set")) {
      
      using namespace ref_impl;
      
      // --- Set project-wide ref
      
      if(!(vm.count("scelname") + vm.count("configname"))) {
        
        if(!fs::exists(chem_ref_path)) {
          // -- Initializing ref
          return initialize_global(chem_ref_path, primclex, set_str, lin_alg_tol);
        }
        else {
          // -- Updating project-wide ref
          return update_global(chem_ref_path, primclex, set_str, lin_alg_tol);
        }
      }
      
      // --- Set config specific ref
      
      else if(vm.count("configname")) {
        
        return update_config(configname, chem_ref_path, primclex, set_str, lin_alg_tol);
      }
      
      // --- Set supercell specific ref
      
      else {
        
        return update_supercell(configname, chem_ref_path, primclex, set_str, lin_alg_tol);
        
      }
      
    }
    
    return 0;
  }

}

