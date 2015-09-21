#include "ref.hh"

#include <cstring>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"

namespace CASM {


  // ///////////////////////////////////////
  // 'ref' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int ref_command(int argc, char *argv[]) {

    po::variables_map vm;
    int choice;
    std::string scellname;
    int refid, configid;

    try {
      po::options_description desc("'casm ref' usage");
      desc.add_options()
      ("help,h", "Write help documentation")
      ("display,d", "Display current reference states")
      ("set", "Set reference state using existing calculated properties")
      ("auto,a", "Set standard reference states")
      ("index,i", po::value< int >(&refid), "Reference state index")
      ("supercell,s", po::value< std::string >(&scellname), "Name of supercell")
      ("configid,c", po::value< int >(&configid), "Configuration ID")
      ("update,u", "Update references based on current reference states");

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
            if(vm.count("index") + vm.count("auto") + vm.count("update") != 1) {
              std::cout << "Error in 'casm ref --set'. Please select one of --auto, --index, \n";
              std::cout << "or --update to use this option." << std::endl;

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
          std::cout << "    For a particular configuration (SCEL, CONFIGID), the following \n";
          std::cout << "    directories are checked for reference state files:\n";
          std::cout << "    1) root/supercells/SCEL/CONFIGID/settings/CURR_CALCTYPE/CURR_REF\n";
          std::cout << "    2) root/supercells/SCEL/settings/CURR_CALCTYPE/CURR_REF\n";
          std::cout << "    3) root/settings/CURR_CALCTYPE/CURR_REF\n";
          std::cout << " \n";
          std::cout << "    This option sets reference states at the project level (#3 above).\n";
          std::cout << "    - Must have used 'casm composition' to select parametric \n";
          std::cout << "      composition axes. \n";
          std::cout << "    - display: prints current reference states \n";
          std::cout << "    - set: set reference state to the calculated properties for \n";
          std::cout << "           configuration CONFIGID of supercell SCEL \n";
          std::cout << "    - set auto: tries to use configurations with extreme param \n";
          std::cout << "                compositions for reference states \n";
          std::cout << "    - set update: updates references based on current reference state files \n";
          std::cout << " \n";
          std::cout << "    Examples:\n";
          std::cout << "      casm ref --display \n";
          std::cout << "      - Print all reference states\n";
          std::cout << "\n";
          std::cout << "      casm ref --set --auto\n";
          std::cout << "      - set all reference states using configurations with\n";
          std::cout << "        extreme param compositions.\n";
          std::cout << "\n";
          std::cout << "      casm ref --set --index 1 --supercell SCEL4_2_2_1_1_0_1 --configid 12\n";
          std::cout << "      - set reference state 1 to be configuration 12 of supercell SCEL4_2_2_1_1_0_1\n";
          std::cout << "\n";
          std::cout << "      casm ref --set --update\n";
          std::cout << "      - updates reference properties based on current reference states files\n";
          std::cout << "      - use this after manually customizing reference states\n";
          std::cout << "\n";

          if(call_help)
            return 1;

          return 0;
        }

        po::notify(vm);

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

    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cout << "Error in 'casm ref': No casm project found." << std::endl;
      return 1;
    }
    fs::current_path(root);

    std::cout << "\n***************************\n" << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    std::string calctype = primclex.settings().calctype();
    std::string ref = primclex.settings().ref();

    auto global_ref = [&](int i) {
      return primclex.dir().ref_state(calctype, ref, i);
    };

    auto supercell_ref = [&](std::string scelname, int i) {
      return primclex.dir().supercell_ref_state(scelname, calctype, ref, i);
    };

    auto configuration_ref = [&](std::string configname, int i) {
      return primclex.dir().configuration_ref_state(configname, calctype, ref, i);
    };

    /// check that references exist:
    if(!primclex.has_composition_axes()) {
      std::cout << "Error in 'casm ref'. Please first use 'casm composition' \n";
      std::cout << "to select your parametric composition axes." << std::endl;
      return 1;
    }

    if(vm.count("display")) {

      std::cout << "\n***************************\n" << std::endl;
      std::cout << "Display reference states:\n\n" << std::endl;

      // Project level reference states
      std::cout << "Project level reference states:" << std::endl;
      for(int i = 0; i < primclex.composition_axes().independent_compositions() + 1; i++) {
        fs::path filepath = global_ref(i);
        if(fs::exists(filepath)) {
          std::cout << "\n Reference: " << filepath << std::endl;
          jsonParser json(filepath);
          std::cout << json << std::endl;
        }
      }
      std::cout << std::endl;

      // Supercell level reference states
      std::cout << "Supercell level reference states:" << std::endl;
      for(int i = 0; i < primclex.get_supercell_list().size(); i++) {

        for(int j = 0; j < primclex.composition_axes().independent_compositions() + 1; j++) {
          fs::path filepath = supercell_ref(primclex.get_supercell(i).get_name(), i);
          if(fs::exists(filepath)) {
            std::cout << "\n Reference: " << filepath << std::endl;
            jsonParser json(filepath);
            std::cout << json << std::endl;
          }
        }
      }
      std::cout << std::endl;

      // Supercell level reference states
      std::cout << "Configuration level reference states:" << std::endl;
      for(int i = 0; i < primclex.get_supercell_list().size(); i++) {
        const Supercell &scel = primclex.get_supercell(i);
        for(int j = 0; j < scel.get_config_list().size(); j++) {
          for(int k = 0; k < primclex.composition_axes().independent_compositions() + 1; k++) {
            fs::path filepath = configuration_ref(scel.get_config(i).name(), i);
            if(fs::exists(filepath)) {
              std::cout << "\n Reference: " << filepath << std::endl;
              jsonParser json(filepath);
              std::cout << json << std::endl;
            }
          }
        }
      }
      std::cout << std::endl;

    }
    else if(vm.count("set")) {

      std::cout << "\n***************************\n" << std::endl;

      if(vm.count("index")) {

        std::cout << "  Set reference state " << refid << "." << std::endl;

        primclex.set_reference_state(refid, primclex.get_supercell(scellname).get_config(configid));

        std::cout << "  Update references... " << std::endl;
        primclex.generate_references();
        std::cout << "    DONE" << std::endl;

        fs::path filepath = primclex.dir().ref_state(calctype, ref, refid);
        std::cout << " Reference: " << filepath << std::endl;
        jsonParser json(filepath);
        std::cout << json << std::endl;

      }
      else if(vm.count("auto")) {

        std::cout << "  Set reference states automatically." << std::endl;

        primclex.set_reference_state_auto();

        for(int i = 0; i < primclex.composition_axes().independent_compositions() + 1; i++) {
          fs::path filepath = global_ref(i);
          std::cout << "\n Reference: " << filepath << std::endl;
          jsonParser json(filepath);
          std::cout << json << std::endl;
        }
      }
      else if(vm.count("update")) {

        std::cout << "  Update references... " << std::endl;
        primclex.generate_references();
        std::cout << "    DONE" << std::endl;

      }
    }

    std::cout << std::endl;

    return 0;
  }

}

