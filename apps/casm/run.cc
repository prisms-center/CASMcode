#include "run.hh"

#include<cstring>
#include<unistd.h>

#include "casm/CASM_classes.hh"
#include "casm_functions.hh"

namespace CASM {

  // ///////////////////////////////////////
  // 'run' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int run_command(int argc, char *argv[]) {
    std::string exec, selection;
    double tol;
    po::variables_map vm;

    /// Set command line options using boost program_options
    po::options_description desc("'casm run' usage");
    desc.add_options()
    ("help,h", "Write help documentation")
    ("write-pos", "Write POS file for each selected configuration before executing the command")
    ("exec,e", po::value<std::string>(&exec)->required(), "Command to execute")
    ("config,c", po::value<std::string>(&selection)->default_value("MASTER"), "Config selection");

    try {
      po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

      /** --help option
       */
      if(vm.count("help")) {
        std::cout << "\n";
        std::cout << desc << std::endl;

        std::cout << "DESCRIPTION\n"
                  << "    Executes the requested command for each selected configuration,\n"
                  << "    with the path to the configuration as an argument.             \n\n"

                  << "    Example: casm run --exec \"vasp.relax\" --write-pos\n"
                  << "    - calls:\n"
                  << "        'vasp.relax $ROOT/training_data/$SCELNAME/$CONFIGID'\n"
                  << "      for each config selected in config_list\n"
                  << "    - The '--write-pos' option makes casm write the POS file  \n"
                  << "      before executing the given command.                     \n\n";




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
    catch(std::exception &e) {
      std::cerr << "ERROR: " << e.what() << ".\n       Exiting..." << std::endl;
      return 1;

    }

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


    try {
      if(!vm.count("config") || (selection == "MASTER")) {
        for(auto it = primclex.selected_config_begin(); it != primclex.selected_config_end(); ++it) {
          if(vm.count("write-pos")) {
            it->write_pos();
          }

          Popen process;

          process.popen(exec + " " + it->get_path().string());

          process.print(std::cout);
        }
      }
      else if(vm.count("config") && fs::exists(fs::path(selection))) {
        ConfigSelection<true> config_select(primclex, selection);
        for(auto it = config_select.selected_config_begin(); it != config_select.selected_config_end(); ++it) {
          if(vm.count("write-pos")) {
            it->write_pos();
          }

          Popen process;

          process.popen(exec + " " + it->get_path().string());

          process.print(std::cout);
        }

      }
      else {
        std::cerr << "ERROR: Invalid input. Option '--config' accepts one argument (either 'MASTER' or a path to a valid configuration selection file)." << std::endl
                  << "       Exiting...\n";
        return 1;
      }
    }
    catch(std::exception &e) {
      std::cerr << "ERROR: Invalid input. Option '--config' accepts one argument (either 'MASTER' or a path to a valid configuration selection file)." << std::endl
                << "       Exiting...\n";
      return 1;
    }




    std::cout << "\n***************************\n" << std::endl;

    std::cout << std::endl;

    return 0;
  };

}


