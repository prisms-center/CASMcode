#include "monte.hh"

#include <string>
#include <iostream>

#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/monte_carlo/MonteIO.hh"
#include "casm/monte_carlo/MonteDriver.hh"
#include "casm_functions.hh"

namespace CASM {

  int monte_command(int argc, char *argv[]) {

    std::string settingsfile;
    fs::path settings_path;
    po::variables_map vm;

    try {

      // Set command line options using boost program_options
      po::options_description desc("'casm monte' usage");
      desc.add_options()
      ("help,h", "Print help message")
      ("settings,s", po::value<std::string>(&settingsfile)->required(), "json file with all the settings you need for your simulation");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;

          std::cout << "Monte is still under construction. Come back later." << std::endl;

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


    fs::path rootpath = find_casmroot(fs::current_path());
    if(rootpath.empty()) {
      std::cout << "Error in 'casm monte': No casm project found." << std::endl;
      return 1;
    }

    std::cout << "\n***************************\n" << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << rootpath << std::endl << std::endl;
    PrimClex primclex(rootpath, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    const DirectoryStructure &dir = primclex.dir();
    ProjectSettings &set = primclex.settings();

    //Get path to settings json file
    if(vm.count("settings")) {
      settings_path = settingsfile;
    }

    //std::cout << "Example settings so far..." << std::endl;
    //jsonParser example_settings = Monte::example_testing_json_settings(primclex);
    //std::ofstream outsettings("monte_settings.json");
    //example_settings.print(outsettings);
    
    MonteSettings monte_settings;
    
    try {
      std::cout << "Reading Monte Carlo settings: " << settings_path << std::endl;
      monte_settings = MonteSettings(settings_path);
      std::cout << "\n-------------------------------\n";
      monte_settings.print(std::cout);
      std::cout << "\n-------------------------------\n\n";
      std::cout << "  DONE." << std::endl << std::endl;
    
    }
    catch(std::exception& e) {
      std::cerr << "ERROR reading Monte Carlo settings.\n\n";
      std::cerr << e.what() << std::endl;
      return 1;
    }
    
    
    if(monte_settings.type() == Monte::TYPE::GrandCanonical) {
    
      try {
        
        std::cout << "Constructing Grand Canonical Monte Carlo driver" << std::endl;
        MonteDriver<GrandCanonical> driver(primclex, GrandCanonicalSettings(settings_path));
        std::cout << "  DONE." << std::endl << std::endl;

        std::cout << "Begin Grand Canonical Monte Carlo runs" << std::endl;
        driver.run(std::cout);
        std::cout << "  DONE." << std::endl << std::endl;
        
      }
      catch(std::exception& e) {
        std::cerr << "ERROR running Grand Canonical Monte Carlo.\n\n";
        std::cerr << e.what() << std::endl;
        return 1;
      }
      
    }
    
    return 0;
  }
}

