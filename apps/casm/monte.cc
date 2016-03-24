#include "monte.hh"

#include <string>
#include <iostream>

#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"
#include "casm/monte_carlo/MonteIO.hh"
#include "casm/monte_carlo/MonteDriver.hh"
#include "casm_functions.hh"

namespace CASM {

  int monte_command(int argc, char *argv[]) {

    fs::path settings_path;
    po::variables_map vm;
    Index condition_index;

    try {

      // Set command line options using boost program_options
      po::options_description desc("'casm monte' usage");
      desc.add_options()
      ("help,h", "Print help message")
      ("settings,s", po::value<fs::path>(&settings_path)->required(), "The Monte Carlo input file. See 'casm format --monte'.")
      ("final-POSCAR", po::value<Index>(&condition_index), "Given the condition index, print a POSCAR for the final state of a monte carlo run.")
      ("traj-POSCAR", po::value<Index>(&condition_index), "Given the condition index, print POSCARs for the state at every sample of monte carlo run. Requires an existing trajectory file.");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

        /** --help option
        */
        if(vm.count("help")) {
          std::cout << "\n";
          std::cout << desc << std::endl;
          
          std::cout << "DESCRIPTION\n" <<
                       "  Perform Monte Carlo calculations.                          \n\n" <<
          
                       "  casm monte --settings input_file.json                      \n" <<
                       "    - Run Monte Carlo calculations given the input file      \n" <<
                       "      settings.                                              \n" <<
                       "    - See 'casm format --monte' for a description of the     \n" <<
                       "      Monte Carlo input file.                                \n\n" <<
                       
                       "  casm monte --settings input_file.json --final-POSCAR 3     \n" <<
                       "    - Write a POSCAR.final file containing the final state of\n" <<
                       "      the Monte Carlo calculation. The argument is a condition\n" <<
                       "      index specifying which run is being requested.\n" <<
                       "    - Written at: output_directory/conditions.3/trajectory/POSCAR.final\n\n" <<
                       
                       "  casm monte --settings input_file.json --traj-POSCAR 5     \n" <<
                       "    - Write the Monte Carlo calculation trajectory as a     \n" <<
                       "      series of POSCAR files containing the state of the    \n" <<
                       "      Monte Carlo calculation every time a sample was taken.\n" <<
                       "      The argument is a condition index specifying which run\n" <<
                       "      is being requested.                                   \n" <<
                       "    - The trajectory file must exist. This is generated when\n" <<
                       "      using input option \"data\"/\"storage\"/\"write_trajectory\" = true  \n" <<
                       "    - Written at: output_directory/conditions.5/trajectory/POSCAR.i,\n" <<
                       "      where i is the sample index.\n\n";
                       
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


    fs::path root = find_casmroot(fs::current_path());
    if(root.empty()) {
      std::cout << "Error in 'casm monte': No casm project found." << std::endl;
      return 1;
    }

    std::cout << "\n***************************\n" << std::endl;

    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;

    const DirectoryStructure &dir = primclex.dir();
    ProjectSettings &set = primclex.settings();
    
    //Get path to settings json file
    settings_path = fs::absolute(settings_path);
  
    //std::cout << "Example settings so far..." << std::endl;
    //jsonParser example_settings = Monte::example_testing_json_settings(primclex);
    //std::ofstream outsettings("monte_settings.json");
    //example_settings.print(outsettings);
    
    MonteSettings monte_settings;
    
    try {
      std::cout << "Reading Monte Carlo settings: " << settings_path << std::endl;
      monte_settings = MonteSettings(settings_path);
      std::cout << "  DONE." << std::endl << std::endl;
    
    }
    catch(std::exception& e) {
      std::cerr << "ERROR reading Monte Carlo settings.\n\n";
      std::cerr << e.what() << std::endl;
      return 1;
    }
    
    if(monte_settings.type() == Monte::TYPE::GrandCanonical) {
      
      if(vm.count("final-POSCAR")) {
        try {
          GrandCanonicalSettings gc_settings(settings_path);
          const GrandCanonical gc(primclex, gc_settings);
          write_POSCAR_final(gc, condition_index);
        }
        catch(std::exception& e) {
          std::cerr << "ERROR printing Grand Canonical Monte Carlo final snapshot for condition: " << condition_index << "\n\n";
          std::cerr << e.what() << std::endl;
          return 1;
        }
      }
      else if(vm.count("traj-POSCAR")) {
        try {
          GrandCanonicalSettings gc_settings(settings_path);
          const GrandCanonical gc(primclex, gc_settings);
          write_POSCAR_trajectory(gc, condition_index);
        }
        catch(std::exception& e) {
          std::cerr << "ERROR printing Grand Canonical Monte Carlo path snapshots for condition: " << condition_index << "\n\n";
          std::cerr << e.what() << std::endl;
          return 1;
        }
      }
      else {
        try {
          
          //std::cout << "\n-------------------------------\n";
          //monte_settings.print(std::cout);
          //std::cout << "\n-------------------------------\n\n";
          
          std::cout << "Constructing Grand Canonical Monte Carlo driver" << std::endl;
          MonteDriver<GrandCanonical> driver(primclex, GrandCanonicalSettings(settings_path));
          std::cout << "  DONE." << std::endl << std::endl;

          std::cout << "Begin Grand Canonical Monte Carlo runs" << std::endl;
          driver.run();
          std::cout << "  DONE." << std::endl << std::endl;
          
        }
        catch(std::exception& e) {
          std::cerr << "ERROR running Grand Canonical Monte Carlo.\n\n";
          std::cerr << e.what() << std::endl;
          return 1;
        }
      }
    }

    
    return 0;
  }
}

