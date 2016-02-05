#include "fit.hh"

#include <string>

#include "casm_functions.hh"
#include "casm/CASM_classes.hh"

namespace CASM {

  int fit_command(int argc, char *argv[]) {

    fs::path selection;
    po::variables_map vm;
    bool force;


    try {

      // Set command line options using boost program_options
      po::options_description desc("'casm fit' usage");
      desc.add_options()
      ("help,h", "Print help message")
      ("config,c", po::value<fs::path>(&selection), "Selected configurations are used as training data for ECI fitting. If not specified, or 'MASTER' given, uses master list selection.")
      ("force,f", "Overrwrite output file");

      try {
        po::store(po::parse_command_line(argc, argv, desc), vm); // can throw

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
      std::cout << "Error in 'casm energy': No casm project found." << std::endl;
      return 1;
    }

    std::cout << "\n***************************\n" << std::endl;


    DirectoryStructure dir(root);
    ProjectSettings set(root);

    fs::path energy_file = dir.energy(set.clex(), set.calctype(), set.ref(), set.bset(), set.eci());
    fs::path eci_in_file = dir.eci_in(set.clex(), set.calctype(), set.ref(), set.bset(), set.eci());
    fs::path corr_in_file = dir.corr_in(set.clex(), set.calctype(), set.ref(), set.bset(), set.eci());

    if(!vm.count("force")) {
      if(fs::exists(energy_file)) {
        std::cerr << "File " << energy_file << " already exists. Use --force to force overwrite." << std::endl;
        return 1;
      }
      if(fs::exists(eci_in_file)) {
        std::cerr << "File " << eci_in_file << " already exists. Use --force to force overwrite." << std::endl;
        return 1;
      }
      if(fs::exists(corr_in_file)) {
        std::cerr << "File " << corr_in_file << " already exists. Use --force to force overwrite." << std::endl;
        return 1;
      }
    }


    // initialize primclex
    std::cout << "Initialize primclex: " << root << std::endl << std::endl;
    PrimClex primclex(root, std::cout);
    std::cout << "  DONE." << std::endl << std::endl;


    if(!fs::exists(dir.clexulator_src(set.name(), set.bset()))) {
      std::cerr << "No basis set found. Please use 'casm bset' first." << std::endl;
      return 1;
    }
    
    Clexulator clexulator(set.global_clexulator(),
                          dir.clexulator_dir(set.bset()),
                          primclex.nlist(),
                          set.compile_options(),
                          set.so_options());

    int N_corr = clexulator.corr_size();

    ConfigSelection<false> config_select;
    if(!vm.count("config") || selection == "MASTER") {
      config_select = ConfigSelection<false>(primclex);
    }
    else {
      config_select = ConfigSelection<false>(primclex, selection);
    }

    // -- write 'energy' file ----
    {
      DataFormatter<Configuration> formatter(ConfigIO::formation_energy(),
                                             ConfigIO::ConstantValue<double>("weight", 1.0),
                                             ConfigIO::configname());

      fs::ofstream sout;
      sout.open(energy_file);
      sout << formatter(config_select.selected_config_cbegin(), config_select.selected_config_cend());
      sout.close();

      std::cout << "Wrote: " << energy_file << "\n\n";
    }

    // -- write 'corr.in' file ----
    {
      DataFormatter<Configuration> formatter;
      formatter.push_back(ConfigIO::Corr(clexulator));

      fs::ofstream sout;
      sout.open(corr_in_file);
      sout << N_corr << " # basis functions\n";
      sout << std::distance(config_select.selected_config_begin(), config_select.selected_config_end()) << " # training values\n";
      sout << "Correlation matrix:\n";
      sout << FormatFlag(sout).print_header(false);
      sout << formatter(config_select.selected_config_cbegin(), config_select.selected_config_cend());
      sout.close();

      std::cout << "Wrote: " << corr_in_file << "\n\n";
    }

    // -- write eci.in ----------------
    {
      primclex.get_global_orbitree().write_eci_in(eci_in_file.string());

      std::cout << "Wrote: " << eci_in_file << "\n\n";
    }

    return 0;
  };

}
