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
    else {
      primclex.read_global_orbitree(dir.clust(set.bset()));
      primclex.generate_full_nlist();
      primclex.generate_supercell_nlists();
    }

    Clexulator clexulator(set.global_clexulator(),
                          dir.clexulator_dir(set.bset()),
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

    // -- count number of selected configurations, and check they all have properties
    bool ok = true;
    int N_values = 0;
    for(auto it = config_select.selected_config_cbegin(); it != config_select.selected_config_cend(); ++it) {
      if(!it->delta_properties().contains("relaxed_energy")) {
        std::cerr << "Configuration: " << it->name() << " does not have a formation energy." << std::endl;
        ok = false;
      }
      N_values++;
    }
    if(!ok) {
      std::cerr << "\nPlease check your configurations and re-try." << std::endl;
      return 1;
    }

    if(N_values == 0) {
      std::cerr << "\nDid not find any selected values. Please update your selection and re-try." << std::endl;
    }

    std::cout << "Calculating convex hull for selected configurations..." << std::endl << std::endl;
    jsonParser hulljson;
    hulljson = update_hull_props(primclex, config_select.selected_config_begin(), config_select.selected_config_end());
    std::cout << "  DONE." << std::endl << std::endl;

    // -- write 'energy' file ----
    {
      DataFormatter<Configuration> formatter(ConfigIO::formation_energy(),
                                             ConfigIO::constant_value("weight", 1.0),
                                             ConfigIO::param_composition(),
                                             ConfigIO::dist_from_hull(),
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
      formatter.push_back(ConfigIO::corr(clexulator));

      fs::ofstream sout;
      sout.open(corr_in_file);
      sout << N_corr << " # basis functions\n";
      sout << N_values << " # training values\n";
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
