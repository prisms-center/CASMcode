#include <string>
#include <iostream>

#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"
#include "casm/monte_carlo/MonteIO.hh"
#include "casm/monte_carlo/MonteDriver.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  void print_monte_help(const po::options_description &desc) {
    std::cout << "\n";
    std::cout << desc << std::endl;

    std::cout << "DESCRIPTION\n" <<
              "  Perform Monte Carlo calculations.                          \n\n" <<

              "  casm monte --settings input_file.json                      \n" <<
              "    - Run Monte Carlo calculations given the input file      \n" <<
              "      settings.                                              \n" <<
              "    - See 'casm format --monte' for a description of the     \n" <<
              "      Monte Carlo input file.                                \n\n" <<

              "  casm monte --settings input_file.json --initial-POSCAR 3     \n" <<
              "    - Write a POSCAR.initial file containing the initial state of\n" <<
              "      the Monte Carlo calculation. The argument is a condition\n" <<
              "      index specifying which run is being requested.\n" <<
              "    - Written at: output_directory/conditions.3/trajectory/POSCAR.initial\n\n" <<

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
              "      where i is the sample index.                          \n\n";

  }

  namespace Completer {

    MonteOption::MonteOption(): OptionHandlerBase("monte") {};

    void MonteOption::initialize() {
      add_help_suboption();
      add_verbosity_suboption();
      add_settings_suboption();

      m_desc.add_options()
      ("initial-POSCAR", po::value<Index>(&m_condition_index), "Given the condition index, print a POSCAR for the initial state of a monte carlo run.")
      ("final-POSCAR", po::value<Index>(&m_condition_index), "Given the condition index, print a POSCAR for the final state of a monte carlo run.")
      ("traj-POSCAR", po::value<Index>(&m_condition_index), "Given the condition index, print POSCARs for the state at every sample of monte carlo run. Requires an existing trajectory file.");
      return;
    }

    Index MonteOption::condition_index() const {
      return m_condition_index;
    };

  }


  int monte_command(const CommandArgs &args) {

    //fs::path settings_path;
    //std::string verbosity_str;
    po::variables_map vm;
    //Index condition_index;

    Completer::MonteOption monte_opt;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, monte_opt.desc()), vm); // can throw

      /** --help option
      */
      if(vm.count("help")) {
        print_monte_help(monte_opt.desc());
        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      if(vm.count("verbosity")) {
        auto res = Log::verbosity_level(monte_opt.verbosity_str());
        if(!res.first) {
          args.err_log.error("--verbosity");
          args.err_log << "Expected: 'none', 'quiet', 'standard', 'verbose', "
                       "'debug', or an integer 0-100 (0: none, 100: all)" << "\n" << std::endl;
          return ERR_INVALID_ARG;
        }
        args.log.set_verbosity(res.second);
      }

    }
    catch(po::error &e) {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      std::cerr << monte_opt.desc() << std::endl;
      return 1;
    }
    catch(std::exception &e) {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return 1;

    }

    const fs::path &root = args.root;
    if(root.empty()) {
      args.err_log.error("No casm project found");
      args.err_log << std::endl;
      return ERR_NO_PROJ;
    }

    // If 'args.primclex', use that, else construct PrimClex in 'uniq_primclex'
    // Then whichever exists, store reference in 'primclex'
    std::unique_ptr<PrimClex> uniq_primclex;
    PrimClex &primclex = make_primclex_if_not(args, uniq_primclex);
    Log &log = args.log;


    const DirectoryStructure &dir = primclex.dir();
    ProjectSettings &set = primclex.settings();

    //Get path to settings json file
    //monte.settings_path() = fs::absolute(monte.settings_path());
    fs::path abs_settings_path = fs::absolute(monte_opt.settings_path());

    //std::cout << "Example settings so far..." << std::endl;
    //jsonParser example_settings = Monte::example_testing_json_settings(primclex);
    //std::ofstream outsettings("monte_settings.json");
    //example_settings.print(outsettings);

    MonteSettings monte_settings;

    try {
      log.read("Monte Carlo settings");
      log << "from: " << abs_settings_path << "\n";
      monte_settings = MonteSettings(abs_settings_path);
    }
    catch(std::exception &e) {
      std::cerr << "ERROR reading Monte Carlo settings.\n\n";
      std::cerr << e.what() << std::endl;
      return 1;
    }
    log << "ensemble: " << monte_settings.ensemble() << "\n";
    log << "method: " << monte_settings.method() << "\n";

    if(args.log.verbosity() == 100) {
      monte_settings.set_debug(true);
    }
    if(monte_settings.debug()) {
      log << "debug: " << monte_settings.debug() << "\n";
    }
    log << std::endl;

    if(monte_settings.ensemble() == Monte::ENSEMBLE::GrandCanonical) {

      if(vm.count("initial-POSCAR")) {
        try {
          GrandCanonicalSettings gc_settings(abs_settings_path);
          const GrandCanonical gc(primclex, gc_settings, log);

          log.write("Initial POSCAR");
          write_POSCAR_initial(gc, monte_opt.condition_index(), log);
          log << std::endl;
        }
        catch(std::exception &e) {
          std::cerr << "ERROR printing Grand Canonical Monte Carlo initial snapshot for condition: " << monte_opt.condition_index() << "\n\n";
          std::cerr << e.what() << std::endl;
          return 1;
        }
      }
      else if(vm.count("final-POSCAR")) {
        try {
          GrandCanonicalSettings gc_settings(abs_settings_path);
          const GrandCanonical gc(primclex, gc_settings, log);

          log.write("Final POSCAR");
          write_POSCAR_final(gc, monte_opt.condition_index(), log);
          log << std::endl;
        }
        catch(std::exception &e) {
          std::cerr << "ERROR printing Grand Canonical Monte Carlo final snapshot for condition: " << monte_opt.condition_index() << "\n\n";
          std::cerr << e.what() << std::endl;
          return 1;
        }
      }
      else if(vm.count("traj-POSCAR")) {
        try {
          GrandCanonicalSettings gc_settings(abs_settings_path);
          const GrandCanonical gc(primclex, gc_settings, log);

          log.write("Trajectory POSCARs");
          write_POSCAR_trajectory(gc, monte_opt.condition_index(), log);
          log << std::endl;
        }
        catch(std::exception &e) {
          std::cerr << "ERROR printing Grand Canonical Monte Carlo path snapshots for condition: " << monte_opt.condition_index() << "\n\n";
          std::cerr << e.what() << std::endl;
          return 1;
        }
      }
      else if(monte_settings.method() == Monte::METHOD::LTE1) {

        try {

          GrandCanonicalSettings gc_settings(abs_settings_path);
          GrandCanonicalDirectoryStructure dir(gc_settings.output_directory());
          if(gc_settings.write_csv()) {
            if(fs::exists(dir.results_csv())) {
              std::cerr << "Existing file at: " << dir.results_csv() << std::endl;
              std::cerr << "  Exiting..." << std::endl;
              return ERR_EXISTING_FILE;
            }
          }
          if(gc_settings.write_json()) {
            if(fs::exists(dir.results_json())) {
              std::cerr << "Existing file at: " << dir.results_json() << std::endl;
              std::cerr << "  Exiting..." << std::endl;
              return ERR_EXISTING_FILE;
            }
          }

          GrandCanonical gc(primclex, gc_settings, log);

          // config, param_potential, T,
          log.custom("LTE Calculation");
          log << "Phi_LTE(1) = potential_energy_gs - kT*ln(Z'(1))/N" << std::endl;
          log << "Z'(1) = sum_i(exp(-dPE_i/kT), summing over ground state and single spin flips" << std::endl;
          log << "dPE_i: (potential_energy_i - potential_energy_gs)*N" << std::endl << std::endl;

          auto init = gc_settings.initial_conditions();
          auto incr = init;
          int num_conditions = 1;

          if(monte_settings.drive_mode() == Monte::DRIVE_MODE::INCREMENTAL) {

            incr = gc_settings.incremental_conditions();
            auto final = gc_settings.final_conditions();
            num_conditions = (final - init) / incr + 1;
          }

          auto cond = init;
          for(int index = 0; index < num_conditions; ++index) {

            if(index != 0) {
              gc.set_conditions(cond);
            }

            if(gc.debug()) {
              const auto &comp_converter = gc.primclex().composition_axes();
              std::cout << "formation_energy: " << std::setprecision(12) << gc.formation_energy() << std::endl;
              std::cout << "  components: " << jsonParser(gc.primclex().composition_axes().components()) << std::endl;
              std::cout << "  comp_n: " << gc.comp_n().transpose() << std::endl;
              std::cout << "  param_chem_pot: " << gc.conditions().param_chem_pot().transpose() << std::endl;
              std::cout << "  comp_x: " << comp_converter.param_composition(gc.comp_n()).transpose() << std::endl;
              std::cout << "potential energy: " << std::setprecision(12) << gc.potential_energy() << std::endl << std::endl;
            }

            double phi_LTE1 = gc.lte_grand_canonical_free_energy();

            write_lte_results(gc_settings, gc, phi_LTE1, log);
            cond += incr;

          }

        }
        catch(std::exception &e) {
          std::cerr << "ERROR calculating single spin flip LTE grand canonical potential.\n\n";
          std::cerr << e.what() << std::endl;
          return 1;
        }
      }
      else if(monte_settings.method() == Monte::METHOD::Metropolis) {

        try {

          //std::cout << "\n-------------------------------\n";
          //monte_settings.print(std::cout);
          //std::cout << "\n-------------------------------\n\n";

          MonteDriver<GrandCanonical> driver(primclex, GrandCanonicalSettings(abs_settings_path), log);
          driver.run();
        }
        catch(std::exception &e) {
          std::cerr << "ERROR running Grand Canonical Monte Carlo.\n\n";
          std::cerr << e.what() << std::endl;
          return 1;
        }
      }
      else {
        std::cerr << "ERROR running Grand Canonical Monte Carlo. No valid option given.\n\n";
        return ERR_INVALID_INPUT_FILE;
      }
    }


    return 0;
  }
}

