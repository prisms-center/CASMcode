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

  void print_monte_help(std::ostream &sout, const po::options_description &desc) {
    sout << "\n";
    sout << desc << std::endl;
  }

  void print_monte_desc(std::ostream &sout, const po::options_description &desc) {
    sout << "\n";
    sout << desc << std::endl;

    sout << "DESCRIPTION\n" <<
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

  template<typename MCType>
  int _initial_POSCAR(PrimClex &primclex, const CommandArgs &args, const Completer::MonteOption &monte_opt);

  template<typename MCType>
  int _final_POSCAR(PrimClex &primclex, const CommandArgs &args, const Completer::MonteOption &monte_opt);

  template<typename MCType>
  int _traj_POSCAR(PrimClex &primclex, const CommandArgs &args, const Completer::MonteOption &monte_opt);

  template<typename MCType>
  int _driver(PrimClex &primclex, const CommandArgs &args, const Completer::MonteOption &monte_opt);

  int _run_GrandCanonical(
    PrimClex &primclex,
    const MonteSettings &monte_settings,
    const CommandArgs &args,
    const Completer::MonteOption &monte_opt);

  int _run_Canonical(
    PrimClex &primclex,
    const MonteSettings &monte_settings,
    const CommandArgs &args,
    const Completer::MonteOption &monte_opt);


  int monte_command(const CommandArgs &args) {

    fs::path settings_path;
    std::string verbosity_str;
    po::variables_map vm;
    Index condition_index;

    // Set command line options using boost program_options
    Completer::MonteOption monte_opt;

    try {
      po::store(po::parse_command_line(args.argc, args.argv, monte_opt.desc()), vm); // can throw

      /** --help option
      */
      if(vm.count("help")) {
        print_monte_help(args.log, monte_opt.desc());
        return 0;
      }

      if(vm.count("desc")) {
        print_monte_desc(args.log, monte_opt.desc());
        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems

      settings_path = monte_opt.settings_path();
      verbosity_str = monte_opt.verbosity_str();
      condition_index = monte_opt.condition_index();

      if(vm.count("verbosity")) {
        auto res = Log::verbosity_level(verbosity_str);
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
      args.err_log << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log << monte_opt.desc() << std::endl;
      return 1;
    }
    catch(std::exception &e) {
      args.err_log << "Unhandled Exception reached the top of main: "
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
    Log &err_log = args.err_log;


    const DirectoryStructure &dir = primclex.dir();
    ProjectSettings &set = primclex.settings();

    //Get path to settings json file
    settings_path = fs::absolute(settings_path);

    //args.log << "Example settings so far..." << std::endl;
    //jsonParser example_settings = Monte::example_testing_json_settings(primclex);
    //std::ofstream outsettings("monte_settings.json");
    //example_settings.print(outsettings);

    MonteSettings monte_settings;

    try {
      log.read("Monte Carlo settings");
      log << "from: " << settings_path << "\n";
      monte_settings = MonteSettings(settings_path);
      log << "ensemble: " << monte_settings.ensemble() << "\n";
      log << "method: " << monte_settings.method() << "\n";
      if(args.log.verbosity() == 100) {
        monte_settings.set_debug(true);
      }
      if(monte_settings.debug()) {
        log << "debug: " << monte_settings.debug() << "\n";
      }
      log << std::endl;

    }
    catch(std::exception &e) {
      args.err_log << "ERROR reading Monte Carlo settings.\n\n";
      args.err_log << e.what() << std::endl;
      return 1;
    }

    if(monte_settings.ensemble() == Monte::ENSEMBLE::GrandCanonical) {
      return _run_GrandCanonical(primclex, monte_settings, args, monte_opt);
    }
    else if(monte_settings.ensemble() == Monte::ENSEMBLE::Canonical) {
      return _run_Canonical(primclex, monte_settings, args, monte_opt);
    }

    return ERR_INVALID_ARG;
  }


  template<typename MCType>
  int _initial_POSCAR(PrimClex &primclex, const CommandArgs &args, const Completer::MonteOption &monte_opt) {
    try {
      MCType::SettingsType mc_settings(primclex, monte_opt.settings_path());
      const MCType mc(primclex, mc_settings, args.log);

      args.log.write("Initial POSCAR");
      write_POSCAR_initial(mc, monte_opt.condition_index(), args.log);
      args.log << std::endl;
      return 0;
    }
    catch(std::exception &e) {
      args.err_log << "ERROR printing " << to_string(MCType::ensemble) <<
                   " Monte Carlo initial snapshot for condition: " << monte_opt.condition_index() << "\n\n";
      args.err_log << e.what() << std::endl;
      return 1;
    }
  }

  template<typename MCType>
  int _final_POSCAR(PrimClex &primclex, const CommandArgs &args, const Completer::MonteOption &monte_opt) {
    try {
      MCType::SettingsType mc_settings(primclex, monte_opt.settings_path());
      const MCType mc(primclex, mc_settings, args.log);

      args.log.write("Final POSCAR");
      write_POSCAR_final(mc, monte_opt.condition_index(), args.log);
      args.log << std::endl;
      return 0;
    }
    catch(std::exception &e) {
      args.err_log << "ERROR printing " << to_string(MCType::ensemble) <<
                   " Monte Carlo initial snapshot for condition: " << monte_opt.condition_index() << "\n\n";
      args.err_log << e.what() << std::endl;
      return 1;
    }
  }

  template<typename MCType>
  int _traj_POSCAR(PrimClex &primclex, const CommandArgs &args, const Completer::MonteOption &monte_opt) {
    try {
      MCType::SettingsType mc_settings(primclex, monte_opt.settings_path());
      const MCType mc(primclex, mc_settings, args.log);

      args.log.write("Trajectory POSCAR");
      write_POSCAR_trajectory(mc, monte_opt.condition_index(), args.log);
      args.log << std::endl;
      return 0;
    }
    catch(std::exception &e) {
      args.err_log << "ERROR printing " << to_string(MCType::ensemble) <<
                   " Monte Carlo path snapshots for condition: " << monte_opt.condition_index() << "\n\n";
      args.err_log << e.what() << std::endl;
      return 1;
    }
  }

  template<typename MCType>
  int _driver(PrimClex &primclex, const CommandArgs &args, const Completer::MonteOption &monte_opt) {
    try {
      MCType::SettingsType mc_settings(primclex, monte_opt.settings_path());
      MonteDriver<MCType> driver(primclex, mc_settings, args.log, args.err_log);
      driver.run();
      return 0;
    }
    catch(std::exception &e) {
      args.err_log << "ERROR running " << to_string(MCType::ensemble) << " Monte Carlo.\n\n";
      args.err_log << e.what() << std::endl;
      return 1;
    }
  }

  int _run_GrandCanonical(
    PrimClex &primclex,
    const MonteSettings &monte_settings,
    const CommandArgs &args,
    const Completer::MonteOption &monte_opt) {

    typedef GrandCanonical MCType;

    if(vm.count("initial-POSCAR")) {
      return _initial_POSCAR<MCType>(primclex, args, monte_opt);
    }
    else if(vm.count("final-POSCAR")) {
      return _final_POSCAR<MCType>(primclex, args, monte_opt);
    }
    else if(vm.count("traj-POSCAR")) {
      return _traj_POSCAR<MCType>(primclex, args, monte_opt);
    }
    else if(monte_settings.method() == Monte::METHOD::LTE1) {

      try {

        GrandCanonicalSettings gc_settings(primclex, settings_path);

        if(gc_settings.dependent_runs()) {
          throw std::invalid_argument("ERROR in LTE1 calculation: dependents_runs must be false");
        }

        bool ok = false;
        if(gc_settings.is_motif_configname() &&
           (gc_settings.motif_configname() == "auto" ||
            gc_settings.motif_configname() == "restricted_auto")) {
          ok = true;
        }

        if(!ok) {
          throw std::invalid_argument("ERROR in LTE1 calculation: must use one of\n"
                                      "  \"driver\"/\"motif\"/\"configname\": \"auto\"\n"
                                      "  \"driver\"/\"motif\"/\"configname\": \"restricted_auto\"");
        }

        GrandCanonicalDirectoryStructure dir(gc_settings.output_directory());
        if(gc_settings.write_csv()) {
          if(fs::exists(dir.results_csv())) {
            args.err_log << "Existing file at: " << dir.results_csv() << std::endl;
            args.err_log << "  Exiting..." << std::endl;
            return ERR_EXISTING_FILE;
          }
        }
        if(gc_settings.write_json()) {
          if(fs::exists(dir.results_json())) {
            args.err_log << "Existing file at: " << dir.results_json() << std::endl;
            args.err_log << "  Exiting..." << std::endl;
            return ERR_EXISTING_FILE;
          }
        }

        GrandCanonical gc(primclex, gc_settings, log);

        // config, param_potential, T,
        log.custom("LTE Calculation");
        log << "Phi_LTE(1) = potential_energy_gs - kT*ln(Z'(1))/N" << std::endl;
        log << "Z'(1) = sum_i(exp(-dPE_i/kT), summing over ground state and single spin flips" << std::endl;
        log << "dPE_i: (potential_energy_i - potential_energy_gs)*N" << "\n\n" << std::endl;

        auto init = gc_settings.initial_conditions();
        auto incr = init;
        int num_conditions = 1;

        if(monte_settings.drive_mode() == Monte::DRIVE_MODE::INCREMENTAL) {

          incr = gc_settings.incremental_conditions();
          auto final = gc_settings.final_conditions();
          num_conditions = (final - init) / incr + 1;
        }

        std::string configname;

        auto cond = init;
        for(int index = 0; index < num_conditions; ++index) {

          configname = gc.set_state(cond, gc_settings).second;

          if(gc.debug()) {
            const auto &comp_converter = gc.primclex().composition_axes();
            args.log << "formation_energy: " << std::setprecision(12) << gc.formation_energy() << std::endl;
            args.log << "  components: " << jsonParser(gc.primclex().composition_axes().components()) << std::endl;
            args.log << "  comp_n: " << gc.comp_n().transpose() << std::endl;
            args.log << "  param_chem_pot: " << gc.conditions().param_chem_pot().transpose() << std::endl;
            args.log << "  comp_x: " << comp_converter.param_composition(gc.comp_n()).transpose() << std::endl;
            args.log << "potential energy: " << std::setprecision(12) << gc.potential_energy() << std::endl << std::endl;
          }

          double phi_LTE1 = gc.lte_grand_canonical_free_energy();

          log.write("Output files");
          write_lte_results(gc_settings, gc, phi_LTE1, configname, log);
          log << std::endl;
          cond += incr;

          log << std::endl;
        }

        return 0;

      }
      catch(std::exception &e) {
        args.err_log << "ERROR calculating single spin flip LTE grand canonical potential.\n\n";
        args.err_log << e.what() << std::endl;
        return 1;
      }
    }
    else if(monte_settings.method() == Monte::METHOD::Metropolis) {
      return _driver<GrandCanonical>(primclex, args, monte_opt);
    }
    else {
      args.err_log << "ERROR running " << to_string(GrandCanonical::ensemble) << " Monte Carlo. No valid option given.\n\n";
      return ERR_INVALID_INPUT_FILE;
    }
  }

  int _run_Canonical(
    PrimClex &primclex,
    const MonteSettings &monte_settings,
    const CommandArgs &args,
    const Completer::MonteOption &monte_opt) {

    typedef Canonical MCType;

    if(vm.count("initial-POSCAR")) {
      return _initial_POSCAR<MCType>(primclex, args, monte_opt);
    }
    else if(vm.count("final-POSCAR")) {
      return _final_POSCAR<MCType>(primclex, args, monte_opt);
    }
    else if(vm.count("traj-POSCAR")) {
      return _traj_POSCAR<MCType>(primclex, args, monte_opt);
    }
    else if(monte_settings.method() == Monte::METHOD::Metropolis) {
      return _driver<MCType>(primclex, args, monte_opt);
    }
    else {
      args.err_log << "ERROR running " << to_string(GrandCanonical::ensemble) << " Monte Carlo. No valid option given.\n\n";
      return ERR_INVALID_INPUT_FILE;
    }
  }
}

