#ifndef CASM_MonteDriver_HH
#define CASM_MonteDriver_HH

#include <string>
#include "casm/external/boost.hh"

#include "casm/monte_carlo/MonteIO.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/MonteCarloEnum_impl.hh"
#include "casm/monte_carlo/MonteSettings.hh"

namespace CASM {

  /**
   * MonteDriver consists of a specialized MonteCarlo object and a list of conditions
   * to equilibrate at. The first condition in the list requires a starting configuration
   * (read from the setting), while subsequent conditions are calculated using the
   * final state of the previous condition.
   *
   * The different kinds of drive modes the user can specify are:
   * INCREMENTAL:   Given a delta in condition values, increment the conditions by the delta after each point
   * CUSTOM:        Calculate for a list of condition values
   */

  template<typename RunType>
  class MonteDriver {

  public:
    typedef typename RunType::CondType CondType;
    typedef typename RunType::SettingsType SettingsType;

    /// \brief Constructor via MonteSettings
    MonteDriver(PrimClex &primclex, const SettingsType &settings, Log &_log, Log &_err_log);

    /// \brief Run everything requested by the MonteSettings
    void run();

  private:

    /// run in debug mode?
    bool debug() const {
      return m_debug;
    }

    ///Return the appropriate std::vector of conditions to visit based from settings. Use for construction.
    std::vector<CondType> make_conditions_list(const PrimClex &primclex, const SettingsType &settings);

    ///Converge the MonteCarlo for conditions 'cond_index'
    void single_run(Index cond_index);

    ///Check for existing calculations to find starting conditions
    Index _find_starting_conditions() const;


    /// target for log messages
    Log &m_log;

    /// target for error messages
    Log &m_err_log;

    ///Copy of initial settings given at construction. Will expand to have MonteCarlo states dumped into it.
    SettingsType m_settings;

    /// describes where to write output
    MonteCarloDirectoryStructure m_dir;

    ///Specifies how to build the conditions list from the settings
    const Monte::DRIVE_MODE m_drive_mode;

    ///Specialized Monte Carlo object to use throughout
    RunType m_mc;

    ///List of specialized conditions to visit in sequential order. Does not include initial conditions.
    const std::vector<CondType> m_conditions_list;

    /// run in debug mode?
    bool m_debug;

    /// Enumerated configurations encountered during Monte Carlo calculations
    notstd::cloneable_ptr<MonteCarloEnum> m_enum;
  };


  /// Perform a single monte carlo step, return true if accepted
  template<typename RunType>
  bool monte_carlo_step(RunType &monte_run);


  template<typename RunType>
  MonteDriver<RunType>::MonteDriver(PrimClex &primclex, const SettingsType &settings, Log &_log, Log &_err_log):
    m_log(_log),
    m_err_log(_err_log),
    m_settings(settings),
    m_dir(m_settings.output_directory()),
    m_drive_mode(m_settings.drive_mode()),
    m_mc(primclex, m_settings, _log),
    m_conditions_list(make_conditions_list(primclex, m_settings)),
    m_debug(m_settings.debug()),
    m_enum(m_settings.is_enumeration() ? new MonteCarloEnum(primclex, settings, _log, m_mc) : nullptr) {
  }

  /// \brief Run calculations for all conditions, outputting data as you finish each one
  ///
  /// - Assumes existing "output_dir/conditions.i/final_state.json" files indicates finished
  ///   calculations that are already included in the results summary "output_dir/results.X",
  ///   and that no other results are written to the results summary.
  /// - If there are existing results, uses "output_dir/conditions.i/final_state.json" as
  ///   the initial state for the next run
  template<typename RunType>
  void MonteDriver<RunType>::run() {

    m_log.check("For existing calculations");

    if(!m_settings.write_json() && !m_settings.write_csv()) {
      throw std::runtime_error(
        std::string("No valid monte carlo output format.\n") +
        "  Expected [\"data\"][\"storage\"][\"output_format\"] to contain a string or array of strings.\n" +
        "  Valid options are 'csv' or 'json'.");
    }

    // Skip any conditions that have already been calculated and saved
    Index start_i = _find_starting_conditions();

    // check if we'll be repeating any calculations that already have files written
    std::vector<Index> repeats;
    for(Index i = start_i; i < m_conditions_list.size(); ++i) {
      if(fs::exists(m_dir.conditions_dir(i))) {
        repeats.push_back(i);
      }
    }

    if(start_i == m_conditions_list.size()) {
      m_log << "calculations already complete." << std::endl;
      return;
    }

    // if existing calculations
    if(start_i > 0 || repeats.size() > 0) {

      m_log << "found existing calculations\n";
      m_log << "will begin with condition " << start_i << "\n";

      if(repeats.size()) {
        jsonParser json;
        to_json(repeats, json);
        m_log << "will overwrite existing results for condition(s): " << json << "\n";
      }
    }
    else {
      m_log << "did not find existing calculations\n";
    }
    m_log << std::endl;

    if(m_settings.dependent_runs()) {

      // if starting from initial condition
      if(start_i == 0) {

        // set intial state
        m_mc.set_state(m_conditions_list[0], m_settings);

        // perform any requested explicit equilibration passes
        if(m_settings.dependent_runs() && m_settings.is_equilibration_passes_first_run()) {

          auto equil_passes = m_settings.equilibration_passes_first_run();

          m_log.write("DoF");
          m_log << "write: " << m_dir.initial_state_firstruneq_json(0) << "\n" << std::endl;

          jsonParser json;
          fs::create_directories(m_dir.conditions_dir(0));
          to_json(m_mc.configdof(), json).write(m_dir.initial_state_firstruneq_json(0));

          m_log.begin("Equilibration passes");
          m_log << equil_passes << " equilibration passes\n" << std::endl;

          MonteCounter equil_counter(m_settings, m_mc.steps_per_pass());
          while(equil_counter.pass() != equil_passes) {
            monte_carlo_step(m_mc);
            equil_counter++;
          }
        }
      }
      else {

        // read end state of previous condition
        ConfigDoF configdof = m_mc.configdof();
        from_json(configdof, jsonParser(m_dir.final_state_json(start_i - 1)));

        m_mc.set_state(
          m_conditions_list[start_i],
          configdof,
          std::string("Using: ") + m_dir.final_state_json(start_i - 1).string());
      }
    }

    // Run for all conditions, outputting data as you finish each one
    for(Index i = start_i; i < m_conditions_list.size(); i++) {
      if(!m_settings.dependent_runs()) {
        m_mc.set_state(m_conditions_list[i], m_settings);
      }
      else if(i != start_i) {

        m_mc.set_conditions(m_conditions_list[i]);

        m_log.custom("Continue with existing DoF");
        m_log << std::endl;
      }

      single_run(i);

      m_log << std::endl;
    }

    return;
  }

  /// \brief Checks existing files to determine where to restart a path
  ///
  /// - Will overwrite or cause to overwrite files in cases where the
  ///   final state or results summary do not exist for some conditions
  ///
  template<typename RunType>
  Index MonteDriver<RunType>::_find_starting_conditions() const {

    Index start_i = 0;
    Index start_max = m_conditions_list.size();
    Index start_json = m_settings.write_json() ? 0 : start_max;
    Index start_csv = m_settings.write_csv() ? 0 : start_max;
    jsonParser json_results;
    fs::ifstream csv_results;
    std::string str;
    std::stringstream ss;

    // can start with condition i+1 if results (i) and final_state.json (i) exist

    // check JSON files
    if(m_settings.write_json() && fs::exists(m_dir.results_json())) {

      json_results.read(m_dir.results_json());

      // can start with i+1 if results[i] and final_state.json (i) exist
      while(json_results.begin()->size() > start_json && fs::exists(m_dir.final_state_json(start_i)) && start_i < start_max) {
        ++start_json;
      }

      start_max = start_json;
    }

    // check CSV files
    if(m_settings.write_csv() && fs::exists(m_dir.results_csv())) {

      csv_results.open(m_dir.results_csv());

      // read header
      std::getline(csv_results, str);
      ss << str << "\n";

      // can start with i+1 if results[i] and final_state.json (i) exist
      while(!csv_results.eof() && fs::exists(m_dir.final_state_json(start_csv)) && start_i < start_max) {
        ++start_csv;
        std::getline(csv_results, str);
        ss << str << "\n";
      }

      start_max = start_csv;
    }
    csv_results.close();

    // use minimum of allowed starting conditions, in case a difference is found
    start_i = std::min(start_json, start_csv);


    // update results summary files to remove any conditions that must be re-calculated

    // for JSON
    if(m_settings.write_json() && fs::exists(m_dir.results_json())) {

      // if results size > start_i, must fix results by removing results that will be re-run
      jsonParser finished_results;
      for(auto it = json_results.begin(); it != json_results.end(); ++it) {
        jsonParser &ref = finished_results[it.name()].put_array();
        for(Index i = 0; i < start_i; ++i) {
          ref.push_back((*it)[i]);
        }
      }
      m_log << "update: " << m_dir.results_json() << "\n";
      finished_results.write(m_dir.results_json());
    }

    // for csv
    if(m_settings.write_csv() && fs::exists(m_dir.results_csv())) {
      m_log << "update: " << m_dir.results_csv() << "\n";
      fs::ofstream out(m_dir.results_csv());
      out << ss.rdbuf();
      out.close();
    }

    return start_i;
  }

  template<typename RunType>
  void MonteDriver<RunType>::single_run(Index cond_index) {

    fs::create_directories(m_dir.conditions_dir(cond_index));

    // perform any requested explicit equilibration passes
    if(m_settings.is_equilibration_passes_each_run()) {

      m_log.write("DoF");
      m_log << "write: " << m_dir.initial_state_runeq_json(cond_index) << "\n" << std::endl;

      jsonParser json;
      to_json(m_mc.configdof(), json).write(m_dir.initial_state_runeq_json(cond_index));
      auto equil_passes = m_settings.equilibration_passes_each_run();

      m_log.begin("Equilibration passes");
      m_log << equil_passes << " equilibration passes\n" << std::endl;

      MonteCounter equil_counter(m_settings, m_mc.steps_per_pass());
      while(equil_counter.pass() != equil_passes) {
        monte_carlo_step(m_mc);
        equil_counter++;
      }
    }

    // initial state (after any equilibriation passes)
    m_log.write("DoF");
    m_log << "write: " << m_dir.initial_state_json(cond_index) << "\n" << std::endl;
    jsonParser json;
    to_json(m_mc.configdof(), json).write(m_dir.initial_state_json(cond_index));

    std::stringstream ss;
    ss << "Conditions " << cond_index;
    m_log.begin(ss.str());
    m_log << std::endl;
    m_log.begin_lap();

    MonteCounter run_counter(m_settings, m_mc.steps_per_pass());
    if(m_enum) {
      m_enum->reset();
    };

    while(true) {

      if(debug()) {
        m_log.custom<Log::debug>("Counter info");
        m_log << "pass: " << run_counter.pass() << "  "
              << "step: " << run_counter.step() << "  "
              << "samples: " << run_counter.samples() << "\n" << std::endl;
      }

      if(m_mc.must_converge()) {

        if(!run_counter.minimums_met()) {

          // keep going, but check for conflicts with maximums
          if(run_counter.maximums_met()) {
            throw std::runtime_error(
              std::string("Error in 'MonteDriver<RunType>::single_run()'\n") +
              "  Conflicting input: Minimum number of passes, steps, or samples not met,\n" +
              "  but maximum number of passes, steps, or samples are met.");
          }
        }
        else {

          if(m_mc.check_convergence_time()) {

            m_log.require<Log::verbose>() << "\n";
            m_log.custom<Log::verbose>("Begin convergence checks");
            m_log << "samples: " << m_mc.sample_times().size() << std::endl;
            m_log << std::endl;

            if(m_mc.is_converged()) {
              break;
            }
          }

          if(run_counter.maximums_met()) {
            break;
          }
        }
      }
      else if(run_counter.is_complete()) {
        // stop
        break;
      }

      bool res = monte_carlo_step(m_mc);

      if(res && m_enum && m_enum->on_accept()) {
        m_enum->insert(m_mc.config());
      }

      run_counter++;

      if(run_counter.sample_time()) {
        if(debug()) {
          m_log.custom<Log::debug>("Sample data");
          m_log << "pass: " << run_counter.pass() << "  "
                << "step: " << run_counter.step() << "  "
                << "take sample " << m_mc.sample_times().size() << "\n" << std::endl;
        }

        m_mc.sample_data(run_counter);
        run_counter.increment_samples();
        if(m_enum && m_enum->on_sample()) {
          m_enum->insert(m_mc.config());
        }
      }
    }
    m_log << std::endl;


    // timing info:
    double s = m_log.lap_time();
    m_log.end(ss.str());
    m_log << "run time: " << s << " (s),  " << s / run_counter.pass() << " (s/pass),  " << s / (run_counter.pass()*run_counter.steps_per_pass() + run_counter.step()) << "(s/step)\n" << std::endl;

    m_log.write("DoF");
    m_log << "write: " << m_dir.final_state_json(cond_index) << "\n" << std::endl;
    to_json(m_mc.configdof(), json).write(m_dir.final_state_json(cond_index));

    m_log.write("Output files");
    m_mc.write_results(cond_index);
    m_log << std::endl;

    if(m_enum) {
      m_enum->save_configs();
    }

    return;
  }

  /**
   * Reads from the settings and constructs an appropriate
   * std::vector of conditions for MonteDriver to visit.
   *
   * Options are:
   * * Single: only visit initial conditions (Returns empty array)
   * * Custom: Provide explicit list of conditions to visit
   * * Incremental: specify initial conditions, final conditions
   * and regular intervals
   */

  template<typename RunType>
  std::vector<typename MonteDriver<RunType>::CondType>
  MonteDriver<RunType>::make_conditions_list(const PrimClex &primclex, const SettingsType &settings) {

    std::vector<CondType> conditions_list;

    switch(m_drive_mode) {

    case Monte::DRIVE_MODE::CUSTOM: {

      // read existing conditions, and check for agreement
      std::vector<CondType> custom_cond(settings.custom_conditions());
      int i = 0;
      while(fs::exists(m_dir.conditions_json(i))) {

        CondType existing;
        jsonParser json(m_dir.conditions_json(i));
        from_json(existing, primclex, json);
        if(existing != custom_cond[i]) {
          m_err_log.error("Conditions mismatch");
          m_err_log << "existing conditions: " << m_dir.conditions_json(i) << "\n";
          m_err_log << existing << "\n\n";
          m_err_log << "specified custom conditions " << i << ":\n";
          m_err_log << custom_cond[i] << "\n" << std::endl;
          throw std::runtime_error("ERROR: custom_conditions list has changed.");
        }
        ++i;

      }
      return custom_cond;
    }

    case Monte::DRIVE_MODE::INCREMENTAL: {

      CondType init_cond(settings.initial_conditions());
      CondType final_cond(settings.final_conditions());
      CondType cond_increment(settings.incremental_conditions());

      CondType incrementing_cond = init_cond;

      int num_increments = 1 + (final_cond - init_cond) / cond_increment;

      for(int i = 0; i < num_increments; i++) {
        conditions_list.push_back(incrementing_cond);
        incrementing_cond += cond_increment;
      }

      int i = 0;
      while(fs::exists(m_dir.conditions_json(i)) && i < conditions_list.size()) {

        CondType existing;
        jsonParser json(m_dir.conditions_json(i));
        from_json(existing, primclex, json);
        if(existing != conditions_list[i]) {
          m_err_log.error("Conditions mismatch");
          m_err_log << "existing conditions: " << m_dir.conditions_json(i) << "\n";
          m_err_log << existing << "\n";
          m_err_log << "incremental conditions " << i << ":\n";
          m_err_log << conditions_list[i] << "\n" << std::endl;
          throw std::runtime_error("ERROR: initial_conditions or incremental_conditions has changed.");
        }
        ++i;

      }

      return conditions_list;
    }

    default: {
      throw std::runtime_error("ERROR: An invalid drive mode was given.");
    }
    }
  }

  template<typename RunType>
  bool monte_carlo_step(RunType &monte_run) {

    typedef typename RunType::EventType EventType;

    const EventType &event = monte_run.propose();

    if(monte_run.check(event)) {
      monte_run.accept(event);
      return true;
    }

    else {
      monte_run.reject(event);
      return false;
    }

  }

}

#endif
