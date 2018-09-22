#include "casm/monte_carlo/MonteSettings.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/MonteCarloEnum.hh"
#include "casm/container/LinearAlgebra.hh"
#include "casm/misc/HallOfFame.hh"

namespace CASM {

  // --- MonteSettings Definitions -------------------------------------------------

  /// \brief Construct MonteSettings by reading a settings JSON file
  ///
  /// - read_path is expected to be within a CASM project directory
  ///
  MonteSettings::MonteSettings(const PrimClex &_primclex, const fs::path &read_path):
    jsonParser(read_path),
    m_primclex(&_primclex) {

    m_root = find_casmroot(fs::absolute(read_path));
    m_output_directory = fs::absolute(read_path).parent_path();

  }


  // --- Project root directory ---------------------------

  fs::path MonteSettings::root() const {
    return m_root;
  }

  const PrimClex &MonteSettings::primclex() const {
    return *m_primclex;
  }


  // --- Type ---------------------------

  /// \brief Return type of Monte Carlo ensemble
  Monte::ENSEMBLE MonteSettings::ensemble() const {
    return _get_setting<Monte::ENSEMBLE>("ensemble", help<Monte::ENSEMBLE>());
  }

  /// \brief Return type of Monte Carlo method
  Monte::METHOD MonteSettings::method() const {
    return _get_setting<Monte::METHOD>("method", help<Monte::METHOD>());
  }

  /// \brief Run in debug mode?
  bool MonteSettings::debug() const {
    auto it = find("debug");
    if(it == end()) {
      return false;
    }
    return it->get<bool>();
  }

  /// \brief Set debug mode
  void MonteSettings::set_debug(bool _debug) {
    (*this)["debug"] = _debug;
  }


  // --- Initialization ---------------------

  /// \brief Returns true if configname of configuration to use as starting motif has been specified
  bool MonteSettings::is_motif_configname() const {
    return _is_setting("driver", "motif", "configname");
  }

  /// \brief Configname of configuration to use as starting motif
  std::string MonteSettings::motif_configname() const {
    std::string help = "string"
                       "  The name of a configuration, of the form \"SCELV_A_B_C_D_E_F/N\"";
    return _get_setting<std::string>("driver", "motif", "configname", help);
  }

  /// \brief Returns true if path to ConfigDoF file to use as starting motif has been specified
  bool MonteSettings::is_motif_configdof() const {
    return _is_setting("driver", "motif", "configdof");
  }

  /// \brief ConfigDoF to use as starting motif
  ConfigDoF MonteSettings::motif_configdof() const {
    std::string help = "string\n"
                       "  Path to file containing DoF, such as an \"final_state.json\" file.";
    fs::path configdof_path = _get_setting<fs::path>("driver", "motif", "configdof", help);
    return jsonParser(configdof_path).get<ConfigDoF>();
  }

  /// \brief Path to ConfigDoF file to use as starting motif
  fs::path MonteSettings::motif_configdof_path() const {
    std::string help = "";
    return _get_setting<fs::path>("driver", "motif", "configdof", help);
  }


  /// \brief Supercell matrix defining the simulation cell
  Eigen::Matrix3i MonteSettings::simulation_cell_matrix() const {
    std::string help = "3x3 transformation matrix, T, such that S = U*T,\n"
                       "  where S, is the supercell lattice vectors,\n"
                       "  and P, is the primitive cell lattice vectors.\n";
    return _get_setting<Eigen::Matrix3i>("supercell", help);
  }


  // --- Driver ---------------------

  /// \brief Given a settings jsonParser figure out the drive mode. Expects drive_mode/incremental,custom
  const Monte::DRIVE_MODE MonteSettings::drive_mode() const {
    return _get_setting<Monte::DRIVE_MODE>("driver", "mode", help<Monte::DRIVE_MODE>());
  }

  /// \brief If dependent runs, start subsequent calculations with the final state
  ///        of the previous calculation. Default true.
  bool MonteSettings::dependent_runs() const {
    if(!_is_setting("driver", "dependent_runs")) {
      return true;
    }
    std::string help = "bool (default=true)\n"
                       "  If true, begin the next calculation with the final DoF from the previous \n"
                       "    calculation.\n"
                       "  If false, begin each calculation with the DoF specified for the \"motif\".\n";
    return _get_setting<bool>("driver", "dependent_runs", help);
  }

  /// \brief Directory where output should go
  const fs::path MonteSettings::output_directory() const {
    return m_output_directory;
  }


  // --- MCData / Sampling ---------------------

  /// \brief Requested confidence level. Default 0.95.
  double MonteSettings::confidence() const {
    std::string help = "number, range (0.0, 1.0), default 0.95)";
    if(_is_setting("data", "confidence")) {
      return _get_setting<double>("data", "confidence", help);
    }
    else {
      return 0.95;
    }
  }

  /// \brief Returns true if snapshots are requested
  bool MonteSettings::write_trajectory() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "write_trajectory";
    std::string help = "bool (default=false)";
    if(!_is_setting(level1, level2, level3)) {
      return false;
    }

    return _get_setting<bool>(level1, level2, level3, help);
  }

  /// \brief Returns true if POSCARs of snapshots are requsted. Requires write_trajectory.
  bool MonteSettings::write_POSCAR_snapshots() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "write_POSCAR_snapshots";
    std::string help = "(bool, default=false)";
    if(!_is_setting(level1, level2, level3)) {
      return false;
    }

    return _get_setting<bool>(level1, level2, level3, help);
  }

  /// \brief Writes all observations
  bool MonteSettings::write_observations() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "write_observations";
    std::string help = "(bool, default=false)";
    if(!_is_setting(level1, level2, level3)) {
      return false;
    }

    return _get_setting<bool>(level1, level2, level3, help);
  }

  /// \brief Write csv versions of files? (csv is the default format if no 'output_format' given)
  bool MonteSettings::write_csv() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "output_format";
    std::string help = "(string or JSON array of string, optional, default='csv')\n"
                       "  Accepts: 'csv' or 'CSV' to write .csv files\n"
                       "           'json' or 'JSON' to write .json files\n"
                       "  Use an array to write in multiple formats.\n";

    if(!_is_setting(level1, level2, level3)) {
      return true;
    }

    const jsonParser &ref = (*this)[level1][level2][level3];

    if(ref.is_array()) {
      std::vector<std::string> formats = _get_setting<std::vector<std::string> >(level1, level2, level3, help);
      for(auto it = formats.begin(); it != formats.end(); ++it) {
        if(*it == "csv" || *it == "CSV") {
          return true;
        }
      }
      return false;
    }
    else {
      std::string input = _get_setting<std::string>(level1, level2, level3, help);
      if(input == "csv" || input == "CSV") {
        return true;
      }
      return false;
    }
  }

  /// \brief Write json versions of files?
  bool MonteSettings::write_json() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "output_format";
    std::string help = "(string or JSON array of string, optional, default='csv')\n"
                       "  Accepts: 'csv' or 'CSV' to write .csv files\n"
                       "           'json' or 'JSON' to write .json files\n"
                       "  Use an array to write in multiple formats.\n";

    if(!_is_setting(level1, level2, level3)) {
      return false;
    }

    const jsonParser &ref = (*this)[level1][level2][level3];

    if(ref.is_array()) {
      std::vector<std::string> formats = _get_setting<std::vector<std::string> >(level1, level2, level3, help);
      for(auto it = formats.begin(); it != formats.end(); ++it) {
        if(*it == "json" || *it == "JSON") {
          return true;
        }
      }
      return false;
    }
    else {
      std::string input = _get_setting<std::string>(level1, level2, level3, help);
      if(input == "json" || input == "JSON") {
        return true;
      }
      return false;
    }
  }


  // --- Enumerating Configurations ---

  /// \brief Returns true if enumeration is requested. (Default false)
  bool MonteSettings::is_enumeration() const {
    if(!_is_setting("data", "enumeration")) {
      return false;
    }
    return true;
  }

  /// \brief Returns 'casm query'-like enumeration metric
  ///
  /// Expects a string containing the Configuration scoring metric. For instance,
  /// 'clex_hull_dist(ALL)' (which is the default).
  ///
  /// Uses boost::lexical_cast<double> on the output to determine the result.
  std::string MonteSettings::enumeration_metric_args() const {

    std::string help = "(string, optional, default=\"clex_hull_dist(ALL)\")\n"
                       "  A 'casm query'-like string that provides a metric for \n"
                       "  ranking configurations as they are encountered during \n"
                       "  a Monte Carlo calculation. The resulting value is used\n"
                       "  to create a hall of fame of 'best' configurations     \n"
                       "  encountered during the calculation. When the          \n"
                       "  calculation is complete configurations in the hall of \n"
                       "  fame are added to the CASM project config list.       \n"
                       "  The 'casm query'-like command should evaluate to a    \n"
                       "  number.";

    if(!_is_setting("data", "enumeration", "metric")) {
      return "clex_hull_dist(ALL)";
    }
    return _get_setting<std::string>("data", "enumeration", "metric", help);
  }

  /// \brief Returns 'casm query'-like enumeration check
  ///
  /// Expects a string that will convert to true if the current Configuration
  /// should be scored for enumeration. For instance, perhaps you only want
  /// to enumerate Configurations that are near or below the currently predicted
  /// convex hull. Then, use 'lt(clex_hull_dist(ALL),0.005)'. Default always
  /// returns true.
  ///
  /// Uses boost::lexical_cast<bool> on the output to determine the result.
  std::string MonteSettings::enumeration_check_args() const {

    std::string help = "(string, optional, default=\"eq(1,1)\")\n"
                       "  A 'casm query'-like string that returns a boolean value \n"
                       "  indicating if (true) a configuration should be considered\n"
                       "  for the enumeration hall of fame.";

    if(!_is_setting("data", "enumeration", "check")) {
      return "eq(1,1)";
    }
    return _get_setting<std::string>("data", "enumeration", "check", help);
  }

  /// \brief Enumeration sample mode (default Monte::ENUM_SAMPLE_MODE::ON_SAMPLE)
  Monte::ENUM_SAMPLE_MODE MonteSettings::enumeration_sample_mode() const {
    if(!_is_setting("data", "enumeration", "sample_mode")) {
      return Monte::ENUM_SAMPLE_MODE::ON_SAMPLE;
    }
    return _get_setting<Monte::ENUM_SAMPLE_MODE>("data", "enumeration", "sample_mode", help<Monte::ENUM_SAMPLE_MODE>());
  }

  /// \brief Insert configurations in their canonical form (default true)
  bool MonteSettings::enumeration_check_existence() const {

    std::string help = "(bool, optional, default=true)\n"
                       "  If true, only configurations that do not already exist\n"
                       "  in the config list are inserted into the enumeration  \n"
                       "  hall of fame.";

    if(!_is_setting("data", "enumeration", "check_existence")) {
      return true;
    }
    return _get_setting<bool>("data", "enumeration", "check_existence", help);
  }

  /// \brief Insert configurations in their canonical form (default true)
  bool MonteSettings::enumeration_insert_canonical() const {

    std::string help = "(bool, optional, default=true)\n"
                       "  If true, configurations are inserted into the         \n"
                       "  enumeration hall of fame in their canonical form. If  \n"
                       "  'check_existence' is true, this must be set to true.";

    bool val;
    if(!_is_setting("data", "enumeration", "insert_canonical")) {
      val = true;
    }
    else {
      val = _get_setting<bool>("data", "enumeration", "insert_canonical", help);
    }

    // if check_existence, insert_canonical must be true
    if(enumeration_check_existence() && !val) {
      throw std::runtime_error(
        "Error in Monte Carlo enumeration in settings: "
        "If 'check_existence' is true, then 'insert_canonical' must be true"
      );
    }
    return val;
  }

  /// \brief Returns enumeration halloffame max size (default 100)
  Index MonteSettings::enumeration_N_halloffame() const {

    std::string help = "(integer, optional, default=100)\n"
                       "  The number of configurations that are allowed in the \n"
                       "  enumeration hall of fame.";

    if(!_is_setting("data", "enumeration", "N_halloffame")) {
      return 100;
    }
    return _get_setting<Index>("data", "enumeration", "N_halloffame", help);
  }

  /// \brief Returns enumeration halloffame tolerance (default 1e-8)
  double MonteSettings::enumeration_tol() const {

    std::string help = "(number, optional, default=1e-8)\n"
                       "  Tolerance used for floating point comparison of         \n"
                       "  configuration scores in the enumeration hall of fame.";

    if(!_is_setting("data", "enumeration", "tolerance")) {
      return 1e-8;
    }
    return _get_setting<double>("data", "enumeration", "tolerance", help);
  }


  /// \brief Returns true if (*this)[level1].contains(level2)
  bool MonteSettings::_is_setting(std::string level1, std::string level2) const {
    try {
      return (*this)[level1].contains(level2);
    }

    catch(std::runtime_error &e) {
      Log &err_log = default_err_log();
      err_log.error<Log::standard>("Reading Monte Carlo settings");
      err_log << "No [\"" << level1 << "\"] setting found.\n" << std::endl;
      throw e;
    }
  }

  /// \brief Returns true if (*this)[level1][level2].contains(level3)
  bool MonteSettings::_is_setting(std::string level1, std::string level2, std::string level3) const {
    try {
      return (*this)[level1][level2].contains(level3);
    }

    catch(std::runtime_error &e) {
      Log &err_log = default_err_log();
      err_log.error<Log::standard>("Reading Monte Carlo settings");
      if(this->contains(level1)) {
        err_log << "No [\"" << level1 << "\"][\"" << level2 << "\"] setting found" << std::endl;
      }
      else {
        err_log << "No [\"" << level1 << "\"] setting found" << std::endl;
      }
      throw;
    }
  }



  // --- EquilibriumMonteSettings Definitions -------------------------------------------------


  // --- MCData / Sampling ---------------------

  /// \brief Sample by pass?
  bool EquilibriumMonteSettings::sample_by_pass() const {
    return Monte::SAMPLE_MODE::PASS == _get_setting<Monte::SAMPLE_MODE>("data", "sample_by", help<Monte::SAMPLE_MODE>());
  }

  /// \brief Sample by step?
  bool EquilibriumMonteSettings::sample_by_step() const {
    return Monte::SAMPLE_MODE::STEP == _get_setting<Monte::SAMPLE_MODE>("data", "sample_by", help<Monte::SAMPLE_MODE>());
  }

  /// \brief Figure out how often to take samples
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::sample_period() const {

    std::string level1 = "data";
    std::string level2 = "sample_period";
    std::string help = "(int, default=1)\n"
                       "  In conjunction with \"sample_by\", determines how often to make observations.";
    if(!_is_setting(level1, level2)) {
      return 1;
    }

    return _get_setting<size_type>(level1, level2, help);
  }

  /// \brief Returns true if explicit equilibration passes for the first run have been specified
  bool EquilibriumMonteSettings::is_equilibration_passes_first_run() const {
    return _is_setting("data", "equilibration_passes_first_run");
  }

  /// \brief Number of explicit equilibration passes requsted for the first run
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::equilibration_passes_first_run() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "equilibration_passes_first_run", help);
  }

  /// \brief Returns true if explicit equilibration passes for each run have been specified
  bool EquilibriumMonteSettings::is_equilibration_passes_each_run() const {
    return _is_setting("data", "equilibration_passes_each_run");
  }

  /// \brief Number of explicit equilibration passes requsted for each run
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::equilibration_passes_each_run() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "equilibration_passes_each_run", help);
  }


  /// \brief Returns true if the number of passes has been specified
  bool EquilibriumMonteSettings::is_N_pass() const {
    return _is_setting("data", "N_pass");
  }

  /// \brief Returns the number of passes requested
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::N_pass() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "N_pass", help);
  }

  /// \brief Returns true if the number of steps has been specified
  bool EquilibriumMonteSettings::is_N_step() const {
    std::string help = "int (optional)";
    return _is_setting("data", "N_step");
  }

  /// \brief Returns the number of steps requested
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::N_step() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "N_step", help);
  }

  /// \brief Returns true if the number of samples has been specified
  bool EquilibriumMonteSettings::is_N_sample() const {
    return _is_setting("data", "N_sample");
  }

  /// \brief Returns the number of samples requested
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::N_sample() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "N_sample", help);
  }


  /// \brief Returns true if a maximum number of passes has been specified
  bool EquilibriumMonteSettings::is_max_pass() const {
    return _is_setting("data", "max_pass");
  }

  /// \brief Maximum number of passes, required if sample by pass
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::max_pass() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "max_pass", help);
  }

  /// \brief Returns true if a minimum number of passes has been specified
  bool EquilibriumMonteSettings::is_min_pass() const {
    return _is_setting("data", "min_pass");
  }

  /// \brief Minimum number of passes
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::min_pass() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "min_pass", help);
  }


  /// \brief Returns true if a maximum number of steps has been specified
  bool EquilibriumMonteSettings::is_max_step() const {
    return _is_setting("data", "max_step");
  }

  /// \brief Maximum number of steps
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::max_step() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "max_step", help);
  }

  /// \brief Returns true if a minimum number of steps has been specified
  bool EquilibriumMonteSettings::is_min_step() const {
    return _is_setting("data", "min_step");

  }

  /// \brief Minimum number of steps
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::min_step() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "min_step", help);
  }


  /// \brief Returns true if a maximum number of samples has been specified
  bool EquilibriumMonteSettings::is_max_sample() const {
    return _is_setting("data", "max_sample");
  }

  /// \brief Maximum number of steps
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::max_sample() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "max_sample", help);
  }

  /// \brief Returns true if a minimum number of sample has been specified
  bool EquilibriumMonteSettings::is_min_sample() const {
    return _is_setting("data", "min_sample");
  }

  /// \brief Minimum number of steps, default 0
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::min_sample() const {
    std::string help = "int (optional)";
    return _get_setting<size_type>("data", "min_sample", help);
  }


  // --- Data ---------------------

  /// \brief Figure out how large data containers should be
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::max_data_length() const {
    if(!_is_setting("data", "sample_by")) {
      return 1024;
    }
    else if(sample_by_pass()) {
      if(is_max_pass()) {
        return (max_pass() / sample_period());
      }
      else if(is_N_pass()) {
        return (N_pass() / sample_period());
      }
      else if(is_N_sample()) {
        return N_sample();
      }
      else {
        return 1024;
      }
    }
    else if(sample_by_step()) {
      if(is_max_step()) {
        return (max_step() / sample_period());
      }
      else if(is_N_step()) {
        return (N_step() / sample_period());
      }
      else if(is_N_sample()) {
        return N_sample();
      }
      else {
        return 1024;
      }
    }
    else {
      Log &err_log = default_err_log();
      err_log.error<Log::standard>("Reading Monte Carlo settings");
      err_log << "Could not determine max data length.\n";
      err_log << "Please check 'sample_by', 'max_pass' or 'max_step', and 'sample_period' in your input file.\n" << std::endl;
      throw std::runtime_error(std::string("Error determining max_data_length"));
    }

  }


}
