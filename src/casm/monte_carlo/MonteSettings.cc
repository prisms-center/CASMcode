#include "casm/monte_carlo/MonteSettings.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/container/LinearAlgebra.hh"
#include <boost/lexical_cast.hpp>

namespace CASM {

  MonteSettings::MonteSettings(const fs::path &read_path):
    jsonParser(read_path) {
    (*this)["parent_path"] = read_path.string();    //TODO: how to make this work with write() cleverly
  }
  
  // --- Type ---------------------------
    
  /// \brief Return type of Monte Carlo calculation
  Monte::TYPE MonteSettings::type() const {
    try {
      return Monte::monte_type((*this)["type"].get<std::string>());
    }
    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::type" << std::endl;
      std::cerr << "Expected [\"type\"]" << std::endl;
      throw e;
    }
  }
  
  
  // --- Initialization ---------------------
    
  /// \brief Configname of configuration to use as starting motif
  std::string MonteSettings::motif_configname() const {
    std::string level1 = "initialization";
    std::string level2 = "motif";
    std::string level3 = "configname";

    try {
      return (*this)[level1][level2][level3].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::motif_config_index" << std::endl;
      std::cerr << "Expected [" << level1 << "][" << level2 << "][" << level3 << "]" << std::endl;
      throw e;
    }
  }
  
  /// \brief Supercell matrix defining the simulation cell
  Eigen::Matrix3i MonteSettings::simulation_cell_matrix() const {
    try {
      return (*this)["initialization"]["matrix"].get<Eigen::Matrix3i>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::simulation_cell_matrix" << std::endl;
      std::cerr << "No matrix given to construct Monte Carlo cell." << std::endl;
      std::cerr << "Expected [\"initialization\"][\"matrix\"]" << std::endl;
      throw e;
    }
  }

  
  // --- Driver ---------------------
    
  /// \brief Given a settings jsonParser figure out the drive mode. Expects drive_mode/single,custom,incremental
  const Monte::DRIVE_MODE MonteSettings::drive_mode() const {
    Monte::DRIVE_MODE dmode;

    std::string dmodestring;

    try {
      dmodestring = (*this)["driver"]["mode"].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::drive_mode" << std::endl;
      std::cerr << "No mode for driver found!" << std::endl;
      std::cerr << "Expected [\"driver\"][\"mode\"]" << std::endl;
      throw e;
    }

    if(dmodestring == "single") {
      dmode = Monte::DRIVE_MODE::SINGLE;
    }

    else if(dmodestring == "incremental") {
      dmode = Monte::DRIVE_MODE::INCREMENTAL;
    }

    else {
      throw std::runtime_error(
        std::string("Error in MonteSettings::drive_mode.\n") +
                    "  Found [\"driver\"][\"mode\"] = \"" + dmodestring + "\", but allowed options are 'single' or 'incremental'.");
    }

    return dmode;
  }

  
  // --- Conditions settings ---------------------
  
  /// \brief Expects initial_conditions
  const MonteSettings &MonteSettings::initial_conditions() const {

    std::string level1 = "driver";
    std::string level2 = "initial_conditions";
    try {
      m_pcondition_subtree = &(*this)[level1][level2];
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::initial_conditions" << std::endl;
      std::cerr << "Tried to point at [\"" << level1 << "\"][\"" << level2 << "\"] but nothing was there!" << std::endl;
      throw e;    //TODO: Make a MonteSettings exception class?
    }

    return *this;
  }

  /// \brief Expects final_conditions
  const MonteSettings &MonteSettings::final_conditions() const {

    std::string level1 = "driver";
    std::string level2 = "final_conditions";
    try {
      m_pcondition_subtree = &(*this)[level1][level2];
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::final_conditions" << std::endl;
      std::cerr << "Tried to point at [\"" << level1 << "\"][\"" << level2 << "\"] but nothing was there!" << std::endl;
      throw e;    //TODO: Make a MonteSettings exception class?
    }

    return *this;
  }

  /// \brief Expects incremental_conditions
  const MonteSettings &MonteSettings::incremental_conditions() const {

    std::string level1 = "driver";
    std::string level2 = "incremental_conditions";
    try {
      m_pcondition_subtree = &(*this)[level1][level2];
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::incremental_conditions" << std::endl;
      std::cerr << "Tried to point at [\"" << level1 << "\"][\"" << level2 << "\"] but nothing was there!" << std::endl;
      throw e;    //TODO: Make a MonteSettings exception class?
    }

    return *this;
  }

  /// \brief Given a settings jsonParser figure out the chemical potential for a particular component. Expects mu/'a'/value
  /**
   * Because there can be many different conditions, this routine expects a subset of
   * the overall settings MonteSettings. You'll have to pass something like mymonteio.initial_conditions()
   * to the constructor of whatever conditions you're making.
   */
  double MonteSettings::chemical_potential(std::string comp_name) const {
    
    try {
      return (*m_pcondition_subtree)["mu"][comp_name].get<double>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::chemical_potential" << std::endl;
      std::cerr << "Looking for chemical_potential for composition: " << comp_name << std::endl;
      std::cerr << "No chemical potential found: \n" << *m_pcondition_subtree << std::endl;
      std::cerr << "Expected [\"mu\"][" << comp_name << "]" << std::endl;
      throw e;
    }
  }

  /// \brief Given a settings jsonParser figure out the temperature of the run. Expects temperature/value
  /**
   * Because there can be many different conditions, this routine expects a subset of
   * the overall settings MonteSettings. You'll have to pass something like mymonteio.initial_conditions()
   * to the constructor of whatever conditions you're making.
   */
  double MonteSettings::temperature() const {
    try {
      return (*m_pcondition_subtree)["temperature"].get<double>();
    }

    catch(std::runtime_error &e) {
      throw std::runtime_error(
        std::string("Error in MonteSettings::temperature.\n") +
                    "  No [\"temperature\"] inside requested conditions.");
    }
  }

  /// \brief Given a settings jsonParser figure out the global tolerance (probably for == operator). Expects tolerance/value
  /**
   * Because there can be many different conditions, this routine expects a subset of
   * the overall settings MonteSettings. You'll have to pass something like mymonteio.initial_conditions()
   * to the constructor of whatever conditions you're making.
   */
  double MonteSettings::tolerance() const {
    try {
      return (*m_pcondition_subtree)["tolerance"].get<double>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "WARNING in Monte::global_tolerance" << std::endl;
      std::cerr << "Expected [\"tolerance\"]" << std::endl;
      std::cerr << "The global tolerance " << TOL << " will be used." << std::endl;
      return TOL;
    }
  }
  

  // --- Project settings ---------------------
  
  /// \brief Given a settings jsonParser figure out what the project clex settings to use are:
  std::string MonteSettings::clex() const {
    std::string level1 = "initialization";
    std::string level2 = "clex";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::clex" << std::endl;
      std::cerr << "No clex settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }
  }

  /// \brief Given a settings jsonParser figure out what the project bset settings to use are:
  std::string MonteSettings::bset() const {
    std::string level1 = "initialization";
    std::string level2 = "bset";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::bset" << std::endl;
      std::cerr << "No bset settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }
  }

  /// \brief Given a settings jsonParser figure out what the project calctype settings to use are:
  std::string MonteSettings::calctype() const {
    std::string level1 = "initialization";
    std::string level2 = "calctype";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::calctype" << std::endl;
      std::cerr << "No calctype settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }
  }

  /// \brief Given a settings jsonParser figure out what the project ref settings to use are:
  std::string MonteSettings::ref() const {
    std::string level1 = "initialization";
    std::string level2 = "ref";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::ref" << std::endl;
      std::cerr << "No ref settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }

  }

  /// \brief Given a settings jsonParser figure out what the project eci settings to use are:
  std::string MonteSettings::eci() const {
    std::string level1 = "initialization";
    std::string level2 = "eci";
    try {
      return (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::eci" << std::endl;
      std::cerr << "No eci settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }

  }

  /// \brief Directory where output should go
  const fs::path MonteSettings::output_directory() const {
    fs::path return_path;
    try {
      return_path = (*this)["data"]["storage"]["output_directory"].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::output_directory" << std::endl;
      std::cerr << "No output path found" << std::endl;
      std::cerr << "Expected [\"data\"][\"storage\"][\"output_directory\"]" << std::endl;
      throw e;
    }

    return return_path;
  }
  
  
  // --- MCData / Sampling ---------------------
  
  /// \brief Sample by pass?
  bool MonteSettings::sample_by_pass() const {

    std::string sample_by;
    std::string level1 = "data";
    std::string level2 = "sample_by";
    try {
      sample_by = (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::sample_by_pass" << std::endl;
      std::cerr << "No 'sample_by' settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }

    if(sample_by == "pass") {
      return true;
    }
    else if(sample_by == "step") {
      return false;
    }
    else {
      std::cerr << "ERROR in Monte::sample_by_pass" << std::endl;
      std::cerr << "Unexpected 'sample_py' setting: '" << sample_by << "'. Please use \"pass\" or \"step\"." << std::endl;
      throw std::runtime_error(std::string("Unexpected 'sample_by' setting: ") + sample_by);
    }
  }

  /// \brief Sample by step?
  bool MonteSettings::sample_by_step() const {

    std::string sample_by;
    std::string level1 = "data";
    std::string level2 = "sample_by";
    try {
      sample_by = (*this)[level1][level2].get<std::string>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::sample_by_step" << std::endl;
      std::cerr << "No 'sample_by' settings found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      throw e;
    }

    if(sample_by == "step") {
      return true;
    }
    else if(sample_by == "pass") {
      return false;
    }
    else {
      std::cerr << "ERROR in Monte::sample_by_pass" << std::endl;
      std::cerr << "Unexpected 'sample_py' setting: '" << sample_by << "'. Please use \"pass\" or \"step\"." << std::endl;
      throw std::runtime_error(std::string("Unexpected 'sample_by' setting: ") + sample_by);
    }
  }

  /// \brief Figure out how often to take samples
  MonteSettings::size_type MonteSettings::sample_period() const {

    size_type value = 1;
    std::string level1 = "data";
    std::string level2 = "sample_period";

    if((*this)[level1].contains(level2)) {
      try {
        value = (*this)[level1]["sample_period"].get<size_type>();
      }

      catch(std::runtime_error &e) {
        throw std::runtime_error(
          std::string("Error in MonteSettings::sample_period\n") +
                      "  Expected [" + level1 + "][" + level2 + "]");
      }
    }

    return value;
  }
  
  
  /// \brief Requested confidence level. Default 0.95.
  double MonteSettings::confidence() const {
    if( _is_setting("data", "confidence")) {
      return _get_setting<double>("data", "confidence");
    }
    else {
      return 0.95;
    }
  }

  
  /// \brief Returns true if (*this)[level1].contains(level2)
  bool MonteSettings::_is_setting(std::string level1, std::string level2) const {
    try {
      return (*this)[level1].contains(level2);
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::is_" << level2 << std::endl;
      std::cerr << "No '" << level1 << "' setting found" << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"]" << std::endl;
      throw e;
    }
  }
  

  /// \brief Returns true if explicit equilibration passes for the first run have been specified
  bool MonteSettings::is_equilibration_passes_first_run() const {
    return _is_setting("data", "equilibration_passes_first_run");
  }
  
  /// \brief Number of explicit equilibration passes requsted for the first run
  MonteSettings::size_type MonteSettings::equilibration_passes_first_run() const {
    return _get_setting<size_type>("data", "equilibration_passes_first_run");
  }
  
  /// \brief Returns true if explicit equilibration passes for each run have been specified
  bool MonteSettings::is_equilibration_passes_each_run() const {
    return _is_setting("data", "equilibration_passes_each_run");
  }
  
  /// \brief Number of explicit equilibration passes requsted for each run
  MonteSettings::size_type MonteSettings::equilibration_passes_each_run() const {
    return _get_setting<size_type>("data", "equilibration_passes_each_run");
  }
  
  
  /// \brief Returns true if the number of passes has been specified
  bool MonteSettings::is_N_pass() const {
    return _is_setting("data", "N_pass");
  }
  
  /// \brief Returns the number of passes requested
  MonteSettings::size_type MonteSettings::N_pass() const {
    return _get_setting<size_type>("data", "N_pass");
  }
  
  /// \brief Returns true if the number of steps has been specified
  bool MonteSettings::is_N_step() const {
    return _is_setting("data", "N_step");
  }

  /// \brief Returns the number of steps requested
  MonteSettings::size_type MonteSettings::N_step() const {
    return _get_setting<size_type>("data", "N_step");
  }
  
  /// \brief Returns true if the number of samples has been specified
  bool MonteSettings::is_N_sample() const {
    return _is_setting("data", "N_sample");
  }

  /// \brief Returns the number of samples requested
  MonteSettings::size_type MonteSettings::N_sample() const {
    return _get_setting<size_type>("data", "N_sample");
  }


  /// \brief Returns true if a maximum number of passes has been specified
  bool MonteSettings::is_max_pass() const {
    return _is_setting("data", "max_pass");
  }
  
  /// \brief Maximum number of passes, required if sample by pass
  MonteSettings::size_type MonteSettings::max_pass() const {
    return _get_setting<size_type>("data", "max_pass");
  }

  /// \brief Returns true if a minimum number of passes has been specified
  bool MonteSettings::is_min_pass() const {
    return _is_setting("data", "min_pass");
  }
  
  /// \brief Minimum number of passes
  MonteSettings::size_type MonteSettings::min_pass() const {
    return _get_setting<size_type>("data", "min_pass");
  }

  
  /// \brief Returns true if a maximum number of steps has been specified
  bool MonteSettings::is_max_step() const {
    return _is_setting("data", "max_step");
  }
  
  /// \brief Maximum number of steps
  MonteSettings::size_type MonteSettings::max_step() const {
    return _get_setting<size_type>("data", "max_step");
  }
  
  /// \brief Returns true if a minimum number of steps has been specified
  bool MonteSettings::is_min_step() const {
    return _is_setting("data", "min_step");
  
  }
  
  /// \brief Minimum number of steps
  MonteSettings::size_type MonteSettings::min_step() const {
    return _get_setting<size_type>("data", "min_step");
  }

  
  /// \brief Returns true if a maximum number of samples has been specified
  bool MonteSettings::is_max_sample() const {
    return _is_setting("data", "max_sample");
  }
  
  /// \brief Maximum number of steps
  MonteSettings::size_type MonteSettings::max_sample() const {
    return _get_setting<size_type>("data", "max_sample");
  }
  
  /// \brief Returns true if a minimum number of sample has been specified
  bool MonteSettings::is_min_sample() const {
    return _is_setting("data", "min_sample");
  }
  
  /// \brief Minimum number of steps, default 0
  MonteSettings::size_type MonteSettings::min_sample() const {
    return _get_setting<size_type>("data", "min_sample");
  }
  
  
  /// \brief Returns true if snapshots are requested
  bool MonteSettings::write_trajectory() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "write_trajectory";
    try {
      
      if(!(*this)[level1][level2].contains(level3)) {
        return false;
      }
      
      return (*this)[level1][level2][level3].get<bool>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::write_trajectory" << std::endl;
      std::cerr << "[\"data\"][\"storage\"] must exist" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"write_trajectory\"] exists, it must be a bool" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"write_trajectory\"] does not exist, default is false (do not write trajectory) " << std::endl;
      throw e;
    }
  }
  
  /// \brief Returns true if POSCARs of snapshots are requsted. Requires write_trajectory.
  bool MonteSettings::write_POSCAR_snapshots() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "write_POSCAR_snapshots";
    try {
      
      if(!(*this)[level1][level2].contains(level3)) {
        return false;
      }
      
      return (*this)[level1][level2][level3].get<bool>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::write_POSCAR_snapshots" << std::endl;
      std::cerr << "[\"data\"][\"storage\"] must exist" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"write_POSCAR_snapshots\"] exists, it must be a bool" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"write_POSCAR_snapshots\"] does not exist, default is false (do not write POSCAR snapshots) " << std::endl;
      throw e;
    }
  }
  
  /// \brief Writes all observations
  bool MonteSettings::write_observations() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "write_observations";
    try {
      
      if(!(*this)[level1][level2].contains(level3)) {
        return false;
      }
      
      return (*this)[level1][level2][level3].get<bool>();
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::write_observations" << std::endl;
      std::cerr << "[\"data\"][\"storage\"] must exist" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"write_observations\"] exists, it must be a bool" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"write_observations\"] does not exist, default is false (do not write observations) " << std::endl;
      throw e;
    }
  }
  
  /// \brief Write csv versions of files? (csv is the default format if no 'output_format' given)
  bool MonteSettings::write_csv() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "output_format";
    try {
      
      if(!(*this)[level1][level2].contains(level3)) {
        return true;
      }
      
      const jsonParser& ref = (*this)[level1][level2][level3];
      
      if(ref.is_string()) {
        std::string input = ref.get<std::string>();
        if(input == "csv" || input == "CSV") {
          return true;
        }
        return false;
      }
      else if(ref.is_array()) {
        for(auto it = ref.cbegin(); it != ref.cend(); ++it) {
          std::string input = it->get<std::string>();
          if(input == "csv" || input == "CSV") {
            return true;
          }
        }
        return false;
      }
      
      throw std::runtime_error(
        std::string("ERROR in 'MonteSettings::write_csv()'\n") +
                    "  Expected [\"data\"][\"storage\"][\"output_format\"] to contain a string or array of strings.");
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::write_csv" << std::endl;
      std::cerr << "[\"data\"][\"storage\"] must exist" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"output_format\"] exists, it may be a string or an array of strings" << std::endl;
      std::cerr << "  if any of those strings is \"csv\" or \"CSV\", return true (do write .csv files)" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"output_format\"] does not exist, default is true (do write .csv files)" << std::endl;
      throw e;
    }
  } 
  
  /// \brief Write json versions of files?
  bool MonteSettings::write_json() const {
    std::string level1 = "data";
    std::string level2 = "storage";
    std::string level3 = "output_format";
    try {
      
      if(!(*this)[level1][level2].contains(level3)) {
        return false;
      }
      
      const jsonParser& ref = (*this)[level1][level2][level3];
      
      if(ref.is_string()) {
        std::string input = ref.get<std::string>();
        if(input == "json" || input == "JSON") {
          return true;
        }
        return false;
      }
      else if(ref.is_array()) {
        for(auto it = ref.cbegin(); it != ref.cend(); ++it) {
          std::string input = it->get<std::string>();
          if(input == "json" || input == "JSON") {
            return true;
          }
        }
        return false;
      }
      
      throw std::runtime_error(
        std::string("ERROR in 'MonteSettings::write_json()'\n") +
                    "  Expected [\"data\"][\"storage\"][\"output_format\"] to contain a string or array of strings.");
    }

    catch(std::runtime_error &e) {
      std::cerr << "ERROR in MonteSettings::write_json" << std::endl;
      std::cerr << "[\"data\"][\"storage\"] must exist" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"output_format\"] exists, it may be a string or an array of strings" << std::endl;
      std::cerr << "  if any of those strings is \"json\" or \"JSON\", return true (do write .json files)" << std::endl;
      std::cerr << "if [\"data\"][\"storage\"][\"output_format\"] does not exist, default is false (do not write .json files)" << std::endl;
      throw e;
    }
  } 
  
  
  // --- Data ---------------------
    
  /// \brief Figure out how large data containers should be
  MonteSettings::size_type MonteSettings::max_data_length() const {
    try {
      if(sample_by_pass()) {
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
        throw std::runtime_error(std::string("Error in MonteSettings::max_data_length()"));
      }
    }
    catch(...) {
      std::cerr << "ERROR: Could not get max data length." << std::endl;
      std::cerr << "Please check 'sample_by', 'max_pass' or 'max_step', and 'sample_period' in your input file" << std::endl;
      throw;
    }

  }

  
  /**
   * This will create a jsonParser object that can be used as a template.
   * I'm writing this mostly for testing purposes but it could come in handy later.
   *
   * I think this routine should get more things appended to it as different classes
   * require different things. Anything not needed in the initialization will get
   * ignored anyway.
   */

  jsonParser example_base_json_settings() {
    
    jsonParser example_settings;
    
    // ---- Initialization settings --------------------
    
    //Seed configuration
    example_settings["initialization"]["motif"]["configname"] = "SCELV_A_B_C_D_E_F/X";

    //Monte Carlo cell
    Eigen::Matrix3i tmat = Eigen::Matrix3i::Zero();
    tmat(0,0) = 10;
    tmat(1,1) = 10;
    tmat(2,2) = 10;
    example_settings["initialization"]["matrix"] = tmat;

    //Cluster expanstion to use, project clex name
    example_settings["initialization"]["clex"] = "formation_energy";

    //Basis set to use, project bset name
    example_settings["initialization"]["bset"] = "default";

    //Calctype settings to use, project calctype name
    example_settings["initialization"]["calctype"] = "default";

    //Calctype settings to use, project ref name
    example_settings["initialization"]["ref"] = "ref";

    //ECI settings name, project eci name
    example_settings["initialization"]["eci"] = "default";
    
    
    // ---- Data settings --------------------
    
    example_settings["data"]["sample_by"] = "pass";
    example_settings["data"]["sample_period"] = 1;
    
    example_settings["data"]["max_pass"] = 10000;
    example_settings["data"]["min_pass"] = 0;
    
    // also allowed:
    //example_settings["data"]["equilibration_passes_first_run"] = 2000;
    //example_settings["data"]["equilibration_passes_each_run"] = 1000;
    
    //example_settings["data"]["max_step"] = 10000;
    //example_settings["data"]["min_step"] = 0;
    
    // also allowed:
    //example_settings["data"]["max_sample"] = 10000;
    //example_settings["data"]["min_sample"] = 0;
    
    // confidence level (default 0.95)
    example_settings["data"]["confidence"] = 0.95;
    
    example_settings["data"]["measurements"] = jsonParser::array(5);
    for(int i=0; i<5; i++) {
      example_settings["data"]["measurements"].push_back(jsonParser::object());
    }
    example_settings["data"]["measurements"][0]["quantity"] = "formation_energy";
    example_settings["data"]["measurements"][1]["quantity"] = "generalized_enthalpy";
    example_settings["data"]["measurements"][2]["quantity"] = "mol_composition";
    example_settings["data"]["measurements"][3]["quantity"] = "formation_energy";
    example_settings["data"]["measurements"][4]["quantity"] = "formation_energy";
    
    // each of the above measurements can also have a requested precision, in which
    // case the monte carlo calculation will run until that precision is reached (<X> +/- prec), 
    // given the confidence level
    // example_settings["data"]["measurements"][0]["precision"] = 0.01
    
    example_settings["data"]["storage"]["output_directory"] = "./mc_results";
    example_settings["data"]["storage"]["write_observations"] = true;
    example_settings["data"]["storage"]["write_trajectory"] = true;
    example_settings["data"]["storage"]["write_POSCAR_snapshots"] = true;
    
    // output file format
    //
    // - may be single string or array of strings
    // - csv is the default format if no 'output_format' given)
    // - options are "csv"/"CSV" or "json"/"JSON"
    example_settings["data"]["storage"]["output_format"] = jsonParser::array();
    example_settings["data"]["storage"]["output_format"].push_back("csv"); // or "CSV"
    example_settings["data"]["storage"]["output_format"].push_back("json"); // or "JSON"
    
    
    
    // ---- Driver settings -------------------
    
    example_settings["driver"]["mode"] = "incremental";
    
    example_settings["driver"]["initial_conditions"]["mu"]["a"] = 2.0;
    example_settings["driver"]["initial_conditions"]["temperature"] = 100.0; // K
    example_settings["driver"]["initial_conditions"]["tolerance"] = 0.001;
    
    example_settings["driver"]["final_conditions"]["mu"]["a"] = 2.0;
    example_settings["driver"]["final_conditions"]["temperature"] = 1000.0; // K
    example_settings["driver"]["final_conditions"]["tolerance"] = 0.001; 
    
    example_settings["driver"]["incremental_conditions"]["mu"]["a"] = 0.0;
    example_settings["driver"]["incremental_conditions"]["temperature"] = 10.0; // K
    example_settings["driver"]["incremental_conditions"]["tolerance"] = 0.001; 
    
    return example_settings;
  }


}



