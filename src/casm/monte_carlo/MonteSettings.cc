#include "casm/monte_carlo/MonteSettings.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/container/LinearAlgebra.hh"

namespace CASM {

  // --- MonteSettings Definitions -------------------------------------------------
  
  /// \brief Construct MonteSettings by reading a settings JSON file
  ///
  /// - read_path is expected to be within a CASM project directory
  ///
  MonteSettings::MonteSettings(const fs::path &read_path):
  jsonParser(read_path) {
    
    m_root = find_casmroot(fs::absolute(read_path));
    m_output_directory = fs::absolute(read_path).parent_path();
    
  }
  
  
  // --- Project root directory ---------------------------
    
  fs::path MonteSettings::root() const {
    return m_root;
  }
  
  
  // --- Type ---------------------------
    
  /// \brief Return type of Monte Carlo ensemble
  Monte::ENSEMBLE MonteSettings::ensemble() const {
    try {
      return Monte::monte_ensemble((*this)["ensemble"].get<std::string>());
    }
    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::ensemble" << std::endl;
      std::cerr << "Expected [\"ensemble\"]" << std::endl;
      throw e;
    }
  }
  
  /// \brief Return type of Monte Carlo method
  Monte::METHOD MonteSettings::method() const {
    try {
      return Monte::monte_method((*this)["method"].get<std::string>());
    }
    catch(std::runtime_error &e) {
      std::cerr << "ERROR in Monte::method" << std::endl;
      std::cerr << "Expected [\"method\"]" << std::endl;
      throw e;
    }
  }
  
  /// \brief Run in debug mode?
  bool MonteSettings::debug() const {
    auto it = find("debug");
    if(it == end()) {
      return false;
    }
    return it->get<bool>();
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

  /// \brief Directory where output should go
  const fs::path MonteSettings::output_directory() const {
    return m_output_directory;
  }
  
  
  // --- MCData / Sampling ---------------------
  
  /// \brief Requested confidence level. Default 0.95.
  double MonteSettings::confidence() const {
    if( _is_setting("data", "confidence")) {
      return _get_setting<double>("data", "confidence");
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
  

  
  // --- EquilibriumMonteSettings Definitions -------------------------------------------------
  
  
  // --- MCData / Sampling ---------------------
  
  /// \brief Sample by pass?
  bool EquilibriumMonteSettings::sample_by_pass() const {

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
  bool EquilibriumMonteSettings::sample_by_step() const {

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
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::sample_period() const {

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
  
  
  /// \brief Returns true if explicit equilibration passes for the first run have been specified
  bool EquilibriumMonteSettings::is_equilibration_passes_first_run() const {
    return _is_setting("data", "equilibration_passes_first_run");
  }
  
  /// \brief Number of explicit equilibration passes requsted for the first run
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::equilibration_passes_first_run() const {
    return _get_setting<size_type>("data", "equilibration_passes_first_run");
  }
  
  /// \brief Returns true if explicit equilibration passes for each run have been specified
  bool EquilibriumMonteSettings::is_equilibration_passes_each_run() const {
    return _is_setting("data", "equilibration_passes_each_run");
  }
  
  /// \brief Number of explicit equilibration passes requsted for each run
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::equilibration_passes_each_run() const {
    return _get_setting<size_type>("data", "equilibration_passes_each_run");
  }
  
  
  /// \brief Returns true if the number of passes has been specified
  bool EquilibriumMonteSettings::is_N_pass() const {
    return _is_setting("data", "N_pass");
  }
  
  /// \brief Returns the number of passes requested
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::N_pass() const {
    return _get_setting<size_type>("data", "N_pass");
  }
  
  /// \brief Returns true if the number of steps has been specified
  bool EquilibriumMonteSettings::is_N_step() const {
    return _is_setting("data", "N_step");
  }

  /// \brief Returns the number of steps requested
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::N_step() const {
    return _get_setting<size_type>("data", "N_step");
  }
  
  /// \brief Returns true if the number of samples has been specified
  bool EquilibriumMonteSettings::is_N_sample() const {
    return _is_setting("data", "N_sample");
  }

  /// \brief Returns the number of samples requested
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::N_sample() const {
    return _get_setting<size_type>("data", "N_sample");
  }


  /// \brief Returns true if a maximum number of passes has been specified
  bool EquilibriumMonteSettings::is_max_pass() const {
    return _is_setting("data", "max_pass");
  }
  
  /// \brief Maximum number of passes, required if sample by pass
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::max_pass() const {
    return _get_setting<size_type>("data", "max_pass");
  }

  /// \brief Returns true if a minimum number of passes has been specified
  bool EquilibriumMonteSettings::is_min_pass() const {
    return _is_setting("data", "min_pass");
  }
  
  /// \brief Minimum number of passes
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::min_pass() const {
    return _get_setting<size_type>("data", "min_pass");
  }

  
  /// \brief Returns true if a maximum number of steps has been specified
  bool EquilibriumMonteSettings::is_max_step() const {
    return _is_setting("data", "max_step");
  }
  
  /// \brief Maximum number of steps
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::max_step() const {
    return _get_setting<size_type>("data", "max_step");
  }
  
  /// \brief Returns true if a minimum number of steps has been specified
  bool EquilibriumMonteSettings::is_min_step() const {
    return _is_setting("data", "min_step");
  
  }
  
  /// \brief Minimum number of steps
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::min_step() const {
    return _get_setting<size_type>("data", "min_step");
  }

  
  /// \brief Returns true if a maximum number of samples has been specified
  bool EquilibriumMonteSettings::is_max_sample() const {
    return _is_setting("data", "max_sample");
  }
  
  /// \brief Maximum number of steps
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::max_sample() const {
    return _get_setting<size_type>("data", "max_sample");
  }
  
  /// \brief Returns true if a minimum number of sample has been specified
  bool EquilibriumMonteSettings::is_min_sample() const {
    return _is_setting("data", "min_sample");
  }
  
  /// \brief Minimum number of steps, default 0
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::min_sample() const {
    return _get_setting<size_type>("data", "min_sample");
  }
  
  
  // --- Data ---------------------
    
  /// \brief Figure out how large data containers should be
  EquilibriumMonteSettings::size_type EquilibriumMonteSettings::max_data_length() const {
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
  

}



