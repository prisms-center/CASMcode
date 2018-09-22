#ifndef CASM_MonteSettings_HH
#define CASM_MonteSettings_HH

#include <string>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/MonteDefinitions.hh"
#include "casm/monte_carlo/MonteSampler.hh"

namespace CASM {

  class MonteCarlo;
  template<typename ObjectType, typename Metric, typename ObjectCompare>
  class HallOfFame;

  class MonteCarloEnumMetric;
  class MonteCarloEnumCheck;

  /*
   * MonteSettings is nothing more than a jsonParser that you can expect
   * to contain a set of values. It's a class meant to hold all the
   * members of the base MonteCarlo class so you can write out all the
   * relevant data and then also use it to reconstruct MonteCarlo.
   *
   * The point of the class is mostly to ensure that when you reconstruct
   * your MonteCarlo you're not giving it some willy nilly jsonParser
   * you came up with. It MUST be something that MonteCarlo itself
   * spit out.
   *
   * Inheritance is private, so you can't use [] as an access operator.
   * Instead use the prescribed accessing routines. This will keep the
   * structure of the settings file consistent and changes will
   * only have to occur through this class.
   *
   * To avoid having millions of access routines, there are jsonParser
   * pointers internally that point to a "current state" of subtrees
   * that may involve the same access routines. For example, instead
   * of having different access routines for initial conditions and
   * final conditions, you just set the current conditions mode, and
   * then the same access routines can be used. From the outside, it
   * looks like you're just requesting a subset of your MonteSettings routine.
   */

  /// \brief Settings for Monte Carlo calculations
  class MonteSettings: protected jsonParser {

  public:

    using jsonParser::print;
    using jsonParser::write;
    using jsonParser::operator==;
    using jsonParser::operator!=;

    typedef Index size_type;

    /// \brief Default constructor
    MonteSettings() {}

    /// \brief Construct MonteSettings by reading a settings JSON file
    MonteSettings(const PrimClex &_primclex, const fs::path &read_path);


    // --- Project root directory ---------------------------

    fs::path root() const;

    const PrimClex &primclex() const;


    // --- Type ---------------------------

    /// \brief Return type of Monte Carlo ensemble
    Monte::ENSEMBLE ensemble() const;

    /// \brief Return type of Monte Carlo method
    Monte::METHOD method() const;

    /// \brief Run in debug mode?
    bool debug() const;

    /// \brief Set debug mode
    void set_debug(bool _debug);

    // --- Initialization ---------------------

    /// \brief Returns true if configname of configuration to use as starting motif has been specified
    bool is_motif_configname() const;

    /// \brief Configname of configuration to use as starting motif
    std::string motif_configname() const;

    /// \brief Returns true if path to ConfigDoF file to use as starting motif has been specified
    bool is_motif_configdof() const;

    /// \brief ConfigDoF to use as starting motif
    ConfigDoF motif_configdof() const;

    /// \brief Path to ConfigDoF file to use as starting motif
    fs::path motif_configdof_path() const;

    /// \brief Supercell matrix defining the simulation cell
    Eigen::Matrix3i simulation_cell_matrix() const;


    // --- Driver ---------------------

    /// \brief Given a settings jsonParser figure out the drive mode. Expects drive_mode/incremental,custom
    virtual const Monte::DRIVE_MODE drive_mode() const;

    /// \brief If dependent runs, start subsequent calculations with the final state
    ///        of the previous calculation. Default true.
    bool dependent_runs() const;


    // --- Sampling -------------------

    /// \brief Given a settings jsonParser figure out the global tolerance (probably for == operator). Expects tolerance/value
    //double tolerance() const;

    /// \brief Requested confidence level. Default 0.95.
    double confidence() const;


    /// \brief Returns true if snapshots are requested
    bool write_trajectory() const;

    /// \brief Returns true if POSCARs of snapshots are requsted. Requires write_trajectory.
    bool write_POSCAR_snapshots() const;

    /// \brief Writes all observations
    bool write_observations() const;

    /// \brief Write csv versions of files? (csv is the default format if no 'output_format' given)
    bool write_csv() const;

    /// \brief Write json versions of files?
    bool write_json() const;

    /// \brief Directory where output should go
    const fs::path output_directory() const;


    // --- Enumerating Configurations ---

    /// \brief Returns true if enumeration is requested. (Default false)
    bool is_enumeration() const;

    /// \brief Returns 'casm query'-like enumeration metric args
    std::string enumeration_metric_args() const;

    /// \brief Returns 'casm query'-like enumeration check args
    std::string enumeration_check_args() const;

    /// \brief Enumeration sample mode (default Monte::ENUM_SAMPLE_MODE::ON_SAMPLE)
    Monte::ENUM_SAMPLE_MODE enumeration_sample_mode() const;

    /// \brief Insert configurations in their canonical form
    bool enumeration_insert_canonical() const;

    /// \brief Only insert configurations that are not already enumerated
    bool enumeration_check_existence() const;

    /// \brief Returns enumeration halloffame max size (default 100)
    Index enumeration_N_halloffame() const;

    /// \brief Returns enumeration halloffame tolerance (default 1e-8)
    double enumeration_tol() const;


  protected:

    /// \brief Returns true if (*this)[level1].contains(level2)
    bool _is_setting(std::string level1, std::string level2) const;

    /// \brief Returns true if (*this)[level1][level2].contains(level3)
    bool _is_setting(std::string level1, std::string level2, std::string level3) const;

    /// \brief Returns (*this)[level1].get<T>();
    template<typename T>
    T _get_setting(std::string level1, std::string msg) const;

    /// \brief Returns (*this)[level1][level2].get<T>();
    template<typename T>
    T _get_setting(std::string level1, std::string level2, std::string msg) const;

    /// \brief Returns (*this)[level1][level2][level3].get<T>();
    template<typename T>
    T _get_setting(std::string level1, std::string level2, std::string level3, std::string msg) const;

    jsonParser _json() const {
      return *this;
    }

  private:

    fs::path m_root;
    fs::path m_output_directory;
    const PrimClex *m_primclex;

  };

  inline bool operator==(const jsonParser &json, const MonteSettings &settings) {
    return settings == json;
  }

  inline bool operator!=(const jsonParser &json, const MonteSettings &settings) {
    return settings != json;
  }


  class EquilibriumMonteSettings : public MonteSettings {


  public:

    /// \brief Default constructor
    EquilibriumMonteSettings() {}

    /// \brief Construct EquilibriumMonteSettings by reading a settings JSON file
    EquilibriumMonteSettings(const PrimClex &_primclex, const fs::path &read_path) :
      MonteSettings(_primclex, read_path) {}


    // --- MCData / Sampling ---------------------

    /// \brief Sample by pass?
    bool sample_by_pass() const;

    /// \brief Sample by step?
    bool sample_by_step() const;

    /// \brief Figure out how often to take samples
    size_type sample_period() const;




    /// \brief Returns true if explicit equilibration passes for the first run have been specified
    bool is_equilibration_passes_first_run() const;

    /// \brief Number of explicit equilibration passes requsted for the first run
    size_type equilibration_passes_first_run() const;

    /// \brief Returns true if explicit equilibration passes for each run have been specified
    bool is_equilibration_passes_each_run() const;

    /// \brief Number of explicit equilibration passes requsted for each run
    size_type equilibration_passes_each_run() const;



    /// \brief Returns true if the number of passes has been specified
    bool is_N_pass() const;

    /// \brief Returns the number of passes requested
    size_type N_pass() const;

    /// \brief Returns true if the number of steps has been specified
    bool is_N_step() const;

    /// \brief Returns the number of steps requested
    size_type N_step() const;

    /// \brief Returns true if the number of samples has been specified
    bool is_N_sample() const;

    /// \brief Returns the number of samples requested
    size_type N_sample() const;



    /// \brief Returns true if a maximum number of passes has been specified
    bool is_max_pass() const;

    /// \brief Maximum number of passes, required if sample by pass
    size_type max_pass() const;

    /// \brief Returns true if a minimum number of passes has been specified
    bool is_min_pass() const;

    /// \brief Minimum number of passes, default 0 if sample by pass
    size_type min_pass() const;


    /// \brief Returns true if a maximum number of steps has been specified
    bool is_max_step() const;

    /// \brief Maximum number of steps, required if sample by step
    size_type max_step() const;

    /// \brief Returns true if a minimum number of steps has been specified
    bool is_min_step() const;

    /// \brief Minimum number of steps, default 0 if sample by step
    size_type min_step() const;


    /// \brief Returns true if a maximum number of samples has been specified
    bool is_max_sample() const;

    /// \brief Maximum number of steps, default std::numeric_limit<size_type>::max()
    size_type max_sample() const;

    /// \brief Returns true if a minimum number of samples has been specified
    bool is_min_sample() const;

    /// \brief Minimum number of steps, default 0
    size_type min_sample() const;


    // --- Data ---------------------

    /// \brief Figure out how large data containers should be
    size_type max_data_length() const;


  };

  /// \brief Returns (*this)[level1][level2].get<T>();
  template<typename T>
  T MonteSettings::_get_setting(std::string level1, std::string msg) const {
    try {
      return (*this)[level1].get<T>();
    }

    catch(std::runtime_error &e) {
      Log &err_log = default_err_log();
      std::stringstream ss;
      ss << "Monte Carlo setting " << "[\"" << level1 << "\"]";

      err_log.error<Log::standard>(ss.str());

      err_log << "Expected " << ss.str() << std::endl;
      err_log << "  Either this was not found, or the type is wrong." << std::endl;
      if(!msg.empty()) {
        err_log << "[\"" << level1 << "\"]: ";
        err_log << msg << std::endl;
      }
      throw;
    }
  }

  /// \brief Returns (*this)[level1][level2].get<T>();
  template<typename T>
  T MonteSettings::_get_setting(std::string level1, std::string level2, std::string msg) const {
    try {
      return (*this)[level1][level2].get<T>();
    }

    catch(std::runtime_error &e) {
      Log &err_log = default_err_log();
      std::stringstream ss;
      ss << "Monte Carlo setting [\"" << level1 << "\"][\"" << level2 << "\"]";

      err_log.error<Log::standard>(ss.str());

      err_log << "ERROR in MonteSettings::" << level2 << std::endl;
      err_log << "Expected " << ss.str() << std::endl;
      err_log << "  Either this was not found, or the type is wrong." << std::endl;
      if(this->contains(level1)) {
        err_log << "Found Settings[\"" << level1 << "\"], but not [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        err_log << "Settings[\"" << level1 << "\"]:\n" << (*this)[level1] << std::endl;
      }
      else {
        err_log << "No Settings[\"" << level1 << "\"] found" << std::endl;
        err_log << "Settings:\n" << static_cast<const jsonParser &>(*this) << std::endl;
      }
      if(!msg.empty()) {
        err_log << "[\"" << level1 << "\"][\"" << level2 << "\"]: ";
        err_log << msg << std::endl;
      }
      throw;
    }
  }

  /// \brief Returns (*this)[level1][level2][level3].get<T>();
  template<typename T>
  T MonteSettings::_get_setting(std::string level1, std::string level2, std::string level3, std::string msg) const {
    try {
      return (*this)[level1][level2][level3].get<T>();
    }

    catch(std::runtime_error &e) {
      Log &err_log = default_err_log();
      std::stringstream ss;
      ss << "Monte Carlo setting [\"" << level1 << "\"][\"" << level2 << "\"][\"" << level3 << "\"]";

      err_log.error<Log::standard>(ss.str());

      err_log << "Expected " << ss.str() << std::endl;
      err_log << "  Either this was not found, or the type is wrong." << std::endl;
      if(this->contains(level1)) {
        if(this->contains(level2)) {
          err_log << "Found Settings[\"" << level1 << "\"][\"" << level2 << "\"], \n"
                  "but not [\"" << level1 << "\"][\"" << level2 << "\"][\"" << level3 << "\"]" << std::endl;
          err_log << "Settings[\"" << level1 << "\"][\"" << level2 << "\"]:\n" << (*this)[level1][level2] << std::endl;
        }
        else {
          err_log << "No Settings[\"" << level1 << "\"][\"" << level2 << "\"] found" << std::endl;
          err_log << "Settings:\n" << (*this)[level1] << std::endl;
        }
      }
      else {
        err_log << "No Settings[\"" << level1 << "\"] found" << std::endl;
        err_log << "Settings:\n" << static_cast<const jsonParser &>(*this) << std::endl;
      }
      if(!msg.empty()) {
        err_log << "[\"" << level1 << "\"][\"" << level2 << "\"][\"" << level3 << "\"]: ";
        err_log << msg << std::endl;
      }
      throw;
    }
  }

}

#endif
