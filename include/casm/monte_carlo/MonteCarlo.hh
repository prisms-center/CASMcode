#ifndef CASM_MonteCarlo_HH
#define CASM_MonteCarlo_HH

#include <string>
#include <vector>

#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/global/definitions.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/monte_carlo/MonteDefinitions.hh"

namespace CASM {

class Log;
class PrimClex;
class SuperNeighborList;

}  // namespace CASM

namespace CASM {
namespace Monte {

class MonteSampler;
class MonteSettings;
class MonteCounter;

struct SamplerNameCompare {
  SamplerNameCompare(){};

  bool operator()(const std::string &A, const std::string &B) const;
};

/// \brief Interface base class for all types of Monte Carlo simulations (not
/// meant to be used polymorphically)
///
/// MonteCarlo serves as a base for all types of Monte Carlo simulations (e.g.
/// GrandCanonical). The base class holds all the data structures to be sampled
/// for every type of derived Monte Carlo class. The derived class is
/// responsible for adding data and using it to perform a calculation.
///
class MonteCarlo {
 public:
  typedef Monte::size_type size_type;

  /// \brief a map of keyname to property value
  ///
  /// - example: m_scalar_property["formation_energy"]
  typedef std::map<std::string, double> ScalarPropertyMap;

  /// \brief a map of keyname to property value Eigen::VectorXd
  ///
  /// - example: m_vector_property["corr"]
  typedef std::map<std::string, Eigen::VectorXd> VectorPropertyMap;

  /// \brief a map of keyname to MonteSampler
  ///
  /// - scalar example: m_sampler["formation_energy"]
  /// - vector example: m_sampler["corr(12)"]
  ///
  typedef std::map<std::string, notstd::cloneable_ptr<MonteSampler>,
                   SamplerNameCompare>
      SamplerMap;

  /// \brief a vector of std::pair(pass, step) indicating when samples were
  /// taken
  typedef std::vector<std::pair<size_type, size_type> > SampleTimes;

  // ---- Accessors -----------------------------

  /// \brief const Access settings used for construction
  const MonteSettings &settings() const { return m_settings; }

  /// \brief const Access the PrimClex that *this is based on
  const PrimClex &primclex() const { return m_primclex; }

  /// \brief const Access the Supercell that *this is based on
  const Supercell &supercell() const { return m_scel; }

  /// \brief const Access current microstate
  const Configuration &config() const { return m_config; }

  /// \brief const Access current microstate
  const ConfigDoF &configdof() const { return m_configdof; }

  // ---- Accessors -----------------------------

  /// \brief Set current microstate and clear samplers
  void reset(const ConfigDoF &dof) {
    _configdof() = dof;
    clear_samples();
  }

  // ---- Properties ----------------

  /// \brief const Access scalar properties map
  const ScalarPropertyMap &scalar_properties() const {
    return m_scalar_property;
  }

  /// \brief const Access a particular scalar property
  const double &scalar_property(std::string property_name) const {
    return m_scalar_property.find(property_name)->second;
  }

  /// \brief const Access vector properties map
  const VectorPropertyMap &vector_properties() const {
    return m_vector_property;
  }

  /// \brief const Access a particular vector property
  const Eigen::VectorXd &vector_property(std::string property_name) const {
    return m_vector_property.find(property_name)->second;
  }

  // ---- Data sampling -------------

  /// \brief Samples all requested property data, and stores pass and step
  /// number sample was taken at
  void sample_data(const MonteCounter &counter);

  /// \brief Clear all data from all samplers
  void clear_samples();

  /// \brief Return true if convergence is requested
  bool must_converge() const { return m_must_converge; }

  /// \brief Returns pair(true, equil_samples) if required equilibration has
  /// occured for all samplers that must converge
  std::pair<bool, size_type> is_equilibrated() const;

  /// \brief Returns true if a convergence check is due
  bool check_convergence_time() const;

  /// \brief Check to see if all the properties required to converge have
  /// converged
  bool is_converged() const;

  /// \brief const Access sampler map
  const SamplerMap &samplers() const { return m_sampler; }

  /// \brief const Access a vector of std::pair<pass, step> indicating when
  /// samples were taken
  const SampleTimes &sample_times() const { return m_sample_time; }

  /// \brief const Access snapshots of the Monte Carlo calculation
  ///
  /// If requested, snapshots are taken at the same time as samples. So examine
  /// sample_times for pass and step information.
  const std::vector<ConfigDoF> &trajectory() const { return m_trajectory; }

  /// \brief return true if running in debug mode
  bool debug() const { return m_debug; }

 protected:
  /// \brief Construct with a starting ConfigDoF as specified the given
  /// MonteSettings and prepare data samplers
  template <typename MonteTypeSettings>
  MonteCarlo(const PrimClex &primclex, const MonteTypeSettings &settings,
             Log &_log);

  /// \brief Access the PrimClex that *this is based on
  const PrimClex &_primclex() const { return m_primclex; }

  /// \brief Access the Supercell that *this is based on
  const Supercell &_supercell() const { return m_scel; }

  /// \brief Access current microstate
  ///
  /// - This can be used by a const member if it undoes any changes
  ///   to the ConfigDoF before returning
  Configuration &_config() const { return m_config; }

  /// \brief Access current microstate
  ///
  /// - This can be used by a const member if it undoes any changes
  ///   to the ConfigDoF before returning
  ConfigDoF &_configdof() const { return m_configdof; }

  Log &_log() const { return m_log; }

  MTRand &_mtrand() { return m_twister; }

  /// \brief Access scalar properties map
  ScalarPropertyMap &_scalar_properties() { return m_scalar_property; }

  /// \brief Access a particular scalar property
  double &_scalar_property(std::string property_name) {
    return m_scalar_property.find(property_name)->second;
  }

  /// \brief const Access vector properties map
  VectorPropertyMap &_vector_properties() { return m_vector_property; }

  /// \brief const Access a particular vector property
  Eigen::VectorXd &_vector_property(std::string property_name) {
    return m_vector_property.find(property_name)->second;
  }

 private:
  /// \brief a map of keyname to property value
  ///
  /// - example: m_scalar_property["formation_energy"]
  ScalarPropertyMap m_scalar_property;

  /// \brief a map of keyname to property value Eigen::VectorXd
  ///
  /// - example: m_vector_property["corr"]
  VectorPropertyMap m_vector_property;

  /// \brief Contains all input settings
  const MonteSettings &m_settings;

  /// \brief PrimClex for this system
  const PrimClex &m_primclex;

  /// \brief Supercell for the calculation.
  Supercell m_scel;

  /// \brief Stores all degrees of freedom of the current microstate
  ///
  /// 'mutable' is used for case where the DoF are modified to calculate
  /// event property values and then reverted within a const function
  mutable Configuration m_config;

  /// Reference to m_config.configdof(), to avoid invalidating id every time
  /// used
  ConfigDoF &m_configdof;

  /// \brief Random number generator
  MTRand m_twister;

  /// \brief Save trajectory?
  bool m_write_trajectory = false;

  /// \brief Snapshots of the Monte Carlo simulation, taken by sample_data() if
  /// m_write_trajectory is true
  std::vector<ConfigDoF> m_trajectory;

  /// \brief True if any MonteSampler must converge
  bool m_must_converge;

  /// \brief Target for messages
  Log &m_log;

  /// \brief Set the next time convergence is due to be checked
  void _set_check_convergence_time() const;

  /// \brief a map of pair<keyname, index> to MonteSampler
  ///
  /// - scalar example: m_sampler[std::make_pair("formation_energy", 0)]
  /// - vector example: m_sampler[std::make_pair("corr", 12)]
  ///
  SamplerMap m_sampler;

  /// \brief a vector of std::pair(pass, step) indicating when samples were
  /// taken
  SampleTimes m_sample_time;

  // Members to remember results of is_equilibrated and is_converged

  mutable std::pair<bool, size_type> m_is_equil;
  mutable bool m_is_equil_uptodate = false;
  mutable bool m_is_converged;
  mutable bool m_is_converged_uptodate = false;

  // which sample to check convergence after
  mutable size_type m_next_convergence_check = 100;
  size_type m_convergence_check_period = 100;

  // in debug mode, allow printing or checking extra things
  bool m_debug;
};

}  // namespace Monte
}  // namespace CASM

#endif
