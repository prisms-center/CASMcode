#ifndef CASM_MonteCarlo_HH
#define CASM_MonteCarlo_HH

#include <vector>
#include "casm/misc/cloneable_ptr.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"
#include "casm/casm_io/Log.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/monte_carlo/MonteSettings.hh"
#include "casm/monte_carlo/MonteSampler.hh"
#include "casm/monte_carlo/MonteCounter.hh"

namespace CASM {

  struct SamplerNameCompare {

    SamplerNameCompare() {};

    bool operator()(const std::string &A, const std::string &B) const;

  };


  /// \brief Interface base class for all types of Monte Carlo simulations (not meant to be used polymorphically)
  ///
  /// MonteCarlo serves as a base for all types of Monte Carlo simulations (e.g. GrandCanonical).
  /// The base class holds all the data structures to be sampled for every type of derived Monte Carlo class.
  /// The derived class is responsible for adding data and using it to perform a calculation.
  ///
  class MonteCarlo {


  public:

    typedef MonteSampler::size_type size_type;

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
    typedef std::map<std::string, notstd::cloneable_ptr<MonteSampler>, SamplerNameCompare> SamplerMap;

    /// \brief a vector of std::pair(pass, step) indicating when samples were taken
    typedef std::vector<std::pair<MonteCounter::size_type, MonteCounter::size_type> > SampleTimes;


    // ---- Accessors -----------------------------

    /// \brief const Access settings used for construction
    const MonteSettings &settings() const {
      return m_settings;
    }

    /// \brief const Access the PrimClex that *this is based on
    const PrimClex &primclex() const {
      return m_primclex;
    }

    /// \brief const Access the Supercell that *this is based on
    const Supercell &supercell() const {
      return m_scel;
    }

    /// \brief Set a pointer to the SuperNeighborList once it is ready
    void set_nlist() {
      m_nlist = &(supercell().nlist());
    }

    /// \brief const Access the SuperNeighborList via pointer stored by 'set_nlist'
    const SuperNeighborList &nlist() const {
      return *m_nlist;
    }

    /// \brief const Access current microstate
    const Configuration &config() const {
      return m_config;
    }

    /// \brief const Access current microstate
    const ConfigDoF &configdof() const {
      return m_configdof;
    }


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

    /// \brief Samples all requested property data, and stores pass and step number sample was taken at
    void sample_data(const MonteCounter &counter);

    /// \brief Clear all data from all samplers
    void clear_samples();

    /// \brief Return true if convergence is requested
    bool must_converge() const {
      return m_must_converge;
    }

    /// \brief Returns pair(true, equil_samples) if required equilibration has occured for all samplers that must converge
    std::pair<bool, MonteSampler::size_type> is_equilibrated() const;

    /// \brief Returns true if a convergence check is due
    bool check_convergence_time() const;

    /// \brief Check to see if all the properties required to converge have converged
    bool is_converged() const;

    /// \brief const Access sampler map
    const SamplerMap &samplers() const {
      return m_sampler;
    }

    /// \brief const Access a vector of std::pair<pass, step> indicating when samples were taken
    const SampleTimes &sample_times() const {
      return m_sample_time;
    }

    /// \brief const Access snapshots of the Monte Carlo calculation
    ///
    /// If requested, snapshots are taken at the same time as samples. So examine sample_times for
    /// pass and step information.
    const std::vector<ConfigDoF> &trajectory() const {
      return m_trajectory;
    }

    /// \brief return true if running in debug mode
    bool debug() const {
      return m_debug;
    }


  protected:

    /// \brief Construct with a starting ConfigDoF as specified the given MonteSettings and prepare data samplers
    template<typename MonteTypeSettings>
    MonteCarlo(PrimClex &primclex, const MonteTypeSettings &settings, Log &_log);

    /// \brief Access the PrimClex that *this is based on
    PrimClex &_primclex() const {
      return m_primclex;
    }

    /// \brief Access the Supercell that *this is based on
    Supercell &_supercell() const {
      return m_scel;
    }

    /// \brief Access current microstate
    ///
    /// - This can be used by a const member if it undoes any changes
    ///   to the ConfigDoF before returning
    Configuration &_config() const {
      return m_config;
    }

    /// \brief Access current microstate
    ///
    /// - This can be used by a const member if it undoes any changes
    ///   to the ConfigDoF before returning
    ConfigDoF &_configdof() const {
      return m_configdof;
    }

    Log &_log() const {
      return m_log;
    }

    MTRand &_mtrand() {
      return m_twister;
    }

    /// \brief Access scalar properties map
    ScalarPropertyMap &_scalar_properties() {
      return m_scalar_property;
    }

    /// \brief Access a particular scalar property
    double &_scalar_property(std::string property_name) {
      return m_scalar_property.find(property_name)->second;
    }

    /// \brief const Access vector properties map
    VectorPropertyMap &_vector_properties() {
      return m_vector_property;
    }

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
    PrimClex &m_primclex;

    /// \brief Supercell for the calculation.
    mutable Supercell m_scel;

    /// \brief Pointer to SuperNeighborList
    const SuperNeighborList *m_nlist;

    /// \brief Stores all degrees of freedom of the current microstate
    ///
    /// 'mutable' is used for case where the DoF are modified to calculate
    /// event property values and then reverted within a const function
    mutable Configuration m_config;

    /// Reference to m_config.configdof(), to avoid invalidating id every time used
    ConfigDoF &m_configdof;

    /// \brief Random number generator
    MTRand m_twister;


    /// \brief Save trajectory?
    bool m_write_trajectory = false;

    /// \brief Snapshots of the Monte Carlo simulation, taken by sample_data() if m_write_trajectory is true
    std::vector<ConfigDoF> m_trajectory;

    /// \brief True if any Sampler must converge
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

    /// \brief a vector of std::pair(pass, step) indicating when samples were taken
    SampleTimes m_sample_time;


    // Members to remember results of is_equilibrated and is_converged

    mutable std::pair<bool, MonteSampler::size_type> m_is_equil;
    mutable bool m_is_equil_uptodate = false;
    mutable bool m_is_converged;
    mutable bool m_is_converged_uptodate = false;

    // which sample to check convergence after
    mutable MonteSampler::size_type m_next_convergence_check = 100;
    MonteSampler::size_type m_convergence_check_period = 100;

    // in debug mode, allow printing or checking extra things
    bool m_debug;

  };

  /// \brief Construct with a starting ConfigDoF as specified the given MonteSettings and prepare data samplers
  template<typename MonteTypeSettings>
  MonteCarlo::MonteCarlo(PrimClex &primclex, const MonteTypeSettings &settings, Log &_log) :
    m_settings(settings),
    m_primclex(primclex),
    m_scel(&primclex, settings.simulation_cell_matrix()),
    m_config(m_scel),
    m_configdof(m_config.configdof()),
    m_write_trajectory(settings.write_trajectory()),
    m_log(_log),
    m_debug(m_settings.debug()) {

    settings.samplers(primclex, std::inserter(m_sampler, m_sampler.begin()));

    m_must_converge = false;
    for(auto it = m_sampler.cbegin(); it != m_sampler.cend(); ++it) {
      if(it->second->must_converge()) {
        m_must_converge = true;
        break;
      }
    }
  }

}
#endif
