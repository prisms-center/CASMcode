#ifndef CASM_GrandCanonical_HH
#define CASM_GrandCanonical_HH

#include "casm/monte_carlo/MonteDefinitions.hh"
#include "casm/monte_carlo/MonteCarlo.hh"
#include "casm/monte_carlo/SiteExchanger.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalEvent.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalConditions.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalSettings.hh"


namespace CASM {

  ///
  /// Derives from base MonteCarlo class, to be used for simulations at constant
  /// temperature and chemical potential.
  ///
  /// As with all the other derived Monte Carlo classes, member functions must
  /// follow a specific naming convention to be used with templated routines currently
  /// defined in MonteDriver.hh:
  ///      -conditions
  ///      -set_conditions
  ///      -propose
  ///      -check
  ///      -accept
  ///      -reject
  ///      -print_run_start_info
  ///      -write_results
  ///
  class GrandCanonical : public MonteCarlo {

  public:

    typedef GrandCanonicalEvent EventType;
    typedef GrandCanonicalConditions CondType;
    typedef GrandCanonicalSettings SettingsType;


    /// \brief Constructs a GrandCanonical object and prepares it for running based on MonteSettings
    GrandCanonical(PrimClex &primclex, const SettingsType &settings, std::ostream &_sout = std::cout);


    /// \brief Return number of steps per pass. Equals number of sites with variable occupation.
    Index steps_per_pass() const;


    /// \brief Return current conditions
    const CondType &conditions() const;

    /// \brief Set conditions and clear previously collected data
    void set_conditions(const CondType &new_conditions);

    /// \brief Set configdof and clear previously collected data
    void set_configdof(const ConfigDoF& configdof);

    /// \brief Propose a new event, calculate delta properties, and return reference to it
    const EventType &propose();

    /// \brief Based on a random number, decide if the change in energy from the proposed event is low enough to be accepted.
    bool check(const EventType &event);

    /// \brief Accept proposed event. Change configuration accordingly and update energies etc.
    void accept(const EventType &event);

    /// \brief Nothing needs to be done to reject a GrandCanonicalEvent
    void reject(const EventType &event);


    /// \brief Print info when a run begins
    void print_run_start_info() const;

    /// \brief Write results to files
    void write_results(Index cond_index) const;


    /// \brief Calculate the single spin flip low temperature expansion of the grand canonical potential
    double lte_grand_canonical_free_energy(std::ostream &sout) const;

    /// \brief Formation energy, normalized per primitive cell
    const double &formation_energy() const {
      return *m_formation_energy;
    }

    /// \brief Potential energy, normalized per primitive cell
    const double &potential_energy() const {
      return *m_potential_energy;
    }

    /// \brief Correlations, normalized per primitive cell
    const Eigen::VectorXd &corr() const {
      return *m_corr;
    }

    /// \brief Number of atoms of each type, normalized per primitive cell
    const Eigen::VectorXd &comp_n() const {
      return *m_comp_n;
    }


  private:

    /// \brief Formation energy, normalized per primitive cell
    double &_formation_energy() {
      return *m_formation_energy;
    }

    /// \brief Potential energy, normalized per primitive cell
    double &_potential_energy() {
      return *m_potential_energy;
    }

    /// \brief Correlations, normalized per primitive cell
    Eigen::VectorXd &_corr() {
      return *m_corr;
    }

    /// \brief Number of atoms of each type, normalized per primitive cell
    Eigen::VectorXd &_comp_n() {
      return *m_comp_n;
    }

    /// \brief Calculate delta properties for an event and update the event with those properties
    void _update_deltas(GrandCanonicalEvent &event,
                        Index mutating_site,
                        int sublat,
                        int current_occupant,
                        int new_occupant) const;

    /// \brief Calculate properties given current conditions
    void _update_properties();

    /// \brief Select initial configdof
    static ConfigDoF _initial_configdof(
      PrimClex &primclex,
      Supercell &scel,
      const GrandCanonicalSettings &settings,
      std::ostream &_sout);


    ///Keeps track of what sites can change to what
    const SiteExchanger m_site_swaps;

    /// Conditions (T, mu). Initially determined by m_settings, but can be changed halfway through the run
    GrandCanonicalConditions m_condition;

    ///Clexulator that is used to calculate the energy from the current simulation state. ConfigIterator breaks const-ness!!!
    mutable Clexulator m_clexulator;

    /// This is the typical eci set that you use to get the energy
    ECIContainer m_formation_energy_eci;

    /// If true, calculate all correlations; if false, calculate correlations with non-zero eci
    bool m_all_correlations;

    /// Event to propose, check, accept/reject:
    EventType m_event;

    /// \brief Change in comp_n due of atom being removed. Equal to -1.0/supercell().volume()
    double m_minus_one_comp_n;

    /// \brief Change in comp_n due of atom being added. Equal to 1.0/supercell().volume()
    double m_plus_one_comp_n;



    // ---- Pointers to properties for faster access

    /// \brief Formation energy, normalized per primitive cell
    double *m_formation_energy;

    /// \brief Potential energy, normalized per primitive cell
    double *m_potential_energy;

    /// \brief Correlations, normalized per primitive cell
    Eigen::VectorXd *m_corr;

    /// \brief Number of atoms of each type, normalized per primitive cell
    Eigen::VectorXd *m_comp_n;

  };
  
}

#endif




