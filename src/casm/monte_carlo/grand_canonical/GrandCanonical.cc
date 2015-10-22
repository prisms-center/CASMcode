
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

namespace CASM {

  
  GrandCanonical::GrandCanonical(PrimClex &primclex, const GrandCanonicalSettings &settings):
    MonteCarlo(primclex, settings), 
    m_site_swaps(supercell()),
    m_condition(settings.initial_conditions()),
    m_clexulator(primclex.global_clexulator()),
    m_formation_energy_eci(primclex.dir().eci_out(settings.clex(), 
                                                  settings.bset(), 
                                                  settings.calctype(), 
                                                  settings.ref(), 
                                                  settings.eci())),
    m_event(primclex.composition_axes().components().size(), m_clexulator.corr_size()),
    m_minus_one_comp_n(-1.0/supercell().volume()),
    m_plus_one_comp_n(1.0/supercell().volume()) {
    
    
    
    /// Prepare for calculating correlations. Maybe this should get put into Clexulator.
    const DirectoryStructure& dir = primclex.dir();
    std::cout << "Read the global orbitree \n";
    if(fs::exists(dir.clexulator_src(primclex.settings().name(), settings.bset()))) {
      primclex.read_global_orbitree(dir.clust(settings.bset()));
    }
    std::cout << "  DONE." << std::endl << std::endl;
    
    primclex.generate_full_nlist();
    supercell().generate_neighbor_list();
    
    // Make sure the simulation is big enough to accommodate the clusters 
    // you're using so that the delta formation energy is calculated accurately
    if(supercell().neighbor_image_overlaps()) {
      throw std::runtime_error(
        std::string("ERROR in 'GrandCanonical(PrimClex &primclex, const MonteSettings &settings)'\n") +
                    "  The simulation cell is too small to fit all the clusters without periodic overlap.\n" +
                    "  This would result in incorrect calculations of dformation_energy.\n" +
                    "  You need a smaller orbitree or a larger simulation cell.");
    }
    
    
    // initialize properties and store pointers to the data strucures
    m_vector_property["corr"] = correlations_vec(m_configdof, supercell(), m_clexulator);
    m_corr = &m_vector_property["corr"];
    
    m_vector_property["comp_n"] = CASM::comp_n(m_configdof, supercell());
    m_comp_n = &m_vector_property["comp_n"];
    
    m_scalar_property["formation_energy"] = m_formation_energy_eci * corr().data();
    m_formation_energy = &m_scalar_property["formation_energy"];
    
    m_scalar_property["potential_energy"] = formation_energy() - comp_n().dot(m_condition.mu());
    m_potential_energy = &m_scalar_property["potential_energy"]; 
    
  }
  
  /// \brief Return number of steps per pass. Equals number of sites with variable occupation.
  Index GrandCanonical::steps_per_pass() const {
    return m_site_swaps.variable_sites().size();
  }
    
  
  /// \brief Return current conditions
  const GrandCanonical::CondType& GrandCanonical::conditions() const {
    return m_condition;
  }
  
  
  /// \brief Set conditions and clear previously collected data
  void GrandCanonical::set_conditions(const GrandCanonicalConditions &new_conditions) {
    m_condition = new_conditions;
    
    clear_samples();
    
    return;
  }
  

  /// \brief Propose a new event, calculate delta properties, and return reference to it
  ///
  /// Randomly picks a site that's allowed more than one occupant, and randomly picks what occupant it
  /// changes to. Then calculates delta properties associated with that change.
  ///
  const GrandCanonical::EventType& GrandCanonical::propose() {
    
    // Randomly pick a site that's allowed more than one occupant
    Index random_variable_site = m_twister.randInt(m_site_swaps.variable_sites().size() - 1);
    
    // Determine what that site's linear index is and what the sublattice index is
    Index mutating_site = m_site_swaps.variable_sites()[random_variable_site];
    Index sublat = m_site_swaps.sublat()[random_variable_site];

    // Determine the current occupant of the mutating site
    int current_occupant = configdof().occ(mutating_site);
    
    // Randomly pick a new occupant for the mutating site
    const std::vector<int> &possible_mutation = m_site_swaps.possible_swap()[sublat][current_occupant];
    int new_occupant = possible_mutation[m_twister.randInt(possible_mutation.size() - 1)];
    
    // Update delta properties in m_event
    _update_deltas(m_event, mutating_site, sublat, current_occupant, new_occupant);
    
    return m_event;
  }
  
  /// \brief Based on a random number, decide if the change in energy from the proposed event is low enough to be accepted.
  bool GrandCanonical::check(const GrandCanonicalEvent &event) {
    
    double rand = m_twister.rand53();
    
    if(event.dpotential_energy() < 0.0 || rand < exp(-event.dpotential_energy() * m_condition.beta() * supercell().volume() )) {
      return true;
    }

    return false;
  }
  
  /// \brief Accept proposed event. Change configuration accordingly and update energies etc.
  ///
  /// Once you've found an event that passes the random number test, you want to accept it. This routine will
  /// apply all the changes in the event to *this (change occupation of one atom) and update the formation energy,
  /// generalized enthalpy, number of species and correlations values.
  ///
  void GrandCanonical::accept(const EventType& event) {

    // First apply changes to configuration (just a single occupant change)
    m_configdof.occ(event.occupational_change().site_index()) = event.occupational_change().to_value();
    
    // Next update all properties that changed from the event
    formation_energy() += event.dformation_energy();
    potential_energy() += event.dpotential_energy();
    corr() += event.dcorr();
    comp_n() += event.dcomp_n();
    
    return;
  }

  /// \brief Nothing needs to be done to reject a GrandCanonicalEvent
  void GrandCanonical::reject(const EventType& event) {
    return;
  }

  /// \brief Calculate the low temperature expansion of the grand canonical free energy
  ///
  /// Returns low temperature expansion estimate of the grand canonical free energy.
  /// Works with the current ConfigDoF as groundstate.
  /// 
  /// Quick derivation:
  /// Z: partition function
  /// boltz(x): exp(-x/kBT)
  /// \Omega: E-SUM(mu*comp_n)
  /// 
  /// The partition function is
  /// Z=SUM(boltz(\Omega_s))    summing over all microstates s
  /// 
  /// \Omega_s can be split into groundstate \Omega_0 and a delta energy D\Omega
  /// \Omega_s=\Omega_0+D\Omega_s
  /// Z=boltz(\Omega_0)*SUM(boltz(D\Omega_s))  summing over all microstates
  /// 
  /// For low temperatures we can approximate Z by truncating the sum after microstates that
  /// only involve point defects and no defects
  /// Z=boltz(\Omega_0)*SUM(boltz(D\Omega_s))  summing over all states with only point defects or no defects
  /// 
  /// The free energy is
  /// Phi=-kB*T*ln(Z)
  /// Phi=-kB*T*(-\Omega_0/kBT+ln(SUM(boltz(D\Omega_s))    Sum is over point defects and no defects (in which case D\Omega_s == 0)
  /// Phi=\Omega_0-kB*T(ln(SUM(boltz(D\Omega_s))            
  /// 
  /// Now for numerical reasons, use ln(SUM(exp(x_s))) = m + ln(SUM(exp(x_s - m))), where m == max(x_s)
  ///
  /// So Phi = \Omega_0-kB*T*(m + ln(SUM(exp(x_s - m))), where x_s = -D\Omega_s / kB / T, and m == max(x_s)
  ///
  double GrandCanonical::lte_grand_canonical_free_energy() const {

    //Now we have \Omega_0 (grand_canonical_energy)
    
    const SiteExchanger& site_exch = m_site_swaps;
    const ConfigDoF& config_dof = m_configdof;
    GrandCanonicalEvent event = m_event;
    
    //This will hold all the possible x_s (x_s = -D\Omega_s / kB / T)  for all point defects and the no defect case
    std::vector<double> x;
    
    // no defect case
    x.push_back(0.0);
    
    //Loop over sites that can change occupants
    for(Index exch_ind = 0; exch_ind < site_exch.variable_sites().size(); exch_ind++) {

      //Transform exchanger index to ConfigDoF index
      Index mutating_site = site_exch.variable_sites()[exch_ind];
      int sublat = site_exch.sublat()[exch_ind];
      int current_occupant = config_dof.occ(mutating_site);

      //Loop over possible occupants for site that can change
      for(Index new_occupant_ind = 0; new_occupant_ind < site_exch.possible_swap()[sublat][current_occupant].size(); new_occupant_ind++) {
        
        int new_occupant = site_exch.sublat_to_mol()[sublat][new_occupant_ind];
        
        _update_deltas(event, mutating_site, sublat, current_occupant, new_occupant);

        //save the result
        x.push_back(-event.dpotential_energy() * m_condition.beta() * supercell().volume());

      }
    }
    
    // So Phi = \Omega_0-kB*T*(m + ln(SUM(exp(x_s - m))), where x_s = -D\Omega_s / kB / T, and m = max(x_s)
    
    // find max of x_s
    double m = *std::max_element(x.cbegin(), x.cend());
    
    // calc ln(SUM(exp(x_s - m))
    double sum_exp = 0.0;
    for(auto it = x.cbegin(); it != x.cend(); ++it) {
      sum_exp += exp(*it-m);
    }
    
    return potential_energy() - (m + std::log(sum_exp)) / m_condition.beta();

  }

  /// \brief Print info when a run begins
  void GrandCanonical::print_run_start_info(std::ostream &sout) const {
    sout << "Begin run.  T = " << m_condition.temperature() << "  ";
    
    Eigen::VectorXd param_mu = primclex().composition_axes().param_mu(m_condition.mu());
    for(int i=0; i<param_mu.size(); i++) {
      sout << "mu_" << ((char) (i + (int) 'a')) << " = " << param_mu(i) << "  ";
    }
    sout << std::endl;
  }
  
  /// \brief Write results to files
  void GrandCanonical::write_results(Index cond_index) const {
    CASM::write_results(settings(), *this);
    write_conditions_json(settings(), *this, cond_index);
    write_observations(settings(), *this, cond_index);
    write_trajectory(settings(), *this, cond_index);
    write_pos_trajectory(settings(), *this, cond_index);
  }
  
  /// \brief Update delta properties in 'event'
  void GrandCanonical::_update_deltas(GrandCanonicalEvent &event,
                                      Index mutating_site,
                                      int sublat,
                                      int current_occupant,
                                      int new_occupant) const {
    
    // ---- set OccMod --------------
    
    event.occupational_change().set(mutating_site, sublat, new_occupant);
    
    // ---- set dspecies --------------
    
    for(int i=0; i<event.dcomp_n().size(); ++i) {
      event.set_dcomp_n(i, 0.0);
    }
    event.set_dcomp_n(m_site_swaps.sublat_to_mol()[sublat][current_occupant], m_minus_one_comp_n);
    event.set_dcomp_n(m_site_swaps.sublat_to_mol()[sublat][new_occupant], m_plus_one_comp_n);
    
    
    // ---- set dcorr --------------
    
    // Point the Clexulator to the right neighborhood 
    m_clexulator.set_nlist(supercell().get_nlist(mutating_site).begin());

    // Calculate the change in correlations due to this event
    m_clexulator.calc_restricted_delta_point_corr(sublat,
                                                  current_occupant,
                                                  new_occupant,
                                                  event.dcorr().data(),
                                                  m_formation_energy_eci.eci_index_list().begin(),
                                                  m_formation_energy_eci.eci_index_list().end());
    
    // ---- set dformation_energy --------------
    
    event.set_dformation_energy((m_formation_energy_eci * event.dcorr().data()) / supercell().volume());
    
    
    // ---- set dpotential_energy --------------
    
    
    event.set_dpotential_energy(event.dformation_energy() - event.dcomp_n().dot(m_condition.mu()));
    
  }
  

  
}

