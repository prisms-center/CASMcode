
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/Norm.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

namespace CASM {

  
  GrandCanonical::GrandCanonical(PrimClex &primclex, const GrandCanonicalSettings &settings, std::ostream& _sout):
    MonteCarlo(primclex, 
               _select_motif(settings.motif_configname(), primclex, settings.initial_conditions(), _sout), 
               settings, 
               _sout), 
    m_site_swaps(supercell()),
    m_condition(settings.initial_conditions()),
    m_clexulator(primclex.global_clexulator()),
    m_formation_energy_eci(
      read_eci(
        primclex.dir().eci(
          settings.clex(), 
          settings.calctype(), 
          settings.ref(), 
          settings.bset(), 
          settings.eci()
        )
      )
    ),
    m_all_correlations(settings.all_correlations()),
    m_event(primclex.composition_axes().components().size(), m_clexulator.corr_size()),
    m_minus_one_comp_n(-1.0/supercell().volume()),
    m_plus_one_comp_n(1.0/supercell().volume()) {
    
    // set the SuperNeighborList... 
    set_nlist();
    
    // Make sure the simulation is big enough to accommodate the clusters 
    // you're using so that the delta formation energy is calculated accurately
    if(nlist().overlaps()) {
      throw std::runtime_error(
        std::string("ERROR in 'GrandCanonical(PrimClex &primclex, const MonteSettings &settings)'\n") +
                    "  The simulation cell is too small to fit all the clusters without periodic overlap.\n" +
                    "  This would result in incorrect calculations of dformation_energy.\n" +
                    "  You need a smaller orbitree or a larger simulation cell.");
    }
    
    _update_properties(); 
  
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
    _update_properties();
    
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
    
    if(debug()) {
      const auto& site_occ = primclex().get_prim().basis[sublat].site_occupant();
      sout << "\n-- Proposed event -- \n\n"
           
           << "  Mutating site (linear index): " << mutating_site << "\n"
           << "  Mutating site (b, i, j, k): " << supercell().uccoord(mutating_site) << "\n"
           << "  Current occupant: " << current_occupant << " (" << site_occ[current_occupant].name << ")\n"
           << "  Proposed occupant: " << new_occupant << " (" << site_occ[new_occupant].name << ")\n\n"
           
           << "  beta: " << m_condition.beta() << "\n"
           << "  T: " << m_condition.temperature() << std::endl;
           
      sout << std::setw(12) << "i" << std::setw(16) << "ECI" << std::setw(16) << "dcorr" << std::endl;
      
      if(m_all_correlations) {
        for(int i=0; i<m_event.dcorr().size(); ++i) {
          
          double eci = 0.0;
          Index index = find_index(m_formation_energy_eci.index(), i);
          if(index != m_formation_energy_eci.index().size()) {
            eci = m_formation_energy_eci.value()[index];
          }
          
          sout << std::setw(12) << i 
               << std::setw(16) << std::setprecision(8) << eci 
               << std::setw(16) << std::setprecision(8) << m_event.dcorr()[i] << std::endl;
        
        }
      }
      else {
        for(int i=0; i<m_formation_energy_eci.value().size(); ++i) {
          sout << std::setw(12) << m_formation_energy_eci.index()[i] 
               << std::setw(16) << std::setprecision(8) << m_formation_energy_eci.value()[i] 
               << std::setw(16) << std::setprecision(8) << m_event.dcorr()[m_formation_energy_eci.index()[i]] << std::endl;
        
        }
      }
      
      auto origin = primclex().composition_axes().origin();
      auto chem_pot = m_condition.chem_pot();
      auto param_chem_pot = m_condition.param_chem_pot();
      auto dcomp_n = m_event.dcomp_n();
      auto dcomp_x = primclex().composition_axes().dparam_composition(dcomp_n);
      auto M = primclex().composition_axes().dparam_dmol();
      auto V = supercell().volume();
           
      sout << "  components: " << jsonParser(primclex().composition_axes().components()) << "\n"
           << "  N*dcomp_n: " << V*dcomp_n.transpose() << "\n"
           << "  chem_pot: " << chem_pot.transpose() << "\n"
           << "  N*dcomp_x: " << V*dcomp_x.transpose() << "\n"
           << "  param_chem_pot: " << param_chem_pot.transpose() << "\n"
           << "   N*param_chem_pot*dcomp_x: " << V*param_chem_pot.dot(dcomp_x) << "\n"
           << "  N*dformation_energy: " << V*m_event.dformation_energy() << "\n"
           << "    N*(dformation_energy - chem_pot*dcomp_n): " << V*(m_event.dformation_energy() - chem_pot.dot(dcomp_n)) << "\n"
           << "    N*(dformation_energy - parm_chem_pot*dcomp_x: " << V*(m_event.dformation_energy() - param_chem_pot.dot(dcomp_x)) << "\n"
           << "  N*dpotential_energy: " << V*m_event.dpotential_energy() << "\n" << std::endl;
           
      
    }
    
    return m_event;
  }
  
  /// \brief Based on a random number, decide if the change in energy from the proposed event is low enough to be accepted.
  bool GrandCanonical::check(const GrandCanonicalEvent &event) {
    
    if(event.dpotential_energy() < 0.0) {
      
      if(debug()) {
        sout << "Probability to accept: 1.0\n" << std::endl;
      }
      return true;
    }
    
    double rand = m_twister.rand53();
    double prob = exp(-event.dpotential_energy() * m_condition.beta() * supercell().volume());
    
    if(debug()) {
      sout << "Probability to accept: " << prob << "\n"
           << "Random number: " << rand << "\n" << std::endl;
    }
    
    return rand < prob;
  }
  
  /// \brief Accept proposed event. Change configuration accordingly and update energies etc.
  ///
  /// Once you've found an event that passes the random number test, you want to accept it. This routine will
  /// apply all the changes in the event to *this (change occupation of one atom) and update the formation energy,
  /// generalized enthalpy, number of species and correlations values.
  ///
  void GrandCanonical::accept(const EventType& event) {
    
    if(debug()) {
      sout << "** Accept event **" << std::endl;
    }
    
    // First apply changes to configuration (just a single occupant change)
    m_configdof.occ(event.occupational_change().site_index()) = event.occupational_change().to_value();
    
    // Next update all properties that changed from the event
    _formation_energy() += event.dformation_energy();
    _potential_energy() += event.dpotential_energy();
    _corr() += event.dcorr();
    _comp_n() += event.dcomp_n();
    
    return;
  }

  /// \brief Nothing needs to be done to reject a GrandCanonicalEvent
  void GrandCanonical::reject(const EventType& event) {
    if(debug()) {
      sout << "** Reject event **" << std::endl;
    }
    return;
  }

  /// \brief Calculate the single spin flip low temperature expansion of the grand canonical potential
  ///
  /// \param sout Stream to print spin flip details 
  ///
  /// Returns low temperature expansion estimate of the grand canonical free energy.
  /// Works with the current ConfigDoF as groundstate.
  /// 
  /// Quick derivation:
  /// Z: partition function
  /// boltz(x): exp(-x/kBT)
  /// \Omega: (E-SUM(chem_pot*comp_n))*N
  /// N: number of unit cells in supercell
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
  /// Phi=(\Omega_0-kB*T(ln(SUM(boltz(D\Omega_s)))))/N            
  /// 
  double GrandCanonical::lte_grand_canonical_free_energy(std::ostream& sout) const {

    const SiteExchanger& site_exch = m_site_swaps;
    const ConfigDoF& config_dof = m_configdof;
    GrandCanonicalEvent event = m_event;
    
    double tol = 1e-12;
    
    auto less = [&](const double& A, const double& B) {
      return A < B - tol;
    };
    
    std::map<double, unsigned long, decltype(less)> hist(less);
    
    // no defect case
    hist[0.0] = 1;
    
    double sum_exp = 0.0;
    
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
        double dpot_nrg = event.dpotential_energy() * supercell().volume();
        
        if(dpot_nrg < 0.0) {
          std::cerr << "Error calculating low temperature expansion: \n"
                    << "  Defect lowered the potential energy. Your motif configuration "
                    << "is not the 0K ground state." << std::endl;
          throw std::runtime_error("Error calculating low temperature expansion. Not in the ground state.");
        }
        
        
        auto it = hist.find(dpot_nrg);
        if(it == hist.end()) {
          hist[dpot_nrg] = 1;
        }
        else {
          it->second++;
        }
      }
    }
    
    
    sout << "\n-- Ground state and point defect potential energy details --\n\n";
    
    sout << "T: " << m_condition.temperature() << std::endl;
    sout << "kT: " << 1.0/m_condition.beta() << std::endl;
    sout << "Beta: " << m_condition.beta() << std::endl << std::endl;
    
    sout << std::setw(16) << "N/unitcell" << " "
         << std::setw(16) << "dPE" << " "
         << std::setw(24) << "N*exp(-dPE_i/kT)" << " "
         << std::setw(16) << "dPhi" << " "
         << std::setw(16) << "Phi" << std::endl;
      
    double tsum = 0.0;
    double phi = 0.0;
    double phi_prev;
    for(auto it=hist.rbegin(); it!=hist.rend(); ++it) {
      phi_prev = phi;
      tsum += it->second*exp(-(it->first)*m_condition.beta());
      phi = std::log(tsum)/m_condition.beta()/supercell().volume();
      
      if(almost_equal(it->first, 0.0, tol)) {
        sout << std::setw(16) << "(gs)" << " ";
      }
      else {
        sout << std::setw(16) << std::setprecision(8) << (1.0*it->second)/supercell().volume() << " ";
      }
      sout << std::setw(16) << std::setprecision(8) << it->first << " "
           << std::setw(24) << std::setprecision(8) << it->second*exp(-it->first*m_condition.beta()) << " "
           << std::setw(16) << std::setprecision(8) << phi - phi_prev << " "
           << std::setw(16) << std::setprecision(8) << potential_energy() - phi << std::endl;
      
    }
    
    sout << "Phi_LTE(1): " << std::setprecision(12) << potential_energy() - phi << std::endl;
    
    return potential_energy() - phi;

  }

  /// \brief Print info when a run begins
  void GrandCanonical::print_run_start_info() const {
    sout << "Begin run.  T = " << m_condition.temperature() << "  ";
    
    Eigen::VectorXd param_chem_pot = primclex().composition_axes().param_chem_pot(m_condition.chem_pot());
    for(int i=0; i<param_chem_pot.size(); i++) {
      sout << "param_chem_pot_" << ((char) (i + (int) 'a')) << " = " << param_chem_pot(i) << "  ";
    }
    sout << std::endl;
  }
  
  /// \brief Write results to files
  void GrandCanonical::write_results(Index cond_index) const {
    CASM::write_results(settings(), *this);
    write_conditions_json(settings(), *this, cond_index);
    write_observations(settings(), *this, cond_index);
    write_trajectory(settings(), *this, cond_index);
    //write_pos_trajectory(settings(), *this, cond_index);
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
    Index curr_species = m_site_swaps.sublat_to_mol()[sublat][current_occupant];
    Index new_species = m_site_swaps.sublat_to_mol()[sublat][new_occupant];
    event.set_dcomp_n(curr_species, m_minus_one_comp_n);
    event.set_dcomp_n(new_species, m_plus_one_comp_n);
    
    
    // ---- set dcorr --------------
    
    // Point the Clexulator to the right neighborhood
    m_clexulator.set_nlist(nlist().sites(nlist().unitcell_index(mutating_site)).data());

    // Calculate the change in correlations due to this event
    if(m_all_correlations) {
      m_clexulator.calc_delta_point_corr(sublat,
                                         current_occupant,
                                         new_occupant,
                                         event.dcorr().data());
    }
    else {
      auto begin = m_formation_energy_eci.index().data();
      auto end = begin + m_formation_energy_eci.index().size();
      m_clexulator.calc_restricted_delta_point_corr(sublat,
                                                    current_occupant,
                                                    new_occupant,
                                                    event.dcorr().data(),
                                                    begin,
                                                    end);
    }
    event.dcorr() /= supercell().volume();
    
    // ---- set dformation_energy --------------
    
    event.set_dformation_energy((m_formation_energy_eci * event.dcorr().data()));
    
    
    // ---- set dpotential_energy --------------
    
    event.set_dpotential_energy(event.dformation_energy() - m_condition.exchange_chem_pot(new_species, curr_species)*m_plus_one_comp_n);
    
  }
  
  /// \brief Calculate properties given current conditions
  void GrandCanonical::_update_properties() {
    
    // initialize properties and store pointers to the data strucures
    m_vector_property["corr"] = correlations_vec(m_configdof, supercell(), m_clexulator);
    m_corr = &m_vector_property["corr"];
    
    m_vector_property["comp_n"] = CASM::comp_n(m_configdof, supercell());
    m_comp_n = &m_vector_property["comp_n"];
    
    m_scalar_property["formation_energy"] = m_formation_energy_eci * corr().data();
    m_formation_energy = &m_scalar_property["formation_energy"];
    
    m_scalar_property["potential_energy"] = formation_energy() - primclex().composition_axes().param_composition(comp_n()).dot(m_condition.param_chem_pot());
    m_potential_energy = &m_scalar_property["potential_energy"]; 
    
    if(debug()) {
      
      sout << std::setw(12) << "i" << std::setw(16) << "ECI" << std::setw(16) << "corr" << std::endl;
      
      if(m_all_correlations) {
        for(int i=0; i<corr().size(); ++i) {
          
          double eci = 0.0;
          Index index = find_index(m_formation_energy_eci.index(), i);
          if(index != m_formation_energy_eci.index().size()) {
            eci = m_formation_energy_eci.value()[index];
          }
          
          sout << std::setw(12) << i 
               << std::setw(16) << std::setprecision(8) << eci 
               << std::setw(16) << std::setprecision(8) << corr()[i] << std::endl;
        
        }
      }
      else {
        for(int i=0; i<m_formation_energy_eci.value().size(); ++i) {
          sout << std::setw(12) << m_formation_energy_eci.index()[i] 
               << std::setw(16) << std::setprecision(8) << m_formation_energy_eci.value()[i] 
               << std::setw(16) << std::setprecision(8) << corr()[m_formation_energy_eci.index()[i]] << std::endl;
        
        }
      }
      
      auto origin = primclex().composition_axes().origin();
      auto chem_pot = m_condition.chem_pot();
      auto exchange_chem_pot = chem_pot;
      exchange_chem_pot.array() -= chem_pot(0);
      auto param_chem_pot = m_condition.param_chem_pot();
      auto comp_x = primclex().composition_axes().param_composition(comp_n());
      auto M = primclex().composition_axes().dparam_dmol();
      
      sout << "\n-- Properties --\n\n"

           << "Semi-grand canonical ensemble: \n"
           << "  Thermodynamic potential (per unitcell), phi = -kT*ln(Z)/N \n"
           << "  Partition function, Z = sum_i exp(-N*potential_energy_i/kT) \n"
           << "  parametric composition, comp_x = M * (comp_n - origin) \n"
           << "  parametric chem potential, param_chem_pot = inv(M).transpose() * chem_pot \n"
           << "  potential_energy_i (per unitcell) = formation_energy_i - param_chem_pot*comp_x_i \n\n"
           
           << "components: " << jsonParser(primclex().composition_axes().components()) << "\n"
           << "M:\n" << M << "\n"
           << "origin:\n" << origin.transpose() << "\n"
           << "comp_n: " << comp_n().transpose() << "\n"
           << "chem_pot: " << chem_pot.transpose() << "\n"
           << "exchange_chem_pot: " << exchange_chem_pot.transpose() << "\n"
           << "  exchange_chem_pot*(comp_n - origin): " << exchange_chem_pot.transpose()*(comp_n() - origin) << "\n"
           << "comp_x: " << comp_x.transpose() << "\n"
           << "param_chem_pot: " << param_chem_pot.transpose() << "\n"
           << "  param_chem_pot*comp_x: " << param_chem_pot.dot(comp_x)  << "\n"
           << "formation_energy: " << formation_energy() << "\n"
           << "  formation_energy - exchange_chem_pot*(comp_n - origin): " << formation_energy() - exchange_chem_pot.transpose()*(comp_n() - origin) << "\n"
           << "  formation_energy - param_chem_pot*(comp_x): " << formation_energy() - param_chem_pot.dot(comp_x) << "\n"
           << "potential_energy: " << potential_energy() << "\n" << std::endl;
    }
    
  }
  
  /// \brief Select initial motif configuration
  ///
  /// \param motif_configname If "auto", use 0K ground state at given mu; else 
  ///        use configuration with given name
  const Configuration& GrandCanonical::_select_motif(
      std::string motif_configname, 
      PrimClex& primclex, 
      const GrandCanonicalConditions& cond,
      std::ostream& _sout) const {
    if(motif_configname == "auto") {
      
      std::cout << "Searching for minimum potential energy motif..." << std::endl;
      
      double tol = 1e-6;
      auto compare = [&](double A, double B) {
        return A < B - tol;
      };
      
      ConfigIO::Clex clex(primclex.global_clexulator(), primclex.global_eci("formation_energy"));
      
      std::multimap<double, const Configuration*, decltype(compare)> configmap(compare);
      for(auto it=primclex.config_begin(); it!=primclex.config_end(); ++it) {
        configmap.insert(std::make_pair(clex(*it) - cond.chem_pot().dot(CASM::comp_n(*it)), &(*it)));
      }
      
      const Configuration& min_config = *(configmap.begin()->second);
      double min_potential_energy = configmap.begin()->first;
      auto eq_range = configmap.equal_range(min_potential_energy);
      if(std::distance(eq_range.first, eq_range.second) > 1) {
        _sout << "Warning: Found degenerate ground states with potential energy: " 
              << std::setprecision(8) << min_potential_energy << std::endl;
        for(auto it=eq_range.first; it!=eq_range.second; ++it) {
          _sout << "  " << it->second->name() << std::endl;
        }
        _sout << "Choosing: " << min_config.name() << std::endl;
      }
      else {
        _sout << "Found: " << min_config.name() << " with potential energy: " 
              << std::setprecision(8) << min_potential_energy << std::endl;
      }
      
      return min_config;
      
      
    }
    else {
      return m_primclex.configuration(motif_configname);
    }
  }
  

  
}

