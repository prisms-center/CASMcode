
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/monte_carlo/MonteCarloEnum_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/Norm.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

namespace CASM {


  GrandCanonical::GrandCanonical(PrimClex &primclex, const GrandCanonicalSettings &settings, Log &log):
    MonteCarloEnum<MonteCarlo>(primclex, settings, log),
    m_site_swaps(supercell()),
    m_formation_energy_clex(primclex, settings.formation_energy(primclex)),
    m_all_correlations(settings.all_correlations()),
    m_event(primclex.composition_axes().components().size(), _clexulator().corr_size()) {

    const auto &desc = m_formation_energy_clex.desc();

    _log().construct("Grand Canonical Monte Carlo");
    _log() << "project: " << this->primclex().get_path() << "\n";
    _log() << "formation_energy cluster expansion: " << desc.name << "\n";
    _log() << std::setw(16) << "property: " << desc.property << "\n";
    _log() << std::setw(16) << "calctype: " << desc.calctype << "\n";
    _log() << std::setw(16) << "ref: " << desc.ref << "\n";
    _log() << std::setw(16) << "bset: " << desc.bset << "\n";
    _log() << std::setw(16) << "eci: " << desc.eci << "\n";
    _log() << "supercell: \n" << supercell().get_transf_mat() << "\n";
    if(enum_configs()) {
      _log() << "configuration enumeration: " << std::boolalpha << true << "\n";
      _log() << "  check: " << enum_check_args() << "\n";
      _log() << "  metric: " << enum_metric_args() << "\n";
    }
    _log() << std::endl;

    // set the SuperNeighborList...
    set_nlist();

    // If the simulation is big enough, use delta cluster functions;
    // else, calculate all cluster functions
    m_use_deltas = !nlist().overlaps();

    reset(_initial_configdof(settings, _log()));

    _log().set("Initial Conditions");
    _log() << settings.initial_conditions() << std::endl << std::endl;
    m_condition = settings.initial_conditions();

    _update_properties();

  }

  /// \brief Return number of steps per pass. Equals number of sites with variable occupation.
  Index GrandCanonical::steps_per_pass() const {
    return m_site_swaps.variable_sites().size();
  }


  /// \brief Return current conditions
  const GrandCanonical::CondType &GrandCanonical::conditions() const {
    return m_condition;
  }


  /// \brief Set conditions and clear previously collected data
  void GrandCanonical::set_conditions(const GrandCanonicalConditions &new_conditions) {
    _log().set("Conditions");
    _log() << new_conditions << std::endl << std::endl;

    m_condition = new_conditions;

    clear_samples();
    _update_properties();

    return;
  }

  /// \brief Set configdof and clear previously collected data
  void GrandCanonical::set_configdof(const ConfigDoF &configdof, const std::string &msg) {
    _log().set("DoF");
    if(!msg.empty()) {
      _log() << msg << "\n";
    }
    _log() << std::endl;

    reset(configdof);
    _update_properties();
  }

  /// \brief Propose a new event, calculate delta properties, and return reference to it
  ///
  /// Randomly picks a site that's allowed more than one occupant, and randomly picks what occupant it
  /// changes to. Then calculates delta properties associated with that change.
  ///
  const GrandCanonical::EventType &GrandCanonical::propose() {

    // Randomly pick a site that's allowed more than one occupant
    Index random_variable_site = _mtrand().randInt(m_site_swaps.variable_sites().size() - 1);

    // Determine what that site's linear index is and what the sublattice index is
    Index mutating_site = m_site_swaps.variable_sites()[random_variable_site];
    Index sublat = m_site_swaps.sublat()[random_variable_site];

    // Determine the current occupant of the mutating site
    int current_occupant = configdof().occ(mutating_site);

    // Randomly pick a new occupant for the mutating site
    const std::vector<int> &possible_mutation = m_site_swaps.possible_swap()[sublat][current_occupant];
    int new_occupant = possible_mutation[_mtrand().randInt(possible_mutation.size() - 1)];

    // Update delta properties in m_event
    _update_deltas(m_event, mutating_site, sublat, current_occupant, new_occupant);

    if(debug()) {
      const auto &site_occ = primclex().get_prim().basis[sublat].site_occupant();
      _log().custom("Propose event");

      _log()  << "  Mutating site (linear index): " << mutating_site << "\n"
              << "  Mutating site (b, i, j, k): " << supercell().uccoord(mutating_site) << "\n"
              << "  Current occupant: " << current_occupant << " (" << site_occ[current_occupant].name << ")\n"
              << "  Proposed occupant: " << new_occupant << " (" << site_occ[new_occupant].name << ")\n\n"

              << "  beta: " << m_condition.beta() << "\n"
              << "  T: " << m_condition.temperature() << std::endl;

      _log() << std::setw(12) << "i" << std::setw(16) << "ECI" << std::setw(16) << "dCorr" << std::endl;

      if(m_all_correlations) {
        for(int i = 0; i < m_event.dCorr().size(); ++i) {

          double eci = 0.0;
          Index index = find_index(_eci().index(), i);
          if(index != _eci().index().size()) {
            eci = _eci().value()[index];
          }

          _log() << std::setw(12) << i
                 << std::setw(16) << std::setprecision(8) << eci
                 << std::setw(16) << std::setprecision(8) << m_event.dCorr()[i] << std::endl;

        }
      }
      else {
        for(int i = 0; i < _eci().value().size(); ++i) {
          _log() << std::setw(12) << _eci().index()[i]
                 << std::setw(16) << std::setprecision(8) << _eci().value()[i]
                 << std::setw(16) << std::setprecision(8) << m_event.dCorr()[_eci().index()[i]] << std::endl;

        }
      }

      auto origin = primclex().composition_axes().origin();
      auto exchange_chem_pot = m_condition.exchange_chem_pot();
      auto param_chem_pot = m_condition.param_chem_pot();
      auto Mpinv = primclex().composition_axes().dparam_dmol();
      auto V = supercell().volume();
      Index curr_species = m_site_swaps.sublat_to_mol()[sublat][current_occupant];
      Index new_species = m_site_swaps.sublat_to_mol()[sublat][new_occupant];

      _log() << "  components: " << jsonParser(primclex().composition_axes().components()) << "\n"
             << "  d(N): " << m_event.dN().transpose() << "\n"
             << "    dx_dn: \n" << Mpinv << "\n"
             << "    param_chem_pot.transpose() * dx_dn: \n" << param_chem_pot.transpose()*Mpinv << "\n"
             << "    param_chem_pot.transpose() * dx_dn * dN: " << param_chem_pot.transpose()*Mpinv *m_event.dN().cast<double>() << "\n"
             << "  d(Nunit * param_chem_pot * x): " << exchange_chem_pot(new_species, curr_species) << "\n"
             << "  d(Ef): " << m_event.dEf() << "\n"
             << "  d(Epot): " << m_event.dEf() - exchange_chem_pot(new_species, curr_species) << "\n"
             << std::endl;


    }

    return m_event;
  }

  /// \brief Based on a random number, decide if the change in energy from the proposed event is low enough to be accepted.
  bool GrandCanonical::check(const GrandCanonicalEvent &event) {

    if(event.dEpot() < 0.0) {

      if(debug()) {
        _log().custom("Check event");
        _log() << "Probability to accept: 1.0\n" << std::endl;
      }
      return true;
    }

    double rand = _mtrand().rand53();
    double prob = exp(-event.dEpot() * m_condition.beta());

    if(debug()) {
      _log().custom("Check event");
      _log() << "Probability to accept: " << prob << "\n"
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
  void GrandCanonical::accept(const EventType &event) {

    if(debug()) {
      _log().custom("Accept Event");
      _log() << std::endl;
    }

    // First apply changes to configuration (just a single occupant change)
    _configdof().occ(event.occupational_change().site_index()) = event.occupational_change().to_value();

    // Next update all properties that changed from the event
    _formation_energy() += event.dEf() / supercell().volume();
    _potential_energy() += event.dEpot() / supercell().volume();
    _corr() += event.dCorr() / supercell().volume();
    _comp_n() += event.dN().cast<double>() / supercell().volume();

    // If enumerating, try to save current configuration
    if(enum_configs()) {
      enum_insert(config());
    }

    return;
  }

  /// \brief Nothing needs to be done to reject a GrandCanonicalEvent
  void GrandCanonical::reject(const EventType &event) {
    if(debug()) {
      _log().custom("Reject Event");
      _log() << std::endl;
    }
    return;
  }

  /// \brief Calculate the single spin flip low temperature expansion of the grand canonical potential
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
  double GrandCanonical::lte_grand_canonical_free_energy() const {

    const SiteExchanger &site_exch = m_site_swaps;
    const ConfigDoF &config_dof = configdof();
    GrandCanonicalEvent event = m_event;

    double tol = 1e-12;

    auto less = [&](const double & A, const double & B) {
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
      const auto &possible = site_exch.possible_swap()[sublat][current_occupant];
      for(auto new_occ_it = possible.begin(); new_occ_it != possible.end(); ++new_occ_it) {

        _update_deltas(event, mutating_site, sublat, current_occupant, *new_occ_it);

        //save the result
        double dpot_nrg = event.dEpot();

        if(dpot_nrg < 0.0) {
          Log &err_log = default_err_log();
          err_log.error<Log::standard>("Calculating low temperature expansion");
          err_log << "  Defect lowered the potential energy. Your motif configuration "
                  << "is not the 0K ground state.\n" << std::endl;
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

    _log().results("Ground state and point defect potential energy details");
    _log() << "T: " << m_condition.temperature() << std::endl;
    _log() << "kT: " << 1.0 / m_condition.beta() << std::endl;
    _log() << "Beta: " << m_condition.beta() << std::endl << std::endl;

    _log() << std::setw(16) << "N/unitcell" << " "
           << std::setw(16) << "dPE" << " "
           << std::setw(24) << "N*exp(-dPE/kT)" << " "
           << std::setw(16) << "dphi" << " "
           << std::setw(16) << "phi" << std::endl;

    double tsum = 0.0;
    double phi = 0.0;
    double phi_prev;
    for(auto it = hist.rbegin(); it != hist.rend(); ++it) {
      phi_prev = phi;
      tsum += it->second * exp(-(it->first) * m_condition.beta());
      phi = std::log(tsum) / m_condition.beta() / supercell().volume();

      if(almost_equal(it->first, 0.0, tol)) {
        _log() << std::setw(16) << "(gs)" << " ";
      }
      else {
        _log() << std::setw(16) << std::setprecision(8) << (1.0 * it->second) / supercell().volume() << " ";
      }
      _log() << std::setw(16) << std::setprecision(8) << it->first << " "
             << std::setw(24) << std::setprecision(8) << it->second *exp(-it->first * m_condition.beta()) << " "
             << std::setw(16) << std::setprecision(8) << phi - phi_prev << " "
             << std::setw(16) << std::setprecision(8) << potential_energy() - phi << std::endl;

    }

    _log() << "phi_LTE(1): " << std::setprecision(12) << potential_energy() - phi << std::endl << std::endl;

    return potential_energy() - phi;

  }

  /// \brief Write results to files
  void GrandCanonical::write_results(Index cond_index) const {
    CASM::write_results(settings(), *this, _log());
    write_conditions_json(settings(), *this, cond_index, _log());
    write_observations(settings(), *this, cond_index, _log());
    write_trajectory(settings(), *this, cond_index, _log());
    //write_pos_trajectory(settings(), *this, cond_index);

    // If enumerating, try to save current configuration
    if(enum_configs()) {
      save_configs();
    }
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

    for(int i = 0; i < event.dN().size(); ++i) {
      event.set_dN(i, 0);
    }
    Index curr_species = m_site_swaps.sublat_to_mol()[sublat][current_occupant];
    Index new_species = m_site_swaps.sublat_to_mol()[sublat][new_occupant];
    event.set_dN(curr_species, -1);
    event.set_dN(new_species, 1);


    // ---- set dcorr --------------

    // Point the Clexulator to the right neighborhood
    _clexulator().set_nlist(nlist().sites(nlist().unitcell_index(mutating_site)).data());

    if(m_use_deltas) {

      // Calculate the change in correlations due to this event
      if(m_all_correlations) {
        _clexulator().calc_delta_point_corr(sublat,
                                            current_occupant,
                                            new_occupant,
                                            event.dCorr().data());
      }
      else {
        auto begin = _eci().index().data();
        auto end = begin + _eci().index().size();
        _clexulator().calc_restricted_delta_point_corr(sublat,
                                                       current_occupant,
                                                       new_occupant,
                                                       event.dCorr().data(),
                                                       begin,
                                                       end);
      }
    }
    else {

      int &occ = _configdof().occ(event.occupational_change().site_index());
      int from_value = occ;

      // Apply changes to configuration (just a single occupant change)
      occ = event.occupational_change().to_value();

      // Calculate the change in correlations due to this event
      event.dCorr() = (correlations_vec(_configdof(), supercell(), _clexulator()) - corr()) * supercell().volume();

      // Unapply changes
      occ = from_value;
    }

    // ---- set dformation_energy --------------

    event.set_dEf(_eci() * event.dCorr().data());


    // ---- set dpotential_energy --------------

    event.set_dEpot(event.dEf() - m_condition.exchange_chem_pot(new_species, curr_species));

  }

  /// \brief Calculate properties given current conditions
  void GrandCanonical::_update_properties() {

    // initialize properties and store pointers to the data strucures
    _vector_properties()["corr"] = correlations_vec(_configdof(), supercell(), _clexulator());
    m_corr = &_vector_property("corr");

    _vector_properties()["comp_n"] = CASM::comp_n(_configdof(), supercell());
    m_comp_n = &_vector_property("comp_n");

    _scalar_properties()["formation_energy"] = _eci() * corr().data();
    m_formation_energy = &_scalar_property("formation_energy");

    _scalar_properties()["potential_energy"] = formation_energy() - primclex().composition_axes().param_composition(comp_n()).dot(m_condition.param_chem_pot());
    m_potential_energy = &_scalar_property("potential_energy");

    if(debug()) {

      _log().custom("Calculate correlations");
      _log() << std::setw(12) << "i" << std::setw(16) << "ECI" << std::setw(16) << "corr" << std::endl;

      if(m_all_correlations) {
        for(int i = 0; i < corr().size(); ++i) {

          double eci = 0.0;
          Index index = find_index(_eci().index(), i);
          if(index != _eci().index().size()) {
            eci = _eci().value()[index];
          }

          _log() << std::setw(12) << i
                 << std::setw(16) << std::setprecision(8) << eci
                 << std::setw(16) << std::setprecision(8) << corr()[i] << std::endl;

        }
      }
      else {
        for(int i = 0; i < _eci().value().size(); ++i) {
          _log() << std::setw(12) << _eci().index()[i]
                 << std::setw(16) << std::setprecision(8) << _eci().value()[i]
                 << std::setw(16) << std::setprecision(8) << corr()[_eci().index()[i]] << std::endl;

        }
      }
      _log() << std::endl;

      auto origin = primclex().composition_axes().origin();
      auto exchange_chem_pot = m_condition.exchange_chem_pot();
      auto param_chem_pot = m_condition.param_chem_pot();
      auto comp_x = primclex().composition_axes().param_composition(comp_n());
      auto M = primclex().composition_axes().dmol_dparam();

      _log().custom("Calculate properties");
      _log() << "Semi-grand canonical ensemble: \n"
             << "  Thermodynamic potential (per unitcell), phi = -kT*ln(Z)/N \n"
             << "  Partition function, Z = sum_i exp(-N*potential_energy_i/kT) \n"
             << "  composition, comp_n = origin + M * comp_x \n"
             << "  parametric chemical potential, param_chem_pot = M.transpose() * chem_pot \n"
             << "  potential_energy (per unitcell) = formation_energy - param_chem_pot*comp_x \n\n"

             << "components: " << jsonParser(primclex().composition_axes().components()) << "\n"
             << "M:\n" << M << "\n"
             << "origin: " << origin.transpose() << "\n"
             << "comp_n: " << comp_n().transpose() << "\n"
             << "comp_x: " << comp_x.transpose() << "\n"
             << "param_chem_pot: " << param_chem_pot.transpose() << "\n"
             << "  param_chem_pot*comp_x: " << param_chem_pot.dot(comp_x)  << "\n"
             << "formation_energy: " << formation_energy() << "\n"
             << "  formation_energy - param_chem_pot*comp_x: " << formation_energy() - param_chem_pot.dot(comp_x) << "\n"
             << "potential_energy: " << potential_energy() << "\n" << std::endl;
    }

  }

  /// \brief Select initial configdof
  ///
  /// - "motif"/"configname": If "auto", use 0K ground state at given mu; else
  ///   use configuration with given name
  /// - "motif"/"configdof": Open configdof file with given name; must be for
  ///   the specified supercell.
  /// - Also prints messages describing what ConfigDoF is being used
  ConfigDoF GrandCanonical::_initial_configdof(
    const GrandCanonicalSettings &settings,
    Log &_log) {

    _log.set("Initial DoF");

    if(settings.is_motif_configname()) {

      std::string motif_configname = settings.motif_configname();
      _log << "motif configname: " << motif_configname << "\n";

      if(motif_configname == "default") {
        _log << "using configuration with default occupation..." << std::endl;
        return Configuration(_supercell(), jsonParser(), Array<int>(_supercell().num_sites(), 0)).configdof();
      }
      else if(motif_configname == "auto" || motif_configname == "restricted_auto") {

        _log << "searching for minimum potential energy motif..." << std::endl;

        double tol = 1e-6;
        auto compare = [&](double A, double B) {
          return A < B - tol;
        };

        GrandCanonicalConditions cond = settings.initial_conditions();

        _log << "using conditions: \n";
        _log << cond << std::endl;

        std::multimap<double, const Configuration *, decltype(compare)> configmap(compare);
        for(auto it = primclex().config_begin(); it != primclex().config_end(); ++it) {
          configmap.insert(std::make_pair(_eci() * correlations(*it, _clexulator()) - cond.param_chem_pot().dot(CASM::comp(*it)), &(*it)));
        }

        if(motif_configname == "restricted_auto") {

          // used to check if configurations can fill the monte carlo supercell
          const Lattice &scel_lat = supercell().get_real_super_lattice();
          auto begin = primclex().get_prim().factor_group().begin();
          auto end = primclex().get_prim().factor_group().end();

          // save iterators pointing to configs that will fill the supercell
          std::vector<decltype(configmap)::const_iterator> allowed;

          // iterate through groups of degenerate configs
          auto next_it = configmap.begin();
          while(next_it != configmap.end()) {
            auto eq_range = configmap.equal_range(next_it->first);

            // save allowed configs
            for(auto it = eq_range.first; it != eq_range.second; ++it) {
              const Lattice &motif_lat = it->second->get_supercell().get_real_super_lattice();
              if(is_supercell(scel_lat, motif_lat, begin, end, TOL).first != end) {
                allowed.push_back(it);
              }
            }

            // if some found, break
            if(allowed.size()) {
              break;
            }

            // else, continue to next group
            next_it = eq_range.second;
          }

          if(!allowed.size()) {
            _log << "Found no enumerated configurations that will fill the supercell\n";
            _log << "using configuration with default occupation..." << std::endl;
            return Configuration(_supercell(), jsonParser(), Array<int>(_supercell().num_sites(), 0)).configdof();
          }

          if(allowed.size() > 1) {
            _log << "Warning: Found degenerate allowed configurations with potential energy: "
                 << std::setprecision(8) << allowed[0]->first << std::endl;
            for(auto it = allowed.begin(); it != allowed.end(); ++it) {
              _log << "  " << (*it)->second->name() << std::endl;
            }
            _log << "using: " << allowed[0]->second->name() << "\n" << std::endl;
          }
          else {
            _log << "using: " << allowed[0]->second->name() << " with potential energy: "
                 << std::setprecision(8) << allowed[0]->first << "\n" << std::endl;
          }

          return fill_supercell(_supercell(), *(allowed[0]->second)).configdof();
        }
        else {
          const Configuration &min_config = *(configmap.begin()->second);
          double min_potential_energy = configmap.begin()->first;
          auto eq_range = configmap.equal_range(min_potential_energy);

          if(std::distance(eq_range.first, eq_range.second) > 1) {
            _log << "Warning: Found degenerate ground states with potential energy: "
                 << std::setprecision(8) << min_potential_energy << std::endl;
            for(auto it = eq_range.first; it != eq_range.second; ++it) {
              _log << "  " << it->second->name() << std::endl;
            }
            _log << "using: " << min_config.name() << "\n" << std::endl;
          }
          else {
            _log << "using: " << min_config.name() << " with potential energy: "
                 << std::setprecision(8) << min_potential_energy << "\n" << std::endl;
          }

          return fill_supercell(_supercell(), min_config).configdof();
        }
      }
      else {
        _log << "using configation: " << motif_configname << "\n" << std::endl;
        return fill_supercell(_supercell(), primclex().configuration(motif_configname)).configdof();
      }

    }
    else if(settings.is_motif_configdof()) {
      _log << "motif configdof: " << settings.motif_configdof_path() << "\n";
      _log << "using configdof: " << settings.motif_configdof_path() << "\n" << std::endl;
      return settings.motif_configdof();
    }
    else {
      throw std::runtime_error("Error: Must specify motif \"configname\" or \"configdof\"");
    }
  }

}

