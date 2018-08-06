
#include "casm/monte_carlo/canonical/Canonical.hh"
#include "casm/monte_carlo/MonteCarloEnum_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/Norm.hh"
#include "casm/monte_carlo/canonical/CanonicalIO.hh"
#include "casm/monte_carlo/MonteIO_impl.hh"

namespace CASM {
  namespace Monte {

    const ENSEMBLE Canonical::ensemble = ENSEMBLE::Canonical;

    /// \brief Constructs a Canonical object and prepares it for running based on MonteSettings
    ///
    /// - Does not set 'state': conditions or ConfigDoF
    Canonical::Canonical(PrimClex &primclex, const CanonicalSettings &settings, Log &log):
      MonteCarlo(primclex, settings, log),
      m_formation_energy_clex(primclex, settings.formation_energy(primclex)),
      m_convert(_supercell()),
      m_cand(m_convert),
      m_all_correlations(settings.all_correlations()),
      m_occ_loc(m_convert, m_cand),
      m_event(primclex.composition_axes().components().size(), _clexulator().corr_size()) {

      const auto &desc = m_formation_energy_clex.desc();

      // set the SuperNeighborList...
      set_nlist();

      // If the simulation is big enough, use delta cluster functions;
      // else, calculate all cluster functions
      m_use_deltas = !nlist().overlaps();

      _log().construct("Canonical Monte Carlo");
      _log() << "project: " << this->primclex().get_path() << "\n";
      _log() << "formation_energy cluster expansion: " << desc.name << "\n";
      _log() << std::setw(16) << "property: " << desc.property << "\n";
      _log() << std::setw(16) << "calctype: " << desc.calctype << "\n";
      _log() << std::setw(16) << "ref: " << desc.ref << "\n";
      _log() << std::setw(16) << "bset: " << desc.bset << "\n";
      _log() << std::setw(16) << "eci: " << desc.eci << "\n";
      _log() << "supercell: \n" << supercell().get_transf_mat() << "\n";
      _log() << "use_deltas: " << std::boolalpha << m_use_deltas << "\n";
      _log() << "\nSampling: \n";
      _log() << std::setw(24) << "quantity" << std::setw(24) << "requested_precision" << "\n";
      for(auto it = samplers().begin(); it != samplers().end(); ++it) {
        _log() << std::setw(24) << it->first;
        if(it->second->must_converge()) {
          _log() << std::setw(24) << it->second->requested_precision() << std::endl;
        }
        else {
          _log() << std::setw(24) << "none" << std::endl;
        }
      }
      _log() << "\nautomatic convergence mode?: " << std::boolalpha << must_converge() << std::endl;
      _log() << std::endl;

      _log() << std::pair<const OccCandidateList &, const Conversions &>(m_cand, m_convert) << std::endl;

    }

    /// \brief Return number of steps per pass. Equals number of sites with variable occupation.
    Index Canonical::steps_per_pass() const {
      return m_occ_loc.size();
    }


    /// \brief Return current conditions
    const Canonical::CondType &Canonical::conditions() const {
      return m_condition;
    }


    /// \brief Set conditions and clear previously collected data
    void Canonical::set_conditions(const CanonicalConditions &new_conditions) {
      _log().set("Conditions");
      _log() << new_conditions << std::endl << std::endl;

      m_condition = new_conditions;

      reset(_enforce_conditions(configdof()));
      _update_properties();

      return;
    }

    /// \brief Set configdof and clear previously collected data
    void Canonical::set_configdof(const ConfigDoF &configdof, const std::string &msg) {
      _log().set("DoF");
      if(!msg.empty()) {
        _log() << msg << "\n";
      }
      _log() << std::endl;

      reset(_enforce_conditions(configdof));
      _update_properties();
    }

    /// \brief Set configdof and conditions and clear previously collected data
    ///
    /// \returns Specified ConfigDoF and configname (or configdof path)
    ///
    std::pair<ConfigDoF, std::string> Canonical::set_state(
      const CanonicalConditions &new_conditions,
      const CanonicalSettings &settings) {

      _log().set("Conditions");
      _log() << new_conditions << std::endl;

      m_condition = new_conditions;

      ConfigDoF configdof;
      std::string configname;

      if(settings.is_motif_configname()) {

        configname = settings.motif_configname();

        if(configname == "default") {
          configdof = _default_motif();
        }
        else if(configname == "auto") {
          std::tie(configdof, configname) = _auto_motif(new_conditions);
        }
        else if(configname == "restricted_auto") {
          std::tie(configdof, configname) = _restricted_auto_motif(new_conditions);
        }
        else {
          configdof = _configname_motif(configname);
        }

      }
      else if(settings.is_motif_configdof()) {
        _log().set("DoF");
        _log() << "motif configdof: " << settings.motif_configdof_path() << "\n";
        _log() << "using configdof: " << settings.motif_configdof_path() << "\n" << std::endl;
        configdof = settings.motif_configdof();
        configname = settings.motif_configdof_path().string();
      }
      else {
        throw std::runtime_error("Error: Must specify motif \"configname\" or \"configdof\"");
      }

      reset(_enforce_conditions(configdof));
      _update_properties();

      return std::make_pair(configdof, configname);
    }

    /// \brief Set configdof and conditions and clear previously collected data
    void Canonical::set_state(const CanonicalConditions &new_conditions,
                              const ConfigDoF &configdof,
                              const std::string &msg) {
      _log().set("Conditions");
      _log() << new_conditions << std::endl << std::endl;

      m_condition = new_conditions;

      _log().set("DoF");
      if(!msg.empty()) {
        _log() << msg << "\n";
      }
      _log() << std::endl;

      reset(_enforce_conditions(configdof));
      _update_properties();

      return;
    }


    /// \brief Propose a new event, calculate delta properties, and return reference to it
    ///
    /// Randomly picks a site that's allowed more than one occupant, and randomly picks what occupant it
    /// changes to. Then calculates delta properties associated with that change.
    ///
    const Canonical::EventType &Canonical::propose() {
      m_occ_loc.propose_canonical(m_event.occ_event(), m_cand.canonical_swap(), _mtrand());
      _update_deltas(m_event);
      return m_event;
    }

    /// \brief Based on a random number, decide if the change in energy from the proposed event is low enough to be accepted.
    bool Canonical::check(const CanonicalEvent &event) {

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
    void Canonical::accept(const EventType &event) {

      if(debug()) {
        _log().custom("Accept Event");
        _log() << std::endl;
      }

      // Apply occ mods && update occ locations table
      m_occ_loc.apply(event.occ_event(), _configdof());

      // Next update all properties that changed from the event
      _formation_energy() += event.dEf() / supercell().volume();
      _potential_energy() += event.dEpot() / supercell().volume();
      _corr() += event.dCorr() / supercell().volume();
      _comp_n() += event.dN().cast<double>() / supercell().volume();

      return;
    }

    /// \brief Nothing needs to be done to reject a CanonicalEvent
    void Canonical::reject(const EventType &event) {
      if(debug()) {
        _log().custom("Reject Event");
        _log() << std::endl;
      }
      return;
    }

    /// \brief Write results to files
    void Canonical::write_results(Index cond_index) const {
      CASM::write_results(settings(), *this, _log());
      write_conditions_json(settings(), *this, cond_index, _log());
      write_observations(settings(), *this, cond_index, _log());
      write_trajectory(settings(), *this, cond_index, _log());
      //write_pos_trajectory(settings(), *this, cond_index);
    }

    /// \brief Get potential energy
    ///
    /// - if(&config == &this->config()) { return potential_energy(); }, else
    ///   calculate potential_energy = formation_energy
    double Canonical::potential_energy(const Configuration &config) const {
      //if(&config == &this->config()) { return potential_energy(); }

      auto corr = correlations(config, _clexulator());
      return _eci() * corr.data();
    }

    void Canonical::_set_nlist(Index l) const {
      _clexulator().set_nlist(nlist().sites(nlist().unitcell_index(l)).data());
    };

    void Canonical::_calc_delta_point_corr(Index l, int new_occ, Eigen::VectorXd &dCorr_comp) const {

      int sublat = _config().get_b(l);
      int curr_occ = _configdof().occ(l);
      _set_nlist(l);

      // Calculate the change in correlations due to this event
      if(m_use_deltas) {
        if(m_all_correlations) {
          _clexulator().calc_delta_point_corr(sublat,
                                              curr_occ,
                                              new_occ,
                                              dCorr_comp.data());
        }
        else {
          auto begin = _eci().index().data();
          auto end = begin + _eci().index().size();
          _clexulator().calc_restricted_delta_point_corr(sublat,
                                                         curr_occ,
                                                         new_occ,
                                                         dCorr_comp.data(),
                                                         begin,
                                                         end);
        }
      }
      else {
        Eigen::VectorXd before { Eigen::VectorXd::Zero(dCorr_comp.size()) };
        Eigen::VectorXd after { Eigen::VectorXd::Zero(dCorr_comp.size()) };

        // Calculate the change in points correlations due to this event
        if(m_all_correlations) {

          // Calculate before
          _clexulator().calc_point_corr(sublat, before.data());

          // Apply change
          _configdof().occ(l) = new_occ;

          // Calculate after
          _clexulator().calc_point_corr(sublat, after.data());
        }
        else {
          auto begin = _eci().index().data();
          auto end = begin + _eci().index().size();

          // Calculate before
          _clexulator().calc_restricted_point_corr(sublat, before.data(), begin, end);

          // Apply change
          _configdof().occ(l) = new_occ;

          // Calculate after
          _clexulator().calc_restricted_point_corr(sublat, after.data(), begin, end);

        }
        dCorr_comp = after - before;

        // Unapply changes
        _configdof().occ(l) = curr_occ;
      }
    }

    /// \brief Calculate delta correlations for an event
    void Canonical::_set_dCorr(CanonicalEvent &event) const {

      const OccEvent &e = event.occ_event();
      const OccTransform &f_a = e.occ_transform[0];
      const OccTransform &f_b = e.occ_transform[1];

      int curr_occ_a = _configdof().occ(f_a.l);
      Index new_occ_a = m_convert.occ_index(f_a.asym, f_a.to_species);
      Index new_occ_b = m_convert.occ_index(f_b.asym, f_b.to_species);

      Eigen::VectorXd dCorr_comp { Eigen::VectorXd::Zero(event.dCorr().size()) };


      // Point the Clexulator to the right neighborhood and right ConfigDoF
      _clexulator().set_config_occ(_configdof().occupation().begin());

      // calc dCorr for first site
      _calc_delta_point_corr(f_a.l, new_occ_a, event.dCorr());

      // change occ on first site
      _configdof().occ(f_a.l) = new_occ_a;

      // calc dCorr for second site
      _calc_delta_point_corr(f_b.l, new_occ_b, dCorr_comp);
      event.dCorr() += dCorr_comp;

      // unchange occ on first site
      _configdof().occ(f_a.l) = curr_occ_a;

      if(debug()) {
        _print_correlations(event.dCorr(), "delta correlations", "dCorr", m_all_correlations);
      }
    }

    /// \brief Print correlations to _log()
    void Canonical::_print_correlations(
      const Eigen::VectorXd &corr,
      std::string title,
      std::string colheader,
      bool all_correlations) const {

      _log().calculate(title);
      _log() << std::setw(12) << "i"
             << std::setw(16) << "ECI"
             << std::setw(16) << colheader
             << std::endl;

      for(int i = 0; i < corr.size(); ++i) {

        double eci = 0.0;
        bool calculated = true;
        Index index = find_index(_eci().index(), i);
        if(index != _eci().index().size()) {
          eci = _eci().value()[index];
        }
        if(!all_correlations && index == _eci().index().size()) {
          calculated = false;
        }

        _log() << std::setw(12) << i
               << std::setw(16) << std::setprecision(8) << eci;
        if(calculated) {
          _log() << std::setw(16) << std::setprecision(8) << corr[i];
        }
        else {
          _log() << std::setw(16) << "unknown";
        }
        _log() << std::endl;

      }
      _log() << std::endl;
    }

    /// \brief Update delta properties in 'event'
    void Canonical::_update_deltas(CanonicalEvent &event) const {

      // ---- set dcorr --------------
      _set_dCorr(event);

      // ---- set dformation_energy --------------
      event.set_dEf(_eci() * event.dCorr().data());
    }

    /// \brief Calculate properties given current conditions
    void Canonical::_update_properties() {

      // initialize properties and store pointers to the data strucures
      _vector_properties()["corr"] = correlations_vec(_configdof(), supercell(), _clexulator());
      m_corr = &_vector_property("corr");

      _vector_properties()["comp_n"] = CASM::comp_n(_configdof(), supercell());
      m_comp_n = &_vector_property("comp_n");

      _scalar_properties()["formation_energy"] = _eci() * corr().data();
      m_formation_energy = &_scalar_property("formation_energy");

      _scalar_properties()["potential_energy"] = formation_energy();
      m_potential_energy = &_scalar_property("potential_energy");

      if(debug()) {

        _print_correlations(corr(), "correlations", "corr", m_all_correlations);

        auto origin = primclex().composition_axes().origin();
        auto comp_x = primclex().composition_axes().param_composition(comp_n());
        auto M = primclex().composition_axes().dmol_dparam();

        _log().custom("Calculate properties");
        _log() << "Canonical ensemble: \n"
               << "  Thermodynamic potential (per unitcell), phi = -kT*ln(Z)/N \n"
               << "  Partition function, Z = sum_i exp(-N*potential_energy_i/kT) \n"
               << "  composition, comp_n = origin + M * comp_x \n"
               << "  potential_energy (per unitcell) = formation_energy \n\n"

               << "components: " << jsonParser(primclex().composition_axes().components()) << "\n"
               << "M:\n" << M << "\n"
               << "origin: " << origin.transpose() << "\n"
               << "comp_n: " << comp_n().transpose() << "\n"
               << "comp_x: " << comp_x.transpose() << "\n"
               << "formation_energy: " << formation_energy() << "\n"
               << "potential_energy: " << potential_energy() << "\n" << std::endl;
      }

    }

    /// \brief Generate supercell filling ConfigDoF from default configuration
    ConfigDoF Canonical::_default_motif() const {
      _log().set("DoF");
      _log() << "motif configname: default\n";
      _log() << "using configuration with default occupation...\n" << std::endl;
      return Configuration(_supercell(), jsonParser(), Array<int>(_supercell().num_sites(), 0)).configdof();
    }

    /// \brief Generate minimum potential energy ConfigDoF
    ///
    /// Raises exception if it doesn't tile the supercell
    std::pair<ConfigDoF, std::string> Canonical::_auto_motif(const CanonicalConditions &cond) const {
      throw std::runtime_error("Canonical Monte Carlo 'auto' motif is not implemented yet");
    }

    /// \brief Generate minimum potential energy ConfigDoF for this supercell
    std::pair<ConfigDoF, std::string> Canonical::_restricted_auto_motif(const CanonicalConditions &cond) const {
      throw std::runtime_error("Canonical Monte Carlo 'restricted_auto' motif is not implemented yet");
    }

    /// \brief Generate supercell filling ConfigDoF from configuration
    ConfigDoF Canonical::_configname_motif(const std::string &configname) const {

      _log().set("DoF");
      _log() << "motif configname: " << configname << "\n";
      _log() << "using configation: " << configname << "\n" << std::endl;

      const Configuration &config = primclex().configuration(configname);
      const SymGroup &g = primclex().get_prim().factor_group();
      return config.fill_supercell(_supercell(), g).configdof();
    }

    /// \brief Find a grand canonical OccSwap to help enforce composition
    ///
    /// - Find the OccSwap that by applying minimizes:
    ///     (comp_n - conditions().mol_composition()).norm()
    /// - If cannot be improved, return end
    /// - If ties, randomly choose weighted by composition
    ///
    std::vector<OccSwap>::const_iterator
    Canonical::_find_grand_canonical_swap(
      const Configuration &config,
      std::vector<OccSwap>::const_iterator begin,
      std::vector<OccSwap>::const_iterator end) {

      double dn = 1. / supercell().volume();
      Eigen::VectorXd target_comp_n = conditions().mol_composition();
      Eigen::VectorXd comp_n = CASM::comp_n(config);

      typedef std::vector<OccSwap>::const_iterator it_type;
      std::vector<std::pair<it_type, double> > best;
      double best_dist = (comp_n - target_comp_n).norm();
      double tol = primclex().settings().lin_alg_tol();

      for(auto it = begin; it != end; ++it) {
        if(m_occ_loc.cand_size(it->cand_a)) {
          Eigen::VectorXd tcomp_n = comp_n;
          tcomp_n[it->cand_a.species_index] -= dn;
          tcomp_n[it->cand_b.species_index] += dn;

          double dist = (tcomp_n - target_comp_n).norm();
          if(dist < best_dist - tol) {
            best.clear();
            best.push_back({it, m_occ_loc.cand_size(it->cand_a)});
            best_dist = dist;
          }
          else if(dist < best_dist + tol) {
            best.push_back({it, m_occ_loc.cand_size(it->cand_a)});
          }
        }
      }

      if(!best.size()) {
        return end;
      }

      // break ties randomly, weighted by number of candidates
      double sum = 0.0;
      for(const auto &val : best) {
        sum += val.second;
      }

      double rand = _mtrand().randExc(sum);
      sum = 0.0;
      int count = 0;
      for(const auto &val : best) {
        sum += val.second;
        if(rand < sum) {
          return val.first;
        }
        ++count;
      }
      throw std::runtime_error("Error enforcing composition");
    }

    /// \brief Enforce composition by repeatedly applying GrandCanonicalSwap
    ///
    /// - Minimizes:
    ///     (comp_n - conditions().mol_composition()).norm()
    ///
    ConfigDoF Canonical::_enforce_conditions(const ConfigDoF &configdof) {

      _log().custom("Enforce composition");
      Configuration tconfig(_supercell(), jsonParser(), configdof);
      m_occ_loc.initialize(tconfig);
      jsonParser json;

      _log() << "    initial comp: " << to_json_array(CASM::comp(tconfig), json) << std::endl;
      _log() << "  initial comp_n: " << to_json_array(CASM::comp_n(tconfig), json) << std::endl;

      Eigen::VectorXd target_comp_n = conditions().mol_composition();
      _log() << "   target comp_n: " << to_json_array(target_comp_n, json) << std::endl;

      int count = 0;
      OccEvent e;
      ConfigDoF &tconfigdof = tconfig.configdof();
      while(true) {
        auto begin = m_cand.grand_canonical_swap().begin();
        auto end = m_cand.grand_canonical_swap().end();
        auto it = _find_grand_canonical_swap(tconfig, begin, end);

        if(it == end) {
          _log() << "   applied swaps: " << count << std::endl;
          break;
        }

        /// apply chosen swap (*it)
        m_occ_loc.propose_grand_canonical(e, *it, _mtrand());
        m_occ_loc.apply(e, tconfigdof);

        ++count;
      }

      _log() << "      final comp: " << to_json_array(CASM::comp(tconfig), json) << std::endl;
      _log() << "    final comp_n: " << to_json_array(CASM::comp_n(tconfig), json) << std::endl << std::endl;

      return tconfig.configdof();
    }
  }
}
