#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"

#include "casm/basis_set/DoF.hh"
#include "casm/clex/ConfigCorrelations.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/Norm.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/monte_carlo/MonteCarloEnum_impl.hh"
#include "casm/monte_carlo/MonteCarlo_impl.hh"
#include "casm/monte_carlo/MonteIO_impl.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonicalSettings_impl.hh"

namespace CASM {
namespace Monte {

template MonteCarlo::MonteCarlo(const PrimClex &primclex,
                                const GrandCanonicalSettings &settings,
                                Log &_log);
template MonteCarloEnum::MonteCarloEnum(const PrimClex &primclex,
                                        const GrandCanonicalSettings &settings,
                                        Log &log, GrandCanonical &mc);

const Monte::ENSEMBLE GrandCanonical::ensemble =
    Monte::ENSEMBLE::GrandCanonical;

// note: the construction process needs refactoring
Clex make_clex(PrimClex const &primclex,
               GrandCanonicalSettings const &settings) {
  ClexDescription const &clex_desc = settings.formation_energy(primclex);
  return Clex{primclex.clexulator(clex_desc.bset), primclex.eci(clex_desc)};
}
// note: the construction process needs refactoring
std::shared_ptr<RandomAlloyCorrCalculator> make_random_alloy_corr_f(
    PrimClex const &primclex, GrandCanonicalSettings const &settings) {
  auto shared_prim = primclex.shared_prim();
  ClexDescription const &clex_desc = settings.formation_energy(primclex);
  std::string basis_set_name = clex_desc.bset;

  std::vector<IntegralCluster> prototypes;
  jsonParser clust_json{primclex.dir().clust(basis_set_name)};
  read_clust(std::back_inserter(prototypes), clust_json, *shared_prim);

  return std::make_shared<RandomAlloyCorrCalculator>(
      shared_prim, primclex.basis_set_specs(basis_set_name), prototypes);
}

/// \brief Constructs a GrandCanonical object and prepares it for running based
/// on MonteSettings
///
/// - Does not set 'state': conditions or ConfigDoF
GrandCanonical::GrandCanonical(const PrimClex &primclex,
                               const GrandCanonicalSettings &settings, Log &log)
    : MonteCarlo(primclex, settings, log),
      m_site_swaps(supercell()),
      m_composition_converter(primclex.composition_axes()),
      m_formation_energy_clex(make_clex(primclex, settings)),
      m_order_parameter(settings.make_order_parameter(primclex)),
      m_random_alloy_corr_f(make_random_alloy_corr_f(primclex, settings)),
      m_convert(_supercell()),
      m_event(m_composition_converter.components().size(),
              _clexulator().corr_size()) {
  const auto &desc = settings.formation_energy(primclex);

  _log().construct("Grand Canonical Monte Carlo");
  _log() << "project: " << this->primclex().dir().root_dir() << "\n";
  _log() << "formation_energy cluster expansion: " << desc.name << "\n";
  _log() << std::setw(16) << "property: " << desc.property << "\n";
  _log() << std::setw(16) << "calctype: " << desc.calctype << "\n";
  _log() << std::setw(16) << "ref: " << desc.ref << "\n";
  _log() << std::setw(16) << "bset: " << desc.bset << "\n";
  _log() << std::setw(16) << "eci: " << desc.eci << "\n";
  if (m_order_parameter != nullptr) {
    jsonParser json;
    to_json(m_order_parameter->dof_space(), json);
    _log().custom("Order Parameter DoF Space");
    _log() << json << std::endl;
    m_order_parameter->update(config());
  }
  _log() << "supercell: \n" << supercell().transf_mat() << "\n";
  _log() << "\nSampling: \n";
  _log() << std::setw(24) << "quantity" << std::setw(24)
         << "requested_precision"
         << "\n";
  for (auto it = samplers().begin(); it != samplers().end(); ++it) {
    _log() << std::setw(24) << it->first;
    if (it->second->must_converge()) {
      _log() << std::setw(24) << it->second->requested_precision() << std::endl;
    } else {
      _log() << std::setw(24) << "none" << std::endl;
    }
  }
  _log() << "\nautomatic convergence mode?: " << std::boolalpha
         << must_converge() << std::endl;
  _log() << std::endl;
}

/// \brief Return number of steps per pass. Equals number of sites with variable
/// occupation.
Index GrandCanonical::steps_per_pass() const {
  return m_site_swaps.variable_sites().size();
}

/// \brief Return current conditions
const GrandCanonical::CondType &GrandCanonical::conditions() const {
  return m_condition;
}

/// \brief Set conditions and clear previously collected data
void GrandCanonical::set_conditions(
    const GrandCanonicalConditions &new_conditions) {
  _log().set("Conditions");
  _log() << new_conditions << std::endl << std::endl;

  m_condition = new_conditions;

  clear_samples();
  _update_properties();

  return;
}

/// \brief Set configdof and clear previously collected data
void GrandCanonical::set_configdof(const ConfigDoF &configdof,
                                   const std::string &msg) {
  _log().set("DoF");
  if (!msg.empty()) {
    _log() << msg << "\n";
  }
  _log() << std::endl;

  reset(configdof);
}

/// \brief Set configdof and conditions and clear previously collected data
///
/// \returns Specified ConfigDoF and configname (or configdof path)
///
std::pair<ConfigDoF, std::string> GrandCanonical::set_state(
    const GrandCanonicalConditions &new_conditions,
    const GrandCanonicalSettings &settings) {
  _log().set("Conditions");
  _log() << new_conditions << std::endl;

  m_condition = new_conditions;

  ConfigDoF configdof = _default_motif();

  std::string configname;

  if (settings.is_motif_configname()) {
    configname = settings.motif_configname();

    if (configname == "default") {
      // configdof = _default_motif();
    } else if (configname == "auto") {
      std::tie(configdof, configname) = _auto_motif(new_conditions);
    } else if (configname == "restricted_auto") {
      std::tie(configdof, configname) = _restricted_auto_motif(new_conditions);
    } else {
      configdof = _configname_motif(configname);
    }

  } else if (settings.is_motif_configdof()) {
    _log().set("DoF");
    _log() << "motif configdof: " << settings.motif_configdof_path() << "\n";
    configdof = settings.motif_configdof(supercell().volume());
    configname = settings.motif_configdof_path().string();
  } else if (settings.is_motif_config()) {
    _log().set("DoF");
    _log() << "motif config: \n" << settings.motif_config_json() << "\n";
    configdof = settings.motif_config(supercell());
    configname = "fromjson";
  } else {
    throw std::runtime_error(
        "Error: Must specify motif \"config\", \"configname\", or "
        "\"configdof\"");
  }

  reset(configdof);
  _update_properties();

  return std::make_pair(configdof, configname);
}

/// \brief Set configdof and conditions and clear previously collected data
void GrandCanonical::set_state(const GrandCanonicalConditions &new_conditions,
                               const ConfigDoF &configdof,
                               const std::string &msg) {
  _log().set("Conditions");
  _log() << new_conditions << std::endl << std::endl;

  m_condition = new_conditions;

  _log().set("DoF");
  if (!msg.empty()) {
    _log() << msg << "\n";
  }
  _log() << std::endl;

  reset(configdof);
  _update_properties();

  return;
}

/// \brief Propose a new event, calculate delta properties, and return reference
/// to it
///
/// Randomly picks a site that's allowed more than one occupant, and randomly
/// picks what occupant it changes to. Then calculates delta properties
/// associated with that change.
///
const GrandCanonical::EventType &GrandCanonical::propose() {
  // Randomly pick a site that's allowed more than one occupant
  Index random_variable_site =
      _mtrand().randInt(m_site_swaps.variable_sites().size() - 1);

  // Determine what that site's linear index is and what the sublattice index is
  Index mutating_site = m_site_swaps.variable_sites()[random_variable_site];
  Index sublat = m_site_swaps.sublattice()[random_variable_site];

  // Determine the current occupant of the mutating site
  int current_occupant = configdof().occ(mutating_site);

  // Randomly pick a new occupant for the mutating site
  const std::vector<int> &possible_mutation =
      m_site_swaps.possible_swap()[sublat][current_occupant];
  int new_occupant =
      possible_mutation[_mtrand().randInt(possible_mutation.size() - 1)];

  if (debug()) {
    const auto &site_occ = primclex().prim().basis()[sublat].occupant_dof();
    _log().custom("Propose event");

    _log() << "  Mutating site (linear index): " << mutating_site << "\n"
           << "  Mutating site (b, i, j, k): "
           << supercell().uccoord(mutating_site) << "\n"
           << "  Current occupant: " << current_occupant << " ("
           << site_occ[current_occupant].name() << ")\n"
           << "  Proposed occupant: " << new_occupant << " ("
           << site_occ[new_occupant].name() << ")\n\n"

           << "  beta: " << m_condition.beta() << "\n"
           << "  T: " << m_condition.temperature() << std::endl
           << std::endl;
  }

  // Update delta properties in m_event
  _update_deltas(m_event, mutating_site, sublat, current_occupant,
                 new_occupant);

  return m_event;
}

/// \brief Based on a random number, decide if the change in energy from the
/// proposed event is low enough to be accepted.
bool GrandCanonical::check(const GrandCanonicalEvent &event) {
  if (event.dEpot() < 0.0) {
    if (debug()) {
      _log().custom("Check event");
      _log() << "Probability to accept: 1.0\n" << std::endl;
    }
    return true;
  }

  double rand = _mtrand().rand53();
  double prob = exp(-event.dEpot() * m_condition.beta());

  if (debug()) {
    _log().custom("Check event");
    _log() << "Probability to accept: " << prob << "\n"
           << "Random number: " << rand << "\n"
           << std::endl;
  }

  return rand < prob;
}

/// \brief Accept proposed event. Change configuration accordingly and update
/// energies etc.
///
/// Once you've found an event that passes the random number test, you want to
/// accept it. This routine will apply all the changes in the event to *this
/// (change occupation of one atom) and update the formation energy, generalized
/// enthalpy, number of species and correlations values.
///
void GrandCanonical::accept(const EventType &event) {
  if (debug()) {
    _log().custom("Accept Event");
    _log() << std::endl;
  }

  // First apply changes to configuration (just a single occupant change)
  _configdof().occ(event.occupational_change().site_index()) =
      event.occupational_change().to_value();

  // Next update all properties that changed from the event
  _formation_energy() += event.dEf() / supercell().volume();
  _potential_energy() += event.dEpot() / supercell().volume();
  _corr() += event.dCorr() / supercell().volume();
  _comp_n() += event.dN().cast<double>() / supercell().volume();
  if (m_order_parameter != nullptr) {
    _eta() += event.deta();
  }
  return;
}

/// \brief Nothing needs to be done to reject a GrandCanonicalEvent
void GrandCanonical::reject(const EventType &event) {
  if (debug()) {
    _log().custom("Reject Event");
    _log() << std::endl;
  }
  return;
}

/// \brief Calculate the single spin flip low temperature expansion of the grand
/// canonical potential
///
/// Returns low temperature expansion estimate of the grand canonical free
/// energy. Works with the current ConfigDoF as groundstate.
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
/// For low temperatures we can approximate Z by truncating the sum after
/// microstates that only involve point defects and no defects
/// Z=boltz(\Omega_0)*SUM(boltz(D\Omega_s))  summing over all states with only
/// point defects or no defects
///
/// The free energy is
/// Phi=-kB*T*ln(Z)
/// Phi=-kB*T*(-\Omega_0/kBT+ln(SUM(boltz(D\Omega_s))    Sum is over point
/// defects and no defects (in which case D\Omega_s == 0)
/// Phi=(\Omega_0-kB*T(ln(SUM(boltz(D\Omega_s)))))/N
///
double GrandCanonical::lte_grand_canonical_free_energy() const {
  const SiteExchanger &site_exch = m_site_swaps;
  const ConfigDoF &config_dof = configdof();
  GrandCanonicalEvent event = m_event;

  double tol = 1e-12;

  auto less = [&](const double &A, const double &B) { return A < B - tol; };

  std::map<double, unsigned long, decltype(less)> hist(less);

  // no defect case
  hist[0.0] = 1;

  // Loop over sites that can change occupants
  for (Index exch_ind = 0; exch_ind < site_exch.variable_sites().size();
       exch_ind++) {
    // Transform exchanger index to ConfigDoF index
    Index mutating_site = site_exch.variable_sites()[exch_ind];
    int sublat = site_exch.sublattice()[exch_ind];
    int current_occupant = config_dof.occ(mutating_site);

    // Loop over possible occupants for site that can change
    const auto &possible = site_exch.possible_swap()[sublat][current_occupant];
    for (auto new_occ_it = possible.begin(); new_occ_it != possible.end();
         ++new_occ_it) {
      _update_deltas(event, mutating_site, sublat, current_occupant,
                     *new_occ_it);

      // save the result
      double dpot_nrg = event.dEpot();
      if (dpot_nrg < 0.0) {
        Log &err_log = CASM::err_log();
        err_log.error<Log::standard>("Calculating low temperature expansion");
        err_log << "  Defect lowered the potential energy. Your motif "
                   "configuration "
                << "is not the 0K ground state.\n"
                << std::endl;
        throw std::runtime_error(
            "Error calculating low temperature expansion. Not in the ground "
            "state.");
      }

      auto it = hist.find(dpot_nrg);
      if (it == hist.end()) {
        hist[dpot_nrg] = 1;
      } else {
        it->second++;
      }
    }
  }

  _log().results("Ground state and point defect potential energy details");
  _log() << "T: " << m_condition.temperature() << std::endl;
  _log() << "kT: " << 1.0 / m_condition.beta() << std::endl;
  _log() << "Beta: " << m_condition.beta() << std::endl << std::endl;

  _log() << std::setw(16) << "N/unitcell"
         << " " << std::setw(16) << "dPE"
         << " " << std::setw(24) << "N*exp(-dPE/kT)"
         << " " << std::setw(16) << "dphi"
         << " " << std::setw(16) << "phi" << std::endl;

  double tsum = 0.0;
  double phi = 0.0;
  double phi_prev;
  for (auto it = hist.rbegin(); it != hist.rend(); ++it) {
    phi_prev = phi;
    tsum += it->second * exp(-(it->first) * m_condition.beta());
    phi = std::log(tsum) / m_condition.beta() / supercell().volume();

    if (CASM::almost_equal(it->first, 0.0, tol)) {
      _log() << std::setw(16) << "(gs)"
             << " ";
    } else {
      _log() << std::setw(16) << std::setprecision(8)
             << (1.0 * it->second) / supercell().volume() << " ";
    }
    _log() << std::setw(16) << std::setprecision(8) << it->first << " "
           << std::setw(24) << std::setprecision(8)
           << it->second * exp(-it->first * m_condition.beta()) << " "
           << std::setw(16) << std::setprecision(8) << phi - phi_prev << " "
           << std::setw(16) << std::setprecision(8) << potential_energy() - phi
           << std::endl;
  }

  _log() << "phi_LTE(1): " << std::setprecision(12) << potential_energy() - phi
         << std::endl
         << std::endl;

  return potential_energy() - phi;
}

/// \brief Write results to files
void GrandCanonical::write_results(Index cond_index) const {
  CASM::Monte::write_results(settings(), *this, _log());
  write_conditions_json(settings(), *this, cond_index, _log());
  write_observations(settings(), *this, cond_index, _log());
  write_trajectory(settings(), *this, cond_index, _log());
  // write_pos_trajectory(settings(), *this, cond_index);
}

/// \brief Get potential energy, normalized per primitive cell
///
/// - if(&config == &this->config()) { return potential_energy(); }, else
///   calculate potential_energy = formation_energy - comp_x.dot(param_chem_pot)
///
/// \throws If config does not have the same supercell as this
double GrandCanonical::potential_energy(const Configuration &config) const {
  if (debug()) {
    _log().calculate("Potential energy");
  }
  // if(&config == &this->config()) { return potential_energy(); }
  auto corr = correlations(config, _clexulator());
  auto comp_n = CASM::comp_n(config);
  auto comp_x = m_composition_converter.param_composition(comp_n);

  if (debug()) {
    _print_correlations(corr, "correlations", "corr");
    _log() << "  comp_n: " << comp_n.transpose() << std::endl;
    _log() << "  comp_x: " << comp_x.transpose() << std::endl;
  }

  double potential_energy = 0;

  if (m_condition.include_formation_energy()) {
    double formation_energy = _eci() * corr.data();
    potential_energy += formation_energy;
    if (debug()) {
      _log() << "  formation_energy: " << formation_energy << std::endl;
      _log() << "  potential_energy (w/ formation_energy): " << potential_energy
             << std::endl;
    }
  }
  if (m_condition.include_param_chem_pot()) {
    potential_energy -= m_condition.param_chem_pot().dot(comp_x);
    if (debug()) {
      _log() << "  potential_energy (w/ param_chem_pot): " << potential_energy
             << std::endl;
    }
  }
  if (m_condition.corr_matching_pot().has_value()) {
    potential_energy +=
        corr_matching_potential(corr, *m_condition.corr_matching_pot());
    if (debug()) {
      _log() << "  potential_energy (w/ corr_matching_pot): "
             << potential_energy << std::endl;
    }
  }
  if (m_condition.random_alloy_corr_matching_pot().has_value()) {
    potential_energy += corr_matching_potential(
        corr, *m_condition.random_alloy_corr_matching_pot());
    if (debug()) {
      _log() << "  potential_energy (w/ random_alloy_corr_matching_pot): "
             << potential_energy << std::endl;
    }
  }
  if (m_condition.param_comp_quad_pot_vector().has_value()) {
    Eigen::VectorXd x = (comp_x - *m_condition.param_comp_quad_pot_target());
    Eigen::VectorXd const &v = *m_condition.param_comp_quad_pot_vector();
    potential_energy += v.dot((x.array() * x.array()).matrix());
    if (debug()) {
      _log() << "  potential_energy (w/ param_comp_quad_pot_vector): "
             << potential_energy << std::endl;
    }
  }
  if (m_condition.param_comp_quad_pot_matrix().has_value()) {
    Eigen::VectorXd x = (comp_x - *m_condition.param_comp_quad_pot_target());
    Eigen::VectorXd const &V = *m_condition.param_comp_quad_pot_matrix();
    potential_energy += x.dot(V * x);
    if (debug()) {
      _log() << "  potential_energy (w/ param_comp_quad_pot_matrix): "
             << potential_energy << std::endl;
    }
  }
  if (m_order_parameter != nullptr) {
    // ---
    // this is a hack to allow GrandCanonical::potential to be
    // used on any Configuration
    Eigen::VectorXd eta;
    if (&this->supercell() != &config.supercell()) {
      OrderParameter tmp = *m_order_parameter;
      eta = tmp(config);
    } else {
      // temporarily reset ConfigDoFValues pointer
      // - requires that the supercell is unchanged
      clexulator::ConfigDoFValues const *prev = m_order_parameter->get();
      m_order_parameter->set(&config.configdof().values());
      eta = m_order_parameter->value();
      m_order_parameter->set(prev);
    }
    // ---
    if (debug()) {
      _log() << "  eta: " << eta.transpose() << std::endl;
    }

    if (m_condition.order_parameter_pot().has_value()) {
      Eigen::VectorXd const &v = *m_condition.order_parameter_pot();
      potential_energy -= v.dot(eta);
      if (debug()) {
        _log() << "  potential_energy (w/ order_parameter_pot): "
               << potential_energy << std::endl;
      }
    }
    if (m_condition.order_parameter_quad_pot_vector().has_value()) {
      Eigen::VectorXd x =
          (eta - *m_condition.order_parameter_quad_pot_target());
      Eigen::VectorXd const &v = *m_condition.order_parameter_quad_pot_vector();
      potential_energy += v.dot((x.array() * x.array()).matrix());
      if (debug()) {
        _log() << "  potential_energy (w/ order_parameter_quad_pot_vector): "
               << potential_energy << std::endl;
      }
    }
    if (m_condition.order_parameter_quad_pot_matrix().has_value()) {
      Eigen::VectorXd x =
          (eta - *m_condition.order_parameter_quad_pot_target());
      Eigen::VectorXd const &V = *m_condition.order_parameter_quad_pot_matrix();
      potential_energy += x.dot(V * x);
      if (debug()) {
        _log() << "  potential_energy (w/ order_parameter_quad_pot_matrix): "
               << potential_energy << std::endl;
      }
    }
  }

  if (debug()) {
    _log() << std::endl;
  }
  return potential_energy;
}

/// \brief Calculate delta correlations for an event
void GrandCanonical::_set_dCorr(Index mutating_site, int new_occupant,
                                Eigen::VectorXd &dCorr) const {
  restricted_delta_corr(dCorr, mutating_site, new_occupant, configdof(),
                        supercell().nlist(), _clexulator(),
                        _eci().index().data(), end_ptr(_eci().index()));

  if (debug()) {
    _print_correlations(dCorr, "delta correlations", "dCorr");
  }
}

/// \brief Print correlations to _log()
void GrandCanonical::_print_correlations(const Eigen::VectorXd &corr,
                                         std::string title,
                                         std::string colheader) const {
  _log().calculate(title);
  _log() << std::setw(12) << "i" << std::setw(16) << "ECI" << std::setw(16)
         << colheader << std::endl;

  for (int i = 0; i < corr.size(); ++i) {
    double eci = 0.0;
    bool calculated = true;
    Index index = find_index(_eci().index(), i);
    if (index != _eci().index().size()) {
      eci = _eci().value()[index];
    }
    if (index == _eci().index().size()) {
      calculated = false;
    }

    _log() << std::setw(12) << i << std::setw(16) << std::setprecision(8)
           << eci;
    if (calculated) {
      _log() << std::setw(16) << std::setprecision(8) << corr[i];
    } else {
      _log() << std::setw(16) << "unknown";
    }
    _log() << std::endl;
  }
  _log() << std::endl;
}

/// \brief Update delta properties in 'event'
void GrandCanonical::_update_deltas(GrandCanonicalEvent &event,
                                    Index mutating_site, int sublat,
                                    int current_occupant,
                                    int new_occupant) const {
  if (debug()) {
    _log().custom("Update deltas");
  }
  // ---- set OccMod --------------
  event.occupational_change().set(mutating_site, sublat, new_occupant);

  // ---- set dspecies --------------
  for (int i = 0; i < event.dN().size(); ++i) {
    event.set_dN(i, 0);
  }
  Index curr_species = m_site_swaps.sublat_to_mol()[sublat][current_occupant];
  Index new_species = m_site_swaps.sublat_to_mol()[sublat][new_occupant];
  event.set_dN(curr_species, -1);
  event.set_dN(new_species, 1);

  // ---- set dcorr --------------
  _set_dCorr(mutating_site, new_occupant, event.dCorr());

  // ---- set dformation_energy --------------
  event.set_dEf(_eci() * event.dCorr().data());

  // ---- set deta -------------
  if (m_order_parameter != nullptr) {
    event.deta() = m_order_parameter->occ_delta(mutating_site, new_occupant);
  }

  if (debug()) {
    _log() << "  dEf: " << event.dEf() << std::endl;
    // note: `dcomp_x` is extensive (includes Nunit),
    //       but `comp_x` is intensive
    double Nunit = supercell().volume();
    Eigen::VectorXd dcomp_x =
        m_composition_converter.dparam_dmol() * event.dN().cast<double>();
    Eigen::VectorXd comp_x =
        m_composition_converter.param_composition(this->comp_n());
    _log() << "  dN: " << event.dN().transpose() << std::endl;
    _log() << "  comp_n: " << this->comp_n().transpose() << std::endl;
    _log() << "  dcomp_n: "
           << (event.dN().cast<double>() / supercell().volume()).transpose()
           << std::endl;
    _log() << "  comp_x: " << comp_x.transpose() << std::endl;
    _log() << "  dcomp_x: " << (dcomp_x / Nunit).transpose() << std::endl;
    if (m_order_parameter != nullptr) {
      _log() << "  eta: " << this->eta().transpose() << std::endl;
      _log() << "  deta: " << event.deta().transpose() << std::endl;
    }
  }

  // ---- set dpotential_energy --------------
  double dEpot = 0;
  if (m_condition.include_formation_energy()) {
    dEpot += event.dEf();
    if (debug()) {
      _log() << "  dEpot (w/ formation_energy): " << dEpot << std::endl;
    }
  }
  if (m_condition.include_param_chem_pot()) {
    dEpot -= m_condition.exchange_chem_pot(new_species, curr_species);
    if (debug()) {
      _log() << "  dEpot (w/ param_chem_pot): " << dEpot << std::endl;
    }
  }
  if (m_condition.corr_matching_pot()) {
    dEpot += delta_corr_matching_potential(this->corr(),
                                           event.dCorr() / supercell().volume(),
                                           *m_condition.corr_matching_pot());
    if (debug()) {
      _log() << "  dEpot (w/ corr_matching_pot): " << dEpot << std::endl;
    }
  }
  if (m_condition.random_alloy_corr_matching_pot()) {
    dEpot += delta_corr_matching_potential(
        this->corr(), event.dCorr() / supercell().volume(),
        *m_condition.random_alloy_corr_matching_pot());
    if (debug()) {
      _log() << "  dEpot (w/ random_alloy_corr_matching_pot): " << dEpot
             << std::endl;
    }
  }
  if (m_condition.param_comp_quad_pot_vector().has_value()) {
    // note: `dcomp_x` is extensive (includes Nunit), but `comp_x` and `target`
    // are intensive
    double Nunit = supercell().volume();
    Eigen::VectorXd dcomp_x =
        m_composition_converter.dparam_dmol() * event.dN().cast<double>();
    Eigen::VectorXd comp_x =
        m_composition_converter.param_composition(this->comp_n());
    Eigen::VectorXd const &target = *m_condition.param_comp_quad_pot_target();
    Eigen::VectorXd const &v = *m_condition.param_comp_quad_pot_vector();

    for (int i = 0; i < dcomp_x.size(); ++i) {
      dEpot += v(i) * (2.0 * (comp_x(i) - target(i)) * dcomp_x(i) +
                       dcomp_x(i) * (dcomp_x(i) / Nunit));
    }
    if (debug()) {
      auto x2 = (comp_x + dcomp_x / Nunit - target);
      auto x1 = (comp_x - target);
      auto expected = Nunit * (v.dot((x2.array() * x2.array()).matrix()) -
                               v.dot((x1.array() * x1.array()).matrix()));

      _log() << "  dEpot (w/ param_comp_quad_pot_vector): " << dEpot << "  ("
             << expected << ")" << std::endl;
    }
  }
  if (m_condition.param_comp_quad_pot_matrix().has_value()) {
    // note: `dcomp_x` is extensive (includes Nunit), but `comp_x` and `target`
    // are intensive
    double Nunit = supercell().volume();
    Eigen::VectorXd dcomp_x =
        m_composition_converter.dparam_dmol() * event.dN().cast<double>();
    Eigen::VectorXd comp_x =
        m_composition_converter.param_composition(this->comp_n());
    Eigen::VectorXd const &target = *m_condition.param_comp_quad_pot_target();
    Eigen::MatrixXd const &V = *m_condition.param_comp_quad_pot_matrix();

    for (int i = 0; i < dcomp_x.size(); ++i) {
      for (int j = 0; j < dcomp_x.size(); ++j) {
        dEpot +=
            Nunit * V(i, j) *
            (comp_x(i) * dcomp_x(j) + comp_x(j) * dcomp_x(i) -
             2.0 * target(i) * dcomp_x(j) + dcomp_x(i) * (dcomp_x(j) / Nunit));
      }
    }
    if (debug()) {
      auto x2 = (comp_x + dcomp_x / Nunit - target);
      auto x1 = (comp_x - target);
      auto expected = Nunit * (x2.dot(V * x2) - x1.dot(V * x1));

      _log() << "  dEpot (w/ param_comp_quad_pot_matrix): " << dEpot << "  ("
             << expected << ")" << std::endl;
    }
  }
  if (m_order_parameter != nullptr) {
    // note: eta and deta are intensive
    double Nunit = supercell().volume();
    Eigen::VectorXd const &eta = this->eta();
    Eigen::VectorXd const &deta = event.deta();
    if (m_condition.order_parameter_pot().has_value()) {
      Eigen::VectorXd const &v = *m_condition.order_parameter_pot();
      dEpot -= Nunit * v.dot(deta);
      if (debug()) {
        auto expected = -Nunit * v.dot(deta);

        _log() << "  dEpot (w/ order_parameter_pot): " << dEpot << "  ("
               << expected << ")" << std::endl;
      }
    }
    if (m_condition.order_parameter_quad_pot_vector().has_value()) {
      Eigen::VectorXd const &target =
          *m_condition.order_parameter_quad_pot_target();
      Eigen::VectorXd const &v = *m_condition.order_parameter_quad_pot_vector();
      for (int i = 0; i < eta.size(); ++i) {
        dEpot += Nunit * v(i) *
                 (2.0 * (eta(i) - target(i)) * deta(i) + deta(i) * deta(i));
      }
      if (debug()) {
        auto x2 = (eta + deta - target);
        auto x1 = (eta - target);
        auto expected = Nunit * (v.dot((x2.array() * x2.array()).matrix()) -
                                 v.dot((x1.array() * x1.array()).matrix()));

        _log() << "  dEpot (w/ order_parameter_quad_pot_vector): " << dEpot
               << "  (" << expected << ")" << std::endl;
      }
    }
    if (m_condition.order_parameter_quad_pot_matrix().has_value()) {
      Eigen::VectorXd const &target =
          *m_condition.order_parameter_quad_pot_target();
      Eigen::MatrixXd const &V = *m_condition.order_parameter_quad_pot_matrix();

      for (int i = 0; i < deta.size(); ++i) {
        for (int j = 0; j < deta.size(); ++j) {
          dEpot += Nunit * V(i, j) *
                   (eta(i) * deta(j) + eta(j) * deta(i) -
                    2.0 * target(i) * deta(j) + deta(i) * deta(j));
        }
      }
      if (debug()) {
        auto x2 = (eta + deta - target);
        auto x1 = (eta - target);
        auto expected = Nunit * (x2.dot(V * x2) - x1.dot(V * x1));

        _log() << "  dEpot (w/ order_parameter_quad_pot_matrix): " << dEpot
               << "  (" << expected << ")" << std::endl;
      }
    }
  }

  if (debug()) {
    _log() << std::endl;
  }
  event.set_dEpot(dEpot);
}

/// \brief Calculate properties given current conditions
void GrandCanonical::_update_properties() {
  // initialize properties and store pointers to the data strucures
  _vector_properties()["corr"] = correlations(_config(), _clexulator());
  m_corr = &_vector_property("corr");

  _vector_properties()["comp_n"] = CASM::comp_n(_configdof(), supercell());
  m_comp_n = &_vector_property("comp_n");

  if (m_order_parameter != nullptr) {
    _vector_properties()["eta"] = m_order_parameter->value();
    m_eta = &_vector_property("eta");
  } else {
    _vector_properties()["eta"] = Eigen::VectorXd(0);
    m_eta = &_vector_property("eta");
  }

  _scalar_properties()["formation_energy"] = _eci() * corr().data();
  m_formation_energy = &_scalar_property("formation_energy");

  _scalar_properties()["potential_energy"] = this->potential_energy(config());
  m_potential_energy = &_scalar_property("potential_energy");
}

/// \brief Generate supercell filling ConfigDoF from default configuration
ConfigDoF GrandCanonical::_default_motif() const {
  _log().set("DoF");
  _log() << "motif configname: default\n";
  _log() << "using configuration with default occupation...\n" << std::endl;

  return Configuration::zeros(_supercell()).configdof();
}

/// \brief Generate minimum potential energy ConfigDoF
///
/// Raises exception if it doesn't tile the supercell
std::pair<ConfigDoF, std::string> GrandCanonical::_auto_motif(
    const GrandCanonicalConditions &cond) const {
  _log().set("DoF");
  _log() << "motif configname: auto\n";
  _log() << "searching for minimum potential energy motif..." << std::endl;

  double tol = 1e-6;
  auto compare = [&](double A, double B) { return A < B - tol; };

  _log() << "using conditions: \n";
  _log() << cond << std::endl;

  std::multimap<double, std::string, decltype(compare)> configmap(compare);
  const auto &db = primclex().db<Configuration>();
  for (const auto &tconfig : db) {
    configmap.insert(
        std::make_pair(_eci() * correlations(tconfig, _clexulator()) -
                           cond.param_chem_pot().dot(CASM::comp(tconfig)),
                       tconfig.name()));
  }

  Configuration min_config = *db.find(configmap.begin()->second);
  double min_potential_energy = configmap.begin()->first;
  auto eq_range = configmap.equal_range(min_potential_energy);

  if (std::distance(eq_range.first, eq_range.second) > 1) {
    _log() << "Warning: Found degenerate ground states with potential energy: "
           << std::setprecision(8) << min_potential_energy << std::endl;
    for (auto it = eq_range.first; it != eq_range.second; ++it) {
      _log() << "  " << db.find(it->second)->name() << std::endl;
    }
    _log() << "using: " << min_config.name() << "\n" << std::endl;
  } else {
    _log() << "using: " << min_config.name()
           << " with potential energy: " << std::setprecision(8)
           << min_potential_energy << "\n"
           << std::endl;
  }

  return std::make_pair(fill_supercell(min_config, _supercell()).configdof(),
                        min_config.name());
}

/// \brief Generate minimum potential energy ConfigDoF for this supercell
std::pair<ConfigDoF, std::string> GrandCanonical::_restricted_auto_motif(
    const GrandCanonicalConditions &cond) const {
  _log().set("DoF");
  _log() << "motif configname: restricted_auto\n";
  _log() << "searching for minimum potential energy motif..." << std::endl;

  double tol = 1e-6;
  auto compare = [&](double A, double B) { return A < B - tol; };

  _log() << "using conditions: \n";
  _log() << cond << std::endl;

  std::multimap<double, std::string, decltype(compare)> configmap(compare);
  const auto &db = primclex().db<Configuration>();
  for (const auto &tconfig : db) {
    configmap.insert(
        std::make_pair(_eci() * correlations(tconfig, _clexulator()) -
                           cond.param_chem_pot().dot(CASM::comp(tconfig)),
                       tconfig.name()));
  }

  // used to check if configurations can fill the monte carlo supercell
  const Lattice &scel_lat = supercell().lattice();
  const SymGroup &g = primclex().prim().factor_group();
  auto begin = g.begin();
  auto end = g.end();

  // save iterators pointing to configs that will fill the supercell
  std::vector<decltype(configmap)::const_iterator> allowed;

  // iterate through groups of degenerate configs
  auto next_it = configmap.begin();
  while (next_it != configmap.end()) {
    auto eq_range = configmap.equal_range(next_it->first);

    // save allowed configs
    for (auto it = eq_range.first; it != eq_range.second; ++it) {
      const Lattice &motif_lat = db.find(it->second)->supercell().lattice();
      if (xtal::is_equivalent_superlattice(scel_lat, motif_lat, begin, end, TOL)
              .first != end) {
        allowed.push_back(it);
      }
    }

    // if some found, break
    if (allowed.size()) {
      break;
    }

    // else, continue to next group
    next_it = eq_range.second;
  }

  if (!allowed.size()) {
    _log()
        << "Found no enumerated configurations that will fill the supercell\n";
    _log() << "using configuration with default occupation..." << std::endl;
    return std::make_pair(Configuration::zeros(_supercell()).configdof(),
                          "default");
  }

  if (allowed.size() > 1) {
    _log() << "Warning: Found degenerate allowed configurations with potential "
              "energy: "
           << std::setprecision(8) << allowed[0]->first << std::endl;
    for (auto it = allowed.begin(); it != allowed.end(); ++it) {
      _log() << "  " << (*it)->second << std::endl;
    }
    _log() << "using: " << allowed[0]->second << "\n" << std::endl;
  } else {
    _log() << "using: " << allowed[0]->second
           << " with potential energy: " << std::setprecision(8)
           << allowed[0]->first << "\n"
           << std::endl;
  }

  Configuration min_config = *db.find(allowed[0]->second);
  return std::make_pair(fill_supercell(min_config, _supercell()).configdof(),
                        min_config.name());
}

/// \brief Generate supercell filling ConfigDoF from configuration
ConfigDoF GrandCanonical::_configname_motif(
    const std::string &configname) const {
  _log().set("DoF");
  _log() << "motif configname: " << configname << "\n";
  _log() << "using configation: " << configname << "\n" << std::endl;

  auto it = primclex().db<Configuration>().find(configname);
  if (it == primclex().db<Configuration>().end()) {
    throw std::runtime_error(std::string("Error: did not find motif '") +
                             configname + "'");
  }
  Configuration config = *it;
  return fill_supercell(config, _supercell()).configdof();
}

}  // namespace Monte
}  // namespace CASM
