#include "casm/monte_carlo/canonical/Canonical.hh"

#include "casm/clex/ConfigCorrelations.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clusterography/io/OrbitPrinter_impl.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/algorithm.hh"
#include "casm/monte_carlo/MonteCarloEnum_impl.hh"
#include "casm/monte_carlo/MonteCarlo_impl.hh"
#include "casm/monte_carlo/MonteCorrelations.hh"
#include "casm/monte_carlo/MonteIO_impl.hh"
#include "casm/monte_carlo/canonical/CanonicalIO.hh"
#include "casm/monte_carlo/canonical/CanonicalSettings_impl.hh"

namespace CASM {
namespace Monte {

template MonteCarlo::MonteCarlo(const PrimClex &primclex,
                                const CanonicalSettings &settings, Log &_log);
template MonteCarloEnum::MonteCarloEnum(const PrimClex &primclex,
                                        const CanonicalSettings &settings,
                                        Log &log, Canonical &mc);

const ENSEMBLE Canonical::ensemble = ENSEMBLE::Canonical;

// note: the construction process needs refactoring
Clex make_clex(PrimClex const &primclex, CanonicalSettings const &settings) {
  ClexDescription const &clex_desc = settings.formation_energy(primclex);
  return Clex{primclex.clexulator(clex_desc.bset), primclex.eci(clex_desc)};
}
// note: the construction process needs refactoring
std::shared_ptr<RandomAlloyCorrCalculator> make_random_alloy_corr_f(
    PrimClex const &primclex, CanonicalSettings const &settings) {
  auto shared_prim = primclex.shared_prim();
  ClexDescription const &clex_desc = settings.formation_energy(primclex);
  std::string basis_set_name = clex_desc.bset;

  std::vector<IntegralCluster> prototypes;
  jsonParser clust_json{primclex.dir().clust(basis_set_name)};
  read_clust(std::back_inserter(prototypes), clust_json, *shared_prim);

  return std::make_shared<RandomAlloyCorrCalculator>(
      shared_prim, primclex.basis_set_specs(basis_set_name), prototypes);
}

/// \brief Constructs a Canonical object and prepares it for running based on
/// MonteSettings
///
/// - Does not set 'state': conditions or ConfigDoF
Canonical::Canonical(const PrimClex &primclex,
                     const CanonicalSettings &settings, Log &log)
    : MonteCarlo(primclex, settings, log),
      m_formation_energy_clex(make_clex(primclex, settings)),
      m_order_parameter(settings.make_order_parameter(primclex)),
      m_order_parameter_subspaces(settings.make_order_parameter_subspaces()),
      m_random_alloy_corr_f(make_random_alloy_corr_f(primclex, settings)),
      m_convert(_supercell()),
      m_cand(m_convert),
      m_occ_loc(m_convert, m_cand),
      m_event(primclex.composition_axes().components().size(),
              _clexulator().corr_size()) {
  const auto &desc = settings.formation_energy(primclex);

  _log().construct("Canonical Monte Carlo");
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
    _log().custom("DoF Space");
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

  _log() << std::pair<const OccCandidateList &, const Conversions &>(m_cand,
                                                                     m_convert)
         << std::endl;
}

/// \brief Return number of steps per pass. Equals number of sites with variable
/// occupation.
Index Canonical::steps_per_pass() const { return m_occ_loc.size(); }

/// \brief Return current conditions
const Canonical::CondType &Canonical::conditions() const { return m_condition; }

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
void Canonical::set_configdof(const ConfigDoF &configdof,
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
std::pair<ConfigDoF, std::string> Canonical::set_state(
    const CanonicalConditions &new_conditions,
    const CanonicalSettings &settings) {
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

  reset(_enforce_conditions(configdof));
  _update_properties();

  return std::make_pair(configdof, configname);
}

/// \brief Set configdof and conditions and clear previously collected data
void Canonical::set_state(const CanonicalConditions &new_conditions,
                          const ConfigDoF &configdof, const std::string &msg) {
  _log().set("Conditions");
  _log() << new_conditions << std::endl << std::endl;

  m_condition = new_conditions;

  _log().set("DoF");
  if (!msg.empty()) {
    _log() << msg << "\n";
  }
  _log() << std::endl;

  reset(_enforce_conditions(configdof));
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
const Canonical::EventType &Canonical::propose() {
  m_occ_loc.propose_canonical(m_event.occ_event(), m_cand.canonical_swap(),
                              _mtrand());

  if (debug()) {
    _log().custom("Propose event");

    for (int i = 0; i < m_event.occ_event().linear_site_index.size(); ++i) {
      Index l = m_event.occ_event().linear_site_index[i];
      Index asym = m_convert.l_to_asym(l);
      Index current_occupant = _configdof().occ(l);
      Index current_species = m_convert.species_index(asym, current_occupant);
      std::string current_species_name =
          m_convert.species_name(current_species);
      Index new_occupant = m_event.occ_event().new_occ[i];
      Index new_species = m_convert.species_index(asym, new_occupant);
      std::string new_species_name = m_convert.species_name(new_species);
      _log() << "- Mutating site (linear index): " << l << "\n"
             << "  Mutating site (b, i, j, k): " << m_convert.l_to_bijk(l)
             << "\n"
             << "  Current occupant: " << current_occupant << " ("
             << current_species_name << ")\n"
             << "  Proposed occupant: " << new_occupant << " ("
             << new_species_name << ")\n";
    }
    _log() << "\n";
    _log() << "  beta: " << m_condition.beta() << "\n"
           << "  T: " << m_condition.temperature() << std::endl
           << std::endl;
  }

  _update_deltas(m_event);
  return m_event;
}

/// \brief Based on a random number, decide if the change in energy from the
/// proposed event is low enough to be accepted.
bool Canonical::check(const CanonicalEvent &event) {
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
void Canonical::accept(const EventType &event) {
  if (debug()) {
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
  if (m_order_parameter != nullptr) {
    _eta() += event.deta();

    if (m_order_parameter_subspaces != nullptr) {
      auto const &subspaces = *m_order_parameter_subspaces;
      for (int i = 0; i < subspaces.size(); ++i) {
        double x = 0.0;
        for (int const &j : subspaces[i]) {
          x += _eta()(j) * _eta()(j);
        }
        _eta_subspace()(i) = sqrt(x);
      }
    }
  }
  return;
}

/// \brief Nothing needs to be done to reject a CanonicalEvent
void Canonical::reject(const EventType &event) {
  if (debug()) {
    _log().custom("Reject Event");
    _log() << std::endl;
  }
  return;
}

/// \brief Write results to files
void Canonical::write_results(Index cond_index) const {
  CASM::Monte::write_results(settings(), *this, _log());
  write_conditions_json(settings(), *this, cond_index, _log());
  write_observations(settings(), *this, cond_index, _log());
  write_trajectory(settings(), *this, cond_index, _log());
  // write_pos_trajectory(settings(), *this, cond_index);
}

/// \brief Get potential energy
///
/// - if(&config == &this->config()) { return potential_energy(); }, else
///   calculate potential_energy = formation_energy
///
/// \throws If config does not have the same supercell as this
double Canonical::potential_energy(const Configuration &config) const {
  if (debug()) {
    _log().calculate("Potential energy");
  }
  // if(&config == &this->config()) { return potential_energy(); }

  auto corr = correlations(config, _clexulator());

  if (debug()) {
    _print_correlations(corr, "correlations", "corr");
  }

  double potential_energy = 0;

  if (m_condition.include_formation_energy()) {
    potential_energy += _eci() * corr.data();
    if (debug()) {
      _log() << "  potential_energy (w/ formation_energy): " << potential_energy
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
  if (m_order_parameter != nullptr) {
    // ---
    // this is a hack to allow Canonical::potential to be
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
      Eigen::VectorXd const &v = *m_condition.order_parameter_quad_pot_matrix();
      potential_energy += x.dot(v * x);
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
void Canonical::_set_dCorr(CanonicalEvent &event) const {
  restricted_delta_corr(event.dCorr(), event.occ_event(), m_convert,
                        configdof(), supercell().nlist(), _clexulator(),
                        _eci().index().data(), end_ptr(_eci().index()));

  if (debug()) {
    _print_correlations(event.dCorr(), "delta correlations", "dCorr");
  }
}

/// \brief Print correlations to _log()
void Canonical::_print_correlations(const Eigen::VectorXd &corr,
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
void Canonical::_update_deltas(CanonicalEvent &event) const {
  if (debug()) {
    _log().custom("Update deltas");
  }
  // ---- set dCorr (extensive) --------------
  _set_dCorr(event);

  // ---- set dformation_energy (extensive) --------------
  event.set_dEf(_eci() * event.dCorr().data());

  // ---- set deta (intensive) -------------
  if (m_order_parameter != nullptr) {
    event.deta() = m_order_parameter->occ_delta(
        event.occ_event().linear_site_index, event.occ_event().new_occ);
  }

  if (debug()) {
    _log() << "  dEf: " << event.dEf() << std::endl;
    if (m_order_parameter != nullptr) {
      _log() << "  eta: " << this->eta().transpose() << std::endl;
      _log() << "  deta: " << event.deta().transpose() << std::endl;
    }
  }

  // ---- set dpotential_energy (extensive) --------------
  double dEpot = 0;
  if (m_condition.include_formation_energy()) {
    dEpot += event.dEf();
    if (debug()) {
      _log() << "  dEpot (w/ formation_energy): " << dEpot << std::endl;
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
void Canonical::_update_properties() {
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

  if (m_order_parameter_subspaces != nullptr) {
    auto const &subspaces = *m_order_parameter_subspaces;
    Eigen::VectorXd eta_subspace = Eigen::VectorXd::Zero(subspaces.size());
    for (int i = 0; i < subspaces.size(); ++i) {
      double x = 0.0;
      for (int const &j : subspaces[i]) {
        if (j < 0 || j >= _eta().size()) {
          throw std::runtime_error("Invalid order_parameter_subspaces");
        }
        x += _eta()(j) * _eta()(j);
      }
      eta_subspace(i) = sqrt(x);
    }
    _vector_properties()["eta_subspace"] = eta_subspace;
    m_eta_subspace = &_vector_property("eta_subspace");
  } else {
    _vector_properties()["eta_subspace"] = Eigen::VectorXd(0);
    m_eta_subspace = &_vector_property("eta_subspace");
  }

  _scalar_properties()["formation_energy"] = _eci() * corr().data();
  m_formation_energy = &_scalar_property("formation_energy");

  _scalar_properties()["potential_energy"] = this->potential_energy(config());
  m_potential_energy = &_scalar_property("potential_energy");
}

/// \brief Generate supercell filling ConfigDoF from default configuration
ConfigDoF Canonical::_default_motif() const {
  _log().set("DoF");
  _log() << "motif configname: default\n";
  _log() << "using configuration with default occupation...\n" << std::endl;

  return Configuration::zeros(_supercell()).configdof();
}

/// \brief Generate minimum potential energy ConfigDoF
///
/// Raises exception if it doesn't tile the supercell
std::pair<ConfigDoF, std::string> Canonical::_auto_motif(
    const CanonicalConditions &cond) const {
  throw std::runtime_error(
      "Canonical Monte Carlo 'auto' motif is not implemented yet");
}

/// \brief Generate minimum potential energy ConfigDoF for this supercell
std::pair<ConfigDoF, std::string> Canonical::_restricted_auto_motif(
    const CanonicalConditions &cond) const {
  throw std::runtime_error(
      "Canonical Monte Carlo 'restricted_auto' motif is not implemented yet");
}

/// \brief Generate supercell filling ConfigDoF from configuration
ConfigDoF Canonical::_configname_motif(const std::string &configname) const {
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

/// \brief Find a grand canonical OccSwap to help enforce composition
///
/// - Find the OccSwap that by applying minimizes:
///     (comp_n - conditions().mol_composition()).norm()
/// - If cannot be improved, return end
/// - If ties, randomly choose weighted by composition
///
std::vector<OccSwap>::const_iterator Canonical::_find_grand_canonical_swap(
    const Configuration &config, std::vector<OccSwap>::const_iterator begin,
    std::vector<OccSwap>::const_iterator end) {
  double dn = 1. / supercell().volume();
  Eigen::VectorXd target_comp_n = conditions().mol_composition();
  Eigen::VectorXd comp_n = CASM::comp_n(config);

  typedef std::vector<OccSwap>::const_iterator it_type;
  std::vector<std::pair<it_type, double> > best;
  double best_dist = (comp_n - target_comp_n).norm();
  double tol = primclex().settings().lin_alg_tol();

  for (auto it = begin; it != end; ++it) {
    if (m_occ_loc.cand_size(it->cand_a)) {
      Eigen::VectorXd tcomp_n = comp_n;
      tcomp_n[it->cand_a.species_index] -= dn;
      tcomp_n[it->cand_b.species_index] += dn;

      double dist = (tcomp_n - target_comp_n).norm();
      if (dist < best_dist - tol) {
        best.clear();
        best.push_back({it, m_occ_loc.cand_size(it->cand_a)});
        best_dist = dist;
      } else if (dist < best_dist + tol) {
        best.push_back({it, m_occ_loc.cand_size(it->cand_a)});
      }
    }
  }

  if (!best.size()) {
    return end;
  }

  // break ties randomly, weighted by number of candidates
  double sum = 0.0;
  for (const auto &val : best) {
    sum += val.second;
  }

  double rand = _mtrand().randExc(sum);
  sum = 0.0;
  int count = 0;
  for (const auto &val : best) {
    sum += val.second;
    if (rand < sum) {
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
  Configuration tconfig(_supercell(), configdof);
  m_occ_loc.initialize(tconfig);
  jsonParser json;

  _log() << "    initial comp: " << to_json_array(CASM::comp(tconfig), json)
         << std::endl;
  _log() << "  initial comp_n: " << to_json_array(CASM::comp_n(tconfig), json)
         << std::endl;

  Eigen::VectorXd target_comp_n = conditions().mol_composition();
  _log() << "   target comp_n: " << to_json_array(target_comp_n, json)
         << std::endl;

  int count = 0;
  OccEvent e;
  ConfigDoF &tconfigdof = tconfig.configdof();
  while (true) {
    auto begin = m_cand.grand_canonical_swap().begin();
    auto end = m_cand.grand_canonical_swap().end();
    auto it = _find_grand_canonical_swap(tconfig, begin, end);

    if (it == end) {
      _log() << "   applied swaps: " << count << std::endl;
      break;
    }

    /// apply chosen swap (*it)
    m_occ_loc.propose_grand_canonical(e, *it, _mtrand());
    m_occ_loc.apply(e, tconfigdof);

    ++count;
  }

  _log() << "      final comp: " << to_json_array(CASM::comp(tconfig), json)
         << std::endl;
  _log() << "    final comp_n: " << to_json_array(CASM::comp_n(tconfig), json)
         << std::endl
         << std::endl;

  return tconfig.configdof();
}

}  // namespace Monte
}  // namespace CASM
