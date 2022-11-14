#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/dataformatter/DataFormatter.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/global/definitions.hh"
#include "casm/monte_carlo/MonteIO_impl.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/monte_carlo/io/json/CorrMatchingPotential_json_io.hh"

namespace CASM {
namespace Monte {

/// \brief Make a results formatter
///
/// Output data is:
/// - N_equil_samples
/// - N_avg_samples
/// - T
/// - Beta
/// - <potential_energy>, prec(<potential_energy>)
/// - <formation_energy>, prec(<formation_energy>)
/// - param_chem_pot(a) ...
/// - <comp(a)> prec(<comp(a)>) ...
/// - <comp_n(a)> prec(<comp_n(a)>) ...
/// - (other properties: correlations, atom_frac, site_frac, etc.)
/// - heat_capacity
/// - susc_x(a,a), ...
/// - susc_x(S,a), ...
/// - susc_n(A,A), ...
/// - susc_n(S,A), ...
///
/// For JSON format is:
/// \code
/// { "key0":[...], "key1":[...], ... }
/// \endcode
///
DataFormatter<ConstMonteCarloPtr> make_results_formatter(
    const GrandCanonical &mc) {
  DataFormatter<ConstMonteCarloPtr> formatter;

  formatter.push_back(MonteCarloIsEquilibratedFormatter());
  formatter.push_back(MonteCarloNEquilSamplesFormatter());
  formatter.push_back(MonteCarloIsConvergedFormatter());
  formatter.push_back(MonteCarloNAvgSamplesFormatter());

  formatter.push_back(MonteCarloTFormatter<GrandCanonical>());

  std::set<std::string> exclude;
  std::string name;

  // always sample Beta, potential_energy, and formation_energy
  {
    formatter.push_back(MonteCarloBetaFormatter<GrandCanonical>());
    name = "potential_energy";
    formatter.push_back(MonteCarloMeanFormatter(name));
    formatter.push_back(MonteCarloPrecFormatter(name));
    exclude.insert(name);

    name = "formation_energy";
    formatter.push_back(MonteCarloMeanFormatter(name));
    formatter.push_back(MonteCarloPrecFormatter(name));
    exclude.insert(name);
  }

  // always sample param_chem_pot, comp
  for (int i = 0;
       i < mc.primclex().composition_axes().independent_compositions(); ++i) {
    formatter.push_back(MonteCarloParamChemPotFormatter<GrandCanonical>(mc, i));
  }

  for (int i = 0;
       i < mc.primclex().composition_axes().independent_compositions(); i++) {
    name = std::string("comp(") + mc.primclex().composition_axes().comp_var(i) +
           ")";
    formatter.push_back(MonteCarloMeanFormatter(name));
    formatter.push_back(MonteCarloPrecFormatter(name));
    exclude.insert(name);
  }

  // always sample comp_n
  auto struc_mol_name = xtal::struc_molecule_name(mc.primclex().prim());
  for (int i = 0; i < struc_mol_name.size(); ++i) {
    name = std::string("comp_n(") + struc_mol_name[i] + ")";
    formatter.push_back(MonteCarloMeanFormatter(name));
    formatter.push_back(MonteCarloPrecFormatter(name));
    exclude.insert(name);
  }

  // if has order parameter, always include order parameter
  if (mc.order_parameter() != nullptr) {
    for (int i = 0; i < mc.order_parameter()->dof_space().subspace_dim(); ++i) {
      name = std::string("order_parameter(") + std::to_string(i) + ")";
      formatter.push_back(MonteCarloMeanFormatter(name));
      formatter.push_back(MonteCarloPrecFormatter(name));
      exclude.insert(name);
    }
  }

  // include mean/prec of other properties
  for (auto it = mc.samplers().cbegin(); it != mc.samplers().cend(); ++it) {
    if (exclude.find(it->first) == exclude.end()) {
      formatter.push_back(MonteCarloMeanFormatter(it->first));
      formatter.push_back(MonteCarloPrecFormatter(it->first));
    }
  }

  // include heat_capacity
  formatter.push_back(MonteCarloHeatCapacityFormatter<GrandCanonical>());

  // include susc_x
  for (int i = 0;
       i < mc.primclex().composition_axes().independent_compositions(); i++) {
    for (int j = i;
         j < mc.primclex().composition_axes().independent_compositions(); j++) {
      auto comp_var_i = mc.primclex().composition_axes().comp_var(i);
      auto comp_var_j = mc.primclex().composition_axes().comp_var(j);
      formatter.push_back(
          MonteCarloSuscXFormatter<GrandCanonical>(comp_var_i, comp_var_j));
    }
  }

  // include thermo-chem susc
  for (int i = 0;
       i < mc.primclex().composition_axes().independent_compositions(); i++) {
    auto comp_var_i = mc.primclex().composition_axes().comp_var(i);
    formatter.push_back(
        MonteCarloThermoChemSuscXFormatter<GrandCanonical>(comp_var_i));
  }

  // include susc_n
  for (int i = 0; i < struc_mol_name.size(); ++i) {
    for (int j = i; j < struc_mol_name.size(); ++j) {
      auto species_i = struc_mol_name[i];
      auto species_j = struc_mol_name[j];
      formatter.push_back(
          MonteCarloSuscNFormatter<GrandCanonical>(species_i, species_j));
    }
  }

  // include thermo-chem susc
  for (int i = 0; i < struc_mol_name.size(); ++i) {
    auto species_i = struc_mol_name[i];
    formatter.push_back(
        MonteCarloThermoChemSuscNFormatter<GrandCanonical>(species_i));
  }

  return formatter;
}

/// \brief Make a LTE results formatter
///
/// Output data is:
/// - T
/// - phi_LTE
/// - Beta
/// - configname
/// - gs_potential_energy
/// - gs_formation_energy
/// - param_chem_pot(a) ...
/// - gs_comp(a) ...
/// - gs_comp_n(a) ...
///
/// For JSON format is:
/// \code
/// { "key0":[...], "key1":[...], ... }
/// \endcode
///
DataFormatter<ConstMonteCarloPtr> make_lte_results_formatter(
    const GrandCanonical &mc, const double &phi_LTE1,
    const std::string &configname) {
  DataFormatter<ConstMonteCarloPtr> formatter;

  bool print_json = true;
  formatter.push_back(ConstantValueFormatter<std::string, ConstMonteCarloPtr>(
      "configname", configname, print_json));
  formatter.push_back(MonteCarloTFormatter<GrandCanonical>());
  formatter.push_back(GrandCanonicalLTEFormatter(phi_LTE1));
  std::set<std::string> exclude;
  std::string name;

  // always sample Beta, potential_energy, and formation_energy
  {
    formatter.push_back(MonteCarloBetaFormatter<GrandCanonical>());
    name = "gs_potential_energy";
    auto evaluator = [=](const ConstMonteCarloPtr &ptr) {
      return static_cast<const GrandCanonical *>(ptr)->potential_energy();
    };
    formatter.push_back(GenericDatumFormatter<double, ConstMonteCarloPtr>(
        name, name, evaluator));
    exclude.insert(name);
  }

  {
    name = "gs_formation_energy";
    auto evaluator = [=](const ConstMonteCarloPtr &ptr) {
      return static_cast<const GrandCanonical *>(ptr)->formation_energy();
    };
    formatter.push_back(GenericDatumFormatter<double, ConstMonteCarloPtr>(
        name, name, evaluator));
    exclude.insert(name);
  }

  // always sample param_chem_pot, comp
  for (int i = 0;
       i < mc.primclex().composition_axes().independent_compositions(); ++i) {
    formatter.push_back(MonteCarloParamChemPotFormatter<GrandCanonical>(mc, i));
  }

  for (int i = 0;
       i < mc.primclex().composition_axes().independent_compositions(); i++) {
    name = std::string("gs_comp(") +
           mc.primclex().composition_axes().comp_var(i) + ")";
    auto evaluator = [=](const ConstMonteCarloPtr &ptr) {
      const GrandCanonical *_ptr = static_cast<const GrandCanonical *>(ptr);
      return _ptr->primclex().composition_axes().param_composition(
          _ptr->comp_n())[i];
    };
    formatter.push_back(GenericDatumFormatter<double, ConstMonteCarloPtr>(
        name, name, evaluator));
    exclude.insert(name);
  }

  // always sample comp_n
  auto struc_mol_name = xtal::struc_molecule_name(mc.primclex().prim());
  for (int i = 0; i < struc_mol_name.size(); ++i) {
    name = std::string("gs_comp_n(") + struc_mol_name[i] + ")";
    auto evaluator = [=](const ConstMonteCarloPtr &ptr) {
      return static_cast<const GrandCanonical *>(ptr)->comp_n()[i];
    };
    formatter.push_back(GenericDatumFormatter<double, ConstMonteCarloPtr>(
        name, name, evaluator));
    exclude.insert(name);
  }

  // if has order parameter, always include order parameter
  if (mc.order_parameter() != nullptr) {
    for (int i = 0; i < mc.order_parameter()->dof_space().subspace_dim(); ++i) {
      name = std::string("order_parameter(") + std::to_string(i) + ")";
      formatter.push_back(MonteCarloMeanFormatter(name));
      formatter.push_back(MonteCarloPrecFormatter(name));
      exclude.insert(name);
    }
  }

  return formatter;
}

/// \brief Store GrandCanonicalConditions in JSON format
///
/// \code
/// {
///   "temperature" : number,
///   "tol: number,
///   "param_chem_pot" : {"a": number, "b": number, ...}
/// }
/// \endcode
///
jsonParser &to_json(const GrandCanonicalConditions &conditions,
                    jsonParser &json) {
  json.put_obj();
  json["temperature"] = conditions.temperature();

  if (conditions.include_param_chem_pot()) {
    json["param_chem_pot"] = jsonParser::object();
    auto param_chem_pot = conditions.param_chem_pot();
    for (int i = 0; i < param_chem_pot.size(); i++) {
      json["param_chem_pot"][CompositionConverter::comp_var(i)] =
          param_chem_pot(i);
    }
  } else {
    json["param_chem_pot"].put_null();
  }

  json["tolerance"] = conditions.tolerance();

  json["include_formation_energy"] = conditions.include_formation_energy();

  to_json(conditions.param_comp_quad_pot_target(),
          json["param_comp_quad_pot_target"], jsonParser::as_array());
  to_json(conditions.param_comp_quad_pot_vector(),
          json["param_comp_quad_pot_vector"], jsonParser::as_array());
  json["param_comp_quad_pot_matrix"] = conditions.param_comp_quad_pot_matrix();

  to_json(conditions.order_parameter_pot(), json["order_parameter_pot"],
          jsonParser::as_array());

  to_json(conditions.order_parameter_quad_pot_target(),
          json["order_parameter_quad_pot_target"], jsonParser::as_array());
  to_json(conditions.order_parameter_quad_pot_vector(),
          json["order_parameter_quad_pot_vector"], jsonParser::as_array());
  json["order_parameter_quad_pot_matrix"] =
      conditions.order_parameter_quad_pot_matrix();

  json["corr_matching_pot"] = conditions.corr_matching_pot();
  json["random_alloy_corr_matching_pot"] =
      conditions.random_alloy_corr_matching_pot();

  return json;
}

/// \brief Read GrandCanonicalConditions from JSON format
///
/// \code
/// {
///   "temperature" : number,
///   "tol: number,
///   // the following are optional and can also be "null"
///   "param_chem_pot" : {"a": number, "b": number, ...}
///   "include_formation_energy": bool // default=true, except default=false for
///                                    //   corr_matching_pot
///   "param_comp_quad_pot_target": [number, number, ...] // potential minimum
///   "param_comp_quad_pot_vector": [number, number, ...] // curvature (vector)
///   "param_comp_quad_pot_matrix": [ // curvature (matrix)
///     [number, number, ...],
///     [number, number, ...],
///     ...]
///   "order_parameter_pot": [number, number, ...] // linear potential
///
///   "order_parameter_quad_pot_target": [number, number, ...] // potential
///   minimum "order_parameter_quad_pot_vector": [number, number, ...] //
///   curvature (vector) "order_parameter_quad_pot_matrix": [ // curvature
///   (matrix)
///     [number, number, ...],
///     [number, number, ...],
///     ...]
///   "corr_matching_pot": {
///     "exact_matching_weight": number // default=0.
///     "tol": number // default=1e-5
///     "targets": [
///       {"index":int, "value":number, "weight":number (default=1.)},
///       ...
///     ]
///   "random_alloy_corr_matching_pot": {
///     "exact_matching_weight": number // default=0.
///     "tol": number // default=1e-5
///     "sublattice_prob": [
///       [number, number, ...], // sublattice 0
///       [number, number, ...], // sublattice 1
///       ...
///     ]
///   }
/// }
/// \endcode
///
void from_json(GrandCanonicalConditions &conditions, const PrimClex &primclex,
               const jsonParser &json, GrandCanonical const &mc) {
  bool include_formation_energy = true;
  bool include_param_chem_pot = false;

  double temp = json["temperature"].get<double>();
  double tol = json["tolerance"].get<double>();

  int Nparam = primclex.composition_axes().independent_compositions();
  Eigen::VectorXd param_chem_pot(Nparam);
  param_chem_pot.setZero();
  if (json.contains("param_chem_pot") && !json["param_chem_pot"].is_null()) {
    include_param_chem_pot = true;
    if (json["param_chem_pot"].is_obj()) {
      for (int i = 0; i < Nparam; i++) {
        param_chem_pot(i) =
            json["param_chem_pot"][CompositionConverter::comp_var(i)]
                .get<double>();
      }
    } else if (json["param_chem_pot"].is_array()) {
      CASM::from_json(param_chem_pot, json["param_chem_pot"]);
      if (param_chem_pot.size() != Nparam) {
        throw std::runtime_error(
            "Error parsing GrandCanonicalConditions JSON: param_chem_pot size "
            "mismatch with composition axes");
      }
    } else {
      throw std::runtime_error(
          "Error parsing GrandCanonicalConditions JSON: param_chem_pot must be "
          "an object or array");
    }
  }
  // --- quadratic potential of parametric composition

  std::optional<Eigen::VectorXd> param_comp_quad_pot_target;
  json.get_if(param_comp_quad_pot_target, "param_comp_quad_pot_target");

  std::optional<Eigen::VectorXd> param_comp_quad_pot_vector;
  json.get_if(param_comp_quad_pot_vector, "param_comp_quad_vector");

  std::optional<Eigen::MatrixXd> param_comp_quad_pot_matrix;
  json.get_if(param_comp_quad_pot_matrix, "param_comp_quad_pot_matrix");

  // --- linear potential of order parameter

  std::optional<Eigen::VectorXd> order_parameter_pot;
  json.get_if(order_parameter_pot, "order_parameter_pot");

  // --- quadratic potential of order parameter

  std::optional<Eigen::VectorXd> order_parameter_quad_pot_target;
  json.get_if(order_parameter_quad_pot_target,
              "order_parameter_quad_pot_target");

  std::optional<Eigen::VectorXd> order_parameter_quad_pot_vector;
  json.get_if(order_parameter_quad_pot_vector,
              "order_parameter_quad_pot_vector");

  std::optional<Eigen::MatrixXd> order_parameter_quad_pot_matrix;
  json.get_if(order_parameter_quad_pot_matrix,
              "order_parameter_quad_pot_matrix");

  // --- correlation matching potential

  std::optional<CorrMatchingParams> corr_matching_pot;
  json.get_if(corr_matching_pot, "corr_matching_pot");
  if (corr_matching_pot.has_value()) {
    include_formation_energy = false;
  }

  std::optional<RandomAlloyCorrMatchingParams> random_alloy_corr_matching_pot;
  if (json.contains("random_alloy_corr_matching_pot") &&
      !json["random_alloy_corr_matching_pot"].is_null()) {
    random_alloy_corr_matching_pot = RandomAlloyCorrMatchingParams();
    random_alloy_corr_matching_pot->random_alloy_corr_f =
        mc.random_alloy_corr_f();
    from_json(*random_alloy_corr_matching_pot,
              json["random_alloy_corr_matching_pot"]);
    include_formation_energy = false;
  }

  // --- override inclusion of formation energy ---
  json.get_if(include_formation_energy, "include_formation_energy");

  conditions = GrandCanonicalConditions(
      primclex, temp, param_chem_pot, tol, include_formation_energy,
      include_param_chem_pot, param_comp_quad_pot_target,
      param_comp_quad_pot_vector, param_comp_quad_pot_matrix,
      order_parameter_pot, order_parameter_quad_pot_target,
      order_parameter_quad_pot_vector, order_parameter_quad_pot_matrix,
      corr_matching_pot, random_alloy_corr_matching_pot);
}

/// \brief Print single spin flip LTE
GenericDatumFormatter<double, ConstMonteCarloPtr> GrandCanonicalLTEFormatter(
    const double &phi_LTE1) {
  auto evaluator = [=](const ConstMonteCarloPtr &mc) { return phi_LTE1; };
  return GenericDatumFormatter<double, ConstMonteCarloPtr>("phi_LTE", "phi_LTE",
                                                           evaluator);
}

/// \brief Will create new file or append to existing results file the results
/// of the latest run
void write_lte_results(const MonteSettings &settings, const GrandCanonical &mc,
                       const double &phi_LTE1, const std::string &configname,
                       Log &_log) {
  try {
    fs::create_directories(settings.output_directory());
    MonteCarloDirectoryStructure dir(settings.output_directory());
    auto formatter = make_lte_results_formatter(mc, phi_LTE1, configname);

    // write csv path results
    if (settings.write_csv()) {
      fs::path file = dir.results_csv();
      _log << "write: " << dir.results_csv() << "\n";
      fs::ofstream sout;

      if (!fs::exists(file)) {
        sout.open(file);
        formatter.print_header(&mc, sout);
      } else {
        sout.open(file, std::ios::app);
      }

      formatter.print(&mc, sout);

      sout.close();
    }

    // write json path results
    if (settings.write_json()) {
      fs::path file = dir.results_json();
      _log << "write: " << dir.results_json() << "\n";

      jsonParser results;
      if (fs::exists(file)) {
        results.read(file);
      } else {
        results = jsonParser::object();
      }

      formatter.to_json_arrays(&mc, results);
      results.write(file);
    }
  } catch (...) {
    std::cerr << "ERROR writing results" << std::endl;
    throw;
  }
}

}  // namespace Monte
}  // namespace CASM
