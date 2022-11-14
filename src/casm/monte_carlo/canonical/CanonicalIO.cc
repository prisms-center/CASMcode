#include "casm/monte_carlo/canonical/CanonicalIO.hh"

#include "casm/casm_io/dataformatter/DataFormatter.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/global/definitions.hh"
#include "casm/monte_carlo/MonteIO_impl.hh"
#include "casm/monte_carlo/canonical/Canonical.hh"
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
/// - <comp(a)> prec(<comp(a)>) ...
/// - <comp_n(a)> prec(<comp_n(a)>) ...
/// - (other properties: correlations, atom_frac, site_frac, etc.)
/// - heat_capacity
///
/// For JSON format is:
/// \code
/// { "key0":[...], "key1":[...], ... }
/// \endcode
///
DataFormatter<ConstMonteCarloPtr> make_results_formatter(const Canonical &mc) {
  DataFormatter<ConstMonteCarloPtr> formatter;

  formatter.push_back(MonteCarloIsEquilibratedFormatter());
  formatter.push_back(MonteCarloNEquilSamplesFormatter());
  formatter.push_back(MonteCarloIsConvergedFormatter());
  formatter.push_back(MonteCarloNAvgSamplesFormatter());

  std::set<std::string> exclude;
  std::string name;

  // always print temperature, Beta
  formatter.push_back(MonteCarloTFormatter<Canonical>());
  formatter.push_back(MonteCarloBetaFormatter<Canonical>());

  // always print comp
  for (int i = 0;
       i < mc.primclex().composition_axes().independent_compositions(); ++i) {
    name = std::string("comp(") + mc.primclex().composition_axes().comp_var(i) +
           ")";
    formatter.push_back(MonteCarloCompFormatter<Canonical>(mc, i));
    exclude.insert(name);
  }

  // always print comp_n
  auto struc_mol_name = xtal::struc_molecule_name(mc.primclex().prim());
  for (int i = 0; i < struc_mol_name.size(); ++i) {
    name = std::string("comp_n(") + struc_mol_name[i] + ")";
    formatter.push_back(MonteCarloCompNFormatter<Canonical>(mc, i));
    exclude.insert(name);
  }

  // always print potential_energy, and formation_energy, though they are the
  // same
  {
    name = "potential_energy";
    formatter.push_back(MonteCarloMeanFormatter(name));
    formatter.push_back(MonteCarloPrecFormatter(name));
    exclude.insert(name);

    name = "formation_energy";
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
  formatter.push_back(MonteCarloHeatCapacityFormatter<Canonical>());

  return formatter;
}

/// \brief Store GrandCanonicalConditions in JSON format
///
/// \code
/// {
///   "temperature" : number,
///   "tol: number,
///   "comp" : {"a": number, "b": number, ...},
///   "comp_n" : {"A": number, "B": number, ...}
/// }
/// \endcode
///
jsonParser &to_json(const CanonicalConditions &conditions, jsonParser &json) {
  json.put_obj();
  json["temperature"] = conditions.temperature();

  json["comp"] = jsonParser::object();
  auto param_composition = conditions.param_composition();
  for (int i = 0; i < param_composition.size(); i++) {
    json["comp"][CompositionConverter::comp_var(i)] = param_composition(i);
  }

  json["comp_n"] = jsonParser::object();
  auto mol_composition = conditions.mol_composition();
  auto components = conditions.primclex().composition_axes().components();
  for (int i = 0; i < components.size(); ++i) {
    json["comp_n"][components[i]] = mol_composition(i);
  }

  json["tolerance"] = conditions.tolerance();

  json["include_formation_energy"] = conditions.include_formation_energy();

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
///   "tolerance": number,
///   "comp" : {"a": number, "b": number, ...} // option 1: comp object
///   "comp" : [number, number, ...] // option 2: comp array
///   "comp_n" : {"A": number, "B": number, ...} // option 3: comp_n object
///   "comp_n" : [number, number, ...] // option 2: comp_n array
///
///   // the following are optional / can also be "null"
///   "include_formation_energy": bool // default=true, except default=false for
///   corr_matching_pot "param_comp_quad_pot_target": [number, number, ...] //
///   option 2: potential minimum "param_comp_quad_pot_vector": [number, number,
///   ...] // curvature (vector) "param_comp_quad_pot_matrix": [ // curvature
///   (matrix)
///     [number, number, ...],
///     [number, number, ...],
///     ...]
///   "order_parameter_pot": [number, number, ...] // linear potential
///
///   // potential minimum
///   "order_parameter_quad_pot_target": [number, number, ...]
///
///   // curvature (vector)
///   "order_parameter_quad_pot_vector": [number, number, ...]
///   "order_parameter_quad_pot_matrix": [ // curvature (matrix)
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
void from_json(CanonicalConditions &conditions, const PrimClex &primclex,
               const jsonParser &json, Canonical const &mc, bool incremental) {
  bool include_formation_energy = true;

  double temp = json["temperature"].get<double>();
  double tol = json["tolerance"].get<double>();

  Eigen::VectorXd comp;

  if (json.contains("comp")) {
    if (json["comp"].is_array()) {
      from_json(comp, json["comp"]);
    } else {
      int Nparam = primclex.composition_axes().independent_compositions();
      comp.resize(Nparam);

      for (int i = 0; i < Nparam; i++) {
        comp(i) = json["comp"][CompositionConverter::comp_var(i)].get<double>();
      }
    }
  } else if (json.contains("comp_n")) {
    Eigen::VectorXd comp_n;

    if (json["comp_n"].is_array()) {
      from_json(comp_n, json["comp_n"]);
    } else {
      auto components = primclex.composition_axes().components();
      comp_n.resize(components.size());

      for (int i = 0; i < components.size(); i++) {
        comp_n(i) = json["comp_n"][components[i]].get<double>();
      }
    }

    if (!incremental) {
      comp_n = comp_n / comp_n.sum() * primclex.prim().basis().size();
      comp = primclex.composition_axes().param_composition(comp_n);
    } else {
      comp = primclex.composition_axes().dparam_composition(comp_n);
    }
  } else {
    throw std::runtime_error(
        "Error reading conditions: No 'comp' or 'comp_n' specified");
  }

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

  conditions = CanonicalConditions(
      primclex, temp, comp, tol, include_formation_energy, order_parameter_pot,
      order_parameter_quad_pot_target, order_parameter_quad_pot_vector,
      order_parameter_quad_pot_matrix, corr_matching_pot,
      random_alloy_corr_matching_pot);
}

}  // namespace Monte
}  // namespace CASM
