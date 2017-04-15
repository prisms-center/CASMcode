#include "casm/monte_carlo/canonical/CanonicalIO.hh"

#include "casm/external/gzstream/gzstream.h"
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/VaspIO.hh"
#include "casm/monte_carlo/canonical/Canonical.hh"
#include "casm/monte_carlo/MonteIO_impl.hh"

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
      for(int i = 0; i < mc.primclex().composition_axes().independent_compositions(); ++i) {
        name = std::string("comp(") + mc.primclex().composition_axes().comp_var(i) + ")";
        formatter.push_back(MonteCarloCompFormatter<Canonical>(mc, i));
        exclude.insert(name);
      }

      // always print comp_n
      auto struc_mol_name = mc.primclex().get_prim().get_struc_molecule_name();
      for(int i = 0; i < struc_mol_name.size(); ++i) {
        name = std::string("comp_n(") + struc_mol_name[i] + ")";
        formatter.push_back(MonteCarloCompNFormatter<Canonical>(mc, i));
        exclude.insert(name);
      }

      // always print potential_energy, and formation_energy, though they are the same
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

      // include mean/prec of other properties
      for(auto it = mc.samplers().cbegin(); it != mc.samplers().cend(); ++it) {
        if(exclude.find(it->first) == exclude.end()) {
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
      for(int i = 0; i < param_composition.size(); i++) {
        json["comp"][CompositionConverter::comp_var(i)] = param_composition(i);
      }

      json["comp_n"] = jsonParser::object();
      auto mol_composition = conditions.mol_composition();
      auto components = conditions.primclex().composition_axes().components();
      for(int i = 0; i < components.size(); ++i) {
        json["comp_n"][components[i]] = mol_composition(i);
      }

      json["tolerance"] = conditions.tolerance();

      return json;
    }

    /// \brief Read GrandCanonicalConditions from JSON format
    ///
    /// \code
    /// {
    ///   "temperature" : number,
    ///   "tol: number,
    ///   "comp" : {"a": number, "b": number, ...} // option 1: comp object
    ///   "comp" : [number, number, ...] // option 2: comp array
    ///   "comp_n" : {"A": number, "B": number, ...} // option 3: comp_n object
    ///   "comp_n" : [number, number, ...] // option 2: comp_n array
    /// }
    /// \endcode
    ///
    void from_json(CanonicalConditions &conditions, const PrimClex &primclex, const jsonParser &json) {

      double temp = json["temperature"].get<double>();
      double tol = json["tolerance"].get<double>();

      Eigen::VectorXd comp;

      if(json.contains("comp")) {

        if(json["comp"].is_array()) {
          from_json(comp, json["comp"]);
        }
        else {
          int Nparam = primclex.composition_axes().independent_compositions();
          comp.resize(Nparam);

          for(int i = 0; i < Nparam; i++) {
            comp(i) = json["comp"][CompositionConverter::comp_var(i)].get<double>();
          }
        }
      }
      else if(json.contains("comp_n")) {

        Eigen::VectorXd comp_n;

        if(json["comp_n"].is_array()) {
          from_json(comp_n, json["comp_n"]);
        }
        else {
          auto components = primclex.composition_axes().components();
          comp_n.resize(components.size());

          for(int i = 0; i < components.size(); i++) {
            comp_n(i) = json["comp_n"][components[i]].get<double>();
          }
        }

        comp_n = comp_n / comp_n.sum() * primclex.get_prim().basis.size();
        comp = primclex.composition_axes().param_composition(comp_n);
      }
      else {
        throw std::runtime_error("Error reading conditions: No 'comp' or 'comp_n' specified");
      }
      conditions = CanonicalConditions(primclex, temp, comp, tol);
    }
  }
}
