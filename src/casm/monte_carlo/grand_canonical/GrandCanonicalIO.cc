#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/monte_carlo/MonteIO_impl.hh"

namespace CASM {

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
  DataFormatter<ConstMonteCarloPtr> make_results_formatter(const GrandCanonical &mc) {

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
    for(int i = 0; i < mc.primclex().composition_axes().independent_compositions(); ++i) {
      formatter.push_back(MonteCarloParamChemPotFormatter<GrandCanonical>(mc, i));
    }

    for(int i = 0; i < mc.primclex().composition_axes().independent_compositions(); i++) {

      name = std::string("comp(") + mc.primclex().composition_axes().comp_var(i) + ")";
      formatter.push_back(MonteCarloMeanFormatter(name));
      formatter.push_back(MonteCarloPrecFormatter(name));
      exclude.insert(name);
    }

    // always sample comp_n
    auto struc_mol_name = mc.primclex().get_prim().get_struc_molecule_name();
    for(int i = 0; i < struc_mol_name.size(); ++i) {
      name = std::string("comp_n(") + struc_mol_name[i] + ")";
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
    formatter.push_back(MonteCarloHeatCapacityFormatter<GrandCanonical>());

    // include susc_x
    for(int i = 0; i < mc.primclex().composition_axes().independent_compositions(); i++) {
      for(int j = i; j < mc.primclex().composition_axes().independent_compositions(); j++) {

        auto comp_var_i = mc.primclex().composition_axes().comp_var(i);
        auto comp_var_j = mc.primclex().composition_axes().comp_var(j);
        formatter.push_back(MonteCarloSuscXFormatter<GrandCanonical>(comp_var_i, comp_var_j));

      }
    }

    // include thermo-chem susc
    for(int i = 0; i < mc.primclex().composition_axes().independent_compositions(); i++) {

      auto comp_var_i = mc.primclex().composition_axes().comp_var(i);
      formatter.push_back(MonteCarloThermoChemSuscXFormatter<GrandCanonical>(comp_var_i));

    }

    // include susc_n
    for(int i = 0; i < struc_mol_name.size(); ++i) {
      for(int j = i; j < struc_mol_name.size(); ++j) {

        auto species_i = struc_mol_name[i];
        auto species_j = struc_mol_name[j];
        formatter.push_back(MonteCarloSuscNFormatter<GrandCanonical>(species_i, species_j));

      }
    }

    // include thermo-chem susc
    for(int i = 0; i < struc_mol_name.size(); ++i) {

      auto species_i = struc_mol_name[i];
      formatter.push_back(MonteCarloThermoChemSuscNFormatter<GrandCanonical>(species_i));

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
  DataFormatter<ConstMonteCarloPtr> make_lte_results_formatter(const GrandCanonical &mc, const double &phi_LTE1, const std::string &configname) {

    DataFormatter<ConstMonteCarloPtr> formatter;

    bool print_json = true;
    formatter.push_back(ConstantValueFormatter<std::string, ConstMonteCarloPtr>("configname", configname, print_json));
    formatter.push_back(MonteCarloTFormatter<GrandCanonical>());
    formatter.push_back(GrandCanonicalLTEFormatter(phi_LTE1));
    std::set<std::string> exclude;
    std::string name;

    // always sample Beta, potential_energy, and formation_energy
    {
      formatter.push_back(MonteCarloBetaFormatter<GrandCanonical>());
      name = "gs_potential_energy";
      auto evaluator = [ = ](const ConstMonteCarloPtr & ptr) {
        return static_cast<const GrandCanonical *>(ptr)->potential_energy();
      };
      formatter.push_back(GenericDatumFormatter<double, ConstMonteCarloPtr>(name, name, evaluator));
      exclude.insert(name);
    }

    {
      name = "gs_formation_energy";
      auto evaluator = [ = ](const ConstMonteCarloPtr & ptr) {
        return static_cast<const GrandCanonical *>(ptr)->formation_energy();
      };
      formatter.push_back(GenericDatumFormatter<double, ConstMonteCarloPtr>(name, name, evaluator));
      exclude.insert(name);
    }



    // always sample param_chem_pot, comp
    for(int i = 0; i < mc.primclex().composition_axes().independent_compositions(); ++i) {
      formatter.push_back(MonteCarloParamChemPotFormatter<GrandCanonical>(mc, i));
    }

    for(int i = 0; i < mc.primclex().composition_axes().independent_compositions(); i++) {

      name = std::string("gs_comp(") + mc.primclex().composition_axes().comp_var(i) + ")";
      auto evaluator = [ = ](const ConstMonteCarloPtr & ptr) {
        const GrandCanonical *_ptr = static_cast<const GrandCanonical *>(ptr);
        return _ptr->primclex().composition_axes().param_composition(_ptr->comp_n())[i];
      };
      formatter.push_back(GenericDatumFormatter<double, ConstMonteCarloPtr>(name, name, evaluator));
      exclude.insert(name);
    }

    // always sample comp_n
    auto struc_mol_name = mc.primclex().get_prim().get_struc_molecule_name();
    for(int i = 0; i < struc_mol_name.size(); ++i) {
      name = std::string("gs_comp_n(") + struc_mol_name[i] + ")";
      auto evaluator = [ = ](const ConstMonteCarloPtr & ptr) {
        return static_cast<const GrandCanonical *>(ptr)->comp_n()[i];
      };
      formatter.push_back(GenericDatumFormatter<double, ConstMonteCarloPtr>(name, name, evaluator));
      exclude.insert(name);
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
  jsonParser &to_json(const GrandCanonicalConditions &conditions, jsonParser &json) {
    json.put_obj();
    json["temperature"] = conditions.temperature();
    json["param_chem_pot"] = jsonParser::object();
    auto param_chem_pot = conditions.param_chem_pot();
    for(int i = 0; i < param_chem_pot.size(); i++) {
      json["param_chem_pot"][CompositionConverter::comp_var(i)] = param_chem_pot(i);
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
  ///   "param_chem_pot" : {"a": number, "b": number, ...}
  /// }
  /// \endcode
  ///
  void from_json(GrandCanonicalConditions &conditions, const PrimClex &primclex, const jsonParser &json) {

    double temp = json["temperature"].get<double>();
    double tol = json["tolerance"].get<double>();

    int Nparam = primclex.composition_axes().independent_compositions();
    Eigen::VectorXd param_chem_pot(Nparam);

    for(int i = 0; i < Nparam; i++) {
      param_chem_pot(i) = json["param_chem_pot"][CompositionConverter::comp_var(i)].get<double>();
    }

    conditions = GrandCanonicalConditions(primclex, temp, param_chem_pot, tol);
  }


  /// \brief Print single spin flip LTE
  GenericDatumFormatter<double, ConstMonteCarloPtr> GrandCanonicalLTEFormatter(const double &phi_LTE1) {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      return phi_LTE1;
    };
    return GenericDatumFormatter<double, ConstMonteCarloPtr>("phi_LTE", "phi_LTE", evaluator);
  }

  /// \brief Will create new file or append to existing results file the results of the latest run
  void write_lte_results(const MonteSettings &settings, const GrandCanonical &mc, const double &phi_LTE1, const std::string &configname, Log &_log) {
    try {

      fs::create_directories(settings.output_directory());
      MonteCarloDirectoryStructure dir(settings.output_directory());
      auto formatter = make_lte_results_formatter(mc, phi_LTE1, configname);

      // write csv path results
      if(settings.write_csv()) {
        fs::path file = dir.results_csv();
        _log << "write: " << dir.results_csv() << "\n";
        fs::ofstream sout;

        if(!fs::exists(file)) {
          sout.open(file);
          formatter.print_header(&mc, sout);
        }
        else {
          sout.open(file, std::ios::app);
        }

        formatter.print(&mc, sout);

        sout.close();
      }

      // write json path results
      if(settings.write_json()) {
        fs::path file = dir.results_json();
        _log << "write: " << dir.results_json() << "\n";

        jsonParser results;
        if(fs::exists(file)) {
          results.read(file);
        }
        else {
          results = jsonParser::object();
        }

        formatter.to_json_arrays(&mc, results);
        results.write(file);
      }
    }
    catch(...) {
      std::cerr << "ERROR writing results" << std::endl;
      throw;
    }

  }
}
