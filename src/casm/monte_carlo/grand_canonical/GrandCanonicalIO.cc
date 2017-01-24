#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

#include "casm/external/gzstream/gzstream.h"
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/casm_io/VaspIO.hh"

namespace CASM {

  /// \brief Print heat capacity, 'heat_capacity'
  ///
  /// - Heat capacity (per unit cell) = cov(potential_energy, potential_energy)*N/(k*T*T)
  GenericDatumFormatter<double, ConstMonteCarloPtr> GrandCanonicalHeatCapacityFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator("potential_energy", "potential_energy");
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const GrandCanonical *>(ptr)->conditions().temperature();
      return cov_evaluator(mc) * N / (KB * T * T);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("heat_capacity");

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Print parametric susceptibility, 'susc_x(a,b)'
  ///
  /// \arg comp_var_i: First parametric composition variable "a", "b", "c", ...
  /// \arg comp_var_j: Second parametric composition variable "a", "b", "c", ...
  ///
  /// - Heat capacity (per unit cell) = cov(x_i, x_j)*N/(k*T)
  GenericDatumFormatter<double, ConstMonteCarloPtr> GrandCanonicalSuscXFormatter(std::string comp_var_i, std::string comp_var_j) {

    std::string name_i = "comp(" + comp_var_i + ")";
    std::string name_j = "comp(" + comp_var_j + ")";

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator(name_i, name_j);
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const GrandCanonical *>(ptr)->conditions().temperature();
      return cov_evaluator(mc) * N / (KB * T);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("susc_x(" + comp_var_i + "," + comp_var_j + ")");

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Print parametric susceptibility, 'susc_x(a,b)'
  ///
  /// \arg species_i: First species "A", "B", "C", ...
  /// \arg species_j: Second species "A", "B", "C", ...
  ///
  /// - Heat capacity (per unit cell) = cov(x_i, x_j)*N/(k*T)
  GenericDatumFormatter<double, ConstMonteCarloPtr> GrandCanonicalSuscNFormatter(std::string species_i, std::string species_j) {

    std::string name_i = "comp_n(" + species_i + ")";
    std::string name_j = "comp_n(" + species_j + ")";

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator(name_i, name_j);
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const GrandCanonical *>(ptr)->conditions().temperature();
      return cov_evaluator(mc) * N / (KB * T);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("susc_n(" + species_i + "," + species_j + ")");

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Print parametric thermo-chemical susceptibility, 'susc_x(S,a)'
  ///
  /// \arg species_i: First species "A", "B", "C", ...
  ///
  /// - thermo-chemical susceptibility (per unit cell) = cov(potential_energy, x_i)*N/(k*T)
  GenericDatumFormatter<double, ConstMonteCarloPtr>
  GrandCanonicalThermoChemSuscXFormatter(std::string comp_var_i) {

    std::string name_i = "comp(" + comp_var_i + ")";

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator("potential_energy", name_i);
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const GrandCanonical *>(ptr)->conditions().temperature();
      return cov_evaluator(mc) * N / (KB * T);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("susc_x(S," + comp_var_i + ")");

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Print thermo-chemical susceptibility, 'susc_n(S,A)'
  GenericDatumFormatter<double, ConstMonteCarloPtr>
  GrandCanonicalThermoChemSuscNFormatter(std::string species_i) {

    std::string name_i = "comp_n(" + species_i + ")";

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator("potential_energy", name_i);
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const GrandCanonical *>(ptr)->conditions().temperature();
      return cov_evaluator(mc) * N / (KB * T);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("susc_n(S," + species_i + ")");

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

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
    auto struc_mol_name = mc.primclex().prim().struc_molecule_name();
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
    formatter.push_back(GrandCanonicalHeatCapacityFormatter());

    // include susc_x
    for(int i = 0; i < mc.primclex().composition_axes().independent_compositions(); i++) {
      for(int j = i; j < mc.primclex().composition_axes().independent_compositions(); j++) {

        auto comp_var_i = mc.primclex().composition_axes().comp_var(i);
        auto comp_var_j = mc.primclex().composition_axes().comp_var(j);
        formatter.push_back(GrandCanonicalSuscXFormatter(comp_var_i, comp_var_j));

      }
    }

    // include thermo-chem susc
    for(int i = 0; i < mc.primclex().composition_axes().independent_compositions(); i++) {

      auto comp_var_i = mc.primclex().composition_axes().comp_var(i);
      formatter.push_back(GrandCanonicalThermoChemSuscXFormatter(comp_var_i));

    }

    // include susc_n
    for(int i = 0; i < struc_mol_name.size(); ++i) {
      for(int j = i; j < struc_mol_name.size(); ++j) {

        auto species_i = struc_mol_name[i];
        auto species_j = struc_mol_name[j];
        formatter.push_back(GrandCanonicalSuscNFormatter(species_i, species_j));

      }
    }

    // include thermo-chem susc
    for(int i = 0; i < struc_mol_name.size(); ++i) {

      auto species_i = struc_mol_name[i];
      formatter.push_back(GrandCanonicalThermoChemSuscNFormatter(species_i));

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

    formatter.push_back(ConstantValueFormatter<std::string, ConstMonteCarloPtr>("configname", configname));
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
    auto struc_mol_name = mc.primclex().prim().struc_molecule_name();
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

  /// \brief Make a observation formatter
  ///
  /// For csv:
  /// \code
  /// # Pass Step X1 X2 ...
  /// \endcode
  ///
  /// For JSON:
  /// \code
  /// {"Pass/Step":[...], "X":[...], ...}
  /// \endcode
  ///
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_observation_formatter(const GrandCanonical &mc) {

    DataFormatter<std::pair<ConstMonteCarloPtr, Index> > formatter;

    formatter.push_back(MonteCarloPassFormatter());
    formatter.push_back(MonteCarloStepFormatter());
    for(auto it = mc.samplers().cbegin(); it != mc.samplers().cend(); ++it) {
      formatter.push_back(MonteCarloObservationFormatter(it->first));
    }
    return formatter;
  }

  /// \brief Make a trajectory formatter
  ///
  /// For csv:
  /// \code
  /// Pass Step occ(0) occ(1) ...
  /// \endcode
  ///
  /// For JSON:
  /// \code
  /// {"Pass" : [...], "Step":[...], "occ":[[...]]}
  /// \endcode
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_trajectory_formatter(const GrandCanonical &mc) {
    DataFormatter<std::pair<ConstMonteCarloPtr, Index> > formatter;
    formatter.push_back(MonteCarloPassFormatter());
    formatter.push_back(MonteCarloStepFormatter());

    // this is probably not the best way...
    for(Index i = 0; i < mc.configdof().occupation().size(); ++i) {
      formatter.push_back(MonteCarloOccFormatter(i));
    }
    return formatter;
  }

  /// \brief Store GrandCanonicalConditions in JSON format
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
  void from_json(GrandCanonicalConditions &conditions, const CompositionConverter &comp_converter, const jsonParser &json) {

    double temp = json["temperature"].get<double>();
    double tol = json["tolerance"].get<double>();

    int Nparam = comp_converter.independent_compositions();
    Eigen::VectorXd param_chem_pot(Nparam);

    for(int i = 0; i < Nparam; i++) {
      param_chem_pot(i) = json["param_chem_pot"][CompositionConverter::comp_var(i)].get<double>();
    }

    conditions = GrandCanonicalConditions(temp, param_chem_pot, comp_converter, tol);
  }

  /// \brief Will create new file or append to existing results file the results of the latest run
  void write_results(const MonteSettings &settings, const GrandCanonical &mc, Log &_log) {
    try {

      fs::create_directories(settings.output_directory());
      GrandCanonicalDirectoryStructure dir(settings.output_directory());
      auto formatter = make_results_formatter(mc);

      // write csv path results
      if(settings.write_csv()) {
        _log << "write: " << dir.results_csv() << "\n";
        fs::path file = dir.results_csv();
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
        _log << "write: " << dir.results_json() << "\n";
        fs::path file = dir.results_json();

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

  /// \brief Write conditions to conditions.cond_index directory
  void write_conditions_json(const MonteSettings &settings, const GrandCanonical &mc, Index cond_index, Log &_log) {
    try {
      GrandCanonicalDirectoryStructure dir(settings.output_directory());
      fs::create_directories(dir.conditions_dir(cond_index));
      jsonParser json;
      to_json(mc.conditions(), json);
      _log << "write: " << dir.conditions_json(cond_index) << "\n";
      json.write(dir.conditions_json(cond_index));
    }
    catch(...) {
      std::cerr << "ERROR writing conditions.json" << std::endl;
      throw;
    }
  }

  /// \brief Will create (and possibly overwrite) new file with all observations from run with conditions.cond_index
  void write_observations(const MonteSettings &settings, const GrandCanonical &mc, Index cond_index, Log &_log) {
    try {
      if(!settings.write_observations()) {
        return;
      }

      GrandCanonicalDirectoryStructure dir(settings.output_directory());
      fs::create_directories(dir.conditions_dir(cond_index));
      auto formatter = make_observation_formatter(mc);

      std::vector<std::pair<ConstMonteCarloPtr, Index> > observations;
      ConstMonteCarloPtr ptr = &mc;
      for(MonteSampler::size_type i = 0; i < mc.sample_times().size(); ++i) {
        observations.push_back(std::make_pair(ptr, i));
      }

      if(settings.write_csv()) {
        gz::ogzstream sout((dir.observations_csv(cond_index).string() + ".gz").c_str());
        _log << "write: " << fs::path(dir.observations_csv(cond_index).string() + ".gz") << "\n";
        sout << formatter(observations.cbegin(), observations.cend());
        sout.close();
      }

      if(settings.write_json()) {
        gz::ogzstream sout((dir.observations_json(cond_index).string() + ".gz").c_str());
        _log << "write: " << fs::path(dir.observations_json(cond_index).string() + ".gz") << "\n";
        jsonParser json = jsonParser::object();
        formatter(observations.cbegin(), observations.cend()).to_json_arrays(json);
        sout << json;
        sout.close();
      }

    }
    catch(...) {
      std::cerr << "ERROR writing observations." << std::endl;
      throw;
    }
  }

  /// \brief Will create (and possibly overwrite) new file with all observations from run with conditions.cond_index
  ///
  /// - Also writes "occupation_key" file giving occupant index -> species for each prim basis site
  ///
  /// For csv:
  /// \code
  /// site_index occ_index_0 occ_index_1 ...
  /// 0          Ni          Al
  /// 1          Ni          -
  /// ...
  /// \endcode
  ///
  /// For JSON:
  /// \code
  /// [["A", "B"],["A" "C"], ... ]
  /// \endcode
  ///
  void write_trajectory(const MonteSettings &settings, const GrandCanonical &mc, Index cond_index, Log &_log) {
    try {

      if(!settings.write_trajectory()) {
        return;
      }

      GrandCanonicalDirectoryStructure dir(settings.output_directory());
      fs::create_directories(dir.conditions_dir(cond_index));
      auto formatter = make_trajectory_formatter(mc);
      const Structure &prim = mc.primclex().prim();

      std::vector<std::pair<ConstMonteCarloPtr, Index> > observations;
      ConstMonteCarloPtr ptr = &mc;
      for(MonteSampler::size_type i = 0; i < mc.sample_times().size(); ++i) {
        observations.push_back(std::make_pair(ptr, i));
      }

      if(settings.write_csv()) {
        gz::ogzstream sout((dir.trajectory_csv(cond_index).string() + ".gz").c_str());
        _log << "write: " << fs::path(dir.trajectory_csv(cond_index).string() + ".gz") << "\n";
        sout << formatter(observations.cbegin(), observations.cend());
        sout.close();

        int max_allowed = 0;
        for(int i = 0; i < prim.basis.size(); i++) {
          if(prim.basis[i].allowed_occupants().size() > max_allowed) {
            max_allowed = prim.basis[i].allowed_occupants().size();
          }
        }

        // --- Write "occupation_key.csv" ------------
        //
        // site_index occ_index_0 occ_index_1 ...
        // 0          Ni          Al
        // 1          Ni          -
        // ...
        fs::ofstream keyout(dir.occupation_key_csv());
        _log << "write: " << dir.occupation_key_csv() << "\n";
        keyout << "site_index";
        for(int i = 0; i < max_allowed; i++) {
          keyout << "\tocc_index_" << i;
        }
        keyout << "\n";

        for(int i = 0; i < prim.basis.size(); i++) {
          keyout << i;
          for(int j = 0; j < max_allowed; j++) {
            if(j < prim.basis[i].allowed_occupants().size()) {
              keyout << "\t" << prim.basis[i].allowed_occupants()[j];
            }
            else {
              keyout << "\t-";
            }
          }
          keyout << "\n";
        }
        keyout.close();

      }

      if(settings.write_json()) {

        jsonParser json = jsonParser::object();
        json["Pass"] = jsonParser::array();
        json["Step"] = jsonParser::array();
        json["DoF"] = jsonParser::array();
        for(auto it = mc.sample_times().cbegin(); it != mc.sample_times().cend(); ++it) {
          json["Pass"].push_back(it->first);
          json["Step"].push_back(it->second);
        }
        for(auto it = mc.trajectory().cbegin(); it != mc.trajectory().cend(); ++it) {
          json["DoF"].push_back(*it);
        }
        gz::ogzstream sout((dir.trajectory_json(cond_index).string() + ".gz").c_str());
        _log << "write: " << fs::path(dir.trajectory_json(cond_index).string() + ".gz") << "\n";
        sout << json;
        sout.close();

        // --- Write "occupation_key.json" ------------
        //
        // [["A", "B"],["A" "C"], ... ]

        jsonParser key = jsonParser::array();
        for(int i = 0; i < prim.basis.size(); i++) {
          key.push_back(prim.basis[i].allowed_occupants());
        }
        key.write(dir.occupation_key_json());
        _log << "write: " << dir.occupation_key_json() << "\n";

      }

    }
    catch(...) {
      std::cerr << "ERROR writing observations." << std::endl;
      throw;
    }

  }

  /// \brief For the initial state, write a POSCAR file.
  ///
  /// The current naming convention is 'POSCAR.initial'
  void write_POSCAR_initial(const GrandCanonical &mc, Index cond_index, Log &_log) {

    GrandCanonicalDirectoryStructure dir(mc.settings().output_directory());
    fs::create_directories(dir.trajectory_dir(cond_index));

    // read initial_state.json
    ConfigDoF config_dof;
    from_json(config_dof, jsonParser(dir.initial_state_json(cond_index)));

    if(!fs::exists(dir.initial_state_json(cond_index))) {
      throw std::runtime_error(
        std::string("ERROR in 'write_POSCAR_initial(const GrandCanonical &mc, Index cond_index)'\n") +
        "  File not found: " + dir.initial_state_json(cond_index).string());
    }

    // write file
    fs::ofstream sout(dir.POSCAR_initial(cond_index));
    _log << "write: " << dir.POSCAR_initial(cond_index) << "\n";
    VaspIO::PrintPOSCAR p(mc.supercell(), config_dof);
    p.sort();
    p.print(sout);
    return;
  }

  /// \brief For the final state, write a POSCAR file.
  ///
  /// The current naming convention is 'POSCAR.final'
  void write_POSCAR_final(const GrandCanonical &mc, Index cond_index, Log &_log) {

    GrandCanonicalDirectoryStructure dir(mc.settings().output_directory());
    fs::create_directories(dir.trajectory_dir(cond_index));

    // read final_state.json
    ConfigDoF config_dof;
    from_json(config_dof, jsonParser(dir.final_state_json(cond_index)));

    if(!fs::exists(dir.final_state_json(cond_index))) {
      throw std::runtime_error(
        std::string("ERROR in 'write_POSCAR_final(const GrandCanonical &mc, Index cond_index)'\n") +
        "  File not found: " + dir.final_state_json(cond_index).string());
    }

    // write file
    fs::ofstream sout(dir.POSCAR_final(cond_index));
    _log << "write: " << dir.POSCAR_final(cond_index) << "\n";
    VaspIO::PrintPOSCAR p(mc.supercell(), config_dof);
    p.sort();
    p.print(sout);
    return;
  }

  /// \brief For every snapshot taken, write a POSCAR file.
  ///
  /// The current naming convention is 'POSCAR.sample'
  /// POSCAR title comment is printed with "Sample: #  Pass: #  Step: #"
  void write_POSCAR_trajectory(const GrandCanonical &mc, Index cond_index, Log &_log) {

    GrandCanonicalDirectoryStructure dir(mc.settings().output_directory());
    fs::create_directories(dir.trajectory_dir(cond_index));

    std::vector<Index> pass;
    std::vector<Index> step;
    std::vector<ConfigDoF> trajectory;

    // create super structure matching supercell
    BasicStructure<Site> primstruc = mc.supercell().prim();
    BasicStructure<Site> superstruc = primstruc.create_superstruc(mc.supercell().real_super_lattice());

    if(mc.settings().write_json()) {

      std::string filename = dir.trajectory_json(cond_index).string() + ".gz";

      if(!fs::exists(filename)) {
        throw std::runtime_error(
          std::string("ERROR in 'write_POSCAR_trajectory(const GrandCanonical &mc, Index cond_index)'\n") +
          "  File not found: " + filename);
      }

      gz::igzstream sin(filename.c_str());
      jsonParser json(sin);
      for(auto it = json["Pass"].cbegin(); it != json["Pass"].cend(); ++it) {
        pass.push_back(it->get<Index>());
      }
      for(auto it = json["Step"].cbegin(); it != json["Step"].cend(); ++it) {
        step.push_back(it->get<Index>());
      }
      ConfigDoF config_dof;
      for(auto it = json["DoF"].cbegin(); it != json["DoF"].cend(); ++it) {
        from_json(config_dof, *it);
        trajectory.push_back(config_dof);
      }

    }
    else if(mc.settings().write_csv()) {

      std::string filename = dir.trajectory_csv(cond_index).string() + ".gz";

      if(!fs::exists(filename)) {
        throw std::runtime_error(
          std::string("ERROR in 'write_POSCAR_trajectory(const GrandCanonical &mc, Index cond_index)'\n") +
          "  File not found: " + filename);
      }

      gz::igzstream sin(filename.c_str());

      // skip the header line
      sin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      Index _pass, _step;
      ConfigDoF config_dof(superstruc.basis.size());

      while(sin) {
        sin >> _pass >> _step;
        pass.push_back(_pass);
        step.push_back(_step);
        for(Index i = 0; i < superstruc.basis.size(); i++) {
          sin >> config_dof.occ(i);
        }
        trajectory.push_back(config_dof);
      }
    }

    for(Index i = 0; i < trajectory.size(); i++) {

      // copy occupation
      for(Index j = 0; j < trajectory[i].size(); j++) {
        superstruc.basis[j].set_occ_value(trajectory[i].occ(j));
      }

      // POSCAR title comment is printed with "Sample: #  Pass: #  Step: #"
      std::stringstream ss;
      ss << "Sample: " << i << "  Pass: " << pass[i] << "  Step: " << step[i];
      superstruc.title = ss.str();

      // write file
      fs::ofstream sout(dir.POSCAR_snapshot(cond_index, i));
      _log << "write: " << dir.POSCAR_snapshot(cond_index, i) << "\n";
      VaspIO::PrintPOSCAR p(mc.supercell(), trajectory[i]);
      p.set_title(ss.str());
      p.sort();
      p.print(sout);
      sout.close();
    }

    return;
  }

  /// \brief Create a jsonParser object that can be used as a template.
  ///
  jsonParser example_grand_canonical_settings() {

    jsonParser example_settings;

    // ---- Initialization settings --------------------

    //Monte Carlo cell
    Eigen::Matrix3i tmat = Eigen::Matrix3i::Zero();
    tmat(0, 0) = 10;
    tmat(1, 1) = 10;
    tmat(2, 2) = 10;
    example_settings["supercell"] = tmat;

    //Cluster expanstion to use, project clex name
    example_settings["model"]["clex"] = "formation_energy";

    //Basis set to use, project bset name
    example_settings["model"]["bset"] = "default";

    //Calctype settings to use, project calctype name
    example_settings["model"]["calctype"] = "default";

    //Calctype settings to use, project ref name
    example_settings["model"]["ref"] = "ref";

    //ECI settings name, project eci name
    example_settings["model"]["eci"] = "default";


    // ---- Data settings --------------------

    example_settings["data"]["sample_by"] = "pass";
    example_settings["data"]["sample_period"] = 1;

    example_settings["data"]["max_pass"] = 10000;
    example_settings["data"]["min_pass"] = 0;

    // also allowed:
    //example_settings["data"]["equilibration_passes_first_run"] = 2000;
    //example_settings["data"]["equilibration_passes_each_run"] = 1000;

    //example_settings["data"]["max_step"] = 10000;
    //example_settings["data"]["min_step"] = 0;

    // also allowed:
    //example_settings["data"]["max_sample"] = 10000;
    //example_settings["data"]["min_sample"] = 0;

    // confidence level (default 0.95)
    example_settings["data"]["confidence"] = 0.95;

    example_settings["data"]["measurements"] = jsonParser::array(5);
    for(int i = 0; i < 5; i++) {
      example_settings["data"]["measurements"].push_back(jsonParser::object());
    }
    example_settings["data"]["measurements"][0]["quantity"] = "formation_energy";
    example_settings["data"]["measurements"][1]["quantity"] = "generalized_enthalpy";
    example_settings["data"]["measurements"][2]["quantity"] = "mol_composition";
    example_settings["data"]["measurements"][3]["quantity"] = "formation_energy";
    example_settings["data"]["measurements"][4]["quantity"] = "formation_energy";

    // each of the above measurements can also have a requested precision, in which
    // case the monte carlo calculation will run until that precision is reached (<X> +/- prec),
    // given the confidence level
    // example_settings["data"]["measurements"][0]["precision"] = 0.01

    example_settings["data"]["storage"]["output_directory"] = "./mc_results";
    example_settings["data"]["storage"]["write_observations"] = true;
    example_settings["data"]["storage"]["write_trajectory"] = true;
    example_settings["data"]["storage"]["write_POSCAR_snapshots"] = true;

    // output file format
    //
    // - may be single string or array of strings
    // - csv is the default format if no 'output_format' given)
    // - options are "csv"/"CSV" or "json"/"JSON"
    example_settings["data"]["storage"]["output_format"] = jsonParser::array();
    example_settings["data"]["storage"]["output_format"].push_back("csv"); // or "CSV"
    example_settings["data"]["storage"]["output_format"].push_back("json"); // or "JSON"



    // ---- Driver settings -------------------

    example_settings["driver"]["mode"] = "incremental";

    example_settings["driver"]["initial_conditions"]["param_chem_pot"]["a"] = 2.0;
    example_settings["driver"]["initial_conditions"]["temperature"] = 100.0; // K
    example_settings["driver"]["initial_conditions"]["tolerance"] = 0.001;

    example_settings["driver"]["final_conditions"]["param_chem_pot"]["a"] = 2.0;
    example_settings["driver"]["final_conditions"]["temperature"] = 1000.0; // K
    example_settings["driver"]["final_conditions"]["tolerance"] = 0.001;

    example_settings["driver"]["incremental_conditions"]["param_chem_pot"]["a"] = 0.0;
    example_settings["driver"]["incremental_conditions"]["temperature"] = 10.0; // K
    example_settings["driver"]["incremental_conditions"]["tolerance"] = 0.001;

    //Seed configuration
    example_settings["driver"]["motif"]["configname"] = "SCELV_A_B_C_D_E_F/X";

    return example_settings;
  }

  /// \brief Print single spin flip LTE
  GenericDatumFormatter<double, ConstMonteCarloPtr> GrandCanonicalLTEFormatter(const double &phi_LTE1) {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      return phi_LTE1;
    };
    return GenericDatumFormatter<double, ConstMonteCarloPtr>("phi_LTE", "phi_LTE", evaluator);
  }

  /// \brief Will create new file or append to existing results file the results of the latest run
  void write_lte_results(const MonteSettings &settings, const GrandCanonical &mc, const double &phi_LTE1, const std::string &configname, Log &_log) {
    try {

      fs::create_directories(settings.output_directory());
      GrandCanonicalDirectoryStructure dir(settings.output_directory());
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

