#ifndef CASM_MonteIO_impl_HH
#define CASM_MonteIO_impl_HH

#include "casm/monte_carlo/MonteIO.hh"

namespace CASM {

  // --- Template definitions ---------------------

  /// \brief Print Temperature for any class MonteType with valid 'double MonteType::conditions().temperature()'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloTFormatter() {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      return static_cast<const MonteType *>(ptr)->conditions().temperature();
    };
    return GenericDatumFormatter<double, ConstMonteCarloPtr>("T", "Temperature", evaluator);
  }

  /// \brief Print Beta for any class MonteType with valid 'double MonteType::conditions().beta()'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloBetaFormatter() {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      return static_cast<const MonteType *>(ptr)->conditions().beta();
    };
    return GenericDatumFormatter<double, ConstMonteCarloPtr>("Beta", "Beta", evaluator);
  }

  /// \brief Print param_chem_pot(x) for any class MonteType with valid 'double MonteType::conditions().param_chem_pot()(index)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloParamChemPotFormatter(const MonteType &mc, int index) {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      return static_cast<const MonteType *>(ptr)->conditions().param_chem_pot()(index);
    };
    std::string header = std::string("param_chem_pot(") + CompositionConverter::comp_var(index) + ")";
    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator);
  }

  /// \brief Print chem_pot(N) for any class MonteType with valid 'double MonteType::conditions().chem_pot(index)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloChemPotFormatter(const MonteType &mc, int index) {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      return static_cast<const MonteType *>(ptr)->conditions().chem_pot(index);
    };
    std::string header = std::string("chem_pot(") + mc.primclex().get_prim().get_struc_molecule_name()[index] + ")";
    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator);
  }

  /// \brief Print comp(x) for any class MonteType with valid 'double MonteType::conditions().param_composition()(index)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloCompFormatter(const MonteType &mc, int index) {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      return static_cast<const MonteType *>(ptr)->conditions().param_composition()(index);
    };
    std::string header = std::string("comp(") + CompositionConverter::comp_var(index) + ")";
    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator);
  }

  /// \brief Print comp(x) for any class MonteType with valid 'double MonteType::conditions().mol_composition()(index)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloCompNFormatter(const MonteType &mc, int index) {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      return static_cast<const MonteType *>(ptr)->conditions().mol_composition()(index);
    };
    std::string header = std::string("comp_n(") + mc.primclex().get_prim().get_struc_molecule_name()[index] + ")";
    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator);
  }

  /// \brief Print site_frac(x) for any class MonteType with valid 'Eigen::VectorXd MonteType::comp_n'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloSiteFracFormatter(const MonteType &mc, int index) {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      const MonteType *mc_ptr = static_cast<const MonteType *>(mc);
      Eigen::VectorXd comp_n = mc_ptr->comp_n();
      Eigen::VectorXd site_frac = comp_n / mc->primclex().get_prim().basis.size();
      return site_frac(index);
    };
    std::string header = std::string("site_frac(") + mc.primclex().get_prim().get_struc_molecule_name()[index] + ")";
    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator);
  }

  /// \brief Print atom_frac(X) for any class MonteType with valid 'Eigen::VectorXd MonteType::comp_n'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloAtomFracFormatter(const MonteType &mc, int index) {
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      const MonteType *mc_ptr = static_cast<const MonteType *>(mc);
      Eigen::VectorXd comp_n = mc_ptr->comp_n();
      if(mc->primclex().vacancy_allowed()) {
        comp_n(mc->primclex().vacancy_index()) = 0.0;
      }
      auto atom_frac = comp_n / comp_n.sum();
      return atom_frac[index];
    };
    std::string header = std::string("atom_frac(") + mc.primclex().get_prim().get_struc_molecule_name()[index] + ")";
    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator);
  }

  /// \brief Print heat capacity, 'heat_capacity'
  ///
  /// - Heat capacity (per unit cell) = cov(potential_energy, potential_energy)*N/(k*T*T)
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloHeatCapacityFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator("potential_energy", "potential_energy");
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const MonteType *>(ptr)->conditions().temperature();
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
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloSuscXFormatter(std::string comp_var_i, std::string comp_var_j) {

    std::string name_i = "comp(" + comp_var_i + ")";
    std::string name_j = "comp(" + comp_var_j + ")";

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator(name_i, name_j);
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const MonteType *>(ptr)->conditions().temperature();
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
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloSuscNFormatter(std::string species_i, std::string species_j) {

    std::string name_i = "comp_n(" + species_i + ")";
    std::string name_j = "comp_n(" + species_j + ")";

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator(name_i, name_j);
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const MonteType *>(ptr)->conditions().temperature();
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
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr>
  MonteCarloThermoChemSuscXFormatter(std::string comp_var_i) {

    std::string name_i = "comp(" + comp_var_i + ")";

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator("potential_energy", name_i);
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const MonteType *>(ptr)->conditions().temperature();
      return cov_evaluator(mc) * N / (KB * T);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("susc_x(S," + comp_var_i + ")");

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Print thermo-chemical susceptibility, 'susc_n(S,A)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr>
  MonteCarloThermoChemSuscNFormatter(std::string species_i) {

    std::string name_i = "comp_n(" + species_i + ")";

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      CovEvaluator cov_evaluator("potential_energy", name_i);
      auto N = mc->supercell().volume();
      ConstMonteCarloPtr ptr = mc;
      auto T = static_cast<const MonteType *>(ptr)->conditions().temperature();
      return cov_evaluator(mc) * N / (KB * T);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("susc_n(S," + species_i + ")");

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Will create new file or append to existing results file the results of the latest run
  template<typename MonteType>
  void write_results(const MonteSettings &settings, const MonteType &mc, Log &_log) {
    try {

      fs::create_directories(settings.output_directory());
      MonteCarloDirectoryStructure dir(settings.output_directory());
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
  template<typename MonteType>
  void write_conditions_json(const MonteSettings &settings, const MonteType &mc, Index cond_index, Log &_log) {
    try {
      MonteCarloDirectoryStructure dir(settings.output_directory());
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

}

#endif
