#ifndef CASM_MonteIO_HH
#define CASM_MonteIO_HH

#include <string>
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/monte_carlo/MonteCarlo.hh"

namespace CASM {

  class MonteCarloDirectoryStructure {

  public:

    MonteCarloDirectoryStructure(fs::path output_dir) :
      m_output_dir(fs::absolute(output_dir)) {}

    /// \brief "output_dir/"
    fs::path output_dir() const {
      return m_output_dir;
    }

    /// \brief Results summary: "output_dir/results.csv"
    fs::path results_csv() const {
      return m_output_dir / "results.csv";
    }

    /// \brief Results summary: "output_dir/results.json"
    fs::path results_json() const {
      return m_output_dir / "results.json";
    }


    /// \brief "output_dir/conditions.cond_index/"
    fs::path conditions_dir(int cond_index) const {
      return m_output_dir / (std::string("conditions.") + std::to_string(cond_index));
    }

    /// \brief "output_dir/conditions.cond_index/conditions.json"
    fs::path conditions_json(int cond_index) const {
      return conditions_dir(cond_index) / "conditions.json";
    }

    /// \brief "output_dir/conditions.cond_index/observations.csv"
    fs::path observations_csv(int cond_index) const {
      return conditions_dir(cond_index) / "observations.csv";
    }

    /// \brief "output_dir/conditions.cond_index/observations.json"
    fs::path observations_json(int cond_index) const {
      return conditions_dir(cond_index) / "observations.json";
    }


    /// \brief "output_dir/conditions.cond_index/trajectory.csv"
    fs::path trajectory_csv(int cond_index) const {
      return conditions_dir(cond_index) / "trajectory.csv";
    }

    /// \brief "output_dir/conditions.cond_index/trajectory.json"
    fs::path trajectory_json(int cond_index) const {
      return conditions_dir(cond_index) / "trajectory.json";
    }

    /// \brief "output_dir/conditions.cond_index/trajectory"
    fs::path trajectory_dir(int cond_index) const {
      return conditions_dir(cond_index) / "trajectory";
    }

    /// \brief "output_dir/conditions.cond_index/trajectory/POSCAR.initial"
    fs::path POSCAR_initial(int cond_index) {
      return trajectory_dir(cond_index) / std::string("POSCAR.initial");
    }

    /// \brief "output_dir/conditions.cond_index/trajectory/POSCAR.final"
    fs::path POSCAR_final(int cond_index) {
      return trajectory_dir(cond_index) / std::string("POSCAR.final");
    }

    /// \brief "output_dir/conditions.cond_index/trajectory/POSCAR.sample"
    fs::path POSCAR_snapshot(int cond_index, MonteSampler::size_type sample_index) {
      return trajectory_dir(cond_index) / (std::string("POSCAR.") + std::to_string(sample_index));
    }

    /// \brief "output_dir/conditions.cond_index/initial_state_firstruneq.json"
    ///
    /// - Initial state before 'first run equilibration'
    fs::path initial_state_firstruneq_json(int cond_index) const {
      return conditions_dir(cond_index) / "initial_state_firstruneq.json";
    }

    /// \brief "output_dir/conditions.cond_index/final_state.json"
    ///
    /// - Initial state before 'each run equilibration'
    fs::path initial_state_runeq_json(int cond_index) const {
      return conditions_dir(cond_index) / "initial_state_runeq.json";
    }

    /// \brief "output_dir/conditions.cond_index/initial_state.json"
    ///
    /// - Initial state before first pass / step
    fs::path initial_state_json(int cond_index) const {
      return conditions_dir(cond_index) / "initial_state.json";
    }

    /// \brief "output_dir/conditions.cond_index/final_state.json"
    fs::path final_state_json(int cond_index) const {
      return conditions_dir(cond_index) / "final_state.json";
    }

    /// \brief "output_dir/occupation_key.csv"
    fs::path occupation_key_csv() const {
      return output_dir() / "occupation_key.csv";
    }

    /// \brief "output_dir/occupation_key.csv"
    fs::path occupation_key_json() const {
      return output_dir() / "occupation_key.json";
    }

  private:

    fs::path m_output_dir;

  };


  /// \brief const pointer to const MonteCarlo
  typedef const MonteCarlo *ConstMonteCarloPtr;

  /// \brief Print mean property values: <prop_name>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloMeanFormatter(std::string prop_name);

  /// \brief Print calculated precision of property values: prec(<prop_name>)
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloPrecFormatter(std::string prop_name);

  /// \brief Functor to help evaluate covariance
  struct CovEvaluator {

    CovEvaluator(std::string _prop_name1, std::string _prop_name2):
      prop_name1(_prop_name1), prop_name2(_prop_name2) {}

    double operator()(const ConstMonteCarloPtr &mc);

    std::string prop_name1;
    std::string prop_name2;

  };

  /// \brief Print covariance: cov(prop_name1, prop_name2)
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloCovFormatter(std::string prop_name1, std::string prop_name2);

  /// \brief Print if equilibrated (not counting explicitly requested equilibration)
  GenericDatumFormatter<bool, ConstMonteCarloPtr> MonteCarloIsEquilibratedFormatter();

  /// \brief Print if converged
  GenericDatumFormatter<bool, ConstMonteCarloPtr> MonteCarloIsConvergedFormatter();

  /// \brief Print Temperature
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloTFormatter();

  /// \brief Print Beta
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloBetaFormatter();

  /// \brief Print param_chem_pot(x)
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloParamChemPotFormatter(const MonteType &mc, int index);

  /// \brief Print chem_pot(N)
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloChemPotFormatter(const MonteType &mc, int index);

  /// \brief Print comp(x)
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloCompFormatter(const MonteType &mc, int index);

  /// \brief Print comp_n(N)
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloCompNFormatter(const MonteType &mc, int index);

  /// \brief Print heat capacity, 'heat_capacity'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloHeatCapacityFormatter();

  /// \brief Print parametric susceptibility, 'susc_x(a,b)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr>
  MonteCarloSuscXFormatter(std::string comp_var_i, std::string comp_var_j);

  /// \brief Print susceptibility, 'susc_n(A,B)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr>
  MonteCarloSuscNFormatter(std::string species_i, std::string species_j);

  /// \brief Print parametric thermo-chemical susceptibility, 'susc_x(S,a)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr>
  MonteCarloThermoChemSuscXFormatter(std::string comp_var_i);

  /// \brief Print thermo-chemical susceptibility, 'susc_n(S,A)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr>
  MonteCarloThermoChemSuscNFormatter(std::string species_i);

  /// \brief Print number of samples used for equilibration (not counting explicitly requested equilibration)
  GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr> MonteCarloNEquilSamplesFormatter();

  /// \brief Print number of samples used in calculating means
  GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr> MonteCarloNAvgSamplesFormatter();

  /// \brief Print Pass number of observation
  GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> > MonteCarloPassFormatter();

  /// \brief Print Step number of observation
  GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> > MonteCarloStepFormatter();

  /// \brief Print value of observation
  GenericDatumFormatter<double, std::pair<ConstMonteCarloPtr, Index> > MonteCarloObservationFormatter(std::string prop_name);

  /// \brief Print value of a particular occupation variable
  GenericDatumFormatter<int, std::pair<ConstMonteCarloPtr, Index> > MonteCarloOccFormatter(Index occ_index);

  /// \brief Make a observation formatter
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_observation_formatter(const MonteCarlo &mc);

  /// \brief Make a trajectory formatter
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_trajectory_formatter(const MonteCarlo &mc);

  /// \brief Will create new file or append to existing file results of the latest run
  template<typename MonteType>
  void write_results(const MonteSettings &settings, const MonteType &mc, Log &_log);

  /// \brief Write conditions to conditions.cond_index directory
  template<typename MonteType>
  void write_conditions_json(const MonteSettings &settings, const MonteType &mc, Index cond_index, Log &_log);

  /// \brief Will create (and possibly overwrite) new file with all observations from run with conditions.cond_index
  void write_observations(const MonteSettings &settings, const MonteCarlo &mc, Index cond_index, Log &_log);

  /// \brief Will create (and possibly overwrite) new file with all observations from run with conditions.cond_index
  void write_trajectory(const MonteSettings &settings, const MonteCarlo &mc, Index cond_index, Log &_log);

  /// \brief For the initial state, write a POSCAR file.
  void write_POSCAR_initial(const MonteCarlo &mc, Index cond_index, Log &_log);

  /// \brief For the final state, write a POSCAR file.
  void write_POSCAR_final(const MonteCarlo &mc, Index cond_index, Log &_log);

  /// \brief For every snapshot taken, write a POSCAR file.
  void write_POSCAR_trajectory(const MonteCarlo &mc, Index cond_index, Log &_log);


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
    typedef const MonteType *ConstMonteTypePtr;
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      return static_cast<const MonteType *>(ptr)->conditions().beta();
    };
    return GenericDatumFormatter<double, ConstMonteCarloPtr>("Beta", "Beta", evaluator);
  }

  /// \brief Print param_chem_pot(x) for any class MonteType with valid 'double MonteType::conditions().param_chem_pot()(index)'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloParamChemPotFormatter(const MonteType &mc, int index) {
    typedef const MonteType *ConstMonteTypePtr;
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
    typedef const MonteType *ConstMonteTypePtr;
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
    typedef const MonteType *ConstMonteTypePtr;
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
    typedef const MonteType *ConstMonteTypePtr;
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      return static_cast<const MonteType *>(ptr)->conditions().mol_composition()(index);
    };
    std::string header = std::string("comp_n(") + CompositionConverter::comp_var(index) + ")";
    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator);
  }

  /// \brief Print site_frac(x) for any class MonteType with valid 'Eigen::VectorXd MonteType::comp_n'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloSiteFracFormatter(const MonteType &mc, int index) {
    typedef const MonteType *ConstMonteTypePtr;
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      Eigen::VectorXd comp_n = static_cast<const MonteType *>(ptr)->comp_n();
      Eigen::VectorXd site_frac = comp_n / mc.primclex().get_prim().basis.size();
      return site_frac(index);
    };
    std::string header = std::string("site_frac(") + mc.primclex().get_prim().get_struc_molecule_name()[index] + ")";
    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator);
  }

  /// \brief Print atom_frac(X) for any class MonteType with valid 'Eigen::VectorXd MonteType::comp_n'
  template<typename MonteType>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloAtomFracFormatter(const MonteType &mc, int index) {
    typedef const MonteType *ConstMonteTypePtr;
    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      ConstMonteCarloPtr ptr = mc;
      Eigen::VectorXd comp_n = static_cast<const MonteType *>(ptr)->comp_n();
      if(mc.primclex().vacancy_allowed()) {
        comp_n(mc.primclex().vacancy_index()) = 0.0;
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
