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
    
    /// \brief Results summary: "output_dir/lte_results.csv"
    fs::path lte_results_csv() const {
      return m_output_dir / "lte_results.csv";
    }
    
    /// \brief Results summary: "output_dir/lte_results.json" 
    fs::path lte_results_json() const {
      return m_output_dir / "lte_results.json";
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

    /// \brief "output_dir/conditions.cond_index/trajectory/POSCAR.final"
    fs::path POSCAR_final(int cond_index) {
      return trajectory_dir(cond_index) / std::string("POSCAR.final");
    }

    /// \brief "output_dir/conditions.cond_index/trajectory/POSCAR.sample"
    fs::path POSCAR_snapshot(int cond_index, MonteSampler::size_type sample_index) {
      return trajectory_dir(cond_index) / (std::string("POSCAR.") + std::to_string(sample_index));
    }

    /// \brief "output_dir/conditions.cond_index/final_state.json"
    fs::path final_state_json(int cond_index) const {
      return conditions_dir(cond_index) / "final_state.json";
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

  /// \brief Print covariance: cov(prop_name1, prop_name2)
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloCovFormatter(std::string prop_name1, std::string prop_name2);

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



}

#endif
