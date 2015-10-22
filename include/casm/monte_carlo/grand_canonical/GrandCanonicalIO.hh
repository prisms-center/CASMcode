#ifndef CASM_GrandCanonicalIO_HH
#define CASM_GrandCanonicalIO_HH

#include <string>
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"
#include "casm/monte_carlo/MonteIO.hh"

namespace CASM {
  
  class GrandCanonicalDirectoryStructure : public MonteCarloDirectoryStructure {
    
    public:
    
    GrandCanonicalDirectoryStructure(fs::path output_dir) :
      MonteCarloDirectoryStructure(fs::absolute(output_dir)) {}
    
    /// \brief "output_dir/occupation_key.csv"
    fs::path occupation_key_csv() const {
      return output_dir() / "occupation_key.csv";
    }
    
    /// \brief "output_dir/occupation_key.csv"
    fs::path occupation_key_json() const {
      return output_dir() / "occupation_key.json";
    }
    
  };
  
  
  template<typename T>
  class DataFormatter;
  
  /// \brief Make a results formatter
  DataFormatter<ConstMonteCarloPtr> make_results_formatter(const GrandCanonical& mc);
  
  /// \brief Make a observation formatter
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_observation_formatter(const GrandCanonical& mc);
  
  /// \brief Make a trajectory formatter
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_trajectory_formatter(const GrandCanonical& mc);
  
  
  /// \brief Store GrandCanonicalConditions in JSON format
  jsonParser &to_json(const GrandCanonicalConditions &conditions, jsonParser &json);  
  
  /// \brief Read GrandCanonicalConditions from JSON format
  void from_json(GrandCanonicalConditions &conditions, const CompositionConverter& comp_converter, const jsonParser &json);  
  
  
  /// \brief Will create new file or append to existing file results of the latest run
  void write_results(const MonteSettings &settings, const GrandCanonical &mc);
  
  /// \brief Write conditions to conditions.cond_index directory
  void write_conditions_json(const MonteSettings &settings, const GrandCanonical &mc, Index cond_index);
  
  /// \brief Will create (and possibly overwrite) new file with all observations from run with conditions.cond_index
  void write_observations(const MonteSettings &settings, const GrandCanonical &mc, Index cond_index);
  
  /// \brief Will create (and possibly overwrite) new file with all observations from run with conditions.cond_index
  void write_trajectory(const MonteSettings &settings, const GrandCanonical &mc, Index cond_index);
  
  /// \brief For every snapshot taken, write a POSCAR file.
  void write_pos_trajectory(const MonteSettings& settings, const GrandCanonical &mc, Index cond_index);
  
  /// \brief Create a jsonParser object that can be used as a template.
  jsonParser example_grand_canonical_settings();
}

#endif
