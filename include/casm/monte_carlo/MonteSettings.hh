#ifndef CASM_MonteSettings_HH
#define CASM_MonteSettings_HH

#include <string>
#include "casm/CASM_global_definitions.hh"
#include "casm/misc/cloneable_ptr.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/clex/ECIContainer.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/monte_carlo/MonteDefinitions.hh"
#include "casm/monte_carlo/MonteSampler.hh"

namespace CASM {

  class MonteCarlo;
  
  /*
   * MonteSettings is nothing more than a jsonParser that you can expect
   * to contain a set of values. It's a class meant to hold all the
   * members of the base MonteCarlo class so you can write out all the
   * relevant data and then also use it to reconstruct MonteCarlo.
   *
   * The point of the class is mostly to ensure that when you reconstruct
   * your MonteCarlo you're not giving it some willy nilly jsonParser
   * you came up with. It MUST be something that MonteCarlo itself
   * spit out.
   *
   * Inheritance is private, so you can't use [] as an access operator.
   * Instead use the prescribed accessing routines. This will keep the
   * structure of the settings file consistent and changes will
   * only have to occur through this class.
   *
   * To avoid having millions of access routines, there are jsonParser
   * pointers internally that point to a "current state" of subtrees
   * that may involve the same access routines. For example, instead
   * of having different access routines for initial conditions and
   * final conditions, you just set the current conditions mode, and
   * then the same access routines can be used. From the outside, it
   * looks like you're just requesting a subset of your MonteSettings routine.
   */

  /// \brief Settings for Monte Carlo calculations
  class MonteSettings: private jsonParser {

  public:
    
    using jsonParser::print;
    using jsonParser::write;
    using jsonParser::operator==;
    using jsonParser::operator!=;

    typedef Index size_type;

    /// \brief Default constructor
    MonteSettings() {}
    
    /// \brief Construct MonteSettings by reading a settings JSON file
    MonteSettings(const fs::path &read_path);
    
    
    // --- Type ---------------------------
    
    /// \brief Return type of Monte Carlo calculation
    Monte::TYPE type() const;
    
    
    // --- Initialization ---------------------
    
    /// \brief Configname of configuration to use as starting motif
    std::string motif_configname() const;
    
    /// \brief Supercell matrix defining the simulation cell
    Eigen::Matrix3i simulation_cell_matrix() const;
    

    // --- Driver ---------------------
    
    /// \brief Given a settings jsonParser figure out the drive mode. Expects drive_mode/single,custom,incremental
    const Monte::DRIVE_MODE drive_mode() const;
    
    
    // --- Conditions settings ---------------------
    
    /// \brief Expects initial_conditions
    const MonteSettings &initial_conditions() const;
    
    /// \brief Expects final_conditions
    const MonteSettings &final_conditions() const;
    
    /// \brief Expects incremental_conditions
    const MonteSettings &incremental_conditions() const;

    /// \brief Given a settings jsonParser figure out the chemical potential for a particular component. Expects mu/'a'/value
    double chemical_potential(std::string comp_name) const;
    
    /// \brief Given a settings jsonParser figure out the temperature of the run. Expects temperature/value
    double temperature() const;
    
    /// \brief Given a settings jsonParser figure out the global tolerance (probably for == operator). Expects tolerance/value
    double tolerance() const;
    
    
    
    
    // --- Project settings ---------------------
    
    /// \brief Given a settings jsonParser figure out what the project clex settings to use are:
    std::string clex() const;

    /// \brief Given a settings jsonParser figure out what the project bset settings to use are:
    std::string bset() const;

    /// \brief Given a settings jsonParser figure out what the project calctype settings to use are:
    std::string calctype() const;

    /// \brief Given a settings jsonParser figure out what the project ref settings to use are:
    std::string ref() const;

    /// \brief Given a settings jsonParser figure out what the project eci settings to use are:
    std::string eci() const;

    /// \brief Directory where output should go
    const fs::path output_directory() const;


    // --- MCData / Sampling ---------------------

    /// \brief Sample by pass?
    bool sample_by_pass() const;
    
    /// \brief Sample by step?
    bool sample_by_step() const;

    /// \brief Figure out how often to take samples
    size_type sample_period() const;
    
    
    /// \brief Requested confidence level. Default 0.95.
    double confidence() const;

    
    /// \brief Returns true if explicit equilibration passes for the first run have been specified
    bool is_equilibration_passes_first_run() const;
    
    /// \brief Number of explicit equilibration passes requsted for the first run
    size_type equilibration_passes_first_run() const;
    
    /// \brief Returns true if explicit equilibration passes for each run have been specified
    bool is_equilibration_passes_each_run() const;
    
    /// \brief Number of explicit equilibration passes requsted for each run
    size_type equilibration_passes_each_run() const;
    
    
    /// \brief Returns true if (*this)[level1].contains(level2)
    bool _is_setting(std::string level1, std::string level2) const;
    
    /// \brief Returns (*this)[level1][level2].get<T>();
    template<typename T>
    T _get_setting(std::string level1, std::string level2) const;
    
    
    /// \brief Returns true if the number of passes has been specified
    bool is_N_pass() const;
    
    /// \brief Returns the number of passes requested
    size_type N_pass() const;
    
    /// \brief Returns true if the number of steps has been specified
    bool is_N_step() const;
  
    /// \brief Returns the number of steps requested
    size_type N_step() const;
    
    /// \brief Returns true if the number of samples has been specified
    bool is_N_sample() const;
  
    /// \brief Returns the number of samples requested
    size_type N_sample() const;
  
    
    
    /// \brief Returns true if a maximum number of passes has been specified
    bool is_max_pass() const;
  
    /// \brief Maximum number of passes, required if sample by pass
    size_type max_pass() const;
    
    /// \brief Returns true if a minimum number of passes has been specified
    bool is_min_pass() const;
    
    /// \brief Minimum number of passes, default 0 if sample by pass
    size_type min_pass() const;

    
    /// \brief Returns true if a maximum number of steps has been specified
    bool is_max_step() const;
    
    /// \brief Maximum number of steps, required if sample by step
    size_type max_step() const;
    
    /// \brief Returns true if a minimum number of steps has been specified
    bool is_min_step() const;
    
    /// \brief Minimum number of steps, default 0 if sample by step
    size_type min_step() const;

    
    /// \brief Returns true if a maximum number of samples has been specified
    bool is_max_sample() const;
    
    /// \brief Maximum number of steps, default std::numeric_limit<size_type>::max()
    size_type max_sample() const;
    
    /// \brief Returns true if a minimum number of samples has been specified
    bool is_min_sample() const;
    
    /// \brief Minimum number of steps, default 0
    size_type min_sample() const;
    
    
    /// \brief Construct MonteSamplers as specified in the MonteSettings
    template<typename SamplerInsertIterator>
    SamplerInsertIterator samplers(const PrimClex &primclex, SamplerInsertIterator result) const;
    

    /// \brief Returns true if snapshots are requested
    bool write_trajectory() const;
    
    /// \brief Returns true if POSCARs of snapshots are requsted. Requires write_trajectory.
    bool write_POSCAR_snapshots() const;
    
    /// \brief Writes all observations
    bool write_observations() const;
    
    /// \brief Write csv versions of files? (csv is the default format if no 'output_format' given)
    bool write_csv() const;
    
    /// \brief Write json versions of files?
    bool write_json() const;
    

    // --- Data ---------------------
    
    /// \brief Figure out how large data containers should be
    size_type max_data_length() const;
    
    
    private:


    /// \brief Pointer into conditions subtree. Any access routine you use involving conditions will grab whatever this is pointing to.
    mutable const jsonParser *m_pcondition_subtree;
    
    /// \brief Pointer into monte carlo subtree. Any access involving Monte Carlo states will grab whatever this is pointing to.
    mutable jsonParser *m_pmontestate_subtree;

    
    
    template<typename jsonParserIteratorType>
    std::tuple<bool, double> _get_precision(jsonParserIteratorType it, std::string input_name) const;
    
    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_comp_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;
    
    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_comp_n_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;
    
    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_site_frac_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;
    
    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_atom_frac_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;
    
    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_all_correlations_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;
    
    template<typename jsonParserIteratorType, typename SamplerInsertIterator>
    SamplerInsertIterator _make_non_zero_eci_correlations_samplers(const PrimClex &primclex, jsonParserIteratorType it, SamplerInsertIterator result) const;
  
  };
  
  inline bool operator==(const jsonParser &json, const MonteSettings &settings) {
    return settings == json;
  }
  
  inline bool operator!=(const jsonParser &json, const MonteSettings &settings) {
    return settings != json;
  }


  /// \brief Return an example json object with dummy values
  jsonParser example_base_json_settings();
  
  
  /// \brief Returns (*this)[level1][level2].get<T>();
  template<typename T>
  T MonteSettings::_get_setting(std::string level1, std::string level2) const {
    try {
      return (*this)[level1][level2].get<T>();
    }

    catch(std::runtime_error &e) {
      T t;
      std::cerr << "ERROR in MonteSettings::" << level2 << std::endl;
      std::cerr << "Expected [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
      std::cerr << "  Either this was not found, or the type is wrong." << std::endl;
      if(this->contains(level1)) {
        std::cerr << "Found Settings[\"" << level1 << "\"], but not [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        std::cerr << "Settings[\"" << level1 << "\"]:\n" << (*this)[level1] << std::endl;
      }
      else {
        std::cerr << "No Settings[\"" << level1 << "\"] found" << std::endl;
        std::cerr << "Settings:\n" << (*this) << std::endl;
      }
      throw;
    }
  }
  
  /// \brief Construct MonteSamplers as specified in the MonteSettings
  ///
  /// The requested MonteSamplers are inserted in 'result' as 
  ///   std::pair<std::string, notstd::cloneable_ptr<MonteSampler> >
  /// 
  template<typename SamplerInsertIterator>
  SamplerInsertIterator MonteSettings::samplers(const PrimClex &primclex, SamplerInsertIterator result) const {
    
    size_type data_maxlength = max_data_length();
    std::string prop_name;
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler*ptr;
    
    
    // construct scalar property samplers
    {
      std::vector<std::string> possible = {
        "formation_energy",
        "potential_energy"};
        
      
      std::string level1 = "data";
      std::string level2 = "measurements";
      try {
        const jsonParser& t_measurements = (*this)[level1][level2];
        for(auto it = t_measurements.cbegin(); it != t_measurements.cend(); it++) {
          prop_name = (*it)["quantity"].get<std::string>();
          
          // check if property found is in list of possible scalar properties
          if(std::find(possible.cbegin(), possible.cend(), prop_name) != possible.cend()) {
            
            std::tie(must_converge, prec) = _get_precision(it, prop_name);
            
            // if 'must converge'
            if(must_converge) {
              ptr = new ScalarMonteSampler(prop_name, prop_name, prec, confidence(), data_maxlength);
            }
            else {
              ptr = new ScalarMonteSampler(prop_name, prop_name, confidence(), data_maxlength);
            }
            
            *result++ = std::make_pair(prop_name, notstd::cloneable_ptr<MonteSampler>(ptr));
            
          }
          
        }
        
      }
      catch(std::runtime_error &e) {
        std::cerr << "ERROR in 'MonteSettings::samplers(const PrimClex &primclex, SamplerInsertIterator result)'\n" << std::endl;
        std::cerr << "Error reading [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        throw;
      }
    }
    
    
    // construct vector properties to sample
    {
      std::vector<std::string> possible = {
        "comp",
        "comp_n",
        "site_frac",
        "atom_frac",
        "all_correlations",
        "non_zero_eci_correlations"};
        
      
      std::string level1 = "data";
      std::string level2 = "measurements";
      try {
        
        const jsonParser& t_measurements = (*this)[level1][level2];
        for(auto it = t_measurements.cbegin(); it != t_measurements.cend(); it++) {
          
          std::string input_name = (*it)["quantity"].get<std::string>();
          
          // check if property found is in list of possible vector properties
          if(std::find(possible.cbegin(), possible.cend(), input_name) != possible.cend()) {
            
            // construct MonteSamplers for 'comp'
            if(input_name == "comp") {
              
               result = _make_comp_samplers(primclex, it, result);
                
            }
            
            // construct MonteSamplers for 'comp_n'
            else if(input_name == "comp_n") {
              
              result = _make_comp_n_samplers(primclex, it, result);
              
            }
            
            // construct MonteSamplers for 'site_frac'
            else if(input_name == "site_frac") {
              
              result = _make_site_frac_samplers(primclex, it, result);
              
            }
            
            // construct MonteSamplers for 'atom_frac'
            else if(input_name == "atom_frac") {
              
              result = _make_atom_frac_samplers(primclex, it, result);
              
            }
            
            // construct MonteSamplers for 'all_correlations'
            else if(input_name == "all_correlations") {
              
              result = _make_all_correlations_samplers(primclex, it, result);
              
            }
            
            // construct MonteSamplers for 'non_zero_eci_correlations'
            else if(input_name == "non_zero_eci_correlations") {
              
              result = _make_non_zero_eci_correlations_samplers(primclex, it, result);
              
            }
          }
        }
      }
      catch(std::runtime_error &e) {
        std::cerr << "ERROR in 'MonteSettings::samplers(const PrimClex &primclex, SamplerInsertIterator result)'\n" << std::endl;
        std::cerr << "Error reading [\"" << level1 << "\"][\"" << level2 << "\"]" << std::endl;
        throw;
      }
    }
    
    return result;
  }
  
  template<typename jsonParserIteratorType>
  std::tuple<bool, double> MonteSettings::_get_precision(jsonParserIteratorType it, std::string input_name) const {
    if( it->contains("precision")) {
      return std::make_tuple(true, (*it)["precision"]. template get<double>());
    }
    else {
      return std::make_tuple(false, 0.0);
    }
  }
  
  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator MonteSettings::_make_comp_samplers(const PrimClex &primclex, 
                                          jsonParserIteratorType it,
                                          SamplerInsertIterator result) const {
    
    size_type data_maxlength = max_data_length();
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler*ptr;
    
    for(size_type i=0; i<primclex.composition_axes().independent_compositions(); i++) {
      
      print_name = std::string("comp(") + std::string(1, (char) (i + ((int) 'a'))) + ")";
      
      std::tie(must_converge, prec) = _get_precision(it, "comp");
  
      // if 'must converge'
      if(must_converge) {
        ptr = new CompMonteSampler(i, primclex.composition_axes(), print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new CompMonteSampler(i, primclex.composition_axes(), print_name, confidence(), data_maxlength);
      }
      
      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));
      
    }
    
    return result;
  }
  
  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator MonteSettings::_make_comp_n_samplers(const PrimClex &primclex, 
                                            jsonParserIteratorType it,
                                            SamplerInsertIterator result) const {
    
    size_type data_maxlength = max_data_length();
    std::string prop_name;
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler*ptr;
    
    for(size_type i=0; i<primclex.composition_axes().components().size(); i++) {
      
      prop_name = "comp_n";
      
      print_name = std::string("comp_n(") + primclex.composition_axes().components()[i] + ")";
      
      std::tie(must_converge, prec) = _get_precision(it, "comp_n");

      // if 'must converge'
      if(must_converge) {
        ptr = new VectorMonteSampler(prop_name, i, print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new VectorMonteSampler(prop_name, i, print_name, confidence(), data_maxlength);
      }
      
      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));
      
    }
    
    return result;
  }
  
  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator MonteSettings::_make_site_frac_samplers(const PrimClex &primclex, 
                                               jsonParserIteratorType it,
                                               SamplerInsertIterator result) const {
    
    size_type data_maxlength = max_data_length();
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler*ptr;
                                               
    for(size_type i=0; i<primclex.composition_axes().components().size(); i++) {
                
      // SiteFracMonteSampler uses 'comp_n' to calculate 'site_frac'
      print_name = std::string("site_frac(") + primclex.composition_axes().components()[i] + ")";
      
      std::tie(must_converge, prec) = _get_precision(it, "site_frac");
  
      // if 'must converge'
      if(must_converge) {
        ptr = new SiteFracMonteSampler(i, primclex.get_prim().basis.size(), print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new SiteFracMonteSampler(i, primclex.get_prim().basis.size(), print_name, confidence(), data_maxlength);
      }
      
      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));
      
    }
    
    return result;
  }
  
  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator MonteSettings::_make_atom_frac_samplers(const PrimClex &primclex, 
                                               jsonParserIteratorType it,
                                               SamplerInsertIterator result) const {
  
    size_type data_maxlength = max_data_length();
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler*ptr;
    
    // find the index for vacancies, if they are allowed, 
    //  if not set to primclex.composition_axes().components().size()
    
    size_type vacancy_index = primclex.composition_axes().components().size();
    for(size_type i=0; i<primclex.composition_axes().components().size(); i++) {
      
      // sample for non-vacancy components
      if(Specie(primclex.composition_axes().components()[i]).is_vacancy()) {
        vacancy_index = i;
        break;
      }
    }
    
    
    for(size_type i=0; i<primclex.composition_axes().components().size(); i++) {
      
      // sample for non-vacancy components
      if(!Specie(primclex.composition_axes().components()[i]).is_vacancy()) {
        
        print_name = std::string("atom_frac(") + primclex.composition_axes().components()[i] + ")";
        
        std::tie(must_converge, prec) = _get_precision(it, "atom_frac");
  
        // if 'must converge'
        if(must_converge) {
          ptr = new AtomFracMonteSampler(i, vacancy_index, print_name, prec, confidence(), data_maxlength);
        }
        else {
          ptr = new AtomFracMonteSampler(i, vacancy_index, print_name, confidence(), data_maxlength);
        }
        
        *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));
        
      }
    }
    
    return result;
  }
  
  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator MonteSettings::_make_all_correlations_samplers(const PrimClex &primclex, 
                                                      jsonParserIteratorType it,
                                                      SamplerInsertIterator result) const {
    
    size_type data_maxlength = max_data_length();
    std::string prop_name;
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler*ptr;
    
    for(size_type i=0; i<primclex.global_clexulator().corr_size(); i++) {
      
      prop_name = "corr";
      print_name = std::string("corr(") + std::to_string(i) + ")";
      
      std::tie(must_converge, prec) = _get_precision(it, "all_correlations");

      // if 'must converge'
      if(must_converge) {
        ptr = new VectorMonteSampler(prop_name, i, print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new VectorMonteSampler(prop_name, i, print_name, confidence(), data_maxlength);
      }
      
      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));
      
    }
    
    return result;
  }
  
  template<typename jsonParserIteratorType, typename SamplerInsertIterator>
  SamplerInsertIterator MonteSettings::_make_non_zero_eci_correlations_samplers(const PrimClex &primclex, 
                                                               jsonParserIteratorType it,
                                                               SamplerInsertIterator result) const {
    
    size_type data_maxlength = max_data_length();
    std::string prop_name;
    std::string print_name;
    bool must_converge;
    double prec;
    MonteSampler*ptr;
    
    ECIContainer _eci = primclex.dir().eci_out(clex(), bset(), calctype(), ref(), eci());
    
    for(size_type ii=0; ii<_eci.eci_index_list().size(); ii++) {
      
      prop_name = "corr";
      
      // store non-zero eci index in i
      size_type i = _eci.eci_index_list()[ii];
      
      print_name = std::string("corr(") + std::to_string(i) + ")";
      
      std::tie(must_converge, prec) = _get_precision(it, "non_zero_eci_correlations");
  
      // if 'must converge'
      if(must_converge) {
        ptr = new VectorMonteSampler(prop_name, i, print_name, prec, confidence(), data_maxlength);
      }
      else {
        ptr = new VectorMonteSampler(prop_name, i, print_name, confidence(), data_maxlength);
      }
      
      *result++ = std::make_pair(print_name, notstd::cloneable_ptr<MonteSampler>(ptr));
      
    }
    
    return result;
  }

  
}

#endif
