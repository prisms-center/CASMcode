#include "casm/monte_carlo/grand_canonical/GrandCanonicalIO.hh"

#include "casm/external/gzstream/gzstream.h"
#include "casm/CASM_global_definitions.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"

namespace CASM {

  /// \brief Make a results formatter
  ///
  /// For csv:
  /// \code
  /// # T Beta mu_a mu_b ... N_equil_samples N_avg_samples <X> prec(<X>) ... ... cov(X,X) ... 
  /// \endcode
  ///
  /// For JSON:
  /// \code
  /// {
  ///   "T":[...], "Beta":[...], "mu_a":[...], "mu_b":[...], ..., 
  ///   "N_equil_samples":[...], "N_avg_samples":[...], 
  ///   "<X>":[...], "prec(<X>)":[...], ..., "cov(X,X)":[...], ... }
  /// \endcode
  ///
  DataFormatter<ConstMonteCarloPtr> make_results_formatter(const GrandCanonical& mc) {
  
  DataFormatter<ConstMonteCarloPtr> formatter;
  
  formatter.push_back(MonteCarloTFormatter<GrandCanonical>());
  formatter.push_back(MonteCarloBetaFormatter<GrandCanonical>());
  for(int i=0; i<mc.primclex().composition_axes().independent_compositions(); ++i) {
    formatter.push_back(MonteCarloMuFormatter<GrandCanonical>(mc, i));
  }
  
  formatter.push_back(MonteCarloNEquilSamplesFormatter());
  formatter.push_back(MonteCarloNAvgSamplesFormatter());
  
  for(auto it=mc.samplers().cbegin(); it != mc.samplers().cend(); ++it) {
    formatter.push_back(MonteCarloMeanFormatter(it->first));
    formatter.push_back(MonteCarloPrecFormatter(it->first));
  }
  
  for(auto it_X=mc.samplers().cbegin(); it_X != mc.samplers().cend(); ++it_X) {
    for(auto it_Y=mc.samplers().cbegin(); it_Y != mc.samplers().cend(); ++it_Y) {
      formatter.push_back(MonteCarloCovFormatter(it_X->first, it_Y->first));
    }
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
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_observation_formatter(const GrandCanonical& mc) {
    
    DataFormatter<std::pair<ConstMonteCarloPtr, Index> > formatter;
    
    formatter.push_back(MonteCarloPassFormatter());
    formatter.push_back(MonteCarloStepFormatter());
    for(auto it=mc.samplers().cbegin(); it != mc.samplers().cend(); ++it) {
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
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_trajectory_formatter(const GrandCanonical& mc) {
    DataFormatter<std::pair<ConstMonteCarloPtr, Index> > formatter;
    formatter.push_back(MonteCarloPassFormatter());
    formatter.push_back(MonteCarloStepFormatter());
    
    // this is probably not the best way...
    for(Index i=0; i<mc.configdof().occupation().size(); ++i) {
      formatter.push_back(MonteCarloOccFormatter(i));
    }
    return formatter;
  }
  
  /// \brief Store GrandCanonicalConditions in JSON format
  jsonParser &to_json(const GrandCanonicalConditions &conditions, jsonParser &json) {

    json["temperature"] = conditions.temperature();
    json["mu"] = jsonParser::object();
    auto param_mu = conditions.param_mu();
    for(int i=0; i<param_mu.size(); i++) {
      json["mu"][std::string(1, (char) (i + (int) 'a'))] = param_mu(i);
    }
    json["tolerance"] = conditions.tolerance();

    return json;
  }
  
  /// \brief Will create new file or append to existing results file the results of the latest run
  void write_results(const MonteSettings &settings, const GrandCanonical &mc) {
    try {
      
      fs::create_directories(settings.output_directory());
      GrandCanonicalDirectoryStructure dir(settings.output_directory());
      auto formatter = make_results_formatter(mc);
      
      // write csv path results
      if(settings.write_csv()) {
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
  void write_conditions_json(const MonteSettings &settings, const GrandCanonical &mc, Index cond_index) {
    try {
      GrandCanonicalDirectoryStructure dir(settings.output_directory());
      fs::create_directories(dir.conditions_dir(cond_index));
      jsonParser json;
      to_json(mc.conditions(), json);
      json.write(dir.conditions_json(cond_index)); 
    } 
    catch(...) {
      std::cerr << "ERROR writing conditions.json" << std::endl;
      throw;
    }
  }
  
  /// \brief Will create (and possibly overwrite) new file with all observations from run with conditions.cond_index
  void write_observations(const MonteSettings &settings, const GrandCanonical &mc, Index cond_index) {
    try {
      if(!settings.write_observations()) {
        return;
      }
      
      GrandCanonicalDirectoryStructure dir(settings.output_directory());
      fs::create_directories(dir.conditions_dir(cond_index));
      auto formatter = make_observation_formatter(mc);
      
      std::vector<std::pair<ConstMonteCarloPtr, Index> > observations;
      ConstMonteCarloPtr ptr = &mc;
      for(MonteSampler::size_type i=0; i<mc.sample_times().size(); ++i) {
        observations.push_back(std::make_pair(ptr, i));
      }
      
      if(settings.write_csv()) {
        gz::ogzstream sout((dir.observations_csv(cond_index).string() + ".gz").c_str());
        sout << formatter(observations.cbegin(), observations.cend());
        sout.close();
      }
      
      if(settings.write_json()) {
        gz::ogzstream sout((dir.observations_json(cond_index).string() + ".gz").c_str());
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
  void write_trajectory(const MonteSettings &settings, const GrandCanonical &mc, Index cond_index) {
    try {
      
      if(!settings.write_trajectory()) {
        return;
      }
      
      GrandCanonicalDirectoryStructure dir(settings.output_directory());
      fs::create_directories(dir.conditions_dir(cond_index));
      auto formatter = make_trajectory_formatter(mc);
      const Structure &prim = mc.primclex().get_prim();
        
      std::vector<std::pair<ConstMonteCarloPtr, Index> > observations;
      ConstMonteCarloPtr ptr = &mc;
      for(MonteSampler::size_type i=0; i<mc.sample_times().size(); ++i) {
        observations.push_back(std::make_pair(ptr, i));
      }
      
      if(settings.write_csv()) {
        gz::ogzstream sout((dir.trajectory_csv(cond_index).string() + ".gz").c_str());
        sout << formatter(observations.cbegin(), observations.cend());
        sout.close();
        
        int max_allowed = 0;
        for(int i=0; i<prim.basis.size(); i++) {
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
        keyout << "site_index";
        for(int i=0; i<max_allowed; i++) {
          keyout << "\tocc_index_" << i;
        }
        keyout << "\n";
        
        for(int i=0; i<prim.basis.size(); i++) {
          keyout << i;
          for(int j=0; j<max_allowed; j++) {
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
        formatter(observations.cbegin(), observations.cend()).to_json_arrays(json);
        gz::ogzstream sout((dir.trajectory_json(cond_index).string() + ".gz").c_str());
        sout << json;
        sout.close();
        
        // --- Write "occupation_key.json" ------------
        //
        // [["A", "B"],["A" "C"], ... ]
        
        jsonParser key = jsonParser::array();
        for(int i=0; i<prim.basis.size(); i++) {
          key.push_back(prim.basis[i].allowed_occupants());
        }
        key.write(dir.occupation_key_json());
        
      }
      
    }
    catch(...) {
      std::cerr << "ERROR writing observations." << std::endl;
      throw;
    }
    
  }
  
  /// \brief For every snapshot taken, write a POSCAR file.
  /// 
  /// The current naming convention is 'POSCAR.sample'
  /// POSCAR title comment is printed with "Sample: #  Pass: #  Step: #"
  void write_pos_trajectory(const MonteSettings& settings, const GrandCanonical &mc, Index cond_index) {
    
    if(!settings.write_POSCAR_snapshots() || !settings.write_trajectory() || mc.trajectory().size() == 0) {
      return;
    }
    
    GrandCanonicalDirectoryStructure dir(settings.output_directory());
    fs::create_directories(dir.trajectory_dir(cond_index));
      
    // create super structure matching supercell
    BasicStructure<Site> primstruc = mc.supercell().get_prim();
    BasicStructure<Site> superstruc = primstruc.create_superstruc(mc.supercell().get_real_super_lattice());
    
    for(Index i = 0; i < mc.trajectory().size(); i++) {
      
      // copy occupation
      for(Index j = 0; j < mc.configdof().size(); j++) {
        superstruc.basis[j].set_occ_value(mc.configdof().occ(j));
      }
      
      // POSCAR title comment is printed with "Sample: #  Pass: #  Step: #" 
      std::stringstream ss;
      ss << "Sample: " << i << "  Pass: " << mc.sample_times()[i].first << "  Step: " << mc.sample_times()[i].second;
      superstruc.title = ss.str();
      
      // write file
      fs::ofstream sout(dir.POSCAR_snapshot(cond_index, i));
      superstruc.print5(sout);
      sout.close();
    }

    return;
  }

}

