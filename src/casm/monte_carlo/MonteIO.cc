#include "casm/monte_carlo/MonteIO.hh"
#include "casm/external/gzstream/gzstream.h"
#include "casm/casm_io/VaspIO.hh"
#include "casm/casm_io/DataFormatter.hh"
#include "casm/monte_carlo/MonteCarlo.hh"

namespace CASM {

  /// \brief Print mean property values: <prop_name>
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloMeanFormatter(std::string prop_name) {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->samplers().find(prop_name)->second->mean(mc->is_equilibrated().second);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("<") + prop_name + ">";

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Print calculated precision of property values: prec(<prop_name>)
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloPrecFormatter(std::string prop_name) {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->samplers().find(prop_name)->second->calculated_precision(mc->is_equilibrated().second);
    };

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("prec(<") + prop_name + ">)";

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  double CovEvaluator::operator()(const ConstMonteCarloPtr &mc) {
    auto equil = mc->is_equilibrated();

    // cov = <X*Y> - <X>*<Y>
    const MonteSampler &sampler1 = *(mc->samplers().find(prop_name1)->second);
    const Eigen::VectorXd &obs1 = sampler1.data().observations();
    const Eigen::VectorXd &X = obs1.segment(equil.second, obs1.size() - equil.second);

    const MonteSampler &sampler2 = *(mc->samplers().find(prop_name2)->second);
    const Eigen::VectorXd &obs2 = sampler2.data().observations();
    const Eigen::VectorXd &Y = obs2.segment(equil.second, obs2.size() - equil.second);

    double Xsum = 0.0;
    double Ysum = 0.0;
    double XYsum = 0.0;
    Index N = X.size();

    for(Index i = 0; i < N; ++i) {
      Xsum += X(i);
      Ysum += Y(i);
      XYsum += X(i) * Y(i);
    }

    //return X.cwiseProduct(Y).mean() - X.mean()*Y.mean();
    return (XYsum - Xsum * Ysum / N) / N;
  }

  /// \brief Print covariance: cov(prop_name1, prop_name2)
  GenericDatumFormatter<double, ConstMonteCarloPtr> MonteCarloCovFormatter(std::string prop_name1, std::string prop_name2) {

    auto evaluator = CovEvaluator(prop_name1, prop_name2);

    auto validator = [ = ](const ConstMonteCarloPtr & mc) {
      return mc->is_equilibrated().first;
    };

    std::string header = std::string("cov(") + prop_name1 + "," + prop_name2 + ")";

    return GenericDatumFormatter<double, ConstMonteCarloPtr>(header, header, evaluator, validator);
  }

  /// \brief Print if equilibrated (not counting explicitly requested equilibration)
  GenericDatumFormatter<bool, ConstMonteCarloPtr> MonteCarloIsEquilibratedFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc)->bool {
      return mc->is_equilibrated().first;
    };

    return GenericDatumFormatter<bool, ConstMonteCarloPtr>("is_equilibrated", "is_equilibrated", evaluator);
  }

  /// \brief Print if converged
  GenericDatumFormatter<bool, ConstMonteCarloPtr> MonteCarloIsConvergedFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc)->bool {
      return mc->is_converged();
    };

    return GenericDatumFormatter<bool, ConstMonteCarloPtr>("is_converged", "is_converged", evaluator);
  }

  /// \brief Print number of samples used for equilibration (not counting explicitly requested equilibration)
  GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr> MonteCarloNEquilSamplesFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc)->unsigned long {
      return mc->is_equilibrated().second;
    };

    return GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr>("N_equil_samples", "N_equil_samples", evaluator);
  }

  /// \brief Print number of samples used in calculating means
  GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr> MonteCarloNAvgSamplesFormatter() {

    auto evaluator = [ = ](const ConstMonteCarloPtr & mc)->unsigned long {
      return mc->sample_times().size() - mc->is_equilibrated().second;
    };

    return GenericDatumFormatter<MonteSampler::size_type, ConstMonteCarloPtr>("N_avg_samples", "N_avg_samples", evaluator);
  }

  /// \brief Print Pass number of observation
  GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> > MonteCarloPassFormatter() {

    auto evaluator = [ = ](const std::pair<ConstMonteCarloPtr, Index> &obs) {
      return obs.first->sample_times()[obs.second].first;
    };

    return GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> >("Pass", "Pass", evaluator);
  }

  /// \brief Print Step number of observation
  GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> > MonteCarloStepFormatter() {

    auto evaluator = [ = ](const std::pair<ConstMonteCarloPtr, Index> &obs) {
      return obs.first->sample_times()[obs.second].second;
    };

    return GenericDatumFormatter<MonteCounter::size_type, std::pair<ConstMonteCarloPtr, Index> >("Step", "Step", evaluator);
  }

  /// \brief Print value of observation
  GenericDatumFormatter<double, std::pair<ConstMonteCarloPtr, Index> > MonteCarloObservationFormatter(std::string prop_name) {

    auto evaluator = [ = ](const std::pair<ConstMonteCarloPtr, Index> &obs) {
      const MonteSampler &sampler = *(obs.first->samplers().find(prop_name)->second);
      return sampler.data().observations()(obs.second);
    };

    return GenericDatumFormatter<double, std::pair<ConstMonteCarloPtr, Index> >(prop_name, prop_name, evaluator);
  }

  /// \brief Print value of a particular occupation variable
  GenericDatumFormatter<int, std::pair<ConstMonteCarloPtr, Index> > MonteCarloOccFormatter(Index occ_index) {

    auto evaluator = [ = ](const std::pair<ConstMonteCarloPtr, Index> &site)->int {
      return site.first->trajectory()[site.second].occ(occ_index);
    };

    std::string header = std::string("occ(") + std::to_string(occ_index) + ")";

    return GenericDatumFormatter<int, std::pair<ConstMonteCarloPtr, Index> >(header, header, evaluator);
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
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_observation_formatter(const MonteCarlo &mc) {

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
  DataFormatter<std::pair<ConstMonteCarloPtr, Index> > make_trajectory_formatter(const MonteCarlo &mc) {
    DataFormatter<std::pair<ConstMonteCarloPtr, Index> > formatter;
    formatter.push_back(MonteCarloPassFormatter());
    formatter.push_back(MonteCarloStepFormatter());

    // this is probably not the best way...
    for(Index i = 0; i < mc.configdof().occupation().size(); ++i) {
      formatter.push_back(MonteCarloOccFormatter(i));
    }
    return formatter;
  }

  /// \brief Will create (and possibly overwrite) new file with all observations from run with conditions.cond_index
  void write_observations(const MonteSettings &settings, const MonteCarlo &mc, Index cond_index, Log &_log) {
    try {
      if(!settings.write_observations()) {
        return;
      }

      MonteCarloDirectoryStructure dir(settings.output_directory());
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
  void write_trajectory(const MonteSettings &settings, const MonteCarlo &mc, Index cond_index, Log &_log) {
    try {

      if(!settings.write_trajectory()) {
        return;
      }

      MonteCarloDirectoryStructure dir(settings.output_directory());
      fs::create_directories(dir.conditions_dir(cond_index));
      auto formatter = make_trajectory_formatter(mc);
      const Structure &prim = mc.primclex().get_prim();

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
  void write_POSCAR_initial(const MonteCarlo &mc, Index cond_index, Log &_log) {

    MonteCarloDirectoryStructure dir(mc.settings().output_directory());
    fs::create_directories(dir.trajectory_dir(cond_index));

    // read initial_state.json
    ConfigDoF config_dof;
    from_json(config_dof, jsonParser(dir.initial_state_json(cond_index)));

    if(!fs::exists(dir.initial_state_json(cond_index))) {
      throw std::runtime_error(
        std::string("ERROR in 'write_POSCAR_initial(const MonteCarlo &mc, Index cond_index)'\n") +
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
  void write_POSCAR_final(const MonteCarlo &mc, Index cond_index, Log &_log) {

    MonteCarloDirectoryStructure dir(mc.settings().output_directory());
    fs::create_directories(dir.trajectory_dir(cond_index));

    // read final_state.json
    ConfigDoF config_dof;
    from_json(config_dof, jsonParser(dir.final_state_json(cond_index)));

    if(!fs::exists(dir.final_state_json(cond_index))) {
      throw std::runtime_error(
        std::string("ERROR in 'write_POSCAR_final(const MonteCarlo &mc, Index cond_index)'\n") +
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
  void write_POSCAR_trajectory(const MonteCarlo &mc, Index cond_index, Log &_log) {

    MonteCarloDirectoryStructure dir(mc.settings().output_directory());
    fs::create_directories(dir.trajectory_dir(cond_index));

    std::vector<Index> pass;
    std::vector<Index> step;
    std::vector<ConfigDoF> trajectory;

    // create super structure matching supercell
    BasicStructure<Site> primstruc = mc.supercell().get_prim();
    BasicStructure<Site> superstruc = primstruc.create_superstruc(mc.supercell().get_real_super_lattice());

    if(mc.settings().write_json()) {

      std::string filename = dir.trajectory_json(cond_index).string() + ".gz";

      if(!fs::exists(filename)) {
        throw std::runtime_error(
          std::string("ERROR in 'write_POSCAR_trajectory(const MonteCarlo &mc, Index cond_index)'\n") +
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
          std::string("ERROR in 'write_POSCAR_trajectory(const MonteCarlo &mc, Index cond_index)'\n") +
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


}
