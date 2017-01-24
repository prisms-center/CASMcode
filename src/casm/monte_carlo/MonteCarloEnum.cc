#include "casm/monte_carlo/MonteCarloEnum_impl.hh"
#include "casm/monte_carlo/MonteCarlo.hh"

namespace CASM {

  /// \brief Insert in hall of fame if 'check' passes
  MonteCarloEnum::HallOfFameType::InsertResult MonteCarloEnum::_insert(const Configuration &config) {
    if(insert_canonical()) {
      return m_halloffame->insert(config.in_canonical_supercell().canonical_form());
    }
    else {
      return m_halloffame->insert(config);
    }
  }

  /// \brief Attempt to insert Configuration into enumeration hall of fame
  ///
  /// \returns Pair of iterator pointing to inserted Configuration (or end), and
  ///          bool indicating if insert was successful
  ///
  /// Configurations are only inserted into hall of fame if:
  /// - 'enum_check' returns true
  /// - the Configuration is not already in the config list
  MonteCarloEnum::HallOfFameType::InsertResult MonteCarloEnum::insert(const Configuration &config) {

    bool check = (*m_enum_check)(config);

    if(!check) {
      if(debug()) {
        _log().custom("Config enumeration");
        _log() << "enum check: " << std::boolalpha << check << std::endl;
        _log() << std::endl;

        print_info();
      }

      return HallOfFameType::InsertResult(
               m_halloffame->end(),
               false,
               std::numeric_limits<double>::quiet_NaN(),
               false,
               m_halloffame->end());
    }

    auto res = _insert(config);

    if(debug()) {
      _log().custom("Config enumeration");
      _log() << "enum check: " << std::boolalpha << check << std::endl;
      _log() << "score: " << res.score << std::endl;
      _log() << "insert config in hall of fame: " << std::boolalpha << res.success << std::endl;
      if(!res.success) {
        if(res.excluded) {
          _log() << "already in config list: " << res.excluded_pos->second.name() << std::endl;
        }
        else if(res.pos != m_halloffame->end()) {
          _log() << "already in hall of fame: #" << std::distance(m_halloffame->begin(), res.pos) << std::endl;
        }
        else {
          _log() << "score not good enough" << std::endl;
        }
      }
      _log() << std::endl;

      print_info();
    }
    return res;
  }

  /// \brief const Access the enumeration hall of fame
  const MonteCarloEnum::HallOfFameType &MonteCarloEnum::halloffame() const {
    if(m_halloffame) {
      return *m_halloffame;
    }
    else {
      throw std::runtime_error("Error accessing Monte Carlo HallOfFame: was not initialized");
    }
  }

  namespace {

    /// \brief Until we write one to get primitive Configurations directly from
    /// non-primitive configurations
    ///
    /// returns configname and is_new
    ///
    std::pair<std::string, bool> _import_primitive(PrimClex &primclex, Supercell &scel, Index config_index) {

      BasicStructure<Site> nonprim = scel.superstructure(config_index);

      BasicStructure<Site> prim;
      nonprim.is_primitive(prim);

      int map_opt = ConfigMapper::none;
      double tol = TOL;
      double vol_tol = 0.25;
      double lattice_weight = 0.5;
      ConfigMapper configmapper(primclex, lattice_weight, vol_tol, map_opt, tol);

      std::string imported_name;
      Eigen::Matrix3d cart_op;
      std::vector<Index> best_assignment;
      jsonParser fullrelax_data;
      bool is_new = configmapper.import_structure_occupation(prim,
                                                             imported_name,
                                                             fullrelax_data,
                                                             best_assignment,
                                                             cart_op,
                                                             true);
      return std::make_pair(imported_name, is_new);
    }
  }

  /// \brief Save configurations in the hall of fame to the config list
  void MonteCarloEnum::save_configs() {

    if(!halloffame().size()) {
      _log().write("Enumerated configurations to master config list");
      _log() << "No configurations in hall of fame\n";
      _log() << std::endl;
      return;
    }

    std::vector<Configuration> output;
    m_data.clear();

    Index config_index;
    Supercell::permute_const_iterator permute_it;
    bool is_new, is_prim;
    std::string configname;

    // transform hall of fame configurations so that they fill the canonical
    // equivalent supercell, and add to project
    for(const auto &val : halloffame()) {

      // get equivalent configuration (not necessarily canonical) in the
      // canonical equivalent supercell stored in the primclex
      Configuration config = val.second.in_canonical_supercell();
      Supercell &canon_scel = config.get_supercell();

      // add config to supercell (the saved config will be canonical)
      is_new = canon_scel.add_config(config, config_index, permute_it);
      Configuration &canon_config = canon_scel.get_config(config_index);

      is_prim = is_primitive(canon_config);
      if(is_new) {
        _halloffame().exclude(canon_config);
      }

      double score = val.first;

      // store config source info
      jsonParser json_src;
      std::stringstream ss;
      ss << std::setprecision(6) << score;
      json_src["monte_carlo_enumeration"]["metric"] = m_metric_args;
      json_src["monte_carlo_enumeration"]["score"] = ss.str();
      canon_config.push_back_source(json_src);

      // store info for printing
      m_data[canon_config.name()] = std::make_pair(is_new, score);
      output.push_back(canon_config);

      // if not primitive, generate and import the primitive configuration
      if(!is_prim) {
        std::tie(configname, is_new) = _import_primitive(this->_primclex(), canon_scel, config_index);
        Configuration &prim_canon_config = _primclex().configuration(configname);
        if(is_new) {
          _halloffame().exclude(prim_canon_config);
        }

        // store config source info
        prim_canon_config.push_back_source(json_src);

        // store info for printing
        m_data[prim_canon_config.name()] = std::make_pair(is_new, score);
        output.push_back(prim_canon_config);
      }
    }
    _primclex().write_config_list();

    auto formatter = m_dict.parse("configname is_primitive is_new score comp potential_energy");
    auto flag = FormatFlag(_log()).print_header(true);

    _log().write("Enumerated configurations to master config list");
    _log() << "configuration enumeration check: " << m_check_args << "\n";
    _log() << "configuration enumeration metric: " << m_metric_args << "\n";
    _log() << flag << formatter(output.begin(), output.end());
    _log() << std::endl;

  }

  void MonteCarloEnum::print_info() const {

    _log().custom("Enumerated configurations hall of fame");
    _log() << "configuration enumeration check: " << m_check_args << "\n";
    _log() << "configuration enumeration metric: " << m_metric_args << "\n";
    _log() << std::setw(16) << "position" << std::setw(16) << "score" << "\n";
    _log() << std::setw(16) << std::string("-", 12) << std::setw(16) << std::string("-", 12) << "\n";

    Index i = 0;
    for(const auto &val : halloffame()) {
      _log() << std::setw(16) << i << std::setw(16) << val.first << "\n";
      i++;
    }
    _log() << std::endl;
  }

  /// \brief Clear hall of fame and reset excluded
  void MonteCarloEnum::reset() {
    m_halloffame->clear();
    if(check_existence()) {
      m_halloffame->clear_excluded();
      m_halloffame->exclude(this->primclex().config_begin(), this->primclex().config_end());
    }
  }

  MonteCarloEnum::HallOfFameType &MonteCarloEnum::_halloffame() {
    if(m_halloffame) {
      return *m_halloffame;
    }
    else {
      throw std::runtime_error("Error accessing Monte Carlo HallOfFame: was not initialized");
    }
  }
}

