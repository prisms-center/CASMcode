#ifndef CASM_MonteCarloEnum_impl
#define CASM_MonteCarloEnum_impl

#include "casm/monte_carlo/MonteCarloEnum.hh"
#include "casm/clex/ConfigIterator.hh"
#include "casm/clex/ConfigMapping.hh"

namespace CASM {

  template<typename Base>
  template<typename MonteTypeSettings>
  MonteCarloEnum<Base>::MonteCarloEnum(PrimClex &primclex, const MonteTypeSettings &settings, Log &_log) :
    Base(primclex, settings, _log) {

    if(settings.is_enumeration()) {
      m_check_args = settings.enumeration_check_args();
      m_metric_args = settings.enumeration_metric_args();

      m_halloffame = settings.enumeration_halloffame(this->primclex().settings());
      m_halloffame->exclude(primclex.config_begin(), primclex.config_end());

      m_enum_check = settings.enumeration_check(this->primclex().settings());
    }
  }

  /// \brief Check if enumeration is requested
  template<typename Base>
  bool MonteCarloEnum<Base>::enum_configs() const {
    return static_cast<bool>(m_halloffame);
  }

  /// \brief Attempt to insert (canonical) Configuration into enumeration hall of fame
  ///
  /// \returns Pair of iterator pointing to inserted Configuration (or end), and
  ///          bool indicating if insert was successful
  ///
  /// Configurations are only inserted into hall of fame if:
  /// - 'enum_check' returns true
  /// - the Configuration is not already in the config list
  template<typename Base>
  typename MonteCarloEnum<Base>::HallOfFameType::InsertResult MonteCarloEnum<Base>::enum_insert(const Configuration &config) {

    Log &log = this->_log();
    bool check = (*m_enum_check)(config);

    if(!check) {
      if(this->debug()) {
        log.custom("Config enumeration");
        log << "enum check: " << std::boolalpha << check << std::endl;
        log << std::endl;

        print_info();
      }

      return HallOfFameType::InsertResult(m_halloffame->end(), false, std::numeric_limits<double>::quiet_NaN(), false);
    }

    /*
    const Supercell& scel = config.get_supercell();
    auto permute_begin = scel.permute_begin();
    auto permute_end = scel.permute_end();
    Supercell::permute_const_iterator permute_it;

    Configuration canon_config = config.canonical_form(permute_begin, permute_end, permute_it);

    auto res = m_halloffame->insert(canon_config);
    */

    auto res = m_halloffame->insert(config);

    if(this->debug()) {
      log.custom("Config enumeration");
      log << "enum check: " << std::boolalpha << check << std::endl;
      log << "score: " << res.score << std::endl;
      log << "insert config in hall of fame: " << std::boolalpha << res.success << std::endl;
      if(!res.success) {
        if(res.pos != m_halloffame->end()) {
          log << "already in hall of fame: #" << std::distance(m_halloffame->begin(), res.pos) << std::endl;
        }
        else {
          log << "score not good enough" << std::endl;
        }
      }
      log << std::endl;

      print_info();
    }
    return res;
  }

  /// \brief const Access the enumeration hall of fame
  template<typename Base>
  const typename MonteCarloEnum<Base>::HallOfFameType &MonteCarloEnum<Base>::enum_halloffame() const {
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
  template<typename Base>
  void MonteCarloEnum<Base>::save_configs() const {
    std::string underline = "------------";

    Log &log = this->_log();

    log.write("Enumerated configurations to master config list");
    log << "configuration enumeration check: " << m_check_args << "\n";
    log << "configuration enumeration metric: " << m_metric_args << "\n";
    log << std::setw(24) << "configname"
        << std::setw(16) << "primitive"
        << std::setw(16) << "new"
        << std::setw(16) << "score"
        << "\n";
    log << std::setw(24) << underline
        << std::setw(16) << underline
        << std::setw(16) << underline
        << std::setw(16) << underline
        << "\n";

    // get canonical equivalent supercell
    Index Nscel = this->primclex().get_supercell_list().size();
    Index scel_index = this->_primclex().add_supercell(this->supercell().get_real_super_lattice());
    Supercell &scel = this->_primclex().get_supercell(scel_index);

    // if this is a new supercell for the project, write SCEL file
    if(Nscel != this->primclex().get_supercell_list().size()) {
      log << "generate supercell: " << scel.get_name() << "\n";
      log << "write: SCEL\n";
      this->primclex().print_supercells();
    }

    Index config_index;
    Supercell::permute_const_iterator permute_it;
    bool is_new, is_prim;
    std::string configname;

    // transform hall of fame configurations so that they fill the canonical
    // equivalent supercell, and add
    for(const auto &val : enum_halloffame()) {
      Configuration config = fill_supercell(scel, val.second);

      jsonParser json_src;
      std::stringstream ss;
      ss << std::setprecision(6) << val.first;
      json_src["monte_carlo_enumeration"]["metric"] = m_metric_args;
      json_src["monte_carlo_enumeration"]["score"] = ss.str();

      is_new = scel.add_config(config, config_index, permute_it);
      scel.get_config(config_index).push_back_source(json_src);
      is_prim = is_primitive(config);
      configname = scel.get_config(config_index).name();

      log << std::setw(24) << configname
          << std::setw(16) << (is_prim ? "primitive" : "non-primitive")
          << std::setw(16) << (is_new ? "new" : "existing")
          << std::setw(16) << val.first << std::endl;

      // if not primitive, generate and import the primitive configuration
      if(!is_prim) {
        auto res = _import_primitive(this->_primclex(), scel, config_index);
        configname = res.first;
        is_new = res.second;
        this->_primclex().configuration(configname).push_back_source(json_src);

        log << std::setw(24) << configname
            << std::setw(16) << "primitive"
            << std::setw(16) << (is_new ? "new" : "existing")
            << std::setw(16) << val.first << std::endl;
      }
    }
    this->_primclex().write_config_list();
    log << std::endl;
  }

  template<typename Base>
  void MonteCarloEnum<Base>::print_info() const {
    Log &log = this->_log();

    log.custom("Enumerated configurations hall of fame");
    log << "configuration enumeration check: " << m_check_args << "\n";
    log << "configuration enumeration metric: " << m_metric_args << "\n";
    log << std::setw(16) << "position" << std::setw(16) << "score" << "\n";
    log << std::setw(16) << std::string("-", 12) << std::setw(16) << std::string("-", 12) << "\n";

    Index i = 0;
    for(const auto &val : enum_halloffame()) {
      log << std::setw(16) << i << std::setw(16) << val.first << "\n";
      i++;
    }
    log << std::endl;
  }

  template<typename Base>
  typename MonteCarloEnum<Base>::HallOfFameType &MonteCarloEnum<Base>::_enum_halloffame() {
    if(m_halloffame) {
      return *m_halloffame;
    }
    else {
      throw std::runtime_error("Error accessing Monte Carlo HallOfFame: was not initialized");
    }
  }

}

#endif

