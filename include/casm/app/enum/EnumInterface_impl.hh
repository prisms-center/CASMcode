#ifndef CASM_EnumInterface_impl
#define CASM_EnumInterface_impl

#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/enum/EnumInterface.hh"
#include "casm/clex/FilteredConfigIterator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/Enumerator.hh"
#include "casm/system/RuntimeLibrary.hh"

namespace CASM {

  // might be useful for other casm commands...
  template<typename OptionType>
  jsonParser make_json_input(const OptionType &opt) {
    if(!opt.settings_path().empty()) {
      return jsonParser {opt.settings_path()};
    }
    else if(!opt.input_str().empty()) {
      return jsonParser::parse(opt.input_str());
    }
    else {
      return jsonParser::object();
    }
  }

  /// \brief Standardizes insertion from enumerators that construct unique
  /// primitive canonical configurations
  ///
  /// \param method Enumerator name, printed to screen
  /// \param primclex PrimClex to add Configurations to
  /// \param begin,end Iterators over canonical Supercell
  /// \param f A function that signature `std::unique_ptr<EnumMethod> f(Supercell&)`
  ///        that returns a std::unique_ptr owning a Configuration enumerator
  /// \param filter_expr An vector of Configuration filtering expressions. No filtering if empty.
  ///
  /// \returns 0 if success, ERR_INVALID_ARG if filter_expr cannot be parsed
  template<typename InConfigIterator, typename ConfigEnumConstructor>
  int insert_unique_canon_configs(
    std::string method,
    const PrimClex &primclex,
    InConfigIterator begin,
    InConfigIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool dry_run) {
    for(; begin != end; ++begin) {
      if(0 != insert_unique_canon_configs(method,
                                          primclex,
                                          ConfigEnumInput(*begin),
                                          f,
                                          filter_expr,
                                          dry_run)) {
        return ERR_INVALID_ARG;
      }
    }
    return 0;
  }

  /// \brief Standardizes insertion from enumerators that construct unique
  /// primitive canonical configurations
  ///
  /// \param method Enumerator name, printed to screen
  /// \param primclex PrimClex to add Configurations to
  /// \param begin,end Iterators over canonical Supercell
  /// \param f A function that signature `std::unique_ptr<EnumMethod> f(Supercell&)`
  ///        that returns a std::unique_ptr owning a Configuration enumerator
  /// \param filter_expr An vector of Configuration filtering expressions. No filtering if empty.
  ///
  /// \returns 0 if success, ERR_INVALID_ARG if filter_expr cannot be parsed
  template<typename ConfigEnumConstructor>
  int insert_unique_canon_configs(
    std::string method,
    const PrimClex &primclex,
    ConfigEnumInput const &in_config,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool dry_run) {

    Log &log = primclex.log();
    auto &db_config = primclex.db<Configuration>();
    auto const &configuration = in_config.configuration();
    auto const &supercell = configuration.supercell();

    std::string dry_run_msg = CASM::dry_run_msg(dry_run);

    Index Ninit = db_config.size();
    log << dry_run_msg << "# configurations in this project: " << Ninit << "\n" << std::endl;

    log.begin(method);


    auto enumerator_ptr = f(in_config);
    auto &enumerator = *enumerator_ptr;
    log << dry_run_msg << "Enumerate configurations for " << configuration.name() << " ...  " << std::flush;

    Index num_before = db_config.scel_range_size(supercell.name());
    if(!filter_expr.empty()) {
      try {
        primclex.db<Configuration>().insert(
          filter_begin(
            enumerator.begin(),
            enumerator.end(),
            filter_expr,
            primclex.settings().query_handler<Configuration>().dict()),
          filter_end(enumerator.end()));
      }
      catch(std::exception &e) {
        primclex.err_log() << "Cannot filter configurations using the expression provided: \n" << e.what() << "\nExiting...\n";
        return ERR_INVALID_ARG;
      }
    }
    else {
      primclex.db<Configuration>().insert(enumerator.begin(), enumerator.end());
    }

    log << (db_config.scel_range_size(supercell.name()) - num_before) << " configs." << std::endl;
    log << dry_run_msg << "  DONE." << std::endl << std::endl;

    Index Nfinal = db_config.size();

    log << dry_run_msg << "# new configurations: " << Nfinal - Ninit << "\n";
    log << dry_run_msg << "# configurations in this project: " << Nfinal << "\n" << std::endl;

    if(!dry_run) {
      log << "Write supercell database..." << std::endl;
      primclex.db<Supercell>().commit();
      log << "  DONE" << std::endl << std::endl;

      log << "Writing configuration database..." << std::endl;
      db_config.commit();
      log << "  DONE" << std::endl;
    }

    return 0;
  }

  /// \brief Standardizes insertion from enumerators that construct configurations
  ///
  /// \param method Enumerator name, printed to screen
  /// \param primclex PrimClex to add Configurations to
  /// \param begin,end Iterators over Supercell
  /// \param f A function that signature `std::unique_ptr<EnumMethod> f(Supercell&)`
  ///        that returns a std::unique_ptr owning a Configuration enumerator
  /// \param filter_expr An vector of Configuration filtering expressions. No filtering if empty.
  /// \param primitive_only If true, only insert primitive Configurations; else
  ///        both primitive and non-primitive.
  ///
  /// \returns 0 if success, ERR_INVALID_ARG if filter_expr cannot be parsed
  ///
  template<typename InConfigIterator, typename ConfigEnumConstructor>
  int insert_configs(
    std::string method,
    const PrimClex &primclex,
    InConfigIterator begin,
    InConfigIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool primitive_only,
    bool dry_run) {
    for(; begin != end; ++begin) {
      if(0 != insert_configs(method,
                             primclex,
                             ConfigEnumInput(*begin),
                             f,
                             filter_expr,
                             primitive_only,
                             dry_run)) {
        return ERR_INVALID_ARG;
      }
    }
    return 0;
  }
  /// \brief Standardizes insertion from enumerators that construct configurations
  ///
  /// \param method Enumerator name, printed to screen
  /// \param primclex PrimClex to add Configurations to
  /// \param begin,end Iterators over Supercell
  /// \param f A function that signature `std::unique_ptr<EnumMethod> f(Supercell&)`
  ///        that returns a std::unique_ptr owning a Configuration enumerator
  /// \param filter_expr An vector of Configuration filtering expressions. No filtering if empty.
  /// \param primitive_only If true, only insert primitive Configurations; else
  ///        both primitive and non-primitive.
  ///
  /// \returns 0 if success, ERR_INVALID_ARG if filter_expr cannot be parsed
  ///
  template<typename ConfigEnumConstructor>
  int insert_configs(
    std::string method,
    const PrimClex &primclex,
    ConfigEnumInput const &in_config,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool primitive_only,
    bool dry_run) {

    Log &log = primclex.log();
    auto &db_config = primclex.db<Configuration>();
    auto const &configuration = in_config.configuration();
    auto const &supercell = configuration.supercell();

    std::string dry_run_msg = CASM::dry_run_msg(dry_run);

    Index Ninit = db_config.size();
    log << dry_run_msg << "# configurations in this project: " << Ninit << "\n" << std::endl;

    log.begin(method);

    log << dry_run_msg << "Enumerate configurations for " << configuration.name() << " ...  " << std::flush;

    auto enumerator_ptr = f(in_config);
    auto &enumerator = *enumerator_ptr;
    Index num_before = db_config.scel_range_size(supercell.name());
    if(!filter_expr.empty()) {
      try {
        auto it = filter_begin(
                    enumerator.begin(),
                    enumerator.end(),
                    filter_expr,
                    primclex.settings().query_handler<Configuration>().dict());
        auto end = filter_end(enumerator.end());
        for(; it != end; ++it) {
          it->insert(primitive_only);
        }
      }
      catch(std::exception &e) {
        primclex.err_log() << "Cannot filter configurations using the expression provided: \n" << e.what() << "\nExiting...\n";
        return ERR_INVALID_ARG;
      }
    }
    else {
      //std::cout << "BEFORE ITER\n";
      auto it = enumerator.begin();
      auto end = enumerator.end();
      for(; it != end; ++it) {
        //std::cout << "ITER: " << it->configdof().global_dof("GLstrain").values().transpose() << "\n";
        //auto result = it->insert(primitive_only);
        it->insert(primitive_only);
        //std::cout << "Iter insert...\n";
        //std::cout << "Iter insert result, primitive as " << result.primitive_it->name() << "\n";
        //std::cout << "Iter insert result, primitive as " << result.canonical_it->name() << "\n";

      }
    }

    log << (db_config.scel_range_size(supercell.name()) - num_before) << " configs." << std::endl;

    log << dry_run_msg << "  DONE." << std::endl << std::endl;

    Index Nfinal = db_config.size();

    log << dry_run_msg << "# new configurations: " << Nfinal - Ninit << "\n";
    log << dry_run_msg << "# configurations in this project: " << Nfinal << "\n" << std::endl;

    if(!dry_run) {
      log << "Write supercell database..." << std::endl;
      primclex.db<Supercell>().commit();
      log << "  DONE" << std::endl << std::endl;

      log << "Writing configuration database..." << std::endl;
      db_config.commit();
      log << "  DONE" << std::endl;
    }

    return 0;
  }

  /// \brief Standardizes insertion from enumerators that construct configurations
  ///
  /// \param method Enumerator name, printed to screen
  /// \param primclex PrimClex to add Configurations to
  /// \param begin,end Iterators over lattice (need not be canonical), to be used for Supercell
  /// \param f A function that signature `std::unique_ptr<EnumMethod> f(Supercell&)`
  ///        that returns a std::unique_ptr owning a Configuration enumerator
  /// \param filter_expr An vector of Configuration filtering expressions. No filtering if empty.
  /// \param primitive_only If true, only insert primitive Configurations; else
  ///        both primitive and non-primitive.
  ///
  /// \returns 0 if success, ERR_INVALID_ARG if filter_expr cannot be parsed
  ///
  template<typename LatticeIterator, typename ConfigEnumConstructor>
  int insert_configs_via_lattice_enum(
    std::string method,
    const PrimClex &primclex,
    LatticeIterator begin,
    LatticeIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool primitive_only,
    bool dry_run) {

    Log &log = primclex.log();
    auto &db_config = primclex.db<Configuration>();
    auto distance = [&](const std::string & scelname) {
      return db_config.scel_range_size(scelname);
    };

    std::string dry_run_msg = CASM::dry_run_msg(dry_run);

    Index Ninit = db_config.size();
    log << dry_run_msg << "# configurations in this project: " << Ninit << "\n" << std::endl;

    log.begin(method);

    for(auto scel_lat_it = begin; scel_lat_it != end; ++scel_lat_it) {
      Supercell scel(&primclex, *scel_lat_it);
      const Supercell &canon_scel = scel.canonical_form();
      log << dry_run_msg << "Enumerate configurations for " << canon_scel.name() << " ...  " << std::flush;

      auto enumerator_ptr = f(ConfigEnumInput(Configuration::zeros(scel)));
      auto &enumerator = *enumerator_ptr;
      Index num_before = distance(canon_scel.name());
      if(!filter_expr.empty()) {
        try {
          auto it = filter_begin(
                      enumerator.begin(),
                      enumerator.end(),
                      filter_expr,
                      primclex.settings().query_handler<Configuration>().dict());
          auto end = filter_end(enumerator.end());
          for(; it != end; ++it) {
            it->insert(primitive_only);
          }
        }
        catch(std::exception &e) {
          primclex.err_log() << "Cannot filter configurations using the expression provided: \n" << e.what() << "\nExiting...\n";
          return ERR_INVALID_ARG;
        }
      }
      else {
        auto it = enumerator.begin();
        auto end = enumerator.end();
        for(; it != end; ++it) {
          it->insert(primitive_only);
        }
      }

      log << (distance(canon_scel.name()) - num_before) << " configs." << std::endl;
    }
    log << dry_run_msg << "  DONE." << std::endl << std::endl;

    Index Nfinal = primclex.db<Configuration>().size();

    log << dry_run_msg << "# new configurations: " << Nfinal - Ninit << "\n";
    log << dry_run_msg << "# configurations in this project: " << Nfinal << "\n" << std::endl;

    if(!dry_run) {
      log << "Write supercell database..." << std::endl;
      primclex.db<Supercell>().commit();
      log << "  DONE" << std::endl << std::endl;

      log << "Writing configuration database..." << std::endl;
      db_config.commit();
      log << "  DONE" << std::endl;
    }

    return 0;
  }

}

#endif
