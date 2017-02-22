#ifndef CASM_Enumerator_impl
#define CASM_Enumerator_impl

#include "casm/system/RuntimeLibrary.hh"
#include "casm/container/Enumerator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

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
  template<typename ScelIterator, typename ConfigEnumConstructor>
  int insert_unique_canon_configs(
    std::string method,
    PrimClex &primclex,
    ScelIterator begin,
    ScelIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr) {

    Log &log = primclex.log();

    Index Ninit = primclex.db<Configuration>().size();
    log << "# configurations in this project: " << Ninit << "\n" << std::endl;

    log.begin(method);

    for(auto scel_it = begin; scel_it != end; ++scel_it) {
      Supercell &scel = *scel_it;
      log << "Enumerate configurations for " << scel.name() << " ...  " << std::flush;

      auto enumerator_ptr = f(scel);
      auto &enumerator = *enumerator_ptr;
      Index num_before = primclex.db<Configuration>().scel_range(scel.name()).size();
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

      log << (primclex.db<Configuration>().scel_range(scel.name()).size() - num_before) << " configs." << std::endl;
    }
    log << "  DONE." << std::endl << std::endl;

    Index Nfinal = primclex.db<Configuration>().size();

    log << "# new configurations: " << Nfinal - Ninit << "\n";
    log << "# configurations in this project: " << Nfinal << "\n" << std::endl;

    log << "Write supercell database..." << std::endl;
    primclex.db<Supercell>().commit();
    log << "  DONE" << std::endl << std::endl;

    log << "Writing configuration database..." << std::endl;
    primclex.db<Configuration>().commit();
    log << "  DONE" << std::endl;

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
  int insert_configs(
    std::string method,
    PrimClex &primclex,
    LatticeIterator begin,
    LatticeIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr,
    bool primitive_only) {

    Log &log = primclex.log();

    Index Ninit = primclex.db<Configuration>().size();
    log << "# configurations in this project: " << Ninit << "\n" << std::endl;

    log.begin(method);

    for(auto scel_lat_it = begin; scel_lat_it != end; ++scel_lat_it) {
      Supercell scel(&primclex, *scel_lat_it);
      Supercell &canon_scel = scel.canonical_form();
      log << "Enumerate configurations for " << canon_scel.name() << " ...  " << std::flush;

      auto enumerator_ptr = f(scel);
      auto &enumerator = *enumerator_ptr;
      Index num_before = primclex.db<Configuration>().scel_range(canon_scel.name()).size();
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

      log << (primclex.db<Configuration>().scel_range(canon_scel.name()).size() - num_before) << " configs." << std::endl;
    }
    log << "  DONE." << std::endl << std::endl;

    Index Nfinal = primclex.db<Configuration>().size();

    log << "# new configurations: " << Nfinal - Ninit << "\n";
    log << "# configurations in this project: " << Nfinal << "\n" << std::endl;

    log << "Write supercell database..." << std::endl;
    primclex.db<Supercell>().commit();
    log << "  DONE" << std::endl << std::endl;

    log << "Writing configuration database..." << std::endl;
    primclex.db<Configuration>().commit();
    log << "  DONE" << std::endl;

    return 0;
  }

}

#endif
