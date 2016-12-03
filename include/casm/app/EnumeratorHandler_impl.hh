#include "casm/app/EnumeratorHandler.hh"
#include "casm/app/ProjectSettings.hh"

#include "casm/container/Enumerator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/casm_functions.hh"
#include "casm/completer/Handlers.hh"

namespace CASM {

  /// \brief Load enumerator plugins from a CASM project
  ///
  /// \param set CASM project settings to read enumerator plugins from
  /// \param enum_it Inserter to EnumeratorMap where the plugins should be stored
  /// \param lib_it Inserter to container of pair<std::string, std::shared_ptr<RuntimeLibrary> >
  ///        holding the libraries with the enumerator plugins
  ///
  /// - Checks the directory `primclex.dir().enumerator_plugins()` for `ENUMNAME.cc` files
  /// - All `*.cc` files found are compiled and linked, and the base of the
  ///   filename (ENUMNAME) is assumed to be the enumerator name
  /// - An extern C function named `make_ENUMNAME_interface` is assumed to exist
  ///   with signature: `CASM::EnumInterfaceBase* make_ENUMNAME_interface()`,
  ///   that returns a pointer which owns a `CASM::EnumInterfaceBase`-derived
  ///   object
  ///
  /// \note The lifetime of the RuntimeLibrary must be longer than the lifetime
  /// of the EnumInterface
  ///
  template<typename EnumeratorMapInserter, typename RuntimeLibInserter>
  std::pair<EnumeratorMapInserter, RuntimeLibInserter>
  load_enumerator_plugins(
    const ProjectSettings &set,
    EnumeratorMapInserter enum_it,
    RuntimeLibInserter lib_it) {

    const DirectoryStructure &dir = set.dir();

    if(dir.root_dir().empty()) {
      return std::make_pair(enum_it, lib_it);
    }

    if(fs::is_directory(dir.enumerator_plugins())) {

      // loop over custom enumerator files *.cc
      for(auto &entry : boost::make_iterator_range(fs::directory_iterator(dir.enumerator_plugins()), {})) {

        fs::path p = entry.path();
        std::string p_s = p.string();
        auto p_size = p_s.size();

        if(fs::is_regular_file(p) && p_s.compare(p_size - 3, p_size, ".cc") == 0) {

          fs::path f = p.filename();
          std::string f_s = f.string();
          auto f_size = f_s.size();

          std::string msg = "compiling new custom enumerator: " + f_s.substr(0, f_size - 3);

          // '-L$CASM_PREFIX/.libs' is a hack so 'make check' works
          auto lib_ptr = std::make_shared<RuntimeLibrary>(
                           p_s.substr(0, p_size - 3),
                           set.compile_options() + " " + include_path(dir.enumerator_plugins()),
                           set.so_options() + " -lcasm " + link_path(set.casm_libdir().first),
                           msg,
                           set);

          auto make_interface = lib_ptr->get_function<EnumInterfaceBase* ()>(
                                  "make_" + f_s.substr(0, f_size - 3) + "_interface");

          std::unique_ptr<EnumInterfaceBase> ptr(make_interface());

          // will clone on insert
          *enum_it++ = *ptr;
          *lib_it++ = std::make_pair(ptr->name(), lib_ptr);
        }
      }
    }

    return std::make_pair(enum_it, lib_it);
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
  template<typename ScelIterator, typename ConfigEnumConstructor>
  int insert_unique_canon_configs(
    std::string method,
    PrimClex &primclex,
    ScelIterator begin,
    ScelIterator end,
    ConfigEnumConstructor f,
    std::vector<std::string> filter_expr) {

    Log &log = primclex.log();

    Index Ninit = std::distance(primclex.config_begin(), primclex.config_end());
    log << "# configurations in this project: " << Ninit << "\n" << std::endl;

    log.begin(method);

    for(auto scel_it = begin; scel_it != end; ++scel_it) {
      Supercell &scel = *scel_it;
      log << "Enumerate configurations for " << scel.get_name() << " ...  " << std::flush;

      auto enumerator_ptr = f(scel);
      auto &enumerator = *enumerator_ptr;
      Index num_before = scel.get_config_list().size();
      if(!filter_expr.empty()) {
        try {
          scel.add_unique_canon_configs(
            filter_begin(
              enumerator.begin(),
              enumerator.end(),
              filter_expr,
              primclex.settings().config_io()),
            filter_end(enumerator.end())
          );
        }
        catch(std::exception &e) {
          primclex.err_log() << "Cannot filter configurations using the expression provided: \n" << e.what() << "\nExiting...\n";
          return ERR_INVALID_ARG;
        }
      }
      else {
        scel.add_unique_canon_configs(enumerator.begin(), enumerator.end());
      }

      log << (scel.get_config_list().size() - num_before) << " configs." << std::endl;
    }
    log << "  DONE." << std::endl << std::endl;

    Index Nfinal = std::distance(primclex.config_begin(), primclex.config_end());

    log << "# new configurations: " << Nfinal - Ninit << "\n";
    log << "# configurations in this project: " << Nfinal << "\n" << std::endl;

    log << "Write SCEL..." << std::endl;
    primclex.print_supercells();
    log << "  DONE" << std::endl << std::endl;

    log << "Writing config_list..." << std::endl;
    primclex.write_config_list();
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

    Index Ninit = std::distance(primclex.config_begin(), primclex.config_end());
    log << "# configurations in this project: " << Ninit << "\n" << std::endl;

    log.begin(method);

    for(auto scel_lat_it = begin; scel_lat_it != end; ++scel_lat_it) {
      Supercell scel(&primclex, *scel_lat_it);
      Supercell &canon_scel = scel.canonical_form();
      log << "Enumerate configurations for " << canon_scel.get_name() << " ...  " << std::flush;

      auto enumerator_ptr = f(scel);
      auto &enumerator = *enumerator_ptr;
      Index num_before = canon_scel.get_config_list().size();
      if(!filter_expr.empty()) {
        try {
          auto it = filter_begin(
                      enumerator.begin(),
                      enumerator.end(),
                      filter_expr,
                      primclex.settings().config_io());
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

      log << (canon_scel.get_config_list().size() - num_before) << " configs." << std::endl;
    }
    log << "  DONE." << std::endl << std::endl;

    Index Nfinal = std::distance(primclex.config_begin(), primclex.config_end());

    log << "# new configurations: " << Nfinal - Ninit << "\n";
    log << "# configurations in this project: " << Nfinal << "\n" << std::endl;

    log << "Write SCEL..." << std::endl;
    primclex.print_supercells();
    log << "  DONE" << std::endl << std::endl;

    log << "Writing config_list..." << std::endl;
    primclex.write_config_list();
    log << "  DONE" << std::endl;

    return 0;
  }

}
