#ifndef CASM_Enumerator_impl
#define CASM_Enumerator_impl

#include "casm/container/Enumerator.hh"
#include "casm/clex/PrimClex.hh"

namespace CASM {

  /// \brief Load enumerator plugins from a CASM project
  ///
  /// \param primclex CASM project to read enumerator plugins from
  /// \param enum_it Inserter to EnumeratorMap where the plugins should be stored
  /// \param lib_it Inserter to container of std::shared_ptr<RuntimeLibrary>
  ///        holding the libraries with the enumerator plugins
  ///
  /// - Checks the directory `primclex.dir().enumerators()` for `ENUMNAME.cc` files
  /// - All `*.cc` files found are compiled and linked, and the base of the
  ///   filename (ENUMNAME) is assumed to be the enumerator name
  /// - An extern C function named `make_ENUMNAME_interface` is assumed to exist
  ///   with signature: `CASM::EnumInterfaceBase* make_ENUMNAME_interface()`,
  ///   that returns a pointer which owns a `CASM::EnumInterfaceBase`-derived
  ///   object
  ///
  /// \note The lifetime of the RuntimeLibrary must be longer than the lifetime
  /// of the EnumeratorMap
  ///
  template<typename EnumeratorMapInserter, typename RuntimeLibInserter>
  std::pair<EnumeratorMapInserter, RuntimeLibInserter>
  load_enumerator_plugins(
    const PrimClex &primclex,
    EnumeratorMapInserter enum_it,
    RuntimeLibInserter lib_it) {

    // If 'args.primclex', use that, else construct PrimClex in 'primclex'
    // Then whichever exists, store reference in 'primclex'
    const DirectoryStructure &dir = primclex.dir();
    const ProjectSettings &set = primclex.settings();

    if(fs::is_directory(dir.enumerators())) {

      // loop over custom enumerator files *.cc
      for(auto &entry : boost::make_iterator_range(fs::directory_iterator(dir.enumerators()), {})) {

        fs::path p = entry.path();
        std::string p_s = p.string();
        auto p_size = p_s.size();

        if(fs::is_regular_file(p) && p_s.compare(p_size - 3, p_size, ".cc") == 0) {

          fs::path f = p.filename();
          std::string f_s = f.string();
          auto f_size = f_s.size();

          std::string msg = "compiling new custom enumerator: " + f_s.substr(0, f_size - 3);

          auto lib_ptr = std::make_shared<RuntimeLibrary>(
                           p_s.substr(0, p_size - 3),
                           set.compile_options() + " -I" + dir.enumerators().string(),
                           set.so_options() + " -lcasm -L" + (set.casm_prefix().first / "lib").string(),
                           primclex.log(),
                           msg);

          auto make_interface = lib_ptr->get_function<EnumInterfaceBase* ()>(
                                  "make_" + f_s.substr(0, f_size - 3) + "_interface");

          std::unique_ptr<EnumInterfaceBase> ptr(make_interface());

          // will clone on insert
          *enum_it++ = *ptr;
          *lib_it++ = lib_ptr;
        }
      }
    }

    return std::make_pair(enum_it, lib_it);
  }

}

#endif
