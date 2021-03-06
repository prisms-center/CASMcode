#include <boost/filesystem.hpp>

#include "casm/app/EnumeratorHandler.hh"
#include "casm/app/LogRuntimeLibrary.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/enum/EnumInterface.hh"
#include "casm/completer/Handlers.hh"
#include "casm/enumerator/Enumerator.hh"
#include "casm/system/RuntimeLibrary.hh"

namespace CASM {

/// \brief Load enumerator plugins from a CASM project
///
/// \param set CASM project settings to read enumerator plugins from
/// \param enum_it Inserter to EnumInterfaceVector where the plugins should be
/// stored \param lib_it Inserter to container of pair<std::string,
/// std::shared_ptr<RuntimeLibrary> >
///        holding the libraries with the enumerator plugins
///
/// - Checks the directory `primclex.dir().enumerator_plugins()` for
/// `ENUMNAME.cc` files
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
template <typename EnumInterfaceVectorInserter, typename RuntimeLibInserter>
std::pair<EnumInterfaceVectorInserter, RuntimeLibInserter>
load_enumerator_plugins(ProjectSettings const &set,
                        EnumInterfaceVectorInserter enum_it,
                        RuntimeLibInserter lib_it) {
  const DirectoryStructure &dir = set.dir();

  if (dir.root_dir().empty()) {
    return std::make_pair(enum_it, lib_it);
  }

  if (fs::is_directory(dir.enumerator_plugins())) {
    // loop over custom enumerator files *.cc
    for (auto &entry : boost::make_iterator_range(
             fs::directory_iterator(dir.enumerator_plugins()), {})) {
      fs::path p = entry.path();
      std::string p_s = p.string();
      auto p_size = p_s.size();

      if (fs::is_regular_file(p) &&
          p_s.compare(p_size - 3, p_size, ".cc") == 0) {
        fs::path f = p.filename();
        std::string f_s = f.string();
        auto f_size = f_s.size();

        std::string msg =
            "compiling new custom enumerator: " + f_s.substr(0, f_size - 3);

        auto lib_ptr = log_make_shared_runtime_lib(
            p_s.substr(0, p_size - 3),
            set.compile_options() + " " +
                include_path(dir.enumerator_plugins()),
            set.so_options() + " -lcasm ", msg);

        auto make_interface = lib_ptr->get_function<EnumInterfaceBase *()>(
            "make_" + f_s.substr(0, f_size - 3) + "_interface");

        notstd::cloneable_ptr<EnumInterfaceBase> ptr(make_interface());

        std::string method_name = ptr->name();
        *enum_it++ = std::move(ptr);
        *lib_it++ = std::make_pair(method_name, lib_ptr);
      }
    }
  }

  return std::make_pair(enum_it, lib_it);
}

}  // namespace CASM
