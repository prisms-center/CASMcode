#ifndef CASM_HamiltonianModules_impl
#define CASM_HamiltonianModules_impl

#include <boost/filesystem.hpp>

#include "casm/app/HamiltonianModules.hh"
#include "casm/app/LogRuntimeLibrary.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/system/RuntimeLibrary.hh"

namespace CASM {
/// \brief Load DoF plugins from a CASM project
template <typename DoFDictInserter, typename RuntimeLibInserter>
std::pair<DoFDictInserter, RuntimeLibInserter> load_dof_plugins(
    const ProjectSettings &set, DoFDictInserter dict_it,
    RuntimeLibInserter lib_it) {
  typedef DoFType::Traits *traits_ptr;
  typedef traits_ptr(signature)();

  const DirectoryStructure &dir = set.dir();

  if (dir.root_dir().empty()) {
    return std::make_pair(dict_it, lib_it);
  }

  if (fs::is_directory(dir.dof_plugins())) {
    // loop over custom query files *.cc
    for (auto &entry : boost::make_iterator_range(
             fs::directory_iterator(dir.dof_plugins()), {})) {
      fs::path p = entry.path();
      std::string p_s = p.string();
      auto p_size = p_s.size();

      if (fs::is_regular_file(p) &&
          p_s.compare(p_size - 3, p_size, ".cc") == 0) {
        fs::path f = p.filename();
        std::string f_s = f.string();
        auto f_size = f_s.size();

        std::string msg = "Compiling new custom DoFTraits for DoF: " +
                          f_s.substr(0, f_size - 3);

        auto lib_ptr = log_make_shared_runtime_lib(
            p_s.substr(0, p_size - 3),
            set.compile_options() + " " + include_path(dir.dof_plugins()),
            set.so_options() + " -lcasm ", msg);

        auto make_dof = lib_ptr->template get_function<signature>(
            "make_" + f_s.substr(0, f_size - 3) + "_dof");

        std::unique_ptr<DoFType::Traits> ptr(make_dof());

        // will clone on insert
        *dict_it++ = *ptr;
        *lib_it++ = std::make_pair(ptr->name(), lib_ptr);
      }
    }
  }

  return std::make_pair(dict_it, lib_it);
}

/// \brief Load SymRepBuilder plugins from a CASM project
template <typename SymRepBuilderDictInserter, typename RuntimeLibInserter>
std::pair<SymRepBuilderDictInserter, RuntimeLibInserter>
load_symrep_builder_plugins(const ProjectSettings &set,
                            SymRepBuilderDictInserter dict_it,
                            RuntimeLibInserter lib_it) {
  typedef SymRepBuilderInterface *bldr_ptr;
  typedef bldr_ptr(signature)();

  const DirectoryStructure &dir = set.dir();

  if (dir.root_dir().empty()) {
    return std::make_pair(dict_it, lib_it);
  }

  if (fs::is_directory(dir.symrep_builder_plugins())) {
    // loop over custom query files *.cc
    for (auto &entry : boost::make_iterator_range(
             fs::directory_iterator(dir.symrep_builder_plugins()), {})) {
      fs::path p = entry.path();
      std::string p_s = p.string();
      auto p_size = p_s.size();

      if (fs::is_regular_file(p) &&
          p_s.compare(p_size - 3, p_size, ".cc") == 0) {
        fs::path f = p.filename();
        std::string f_s = f.string();
        auto f_size = f_s.size();

        std::string msg = "Compiling new custom SymRepBuilderInterface: " +
                          f_s.substr(0, f_size - 3);

        // '-L$CASM_PREFIX/.libs' is a hack so 'make check' works
        auto lib_ptr = log_make_shared_runtime_lib(
            p_s.substr(0, p_size - 3),
            set.compile_options() + " " +
                include_path(dir.symrep_builder_plugins()),
            set.so_options() + " -lcasm ", msg);

        auto make_traits = lib_ptr->template get_function<signature>(
            "make_" + f_s.substr(0, f_size - 3) + "_symrep_builder");

        std::unique_ptr<SymRepBuilderInterface> ptr(make_traits());
        // will clone on insert
        *dict_it++ = *ptr;
        *lib_it++ = std::make_pair(ptr->name(), lib_ptr);
      }
    }
  }

  return std::make_pair(dict_it, lib_it);
}
}  // namespace CASM

#endif
