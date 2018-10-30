#ifndef CASM_HamiltonianModules_impl
#define CASM_HamiltonianModules_impl
#include <boost/filesystem.hpp>
#include "casm/app/HamiltonianModules.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/basis_set/DoF.hh"
#include "casm/crystallography/MoleculeAttribute.hh"
namespace CASM {
  /// \brief Load DoF plugins from a CASM project
  template<typename DoFDictInserter, typename RuntimeLibInserter>
  std::pair<DoFDictInserter, RuntimeLibInserter>
  load_dof_plugins(
    const ProjectSettings &set,
    DoFDictInserter dict_it,
    RuntimeLibInserter lib_it) {
    typedef DoF::BasicTraits *traits_ptr;
    typedef traits_ptr(signature)();

    const DirectoryStructure &dir = set.dir();

    if(dir.root_dir().empty()) {
      return std::make_pair(dict_it, lib_it);
    }

    if(fs::is_directory(dir.dof_plugins())) {

      // loop over custom query files *.cc
      for(auto &entry : boost::make_iterator_range(fs::directory_iterator(dir.dof_plugins()), {})) {

        fs::path p = entry.path();
        std::string p_s = p.string();
        auto p_size = p_s.size();

        if(fs::is_regular_file(p) && p_s.compare(p_size - 3, p_size, ".cc") == 0) {

          fs::path f = p.filename();
          std::string f_s = f.string();
          auto f_size = f_s.size();

          std::string msg = "compiling new custom dof: " + f_s.substr(0, f_size - 3);

          // '-L$CASM_PREFIX/.libs' is a hack so 'make check' works
          auto lib_ptr = std::make_shared<RuntimeLibrary>(
                           p_s.substr(0, p_size - 3),
                           set.compile_options() + " " + include_path(dir.dof_plugins()),
                           set.so_options() + " -lcasm ",
                           msg,
                           set);

          auto make_dof = lib_ptr->template get_function<signature>(
            "make_" + f_s.substr(0, f_size - 3) + "_dof");

          std::unique_ptr<DoF::BasicTraits> ptr(make_dof());

          // will clone on insert
          *dict_it++ = *ptr;
          *lib_it++ = std::make_pair(ptr->name(), lib_ptr);
        }
      }
    }

    return std::make_pair(dict_it, lib_it);
  }


  /// \brief Load MoleculeAttribute plugins from a CASM project
  template<typename AttributeDictInserter, typename RuntimeLibInserter>
  std::pair<AttributeDictInserter, RuntimeLibInserter>
  load_mol_attribute_plugins(
    const ProjectSettings &set,
    AttributeDictInserter dict_it,
    RuntimeLibInserter lib_it) {
    typedef MoleculeAttribute_impl::BasicTraits *traits_ptr;
    typedef traits_ptr(signature)();

    const DirectoryStructure &dir = set.dir();

    if(dir.root_dir().empty()) {
      return std::make_pair(dict_it, lib_it);
    }

    if(fs::is_directory(dir.molecule_traits_plugins())) {

      // loop over custom query files *.cc
      for(auto &entry : boost::make_iterator_range(fs::directory_iterator(dir.molecule_traits_plugins()), {})) {

        fs::path p = entry.path();
        std::string p_s = p.string();
        auto p_size = p_s.size();

        if(fs::is_regular_file(p) && p_s.compare(p_size - 3, p_size, ".cc") == 0) {

          fs::path f = p.filename();
          std::string f_s = f.string();
          auto f_size = f_s.size();

          std::string msg = "compiling new custom molecule attribute: " + f_s.substr(0, f_size - 3);

          // '-L$CASM_PREFIX/.libs' is a hack so 'make check' works
          auto lib_ptr = std::make_shared<RuntimeLibrary>(
                           p_s.substr(0, p_size - 3),
                           set.compile_options() + " " + include_path(dir.molecule_traits_plugins()),
                           set.so_options() + " -lcasm ",
                           msg,
                           set);

          auto make_traits = lib_ptr->template get_function<signature>(
            "make_" + f_s.substr(0, f_size - 3) + "_molecule_traits");

          std::unique_ptr<MoleculeAttribute_impl::BasicTraits> ptr(make_traits());

          // will clone on insert
          *dict_it++ = *ptr;
          *lib_it++ = std::make_pair(ptr->name(), lib_ptr);
        }
      }
    }

    return std::make_pair(dict_it, lib_it);
  }
}

#endif
