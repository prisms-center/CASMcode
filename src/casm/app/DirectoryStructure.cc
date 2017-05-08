#include "casm/app/DirectoryStructure.hh"

#include <boost/filesystem.hpp>
#include "casm/clex/Configuration.hh"
#include "casm/casm_io/Log.hh"

namespace CASM {

  /// return path to current or parent directory containing ".casm" directory
  ///   if none found, return empty path
  fs::path find_casmroot(const fs::path &cwd) {
    fs::path dir(cwd);
    fs::path casmroot(".casm");

    while(!dir.empty()) {
      // if cwd contains ".casm", return cwd
      if(fs::is_directory(dir / casmroot))
        return dir;
      dir = dir.parent_path();
    }
    return dir;
  };

  /// return relative path to current or parent directory containing ".casm" directory
  ///   if none found, return empty path
  fs::path relative_casmroot(const fs::path &cwd) {
    fs::path dir(cwd);
    fs::path casmroot = find_casmroot(cwd);
    fs::path relpath("");

    while(dir != casmroot) {
      dir = dir.parent_path();
      relpath /= "..";
    }
    return relpath;
  };

  /// \brief Remove files recursively
  void recurs_rm_files(fs::path p, bool dry_run, Log &log) {
    if(!fs::exists(p)) {
      return;
    }

    auto it = fs::directory_iterator(p);
    auto end = fs::directory_iterator();

    for(; it != end; ++it) {
      if(fs::is_regular_file(*it)) {
        log << "rm " << *it << std::endl;
        if(!dry_run) {
          fs::remove(*it);
        }
      }
      else {
        recurs_rm_files(*it, dry_run, log);
        log << "rm " << *it << std::endl;
      }
    }

    log << "rm " << p << std::endl;
    if(!dry_run) {
      fs::remove(p);
    }
  }

  namespace {
    /// \brief Copy files recursively, and returns a count of copied files
    Index recurs_cp_files_impl(const fs::path &from_dir, const fs::path &to_dir, bool dry_run, Index count, Log &log) {
      auto it = fs::directory_iterator(from_dir);
      auto end = fs::directory_iterator();
      for(; it != end; ++it) {
        if(fs::is_regular_file(*it)) {
          log << "cp " << *it << " " << to_dir << std::endl;
          count++;
          if(!dry_run) {
            fs::copy_file(*it, to_dir / it->path().filename());
          }
        }
        else {
          fs::path new_to_dir = to_dir / it->path().filename();
          if(!dry_run) {
            fs::create_directories(new_to_dir);
          }
          count = recurs_cp_files_impl(*it, new_to_dir, dry_run, count, log);
        }
      }
      return count;
    }
  }

  /// \brief Copy files recursively, and returns a count of copied files
  Index recurs_cp_files(const fs::path &from_dir, const fs::path &to_dir, bool dry_run, Log &log) {
    Index count = 0;
    recurs_cp_files_impl(from_dir, to_dir, dry_run, count, log);
    return count;
  }

  DirectoryStructure::DirectoryStructure(const fs::path _root) {
    _init(fs::absolute(_root));
  }


  // ** Query filesystem **

  /// \brief Check filesystem directory structure and return list of all basis set names
  std::vector<std::string> DirectoryStructure::all_bset() const {
    return _all_settings("bset", m_root / m_bset_dir);
  }

  /// \brief Check filesystem directory structure and return list of all calctype names
  std::vector<std::string> DirectoryStructure::all_calctype() const {
    return _all_settings("calctype", m_root / m_calc_dir / m_set_dir);
  }

  /// \brief Check filesystem directory structure and return list of all ref names for a given calctype
  std::vector<std::string> DirectoryStructure::all_ref(std::string calctype) const {
    return _all_settings("ref", calc_settings_dir(calctype));
  }

  /// \brief Check filesystem directory structure and return list of all property names
  std::vector<std::string> DirectoryStructure::all_property() const {
    return _all_settings("clex", m_root / m_clex_dir);
  }

  /// \brief Check filesystem directory structure and return list of all eci names
  std::vector<std::string> DirectoryStructure::all_eci(std::string property, std::string calctype, std::string ref, std::string bset) const {
    return _all_settings("eci", m_root / m_clex_dir / _property(property) / _calctype(calctype) / _ref(ref) / _bset(bset));
  }


  // ** File and Directory paths **


  // -- Project directory --------

  /// \brief Return casm project directory path
  fs::path DirectoryStructure::root_dir() const {
    return m_root;
  }

  /// \brief Return prim.json path
  fs::path DirectoryStructure::prim() const {
    return m_root / "prim.json";
  }

  /// \brief Return PRIM path
  fs::path DirectoryStructure::PRIM() const {
    return m_root / "PRIM";
  }


  // -- Hidden .casm directory --------

  /// \brief Return hidden .casm dir path
  fs::path DirectoryStructure::casm_dir() const {
    return m_root / m_casm_dir;
  }

  /// \brief Return project_settings.json path
  fs::path DirectoryStructure::project_settings() const {
    return m_root / m_casm_dir / "project_settings.json";
  }

  /// \brief Return master scel_list.json path
  fs::path DirectoryStructure::scel_list() const {
    return m_root / m_casm_dir / "scel_list.json";
  }

  /// \brief Return master config_list.json file path
  fs::path DirectoryStructure::config_list() const {
    return m_root / m_casm_dir / "config_list.json";
  }

  /// \brief Return enumerators plugin dir
  fs::path DirectoryStructure::enumerator_plugins() const {
    return m_root / m_casm_dir / "enumerators";
  }

  /// \brief Return enumerators plugin dir
  template<typename DataObject>
  fs::path DirectoryStructure::query_plugins() const {
    return m_root / m_casm_dir / "query" / traits<DataObject>::name;
  }

  template<typename DataObject>
  fs::path DirectoryStructure::master_selection() const {
    return query_plugins<DataObject>() / "master_selection";
  }

  /// \brief File containing DataObject name aliases (not query function aliases)
  template<typename DataObject>
  fs::path DirectoryStructure::aliases() const {
    return query_plugins<DataObject>() / "aliases.json";
  }

  // -- Symmetry --------

  /// \brief Return symmetry directory path
  fs::path DirectoryStructure::symmetry_dir() const {
    return m_root / m_sym_dir;
  }

  /// \brief Return lattice_point_group.json path
  fs::path DirectoryStructure::lattice_point_group() const {
    return m_root / m_sym_dir / "lattice_point_group.json";
  }

  /// \brief Return factor_group.json path
  fs::path DirectoryStructure::factor_group() const {
    return m_root / m_sym_dir / "factor_group.json";
  }

  /// \brief Return crystal_point_group.json path
  fs::path DirectoryStructure::crystal_point_group() const {
    return m_root / m_sym_dir / "crystal_point_group.json";
  }


  // -- Basis sets --------

  /// \brief Return path to directory contain basis set info
  fs::path DirectoryStructure::bset_dir(std::string bset) const {
    return m_root / m_bset_dir / _bset(bset);
  }

  /// \brief Return basis function specs (bspecs.json) file path
  fs::path DirectoryStructure::bspecs(std::string bset) const {
    return bset_dir(bset) / "bspecs.json";
  }

  // \brief Returns path to the clust.json file
  fs::path DirectoryStructure::clust(std::string bset) const {
    return bset_dir(bset) / "clust.json";
  }

  // \brief Returns path to the basis.json file
  fs::path DirectoryStructure::basis(std::string bset) const {
    return bset_dir(bset) / "basis.json";
  }

  /// \brief Returns path to directory containing global clexulator
  fs::path DirectoryStructure::clexulator_dir(std::string bset) const {
    return bset_dir(bset);
  }

  /// \brief Returns path to global clexulator source file
  fs::path DirectoryStructure::clexulator_src(std::string project, std::string bset) const {
    return bset_dir(bset) / (project + "_Clexulator.cc");
  }

  /// \brief Returns path to global clexulator o file
  fs::path DirectoryStructure::clexulator_o(std::string project, std::string bset) const {
    return bset_dir(bset) / (project + "_Clexulator.o");
  }

  /// \brief Returns path to global clexulator so file
  fs::path DirectoryStructure::clexulator_so(std::string project, std::string bset) const {
    return bset_dir(bset) / (project + "_Clexulator.so");
  }

  /// \brief Returns path to eci.in, in bset directory
  fs::path DirectoryStructure::eci_in(std::string bset) const {
    return bset_dir(bset) / "eci.in";
  }

  /// \brief Returns path to corr.in, in bset directory
  fs::path DirectoryStructure::corr_in(std::string bset) const {
    return bset_dir(bset) / "corr.in";
  }


  // -- Calculations and reference --------

  /// \brief Return 'training_data' directorty path
  fs::path DirectoryStructure::training_data() const {
    return m_root / m_calc_dir;
  }

  /// \brief Return SCEL path
  fs::path DirectoryStructure::SCEL() const {
    return m_root / m_calc_dir / "SCEL";
  }

  /// \brief Return supercell directory path (scelname has format SCELV_A_B_C_D_E_F)
  fs::path DirectoryStructure::supercell_dir(std::string scelname) const {
    return m_root / m_calc_dir / scelname;
  }

  /// \brief Return supercell LAT file path (scelname has format SCELV_A_B_C_D_E_F)
  fs::path DirectoryStructure::LAT(std::string scelname) const {
    return m_root / m_calc_dir / scelname / "LAT";
  }

  /// \brief Return configuration directory path (configname has format SCELV_A_B_C_D_E_F/I)
  fs::path DirectoryStructure::configuration_dir(std::string configname) const {
    return m_root / m_calc_dir / configname;
  }

  /// \brief Return path to POS file
  fs::path DirectoryStructure::POS(std::string configname) const {
    return configuration_dir(configname) / "POS";
  }

  /// \brief Return calculation settings directory path, for global settings
  fs::path DirectoryStructure::calc_settings_dir(std::string calctype) const {
    return m_root / m_calc_dir / m_set_dir / _calctype(calctype);
  }

  /// \brief Return calculation settings directory path, for supercell specific settings
  fs::path DirectoryStructure::supercell_calc_settings_dir(std::string scelname, std::string calctype) const {
    return supercell_dir(scelname) / m_set_dir / _calctype(calctype);
  }

  /// \brief Return calculation settings directory path, for configuration specific settings
  fs::path DirectoryStructure::configuration_calc_settings_dir(std::string configname, std::string calctype) const {
    return configuration_dir(configname) / m_set_dir / _calctype(calctype);
  }

  /// \brief Return directory containing properties.calc.json
  fs::path DirectoryStructure::configuration_calc_dir(std::string configname, std::string calctype) const {
    return configuration_dir(configname) / _calctype(calctype);
  }

  /// \brief Return properties.calc.json file path
  fs::path DirectoryStructure::calculated_properties(std::string configname, std::string calctype) const {
    return configuration_calc_dir(configname, calctype) / "properties.calc.json";
  }

  /// \brief Return calculation status file path
  fs::path DirectoryStructure::calc_status(std::string configname, std::string calctype) const {
    return configuration_calc_dir(configname, calctype) / "status.json";
  }


  /// \brief Return calculation reference settings directory path, for global settings
  fs::path DirectoryStructure::ref_dir(std::string calctype, std::string ref) const {
    return calc_settings_dir(calctype) / _ref(ref);
  }

  /// \brief Return composition axes file path
  fs::path DirectoryStructure::composition_axes() const {
    return casm_dir() / "composition_axes.json";
  }

  /// \brief Return chemical reference file path
  fs::path DirectoryStructure::chemical_reference(std::string calctype, std::string ref) const {
    return ref_dir(calctype, ref) / "chemical_reference.json";
  }


  // -- Cluster expansions --------

  /// \brief Returns path to eci directory
  fs::path DirectoryStructure::clex_dir(std::string property) const {
    return m_root / m_clex_dir / _property(property);
  }

  /// \brief Returns path to eci directory
  fs::path DirectoryStructure::eci_dir(std::string property, std::string calctype, std::string ref, std::string bset, std::string eci) const {
    return clex_dir(property) / _calctype(calctype) / _ref(ref) / _bset(bset) / _eci(eci);
  }

  /// \brief Returns path to eci.json
  fs::path DirectoryStructure::eci(std::string property, std::string calctype, std::string ref, std::string bset, std::string eci) const {
    return eci_dir(property, calctype, ref, bset, eci) / "eci.json";
  }


  // -- other maybe temporary --------------------------

  /// \brief Return cluster specs (CSPECS) file path
  fs::path DirectoryStructure::CSPECS(std::string bset) const {
    return bset_dir(bset) / "CSPECS";
  }

  // \brief Returns path to the clust.json file
  fs::path DirectoryStructure::FCLUST(std::string bset) const {
    return bset_dir(bset) / "FCLUST.json";
  }

  // -- deprecated ------------------------------------

  /// \brief Returns path to eci.out
  fs::path DirectoryStructure::eci_out(std::string property, std::string calctype, std::string ref, std::string bset, std::string eci) const {
    return eci_dir(property, calctype, ref, bset, eci) / "eci.out";
  }

  /// \brief Query aliases file
  fs::path DirectoryStructure::query_alias() const {
    return m_root / m_casm_dir / "query_alias.json";
  }

  std::string DirectoryStructure::_bset(std::string bset) const {
    return std::string("bset.") + bset;
  }

  std::string DirectoryStructure::_calctype(std::string calctype) const {
    return std::string("calctype.") + calctype;
  }

  std::string DirectoryStructure::_ref(std::string ref) const {
    return std::string("ref.") + ref;
  }

  std::string DirectoryStructure::_property(std::string property) const {
    return std::string("clex.") + property;
  }

  std::string DirectoryStructure::_eci(std::string eci) const {
    return std::string("eci.") + eci;
  }

  std::string DirectoryStructure::_ref_state(int index) const {
    return std::string("properties.ref_state.") + std::to_string(index) + ".json";
  }


  void DirectoryStructure::_init(const fs::path &_root) {
    m_root = _root;
    m_casm_dir = ".casm";
    m_bset_dir = "basis_sets";
    m_calc_dir = "training_data";
    m_set_dir = "settings";
    m_sym_dir = "symmetry";
    m_clex_dir = "cluster_expansions";
  }

  /// \brief Find all directories at 'location' that match 'pattern.something'
  ///        and return a std::vector of the 'something'
  std::vector<std::string> DirectoryStructure::_all_settings(std::string pattern, fs::path location) const {

    std::vector<std::string> all;
    std::string dir;
    pattern += ".";

    // get all
    if(!fs::exists(location)) {
      return all;
    }
    fs::directory_iterator it(location);
    fs::directory_iterator end_it;
    for(; it != end_it; ++it) {
      if(fs::is_directory(*it)) {
        dir = it->path().filename().string();
        if(dir.substr(0, pattern.size()) == pattern) {
          all.push_back(dir.substr(pattern.size(), dir.size()));
        }
      }
    }

    std::sort(all.begin(), all.end());

    return all;
  }

}

#include "casm/database/DatabaseTypeTraits.hh"

// explicit instantiations
#define INST_DirectoryStructure_all(r, data, type) \
template fs::path DirectoryStructure::query_plugins<type>() const; \
template fs::path DirectoryStructure::master_selection<type>() const; \
template fs::path DirectoryStructure::aliases<type>() const;

namespace CASM {
  BOOST_PP_SEQ_FOR_EACH(INST_DirectoryStructure_all, _, CASM_DB_TYPES)
}

