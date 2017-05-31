#ifndef CASM_DirectoryStructure
#define CASM_DirectoryStructure

#include <string>
#include <vector>

#include <boost/filesystem/path.hpp>
#include "casm/CASM_global_definitions.hh"


namespace CASM {

  /**
   * \ingroup Project
   *
   * @{
   */

  class Log;

  /// return path to current or parent directory containing ".casm" directory
  ///   if none found, return empty path
  fs::path find_casmroot(const fs::path &cwd);

  /// return relative path to current or parent directory containing ".casm" directory
  ///   if none found, return empty path
  fs::path relative_casmroot(const fs::path &cwd);

  /// \brief Remove files recursively
  void recurs_rm_files(fs::path p, bool dry_run, Log &log);

  /// \brief Copy files recursively, and returns a count of copied files
  Index recurs_cp_files(const fs::path &from_dir, const fs::path &to_dir, bool dry_run, Log &log);

  /// \brief Specification of CASM project directory structure
  class DirectoryStructure {

  public:

    DirectoryStructure() {}

    DirectoryStructure(const fs::path _root);


    // ** Query filesystem **

    /// \brief Check filesystem directory structure and return list of all basis set names
    std::vector<std::string> all_bset() const;

    /// \brief Check filesystem directory structure and return list of all calctype names
    std::vector<std::string> all_calctype() const;

    /// \brief Check filesystem directory structure and return list of all ref names for a given calctype
    std::vector<std::string> all_ref(std::string calctype) const;

    /// \brief Check filesystem directory structure and return list of all property names
    std::vector<std::string> all_property() const;

    /// \brief Check filesystem directory structure and return list of all eci names
    std::vector<std::string> all_eci(std::string property, std::string calctype, std::string ref, std::string bset) const;


    // ** File and Directory paths **


    // -- Project directory --------

    /// \brief Return casm project directory path
    fs::path root_dir() const;

    /// \brief Return prim.json path
    fs::path prim() const;

    /// \brief Return PRIM path
    fs::path PRIM() const;


    // -- Hidden .casm directory --------

    /// \brief Return hidden .casm dir path
    fs::path casm_dir() const;

    /// \brief Return project_settings.json path
    fs::path project_settings() const;

    /// \brief Return master scel_list.json path
    fs::path scel_list() const;

    /// \brief Return master config_list.json file path
    fs::path config_list() const;

    /// \brief Return enumerators plugin dir
    fs::path enumerator_plugins() const;

    /// \brief Return enumerators plugin dir
    template<typename DataObject>
    fs::path query_plugins() const;

    template<typename DataObject>
    fs::path master_selection() const;

    /// \brief File containing DataObject name aliases (not query function aliases)
    template<typename DataObject>
    fs::path aliases() const;

    // -- Symmetry --------

    /// \brief Return symmetry directory path
    fs::path symmetry_dir() const;

    /// \brief Return lattice_point_group.json path
    fs::path lattice_point_group() const;

    /// \brief Return factor_group.json path
    fs::path factor_group() const;

    /// \brief Return crystal_point_group.json path
    fs::path crystal_point_group() const;


    // -- Basis sets --------

    /// \brief Return path to directory contain basis set info
    fs::path bset_dir(std::string bset) const;

    /// \brief Return basis function specs (bspecs.json) file path
    fs::path bspecs(std::string bset) const;

    // \brief Returns path to the clust.json file
    fs::path clust(std::string bset) const;

    // \brief Returns path to the basis.json file
    fs::path basis(std::string bset) const;

    /// \brief Returns path to directory containing global clexulator
    fs::path clexulator_dir(std::string bset) const;

    /// \brief Returns path to global clexulator source file
    fs::path clexulator_src(std::string project, std::string bset) const;

    /// \brief Returns path to global clexulator o file
    fs::path clexulator_o(std::string project, std::string bset) const;

    /// \brief Returns path to global clexulator so file
    fs::path clexulator_so(std::string project, std::string bset) const;

    /// \brief Returns path to eci.in, in bset directory
    fs::path eci_in(std::string bset) const;

    /// \brief Returns path to corr.in, in bset directory
    fs::path corr_in(std::string bset) const;


    // -- Calculations and reference --------

    /// \brief Return 'training_data' directorty path
    fs::path training_data() const;

    /// \brief Return SCEL path
    fs::path SCEL() const;

    /// \brief Return supercell directory path (scelname has format SCELV_A_B_C_D_E_F)
    fs::path supercell_dir(std::string scelname) const;

    /// \brief Return supercell LAT file path (scelname has format SCELV_A_B_C_D_E_F)
    fs::path LAT(std::string scelname) const;

    /// \brief Return configuration directory path (configname has format SCELV_A_B_C_D_E_F/I)
    fs::path configuration_dir(std::string configname) const;

    /// \brief Return path to POS file
    fs::path POS(std::string configname) const;

    /// \brief Return calculation settings directory path, for global settings
    fs::path calc_settings_dir(std::string calctype) const;

    /// \brief Return calculation settings directory path, for supercell specific settings
    fs::path supercell_calc_settings_dir(std::string scelname, std::string calctype) const;

    /// \brief Return calculation settings directory path, for configuration specific settings
    fs::path configuration_calc_settings_dir(std::string configname, std::string calctype) const;

    /// \brief Return directory containing properties.calc.json
    fs::path configuration_calc_dir(std::string configname, std::string calctype) const;

    /// \brief Return properties.calc.json file path
    fs::path calculated_properties(std::string configname, std::string calctype) const;

    /// \brief Return calculation status file path
    fs::path calc_status(std::string configname, std::string calctype) const;


    /// \brief Return calculation reference settings directory path, for global settings
    fs::path ref_dir(std::string calctype, std::string ref) const;

    /// \brief Return composition axes file path
    fs::path composition_axes() const;

    /// \brief Return chemical reference file path
    fs::path chemical_reference(std::string calctype, std::string ref) const;


    // -- Cluster expansions --------

    /// \brief Returns path to eci directory
    fs::path clex_dir(std::string property) const;

    /// \brief Returns path to eci directory
    fs::path eci_dir(std::string property, std::string calctype, std::string ref, std::string bset, std::string eci) const;

    /// \brief Returns path to eci.json
    fs::path eci(std::string property, std::string calctype, std::string ref, std::string bset, std::string eci) const;


    // -- other maybe temporary --------------------------

    /// \brief Return cluster specs (CSPECS) file path
    fs::path CSPECS(std::string bset) const;

    // \brief Returns path to the clust.json file
    fs::path FCLUST(std::string bset) const;

    // -- deprecated ------------------------------------

    /// \brief Returns path to eci.out
    fs::path eci_out(std::string property, std::string calctype, std::string ref, std::string bset, std::string eci) const;

    /// \brief Query aliases file (deprecated: now stored in project_settings.json)
    fs::path query_alias() const;


  private:

    std::string _bset(std::string bset) const;

    std::string _calctype(std::string calctype) const;

    std::string _ref(std::string ref) const;

    std::string _property(std::string property) const;

    std::string _eci(std::string eci) const;

    std::string _ref_state(int index) const;


    void _init(const fs::path &_root);

    /// \brief Find all directories at 'location' that match 'pattern.something'
    ///        and return a std::vector of the 'something'
    std::vector<std::string> _all_settings(std::string pattern, fs::path location) const;

    fs::path m_root;
    std::string m_casm_dir;
    std::string m_bset_dir;
    std::string m_calc_dir;
    std::string m_set_dir;
    std::string m_sym_dir;
    std::string m_clex_dir;

  };

  /** @} */
}

#endif
