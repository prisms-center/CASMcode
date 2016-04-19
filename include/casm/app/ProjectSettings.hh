#ifndef CASM_ProjectSettings
#define CASM_ProjectSettings

#include <string>
#include <vector>

#include "casm/casm_io/SafeOfstream.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/app/DirectoryStructure.hh"

#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  /// \brief Read/modify settings of an already existing CASM project
  ///
  /// - Use ProjectBuilder to create a new CASM project
  /// - Only allows modifying settings if the appropriate directories exist
  ///
  class ProjectSettings {

  public:

    /// \brief Default constructor
    ProjectSettings() {}

    /// \brief Construct CASM project settings for a new project
    ///
    /// \param root Path to new CASM project directory
    /// \param name Name of new CASM project. Use a short title suitable for prepending to file names.
    ///
    explicit ProjectSettings(fs::path root, std::string name);

    /// \brief Construct CASM project settings from existing project
    ///
    /// \param root Path to existing CASM project directory. Project settings will be read.
    ///
    explicit ProjectSettings(fs::path root);


    /// \brief Get project name
    std::string name() const;

    /// \brief Get current properties
    const std::vector<std::string> &properties() const;

    /// \brief Get current basis set name
    std::string bset() const;

    /// \brief Get current calctype name
    std::string calctype() const;

    /// \brief Get current ref name
    std::string ref() const;

    /// \brief Get current cluster expansion name
    std::string clex() const;

    /// \brief Get current eci name
    std::string eci() const;

    /// \brief Get neighbor list weight matrix
    Eigen::Matrix3l nlist_weight_matrix() const;

    /// \brief Get set of sublattice indices to include in neighbor lists
    const std::set<int> &nlist_sublat_indices() const;

    /// \brief Get current compilation options string
    std::string compile_options() const;

    /// \brief Get current shared library options string
    std::string so_options() const;

    /// \brief Get current command used by 'casm view'
    std::string view_command() const;

    /// \brief Get current project tol
    double tol() const;


    // ** Clexulator names **

    std::string global_clexulator() const;


    // ** Add directories for additional project data **

    /// \brief Create new project data directory
    bool new_casm_dir() const;

    /// \brief Create new symmetry directory
    bool new_symmetry_dir() const;

    /// \brief Add a basis set directory
    bool new_bset_dir(std::string bset) const;

    /// \brief Add a cluster expansion directory
    bool new_clex_dir(std::string clex) const;


    /// \brief Add calculation settings directory path
    bool new_calc_settings_dir(std::string calctype) const;

    /// \brief Add calculation settings directory path, for supercell specific settings
    bool new_supercell_calc_settings_dir(std::string scelname, std::string calctype) const;

    /// \brief Add calculation settings directory path, for configuration specific settings
    bool new_configuration_calc_settings_dir(std::string configname, std::string calctype) const;


    /// \brief Add a ref directory
    bool new_ref_dir(std::string calctype, std::string ref) const;


    /// \brief Add an eci directory
    bool new_eci_dir(std::string clex, std::string calctype, std::string ref, std::string bset, std::string eci) const;


    // ** Change current settings **

    /// \brief Access current properties
    std::vector<std::string> &properties();

    /// \brief Set current basis set to 'bset', if 'bset' exists
    bool set_bset(std::string bset);

    /// \brief Set current calctype to 'calctype', if 'calctype' exists
    bool set_calctype(std::string calctype);

    /// \brief Set current calculation reference to 'ref', if 'ref' exists
    bool set_ref(std::string calctype, std::string ref);

    /// \brief Set current cluster expansion to 'clex', if 'clex' exists
    bool set_clex(std::string clex);

    /// \brief Set current eci to 'eci', if 'eci' exists
    bool set_eci(std::string clex, std::string calctype, std::string ref, std::string bset, std::string eci);

    /// \brief Set neighbor list weight matrix (will delete existing Clexulator
    /// source and compiled code)
    bool set_nlist_weight_matrix(Eigen::Matrix3l M);

    /// \brief Set range of sublattice indices to include in neighbor lists (will
    /// delete existing Clexulator source and compiled code)
    template<typename SublatIterator>
    bool set_nlist_sublat_indices(SublatIterator begin, SublatIterator end);

    /// \brief Set compile options to 'opt'
    bool set_compile_options(std::string opt);

    /// \brief Set shared library options to 'opt'
    bool set_so_options(std::string opt);

    /// \brief Set command used by 'casm view'
    bool set_view_command(std::string opt);

    /// \brief Set shared library options to 'opt'
    bool set_tol(double _tol);


    /// \brief Save settings to project settings file
    void commit() const;


  private:

    /// \brief Changing the neighbor list properties requires updating Clexulator source code
    void _reset_clexulators();

    DirectoryStructure m_dir;

    std::string m_name;

    // CASM project current settings: used to determine where to write things
    std::string m_bset;
    std::string m_calctype;
    std::string m_ref;
    std::string m_clex;
    std::string m_eci;

    // neighbor list settings
    Eigen::Matrix3l m_nlist_weight_matrix;
    std::set<int> m_nlist_sublat_indices;

    // Properties to read from calculations
    std::vector<std::string> m_properties;

    // Runtime library compilation settings: compilation options
    std::string m_compile_options;
    std::string m_so_options;

    // Command executed by 'casm view'
    std::string m_view_command;

    // Default tolerance
    double m_tol;

  };

  jsonParser &to_json(const ProjectSettings &set, jsonParser &json);


  /// \brief Set range of sublattice indices to include in neighbor lists (will
  /// delete existing Clexulator source and compiled code)
  template<typename SublatIterator>
  bool ProjectSettings::set_nlist_sublat_indices(SublatIterator begin, SublatIterator end) {
    _reset_clexulators();
    m_nlist_sublat_indices = std::set<int>(begin, end);
    return true;
  }

}

#endif
