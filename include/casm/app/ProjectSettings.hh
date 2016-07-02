#ifndef CASM_ProjectSettings
#define CASM_ProjectSettings

#include <string>
#include <vector>

#include "casm/casm_io/SafeOfstream.hh"
#include "casm/system/RuntimeLibrary.hh"
#include "casm/app/DirectoryStructure.hh"

#include "casm/casm_io/DataFormatter.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/container.hh"

namespace CASM {

  class Configuration;

  template <bool IsConst>
  class ConfigSelection;
  typedef ConfigSelection<true> ConstConfigSelection;

  namespace ConfigIO {
    class Selected;
  }

  struct ClexDescription {
    ClexDescription() {}

    ClexDescription(std::string _name,
                    std::string _property,
                    std::string _calctype,
                    std::string _ref,
                    std::string _bset,
                    std::string _eci) :
      name(_name), property(_property), calctype(_calctype), ref(_ref), bset(_bset), eci(_eci) {}

    std::string name;
    std::string property;
    std::string calctype;
    std::string ref;
    std::string bset;
    std::string eci;
  };

  /// \brief Compare using name strings: A.name < B.name
  bool operator<(const ClexDescription &A, const ClexDescription &B);

  jsonParser &to_json(const ClexDescription &desc, jsonParser &json);

  void from_json(ClexDescription &desc, const jsonParser &json);

  bool clex_exists(const DirectoryStructure &dir, const ClexDescription &desc);


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


    const std::map<std::string, ClexDescription> &cluster_expansions() const;

    bool has_clex(std::string name) const;

    const ClexDescription &clex(std::string name) const;

    const ClexDescription &default_clex() const;

    bool new_clex(const ClexDescription &desc);

    bool erase_clex(const ClexDescription &desc);

    bool set_default_clex(const std::string &clex_name);

    bool set_default_clex(const ClexDescription &desc);


    /// \brief Get neighbor list weight matrix
    Eigen::Matrix3l nlist_weight_matrix() const;

    /// \brief Get set of sublattice indices to include in neighbor lists
    const std::set<int> &nlist_sublat_indices() const;

    /// \brief Get c++ compiler
    std::pair<std::string, std::string> cxx() const;

    /// \brief Get c++ compiler options
    std::pair<std::string, std::string> cxxflags() const;

    /// \brief Get shared object options
    std::pair<std::string, std::string> soflags() const;

    /// \brief Get casm prefix
    std::pair<fs::path, std::string> casm_prefix() const;

    /// \brief Get boost prefix
    std::pair<fs::path, std::string> boost_prefix() const;

    /// \brief Get current compilation options string
    std::string compile_options() const;

    /// \brief Get current shared library options string
    std::string so_options() const;

    /// \brief Get current command used by 'casm view'
    std::string view_command() const;

    /// \brief Get current project crystallography tolerance
    double crystallography_tol() const;

    /// \brief Get current project linear algebra tolerance
    double lin_alg_tol() const;


    // ** Configuration properties **

    const DataFormatterDictionary<Configuration> &config_io() const;

    /// \brief Set the selection to be used for the 'selected' column
    void set_selected(const ConfigIO::Selected &selection);

    /// \brief Set the selection to be used for the 'selected' column
    void set_selected(const ConstConfigSelection &selection);

    /// \brief Add user-defined query alias
    void add_alias(const std::string &alias_name, const std::string &alias_command, std::ostream &serr);

    /// \brief Return map containing aliases
    const std::map<std::string, std::string> &aliases() const;


    // ** Clexulator names **

    std::string clexulator() const;


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


    /// \brief Set neighbor list weight matrix (will delete existing Clexulator
    /// source and compiled code)
    bool set_nlist_weight_matrix(Eigen::Matrix3l M);

    /// \brief Set range of sublattice indices to include in neighbor lists (will
    /// delete existing Clexulator source and compiled code)
    template<typename SublatIterator>
    bool set_nlist_sublat_indices(SublatIterator begin, SublatIterator end);


    /// \brief Set c++ compiler (empty string to use default)
    bool set_cxx(std::string opt);

    /// \brief Set c++ compiler options (empty string to use default)
    bool set_cxxflags(std::string opt);

    /// \brief Set shared object options (empty string to use default)
    bool set_soflags(std::string opt);

    /// \brief Set casm prefix (empty string to use default)
    bool set_casm_prefix(fs::path prefix);

    /// \brief Set boost prefix (empty string to use default)
    bool set_boost_prefix(fs::path prefix);


    /// \brief Set command used by 'casm view'
    bool set_view_command(std::string opt);

    /// \brief Set crystallography tolerance
    bool set_crystallography_tol(double _tol);

    /// \brief Set linear algebra tolerance
    bool set_lin_alg_tol(double _tol);


    /// \brief Save settings to project settings file
    void commit() const;

    /// \brief Output as JSON
    jsonParser &to_json(jsonParser &json) const;

    /// \brief Print summary of ProjectSettings, as for 'casm settings -l'
    void print_summary(Log &log) const;

    /// \brief (deprecated) Set compile options to 'opt' (empty string to use default)
    bool set_compile_options(std::string opt);

    /// \brief (deprecated) Set shared library options to 'opt' (empty string to use default)
    bool set_so_options(std::string opt);



  private:

    /// \brief Changing the neighbor list properties requires updating Clexulator source code
    void _reset_clexulators();

    /// \brief initialize default compiler options
    void _load_default_options();


    DirectoryStructure m_dir;

    std::string m_name;

    // CASM project current settings

    // name : ClexDescription map
    std::map<std::string, ClexDescription> m_clex;

    // name of default cluster expansion
    std::string m_default_clex;

    // neighbor list settings
    Eigen::Matrix3l m_nlist_weight_matrix;
    std::set<int> m_nlist_sublat_indices;

    // Properties to read from calculations
    std::vector<std::string> m_properties;

    // Runtime library compilation settings: compilation options
    std::pair<std::string, std::string> m_cxx;
    std::pair<std::string, std::string> m_cxxflags;
    std::pair<std::string, std::string> m_soflags;
    std::pair<fs::path, std::string> m_casm_prefix;
    std::pair<fs::path, std::string> m_boost_prefix;

    // deprecated reading exactly from settings file
    std::string m_depr_compile_options;
    // deprecated reading exactly from settings file
    std::string m_depr_so_options;

    // Command executed by 'casm view'
    std::string m_view_command;

    // Crystallography tolerance
    double m_crystallography_tol;

    // Linear algebra tolerance
    double m_lin_alg_tol;

    // ConfigIO
    DataFormatterDictionary<Configuration> m_config_io_dict;

    // ConfigIO aliases
    std::map<std::string, std::string> m_aliases;

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
