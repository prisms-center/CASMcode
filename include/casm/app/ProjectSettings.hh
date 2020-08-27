#ifndef CASM_ProjectSettings
#define CASM_ProjectSettings

#include <string>
#include <vector>
#include <map>
#include <set>

#include "casm/global/eigen.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {

  class Configuration;
  class jsonParser;
  class EnumeratorHandler;
  template<typename T> class QueryHandler;
  class HamiltonianModules;

  /** \defgroup Project
   *
   *  Relates to CASM project settings, directory structure, etc.
   *
   *  A CASM project encompasses all the settings, calculations, cluster
   *  expansions, Monte Carlo results, etc. related to a single parent
   *  crystal structure (known as the 'prim').
   *
   *  All the results and data related to a CASM project may be stored in a
   *  directory tree with structure defined by the DirectoryStructure class, and
   *  accessible through the top-level data structure PrimClex.
   *
   *  @{
  */

  /// Return true if project name is valid
  ///
  /// Notes:
  /// - Should be a short descriptive name such as "ZrO" or "NiAl"
  /// - Must be a valid C++ identifier:
  ///   - only alphanumeric characters and underscores allowed
  ///   - cannot start with a number
  bool is_valid_project_name(std::string project_name);

  /// Throw if project name is invalid
  void throw_if_project_name_is_not_valid(std::string project_name);

  /// CASM project settings
  ///
  /// Contains:
  /// - The project name
  /// - Parameters that control which basis sites are included in the neighborlist and the order of
  ///   sites in the neighborlist
  /// - Compiler and linker settings for Clexulator
  /// - Crystallography and linear algebra default tolerances
  /// - Names of properties that must be calculated, by object type, and calctype
  /// - List of cluster expansions (via a vector of ClexDescription, which is key type that allows
  ///   looking up cluster expansion data)
  /// - A default cluster expansion (default ClexDescription)
  /// - Subcommands to call visualization programs
  /// - Default database name
  /// - Query aliases
  ///
  /// Also contains:
  /// - EnumeratorHandler: a dictionary-like class for holding standard enumeration methods and plugins
  /// - HamiltonianModules: a dictionary-like class for holding standard and (todo) custom
  ///    traits for anisotropic value types, dof type, and symmetry representations
  /// - QueryHandler: a dictionary-like class for holding query functors, plugins, and aliases for
  ///    querying object properties and formatted output
  ///
  /// Optionally contains:
  /// - The associated project root directory path and a DirectoryStructure object which,
  ///   if present, allows creating project directories
  ///
  /// Compiler and linker settings:
  /// - For easier debugging, each setting is returned as a pair of value and the source where the
  /// value was obtained from. This could be "project_settings" for values set explicitly in this
  /// object or read from the project settings JSON file, or it could be an environment variable,
  /// or it could be a default value.
  /// - The environment variables checked and default values are: (using implementation in RuntimeLibrary.cc)
  ///   - cxx: default="g++", check for environment variables: `CASM_CXX`, `CXX`
  ///   - cxxflags: default="-O3 -Wall -fPIC --std=c++11", check for environment variable: `CASM_CXXFLAGS`
  ///   - soflags: default="-shared -lboost_system", check for environment variable: `CASM_SOFLAGS`
  ///   - casm_includedir: (where to find the "casm" headers directory tree)
  ///     - default=(attempts to find in standard locations relative to the `ccasm` executable),
  ///       - for example: if <prefix>/bin/ccasm exists then check <prefix>/include/casm
  ///     - check for environment variables:
  ///       - `CASM_INCLUDEDIR`
  ///       - `CASM_PREFIX` (to indicate `CASM_PREFIX/include`)
  ///   - casm_libdir: (where to find libcasm)
  ///     - default=(attempts to find in standard locations relative to the `ccasm` executable),
  ///       - for example: if <prefix>/bin/ccasm exists then check <prefix>/{lib,lib64,lib/x86_64-linux-gnu}libcasm.{dylib,so}
  ///     - check for environment variables:
  ///       - `CASM_LIBDIR`
  ///       - `CASM_PREFIX` (to indicate `CASM_PREFIX/lib`)
  ///   - boost_includedir: (where to find the "boost" headers directory tree)
  ///     - default=(attempts to find in standard locations relative to the `ccasm` executable),
  ///     - check for environment variables:
  ///       - `CASM_BOOST_INCLUDEDIR`
  ///       - `CASM_BOOST_PREFIX` (to indicate `CASM_BOOST_PREFIX/include`)
  ///   - boost_libdir: (where to find boost libraries)
  ///     - default=(attempts to find in standard locations relative to the `ccasm` executable),
  ///     - check for environment variables:
  ///       - `CASM_BOOST_LIBDIR`
  ///       - `CASM_BOOST_PREFIX` (to indicate `CASM_BOOST_PREFIX/lib`)
  ///
  /// Notes on previous versions:
  /// - Prior to v0.4: JSON file attribute "curr_properties" was expected to be an array of string
  ///   containing properties that must exist in properties.calc.json for a configuration to be
  //    considered calculated.
  class ProjectSettings {

  public:

    /// Construct CASM project settings for a new project (memory only)
    ///
    /// \param project_name Name of new CASM project. Use a short title suitable for prepending to file names.
    ///
    explicit ProjectSettings(std::string project_name);

    /// Construct CASM project settings for a new project (to be persisted on disk)
    ///
    /// \param project_name Name of new CASM project. Use a short title suitable for prepending to file names.
    /// \param root Path to new CASM project directory
    ///
    explicit ProjectSettings(std::string project_name, fs::path root);

    ~ProjectSettings();


    /// Get project name
    std::string project_name() const;


    // --- optional: root_dir & DirectoryStructure ---

    /// Check if DirectoryStructure exists
    bool has_dir() const;

    /// Access DirectoryStructure object. Throw if not set.
    DirectoryStructure const &dir() const;

    /// Set DirectoryStructure
    bool set_root_dir(fs::path root);

    /// Access dir().root_dir(). Throw if not set.
    fs::path root_dir() const;


    // --- required properties settings ---

    typedef std::string ObjectTypeName;
    typedef std::string CalcTypeName;

    /// required_properties is a map of properties names required for a calculation to be complete:
    ///   traits<Configuration>::name -> calctype -> {propname1, propname2, ...}
    typedef std::map<ObjectTypeName, std::map<CalcTypeName, std::vector<std::string>>>
    required_properties_map_type;

    /// Access properties required for an object to be considered calculated
    ///
    /// Note:
    /// - required_properties is a map of properties names required for a calculation to be complete:
    ///   traits<ObjectType>::name -> calctype -> {propname1, propname2, ...}
    required_properties_map_type const &required_properties() const;

    /// Access properties required for an object to be considered calculated
    ///
    /// Note:
    /// - required_properties is a map of properties names required for a calculation to be complete:
    ///   traits<ObjectType>::name -> calctype -> {propname1, propname2, ...}
    void set_required_properties(required_properties_map_type const &_required_properties);

    /// const Access properties required for an object to be considered calculated
    ///
    /// \param type_name string, i.e. traits<Configuration>::name
    /// \param calctype string
    std::vector<std::string> const &required_properties(std::string type_name, std::string calctype) const;

    /// Access properties required for an object to be considered calculated
    ///
    /// \param type_name string, i.e. traits<Configuration>::name
    /// \param calctype string
    void set_required_properties(std::string type_name, std::string calctype, std::vector<std::string> const &_required_properties);

    // --- cluster expansion settings: ClexDescription, by ClexDescription.name ---

    /// Const access map of all ClexDescription
    std::map<std::string, ClexDescription> const &cluster_expansions() const;

    /// Check if a ClexDescription exists
    bool has_clex(std::string clex_name) const;

    /// Get a ClexDescription by name
    ClexDescription const &clex(std::string clex_name) const;

    /// Insert a ClexDescription into ProjectSettings. Fails if existing ClexDescription with same name
    bool insert_clex(ClexDescription const &desc);

    /// Erase existing ClexDescription from ProjectSettings
    ///
    /// If `desc` is the default clex, the default clex is cleared.
    bool erase_clex(ClexDescription const &desc);

    /// Get default ClexDescription name
    std::string default_clex_name() const;

    /// Set default ClexDescription by name
    bool set_default_clex_name(std::string const &clex_name);

    /// Get default ClexDescription
    ClexDescription const &default_clex() const;

    /// Set default ClexDescription. Will overwrite existing ClexDescription with same name
    bool set_default_clex(ClexDescription const &desc);


    // --- compiler and linker settings ---

    /// Get c++ compiler (pair of value and source for the value)
    std::pair<std::string, std::string> cxx() const;

    /// Get c++ compiler options (pair of value and source for the value)
    std::pair<std::string, std::string> cxxflags() const;

    /// Get shared object options (pair of value and source for the value)
    std::pair<std::string, std::string> soflags() const;

    /// Get casm includedir (pair of value and source for the value)
    std::pair<fs::path, std::string> casm_includedir() const;

    /// Get casm libdir (pair of value and source for the value)
    std::pair<fs::path, std::string> casm_libdir() const;

    /// Get boost includedir (pair of value and source for the value)
    std::pair<fs::path, std::string> boost_includedir() const;

    /// Get boost libdir (pair of value and source for the value)
    std::pair<fs::path, std::string> boost_libdir() const;

    /// Combines cxx, cxxflags, casm_includedir, and boost_includedir to make compiler options string
    std::string compile_options() const;

    /// Combines cxx, soflags, casm_libdir, and boost_libdir to make shared library options string
    std::string so_options() const;


    /// Get current command used by 'casm view'
    std::string view_command() const;

    /// Get current video viewing command used by 'casm view'
    std::string view_command_video() const;


    /// Get current project crystallography tolerance
    double crystallography_tol() const;

    /// Get current project linear algebra tolerance
    double lin_alg_tol() const;


    // ** Enumerators **

    EnumeratorHandler &enumerator_handler();

    EnumeratorHandler const &enumerator_handler() const;


    // ** Database **

    /// Set default database type name (for future)
    void set_default_database_name(std::string _default_database_name);

    /// Get default database type name
    std::string default_database_name() const;


    // ** Querie aliases **

    typedef std::string QueryAliasName;
    typedef std::string QueryAliasValue;

    /// query_alias is a map by object type name of query alias {name,value} pairs:
    ///   traits<Configuration>::name -> {query alias name, query alias value}
    typedef std::map<ObjectTypeName, std::map<QueryAliasName, QueryAliasValue>> query_alias_map_type;

    /// Access query aliases stored in the project settings
    ///
    /// Note:
    /// - query_alias is a map, by object type name, of query alias {name,value} pairs:
    ///   traits<Configuration>::name -> {query alias name, query alias value}
    query_alias_map_type const &query_alias() const;

    /// Set query aliases to be stored in the project settings
    ///
    /// Note:
    /// - query_alias is a map, by object type name, of query alias {name,value} pairs:
    ///   traits<Configuration>::name -> {query alias name, query alias value}
    void set_query_alias(query_alias_map_type const &_query_alias);

    /// Set query alias. Invalidates query_handler references.
    ///
    /// Note:
    /// - query_alias is a map, by object type name, of query alias {name,value} pairs:
    ///   traits<Configuration>::name -> {query alias name, query alias value}
    void set_query_alias(std::string type_name, std::string alias_name, std::string alias_value);

    template<typename DataObject>
    QueryHandler<DataObject> &query_handler();

    template<typename DataObject>
    QueryHandler<DataObject> const &query_handler() const;



    // ** Hamiltonian Modules **

    HamiltonianModules &hamiltonian_modules();

    HamiltonianModules const &hamiltonian_modules()const;


    // ** Clexulator names **

    /// Name to use for clexulator printing
    std::string global_clexulator_name() const;


    /// Check if neighbor list weight matrix exists
    bool has_m_nlist_weight_matrix() const;

    /// Get neighbor list weight matrix
    Eigen::Matrix3l nlist_weight_matrix() const;

    /// Set neighbor list weight matrix (will delete existing Clexulator source and compiled code)
    ///
    /// Note:
    /// - Clexulators are dependent on the neighbor list parameters. If you changes these parameters
    ///   you must regenerate all clexulators. This is no longer done by this function.
    bool set_nlist_weight_matrix(Eigen::Matrix3l M);

    /// Check if set of sublattice indices to include in neighbor lists exists
    bool has_nlist_sublat_indices() const;

    /// Get set of sublattice indices to include in neighbor lists
    std::set<int> const &nlist_sublat_indices() const;

    /// Set range of sublattice indices to include in neighbor lists
    ///
    /// Note:
    /// - Clexulators are dependent on the neighbor list parameters. If you changes these parameters
    ///   you must regenerate all clexulators. This is no longer done by this function.
    bool set_nlist_sublat_indices(std::set<int> value);


    /// Set c++ compiler (empty string to use default)
    bool set_cxx(std::string opt);

    /// Set c++ compiler options (empty string to use default)
    bool set_cxxflags(std::string opt);

    /// Set shared object options (empty string to use default)
    bool set_soflags(std::string opt);


    /// Set casm prefix (empty string to use default)
    bool set_casm_prefix(fs::path dir);

    /// Set casm includedir (empty string to use default)
    bool set_casm_includedir(fs::path dir);

    /// Set casm libdir (empty string to use default)
    bool set_casm_libdir(fs::path dir);


    /// Set boost prefix (empty string to use default)
    bool set_boost_prefix(fs::path dir);

    /// Set boost includedir (empty string to use default)
    bool set_boost_includedir(fs::path dir);

    /// Set boost libdir (empty string to use default)
    bool set_boost_libdir(fs::path dir);


    /// Set command used by 'casm view'
    bool set_view_command(std::string opt);

    /// Set video viewing command used by 'casm view'
    bool set_view_command_video(std::string opt);


    /// Set crystallography tolerance
    bool set_crystallography_tol(double _tol);

    /// Set linear algebra tolerance
    bool set_lin_alg_tol(double _tol);


  private:

    // project name
    std::string m_project_name;

    // directory structure (may be empty)
    notstd::cloneable_ptr<DirectoryStructure> m_dir;

    // EnumeratorHandler (type erased)
    notstd::cloneable_ptr<notstd::Cloneable> m_enumerator_handler;

    // Datatype name : QueryHandler<DataType> map (type erased)
    std::map<std::string, notstd::cloneable_ptr<notstd::Cloneable> > m_query_handler;

    // traits<ObjecType>::name -> std::pair{alias name, alias value}
    query_alias_map_type m_query_alias;

    // HamiltonianModules (type erased)
    mutable notstd::cloneable_ptr<notstd::Cloneable> m_hamiltonian_modules;

    // CASM project current settings

    // name : ClexDescription map
    std::map<std::string, ClexDescription> m_clex;

    // name of default cluster expansion
    std::string m_default_clex_name;

    // neighbor list parameters (may be empty)
    notstd::cloneable_ptr<Eigen::Matrix3l> m_nlist_weight_matrix;
    notstd::cloneable_ptr<std::set<int>> m_nlist_sublat_indices;

    // Properties required to be read from calculations
    // traits<ObjecType>::name -> calctype -> {prop1, prop2, ...}
    required_properties_map_type m_required_properties;

    // Runtime library compilation settings: compilation options
    std::pair<std::string, std::string> m_cxx;
    std::pair<std::string, std::string> m_cxxflags;
    std::pair<std::string, std::string> m_soflags;
    std::pair<fs::path, std::string> m_casm_includedir;
    std::pair<fs::path, std::string> m_casm_libdir;
    std::pair<fs::path, std::string> m_boost_includedir;
    std::pair<fs::path, std::string> m_boost_libdir;

    // Command executed by 'casm view'
    std::string m_view_command;

    // Command executed by 'casm view'
    std::string m_view_command_video;

    // Crystallography tolerance
    double m_crystallography_tol;

    // Linear algebra tolerance
    double m_lin_alg_tol;

    // Database
    std::string m_default_database_name;
  };


  /// Add directories for all cluster expansions
  ///
  /// For all ClexDescription currently saved in `this->cluster_expansions()`,
  ///   create bset, calc_settings, ref, and eci directories
  bool create_all_directories(ProjectSettings const &set);

  /// Print summary of compiler settings, as for 'casm settings -l'
  void print_compiler_settings_summary(ProjectSettings const &set, Log &log);

  /// Print summary of ProjectSettings, as for 'casm settings -l'
  void print_summary(ProjectSettings const &set, Log &log);

  /// Serialize ProjectSettings to JSON
  jsonParser &to_json(const ProjectSettings &set, jsonParser &json);

  template<typename T> struct jsonConstructor;

  template<>
  struct jsonConstructor<ProjectSettings> {

    /// Construct ProjectSettings from JSON
    ///
    /// \param json Input JSON
    ///
    /// Note:
    /// - After construction, root dir is not set (`ProjectSettings::has_dir() == false`)
    static ProjectSettings from_json(jsonParser const &json);
  };

  /// Write ProjectSettings to JSON file
  ///
  /// \param project_settings_path Location of ProjectSettings JSON file
  /// \param logging For showing errors and warnings
  void write_project_settings(ProjectSettings const &set, fs::path project_settings_path);

  /// Read ProjectSettings from JSON file
  ///
  /// \param project_settings_path Location of ProjectSettings JSON file
  /// \param logging For showing errors and warnings
  ///
  /// Note:
  /// - After construction, root dir is not set (`ProjectSettings::has_dir() == false`)
  ProjectSettings read_project_settings(fs::path project_settings_path);

  /// Write project settings file for an existing project
  ///
  /// Equivalent to `write_project_settings(set, set.dir().project_settings())`
  void commit(ProjectSettings const &set);

  /// Open the project settings JSON file of an existing project
  ///
  /// Note:
  /// - After construction, root dir is set (`ProjectSettings::has_dir() == true`)
  /// - `path_in_project` does not need to be a root dir exactly, it may be anywere inside a CASM project
  /// - Throws if `path_in_project` is not inside any CASM project or reading fails
  ProjectSettings open_project_settings(fs::path path_in_project);

  /** @} */
}

#endif
