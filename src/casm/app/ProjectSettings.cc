#include "casm/app/ProjectSettings.hh"

#include <tuple>
#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/EnumeratorHandler.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/app/HamiltonianModules.hh"
//#include "casm/basis_set/DoF.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/misc/algorithm.hh"
#include "casm/casm_io/SafeOfstream.hh"
#include "casm/database/DatabaseTypes.hh"

namespace CASM {

  /// \brief Default weight matrix for approximately spherical neighborhood in Cartesian coordinates
  ///
  /// Equivalent to:
  /// \code
  /// PrimNeighborList::make_weight_matrix(prim.lattice().lat_column_mat(), 10, tol());
  /// \endcode
  Eigen::Matrix3l _default_nlist_weight_matrix(const Structure &prim, double tol) {
    return PrimNeighborList::make_weight_matrix(prim.lattice().lat_column_mat(), 10, tol);
  }

  /// \brief Default includes sublattices with >= 2 components
  std::set<int> _default_nlist_sublat_indices(const Structure &prim) {
    // for now, include sublattices with >= 2 components
    std::set<int> sublat_indices;
    for(int b = 0; b < prim.basis().size(); ++b) {
      if(prim.basis()[b].occupant_dof().size() >= 2 || prim.basis()[b].dof_size() > 0) {
        sublat_indices.insert(b);
      }
    }
    return sublat_indices;
  }

  namespace {
    struct AddAliasesFromJSON {
      AddAliasesFromJSON(ProjectSettings &_set, jsonParser &_json) :
        set(_set), json(_json) {}

      template<typename T>
      void eval() {
        if(json.contains(traits<T>::name)) {
          auto &q = set.query_handler<T>();
          auto res = json.find(traits<T>::name);
          auto it = res->begin();
          auto end = res->end();
          for(; it != end; ++it) {
            q.add_alias(it.name(), it->template get<std::string>());
          }
        }
      }

      ProjectSettings &set;
      jsonParser &json;
    };

    struct AddAliasesToJSON {
      AddAliasesToJSON(const ProjectSettings &_set, jsonParser &_json) : json(_json), set(_set) {}

      template<typename T> void eval() {
        json[traits<T>::name] = set.query_handler<T>().aliases();
      }

      jsonParser &json;
      const ProjectSettings &set;
    };
  }

  /// \brief Default constructor
  ProjectSettings::ProjectSettings() :
    m_crystallography_tol(CASM::TOL),
    m_lin_alg_tol(1e-10) {}


  /// \brief Construct CASM project settings for a new project
  ///
  /// \param root Path to new CASM project directory
  /// \param name Name of new CASM project. Use a short title suitable for prepending to file names.
  ///
  ProjectSettings::ProjectSettings(fs::path root, std::string name, const Logging &logging) :
    Logging(logging),
    m_dir(root),
    m_name(name),
    m_db_name("jsonDB") {

    if(fs::exists(m_dir.casm_dir())) {
      throw std::runtime_error(
        std::string("Error in 'ProjectSettings(fs::path root, std::string name)'.\n") +
        "  A CASM project already exists at this location: '" + root.string() + "'");
    }

    // check for a prim.json
    if(!fs::is_regular_file(m_dir.prim())) {
      throw std::runtime_error(
        std::string("Error in 'ProjectSettings(fs::path root, std::string name)'.\n") +
        "  No prim.json file found at: " + m_dir.prim().string());
    }

    // generate default nlist settings
    Structure prim(read_prim(m_dir.prim(), hamiltonian_modules(), TOL));
    m_nlist_weight_matrix = _default_nlist_weight_matrix(prim, TOL);
    m_nlist_sublat_indices = _default_nlist_sublat_indices(prim);
  }

  /// \brief Construct CASM project settings from existing project
  ///
  /// \param root Path to existing CASM project directory. Project settings will be read.
  ///
  ProjectSettings::ProjectSettings(fs::path root, const Logging &logging) :
    Logging(logging),
    m_dir(root),
    m_db_name("jsonDB") {

    if(fs::exists(m_dir.casm_dir())) {

      try {

        // read .casmroot current settings
        fs::ifstream file(m_dir.project_settings());
        jsonParser settings(file);

        if(settings.contains("curr_properties")) { // deprecated after v0.2.X
          from_json(m_properties, settings["curr_properties"]);
        }
        else { // v0.3+
          from_json(m_properties, settings["properties"]);
        }

        if(settings.contains("cluster_expansions") && settings["cluster_expansions"].size()) {
          from_json(m_clex, settings["cluster_expansions"]);

          // if no 'default_clex', use 'formation_energy' if that exists, else use first in clex list
          if(!settings.get_if(m_default_clex, "default_clex")) {
            if(m_clex.find("formation_energy") != m_clex.end()) {
              m_default_clex = "formation_energy";
            }
            else {
              m_default_clex = m_clex.begin()->first;
            }
          }
        }
        else {
          ClexDescription desc("formation_energy", "formation_energy", "default", "default", "default", "default");
          m_clex[desc.name] = desc;
          m_default_clex = desc.name;
        }

        auto _read_if = [&](std::pair<std::string, std::string> &opt, std::string name) {
          if(settings.get_if(opt.first, name)) {
            opt.second = "project_settings";
          }
        };

        auto _read_path_if = [&](std::pair<fs::path, std::string> &opt, std::string name) {
          if(settings.get_if(opt.first, name)) {
            opt.second = "project_settings";
          }
        };

        _read_if(m_cxx, "cxx");
        _read_if(m_cxxflags, "cxxflags");
        _read_if(m_soflags, "soflags");

        fs::path tmp;
        if(settings.get_if(tmp, "casm_prefix")) {
          m_casm_includedir.first = tmp / "include";
          m_casm_includedir.second = "project_settings";
          m_casm_libdir.first = tmp / "lib";
          m_casm_libdir.second = "project_settings";
        }
        if(settings.get_if(tmp, "boost_prefix")) {
          m_boost_includedir.first = tmp / "include";
          m_boost_includedir.second = "project_settings";
          m_boost_libdir.first = tmp / "lib";
          m_boost_libdir.second = "project_settings";
        }

        _read_path_if(m_casm_includedir, "casm_includedir");
        _read_path_if(m_casm_libdir, "casm_libdir");
        _read_path_if(m_boost_includedir, "boost_includedir");
        _read_path_if(m_boost_libdir, "boost_libdir");

        settings.get_if(m_depr_compile_options, "compile_options");
        settings.get_if(m_depr_so_options, "so_options");


        // other options
        settings.get_if(m_view_command, "view_command");
        settings.get_if(m_view_command_video, "view_command_video");
        from_json(m_name, settings["name"]);

        // precision options
        settings.get_else(m_crystallography_tol, "tol", TOL); // deprecated 'tol'
        settings.get_if(m_crystallography_tol, "crystallography_tol");

        settings.get_else(m_lin_alg_tol, "lin_alg_tol", 1e-10);

        // read nlist settings, or generate defaults
        Structure prim;
        bool and_commit = false;
        if(!settings.contains("nlist_weight_matrix") || !settings.contains("nlist_sublat_indices")) {
          _reset_clexulators();
          prim = Structure(read_prim(m_dir.prim(), hamiltonian_modules(), crystallography_tol()));
          and_commit = true;
        }

        if(settings.contains("nlist_weight_matrix")) {
          from_json(m_nlist_weight_matrix, settings["nlist_weight_matrix"]);
        }
        else {
          m_nlist_weight_matrix = _default_nlist_weight_matrix(prim, crystallography_tol());
        }

        if(settings.contains("nlist_sublat_indices")) {
          from_json(m_nlist_sublat_indices, settings["nlist_sublat_indices"]);
          //if(m_nlist_sublat_indices.size() && (*m_nlist_sublat_indices.rbegin())>=prim.basis().size()){
          //std::cout << "ABOUT TO BAIL: " <<  *m_nlist_sublat_indices.rbegin() << " vs " << prim.basis().size() << std::endl;
          //throw std::runtime_error("JSON settings field \"nlist_sublat_indices\" is incompatible with prim.json. This may occur if prim.json is edited without re-initializing project.");
          //}
        }
        else {
          m_nlist_sublat_indices = _default_nlist_sublat_indices(prim);
        }

        // add aliases to dictionary
        DB::for_each_type(AddAliasesFromJSON(*this, settings["query_alias"]));

        if(and_commit) {
          commit();
        }

      }
      catch(std::exception &e) {
        std::cerr << "Error in ProjectSettings::ProjectSettings(const fs::path root).\n" <<
                  "  Error reading " << m_dir.project_settings() << std::endl;
        throw e;
      }
    }
    else if(!fs::exists(m_dir.casm_dir())) {
      throw std::runtime_error(
        std::string("Error in 'ProjectSettings(fs::path root, std::string name)'.\n") +
        "  A CASM project does not exist at this location: '" + root.string() + "'");
    }

  }

  ProjectSettings::~ProjectSettings() {}

  /// \brief Get project name
  std::string ProjectSettings::name() const {
    return m_name;
  }

  /// \brief Access current properties
  template<typename DataObject>
  std::vector<std::string> &ProjectSettings::properties() {
    std::string config_type_name = traits<DataObject>::name;
    if(!DB::config_types().count(config_type_name)) {
      std::stringstream msg;
      msg << "Unrecognized config type: " << config_type_name;
      throw std::runtime_error(msg.str());
    }
    return m_properties[config_type_name];
  }


  const std::map<std::string, ClexDescription> &ProjectSettings::cluster_expansions() const {
    return m_clex;
  }

  bool ProjectSettings::has_clex(std::string name) const {
    return m_clex.find(name) != m_clex.end();
  }

  const ClexDescription &ProjectSettings::clex(std::string name) const {
    return m_clex.find(name)->second;
  }

  bool ProjectSettings::new_clex(const ClexDescription &desc) {
    if(m_clex.find(desc.name) != m_clex.end()) {
      return false;
    }
    m_clex[desc.name] = desc;
    return true;
  }

  bool ProjectSettings::erase_clex(const ClexDescription &desc) {
    if(cluster_expansions().size() == 1) {
      return false;
    }

    if(m_default_clex == desc.name) {
      for(auto it = m_clex.begin(); it != m_clex.end(); ++it) {
        if(it->first != desc.name) {
          m_default_clex = it->first;
          break;
        }
      }
    }

    return m_clex.erase(desc.name);
  }

  const ClexDescription &ProjectSettings::default_clex() const {
    return m_clex.find(m_default_clex)->second;
  }

  bool ProjectSettings::set_default_clex(const std::string &clex_name) {
    if(m_clex.find(clex_name) != m_clex.end()) {
      m_default_clex = clex_name;
      return true;
    }
    return false;
  }

  /// \brief Will overwrite existing ClexDescription with same name
  bool ProjectSettings::set_default_clex(const ClexDescription &desc) {
    m_clex[desc.name] = desc;
    m_default_clex = desc.name;
    return true;
  }


  /// \brief Get neighbor list weight matrix
  Eigen::Matrix3l ProjectSettings::nlist_weight_matrix() const {
    return m_nlist_weight_matrix;
  }

  /// \brief Get set of sublattice indices to include in neighbor lists
  const std::set<int> &ProjectSettings::nlist_sublat_indices() const {
    return m_nlist_sublat_indices;
  }

  /// \brief Get c++ compiler
  std::pair<std::string, std::string> ProjectSettings::cxx() const {
    return m_cxx.first.empty() ? RuntimeLibrary::default_cxx() : m_cxx;
  }

  /// \brief Get c++ compiler options
  std::pair<std::string, std::string> ProjectSettings::cxxflags() const {
    return m_cxxflags.first.empty() ? RuntimeLibrary::default_cxxflags() : m_cxxflags;
  }

  /// \brief Get shared object options
  std::pair<std::string, std::string> ProjectSettings::soflags() const {
    return m_soflags.first.empty() ? RuntimeLibrary::default_soflags() : m_soflags;
  }

  /// \brief Get casm includedir
  std::pair<fs::path, std::string> ProjectSettings::casm_includedir() const {
    return m_casm_includedir.first.empty() ? RuntimeLibrary::default_casm_includedir() : m_casm_includedir;
  }

  /// \brief Get casm libdir
  std::pair<fs::path, std::string> ProjectSettings::casm_libdir() const {
    return m_casm_libdir.first.empty() ? RuntimeLibrary::default_casm_libdir() : m_casm_libdir;
  }

  /// \brief Get boost includedir
  std::pair<fs::path, std::string> ProjectSettings::boost_includedir() const {
    return m_boost_includedir.first.empty() ? RuntimeLibrary::default_boost_includedir() : m_boost_includedir;
  }

  /// \brief Get boost libdir
  std::pair<fs::path, std::string> ProjectSettings::boost_libdir() const {
    return m_boost_libdir.first.empty() ? RuntimeLibrary::default_boost_libdir() : m_boost_libdir;
  }

  /// \brief Get current compilation options string
  std::string ProjectSettings::compile_options() const {
    if(!m_depr_compile_options.empty()) {
      return m_depr_compile_options;
    }
    else {
      // else construct from pieces
      return cxx().first + " " + cxxflags().first + " " +
             include_path(casm_includedir().first) + " " +
             include_path(boost_includedir().first);
    }
  }

  /// \brief Get current shared library options string
  std::string ProjectSettings::so_options() const {
    // default to read deprecated 'so_options'
    if(!m_depr_so_options.empty()) {
      return m_depr_so_options;
    }
    else {
      // else construct from pieces
      return cxx().first + " " + soflags().first + " " +
             link_path(boost_libdir().first) + " " +
             link_path(casm_libdir().first);
    }
  }

  /// \brief Get current command used by 'casm view'
  std::string ProjectSettings::view_command() const {
    return m_view_command;
  }

  /// \brief Get current command used by 'casm view'
  std::string ProjectSettings::view_command_video() const {
    return m_view_command_video;
  }

  /// \brief Get current project crystallography tolerance
  double ProjectSettings::crystallography_tol() const {
    return m_crystallography_tol;
  }

  /// \brief Get current project linear algebra tolerance
  double ProjectSettings::lin_alg_tol() const {
    return m_lin_alg_tol;
  }


  // ** Enumerators **

  EnumeratorHandler &ProjectSettings::enumerator_handler() {
    if(!m_enumerator_handler) {
      m_enumerator_handler = notstd::make_cloneable<EnumeratorHandler>(*this);
    }
    return *m_enumerator_handler;
  }

  const EnumeratorHandler &ProjectSettings::enumerator_handler() const {
    return const_cast<ProjectSettings &>(*this).enumerator_handler();
  }


  // ** Database **

  void ProjectSettings::set_db_name(std::string _db_name) {
    m_db_name = _db_name;
  }

  std::string ProjectSettings::db_name() const {
    return m_db_name;
  }

  // ** Queries **

  template<typename DataObject>
  QueryHandler<DataObject> &ProjectSettings::query_handler() {
    auto res = m_query_handler.find(traits<DataObject>::name);
    if(res == m_query_handler.end()) {
      res = m_query_handler.insert(
              std::make_pair(
                traits<DataObject>::name,
                notstd::cloneable_ptr<notstd::Cloneable>(new QueryHandler<DataObject>(*this))
              )
            ).first;
    }
    return static_cast<QueryHandler<DataObject>& >(*res->second);
  }

  template<typename DataObject>
  const QueryHandler<DataObject> &ProjectSettings::query_handler() const {
    return const_cast<ProjectSettings &>(*this).query_handler<DataObject>();
  }

  // ** Hamiltonian Modules **

  HamiltonianModules &ProjectSettings::hamiltonian_modules() {
    if(!m_hamiltonian_modules) {
      m_hamiltonian_modules = notstd::make_cloneable<HamiltonianModules>(this);
    }
    return *m_hamiltonian_modules;
  }

  HamiltonianModules const &ProjectSettings::hamiltonian_modules()const {
    if(!m_hamiltonian_modules) {
      m_hamiltonian_modules = notstd::make_cloneable<HamiltonianModules>(this);
    }
    return *m_hamiltonian_modules;
  }

  // ** Clexulator names **

  std::string ProjectSettings::global_clexulator_name() const {
    return name() + "_Clexulator";
  }


  // ** Add directories for additional project data **

  /// \brief Create new project data directory
  bool ProjectSettings::new_casm_dir() const {
    return fs::create_directory(m_dir.casm_dir());
  }

  /// \brief Create new symmetry directory
  bool ProjectSettings::new_symmetry_dir() const {
    return fs::create_directory(m_dir.symmetry_dir());
  }

  /// \brief Create new reports directory
  bool ProjectSettings::new_reports_dir() const {
    return fs::create_directory(m_dir.reports_dir());
  }
  /// \brief Add a basis set directory
  bool ProjectSettings::new_bset_dir(std::string bset) const {
    return fs::create_directories(m_dir.bset_dir(bset));
  }

  /// \brief Add a cluster expansion directory
  bool ProjectSettings::new_clex_dir(std::string clex) const {
    return fs::create_directories(m_dir.clex_dir(clex));
  }


  /// \brief Add calculation settings directory path
  bool ProjectSettings::new_calc_settings_dir(std::string calctype) const {
    return fs::create_directories(m_dir.calc_settings_dir(calctype));
  }

  /// \brief Add calculation settings directory path, for supercell specific settings
  bool ProjectSettings::new_supercell_calc_settings_dir(std::string scelname, std::string calctype) const {
    return fs::create_directories(m_dir.supercell_calc_settings_dir(scelname, calctype));
  }

  /// \brief Add calculation settings directory path, for configuration specific settings
  bool ProjectSettings::new_configuration_calc_settings_dir(std::string configname, std::string calctype) const {
    return fs::create_directories(m_dir.configuration_calc_settings_dir(configname, calctype));
  }


  /// \brief Add a ref directory
  bool ProjectSettings::new_ref_dir(std::string calctype, std::string ref) const {
    return fs::create_directories(m_dir.ref_dir(calctype, ref));
  }

  /// \brief Add an eci directory
  bool ProjectSettings::new_eci_dir(std::string property, std::string calctype, std::string ref, std::string bset, std::string eci) const {
    return fs::create_directories(m_dir.eci_dir(property, calctype, ref, bset, eci));
  }


  // ** Change current settings **

  /// \brief const Access current properties
  template<typename DataObject>
  const std::vector<std::string> &ProjectSettings::properties() const {
    return const_cast<ProjectSettings &>(*this).properties<DataObject>();
  }


  /// \brief Set neighbor list weight matrix (will delete existing Clexulator
  /// source and compiled code)
  bool ProjectSettings::set_nlist_weight_matrix(Eigen::Matrix3l M) {

    // changing the neighbor list properties requires updating Clexulator source code
    // for now we just remove existing source/compiled files
    _reset_clexulators();

    m_nlist_weight_matrix = M;
    return true;
  }


  /// \brief Set c++ compiler (empty string to use default)
  bool ProjectSettings::set_cxx(std::string opt) {
    m_cxx = std::make_pair(opt, "project_settings");
    return true;
  }

  /// \brief Set c++ compiler options (empty string to use default)
  bool ProjectSettings::set_cxxflags(std::string opt)  {
    m_cxxflags = std::make_pair(opt, "project_settings");
    return true;
  }

  /// \brief Set shared object options (empty string to use default)
  bool ProjectSettings::set_soflags(std::string opt)  {
    m_soflags = std::make_pair(opt, "project_settings");
    return true;
  }

  /// \brief Set casm prefix (empty string to use default)
  bool ProjectSettings::set_casm_prefix(fs::path prefix)  {
    m_casm_includedir = std::make_pair(prefix / "include", "project_settings");
    m_casm_libdir = std::make_pair(prefix / "lib", "project_settings");
    return true;
  }

  /// \brief Set casm includedir (empty string to use default)
  bool ProjectSettings::set_casm_includedir(fs::path dir)  {
    m_casm_includedir = std::make_pair(dir, "project_settings");
    return true;
  }

  /// \brief Set casm libdir (empty string to use default)
  bool ProjectSettings::set_casm_libdir(fs::path dir)  {
    m_casm_libdir = std::make_pair(dir, "project_settings");
    return true;
  }

  /// \brief Set boost prefix (empty string to use default)
  bool ProjectSettings::set_boost_prefix(fs::path prefix)  {
    m_boost_includedir = std::make_pair(prefix / "include", "project_settings");
    m_boost_libdir = std::make_pair(prefix / "lib", "project_settings");
    return true;
  }

  /// \brief Set boost includedir (empty string to use default)
  bool ProjectSettings::set_boost_includedir(fs::path dir)  {
    m_boost_includedir = std::make_pair(dir, "project_settings");
    return true;
  }

  /// \brief Set boost libdir (empty string to use default)
  bool ProjectSettings::set_boost_libdir(fs::path dir)  {
    m_boost_libdir = std::make_pair(dir, "project_settings");
    return true;
  }

  /// \brief (deprecated) Set compile options to 'opt' (empty string to use default)
  bool ProjectSettings::set_compile_options(std::string opt) {
    m_depr_compile_options = opt;
    return true;
  }

  /// \brief (deprecated) Set shared library options to 'opt' (empty string to use default)
  bool ProjectSettings::set_so_options(std::string opt) {
    m_depr_so_options = opt;
    return true;
  }


  /// \brief Set command used by 'casm view'
  bool ProjectSettings::set_view_command(std::string opt) {
    m_view_command = opt;
    return true;
  }

  /// \brief Set command used by 'casm view'
  bool ProjectSettings::set_view_command_video(std::string opt) {
    m_view_command_video = opt;
    return true;
  }

  /// \brief Set crystallography tolerance
  bool ProjectSettings::set_crystallography_tol(double _tol) {
    m_crystallography_tol = _tol;
    return true;
  }

  /// \brief Set linear algebra tolerance
  bool ProjectSettings::set_lin_alg_tol(double _tol) {
    m_lin_alg_tol = _tol;
    return true;
  }

  /// \brief Save settings to file
  void ProjectSettings::commit() const {

    try {
      SafeOfstream file;
      file.open(m_dir.project_settings());
      jsonParser json;
      to_json(json).print(file.ofstream(), 2, 18);
      file.close();
    }
    catch(...) {
      std::cerr << "Uncaught exception in ProjectSettings::commit()" << std::endl;
      /// re-throw exceptions
      throw;
    }

  }

  /// \brief Changing the neighbor list properties requires updating Clexulator
  /// source code. This will remove existing source/compiled files
  void ProjectSettings::_reset_clexulators() {
    auto all_bset = m_dir.all_bset();
    for(auto it = all_bset.begin(); it != all_bset.end(); ++it) {
      fs::remove(m_dir.clexulator_src(name(), *it));
      fs::remove(m_dir.clexulator_o(name(), *it));
      fs::remove(m_dir.clexulator_so(name(), *it));
    }
  }

  void ProjectSettings::_load_default_options() {
    m_cxx = RuntimeLibrary::default_cxx();
    m_cxxflags = RuntimeLibrary::default_cxxflags();
    m_casm_includedir = RuntimeLibrary::default_casm_includedir();
    m_casm_libdir = RuntimeLibrary::default_casm_libdir();
    m_boost_includedir = RuntimeLibrary::default_boost_includedir();
    m_boost_libdir = RuntimeLibrary::default_boost_libdir();
    m_soflags = RuntimeLibrary::default_soflags();
  }

  jsonParser &ProjectSettings::to_json(jsonParser &json) const {

    json = jsonParser::object();

    json["name"] = name();
    json["cluster_expansions"] = cluster_expansions();
    json["properties"] = m_properties;
    json["default_clex"] = m_default_clex;
    json["nlist_weight_matrix"] = nlist_weight_matrix();
    json["nlist_sublat_indices"] = nlist_sublat_indices();


    auto _write_if = [&](std::string name, std::string val) {
      if(!val.empty()) {
        json[name] = val;
      };
    };

    _write_if("cxx", m_cxx.first);
    _write_if("cxxflags", m_cxxflags.first);
    _write_if("soflags", m_soflags.first);
    _write_if("casm_includedir", m_casm_includedir.first.string());
    _write_if("casm_libdir", m_casm_libdir.first.string());
    _write_if("boost_includedir", m_boost_includedir.first.string());
    _write_if("boost_libdir", m_boost_libdir.first.string());
    _write_if("compile_options", m_depr_compile_options);
    _write_if("so_options", m_depr_so_options);

    json["view_command"] = view_command();
    json["view_command_video"] = view_command_video();
    json["crystallography_tol"] = crystallography_tol();
    json["crystallography_tol"].set_scientific();
    json["lin_alg_tol"] = lin_alg_tol();
    json["lin_alg_tol"].set_scientific();

    DB::for_each_type(AddAliasesToJSON(*this, json["query_alias"]));

    return json;
  }

  namespace {
    std::string _wdefaultval(std::string name, std::pair<std::string, std::string> val) {
      return name + ": '" + val.first + "' (" + val.second + ")\n";
    }

    std::string _wdefaultval(std::string name, std::pair<fs::path, std::string> val) {
      return _wdefaultval(name, std::make_pair(val.first.string(), val.second));
    }
  }

  void ProjectSettings::print_compiler_settings_summary(Log &log) const {
    log.custom<Log::standard>("Compiler settings");
    log << _wdefaultval("cxx", cxx())
        << _wdefaultval("cxxflags", cxxflags())
        << _wdefaultval("soflags", soflags())
        << _wdefaultval("casm_includedir", casm_includedir())
        << _wdefaultval("casm_libdir", casm_libdir())
        << _wdefaultval("boost_includedir", boost_includedir())
        << _wdefaultval("boost_libdir", boost_libdir()) << std::endl;

    if(!m_depr_compile_options.empty()) {
      log << "Note: using deprecated 'compile_options' value from .casm/project_settings.json \n"
          "explicitly instead of individual compiler settings (cxx, cxxflags, casm_includedir,\n"
          "boost_includedir).\n"
          "Delete 'compile_options' from .casm/project_settings.json manually \n"
          "to use begin using the individually set settings.\n";
    }
    log << "compile command: '" << compile_options() << "'\n\n";

    if(!m_depr_so_options.empty()) {
      log << "Note: using deprecated 'so_options' value from .casm/project_settings.json \n"
          "explicitly instead of individual compiler settings (cxx, cxxflags, boost_libdir).\n"
          "Delete 'so_options' from .casm/project_settings.json manually \n"
          "to use begin using the individually set settings.\n";
    }
    log << "so command: '" << so_options() << "'\n\n";
  }

  /// \brief Print summary of ProjectSettings, as for 'casm settings -l'
  void ProjectSettings::print_summary(Log &log) const {

    log.custom<Log::standard>("Cluster expansions");
    for(auto it = m_clex.begin(); it != m_clex.end(); ++it) {

      const ClexDescription &desc = it->second;
      bool is_default = (desc.name == m_default_clex);
      int indent = 0;
      desc.print(log, is_default, indent);
    }
    log << std::endl;

    const ClexDescription &default_desc = default_clex();

    std::vector<std::string> all = m_dir.all_bset();
    log.custom<Log::standard>("Basis sets");
    for(int i = 0; i < all.size(); i++) {
      log << all[i];
      if(all[i] == default_desc.bset) {
        log << "*";
      }
      log << "\n";
    }
    log << std::endl;

    all = m_dir.all_calctype();
    log.custom<Log::standard>("Training Data & References");
    for(int i = 0; i < all.size(); i++) {
      std::vector<std::string> all_ref = m_dir.all_ref(all[i]);
      for(int j = 0; j < all_ref.size(); j++) {
        log << all[i] << " / " << all_ref[j];
        if(all[i] == default_desc.calctype && all_ref[j] == default_desc.ref) {
          log << "*";
        }
        log << "\n";
      }
      log << std::endl;
    }

    all = m_dir.all_eci(default_desc.property, default_desc.calctype, default_desc.ref, default_desc.bset);
    log.custom<Log::standard>("ECI for current cluster expansion settings group");
    for(int i = 0; i < all.size(); i++) {
      log << all[i];
      if(all[i] == default_desc.eci) {
        log << "*";
      }
      log << "\n";
    }
    log << std::endl;

    log << "*: indicates the default settings used by CASM whenever particular \n"
        "settings are not explicitly specified (i.e. the basis set to evaluate \n"
        "for 'casm query -k corr')\n\n";

    print_compiler_settings_summary(log);

    log.custom<Log::standard>("'casm view'");
    log << "command: '" << view_command() << "'\n\n";
    log << "video command: '" << view_command_video() << "'\n\n";

  }

  jsonParser &to_json(const ProjectSettings &set, jsonParser &json) {
    return set.to_json(json);
  }

}

// explicit instantiations

#define INST_ProjectSettings_all(r, data, type) \
template QueryHandler<type> &ProjectSettings::query_handler<type>(); \
template const QueryHandler<type> &ProjectSettings::query_handler<type>() const;

#define INST_ProjectSettings_config(r, data, type) \
template std::vector<std::string> &ProjectSettings::properties<type>(); \
template const std::vector<std::string> &ProjectSettings::properties<type>() const;

namespace CASM {
  BOOST_PP_SEQ_FOR_EACH(INST_ProjectSettings_all, _, CASM_DB_TYPES)
  BOOST_PP_SEQ_FOR_EACH(INST_ProjectSettings_config, _, CASM_DB_CONFIG_TYPES)
}
