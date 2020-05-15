#include "casm/app/ProjectSettings.hh"

#include <regex>
#include <tuple>
#include "casm/system/RuntimeLibrary.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/app/EnumeratorHandler.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/app/HamiltonianModules.hh"
#include "casm/casm_io/SafeOfstream.hh"

namespace CASM {

  bool is_valid_project_name(std::string name) {
    /// check if m_name is suitable:
    if(!std::regex_match(name, std::regex(R"([_a-zA-Z]\w*)"))) {
      return false;
    }
    return true;
  }

  void check_project_name(std::string name) {
    /// check if m_name is suitable:
    if(!is_valid_project_name(name)) {
      throw std::runtime_error(
        std::string("  Invalid Project name: '") + name + "'\n"
        "  Must be a valid C++ identifier: \n"
        "  - only alphanumeric characters and underscores allowed\n"
        "  - cannot start with a number");
    }
  }

  ProjectSettings::ProjectSettings(std::string name, const Logging &logging) :
    Logging(logging),
    m_name(name),
    m_crystallography_tol(CASM::TOL),
    m_lin_alg_tol(1e-10),
    m_db_name("jsonDB") {

    check_project_name(name);
  }

  ProjectSettings::ProjectSettings(std::string name, fs::path root, const Logging &logging) :
    Logging(logging),
    m_name(name),
    m_dir(notstd::make_unique<DirectoryStructure>(root)),
    m_crystallography_tol(CASM::TOL),
    m_lin_alg_tol(1e-10),
    m_db_name("jsonDB") {

    check_project_name(name);

  }

  ProjectSettings::~ProjectSettings() {}

  std::string ProjectSettings::name() const {
    return m_name;
  }

  bool ProjectSettings::has_dir() const {
    return bool {m_dir};
  }

  DirectoryStructure const &ProjectSettings::dir() const {
    if(!m_dir) {
      throw std::runtime_error("Error accessing DirectoryStructure from ProjectSettings: Does not exist.");
    }
    return *m_dir;
  }

  bool ProjectSettings::set_root_dir(fs::path root) {
    fs::path checkroot = find_casmroot(root);
    if(checkroot == root) {
      throw std::runtime_error(std::string("Error in ProjectSettings::set_root_dir: Project already exists at '") + root.string() + "'");
    }
    m_dir = notstd::make_cloneable<DirectoryStructure>(root);
    return bool {m_dir};
  }

  fs::path ProjectSettings::root_dir() const {
    return dir().root_dir();
  }


  ProjectSettings::required_properties_map_type const &ProjectSettings::required_properties() const {
    return m_required_properties;
  }

  void ProjectSettings::set_required_properties(required_properties_map_type const &_required_properties) {
    m_required_properties = _required_properties;
  }

  std::vector<std::string> const &ProjectSettings::required_properties(std::string type_name, std::string calctype) const {
    return m_required_properties.find(type_name)->second.find(calctype)->second;
  }

  std::vector<std::string> &ProjectSettings::required_properties(std::string type_name, std::string calctype) {
    return m_required_properties[type_name][calctype];
  }


  std::map<std::string, ClexDescription> const &ProjectSettings::cluster_expansions() const {
    return m_clex;
  }

  bool ProjectSettings::has_clex(std::string clex_name) const {
    return m_clex.find(clex_name) != m_clex.end();
  }

  ClexDescription const &ProjectSettings::clex(std::string clex_name) const {
    return m_clex.find(clex_name)->second;
  }

  bool ProjectSettings::insert_clex(const ClexDescription &desc) {
    if(m_clex.find(desc.name) != m_clex.end()) {
      return false;
    }
    m_clex[desc.name] = desc;
    return true;
  }

  bool ProjectSettings::erase_clex(const ClexDescription &desc) {
    if(m_default_clex_name == desc.name) {
      m_default_clex_name.clear();
    }
    return m_clex.erase(desc.name);
  }

  std::string ProjectSettings::default_clex_name() const {
    return m_default_clex_name;
  }

  bool ProjectSettings::set_default_clex_name(const std::string &clex_name) {
    if(m_clex.find(clex_name) != m_clex.end()) {
      m_default_clex_name = clex_name;
      return true;
    }
    return false;
  }

  ClexDescription const &ProjectSettings::default_clex() const {
    return m_clex.find(m_default_clex_name)->second;
  }

  bool ProjectSettings::set_default_clex(const ClexDescription &desc) {
    m_clex[desc.name] = desc;
    m_default_clex_name = desc.name;
    return true;
  }


  std::pair<std::string, std::string> ProjectSettings::cxx() const {
    return m_cxx.first.empty() ? RuntimeLibrary::default_cxx() : m_cxx;
  }

  std::pair<std::string, std::string> ProjectSettings::cxxflags() const {
    return m_cxxflags.first.empty() ? RuntimeLibrary::default_cxxflags() : m_cxxflags;
  }

  std::pair<std::string, std::string> ProjectSettings::soflags() const {
    return m_soflags.first.empty() ? RuntimeLibrary::default_soflags() : m_soflags;
  }

  std::pair<fs::path, std::string> ProjectSettings::casm_includedir() const {
    return m_casm_includedir.first.empty() ? RuntimeLibrary::default_casm_includedir() : m_casm_includedir;
  }

  std::pair<fs::path, std::string> ProjectSettings::casm_libdir() const {
    return m_casm_libdir.first.empty() ? RuntimeLibrary::default_casm_libdir() : m_casm_libdir;
  }

  std::pair<fs::path, std::string> ProjectSettings::boost_includedir() const {
    return m_boost_includedir.first.empty() ? RuntimeLibrary::default_boost_includedir() : m_boost_includedir;
  }

  std::pair<fs::path, std::string> ProjectSettings::boost_libdir() const {
    return m_boost_libdir.first.empty() ? RuntimeLibrary::default_boost_libdir() : m_boost_libdir;
  }

  std::string ProjectSettings::compile_options() const {
    // construct from pieces
    return cxx().first + " " + cxxflags().first + " " +
           include_path(casm_includedir().first) + " " +
           include_path(boost_includedir().first);
  }

  std::string ProjectSettings::so_options() const {
    // construct from pieces
    return cxx().first + " " + soflags().first + " " +
           link_path(boost_libdir().first) + " " +
           link_path(casm_libdir().first);
  }

  std::string ProjectSettings::view_command() const {
    return m_view_command;
  }

  std::string ProjectSettings::view_command_video() const {
    return m_view_command_video;
  }

  double ProjectSettings::crystallography_tol() const {
    return m_crystallography_tol;
  }

  double ProjectSettings::lin_alg_tol() const {
    return m_lin_alg_tol;
  }


  EnumeratorHandler &ProjectSettings::enumerator_handler() {
    if(!m_enumerator_handler) {
      m_enumerator_handler = notstd::make_cloneable<EnumeratorHandler>(*this);
    }
    return *m_enumerator_handler;
  }

  EnumeratorHandler const &ProjectSettings::enumerator_handler() const {
    return const_cast<ProjectSettings &>(*this).enumerator_handler();
  }


  void ProjectSettings::set_db_name(std::string _db_name) {
    m_db_name = _db_name;
  }

  std::string ProjectSettings::db_name() const {
    return m_db_name;
  }

  ProjectSettings::query_alias_map_type const &ProjectSettings::query_alias() const {
    return m_query_alias;
  }

  void ProjectSettings::set_query_alias(query_alias_map_type const &_query_alias) {
    m_query_handler.clear();
    m_query_alias = _query_alias;
  }

  void ProjectSettings::set_query_alias(std::string type_name, std::string alias_name, std::string alias_value) {
    m_query_handler.erase(type_name);
    m_query_alias[type_name][alias_name] = alias_value;
  }

  template<typename DataObject>
  QueryHandler<DataObject> &ProjectSettings::query_handler() {
    auto res = m_query_handler.find(traits<DataObject>::name);
    if(res != m_query_handler.end()) {
      return static_cast<QueryHandler<DataObject>& >(*res->second);
    }

    res = m_query_handler.insert(
            std::make_pair(
              traits<DataObject>::name,
              notstd::cloneable_ptr<notstd::Cloneable>(new QueryHandler<DataObject>(*this))
            )
          ).first;

    QueryHandler<DataObject> &q = static_cast<QueryHandler<DataObject>& >(*res->second);

    // if query aliases exist, add all for this type
    auto type_it = m_query_alias.find(traits<DataObject>::name);
    if(type_it != m_query_alias.end()) {
      for(auto const &value : type_it->second) {
        q.add_alias(value.first, value.second);
      }
    }
    return q;
  }

  template<typename DataObject>
  QueryHandler<DataObject> const &ProjectSettings::query_handler() const {
    return const_cast<ProjectSettings &>(*this).query_handler<DataObject>();
  }


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


  std::string ProjectSettings::global_clexulator_name() const {
    return name() + "_Clexulator";
  }


  bool ProjectSettings::has_m_nlist_weight_matrix() const {
    return bool {m_nlist_weight_matrix};
  }

  Eigen::Matrix3l ProjectSettings::nlist_weight_matrix() const {
    return *m_nlist_weight_matrix;
  }

  bool ProjectSettings::set_nlist_weight_matrix(Eigen::Matrix3l M) {
    m_nlist_weight_matrix = notstd::clone(M);
    return bool {m_nlist_weight_matrix};
  }

  bool ProjectSettings::has_nlist_sublat_indices() const {
    return bool {m_nlist_sublat_indices};
  }

  const std::set<int> &ProjectSettings::nlist_sublat_indices() const {
    return *m_nlist_sublat_indices;
  }

  bool ProjectSettings::set_nlist_sublat_indices(std::set<int> value) {
    m_nlist_sublat_indices = notstd::clone(value);
    return bool {m_nlist_sublat_indices};
  }


  bool ProjectSettings::set_cxx(std::string opt) {
    m_cxx = std::make_pair(opt, "project_settings");
    return true;
  }

  bool ProjectSettings::set_cxxflags(std::string opt)  {
    m_cxxflags = std::make_pair(opt, "project_settings");
    return true;
  }

  bool ProjectSettings::set_soflags(std::string opt)  {
    m_soflags = std::make_pair(opt, "project_settings");
    return true;
  }

  bool ProjectSettings::set_casm_prefix(fs::path prefix)  {
    m_casm_includedir = std::make_pair(prefix / "include", "project_settings");
    m_casm_libdir = std::make_pair(prefix / "lib", "project_settings");
    return true;
  }

  bool ProjectSettings::set_casm_includedir(fs::path dir)  {
    m_casm_includedir = std::make_pair(dir, "project_settings");
    return true;
  }

  bool ProjectSettings::set_casm_libdir(fs::path dir)  {
    m_casm_libdir = std::make_pair(dir, "project_settings");
    return true;
  }

  bool ProjectSettings::set_boost_prefix(fs::path prefix)  {
    m_boost_includedir = std::make_pair(prefix / "include", "project_settings");
    m_boost_libdir = std::make_pair(prefix / "lib", "project_settings");
    return true;
  }

  bool ProjectSettings::set_boost_includedir(fs::path dir)  {
    m_boost_includedir = std::make_pair(dir, "project_settings");
    return true;
  }

  bool ProjectSettings::set_boost_libdir(fs::path dir)  {
    m_boost_libdir = std::make_pair(dir, "project_settings");
    return true;
  }


  bool ProjectSettings::set_view_command(std::string opt) {
    m_view_command = opt;
    return true;
  }

  bool ProjectSettings::set_view_command_video(std::string opt) {
    m_view_command_video = opt;
    return true;
  }

  bool ProjectSettings::set_crystallography_tol(double _tol) {
    m_crystallography_tol = _tol;
    return true;
  }

  bool ProjectSettings::set_lin_alg_tol(double _tol) {
    m_lin_alg_tol = _tol;
    return true;
  }


  bool create_all_directories(ProjectSettings const &set) {
    auto const &dir = set.dir();
    auto const &cluster_expansions = set.cluster_expansions();

    bool result {true};
    result &= dir.new_casm_dir();
    result &= dir.new_symmetry_dir();
    result &= dir.new_reports_dir();
    for(auto const &pair : cluster_expansions) {
      ClexDescription const &desc = pair.second;
      result &= new_dir(dir, desc);
    }
    return result;
  }

  namespace {
    std::string _wdefaultval(std::string name, std::pair<std::string, std::string> val) {
      return name + ": '" + val.first + "' (" + val.second + ")\n";
    }

    std::string _wdefaultval(std::string name, std::pair<fs::path, std::string> val) {
      return _wdefaultval(name, std::make_pair(val.first.string(), val.second));
    }
  }

  void print_compiler_settings_summary(ProjectSettings const &set, Log &log) {
    log.custom<Log::standard>("Compiler settings");
    log << _wdefaultval("cxx", set.cxx())
        << _wdefaultval("cxxflags", set.cxxflags())
        << _wdefaultval("soflags", set.soflags())
        << _wdefaultval("casm_includedir", set.casm_includedir())
        << _wdefaultval("casm_libdir", set.casm_libdir())
        << _wdefaultval("boost_includedir", set.boost_includedir())
        << _wdefaultval("boost_libdir", set.boost_libdir()) << std::endl;
    log << "compile command: '" << set.compile_options() << "'\n\n";
    log << "so command: '" << set.so_options() << "'\n\n";
  }

  void print_summary(ProjectSettings const &set, Log &log) {

    log.custom<Log::standard>("Cluster expansions");
    auto const &clex = set.cluster_expansions();
    auto const &default_clex_name = set.default_clex_name();
    for(auto it = clex.begin(); it != clex.end(); ++it) {

      const ClexDescription &desc = it->second;
      bool is_default = (desc.name == default_clex_name);
      int indent = 0;
      desc.print(log, is_default, indent);
    }
    log << std::endl;

    const ClexDescription &default_desc = set.default_clex();

    // print basis sets
    std::vector<std::string> all_bset = set.dir().all_bset();
    log.custom<Log::standard>("Basis sets");
    for(int i = 0; i < all_bset.size(); i++) {
      log << all_bset[i];
      if(all_bset[i] == default_desc.bset) {
        log << "*";
      }
      log << "\n";
    }
    log << std::endl;

    // print training data (calctype) & references
    std::vector<std::string> all_calctype = set.dir().all_calctype();
    log.custom<Log::standard>("Training Data & References");
    for(int i = 0; i < all_calctype.size(); i++) {
      std::vector<std::string> all_ref = set.dir().all_ref(all_calctype[i]);
      for(int j = 0; j < all_ref.size(); j++) {
        log << all_calctype[i] << " / " << all_ref[j];
        if(all_calctype[i] == default_desc.calctype && all_ref[j] == default_desc.ref) {
          log << "*";
        }
        log << "\n";
      }
      log << std::endl;
    }

    // print training data (calctype) & references
    std::vector<std::string> all_eci = set.dir().all_eci(default_desc.property, default_desc.calctype, default_desc.ref, default_desc.bset);
    log.custom<Log::standard>("ECI for current cluster expansion settings group");
    for(int i = 0; i < all_eci.size(); i++) {
      log << all_eci[i];
      if(all_eci[i] == default_desc.eci) {
        log << "*";
      }
      log << "\n";
    }
    log << std::endl;

    log << "*: indicates the default settings used by CASM whenever particular \n"
        "settings are not explicitly specified (i.e. the basis set to evaluate \n"
        "for 'casm query -k corr')\n\n";

    print_compiler_settings_summary(set, log);

    log.custom<Log::standard>("'casm view'");
    log << "command: '" << set.view_command() << "'\n\n";
    log << "video command: '" << set.view_command_video() << "'\n\n";

  }

  jsonParser &to_json(const ProjectSettings &set, jsonParser &json) {
    json = jsonParser::object();

    json["name"] = set.name();
    json["cluster_expansions"] = set.cluster_expansions();
    json["required_properties"] = set.required_properties();
    json["default_clex"] = set.default_clex_name();
    json["nlist_weight_matrix"] = set.nlist_weight_matrix();
    json["nlist_sublat_indices"] = set.nlist_sublat_indices();

    // only write compiler & linker settings if stored in project settings, not if obtained from
    // default values or environment variables

    auto _write_str_if = [&](std::string name, std::pair<std::string, std::string> value_and_source) {
      if(value_and_source.second == "project_settings") {
        json[name] = value_and_source.first;
      };
    };

    auto _write_path_if = [&](std::string name, std::pair<fs::path, std::string> value_and_source) {
      if(value_and_source.second == "project_settings") {
        json[name] = value_and_source.first.string();
      };
    };

    _write_str_if("cxx", set.cxx());
    _write_str_if("cxxflags", set.cxxflags());
    _write_str_if("soflags", set.soflags());
    _write_path_if("casm_includedir", set.casm_includedir());
    _write_path_if("casm_libdir", set.casm_libdir());
    _write_path_if("boost_includedir", set.boost_includedir());
    _write_path_if("boost_libdir", set.boost_libdir());

    json["view_command"] = set.view_command();
    json["view_command_video"] = set.view_command_video();
    json["crystallography_tol"] = set.crystallography_tol();
    json["crystallography_tol"].set_scientific();
    json["lin_alg_tol"] = set.lin_alg_tol();
    json["lin_alg_tol"].set_scientific();
    json["query_alias"] = set.query_alias();

    return json;
  }

  ProjectSettings jsonConstructor<ProjectSettings>::from_json(
    jsonParser const &json,
    Logging const &logging) {
    try {

      ProjectSettings settings {json["name"].get<std::string>(), logging};

      // read required_properties map
      if(json.contains("required_properties") && json["required_properties"].size()) {
        auto it = json.find("required_properties");
        ProjectSettings::required_properties_map_type tmp;
        CASM::from_json(tmp, *it);
        settings.set_required_properties(tmp);
      }

      // read cluster_expansions list and set default_clex
      if(json.contains("cluster_expansions") && json["cluster_expansions"].size()) {
        auto it = json.find("cluster_expansions");
        auto clex_it = it->begin();
        auto clex_end = it->end();
        for(; clex_it != clex_end; ++clex_it) {
          settings.insert_clex(clex_it->get<ClexDescription>());
        }

        if(json.contains("default_clex")) {
          settings.set_default_clex_name(json["default_clex"].get<std::string>());
        }
        else {
          if(settings.has_clex("formation_energy")) {
            settings.set_default_clex_name("formation_energy");
          }
          else {
            settings.set_default_clex_name(settings.cluster_expansions().begin()->first);
          }
        }
      }
      // if no cluster_expansions list, create with default_configuration_clex
      else {
        ClexDescription desc = default_configuration_clex();
        settings.insert_clex(desc);
        settings.set_default_clex_name(desc.name);
      }

      // read compiler & linker options

      std::string tmp_str;
      if(json.get_if(tmp_str, "cxx")) {
        settings.set_cxx(tmp_str);
      }
      if(json.get_if(tmp_str, "cxxflags")) {
        settings.set_cxxflags(tmp_str);
      }
      if(json.get_if(tmp_str, "soflags")) {
        settings.set_soflags(tmp_str);
      }

      fs::path tmp_path;
      if(json.get_if(tmp_path, "casm_prefix")) {
        settings.set_casm_prefix(tmp_path); // deprecated
      }
      if(json.get_if(tmp_path, "boost_prefix")) {
        settings.set_boost_prefix(tmp_path); // deprecated
      }
      if(json.get_if(tmp_path, "casm_includedir")) {
        settings.set_casm_includedir(tmp_path);
      }
      if(json.get_if(tmp_path, "casm_libdir")) {
        settings.set_casm_libdir(tmp_path);
      }
      if(json.get_if(tmp_path, "boost_includedir")) {
        settings.set_boost_includedir(tmp_path);
      }
      if(json.get_if(tmp_path, "boost_libdir")) {
        settings.set_boost_libdir(tmp_path);
      }

      // other options
      if(json.get_if(tmp_str, "view_command")) {
        settings.set_view_command(tmp_str);
      }
      if(json.get_if(tmp_str, "view_command_video")) {
        settings.set_view_command_video(tmp_str);
      }

      // precision options -- always set
      double tmp_double = TOL;
      json.get_if(tmp_double, "crystallography_tol");
      settings.set_crystallography_tol(tmp_double);

      json.get_else(tmp_double, "lin_alg_tol", 1e-10);
      settings.set_lin_alg_tol(tmp_double);

      // read nlist json -- no longer generate defaults
      if(json.contains("nlist_weight_matrix")) {
        settings.set_nlist_weight_matrix(json["nlist_weight_matrix"].get<Eigen::Matrix3l>());
      }
      if(json.contains("nlist_sublat_indices")) {
        settings.set_nlist_sublat_indices(json["nlist_sublat_indices"].get<std::set<int>>());
      }

      // read query aliases
      if(json.contains("query_alias") && json["query_alias"].size()) {
        auto it = json.find("query_alias");
        ProjectSettings::query_alias_map_type tmp;
        CASM::from_json(tmp, *it);
        settings.set_query_alias(tmp);
      }

      return settings;

    }
    catch(std::exception &e) {
      logging.err_log() << "Error reading ProjectSettings from JSON." << std::endl;
      throw e;
    }
  }

  void commit(ProjectSettings const &set) {
    jsonParser json;
    to_json(set, json);

    SafeOfstream file;
    file.open(set.dir().project_settings());
    json.print(file.ofstream(), 2, 18);
    file.close();
  }

  ProjectSettings read_project_settings(fs::path path, Logging const &logging) {
    try {
      fs::ifstream file {path};
      jsonParser json {file};
      return jsonConstructor<ProjectSettings>::from_json(json, logging);
    }
    catch(std::exception &e) {
      logging.err_log() << "Error reading ProjectSettings from: '" << path << "'" << std::endl;
      throw e;
    }
  }

}
