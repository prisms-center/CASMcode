#ifndef CASM_ProjectSettings
#define CASM_ProjectSettings

#include <string>
#include <vector>
#include <boost/filesystem.hpp>

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
    explicit ProjectSettings(fs::path root, std::string name) :
      m_dir(root),
      m_name(name) {

      if(fs::exists(m_dir.casm_dir())) {
        throw std::runtime_error(
          std::string("Error in 'ProjectSettings(fs::path root, std::string name)'.\n") +
          "  A CASM project already exists at this location: '" + root.string() + "'");
      }

    }

    /// \brief Construct CASM project settings from existing project
    ///
    /// \param root Path to existing CASM project directory. Project settings will be read.
    ///
    explicit ProjectSettings(fs::path root) :
      m_dir(root) {

      if(fs::exists(m_dir.casm_dir())) {
        try {

          // read .casmroot current settings
          fs::ifstream file(m_dir.project_settings());
          jsonParser settings(file);

          from_json(m_properties, settings["curr_properties"]);
          from_json(m_bset, settings["curr_bset"]);
          from_json(m_calctype, settings["curr_calctype"]);
          from_json(m_ref, settings["curr_ref"]);
          from_json(m_clex, settings["curr_clex"]);
          from_json(m_eci, settings["curr_eci"]);
          settings.get_else(m_compile_options, "compile_options", RuntimeLibrary::default_compile_options());
          settings.get_else(m_so_options, "so_options", RuntimeLibrary::default_so_options());
          from_json(m_name, settings["name"]);
          from_json(m_tol, settings["tol"]);
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


    /// \brief Get project name
    std::string name() const {
      return m_name;
    }

    /// \brief Get current properties
    const std::vector<std::string> &properties() const {
      return m_properties;
    }

    /// \brief Get current basis set name
    std::string bset() const {
      return m_bset;
    }

    /// \brief Get current calctype name
    std::string calctype() const {
      return m_calctype;
    }

    /// \brief Get current ref name
    std::string ref() const {
      return m_ref;
    }

    /// \brief Get current cluster expansion name
    std::string clex() const {
      return m_clex;
    }

    /// \brief Get current eci name
    std::string eci() const {
      return m_eci;
    }

    /// \brief Get current compilation options string
    std::string compile_options() const {
      return m_compile_options;
    }

    /// \brief Get current shared library options string
    std::string so_options() const {
      return m_so_options;
    }

    /// \brief Get current project tol
    double tol() const {
      return m_tol;
    }


    // ** Clexulator names **

    std::string global_clexulator() const {
      return name() + "_Clexulator";
    }


    // ** Add directories for additional project data **

    /// \brief Create new project data directory
    bool new_casm_dir() const {
      return fs::create_directory(m_dir.casm_dir());
    }

    /// \brief Create new symmetry directory
    bool new_symmetry_dir() const {
      return fs::create_directory(m_dir.symmetry_dir());
    }

    /// \brief Add a basis set directory
    bool new_bset_dir(std::string bset) const {
      return fs::create_directories(m_dir.bset_dir(bset));
    }

    /// \brief Add a cluster expansion directory
    bool new_clex_dir(std::string clex) const {
      return fs::create_directories(m_dir.clex_dir(clex));
    }


    /// \brief Add calculation settings directory path
    bool new_calc_settings_dir(std::string calctype) const {
      return fs::create_directories(m_dir.calc_settings_dir(calctype));
    }

    /// \brief Add calculation settings directory path, for supercell specific settings
    bool new_supercell_calc_settings_dir(std::string scelname, std::string calctype) const {
      return fs::create_directories(m_dir.supercell_calc_settings_dir(scelname, calctype));
    }

    /// \brief Add calculation settings directory path, for configuration specific settings
    bool new_configuration_calc_settings_dir(std::string configname, std::string calctype) const {
      return fs::create_directories(m_dir.configuration_calc_settings_dir(configname, calctype));
    }


    /// \brief Add a ref directory
    bool new_ref_dir(std::string calctype, std::string ref) const {
      return fs::create_directories(m_dir.ref_dir(calctype, ref));
    }

    /// \brief Add calculation settings directory path, for supercell specific settings
    bool new_supercell_ref_dir(std::string scelname, std::string calctype, std::string ref) const {
      return fs::create_directories(m_dir.supercell_ref_dir(scelname, calctype, ref));
    }

    /// \brief Add calculation settings directory path, for configuration specific settings
    bool new_configuration_ref_dir(std::string configname, std::string calctype, std::string ref) const {
      return fs::create_directories(m_dir.configuration_ref_dir(configname, calctype, ref));
    }


    /// \brief Add an eci directory
    bool new_eci_dir(std::string clex, std::string calctype, std::string ref, std::string bset, std::string eci) const {
      return fs::create_directories(m_dir.eci_dir(clex, calctype, ref, bset, eci));
    }


    // ** Change current settings **

    /// \brief Access current properties
    std::vector<std::string> &properties() {
      return m_properties;
    }

    /// \brief Set current basis set to 'bset', if 'bset' exists
    bool set_bset(std::string bset) {
      auto all = m_dir.all_bset();
      if(std::find(all.begin(), all.end(), bset) != all.end()) {
        m_bset = bset;
        return true;
      }
      return false;
    }

    /// \brief Set current calctype to 'calctype', if 'calctype' exists
    bool set_calctype(std::string calctype) {
      auto all = m_dir.all_calctype();
      if(std::find(all.begin(), all.end(), calctype) != all.end()) {
        m_calctype = calctype;
        return true;
      }
      return false;
    }

    /// \brief Set current calculation reference to 'ref', if 'ref' exists
    bool set_ref(std::string calctype, std::string ref) {
      auto all = m_dir.all_ref(calctype);
      if(std::find(all.begin(), all.end(), ref) != all.end()) {
        m_ref = ref;
        return true;
      }
      return false;
    }

    /// \brief Set current cluster expansion to 'clex', if 'clex' exists
    bool set_clex(std::string clex) {
      auto all = m_dir.all_clex();
      if(std::find(all.begin(), all.end(), clex) != all.end()) {
        m_clex = clex;
        return true;
      }
      return false;
    }

    /// \brief Set current eci to 'eci', if 'eci' exists
    bool set_eci(std::string clex, std::string calctype, std::string ref, std::string bset, std::string eci) {
      auto all = m_dir.all_eci(clex, calctype, ref, bset);
      if(std::find(all.begin(), all.end(), eci) != all.end()) {
        m_eci = eci;
        return true;
      }
      return false;
    }

    /// \brief Set compile options to 'opt'
    bool set_compile_options(std::string opt) {
      m_compile_options = opt;
      return true;
    }

    /// \brief Set shared library options to 'opt'
    bool set_so_options(std::string opt) {
      m_so_options = opt;
      return true;
    }

    /// \brief Set shared library options to 'opt'
    bool set_tol(double _tol) {
      m_tol = _tol;
      return true;
    }


    /// \brief Save settings to project settings file
    void commit() const;


  private:

    DirectoryStructure m_dir;

    std::string m_name;

    // CASM project current settings: used to determine where to write things
    std::string m_bset;
    std::string m_calctype;
    std::string m_ref;
    std::string m_clex;
    std::string m_eci;

    // Properties to read from calculations
    std::vector<std::string> m_properties;

    // Runtime library compilation settings: compilation options
    std::string m_compile_options;
    std::string m_so_options;

    // Default tolerance
    double m_tol;

  };

  inline jsonParser &to_json(const ProjectSettings &set, jsonParser &json) {

    json = jsonParser::object();

    json["name"] = set.name();
    json["curr_properties"] = set.properties();
    json["curr_clex"] = set.clex();
    json["curr_calctype"] = set.calctype();
    json["curr_ref"] = set.ref();
    json["curr_bset"] = set.bset();
    json["curr_eci"] = set.eci();
    json["compile_options"] = set.compile_options();
    json["so_options"] = set.so_options();
    json["tol"] = set.tol();

    return json;
  }


  /// \brief Save settings to file
  inline void ProjectSettings::commit() const {

    try {
      SafeOfstream file;
      file.open(m_dir.project_settings());

      jsonParser json;
      to_json(*this, json);

      json.print(file.ofstream());
      file.close();
    }
    catch(...) {
      std::cerr << "Uncaught exception in ProjectSettings::commit()" << std::endl;
      /// re-throw exceptions
      throw;
    }

  }

}

#endif
