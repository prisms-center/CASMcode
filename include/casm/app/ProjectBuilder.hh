#ifndef CASM_ProjectBuilder
#define CASM_ProjectBuilder

#include <string>
#include <vector>

#include "casm/casm_io/jsonParser.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/app/AppIO.hh"

namespace CASM {

  /// \brief Sets up directories and files for a new CASM project
  class ProjectBuilder {

  public:

    /// \brief Construct a CASM ProjectBuilder
    ///
    /// \param _root The directory where a new CASM project should be created.
    /// \param _name The name of the CASM project. Should be a short name suitable for prepending to files.
    /// \param _cluster_expansion The name of the default cluster expansion, i.e. "formation_energy"
    ///
    ProjectBuilder(fs::path _root, std::string _name, std::string _cluster_expansion) :
      m_root(_root),
      m_name(_name),
      m_clex(_cluster_expansion) {

      /// check if m_name is suitable:
      if(!boost::regex_match(m_name, boost::regex(R"([_a-zA-Z]\w*)"))) {
        throw std::runtime_error(
          std::string("Error constructing ProjectBuilder.\n") +
          "  Invalid Project name: '" + m_name + "'\n"
          "  Must be a valid C++ identifier: \n"
          "  - only alphanumeric characters and underscores allowed\n"
          "  - cannot start with a number");
      }

    }


    ProjectBuilder &set_bset(std::string _bset) {
      m_bset = _bset;
      return *this;
    }

    ProjectBuilder &set_calctype(std::string _calctype) {
      m_calctype = _calctype;
      return *this;
    }

    ProjectBuilder &set_ref(std::string _ref) {
      m_ref = _ref;
      return *this;
    }

    ProjectBuilder &set_eci(std::string _eci) {
      m_eci = _eci;
      return *this;
    }

    ProjectBuilder &set_compile_options(std::string _compile_options) {
      m_compile_options = _compile_options;
      return *this;
    }

    ProjectBuilder &set_so_options(std::string _so_options) {
      m_so_options = _so_options;
      return *this;
    }

    ProjectBuilder &set_crystallography_tol(double _tol) {
      m_crystallography_tol = _tol;
      return *this;
    }
    
    ProjectBuilder &set_lin_alg_tol(double _tol) {
      m_lin_alg_tol = _tol;
      return *this;
    }

    /// \brief Builds a new CASM project
    void build() const;


  private:

    // require user initialization:

    fs::path m_root;
    std::string m_name;
    std::string m_clex;

    // allow default initialization:

    std::vector<std::string> m_properties {"relaxed_energy"};
    std::string m_bset = "default";
    std::string m_calctype = "default";
    std::string m_ref = "default";
    std::string m_eci = "default";
    std::string m_compile_options = RuntimeLibrary::default_compile_options();
    std::string m_so_options = RuntimeLibrary::default_so_options() + " -lboost_system";
    double m_crystallography_tol = CASM::TOL;
    double m_lin_alg_tol = 1e-10;

  };

}

#endif
