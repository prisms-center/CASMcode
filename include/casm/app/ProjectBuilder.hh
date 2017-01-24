#ifndef CASM_ProjectBuilder
#define CASM_ProjectBuilder

#include <string>
#include <vector>

#include "casm/casm_io/jsonParser.hh"

#include "casm/app/ProjectSettings.hh"
#include "casm/app/AppIO.hh"

namespace CASM {

  /// \brief Sets up directories and files for a new CASM project
  ///
  /// \ingroup Project
  ///
  class ProjectBuilder {

  public:

    /// \brief Construct a CASM ProjectBuilder
    ///
    /// \param _root The directory where a new CASM project should be created.
    /// \param _name The name of the CASM project. Should be a short name suitable for prepending to files.
    /// \param _property The name of the default cluster expansion property, i.e. "formation_energy"
    ///
    ProjectBuilder(fs::path _root, std::string _name, std::string _property) :
      m_root(_root),
      m_name(_name),
      m_property(_property) {

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
    std::string m_property;

    // allow default initialization:

    std::vector<std::string> m_properties {"relaxed_energy"};
    double m_crystallography_tol = CASM::TOL;
    double m_lin_alg_tol = 1e-10;

  };

}

#endif
