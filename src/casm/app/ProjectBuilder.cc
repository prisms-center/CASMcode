#include "casm/app/ProjectBuilder.hh"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <regex>
#include "casm/app/AppIO.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ClexDescription.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/symmetry/SymGroup.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/database/DatabaseTypes.hh"

namespace CASM {

  namespace {
    struct SetPropsImpl {
      typedef std::map<std::string, std::vector<std::string> > map_type;
      ProjectSettings &set;
      map_type props;
      SetPropsImpl(ProjectSettings &_set, const map_type &_props):
        set(_set), props(_props) {}
      template<typename T> void eval() {
        set.properties<T>() = props[traits<T>::name];
      }
    };
  }

  bool ProjectBuilder::valid_title(std::string const &title) {
    return std::regex_match(title, std::regex(R"([_a-zA-Z]\w*)"));
  }

  /// \brief Construct a CASM ProjectBuilder
  ///
  /// \param _root The directory where a new CASM project should be created.
  /// \param _title The title of the CASM project. Should be a short title suitable for prepending to files.
  /// \param _property The name of the default cluster expansion property, i.e. "formation_energy"
  ///
  ProjectBuilder::ProjectBuilder(fs::path _root, std::string _title, std::string _property) :
    m_root(_root),
    m_title(_title),
    m_property(_property) {

    /// check if m_title is suitable:
    if(!valid_title(m_title)) {
      throw std::runtime_error(
        std::string("Error constructing ProjectBuilder.\n") +
        "  Invalid Project title: '" + m_title + "'\n"
        "  Must be a valid C++ identifier: \n"
        "  - only alphanumeric characters and underscores allowed\n"
        "  - cannot start with a number");
    }

    fs::path test_root = find_casmroot(m_root);
    if(test_root == m_root) {
      throw std::runtime_error(
        std::string("Error in 'ProjectBuilder::build()'.\n") +
        "  Already in a casm project: " + test_root.string());
    }

    // check for a prim.json
    fs::path prim_path = m_root / "prim.json";
    if(!fs::is_regular_file(prim_path)) {
      throw std::runtime_error(std::string("Error constructing 'ProjectBuilder'.\n") +
                               "  No prim.json file found at: " + prim_path.string());
    }

    // Read prim
    // HamiltonianModules needs to know where to look for plug-ins, but it's unclear where that should be
    HamiltonianModules hamiltonian_modules;
    m_prim = read_prim(prim_path, m_crystallography_tol, &hamiltonian_modules);
  }

  /// \brief Construct a CASM ProjectBuilder
  ///
  /// \param _root The directory where a new CASM project should be created.
  /// \param _title The title of the CASM project. Should be a short title suitable for prepending to files.
  /// \param _property The name of the default cluster expansion property, i.e. "formation_energy"
  ///
  ProjectBuilder::ProjectBuilder(BasicStructure const &_prim, fs::path _root, std::string _title, std::string _property) :
    m_prim(_prim),
    m_root(_root),
    m_title(_title),
    m_property(_property) {

    /// check if m_title is suitable:
    if(!valid_title(m_title)) {
      throw std::runtime_error(
        std::string("Error constructing ProjectBuilder.\n") +
        "  Invalid Project title: '" + m_title + "'\n"
        "  Must be a valid C++ identifier: \n"
        "  - only alphanumeric characters and underscores allowed\n"
        "  - cannot start with a number");
    }

    fs::path test_root = find_casmroot(m_root);
    if(test_root == m_root) {
      throw std::runtime_error(
        std::string("Error in 'ProjectBuilder::build()'.\n") +
        "  Already in a casm project: " + test_root.string());
    }

  }

  /// \brief Builds a new CASM project
  ///
  /// - requires a prim.json in the project directory
  ///   - does not perform any checks if structure is primitive, right-handed, reduced, etc.
  /// - creates directories
  /// - create project_settings, symmetry, and standard composition axes files
  ///
  ///
  void ProjectBuilder::build() const {

    try {

      DirectoryStructure dir(m_root);
      fs::create_directories(dir.casm_dir());
      write_prim(m_prim, dir.prim(), FRAC);

      ProjectSettings set(m_root, m_title);

      // create basic directories
      set.new_casm_dir();
      set.new_symmetry_dir();
      set.new_reports_dir();
      set.new_bset_dir("default");
      set.new_calc_settings_dir("default");
      set.new_ref_dir("default", "default");
      set.new_eci_dir(m_property, "default", "default", "default", "default");

      // set project settings

      auto exc = [ = ](std::string type) {
        throw std::runtime_error(
          std::string("Error in 'ProjectBuilder::build()'.\n") +
          "  Could not set " + type);
      };

      DB::for_each_config_type(SetPropsImpl(set, m_properties));

      ClexDescription desc("formation_energy", "formation_energy", "default", "default", "default", "default");
      set.new_clex(desc);
      set.set_default_clex(desc.name);

      if(!set.set_crystallography_tol(m_crystallography_tol)) {
        exc("crystallography_tol");
      }
      if(!set.set_lin_alg_tol(m_lin_alg_tol)) {
        exc("lin_alg_tol");
      }

      set.commit();

      Structure prim(m_prim);

      // Calculate symmetry  --------------------
      // get lattice point group and character table
      SymGroup lattice_point_grp(SymGroup::lattice_point_group(prim.lattice()));
      lattice_point_grp.character_table();

      // get factor group
      /* prim.generate_factor_group(); */
      /* prim.set_site_internals(); */

      // Write symmetry info files

      // Write lattice point group
      {
        fs::ofstream outfile;
        jsonParser json;
        outfile.open(dir.lattice_point_group());
        write_symgroup(lattice_point_grp, json);
        json.print(outfile);
        outfile.close();
      }

      // Write factor group
      {
        fs::ofstream outfile;
        jsonParser json;
        outfile.open(dir.factor_group());
        write_symgroup(prim.factor_group(), json);
        json.print(outfile);
        outfile.close();
      }

      // Write crystal point group
      {
        fs::ofstream outfile;
        jsonParser json;
        outfile.open(dir.crystal_point_group());
        write_symgroup(prim.point_group(), json);
        json.print(outfile);
        outfile.close();
      }


      // Generate empty composition_axes.json --------------------
      CompositionAxes().write(dir.composition_axes());

    }
    catch(...) {
      std::cerr << "Uncaught exception in ProjectBuilder::build()" << std::endl;
      /// re-throw exceptions
      throw;
    }
  }
}

