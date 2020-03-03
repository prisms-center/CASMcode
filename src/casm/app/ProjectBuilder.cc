#include "casm/app/ProjectBuilder.hh"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/regex.hpp>
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


  /// \brief Construct a CASM ProjectBuilder
  ///
  /// \param _root The directory where a new CASM project should be created.
  /// \param _name The name of the CASM project. Should be a short name suitable for prepending to files.
  /// \param _property The name of the default cluster expansion property, i.e. "formation_energy"
  ///
  ProjectBuilder::ProjectBuilder(fs::path _root, std::string _name, std::string _property) :
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

      fs::path test_root = find_casmroot(m_root);
      if(test_root == m_root) {
        std::cerr << "Attempting to create a casm project here: " << m_root << std::endl;
        std::cerr << "Found existing casm project here: " << test_root << std::endl;
        throw std::runtime_error(
          std::string("Error in 'ProjectBuilder::build()'.\n") +
          "  Already in a casm project: " + test_root.string());
      }

      DirectoryStructure dir(m_root);

      // check for a prim.json
      if(!fs::is_regular_file(dir.prim())) {
        throw std::runtime_error(
          std::string("Error in 'ProjectBuilder::build()'.\n") +
          "  No prim.json file found at: " + dir.prim().string());
      }

      ProjectSettings set(m_root, m_name);

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


      // Read prim
      fs::ifstream primfile(dir.prim());
      Structure prim(read_prim(jsonParser(primfile), m_crystallography_tol, &(set.hamiltonian_modules())));
      primfile.close();


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

