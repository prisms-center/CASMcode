#ifndef CASM_ProjectBuilder
#define CASM_ProjectBuilder

#include <string>
#include <vector>
#include <boost/filesystem.hpp>

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
      m_clex(_cluster_expansion) {}


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

    ProjectBuilder &set_tol(double _tol) {
      m_tol = _tol;
      return *this;
    }

    /// \brief Builds a new CASM project
    ///
    /// - requires a prim.json in the project directory
    ///   - does not perform any checks if structure is primitive, right-handed, reduced, etc.
    /// - creates directories
    /// - create project_settings, symmetry, and standard composition axes files
    ///
    ///
    void build() const {

      try {

        fs::path test_root = find_casmroot(m_root);
        if(test_root == m_root) {
          throw std::runtime_error(
            std::string("Error in 'ProjectBuilder::build()'.\n") +
            "  Already in a casm project: " + test_root.string());
        }

        ProjectSettings set(m_root, m_name);
        DirectoryStructure dir(m_root);

        // check for a prim.json
        if(!fs::is_regular_file(dir.prim())) {
          throw std::runtime_error(
            std::string("Error in 'ProjectBuilder::build()'.\n") +
            "  No prim.json file found at: " + dir.prim().string());
        }

        // create basic directories
        set.new_casm_dir();
        set.new_symmetry_dir();
        set.new_bset_dir(m_bset);
        set.new_calc_settings_dir(m_calctype);
        set.new_ref_dir(m_calctype, m_ref);
        set.new_eci_dir(m_clex, m_calctype, m_ref, m_bset, m_eci);

        // set project settings

        auto exc = [ = ](std::string type) {
          throw std::runtime_error(
            std::string("Error in 'ProjectBuilder::build()'.\n") +
            "  Could not set " + type);
        };

        set.properties() = m_properties;

        if(!set.set_bset(m_bset)) {
          exc("basis_set");
        }
        if(!set.set_calctype(m_calctype)) {
          exc("calctype");
        }
        if(!set.set_ref(m_calctype, m_ref)) {
          exc("ref");
        }
        if(!set.set_clex(m_clex)) {
          exc("clex");
        }
        if(!set.set_eci(m_clex, m_calctype, m_ref, m_bset, m_eci)) {
          exc("eci");
        }
        if(!set.set_compile_options(m_compile_options)) {
          exc("compile options");
        }
        if(!set.set_so_options(m_so_options)) {
          exc("so options");
        }
        if(!set.set_tol(m_tol)) {
          exc("tol");
        }

        set.commit();

        // Calculate symmetry and default composition axes
        // Generated project objects
        Structure prim;
        SymGroup lattice_point_grp;

        // read prim
        fs::ifstream primfile(dir.prim());
        prim = Structure(read_prim(jsonParser(primfile)));
        primfile.close();

        // get lattice point group and character table
        prim.lattice().generate_point_group(lattice_point_grp, m_tol);
        lattice_point_grp.get_character_table();

        // get factor group
        prim.generate_factor_group(m_tol);
        prim.set_site_internals();


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

        // Generate standard composition axes...
        CompositionAxes opt;

        opt.standard.clear();
        std::vector<CompositionConverter> v;
        standard_composition_axes(prim, std::back_inserter(v));
        for(int i = 0; i < v.size(); i++) {
          opt.standard[std::to_string(i)] = v[i];
        }

        opt.write(dir.composition_axes(set.calctype(), set.ref()));

      }
      catch(...) {
        std::cerr << "Uncaught exception in ProjectBuilder::build()" << std::endl;
        /// re-throw exceptions
        throw;
      }
    }


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
    double m_tol = CASM::TOL;

  };

}

#endif
