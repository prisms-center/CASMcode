#include "casm/app/ProjectBuilder.hh"
#include "casm/clex/NeighborList.hh"

namespace CASM {

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

      set.properties() = m_properties;

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
      Structure prim;
      fs::ifstream primfile(dir.prim());
      prim = Structure(read_prim(jsonParser(primfile)));
      primfile.close();


      // Calculate symmetry  --------------------
      // get lattice point group and character table
      SymGroup lattice_point_grp;
      prim.lattice().generate_point_group(lattice_point_grp, m_crystallography_tol);
      lattice_point_grp.character_table();

      // get factor group
      prim.generate_factor_group(m_crystallography_tol);
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


      // Generate standard composition axes --------------------
      CompositionAxes opt;

      opt.standard.clear();
      std::vector<CompositionConverter> v;
      standard_composition_axes(prim, std::back_inserter(v));
      for(int i = 0; i < v.size(); i++) {
        opt.standard[std::to_string(i)] = v[i];
      }

      opt.write(dir.composition_axes());

    }
    catch(...) {
      std::cerr << "Uncaught exception in ProjectBuilder::build()" << std::endl;
      /// re-throw exceptions
      throw;
    }
  }
}

