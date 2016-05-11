#include "Common.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/app/ProjectBuilder.hh"

namespace test {

  /// \brief Build a CASM project at 'proj_dir/title' using the prim
  void make_project(const Proj &proj) {

    jsonParser json;
    write_prim(proj.prim, json, FRAC);
    json["description"] = proj.desc;

    fs::create_directory(proj.dir);

    json.write(proj.dir / "prim.json");

    // build a project
    ProjectBuilder builder(proj.dir, proj.title, "formation_energy");
    builder.build();

  }

  /// \brief Remove a CASM project, checking first that there is a '.casm' dir
  ///
  /// Be careful! This does a recursive remove of the entire proj_dir!
  void rm_project(const Proj &proj) {

    if(fs::exists(proj.dir / ".casm")) {
      fs::remove_all(proj.dir);
    }

  }

  /// \brief Check some aspects of a SymGroup json
  void check_symgroup(const jsonParser &json, int N_op, int N_class) {
    BOOST_CHECK_EQUAL(json["character_table"].size(), N_class);
    BOOST_CHECK_EQUAL(json["conjugacy_class"].size(), N_class);

    BOOST_CHECK_EQUAL(json["symop"].size(), N_op);
    BOOST_CHECK_EQUAL(json["symop"][0]["type"].get<std::string>(), "identity");

    BOOST_CHECK_EQUAL(json["inverse"].size(), N_op);
    BOOST_CHECK_EQUAL(json["multiplication_table"].size(), N_op);
    for(auto i = 0; i < json["multiplication_table"].size(); ++i) {
      BOOST_CHECK_EQUAL(json["multiplication_table"][i].size(), N_op);
    }
  }

  /// \brief Check that 'casm init' runs without error and expected files are created
  void Proj::check_init() {

    test::rm_project(*this);

    test::make_project(*this);

    m_set = ProjectSettings(dir);

    BOOST_CHECK_EQUAL(true, fs::exists(dir));

    // prim and settings
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.prim()));
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.project_settings()));

    // symmetry
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.crystal_point_group()));
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.factor_group()));
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.lattice_point_group()));

    // composition axes
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.composition_axes(m_set.calctype(), m_set.ref())));
  }

  /// \brief Check that 'casm sym' runs without error
  void Proj::check_symmetry() {
    m_p.popen(cd_and() + "casm sym");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }

  /// \brief Check number of symmetry operations and classes found
  void Proj::_check_symmetry(int lat_pg_op, int lat_pg_class,
                             int xtal_pg_op, int xtal_pg_class,
                             int fg_op, int fg_class,
                             std::string lat_pg_name, std::string xtal_pg_name) {


    check_symgroup(jsonParser(m_dirs.lattice_point_group()), lat_pg_op, lat_pg_class);
    check_symgroup(jsonParser(m_dirs.crystal_point_group()), xtal_pg_op, xtal_pg_class);
    check_symgroup(jsonParser(m_dirs.factor_group()), fg_op, fg_class);

    m_p.popen(cd_and() + "casm sym");

    BOOST_CHECK_EQUAL(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Lattice point group is:\s+)" + lat_pg_name)), true);
    BOOST_CHECK_EQUAL(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Crystal point group is:\s+)" + xtal_pg_name)), true);
  }

  /// \brief Default checks '-d' runs without error
  void Proj::check_composition() {
    m_p.popen(cd_and() + "casm composition -d");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }

  /// \brief Default checks that enumerating supercells and configurations can
  /// be run for '--max 2' without error, but doesn't check results
  void Proj::check_enum() {
    m_p.popen(cd_and() + "casm enum --supercells --max 2");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
    m_p.popen(cd_and() + "casm enum --configs --max 2");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }

  /// \brief Default checks that several options run without error
  void Proj::check_select() {
    m_p.popen(cd_and() + "casm select --set-on");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);

    m_p.popen(cd_and() + "casm select --set-off");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }

  /// \brief Default checks that several options run without error
  void Proj::check_query() {
    m_p.popen(cd_and() + "casm query --columns comp");
    BOOST_CHECK_MESSAGE(m_p.exit_code() == 0, m_p.gets());

  }

}

