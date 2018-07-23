#include "Common.hh"

#include <thread>
#include <chrono>
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/app/ProjectBuilder.hh"

namespace test {

  /// \brief Check expected JSON vs calculated JSON using BOOST_CHECK_EQUAL
  ///
  /// Checks:
  /// \code
  /// if(expected.contains(test)) {
  ///   BOOST_CHECK(expected[test].almost_equal, calculated);
  /// }
  /// \endcode
  ///
  /// If \code !expected.contains(test) && !quiet \endcode, print the calculated
  /// JSON so that it can be added to the test data.
  bool check(std::string test,
             const jsonParser &expected,
             const jsonParser &calculated,
             fs::path test_cases_path,
             bool quiet,
             double tol) {
    if(!expected.contains(test) && !quiet) {
      std::cout << "Test case: " << expected["title"] << " has no \"" << test << "\" test data." << std::endl;
      std::cout << "To use the current CASM results, add the following to the " << expected["title"]
                << " test case in " << test_cases_path << std::endl;
      jsonParser j = jsonParser::object();
      j[test] = calculated;
      std::cout << j << std::endl;
    }

    bool ok;
    fs::path diff_path;
    if(tol == 0.0) {
      ok = (expected[test] == calculated);
      if(!ok) {
        diff_path = find_diff(expected[test], calculated);
      }
    }
    else {
      ok = expected[test].almost_equal(calculated, tol);
      if(!ok) {
        fs::path diff_path = find_diff(expected[test], calculated, tol);
      }
    }

    if(!ok) {
      std::cout << "Difference at: " << diff_path << std::endl;
      std::cout << "Expected: \n" << expected[test].at(diff_path) << "\n"
                << "Found: \n" << calculated.at(diff_path) << std::endl;
    }
    BOOST_CHECK(ok);

    return ok;
  }

  /// \brief Create a new project directory, appending ".(#)" to ensure
  /// it is a new project
  fs::path proj_dir(fs::path init) {

    fs::path result = fs::absolute(init);
    int index = 0;
    std::string dot = ".";
    while(!fs::create_directories(result)) {
      result = fs::path(init.string() + dot + std::to_string(index));
      ++index;
    }

    return result;
  }

  /// \brief Build a CASM project at 'proj_dir/title' using the prim
  void Proj::make() {

    if(!fs::exists(dir / ".casm")) {
      jsonParser json;
      write_prim(prim, json, FRAC);
      json["description"] = desc;

      json.write(dir / "prim.json");

      // build a project
      ProjectBuilder builder(dir, title, "formation_energy");
      builder.build();
    }

    // (re)load ProjectSettings
    m_set = notstd::make_cloneable<ProjectSettings>(dir);
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

    make();

    m_set->set_casm_prefix(fs::current_path());

    // handle scons and autotools
    if(!fs::exists(m_set->casm_libdir().first / "libcasm.dylib") &&
       !fs::exists(m_set->casm_libdir().first / "libcasm.so")) {
      m_set->set_casm_libdir(fs::current_path() / ".libs");
    }

    m_set->commit();

    BOOST_CHECK_EQUAL(true, fs::exists(dir));

    // prim and settings
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.prim()));
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.project_settings()));

    // symmetry
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.crystal_point_group()));
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.factor_group()));
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.lattice_point_group()));

    // composition axes
    BOOST_CHECK_EQUAL(true, fs::exists(m_dirs.composition_axes()));
  }

  /// \brief Check that 'casm sym' runs without error
  void Proj::check_symmetry() {
    m_p.popen(cd_and() + "ccasm sym");
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

    m_p.popen(cd_and() + "ccasm sym");

    BOOST_CHECK_EQUAL(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Lattice point group is:\s+)" + lat_pg_name)), true);
    BOOST_CHECK_EQUAL(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Crystal point group is:\s+)" + xtal_pg_name)), true);
  }

  /// \brief Default checks '-d' runs without error
  void Proj::check_composition() {
    m_p.popen(cd_and() + "ccasm composition -d");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }

  /// \brief Default checks that enumerating supercells and configurations can
  /// be run for '--max 2' without error, but doesn't check results
  void Proj::check_enum() {
    m_p.popen(cd_and() + "ccasm enum --method ScelEnum --max 2");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);

    m_p.popen(cd_and() + "ccasm enum --method ConfigEnumAllOccupations --all");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }

  /// \brief Default checks that several options run without error
  void Proj::check_select() {
    m_p.popen(cd_and() + "ccasm select --set-on");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);

    m_p.popen(cd_and() + "ccasm select --set-off");
    BOOST_CHECK_EQUAL(m_p.exit_code(), 0);
  }

  /// \brief Default checks that several options run without error
  void Proj::check_query() {
    m_p.popen(cd_and() + "ccasm query --columns comp");
    BOOST_CHECK_MESSAGE(m_p.exit_code() == 0, m_p.gets());

  }

}
