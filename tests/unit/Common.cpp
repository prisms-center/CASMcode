#include "Common.hh"
#include "gtest/gtest.h"
#include "autotools.hh"

#include <thread>
#include <chrono>
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include "Proj.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/AppIO.hh"

//TODO: I'm pretty sure having these macros outside of the testing blocks is
//causing undefined behavior
//

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
  void check(std::string test,
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

    EXPECT_EQ(ok, true);

    return;
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

    if(find_casmroot(dir) != dir) {
      // build a project
      build_project(make_default_project_settings(prim, title, dir), Structure {prim});
    }

    // (re)load ProjectSettings
    m_set = notstd::clone(open_project_settings(dir));
  }

  /// \brief Check some aspects of a SymGroup json
  void check_symgroup(const jsonParser &json, int N_op, int N_class) {
    // temporarily disabled while transiting character table determination
    //EXPECT_EQ(json["character_table"].size(), N_class);

    //EXPECT_EQ(json["conjugacy_class"].size(), N_class);

    EXPECT_EQ(json["group_operations"].size(), N_op);
    EXPECT_EQ((*json["group_operations"].begin())["info"]["type"].get<std::string>(), "identity");

    EXPECT_EQ(json["group_structure"]["multiplication_table"].size(), N_op);
    for(auto i = 0; i < json["group_structure"]["multiplication_table"].size(); ++i) {
      EXPECT_EQ(json["group_structure"]["multiplication_table"][i].size(), N_op);
    }
  }

  /// \brief Check that 'casm init' runs without error and expected files are created
  void Proj::check_init() {

    make();

    m_set->set_casm_libdir(autotools::abs_libdir());
    m_set->set_casm_includedir(autotools::abs_includedir());

    commit(*m_set);

    EXPECT_EQ(true, fs::exists(dir));

    // prim and settings
    EXPECT_EQ(true, fs::exists(m_dirs.prim()));
    EXPECT_EQ(true, fs::exists(m_dirs.project_settings()));

    // symmetry
    EXPECT_EQ(true, fs::exists(m_dirs.crystal_point_group()));
    EXPECT_EQ(true, fs::exists(m_dirs.factor_group()));
    EXPECT_EQ(true, fs::exists(m_dirs.lattice_point_group()));

    // composition axes
    EXPECT_EQ(true, fs::exists(m_dirs.composition_axes()));
  }

  /// \brief Check that 'casm sym' runs without error
  void Proj::check_symmetry() {
    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " sym");
    EXPECT_EQ(m_p.exit_code(), 0);
  }

  /// \brief Check number of symmetry operations and classes found
  void Proj::_check_symmetry(int lat_pg_op, int lat_pg_class,
                             int xtal_pg_op, int xtal_pg_class,
                             int fg_op, int fg_class,
                             std::string lat_pg_name, std::string xtal_pg_name) {


    check_symgroup(jsonParser(m_dirs.lattice_point_group()), lat_pg_op, lat_pg_class);
    check_symgroup(jsonParser(m_dirs.crystal_point_group()), xtal_pg_op, xtal_pg_class);
    check_symgroup(jsonParser(m_dirs.factor_group()), fg_op, fg_class);

    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " sym");

    EXPECT_EQ(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Lattice point group is:\s+)" + lat_pg_name)), true);
    EXPECT_EQ(boost::regex_search(m_p.gets(), m_match, boost::regex(R"(Crystal point group is:\s+)" + xtal_pg_name)), true);
  }

  /// \brief Default checks '-d' runs without error
  void Proj::check_composition() {
    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " composition -d");
    EXPECT_EQ(m_p.exit_code(), 0);
  }


  /// \brief Default checks that enumerating supercells and configurations can
  /// be run for '--max 2' without error, but doesn't check results
  void Proj::check_enum() {
    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " enum --method ScelEnum --max 2");
    EXPECT_EQ(m_p.exit_code(), 0);

    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " enum --method ConfigEnumAllOccupations --all");
    EXPECT_EQ(m_p.exit_code(), 0);
  }

  /// \brief Default checks that several options run without error
  void Proj::check_select() {
    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " select --set-on");
    EXPECT_EQ(m_p.exit_code(), 0);

    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " select --set-off");
    EXPECT_EQ(m_p.exit_code(), 0);
  }

  /// \brief Default checks that several options run without error
  void Proj::check_query() {
    m_p.popen(cd_and() + autotools::abs_ccasm_path() + " query --columns comp");
    EXPECT_EQ(m_p.exit_code(), 0) << m_p.gets();
  }

}
