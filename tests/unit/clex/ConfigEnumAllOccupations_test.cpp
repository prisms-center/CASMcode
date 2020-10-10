#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/clex/ConfigEnumAllOccupations.hh"

/// What is being used to test it:
#include "casm/crystallography/Structure.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/database/Database.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/enum.hh"
#include "casm/app/enum/methods/ConfigEnumAllOccupationsInterface.hh"
#include "casm/app/enum/methods/ScelEnumInterface.hh"
#include "casm/app/CLIParse.hh"

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "clex/TestClexEnumerators.hh"

using namespace CASM;
using namespace test;

TEST(ConfigEnumTest, ConfigEnumAllOccupationsTest) {

  // tests ConfigEnumAllOccupations and ConfigEnumEquivalents
  ScopedNullLogging logging;

  // read test file
  fs::path test_cases_path(autotools::abs_srcdir() + "/tests/unit/clex/ConfigEnumAllOccupations_test_cases.json");
  jsonParser tests(test_cases_path);
  double tol = TOL;
  fs::path test_proj_dir(autotools::abs_srcdir() + "/tests/unit/clex/test_proj");

  for(auto test_it = tests.begin(); test_it != tests.end(); ++test_it) {

    // input and expected output data
    jsonParser &j = *test_it;

    // if false: print calculated results if no test data; if true: suppress
    bool quiet = false;
    j.get_else(quiet, "quiet", false);

    EXPECT_TRUE(j.contains("title")) << "test case 'title' is required";
    EXPECT_TRUE(j.contains("prim")) << "test case 'prim' is required";
    EXPECT_TRUE(j.contains("min_vol")) << "test case 'min_vol' is required";
    EXPECT_TRUE(j.contains("max_vol")) << "test case 'max_vol' is required";

    // generate prim
    Structure prim(read_prim(j["prim"], tol));

    // clean up test proj
    if(fs::exists(test_proj_dir / ".casm")) {
      fs::remove_all(test_proj_dir);
    }

    fs::create_directory(test_proj_dir);

    // build a project
    auto project_name = j["title"].get<std::string>();
    auto project_settings = make_default_project_settings(prim, project_name, test_proj_dir);
    build_project(project_settings, prim);

    // read primclex
    PrimClex primclex(test_proj_dir);

    // generate supercells
    xtal::ScelEnumProps enum_props(j["min_vol"].get<int>(), j["max_vol"].get<int>() + 1);
    ScelEnumByProps scel_enum(primclex.shared_prim(), enum_props);
    for(const auto &scel : scel_enum) {
      (void) scel;
    }

    // run checks:
    jsonParser json_scel;
    json_scel = (Index) primclex.generic_db<Supercell>().size();
    check("Nscel", j, json_scel, test_cases_path, quiet);

    // generate configurations
    jsonParser json = jsonParser::array();
    for(auto &scel : primclex.generic_db<Supercell>()) {
      ConfigEnumAllOccupations e(scel);
      json.push_back(std::distance(e.begin(), e.end()));
    }

    // run checks:
    check("Nconfigs", j, json, test_cases_path, quiet);

    // ... add more here ...

    // clean up test proj
    if(fs::exists(test_proj_dir / ".casm")) {
      fs::remove_all(test_proj_dir);
    }
  }
}

TEST(ConfigEnumTest, ConfigEnumAllOccupationsInterfaceTest) {

  // create a project
  test::FCCTernaryProj proj;
  proj.check_init();

  // construct PrimClex
  ScopedNullLogging logging;
  PrimClex primclex(proj.dir);

  // --dry-run test
  {
    std::string cli_str = "casm enum --method ScelEnum --max 4 --dry-run";
    run_enum_interface<ScelEnumInterface>(cli_str, primclex);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  primclex.generic_db<Supercell>().close();
  primclex.generic_db<Supercell>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 0);


  {
    std::string cli_str = "casm enum --method ScelEnum --max 4";
    run_enum_interface<ScelEnumInterface>(cli_str, primclex);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  primclex.generic_db<Supercell>().close();
  primclex.generic_db<Supercell>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);

  // --dry-run test
  {
    std::string cli_str = "casm enum --method ConfigEnumAllOccupations -a --dry-run";
    run_enum_interface<ConfigEnumAllOccupationsInterface>(cli_str, primclex);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 0);

  {
    std::string cli_str = "casm enum --method ConfigEnumAllOccupations -a";
    run_enum_interface<ConfigEnumAllOccupationsInterface>(cli_str, primclex);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);

}
