#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/ConfigEnumAllOccupations.hh"

/// What is being used to test it:
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/app/AppIO.hh"
#include "casm/app/ProjectBuilder.hh"
#include "Common.hh"

using namespace CASM;
using namespace test;

BOOST_AUTO_TEST_SUITE(ConfigEnumTest)

BOOST_AUTO_TEST_CASE(ConfigEnumAllOccupationsTest) {

  // tests ConfigEnumAllOccupations and ConfigEnumEquivalents

  // read test file
  fs::path test_cases_path("tests/unit/clex/ConfigEnumAllOccupations_test_cases.json");
  jsonParser tests(test_cases_path);
  // double tol = TOL;
  fs::path test_proj_dir("tests/unit/clex/test_proj");

  for(auto test_it = tests.begin(); test_it != tests.end(); ++test_it) {

    // input and expected output data
    jsonParser &j = *test_it;

    // if false: print calculated results if no test data; if true: suppress
    bool quiet = false;
    j.get_else(quiet, "quiet", false);

    BOOST_CHECK_MESSAGE(j.contains("title"), "test case 'title' is required");
    BOOST_CHECK_MESSAGE(j.contains("prim"), "test case 'prim' is required");
    BOOST_CHECK_MESSAGE(j.contains("min_vol"), "test case 'min_vol' is required");
    BOOST_CHECK_MESSAGE(j.contains("max_vol"), "test case 'max_vol' is required");

    // generate prim
    Structure prim(read_prim(j["prim"]));

    // clean up test proj
    if(fs::exists(test_proj_dir / ".casm")) {
      fs::remove_all(test_proj_dir);
    }

    fs::create_directory(test_proj_dir);

    j["prim"].write(test_proj_dir / "prim.json");

    // build a project
    ProjectBuilder builder(test_proj_dir, j["title"].get<std::string>(), "formation_energy");
    builder.build();

    // read primclex
    PrimClex primclex(test_proj_dir, null_log());

    // generate supercells
    ScelEnumProps enum_props(j["min_vol"].get<int>(), j["max_vol"].get<int>() + 1);
    primclex.generate_supercells(enum_props);

    // run checks:
    jsonParser json_scel;
    json_scel = (Index) primclex.get_supercell_list().size();
    check("Nscel", j, json_scel, test_cases_path, quiet);

    // generate configurations
    jsonParser json = jsonParser::array();
    for(auto &scel : primclex.get_supercell_list()) {
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

BOOST_AUTO_TEST_SUITE_END()
