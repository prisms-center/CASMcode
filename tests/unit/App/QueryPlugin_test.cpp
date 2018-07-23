#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
///   the command line executable

/// What is being used to test it:
#include <boost/filesystem.hpp>

#include "Common.hh"
#include "casm/app/casm_functions.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(QueryPluginTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();
  proj.check_enum();

  PrimClex primclex(proj.dir, Logging(null_log(), null_log(), null_log()));

  auto cp = [&](std::string _filename) {

    fs::path filename(_filename);
    fs::path src = "tests/unit/App" / filename;
    BOOST_REQUIRE(fs::exists(src));

    fs::path dest = primclex.dir().query_plugins<Configuration>();
    fs::create_directories(dest);
    BOOST_REQUIRE(fs::exists(dest));

    fs::copy_file(src, dest / filename, fs::copy_option::overwrite_if_exists);
    BOOST_REQUIRE(fs::exists(dest / filename));

  };

  // functor formatter
  cp("TestCompN.hh");
  cp("TestCompN.cc");

  // generic formatter
  cp("TestConfigname.hh");
  cp("TestConfigname.cc");

  // refresh to load plugins
  primclex.refresh(true);

  auto check = [&](std::string str) {
    CommandArgs args(str, &primclex, primclex.dir().root_dir(), Logging::null());
    return !casm_api(args);
  };

  BOOST_CHECK(check(R"(ccasm select -h)"));

  BOOST_CHECK(check(R"(ccasm select --set-on)"));

  BOOST_CHECK(check(R"(ccasm query -h p)"));

  BOOST_CHECK(check(R"(ccasm query -k test_comp_n)"));

  BOOST_CHECK(check(R"(ccasm query -k 'test_comp_n(Zr)')"));

  BOOST_CHECK(check(R"(ccasm query -k 'test_comp_n(O)')"));

  BOOST_CHECK(check(R"(ccasm query -k 'test_configname')"));

}

BOOST_AUTO_TEST_SUITE_END()
