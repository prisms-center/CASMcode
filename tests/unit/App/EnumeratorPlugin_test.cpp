#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
///   the command line executable

/// What is being used to test it:
#include <boost/filesystem.hpp>

#include "Common.hh"
#include "casm/app/casm_functions.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(EnumeratorPluginTest)

BOOST_AUTO_TEST_CASE(Test1) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  PrimClex primclex(proj.dir, null_log());

  auto cp = [&](std::string _filename) {

    fs::path filename(_filename);
    fs::path src = "tests/unit/App" / filename;
    BOOST_REQUIRE(fs::exists(src));

    fs::path dest = primclex.dir().enumerator_plugins();
    fs::create_directories(dest);
    BOOST_REQUIRE(fs::exists(dest));

    fs::copy_file(src, dest / filename, fs::copy_option::overwrite_if_exists);
    BOOST_REQUIRE(fs::exists(dest / filename));

  };

  cp("TestEnum.hh");
  cp("TestEnum.cc");

  // refresh to load plugins
  primclex.refresh(true);

  auto check = [&](std::string str) {
    CommandArgs args(str, &primclex, primclex.dir().root_dir(), primclex);
    BOOST_CHECK(!enum_command(args));
  };

  check(R"(enum -h)");

  check(R"(enum --desc TestEnum)");

  check(R"(enum --method TestEnum -i '{"supercells": {"max": 4}}')");

  BOOST_CHECK_EQUAL(std::distance(primclex.config_begin(), primclex.config_end()), 336);

}

BOOST_AUTO_TEST_SUITE_END()
