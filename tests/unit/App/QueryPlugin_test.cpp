#include "gtest/gtest.h"

/// What is being tested:
///   the command line executable

/// What is being used to test it:
#include <boost/filesystem.hpp>

#include "Common.hh"
#include "casm/app/casm_functions.hh"

using namespace CASM;

TEST(QueryPlugin, Test1) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();
  proj.check_enum();

  PrimClex primclex(proj.dir, Logging(null_log(), null_log(), null_log()));

  //TODO: This is more code duplication
  auto cp = [&](std::string _filename) {

    fs::path filename(_filename);
    fs::path src = "tests/unit/App" / filename;
    ASSERT_TRUE(fs::exists(src));

    fs::path dest = primclex.dir().query_plugins<Configuration>();
    fs::create_directories(dest);
    ASSERT_TRUE(fs::exists(dest));

    fs::copy_file(src, dest / filename, fs::copy_option::overwrite_if_exists);
    ASSERT_TRUE(fs::exists(dest / filename));

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

  EXPECT_TRUE(check(R"(ccasm select -h)"));

  EXPECT_TRUE(check(R"(ccasm select --set-on)"));

  EXPECT_TRUE(check(R"(ccasm query -h p)"));

  EXPECT_TRUE(check(R"(ccasm query -k test_comp_n)"));

  EXPECT_TRUE(check(R"(ccasm query -k 'test_comp_n(Zr)')"));

  EXPECT_TRUE(check(R"(ccasm query -k 'test_comp_n(O)')"));

  EXPECT_TRUE(check(R"(ccasm query -k 'test_configname')"));

}

