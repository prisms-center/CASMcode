/// What is being tested:
///   the command line executable

#include "gtest/gtest.h"

/// What is being used to test it:
#include <boost/filesystem.hpp>

#include "Common.hh"
#include "ZrOProj.hh"
#include "casm/app/APICommand.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/enum.hh"

using namespace CASM;

TEST(EnumeratorPlugin, Test1) {
  // ScopedNullLogging logging;
  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  PrimClex primclex(proj.dir);

  auto cp = [&](std::string _filename) {
    fs::path src = test::data_file("App", _filename);
    ASSERT_TRUE(fs::exists(src));

    fs::path dest = primclex.dir().enumerator_plugins();
    fs::create_directories(dest);
    ASSERT_TRUE(fs::exists(dest));

    fs::copy_file(src, dest / _filename, fs::copy_option::overwrite_if_exists);
    ASSERT_TRUE(fs::exists(dest / _filename));
  };

  cp("TestEnum.cc");

  // refresh to load plugins
  primclex.settings().set_cxxflags(
      "-O3 -Wall -fPIC --std=c++17 -DGZSTREAM_NAMESPACE=gz");
  primclex.settings().set_soflags(
      "-shared -lboost_system -lboost_filesystem -lz");
  commit(primclex.settings());
  primclex.refresh(true);

  auto check = [&](std::string str) {
    CommandArgs args(str, &primclex, primclex.dir().root_dir());
    ASSERT_TRUE(!run_api_command<EnumCommand>(args));
  };

  check(R"(enum -h)");

  check(R"(enum --desc TestEnum)");

  check(R"(enum --method TestEnum -i '{"supercells": {"max": 4}}')");

  ASSERT_EQ(primclex.generic_db<Configuration>().size(), 336);
}
