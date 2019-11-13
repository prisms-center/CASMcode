/// What is being tested:
///   the command line executable

#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being used to test it:
#include <boost/filesystem.hpp>

#include "Common.hh"
#include "ZrOProj.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/APICommand.hh"
#include "casm/app/enum.hh"


using namespace CASM;

TEST(EnumeratorPlugin, Test1) {

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  PrimClex primclex(proj.dir, null_log());

  auto cp = [&](std::string _filename) {

    fs::path filename(_filename);
    fs::path src = std::string(autotools::abs_srcdir() + "/tests/unit/App") / filename;
    ASSERT_TRUE(fs::exists(src));

    fs::path dest = primclex.dir().enumerator_plugins();
    fs::create_directories(dest);
    ASSERT_TRUE(fs::exists(dest));

    fs::copy_file(src, dest / filename, fs::copy_option::overwrite_if_exists);
    ASSERT_TRUE(fs::exists(dest / filename));

  };

  cp("TestEnum.hh");
  cp("TestEnum_impl.hh");
  cp("TestEnum.cc");

  // refresh to load plugins
  primclex.refresh(true);

  auto check = [&](std::string str) {
    CommandArgs args(str, &primclex, primclex.dir().root_dir(), primclex);
    ASSERT_TRUE(!run_api_command<EnumCommand>(args));
  };

  check(R"(enum -h)");

  check(R"(enum --desc TestEnum)");

  check(R"(enum --method TestEnum -i '{"supercells": {"max": 4, "existing_only" : false}}')");

  ASSERT_EQ(primclex.generic_db<Configuration>().size(), 336);
}
