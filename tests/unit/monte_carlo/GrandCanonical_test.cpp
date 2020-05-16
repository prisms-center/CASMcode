#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
///   the command line executable

/// What is being used to test it:
#include <boost/filesystem.hpp>

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "casm/app/casm_functions.hh"
#include "casm/monte_carlo/grand_canonical/GrandCanonical.hh"

using namespace CASM;


TEST(GrandCanonicalTest, Test0) {

  std::cout << "skipping GrandCanonical_test" << std::endl;
  if(true) {
    return;
  }

  test::ZrOProj proj;
  proj.check_init();
  proj.check_composition();

  Logging logging = Logging::null();
  //Logging logging;
  PrimClex primclex(proj.dir, logging);

  fs::path eci_src = autotools::abs_srcdir() + "/tests/unit/monte_carlo/eci_0.json";
  fs::path eci_dest = primclex.dir().eci("formation_energy", "default", "default", "default", "default");
  fs::copy_file(eci_src, eci_dest, fs::copy_option::overwrite_if_exists);

  fs::path bspecs_src = autotools::abs_srcdir() + "/tests/unit/monte_carlo/bspecs_0.json";
  fs::path bspecs_dest = primclex.dir().bspecs("default");
  fs::copy_file(bspecs_src, bspecs_dest, fs::copy_option::overwrite_if_exists);

  fs::path settings_src = autotools::abs_srcdir() + "/tests/unit/monte_carlo/metropolis_grand_canonical_0.json";
  fs::path mc_dir = primclex.dir().root_dir() / "mc_0";
  fs::create_directory(mc_dir);
  fs::path settings_dest = mc_dir / settings_src.filename();
  fs::copy_file(settings_src, settings_dest, fs::copy_option::overwrite_if_exists);


  // for autotools
  primclex.settings().set_casm_libdir(fs::current_path() / ".libs");
  commit(primclex.settings());

  auto check = [&](std::string str) {
    CommandArgs args(str, &primclex, primclex.dir().root_dir(), Logging::null());
    return !casm_api(args);
  };

  EXPECT_TRUE(check(R"(casm bset -u)"));

  EXPECT_TRUE(check(std::string("casm monte -s ") + settings_dest.string()));

}
