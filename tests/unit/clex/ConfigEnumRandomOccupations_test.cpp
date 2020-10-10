#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/ConfigEnumRandomOccupations.hh"

/// What is being used to test it:
#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "clex/TestClexEnumerators.hh"
#include "casm/app/enum.hh"
#include "casm/app/enum/methods/ConfigEnumRandomOccupationsInterface.hh"
#include "casm/app/enum/methods/ScelEnumInterface.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/Database.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/external/MersenneTwister/MersenneTwister.h"

using namespace CASM;

TEST(ConfigEnumRandomOccupationsTest, Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  ScopedNullLogging logging;
  PrimClex primclex(proj.dir);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = primclex.prim().lattice().vectors();

  Supercell scel = Supercell(&primclex, Lattice {2.*a, 2.*b, 2.*c}).canonical_form();

  Index n_configs = 100;
  MTRand mtrand;

  {
    ConfigEnumRandomOccupations e(scel, n_configs, mtrand);
    EXPECT_EQ(n_configs, std::distance(e.begin(), e.end()));
  }


  {
    xtal::ScelEnumProps enum_props(1, 5);
    ScelEnumByProps e(primclex.shared_prim(), enum_props);
    for(const auto &scel : e) {
      (void) scel;
    }
    EXPECT_EQ(primclex.generic_db<Supercell>().size(), 14);
  }

  {
    jsonParser json_options;
    json_options["n_config"] = 200;
    std::string cli_str = "casm enum --method ConfigEnumRandomOccupations -a";
    test::run_enum_interface<ConfigEnumRandomOccupationsInterface>(cli_str, primclex, json_options);
    EXPECT_GT(primclex.generic_db<Configuration>().size(), 150);
  }

}

TEST(ConfigEnumRandomOccupationsTest, ConfigEnumRandomOccupationsRunTest) {

  // create a project
  test::FCCTernaryProj proj;
  proj.check_init();

  // construct PrimClex
  ScopedNullLogging logging;
  PrimClex primclex(proj.dir);

  {
    std::string cli_str = "casm enum --method ScelEnum --max 4";
    test::run_enum_interface<ScelEnumInterface>(cli_str, primclex);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  primclex.generic_db<Supercell>().close();
  primclex.generic_db<Supercell>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);

  std::string cmd = R"(casm enum --method ConfigEnumRandomOccupations -a)";
  jsonParser json_options;
  json_options["n_config"] = 10;

  // --dry-run test
  {
    std::string cli_str = cmd + " --dry-run";
    test::run_enum_interface<ConfigEnumRandomOccupationsInterface>(cli_str, primclex, json_options);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_GT(primclex.generic_db<Configuration>().size(), 0);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 0);

  {
    std::string cli_str = cmd;
    test::run_enum_interface<ConfigEnumRandomOccupationsInterface>(cli_str, primclex, json_options);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_GT(primclex.generic_db<Configuration>().size(), 0);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_GT(primclex.generic_db<Configuration>().size(), 0);

}
