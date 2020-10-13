#include "gtest/gtest.h"

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "App/TestEnumeratorInterface.hh"

#include "casm/app/enum.hh"
#include "casm/app/enum/methods/ConfigEnumRandomOccupationsInterface.hh"
#include "casm/app/enum/methods/ScelEnumInterface.hh"
#include "casm/clex/ConfigEnumRandomOccupations.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/Database.hh"
#include "casm/enumerator/ConfigEnumInput.hh"

TEST(ConfigEnumRandomOccupationsTest, ConfigEnumRandomOccupationsRunTest) {

  ScopedNullLogging logging;
  test::FCCTernaryProj proj;
  proj.check_init();
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
