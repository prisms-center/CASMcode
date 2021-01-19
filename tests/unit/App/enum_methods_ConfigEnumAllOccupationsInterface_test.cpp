#include "App/TestEnumeratorInterface.hh"
#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/app/enum.hh"
#include "casm/app/enum/methods/ConfigEnumAllOccupationsInterface.hh"
#include "casm/app/enum/methods/ScelEnumInterface.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/completer/Handlers.hh"
#include "casm/database/Database.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "gtest/gtest.h"

TEST(enum_methods_ConfigEnumAllOccupationsInterfaceTest, Test1) {
  ScopedNullLogging logging;
  test::FCCTernaryProj proj;
  proj.check_init();
  PrimClex primclex(proj.dir);

  // --dry-run test
  {
    std::string cli_str = "casm enum --method ScelEnum --max 4 --dry-run";
    test::run_enum_interface<ScelEnumInterface>(cli_str, primclex);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  primclex.generic_db<Supercell>().close();
  primclex.generic_db<Supercell>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 0);

  {
    std::string cli_str = "casm enum --method ScelEnum --max 4";
    test::run_enum_interface<ScelEnumInterface>(cli_str, primclex);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  primclex.generic_db<Supercell>().close();
  primclex.generic_db<Supercell>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);

  // --dry-run test
  {
    std::string cli_str =
        "casm enum --method ConfigEnumAllOccupations -a --dry-run";
    test::run_enum_interface<ConfigEnumAllOccupationsInterface>(cli_str,
                                                                primclex);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 0);

  {
    std::string cli_str = "casm enum --method ConfigEnumAllOccupations -a";
    test::run_enum_interface<ConfigEnumAllOccupationsInterface>(cli_str,
                                                                primclex);
  }
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);
  primclex.generic_db<Configuration>().close();
  primclex.generic_db<Configuration>().open();
  EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
  EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);
}
