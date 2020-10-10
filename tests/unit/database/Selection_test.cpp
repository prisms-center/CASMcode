#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/database/Selection.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/enumerator/Enumerator.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/app/enum.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/ConfigDatabaseTools_impl.hh"
#include "casm/database/Database.hh"
#include "casm/database/ScelDatabase.hh"

// template definitions
#include "casm/app/QueryHandler_impl.hh"
#include "casm/app/enum/EnumInterface_impl.hh"

using namespace CASM;

TEST(Selection_Test, Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  ScopedNullLogging logging;
  PrimClex primclex(proj.dir);
  //const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  auto &supercell_db = primclex.db<Supercell>();
  auto &configuration_db = primclex.db<Configuration>();

  Completer::EnumOption enum_opt;
  enum_opt.desc();

  // -- Generate Supercell & Configuration --

  ScelEnumByProps enum_scel {primclex.shared_prim(), xtal::ScelEnumProps(1, 5)};
  EXPECT_EQ(true, true);

  bool primitive_only = true;
  for(auto const &scel : enum_scel) {
    ConfigEnumAllOccupations enum_config {scel};
    for(auto const &config : enum_config) {
      make_canonical_and_insert(enum_config, config, supercell_db, configuration_db, primitive_only);
    }
  }
  EXPECT_EQ(true, true);


  // Test Supercell Selection
  {
    //auto &dict = primclex.settings().query_handler<Supercell>().dict();
    EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
    DB::Selection<Supercell> selection(primclex);
    EXPECT_EQ(selection.size(), 13);
    EXPECT_EQ(selection.selected_size(), 0);
  }

  // Test Configuration Selection
  {
    EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);

    DB::Selection<Configuration> selection(primclex);
    EXPECT_EQ(selection.size(), 126);
    EXPECT_EQ(selection.selected_size(), 0);

    auto &dict = primclex.settings().query_handler<Configuration>().dict();
    selection.set(dict, "lt(scel_size,3)");
    EXPECT_EQ(selection.selected_size(), 9);
  }
}
