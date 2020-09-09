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
#include "casm/database/Database.hh"

// template definitions
#include "casm/app/QueryHandler_impl.hh"
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/app/enum/EnumInterface_impl.hh"

using namespace CASM;

TEST(Selection_Test, Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  //const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);

  Completer::EnumOption enum_opt;
  enum_opt.desc();

  // -- Generate Supercell & Configuration --

  ScelEnumByProps enum_scel(primclex, xtal::ScelEnumProps(1, 5));
  EXPECT_EQ(true, true);

  ConfigEnumAllOccupations::run(primclex, enum_scel.begin(), enum_scel.end());
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
