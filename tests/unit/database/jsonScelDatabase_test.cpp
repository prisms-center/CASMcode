#include "gtest/gtest.h"

/// What is being tested:
#include "casm/database/ScelDatabase.hh"
#include "casm/database/json/jsonDatabase.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"

using namespace CASM;


TEST(jsonScelDatabase_Test, Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  EXPECT_EQ(fs::equivalent(proj.dir, primclex.dir().root_dir()), true);

  DB::jsonDatabase<Supercell> db_scel(primclex);

  db_scel.open();
  EXPECT_EQ(db_scel.size(), 0);

  int minvol = 1;
  int maxvol = 10;
  ScelEnumProps enum_props(minvol, maxvol + 1);
  SuperlatticeEnumerator lat_enum(prim.lattice(), Adapter::symop_to_matrix(prim.factor_group()), enum_props);
  for(auto it = lat_enum.begin(); it != lat_enum.end(); ++it) {
    db_scel.emplace(&primclex, it->canonical_form(prim.point_group()));
  }
  EXPECT_EQ(db_scel.size(), 87);

  db_scel.commit();
  //fs::ifstream file(primclex.dir().scel_list());
  //std::cout << file.rdbuf() << std::endl;

  db_scel.close();
  EXPECT_EQ(db_scel.size(), 0);

  db_scel.open();
  EXPECT_EQ(db_scel.size(), 87);

}

