#include "gtest/gtest.h"

/// What is being tested:
#include "casm/database/ScelDatabase.hh"
#include "casm/database/json/jsonDatabase.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/Niggli.hh"

using namespace CASM;

TEST(ScelDatabase_Test, Test1) {

  EXPECT_EQ(1, 1);
  test::FCCTernaryProj proj;
  proj.check_init();
  EXPECT_EQ(1, 1);

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  EXPECT_EQ(fs::equivalent(proj.dir, primclex.dir().root_dir()), true);

  DB::jsonDatabase<Supercell> db_scel(primclex);
  EXPECT_EQ(1, 1);

  db_scel.open();
  EXPECT_EQ(1, 1);
  EXPECT_EQ(db_scel.size(), 0);

  Supercell scel(&primclex, Eigen::Matrix3i::Identity());
  db_scel.insert(scel);
  EXPECT_EQ(db_scel.size(), 1);

  db_scel.erase(scel.name());
  EXPECT_EQ(db_scel.size(), 0);

  db_scel.emplace(&primclex, Eigen::Matrix3i::Identity());
  EXPECT_EQ(db_scel.size(), 1);

  db_scel.insert(scel);
  EXPECT_EQ(db_scel.size(), 1);

  db_scel.emplace(&primclex, Eigen::Matrix3i::Identity());
  EXPECT_EQ(db_scel.size(), 1);

  EXPECT_EQ(db_scel.name(scel.name()), scel.name());
  EXPECT_EQ(db_scel.alias(scel.name()).empty(), true);

  std::string other("other_name");

  EXPECT_EQ(db_scel.name(other), other);
  EXPECT_EQ(db_scel.alias(other).empty(), true);

  auto it = db_scel.find(other);
  EXPECT_EQ(it == db_scel.end(), true);
  EXPECT_EQ(it != db_scel.end(), false);

  it = db_scel.find(scel.name());
  EXPECT_EQ(it != db_scel.end(), true);
  EXPECT_EQ(it == db_scel.end(), false);

  EXPECT_EQ(db_scel.count(scel.name()), 1);
  EXPECT_EQ(db_scel.count(other), 0);

  db_scel.erase(it);
  EXPECT_EQ(db_scel.size(), 0);

  int minvol = 1;
  int maxvol = 10;
  ScelEnumProps enum_props(minvol, maxvol + 1);
  auto fg = prim.factor_group();
  SuperlatticeEnumerator lat_enum(fg.begin(), fg.end(), prim.lattice(), enum_props);
  EXPECT_EQ(true, true);
  EXPECT_EQ(std::distance(lat_enum.begin(), lat_enum.end()), 87);
  for(auto it = lat_enum.begin(); it != lat_enum.end(); ++it) {
    db_scel.emplace(&primclex, xtal::canonical::equivalent(*it, prim.point_group()));
  }
  EXPECT_EQ(db_scel.size(), 87);

}
