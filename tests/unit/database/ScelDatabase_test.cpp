#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/database/ScelDatabase.hh"
#include "casm/database/json/jsonDatabase.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SupercellEnumerator.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(jsonDatabase_Supercell_Test)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);

  BOOST_CHECK_EQUAL(proj.dir, primclex.dir().root_dir());

  DB::jsonDatabase<Supercell> db_scel(primclex);
  BOOST_CHECK_EQUAL(1, 1);

  db_scel.open();
  BOOST_CHECK_EQUAL(1, 1);
  BOOST_CHECK_EQUAL(db_scel.size(), 0);

  Supercell scel(&primclex, Eigen::Matrix3i::Identity());
  db_scel.insert(scel);
  BOOST_CHECK_EQUAL(db_scel.size(), 1);

  db_scel.erase(scel.name());
  BOOST_CHECK_EQUAL(db_scel.size(), 0);

  db_scel.emplace(&primclex, Eigen::Matrix3i::Identity());
  BOOST_CHECK_EQUAL(db_scel.size(), 1);

  db_scel.insert(scel);
  BOOST_CHECK_EQUAL(db_scel.size(), 1);

  db_scel.emplace(&primclex, Eigen::Matrix3i::Identity());
  BOOST_CHECK_EQUAL(db_scel.size(), 1);

  BOOST_CHECK_EQUAL(db_scel.name(scel.name()), scel.name());
  BOOST_CHECK_EQUAL(db_scel.alias(scel.name()).empty(), true);

  std::string other("other_name");

  BOOST_CHECK_EQUAL(db_scel.name(other), other);
  BOOST_CHECK_EQUAL(db_scel.alias(other).empty(), true);

  auto it = db_scel.find(other);
  BOOST_CHECK_EQUAL(it == db_scel.end(), true);
  BOOST_CHECK_EQUAL(it != db_scel.end(), false);

  it = db_scel.find(scel.name());
  BOOST_CHECK_EQUAL(it != db_scel.end(), true);
  BOOST_CHECK_EQUAL(it == db_scel.end(), false);

  BOOST_CHECK_EQUAL(db_scel.count(scel.name()), 1);
  BOOST_CHECK_EQUAL(db_scel.count(other), 0);

  db_scel.erase(it);
  BOOST_CHECK_EQUAL(db_scel.size(), 0);

  int minvol = 1;
  int maxvol = 10;
  ScelEnumProps enum_props(minvol, maxvol + 1);
  SupercellEnumerator<Lattice> lat_enum(prim.lattice(), prim.factor_group(), enum_props);
  for(auto it = lat_enum.begin(); it != lat_enum.end(); ++it) {
    db_scel.emplace(&primclex, it.matrix());
  }
  BOOST_CHECK_EQUAL(db_scel.size(), 87);

}

BOOST_AUTO_TEST_SUITE_END()
