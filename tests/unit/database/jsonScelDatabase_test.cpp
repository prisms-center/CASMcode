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

BOOST_AUTO_TEST_SUITE(jsonScelDatabase_Test)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  BOOST_CHECK_EQUAL(fs::equivalent(proj.dir, primclex.dir().root_dir()), true);

  DB::jsonDatabase<Supercell> db_scel(primclex);

  db_scel.open();
  BOOST_CHECK_EQUAL(db_scel.size(), 0);

  int minvol = 1;
  int maxvol = 10;
  ScelEnumProps enum_props(minvol, maxvol + 1);
  SupercellEnumerator<Lattice> lat_enum(prim.lattice(), prim.factor_group(), enum_props);
  for(auto it = lat_enum.begin(); it != lat_enum.end(); ++it) {
    Supercell scel(&primclex, it.matrix());
    db_scel.insert(scel.canonical_form());
  }
  BOOST_CHECK_EQUAL(db_scel.size(), 87);

  db_scel.commit();
  //fs::ifstream file(primclex.dir().scel_list());
  //std::cout << file.rdbuf() << std::endl;

  db_scel.close();
  BOOST_CHECK_EQUAL(db_scel.size(), 0);

  db_scel.open();
  BOOST_CHECK_EQUAL(db_scel.size(), 87);

}

BOOST_AUTO_TEST_SUITE_END()
