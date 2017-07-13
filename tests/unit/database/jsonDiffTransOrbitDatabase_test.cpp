#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/database/DiffTransOrbitDatabase.hh"
#include "casm/database/json/jsonDatabase.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
//#include "casm/app/AppIO.hh"
//#include "casm/app/AppIO_impl.hh"
#include "casm/crystallography/Site.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clex/PrimClex.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(jsonDiffTransOrbitDatabase_Test)

BOOST_AUTO_TEST_CASE(Test1) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  BOOST_CHECK_EQUAL(fs::equivalent(proj.dir, primclex.dir().root_dir()), true);

  DB::jsonDatabase<Kinetics::PrimPeriodicDiffTransOrbit> db_diff_trans(primclex);

  db_diff_trans.open();
  BOOST_CHECK_EQUAL(db_diff_trans.size(), 0);

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  fs::path difftrans_path = "tests/unit/kinetics/diff_trans.json";
  jsonParser diff_trans_json {difftrans_path};
  auto end = make_prim_periodic_orbits(
               primclex.prim(), diff_trans_json["bspecs"], alloy_sites_filter, primclex.crystallography_tol(), std::back_inserter(orbits), primclex.log());

  std::vector< PrimPeriodicDiffTransOrbit > diff_trans_orbits;
  auto end2 = make_prim_periodic_diff_trans_orbits(
                orbits.begin(), orbits.end(), primclex.crystallography_tol(), std::back_inserter(diff_trans_orbits));
  for(auto it = diff_trans_orbits.begin(); it != diff_trans_orbits.end(); ++it) {
    db_diff_trans.insert(*it);
  }
  BOOST_CHECK_EQUAL(db_diff_trans.size(), 28);

  db_diff_trans.commit();

  db_diff_trans.close();
  BOOST_CHECK_EQUAL(db_diff_trans.size(), 0);

  db_diff_trans.open();
  BOOST_CHECK_EQUAL(db_diff_trans.size(), 28);


}

BOOST_AUTO_TEST_SUITE_END()
