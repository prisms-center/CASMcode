#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/database/DiffTransConfigDatabase.hh"
#include "casm/database/json/jsonDatabase.hh"


/// What is being used to test it:

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/casm_io/stream_io/container.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffTransConfigEnumPerturbations.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"


using namespace CASM;

BOOST_AUTO_TEST_SUITE(DiffTransConfigDatabase_Test)

BOOST_AUTO_TEST_CASE(Test1) {
  /// THIS TEST IS JUST A TEMPLATE NEEDS TO BE FULLY FILLED OUT
  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  BOOST_CHECK_EQUAL(fs::equivalent(proj.dir, primclex.dir().root_dir()), true);

  DB::jsonDatabase<Kinetics::DiffTransConfiguration> db_diff_trans_config(primclex);

  fs::path bspecs_path = "tests/unit/kinetics/bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  fs::path diffperturb_path = "tests/unit/kinetics/diff_perturb.json";
  jsonParser diff_perturb_json {diffperturb_path};

  // Make PrimPeriodicIntegralClusterOrbit
  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());
  // Make PrimPeriodicDiffTransOrbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 2,
    orbits.begin() + 4,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits));


  db_diff_trans_config.open();
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 0);

  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = prim.lattice().vectors();

  Supercell tscel(
    &primclex,
    Lattice(2.*a, 2.*b, c));
  const Supercell &scel = *tscel.insert().first;
  //construct an arbitrary DiffTransConfiguration here
  Configuration config(scel, jsonParser(), ConfigDoF({0, 0, 0, 0}));
  /*
  DiffusionTransformation diff_trans;
  Kinetics::DiffTransConfiguration diff_trans_config(config, diff_trans);
  auto res = db_diff_trans_config.insert(diff_trans_config);
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 1);
  BOOST_CHECK_EQUAL(db_diff_trans_config.begin()->id(), "0");

  db_diff_trans_config.erase(res.first);
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 0);
  */
  Kinetics::DiffTransConfigEnumPerturbations enum_diff_trans_config(config, diff_trans_orbits[0], diff_perturb_json["local_bspecs"]);
  for(const auto &diff_trans_config : enum_diff_trans_config) {
    db_diff_trans_config.insert(diff_trans_config);
  }
  db_diff_trans_config.commit();
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 56);
  /*
  for(const auto &diff_trans_config : db_diff_trans_config) {
    //std::cout << "id: " << config.id() << "  occ: " << config.occupation() << std::endl;
    BOOST_CHECK_EQUAL(diff_trans_config.cache().contains("multiplicity"), false);
    BOOST_CHECK_EQUAL(diff_trans_config.multiplicity() != 0, true);
    BOOST_CHECK_EQUAL(diff_trans_config.cache().contains("multiplicity"), true);
  }
  db_diff_trans_config.commit();
  */
  {
    auto next = db_diff_trans_config.begin();
    auto it = next++;
    auto end = db_diff_trans_config.end();
    for(; next != end; ++it, ++next) {
      BOOST_CHECK_EQUAL(*it < *next, true);
    }
  }

  db_diff_trans_config.close();
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 0);

  db_diff_trans_config.open();
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 56);
  /*for(const auto &diff_trans_config : db_diff_trans_config) {
    //std::cout << "id: " << config.id() << "  occ: " << config.occupation() << std::endl;
    BOOST_CHECK_EQUAL(diff_trans_config.cache().contains("multiplicity"), true);
  }*/
  {
    auto next = db_diff_trans_config.begin();
    auto it = next++;
    auto end = db_diff_trans_config.end();
    for(; next != end; ++it, ++next) {
      BOOST_CHECK_EQUAL(*it < *next, true);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
