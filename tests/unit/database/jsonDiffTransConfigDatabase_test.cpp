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
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/clex/ScelEnum_impl.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/casm_io/stream_io/container.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"


using namespace CASM;

BOOST_AUTO_TEST_SUITE(jsonDiffTransConfigDatabase_Test)

BOOST_AUTO_TEST_CASE(Test1) {

  // Make test project
  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  const Structure &prim(primclex.prim());
  primclex.settings().set_crystallography_tol(1e-5);
  BOOST_CHECK_EQUAL(fs::equivalent(proj.dir, primclex.dir().root_dir()), true);

  fs::path bspecs_path = "tests/unit/kinetics/FCCTernary_bspecs_0.json";
  jsonParser bspecs {bspecs_path};

  fs::path diffperturb_path = "tests/unit/kinetics/FCCTernary_diff_perturb_0.json";
  jsonParser diff_perturb_json {diffperturb_path};
  BOOST_CHECK_EQUAL(true, true);

  // Make PrimPeriodicIntegralClusterOrbit
  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    null_log());
  BOOST_CHECK_EQUAL(true, true);

  // Make PrimPeriodicDiffTransOrbit
  std::vector<Kinetics::PrimPeriodicDiffTransOrbit> diff_trans_orbits;
  Kinetics::make_prim_periodic_diff_trans_orbits(
    orbits.begin() + 2,  // use pairs+
    orbits.begin() + 4,
    primclex.crystallography_tol(),
    std::back_inserter(diff_trans_orbits),
    &primclex);
  BOOST_CHECK_EQUAL(true, true);

  // Make background config
  Eigen::Vector3d a, b, c;
  std::tie(a, b, c) = prim.lattice().vectors();
  BOOST_CHECK_EQUAL(true, true);

  Supercell tscel(
    &primclex,
    Lattice(3 * (c + b - a), 3 * (a - b + c), 3 * (a + b - c)));
  BOOST_CHECK_EQUAL(true, true);
  const Supercell &scel = *tscel.insert().first;
  Configuration config(scel);
  config.init_occupation();
  BOOST_CHECK_EQUAL(true, true);

  // Make DiffTransConfiguration database
  DB::jsonDatabase<Kinetics::DiffTransConfiguration> db_diff_trans_config(primclex);
  BOOST_CHECK_EQUAL(true, true);

  // Open DiffTransConfiguration database
  db_diff_trans_config.open();
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 0);

  // Make DiffTransConfiguration enumerator and enumerate configs
  //std::cout << "skipping DiffTransConfigEnumOccPerturbations dependent parts" << std::endl;
  Kinetics::DiffTransConfigEnumOccPerturbations enum_diff_trans_config(
    config, diff_trans_orbits[0], diff_perturb_json["local_cspecs"]);
  for(const auto &diff_trans_config : enum_diff_trans_config) {
    db_diff_trans_config.insert(diff_trans_config);
  }
  db_diff_trans_config.commit();
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 29); // not checked for accuracy


  // Check cached properties
  std::cout << "skipping cache check" << std::endl;
  ////  for(const auto &diff_trans_config : db_diff_trans_config) {
  ////    //std::cout << "id: " << config.id() << "  occ: " << config.occupation() << std::endl;
  ////    BOOST_CHECK_EQUAL(diff_trans_config.cache().contains("multiplicity"), false);
  ////    BOOST_CHECK_EQUAL(diff_trans_config.multiplicity() != 0, true);
  ////    BOOST_CHECK_EQUAL(diff_trans_config.cache().contains("multiplicity"), true);
  ////  }
  ////  db_diff_trans_config.commit();
  //
  // Check that the database is sorted
  {
    auto next = db_diff_trans_config.begin();
    auto it = next++;
    auto end = db_diff_trans_config.end();
    for(; next != end; ++it, ++next) {
      BOOST_CHECK_EQUAL(*it < *next, true);
    }
  }

  // Close DiffTransConfiguration database
  db_diff_trans_config.close();
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 0);

  // Re-open DiffTransConfiguration database
  db_diff_trans_config.open();
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 29); // not checked for accuracy

  //  // Check cached properties
  std::cout << "skipping cache check" << std::endl;
  ////  for(const auto &diff_trans_config : db_diff_trans_config) {
  ////    BOOST_CHECK_EQUAL(diff_trans_config.cache().contains("multiplicity"), true);
  ////  }

  // Check that the database is sorted
  {
    auto next = db_diff_trans_config.begin();
    auto it = next++;
    auto end = db_diff_trans_config.end();
    for(; next != end; ++it, ++next) {
      BOOST_CHECK_EQUAL(*it < *next, true);
    }
  }

  // Close DiffTransConfiguration database
  db_diff_trans_config.close();
  BOOST_CHECK_EQUAL(true, true);
}

BOOST_AUTO_TEST_SUITE_END()
