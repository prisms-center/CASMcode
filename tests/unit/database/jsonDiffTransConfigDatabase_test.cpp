#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/database/DiffTransConfigDatabase.hh"
#include "casm/database/json/jsonDatabase.hh"


/// What is being used to test it:

#include "Common.hh"
#include "TestConfiguration.hh"
#include "TestOrbits.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SupercellEnumerator.hh"
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/clex/ScelEnum_impl.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/casm_io/jsonFile.hh"
#include "casm/casm_io/stream_io/container.hh"
#include "casm/kinetics/DiffusionTransformationEnum.hh"
#include "casm/kinetics/DiffTransConfigEnumOccPerturbations.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/app/AppIO_impl.hh"
#include "Common.hh"


using namespace CASM;

namespace {
  struct TestConfig0 : test::TestConfiguration {

    TestConfig0(const PrimClex &primclex) :
      TestConfiguration(
        primclex,
        lattice(primclex.prim()),
        std::vector<int>(108, 0)) {}

    static Lattice lattice(const Structure &prim) {
      Eigen::Vector3d a, b, c;
      std::tie(a, b, c) = prim.lattice().vectors();
      return Lattice(3 * (c + b - a), 3 * (a - b + c), 3 * (a + b - c));
    }

  };

  struct TestOrbits0 : test::TestPrimPeriodicDiffusionTransformationOrbits {
    TestOrbits0(const PrimClex &primclex) :
      test::TestPrimPeriodicDiffusionTransformationOrbits(
        primclex,
        jsonFile("tests/unit/kinetics/FCCTernary_bspecs_0.json"),
        2, 4) {}
  };

  struct TestEnumerator0 {

    TestEnumerator0() :
      diff_perturb_specs("tests/unit/kinetics/FCCTernary_diff_perturb_0.json") {}

    typedef Kinetics::DiffTransConfigEnumOccPerturbations EnumType;
    typedef std::unique_ptr<EnumType> EnumPtr;

    EnumPtr enumerator(
      const Configuration &bg_config,
      const Kinetics::PrimPeriodicDiffTransOrbit &dt_orbit) const {
      return notstd::make_unique<EnumType>(bg_config, dt_orbit, diff_perturb_specs["local_cspecs"]);
    }

    jsonFile diff_perturb_specs;
  };
}

BOOST_AUTO_TEST_SUITE(jsonDiffTransConfigDatabase_Test)

BOOST_AUTO_TEST_CASE(Test1) {

  // Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  PrimClex primclex(proj.dir, null_log());

  TestOrbits0 to(primclex);
  TestEnumerator0 te;
  TestConfig0 tc(primclex);

  // Make DiffTransConfiguration database
  DB::jsonDatabase<Kinetics::DiffTransConfiguration> db_diff_trans_config(primclex);
  BOOST_CHECK_EQUAL(true, true);

  // Open DiffTransConfiguration database
  db_diff_trans_config.open();
  BOOST_CHECK_EQUAL(db_diff_trans_config.size(), 0);

  // Make DiffTransConfiguration enumerator and enumerate configs
  auto enum_ptr = te.enumerator(tc.config, to.diff_trans_orbits[0]);

  //std::cout << "skipping DiffTransConfigEnumOccPerturbations dependent parts" << std::endl;
  for(const auto &diff_trans_config : *enum_ptr) {
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
