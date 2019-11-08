#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/database/DiffTransConfigDatabase.hh"
#include "casm/database/json/jsonDatabase.hh"


/// What is being used to test it:

#include "Common.hh"
#include "TestConfiguration.hh"
#include "TestOrbits.hh"
#include "FCCTernaryProj.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/clex/ConfigEnumAllOccupations_impl.hh"
#include "casm/clex/ScelEnum_impl.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/casm_io/json/jsonFile.hh"
#include "casm/casm_io/container/stream_io.hh"
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
        jsonFile(autotools::abs_srcdir() + "/tests/unit/kinetics/FCCTernary_bspecs_0.json"),
        2, 4) {}
  };

  struct TestEnumerator0 {

    TestEnumerator0() :
      diff_perturb_specs(autotools::abs_srcdir() + "/tests/unit/kinetics/FCCTernary_diff_perturb_0.json") {}

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


TEST(jsonDiffTransConfigDatabase_Test, Test1) {

  // Make test project
  test::FCCTernaryProj proj;
  proj.check_init();
  PrimClex primclex(proj.dir, null_log());

  TestOrbits0 to(primclex);
  TestEnumerator0 te;
  TestConfig0 tc(primclex);

  std::cout << "point 1\n";
  // Make DiffTransConfiguration database
  DB::jsonDatabase<Kinetics::DiffTransConfiguration> db_diff_trans_config(primclex);
  EXPECT_EQ(true, true);

  std::cout << "point 2\n";
  // Open DiffTransConfiguration database
  db_diff_trans_config.open();
  EXPECT_EQ(db_diff_trans_config.size(), 0);

  std::cout << "point 3\n";
  // Make DiffTransConfiguration enumerator and enumerate configs
  auto enum_ptr = te.enumerator(tc.config, to.diff_trans_orbits[0]);

  //std::cout << "skipping DiffTransConfigEnumOccPerturbations dependent parts" << std::endl;
  for(const auto &diff_trans_config : *enum_ptr) {
    db_diff_trans_config.insert(diff_trans_config);
  }
  db_diff_trans_config.commit();
  EXPECT_EQ(db_diff_trans_config.size(), 29); // not checked for accuracy
  std::cout << "point 4\n";

  // Check cached properties
  std::cout << "skipping cache check" << std::endl;
  ////  for(const auto &diff_trans_config : db_diff_trans_config) {
  ////    //std::cout << "id: " << config.id() << "  occ: " << config.occupation() << std::endl;
  ////    EXPECT_EQ(diff_trans_config.cache().contains("multiplicity"), false);
  ////    EXPECT_EQ(diff_trans_config.multiplicity() != 0, true);
  ////    EXPECT_EQ(diff_trans_config.cache().contains("multiplicity"), true);
  ////  }
  ////  db_diff_trans_config.commit();
  //
  // Check that the database is sorted
  {
    auto next = db_diff_trans_config.begin();
    auto it = next++;
    auto end = db_diff_trans_config.end();
    for(; next != end; ++it, ++next) {
      EXPECT_EQ(*it < *next, true);
    }
  }
  std::cout << "point 5\n";
  // Close DiffTransConfiguration database
  db_diff_trans_config.close();
  EXPECT_EQ(db_diff_trans_config.size(), 0);

  // Re-open DiffTransConfiguration database
  db_diff_trans_config.open();
  EXPECT_EQ(db_diff_trans_config.size(), 29); // not checked for accuracy
  std::cout << "point 6\n";
  //  // Check cached properties
  std::cout << "skipping cache check" << std::endl;
  ////  for(const auto &diff_trans_config : db_diff_trans_config) {
  ////    EXPECT_EQ(diff_trans_config.cache().contains("multiplicity"), true);
  ////  }

  // Check that the database is sorted
  {
    auto next = db_diff_trans_config.begin();
    auto it = next++;
    auto end = db_diff_trans_config.end();
    for(; next != end; ++it, ++next) {
      EXPECT_EQ(*it < *next, true);
    }
  }

  // Close DiffTransConfiguration database
  db_diff_trans_config.close();
  EXPECT_EQ(true, true);
}
