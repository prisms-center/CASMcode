#include "gtest/gtest.h"

/// What is being tested:
#include "casm/symmetry/InvariantSubgroup.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"
#include "casm/symmetry/Orbit.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/symmetry/SymBasisPermute.hh"

/// What is being used to test it:
#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "TestConfiguration.hh"
#include "casm/casm_io/json_io/jsonFile.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/symmetry/Orbit_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/kinetics/DiffusionTransformation_impl.hh"
#include "casm/symmetry/ConfigSubOrbits_impl.hh"

using namespace CASM;

namespace {
  struct TestConfig0 : test::TestConfiguration {

    TestConfig0(const PrimClex &primclex) :
      TestConfiguration(
        primclex,
        Eigen::Vector3i(2, 1, 1).asDiagonal(), {
      0, 0,  0, 0,  1, 1,  0, 0
    }) {

      EXPECT_EQ(this->scel_fg().size(), 16);
      EXPECT_EQ(this->config_sym_fg().size(), 8);
    }

  };
}


TEST(OrbitTest, Test0) {
  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();

  const auto &g = prim.factor_group();

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 0, 0, 0, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    EXPECT_EQ(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 1, 0, 0, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    EXPECT_EQ(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 2, 0, 0, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    EXPECT_EQ(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 2, 0, 0, 0));
    generating_element.elements().push_back(UnitCellCoord(prim, 3, 0, 0, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    EXPECT_EQ(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().push_back(UnitCellCoord(prim, 0, 0, 0, 0));
    generating_element.elements().push_back(UnitCellCoord(prim, 1, 0, 1, 0));
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(primclex);
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g, sym_compare);
    EXPECT_EQ(orbit.size(), 6);
  }


  // --- DiffTrans tests ---

  {
    using namespace Kinetics;

    // Construct
    DiffusionTransformation diff_trans(prim);
    EXPECT_EQ(true, true);
    EXPECT_EQ(diff_trans.occ_transform().size(), 0);

    UnitCellCoord uccoordA(prim, 2, 0, 0, 0);
    UnitCellCoord uccoordB(prim, 3, 0, 0, 0);
    Index iVa = 0;
    Index iO = 1;

    // Add transform (so that it's not sorted as constructed)
    diff_trans.occ_transform().emplace_back(uccoordB, iO, iVa);
    diff_trans.occ_transform().emplace_back(uccoordA, iVa, iO);
    EXPECT_EQ(true, true);
    EXPECT_EQ(diff_trans.is_valid_occ_transform(), true);

    // Add trajectory
    diff_trans.species_traj().emplace_back(SpeciesLocation(uccoordA, iVa, 0), SpeciesLocation(uccoordB, iVa, 0));
    diff_trans.species_traj().emplace_back(SpeciesLocation(uccoordB, iO, 0), SpeciesLocation(uccoordA, iO, 0));
    EXPECT_EQ(true, true);
    EXPECT_EQ(diff_trans.is_valid_species_traj(), true);
    EXPECT_EQ(diff_trans.is_valid(), true);

    PrimPeriodicDiffTransSymCompare sym_compare(primclex);
    PrimPeriodicDiffTransOrbit orbit(diff_trans, g, sym_compare, &primclex);
    EXPECT_EQ(orbit.size(), 2);
  }
}

TEST(OrbitTest, Test1) {

  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const Lattice &lat = prim.lattice();
  Supercell prim_scel(&primclex, Eigen::Matrix3i::Identity());
  EXPECT_EQ(true, true);

  // Make PrimPeriodicIntegralClusterOrbit
  jsonFile bspecs {autotools::abs_srcdir() + "/tests/unit/kinetics/ZrO_bspecs_0.json"};

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());
  EXPECT_EQ(true, true);

  EXPECT_EQ(orbits.size(), 74);
}
