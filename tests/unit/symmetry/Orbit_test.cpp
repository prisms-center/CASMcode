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
#include "TestConfiguration.hh"
#include "ZrOProj.hh"
#include "casm/casm_io/json/jsonFile.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/symmetry/ConfigSubOrbits_impl.hh"
#include "casm/symmetry/Orbit_impl.hh"

using namespace CASM;

namespace {
struct TestConfig0 : test::TestConfiguration {
  TestConfig0(const PrimClex &primclex)
      : TestConfiguration(primclex, Eigen::Vector3l(2, 1, 1).asDiagonal(),
                          test::eigen_vector<int>({0, 0, 0, 0, 1, 1, 0, 0})) {
    EXPECT_EQ(this->scel_fg().size(), 16);
    EXPECT_EQ(this->config_sym_fg().size(), 8);
  }
};
}  // namespace

TEST(OrbitTest, Test0) {
  test::ZrOProj proj;
  proj.check_init();

  ScopedNullLogging logging;
  PrimClex primclex(proj.dir);
  const Structure &prim = primclex.prim();

  const auto &g = prim.factor_group();

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().emplace_back(0, 0, 0, 0);
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(
        primclex.shared_prim(), primclex.crystallography_tol());
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g,
                                             sym_compare);
    EXPECT_EQ(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().emplace_back(1, 0, 0, 0);
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(
        primclex.shared_prim(), primclex.crystallography_tol());
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g,
                                             sym_compare);
    EXPECT_EQ(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().emplace_back(2, 0, 0, 0);
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(
        primclex.shared_prim(), primclex.crystallography_tol());
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g,
                                             sym_compare);
    EXPECT_EQ(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().emplace_back(2, 0, 0, 0);
    generating_element.elements().emplace_back(3, 0, 0, 0);
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(
        primclex.shared_prim(), primclex.crystallography_tol());
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g,
                                             sym_compare);
    EXPECT_EQ(orbit.size(), 2);
  }

  {
    IntegralCluster generating_element(prim);
    generating_element.elements().emplace_back(0, 0, 0, 0);
    generating_element.elements().emplace_back(1, 0, 1, 0);
    PrimPeriodicSymCompare<IntegralCluster> sym_compare(
        primclex.shared_prim(), primclex.crystallography_tol());
    PrimPeriodicOrbit<IntegralCluster> orbit(generating_element, g,
                                             sym_compare);
    EXPECT_EQ(orbit.size(), 6);
  }
}
