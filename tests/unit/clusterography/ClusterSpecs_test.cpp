#include "gtest/gtest.h"
#include "autotools.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/crystallography/Structure.hh"
#include "crystallography/TestStructures.hh"

using namespace CASM;
using namespace test;

namespace {
  Eigen::Matrix3i _fcc_conventional_transf_mat() {
    Eigen::Matrix3i transf_mat;
    transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    return transf_mat;
  }
}

class ClusterSpecsTest : public testing::Test {
protected:

  std::shared_ptr<CASM::Structure const> shared_prim;

  ClusterSpecsTest():
    shared_prim(std::make_shared<CASM::Structure const>(test::FCC_ternary_prim())) {}

};

TEST_F(ClusterSpecsTest, PeriodicTest) {

  // primitive FCC clusters

  double a = shared_prim->lattice()[0].norm();

  std::vector<double> max_length {
    0, // null cluster
    0, // point clusters
    a + TOL // 1NN pair clusters only
  };

  PeriodicMaxLengthClusterSpecs cluster_specs {
    shared_prim,
    shared_prim->factor_group(),
    alloy_sites_filter,
    max_length};

  auto orbits = cluster_specs.make_periodic_orbits(null_log());
  EXPECT_EQ(orbits.size(), 3);
  EXPECT_EQ(orbits[0].prototype().size(), 0); // 1 null cluster
  EXPECT_EQ(orbits[0].size(), 1);
  EXPECT_EQ(orbits[1].prototype().size(), 1); // 1 point cluster
  EXPECT_EQ(orbits[1].size(), 1);
  EXPECT_EQ(orbits[2].prototype().size(), 2); // 6 1NN pair clusters (per site, no double counting)
  EXPECT_EQ(orbits[2].size(), 6);
}

TEST_F(ClusterSpecsTest, WithinScelTest) {

  // conventional FCC (4 atom unit cell) clusters "within_scel"
  Supercell supercell = std::make_shared<Supercell>(shared_prim, _fcc_conventional_transf_mat());

  double a = shared_prim->lattice()[0].norm();

  std::vector<double> max_length {
    0, // null cluster
    0, // point clusters
    a + TOL // 1NN pair clusters only
  };

  WithinScelMaxLengthClusterSpecs cluster_specs {
    shared_prim,
    supercell.sym_info().transformation_matrix_to_super(),
    make_generating_group(supercell),
    alloy_sites_filter,
    max_length};

  auto orbits = cluster_specs.make_within_scel_orbits(log());

  // Construct OrbitPrinterOptions:
  OrbitPrinterOptions printer_options;
  printer_options.coord_type = INTEGRAL;
  printer_options.print_equivalence_map = true;
  printer_options.print_invariant_group = true;

  // Print prototype clusters and all equivalent clusters
  printer_options.orbit_print_mode = ORBIT_PRINT_MODE::FULL;
  print_clust(orbits.begin(), orbits.end(), log(), printer_options);

  EXPECT_EQ(orbits.size(), 3);
  EXPECT_EQ(orbits[0].prototype().size(), 0);
  EXPECT_EQ(orbits[0].size(), 1);
  EXPECT_EQ(orbits[1].prototype().size(), 1);
  EXPECT_EQ(orbits[1].size(), 1);
  EXPECT_EQ(orbits[2].prototype().size(), 2);
  EXPECT_EQ(orbits[2].size(), 6);

}
