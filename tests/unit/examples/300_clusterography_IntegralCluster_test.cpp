#include "gtest/gtest.h"

#include "casm/clusterography/IntegralCluster.hh"
#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Superlattice.hh"

#include "crystallography/TestStructures.hh" // for test::ZrO_prim

TEST(ExampleClusterographyIntegralCluster, IntegralClusterConstructor) {

  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // Initialize an empty IntegralCluster with a reference to the prim structure
  CASM::IntegralCluster cluster {*shared_prim};

  // UnitCellCoord can be emplaced back into the cluster elements vector
  EXPECT_EQ(cluster.elements().size(), 0);
  cluster.elements().emplace_back(0, 0, 0, 0); // UnitCellCoord {0, 0, 0, 0}
  cluster.elements().emplace_back(0, 1, 0, 0); // UnitCellCoord {0, 1, 0, 0}
  EXPECT_EQ(cluster.elements().size(), 2);

  // IntegralCluster site coordinates can be obtained
  EXPECT_TRUE(cluster.coordinate(0).almost_equal(cluster.elements()[0].coordinate(*shared_prim)));


  // IntegralCluster invariants are cluster size and site-to-site distances
  ClusterInvariants invariants {cluster};
  EXPECT_EQ(invariants.size(), 2);
  EXPECT_EQ(invariants.displacement().size(), 1);

  EXPECT_EQ(
    invariants.displacement().front(),
    (cluster.coordinate(0) - cluster.coordinate(1)).as_vec(CASM::CART).norm());
  EXPECT_EQ(invariants.displacement().front(), cluster.coordinate(0).dist(cluster.coordinate(1)));

  // There is also a variant that calculates site-to-site distances using the minimum periodic
  //   boundary distance between any site images.
  CASM::xtal::Superlattice superlattice {shared_prim->lattice(), Eigen::Matrix3l::Identity() * 3};
  WithinScelClusterInvariants within_scel_invariants {cluster, superlattice.transformation_matrix_to_super()};

}
