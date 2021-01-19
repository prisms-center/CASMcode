#include "casm/clusterography/ClusterInvariants.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Superlattice.hh"
#include "crystallography/TestStructures.hh"  // for test::ZrO_prim
#include "gtest/gtest.h"

TEST(ExampleClusterographyIntegralCluster, IntegralClusterConstructor) {
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // Initialize an empty IntegralCluster with a reference to the prim structure
  CASM::IntegralCluster cluster{*shared_prim};

  // UnitCellCoord can be emplaced back into the cluster elements vector
  EXPECT_EQ(cluster.elements().size(), 0);
  cluster.elements().emplace_back(0, 0, 0, 0);  // UnitCellCoord {0, 0, 0, 0}
  cluster.elements().emplace_back(0, 1, 0, 0);  // UnitCellCoord {0, 1, 0, 0}

  // Get cluster size
  EXPECT_EQ(cluster.elements().size(), 2);

  // Get cluster site coordinates
  EXPECT_TRUE(cluster.coordinate(0).almost_equal(
      cluster.elements()[0].coordinate(*shared_prim)));
}

TEST(ExampleClusterographyIntegralCluster, ClusterInvariants) {
  // Cluster invariants, properties of the cluster that invariant to application
  // of symmetry
  //   are useful for understanding, comparing, filtering, etc. The particular
  //   properties that are invariant may be different in different contexts. For
  //   instance, when considering a system with or without periodic boundary
  //   conditions, or clusters in the infinite periodic crystal vs clusters
  //   around a particular local site, the number of sites in the cluster is
  //   always invariant, but the relevant site-to-site distance may be
  //   different.

  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // Initialize an empty IntegralCluster with a reference to the prim structure
  CASM::IntegralCluster cluster{*shared_prim};

  // UnitCellCoord can be emplaced back into the cluster elements vector
  EXPECT_EQ(cluster.elements().size(), 0);
  cluster.elements().emplace_back(0, 0, 0, 0);  // UnitCellCoord {0, 0, 0, 0}
  cluster.elements().emplace_back(0, 1, 0, 0);  // UnitCellCoord {0, 1, 0, 0}

  // The basic IntegralCluster invariants are cluster size and site-to-site
  // distances
  CASM::ClusterInvariants invariants{cluster};
  EXPECT_EQ(invariants.size(), 2);
  EXPECT_EQ(invariants.displacement().size(), 1);

  EXPECT_EQ(invariants.displacement().front(),
            (cluster.coordinate(0) - cluster.coordinate(1))
                .as_vec(CASM::CART)
                .norm());
  EXPECT_EQ(invariants.displacement().front(),
            cluster.coordinate(0).dist(cluster.coordinate(1)));

  // There is also a variant that calculates site-to-site distances using the
  // minimum periodic
  //   boundary distance between any site images.
  CASM::xtal::Superlattice superlattice{shared_prim->lattice(),
                                        Eigen::Matrix3l::Identity() * 3};
  CASM::WithinScelClusterInvariants within_scel_invariants{
      cluster, superlattice.transformation_matrix_to_super()};
}
