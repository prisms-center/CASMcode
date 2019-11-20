#include "gtest/gtest.h"

/// What is being tested:
#include "casm/clex/NeighborList.hh"

/// What is being used to test it:

#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/crystallography/Structure.hh"
#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/IntegralCluster_impl.hh"

using namespace CASM;

TEST(NeighborListTest, PrimNeighborListBasics) {
  Structure prim(test::FCC_ternary_prim());
  std::set<int> sublat_indices;
  for(int i = 0; i < prim.basis().size(); i++) {
    sublat_indices.insert(i);
  }

  // construct
  PrimNeighborList nlist(
    PrimNeighborList::make_weight_matrix(prim.lattice().lat_column_mat(), 10, TOL),
    sublat_indices.begin(),
    sublat_indices.end()
  );

  // weight matrix
  Eigen::Matrix3l W;
  W << 2, 1, 1,
  1, 2, 1,
  1, 1, 2;
  EXPECT_EQ(nlist.weight_matrix() == W, true);

  // expand
  std::set<UnitCellCoord> nbors;
  nbors.emplace(0, UnitCell(3, 0, 0));
  nlist.expand(nbors.begin(), nbors.end());

  // size
  EXPECT_EQ(nlist.size(), 177);

  // copy
  PrimNeighborList nlist2 = nlist;
  EXPECT_EQ(nlist2.size(), 177);

  // clone
  notstd::cloneable_ptr<PrimNeighborList> ptr1(nlist.clone());
  EXPECT_EQ(ptr1->size(), 177);
  std::unique_ptr<PrimNeighborList> ptr2 = nlist.clone();
  EXPECT_EQ(ptr2->size(), 177);

}

TEST(NeighborListTest, SuperNeighborListBasics) {
  Structure prim(test::FCC_ternary_prim());
  std::set<int> sublat_indices;
  for(int i = 0; i < prim.basis().size(); i++) {
    sublat_indices.insert(i);
  }

  // construct
  PrimNeighborList nlist(
    PrimNeighborList::make_weight_matrix(prim.lattice().lat_column_mat(), 10, TOL),
    sublat_indices.begin(),
    sublat_indices.end()
  );

  // expand
  std::set<UnitCellCoord> nbors;
  nbors.emplace(0, UnitCell(3, 0, 0));
  nlist.expand(nbors.begin(), nbors.end());

  // size
  EXPECT_EQ(nlist.size(), 177);

  // construct SuperNeighborList
  Eigen::Matrix3i T;
  T << 2, 0, 0,
  0, 2, 0,
  0, 0, 2;
  Lattice super_lat = make_superlattice(prim.lattice(), T);
  PrimGrid grid(prim.lattice(), super_lat);
  SuperNeighborList super_nlist(grid, nlist);

  // size
  for(int i = 0; i < grid.size(); ++i) {
    EXPECT_EQ(super_nlist.sites(i).size(), 177);
  }

  for(int i = 0; i < grid.size(); ++i) {
    EXPECT_EQ(super_nlist.unitcells(i).size(), 177);
  }

  // overlaps
  EXPECT_EQ(super_nlist.overlaps(), true);

  // copy
  SuperNeighborList super_nlist2 = super_nlist;

  // clone
  notstd::cloneable_ptr<SuperNeighborList> ptr1(super_nlist.clone());
  std::unique_ptr<SuperNeighborList> ptr2 = super_nlist.clone();
  notstd::cloneable_ptr<SuperNeighborList> ptr3 = ptr1;

}

TEST(NeighborListTest, NeighborListTestLatticeTests) {

  Eigen::Matrix3d latvec;
  latvec.col(0) <<  2.955270000000, 0.000000000000, 0.000000000000;
  latvec.col(1) <<  1.477635000000, 2.559338895042, 0.000000000000;
  latvec.col(2) << 1.477635000000, 0.853112965014, 11.758280000000;

  PrimNeighborList::Matrix3Type W = PrimNeighborList::make_weight_matrix(latvec, 10, TOL);

  PrimNeighborList::Matrix3Type W_check;
  W_check.row(0) << 2, 1, 1;
  W_check.row(1) << 1, 2, 1;
  W_check.row(2) << 1, 1, 32;

  EXPECT_EQ(W == W_check, true);
}

TEST(NeighborListTest, Proj) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  primclex.settings().set_crystallography_tol(TOL);
  double tol = primclex.crystallography_tol();
  Structure prim = primclex.prim();
  const ProjectSettings &set = primclex.settings();

  // initialize nlist
  PrimNeighborList nlist(
    set.nlist_weight_matrix(),
    set.nlist_sublat_indices().begin(),
    set.nlist_sublat_indices().end()
  );

  // generate orbitree
  jsonParser bspecs_json(proj.bspecs());
  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  double crystallography_tol = TOL;
  make_prim_periodic_orbits(prim,
                            bspecs_json,
                            alloy_sites_filter,
                            crystallography_tol,
                            std::back_inserter(orbits),
                            null_log());

  // expand the nlist to contain 'tree'
  std::set<UnitCellCoord> nbors;
  prim_periodic_neighborhood(orbits.begin(), orbits.end(), std::inserter(nbors, nbors.begin()));

  //std::cout << "expand nlist" << std::endl;
  nlist.expand(nbors.begin(), nbors.end());
  EXPECT_EQ(nlist.size(), 177);

  //std::cout << "expand nlist again" << std::endl;
  nbors.clear();
  nbors.emplace(0, UnitCell(4, 0, 0));
  nlist.expand(nbors.begin(), nbors.end());
  EXPECT_EQ(nlist.size(), 381);

}
