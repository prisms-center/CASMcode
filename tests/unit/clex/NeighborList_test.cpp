#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/clex/NeighborList.hh"

/// What is being used to test it:

#include "Common.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(NeighborListTest)

BOOST_AUTO_TEST_CASE(PrimNeighborListBasics) {
  Structure prim(test::FCC_ternary_prim());
  std::set<int> sublat_indices;
  for(int i = 0; i < prim.basis.size(); i++) {
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
  BOOST_CHECK_EQUAL(nlist.weight_matrix() == W, true);

  // expand
  std::set<UnitCellCoord> nbors;
  nbors.insert(UnitCellCoord(0, UnitCell(3, 0, 0)));
  nlist.expand(nbors.begin(), nbors.end());

  // size
  BOOST_CHECK_EQUAL(nlist.size(), 177);

  // copy
  PrimNeighborList nlist2 = nlist;
  BOOST_CHECK_EQUAL(nlist2.size(), 177);

  // clone
  notstd::cloneable_ptr<PrimNeighborList> ptr1(nlist);
  BOOST_CHECK_EQUAL(ptr1->size(), 177);
  std::unique_ptr<PrimNeighborList> ptr2 = nlist.clone();
  BOOST_CHECK_EQUAL(ptr2->size(), 177);

}

BOOST_AUTO_TEST_CASE(SuperNeighborListBasics) {
  Structure prim(test::FCC_ternary_prim());
  std::set<int> sublat_indices;
  for(int i = 0; i < prim.basis.size(); i++) {
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
  nbors.insert(UnitCellCoord(0, UnitCell(3, 0, 0)));
  nlist.expand(nbors.begin(), nbors.end());

  // size
  BOOST_CHECK_EQUAL(nlist.size(), 177);

  // construct SuperNeighborList
  Eigen::Matrix3i T;
  T << 2, 0, 0,
  0, 2, 0,
  0, 0, 2;
  Lattice super_lat = make_supercell(prim.lattice(), T);
  PrimGrid grid(prim.lattice(), super_lat);
  SuperNeighborList super_nlist(grid, nlist);

  // size
  for(int i = 0; i < grid.size(); ++i) {
    BOOST_CHECK_EQUAL(super_nlist.sites(i).size(), 177);
  }

  for(int i = 0; i < grid.size(); ++i) {
    BOOST_CHECK_EQUAL(super_nlist.unitcells(i).size(), 177);
  }

  // overlaps
  BOOST_CHECK_EQUAL(super_nlist.overlaps(), true);

  // copy
  SuperNeighborList super_nlist2 = super_nlist;

  // clone
  notstd::cloneable_ptr<SuperNeighborList> ptr1(super_nlist);
  std::unique_ptr<SuperNeighborList> ptr2 = super_nlist.clone();
  notstd::cloneable_ptr<SuperNeighborList> ptr3 = ptr1;

}

BOOST_AUTO_TEST_CASE(NeighborListTestLatticeTests) {

  Eigen::Matrix3d latvec;
  latvec.col(0) <<  2.955270000000, 0.000000000000, 0.000000000000;
  latvec.col(1) <<  1.477635000000, 2.559338895042, 0.000000000000;
  latvec.col(2) << 1.477635000000, 0.853112965014, 11.758280000000;

  PrimNeighborList::Matrix3Type W = PrimNeighborList::make_weight_matrix(latvec, 10, TOL);

  PrimNeighborList::Matrix3Type W_check;
  W_check.row(0) << 2, 1, 1;
  W_check.row(1) << 1, 2, 1;
  W_check.row(2) << 1, 1, 32;

  BOOST_CHECK_EQUAL(W == W_check, true);
}

BOOST_AUTO_TEST_CASE(Proj) {

  test::FCCTernaryProj proj;
  proj.check_init();

  PrimClex primclex(proj.dir, null_log());
  primclex.settings().set_crystallography_tol(TOL);
  double tol = primclex.crystallography_tol();
  Structure prim = primclex.get_prim();

  const ProjectSettings &set = primclex.settings();

  // initialize nlist
  PrimNeighborList nlist(
    set.nlist_weight_matrix(),
    set.nlist_sublat_indices().begin(),
    set.nlist_sublat_indices().end()
  );

  // generate orbitree
  SiteOrbitree tree(prim.lattice(), tol);
  jsonParser bspecs_json(proj.bspecs());
  tree = make_orbitree(prim, bspecs_json, tol);

  // expand the nlist to contain 'tree'
  std::set<UnitCellCoord> nbors;
  neighborhood(std::inserter(nbors, nbors.begin()), tree, prim, TOL);

  //std::cout << "expand nlist" << std::endl;
  nlist.expand(nbors.begin(), nbors.end());
  BOOST_CHECK_EQUAL(nlist.size(), 177);

  //std::cout << "expand nlist again" << std::endl;
  nbors.clear();
  nbors.insert(UnitCellCoord(0, UnitCell(4, 0, 0)));
  nlist.expand(nbors.begin(), nbors.end());
  BOOST_CHECK_EQUAL(nlist.size(), 381);

}

BOOST_AUTO_TEST_SUITE_END()
