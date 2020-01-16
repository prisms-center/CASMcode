#include "gtest/gtest.h"
#include <stdexcept>

/// What is being tested:
#include "casm/crystallography/UnitCellCoord.hh"

/// What is being used to test it:
#include "FCCTernaryProj.hh"
#include "casm/external/Eigen/src/Core/Matrix.h"
#include "casm/misc/CASM_Eigen_math.hh"

using namespace CASM;

namespace {
  Eigen::Matrix3i transformation_matrix() {
    Eigen::Matrix3i trans_mat;
    trans_mat << -1, 1, 2, 1, -1, 2, 1, 1, -2;

    return trans_mat;
  }

  double fcc_lattice_parameter() {
    return 1.75;
  }

  xtal::Lattice fcc_lattice() {
    double a = fcc_lattice_parameter();
    Eigen::Matrix3d fcc_column_matrix;
    fcc_column_matrix << 0, a, a, a, 0, a, a, a, 0;
    return xtal::Lattice(fcc_column_matrix);
  }

  xtal::Lattice fcc_superlattice() {
    auto fcc_superlattice = xtal::make_superlattice(::fcc_lattice(), ::transformation_matrix());
    return fcc_superlattice;
  }

  xtal::Lattice fcc_conventional() {
    Eigen::Matrix3i trans_mat;
    trans_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    auto fcc_superlattice = xtal::make_superlattice(::fcc_lattice(), trans_mat);
    return fcc_superlattice;
  }

} // namespace

TEST(UnitCellCoordTest, Test1) {

  BasicStructure<Site> prim = test::FCC_ternary_prim();
  {
    UnitCellCoord uccoord(0, -1, 1, 1);
    Eigen::Vector3d vec = uccoord.site(prim).cart();
    EXPECT_EQ(almost_equal(vec, Eigen::Vector3d(4., 0., 0.)), true);
  }
}

TEST(UnitCellCoordTest, construct_from_coordinate) {
  Eigen::Vector3d exact_cart = fcc_lattice().lat_column_mat() * Eigen::Vector3d(2, 0, 5);
  Coordinate exact_coord(exact_cart, fcc_lattice(), CART);
  EXPECT_EQ(UnitCell::from_coordinate(exact_coord), UnitCell(2, 0, 5)) << "Failed with exact coordinate";

  Eigen::Vector3d fuzzy_cart = exact_cart - Eigen::Vector3d(0.00001, -0.00002, -0.000005);
  Coordinate fuzzy_coord(fuzzy_cart, fcc_lattice(), CART);
  EXPECT_EQ(UnitCell::from_coordinate(fuzzy_coord), UnitCell(2, 0, 5)) << "Failed with fuzzy coordinate";

  Eigen::Vector3d incompatible_cart = exact_cart * 1.2;
  Coordinate incompatible_coord(incompatible_cart, fcc_lattice(), CART);

  bool good_catch = false;
  try {
    UnitCell impossible_uc = UnitCell::from_coordinate(incompatible_coord);
  }
  catch(const std::runtime_error &e) {
    good_catch = true;
  }
  EXPECT_TRUE(good_catch) << "Exception misbehaved";
}

TEST(UnitCellCoordTest, construct_from_cartesian) {
  Lattice fcc_lat =::fcc_lattice();
  Eigen::Vector3d cart_coord = 2 * fcc_lat[0] - 5 * fcc_lat[1] + 9 * fcc_lat[2];
  UnitCell uc_from_cart = UnitCell::from_cartesian(cart_coord, fcc_lat);
  EXPECT_EQ(uc_from_cart, UnitCell(2, -5, 9));
}

TEST(UnitCellCoordTest, recover_coordinate) {
  double a = fcc_lattice_parameter();

  UnitCell ucc(2, 0, 5);
  Eigen::Vector3d exact_cart = fcc_lattice().lat_column_mat() * Eigen::Vector3d(2, 0, 5);
  Coordinate exact_coord(exact_cart, fcc_lattice(), CART);
  Coordinate recovered_coord = ucc.coordinate(fcc_lattice());

  EXPECT_TRUE(exact_coord.almost_equal(recovered_coord));
}

TEST(UnitCellCoordTest, reset_tiling_unit) {
  Lattice conventional_lattice = fcc_conventional();
  Lattice stacked_conventional_lattice = fcc_superlattice();

  UnitCell relative_to_conventional_uc(4, 4, 4);
  UnitCell relative_to_stacked_uc = relative_to_conventional_uc.reset_tiling_unit(conventional_lattice, stacked_conventional_lattice);
  EXPECT_EQ(relative_to_stacked_uc, UnitCell(4, 4, 2));

  relative_to_stacked_uc = UnitCell(4, 4, 4);
  relative_to_conventional_uc = relative_to_stacked_uc.reset_tiling_unit(stacked_conventional_lattice, conventional_lattice);
  EXPECT_EQ(relative_to_conventional_uc, UnitCell(4, 4, 8));

  bool good_catch = false;
  relative_to_conventional_uc = UnitCell(4, 4, 3);
  try {
    relative_to_stacked_uc = relative_to_conventional_uc.reset_tiling_unit(conventional_lattice, stacked_conventional_lattice);
  }
  catch(const std::runtime_error &e) {
    good_catch = true;
  }
  EXPECT_TRUE(good_catch);
}
