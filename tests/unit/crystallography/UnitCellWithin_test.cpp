#include "gtest/gtest.h"
#include <stdexcept>

#include "casm/crystallography/UnitCellWithin.hh"

#include "casm/crystallography/Lattice.hh"
#include "casm/external/Eigen/src/Core/Matrix.h"

using namespace CASM;

namespace {
  xtal::Lattice fcc_lattice() {
    double a = 1.75;
    Eigen::Matrix3d fcc_column_matrix;
    fcc_column_matrix << 0, a, a, a, 0, a, a, a, 0;
    return xtal::Lattice(fcc_column_matrix);
  }

  xtal::Lattice bcc_lattice() {
    double a = 2.8;
    Eigen::Matrix3d bcc_column_matrix;
    bcc_column_matrix << a, 0, 0, 0, a, 0, 0, 0, a;
    return xtal::Lattice(bcc_column_matrix);
  }

  xtal::Lattice hcp_lattice();
} // namespace

TEST(UnitCellWithinTest, construct_via_int_transformation) {
  Eigen::Matrix3i transformation_matrix;
  transformation_matrix << 1, 0, 3, 1, 1, -2, 1, 2, 0;

  xtal::UnitCellWithin bring_within(transformation_matrix);
}

TEST(UnitCellWithinTest, construct_via_transformation) {
  Eigen::Matrix3l transformation_matrix;
  transformation_matrix << 1, 0, 3, 1, 1, -2, 1, 2, 0;

  xtal::UnitCellWithin bring_within(transformation_matrix);
}

TEST(UnitCellWithinTest, construct_via_bad_transformation) {
  Eigen::Matrix3l transformation_matrix;
  transformation_matrix << 0, 0, 0, 1, 1, -2, 1, 2, 0;

  bool good_catch = false;
  try {
    xtal::UnitCellWithin bring_within(transformation_matrix);
  }
  catch(const std::runtime_error &e) {
    good_catch = true;
  }

  EXPECT_TRUE(good_catch) << "Determinant is " << transformation_matrix.determinant();
}

TEST(UnitCellWithinTest, construct_via_superlattice) {
  Eigen::Matrix3l transformation_matrix;
  transformation_matrix << 1, 0, 3, 1, 1, -2, 1, 2, 0;

  auto fcc_superlattice = xtal::make_superlattice(::fcc_lattice(), transformation_matrix.cast<int>());

  xtal::UnitCellWithin bring_within(::fcc_lattice(), fcc_superlattice);
}

TEST(UnitCellWithinTest, construct_via_bad_superlattice) {
  Eigen::Matrix3l transformation_matrix;
  transformation_matrix << 1, 0, 3, 1, 1, -2, 1, 2, 0;

  auto fcc_superlattice = xtal::make_superlattice(::fcc_lattice(), transformation_matrix.cast<int>());

  bool good_catch = false;
  try {
    xtal::UnitCellWithin bring_within(::bcc_lattice(), fcc_superlattice);
  }
  catch(const std::runtime_error &e) {
    good_catch = true;
  }

  EXPECT_TRUE(good_catch);
}

