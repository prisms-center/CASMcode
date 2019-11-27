#include "gtest/gtest.h"

#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/crystallography/Lattice.hh"

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

  xtal::LinearIndexConverter::matrix_type transformation_matrix() {
    Eigen::Matrix3l transformation_matrix;
    transformation_matrix << 1, 9, 3, 5, 1, -2, 1, 22, 50;
    return transformation_matrix;
  }

} // namespace

TEST(LinearIndexConverterTest, construct_via_transformation) {
  auto trans_mat = transformation_matrix();
  xtal::LinearIndexConverter ix_bijk_converter(trans_mat, 10);
}

TEST(LinearIndexConverterTest, construct_via_superlattice) {
  auto trans_mat = transformation_matrix();
  auto fcc_superlattice = xtal::make_superlattice(::fcc_lattice(), trans_mat.cast<int>());
  xtal::LinearIndexConverter ix_bijk_converter(::fcc_lattice(), fcc_superlattice, 3);
}

TEST(LinearIndexConverterTest, construct_via_bad_basis_sites) {
  auto trans_mat = transformation_matrix();
  bool good_catch = false;
  int total_basis_atoms = 0;

  try {
    xtal::LinearIndexConverter ix_bijk_converter(trans_mat, total_basis_atoms);
  }

  catch(const std::runtime_error &e) {
    good_catch = true;
  }

  EXPECT_TRUE(good_catch) << "Constructed with " << total_basis_atoms << " basis atoms";
}

TEST(LinearIndexConverterTest, get_index) {
  auto trans_mat = transformation_matrix();
  xtal::LinearIndexConverter ix_bijk_converter(trans_mat, 10);

  auto ix = ix_bijk_converter[0];
  EXPECT_EQ(ix, xtal::UnitCellCoord(0, 0, 0, 0));
}
