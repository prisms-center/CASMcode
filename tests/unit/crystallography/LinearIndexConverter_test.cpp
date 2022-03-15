#include "casm/crystallography/LinearIndexConverter.hh"

#include "casm/crystallography/Lattice.hh"
#include "gtest/gtest.h"

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

xtal::UnitCellCoordIndexConverter::matrix_type transformation_matrix() {
  Eigen::Matrix3l transformation_matrix;
  transformation_matrix << 1, 9, 3, 5, 1, -2, 1, 22, 50;
  return transformation_matrix;
}

}  // namespace

TEST(UnitCellCoordIndexConverterTest, construct_via_transformation) {
  auto trans_mat = transformation_matrix();
  xtal::UnitCellCoordIndexConverter ix_bijk_converter(trans_mat, 10);
}

TEST(UnitCellCoordIndexConverterTest, construct_via_bad_basis_sites) {
  auto trans_mat = transformation_matrix();
  bool good_catch = false;
  int total_basis_atoms = 0;

  try {
    xtal::UnitCellCoordIndexConverter ix_bijk_converter(trans_mat,
                                                        total_basis_atoms);
  }

  catch (const std::runtime_error &e) {
    good_catch = true;
  }

  EXPECT_TRUE(good_catch) << "Constructed with " << total_basis_atoms
                          << " basis atoms";
}

TEST(UnitCellCoordIndexConverterTest, get_index) {
  Eigen::Matrix3l trans_mat;
  trans_mat << 10, 0, 0, 0, 10, 0, 0, 0, 10;
  xtal::UnitCellCoordIndexConverter ix_bijk_converter(trans_mat, 10);

  auto ix = ix_bijk_converter(0);
  EXPECT_EQ(ix, xtal::UnitCellCoord(0, 0, 0, 0));

  ix = ix_bijk_converter(trans_mat.determinant() * 10 - 1);
  EXPECT_EQ(ix, xtal::UnitCellCoord(9, 9, 9, 9));
}

TEST(UnitCellCoordIndexConverterTest, self_consistency) {
  auto trans_mat = transformation_matrix();
  int total_basis_atoms = 10;
  xtal::UnitCellCoordIndexConverter ix_bijk_converter(trans_mat,
                                                      total_basis_atoms);

  for (int i = 0; i < std::abs(trans_mat.determinant()) * total_basis_atoms;
       ++i) {
    auto bijk = ix_bijk_converter(i);
    auto ix_recover = ix_bijk_converter(bijk);
    auto bijk_recover = ix_bijk_converter(ix_recover);
    EXPECT_EQ(bijk, bijk_recover);
  }
}

TEST(UnitCellCoordIndexConverterTest, caching_consistency) {
  auto trans_mat = transformation_matrix();
  int total_basis_atoms = 10;
  xtal::UnitCellCoordIndexConverter ix_bijk_converter(trans_mat,
                                                      total_basis_atoms);

  // This guy falls outside the superlattice
  xtal::UnitCellCoord bijk_outside(9, 1000, 1000, 1000);

  auto ix_within = ix_bijk_converter(bijk_outside);
  auto bijk_within = ix_bijk_converter(ix_within);
  auto ix_within_recover = ix_bijk_converter(bijk_within);

  EXPECT_EQ(ix_within, ix_within_recover)
      << ix_within << " vs " << ix_within_recover;
}

TEST(UnitCellCoordIndexConverterTest, bad_basis_site_index_query) {
  auto trans_mat = transformation_matrix();
  bool good_catch = false;
  int total_basis_atoms = 10;

  xtal::UnitCellCoordIndexConverter ix_bijk_converter(trans_mat,
                                                      total_basis_atoms);

  try {
    auto ix = ix_bijk_converter(xtal::UnitCellCoord(10, 0, 0, 0));
  } catch (const std::runtime_error &e) {
    good_catch = true;
  }

  EXPECT_TRUE(good_catch) << "Constructed with " << total_basis_atoms
                          << " basis atoms";
}

TEST(UnitCellCoordIndexConverterTest, bad_index_query) {
  auto trans_mat = transformation_matrix();
  bool good_catch = false;
  int total_basis_atoms = 10;

  xtal::UnitCellCoordIndexConverter ix_bijk_converter(trans_mat,
                                                      total_basis_atoms);

  try {
    auto bijk = ix_bijk_converter(1000000);
  } catch (const std::runtime_error &e) {
    good_catch = true;
  }

  EXPECT_TRUE(good_catch) << "Constructed with " << total_basis_atoms
                          << " basis atoms";
}

TEST(UnitCellCoordIndexConverterTest, bad_fall_outside_superlattice_query) {
  auto trans_mat = transformation_matrix();
  bool good_catch = false;
  int total_basis_atoms = 10;

  xtal::UnitCellCoordIndexConverter ix_bijk_converter(trans_mat,
                                                      total_basis_atoms);
  ix_bijk_converter.never_bring_within();

  try {
    auto ix = ix_bijk_converter(xtal::UnitCellCoord(1000, 1000, 1000, 0));
  } catch (const std::runtime_error &e) {
    good_catch = true;
  }

  EXPECT_TRUE(good_catch) << "Constructed with " << total_basis_atoms
                          << " basis atoms";
}
