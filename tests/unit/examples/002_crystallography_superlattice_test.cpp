#include "gtest/gtest.h"

#include "casm/crystallography/Superlattice.hh"
#include "casm/misc/CASM_Eigen_math.hh" // for almost_equal of Eigen types
#include "casm/misc/CASM_math.hh" // for almost_equal of other types

TEST(ExampleCrystallographySuperlattice, SuperlatticeConstruction) {

  CASM::xtal::Lattice unit_lattice {Eigen::Matrix3d::Identity()};

  // make a 2x2x2 supercell of the unit lattice
  Eigen::Matrix3l T;
  T << 2, 0, 0,
  0, 2, 0,
  0, 0, 2;

  CASM::xtal::Lattice basic_superlattice = make_superlattice(unit_lattice, T);
  EXPECT_TRUE(almost_equal(basic_superlattice.lat_column_mat(), unit_lattice.lat_column_mat() * T.cast<double>()));

  // There is also a class, Superlattice, which can be used to hold the unit lattice,
  //   the super lattice, and the transformation matrix
  CASM::xtal::Superlattice superlattice_1 {unit_lattice, T};
  EXPECT_EQ(superlattice_1.transformation_matrix_to_super(), T);
  EXPECT_TRUE(almost_equal(superlattice_1.superlattice().lat_column_mat(), basic_superlattice.lat_column_mat()));

  // The 'size' of a Superlattice is the number of unit lattices that fit inside the super lattice
  EXPECT_EQ(superlattice_1.size(), T.determinant());

  // A Superlattice, can also be constructed with a unit lattice and a super lattice,
  //   in which case the transformation matrix is calculated
  CASM::xtal::Superlattice superlattice_2 {unit_lattice, basic_superlattice};
  EXPECT_EQ(superlattice_1.transformation_matrix_to_super(), T);

  // The Superlattice constructor will throw if the provided super lattice is not an integer
  //   multiple of the unit lattice
  CASM::xtal::Lattice not_a_superlattice = basic_superlattice;
  not_a_superlattice[0] = basic_superlattice[0] * 1.1;
  ASSERT_THROW((CASM::xtal::Superlattice {unit_lattice, not_a_superlattice}), std::runtime_error);

}
