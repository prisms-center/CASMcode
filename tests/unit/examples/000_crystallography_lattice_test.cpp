#include "gtest/gtest.h"

#include "casm/crystallography/Lattice.hh"
#include "casm/misc/CASM_Eigen_math.hh" // for almost_equal of Eigen types
#include "casm/misc/CASM_math.hh" // for almost_equal of other types

TEST(ExampleCrystallographyLattice, LatticeVectors) {

  // Construct a lattice with the lattice vectors
  Eigen::Vector3d a {1., 0., 1.};
  Eigen::Vector3d b {0., 2., 0.};
  Eigen::Vector3d c {0., 0., 3.};
  CASM::xtal::Lattice lattice {a, b, c, };

  // Check floating point equivalence with CASM::almost_equal, using default tolerance
  //   CASM::TOL = 1e-5, which is usually ok for crystallographic values
  double tol = CASM::TOL;

  // Access the lattice vectors via 'Lattice::operator[](Index i) const'
  EXPECT_TRUE(almost_equal(lattice[0], a, tol));
  EXPECT_TRUE(almost_equal(lattice[1], b, tol));
  EXPECT_TRUE(almost_equal(lattice[2], c, tol));

  // Access the lattice vectors via 'Lattice::vectors() const'
  Eigen::Vector3d a_check, b_check, c_check;
  std::tie(a_check, b_check, c_check) = lattice.vectors();
  EXPECT_TRUE(almost_equal(a_check, a, tol));
  EXPECT_TRUE(almost_equal(b_check, b, tol));
  EXPECT_TRUE(almost_equal(c_check, c, tol));

  // Access the lattice vectors as a column matrix via 'Lattice::lat_column_mat() const'
  Eigen::Matrix3d lat_column_mat;
  lat_column_mat = lattice.lat_column_mat();
  EXPECT_TRUE(almost_equal(lat_column_mat.col(0), a, tol));
  EXPECT_TRUE(almost_equal(lat_column_mat.col(1), b, tol));
  EXPECT_TRUE(almost_equal(lat_column_mat.col(2), c, tol));

}

TEST(ExampleCrystallographyLattice, LatticeConstructors) {

  // Construct a lattice with the lattice vectors
  Eigen::Vector3d a {1., 0., 1.};
  Eigen::Vector3d b {0., 2., 0.};
  Eigen::Vector3d c {0., 0., 3.};
  CASM::xtal::Lattice lattice_1 {a, b, c};

  // A tolerance used for crystallographic comparisons is included in Lattice. By default the
  // tolerance is set to CASM::TOL.
  EXPECT_EQ(lattice_1.tol(), CASM::TOL);

  // Construct a lattice with the lattice vectors, and a custom tolerance
  double custom_tol = CASM::TOL * 0.1;
  CASM::xtal::Lattice lattice_2 {a, b, c, custom_tol};
  EXPECT_EQ(lattice_2.tol(), custom_tol);

  // You can also construct a lattice with the lattice vectors as a column matrix
  Eigen::Matrix3d M;
  M << 1., 0., 0.,
  0., 2., 0.,
  1., 0., 3.;
  CASM::xtal::Lattice lattice_3 {M};
  EXPECT_TRUE(almost_equal(lattice_3[0], M.col(0)));
  EXPECT_TRUE(almost_equal(lattice_3[1], M.col(1)));
  EXPECT_TRUE(almost_equal(lattice_3[2], M.col(2)));

}

TEST(ExampleCrystallographyLattice, LatticeProperties) {

  Eigen::Vector3d vec1 {1., 1., 0.};
  Eigen::Vector3d vec2 {0., 1., 0.};
  Eigen::Vector3d vec3 {0., 0., 2.};
  CASM::xtal::Lattice lattice {vec1, vec2, vec3};

  // get lattice vector angles in degrees: lattice.angle(i) is the angle between vectors (i+1)%3, (i+2)%3
  EXPECT_TRUE(CASM::almost_equal(lattice.angle(0), 90.0));
  EXPECT_TRUE(CASM::almost_equal(lattice.angle(1), 90.0));
  EXPECT_TRUE(CASM::almost_equal(lattice.angle(2), 45.0));

  // get lattice vector lengths:
  EXPECT_TRUE(CASM::almost_equal(lattice.length(0), vec1.norm()));
  EXPECT_TRUE(CASM::almost_equal(lattice.length(1), vec2.norm()));
  EXPECT_TRUE(CASM::almost_equal(lattice.length(2), vec3.norm()));

  // get lattice volume:
  EXPECT_TRUE(CASM::almost_equal(lattice.volume(), vec1.dot(vec2.cross(vec3))));

  // get the reciprocal lattice:
  CASM::xtal::Lattice reciprocal_lattice = lattice.reciprocal();
  EXPECT_TRUE(almost_equal(reciprocal_lattice[0], 2 * M_PI * vec2.cross(vec3) / lattice.volume()));
  EXPECT_TRUE(almost_equal(reciprocal_lattice[1], 2 * M_PI * vec3.cross(vec1) / lattice.volume()));
  EXPECT_TRUE(almost_equal(reciprocal_lattice[2], 2 * M_PI * vec1.cross(vec2) / lattice.volume()));

  CASM::xtal::Lattice lattice_check = reciprocal_lattice.reciprocal();
  EXPECT_TRUE(almost_equal(lattice_check.lat_column_mat(), lattice.lat_column_mat()));
}

TEST(ExampleCrystallographyLattice, ReciprocalLattice) {

  // A lattice with identity lat_column_mat is right-handed
  CASM::xtal::Lattice lattice_1 {Eigen::Matrix3d::Identity()};
  EXPECT_TRUE(lattice_1.is_right_handed());

  // Make a left-handed matrix by inverting lat_column_mat
  CASM::xtal::Lattice lattice_2 {-lattice_1.lat_column_mat()};
  EXPECT_TRUE(!lattice_2.is_right_handed());

  // Make it right-handed (which inverts lat_column_mat if it is left-handed)
  lattice_2.make_right_handed();
  EXPECT_TRUE(lattice_2.is_right_handed());
  EXPECT_TRUE(almost_equal(lattice_2.lat_column_mat(), lattice_1.lat_column_mat()));
}

// TODO: Other useful examples
// - reduced form vs Niggli form vs canonical form
