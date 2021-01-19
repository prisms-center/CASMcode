#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/misc/CASM_Eigen_math.hh"  // for almost_equal of Eigen types
#include "gtest/gtest.h"

TEST(ExampleCrystallographyCoordinate, CoordinateConstructor) {
  // First, construct a Lattice
  Eigen::Vector3d a{1.0, 0.0, 0.0};
  Eigen::Vector3d b{0.0, 2.0, 0.0};
  Eigen::Vector3d c{0.0, 0.0, 3.0};
  CASM::xtal::Lattice lattice{a, b, c};

  // The CASM::xtal::Coordinate is always referenced to a CASM::xtal::Lattice
  // When it is constructed with coordinate values, the constructor must also be
  // told whether those values are Cartesian coordinates (using CART) or
  // fractional coordinates w.r.t the lattice vectors (using FRAC).

  // These two coordinates are equivalent:
  CASM::xtal::Coordinate coord_1{Eigen::Vector3d{1.0, 2.0, 3.0}, lattice,
                                 CASM::CART};
  CASM::xtal::Coordinate coord_2{Eigen::Vector3d{1.0, 1.0, 1.0}, lattice,
                                 CASM::FRAC};

  EXPECT_TRUE(
      almost_equal(coord_1.as_vec(CASM::CART), Eigen::Vector3d{1.0, 2.0, 3.0}));
  EXPECT_TRUE(
      almost_equal(coord_2.as_vec(CASM::CART), Eigen::Vector3d{1.0, 2.0, 3.0}));

  EXPECT_TRUE(
      almost_equal(coord_1.as_vec(CASM::FRAC), Eigen::Vector3d{1.0, 1.0, 1.0}));
  EXPECT_TRUE(
      almost_equal(coord_2.as_vec(CASM::FRAC), Eigen::Vector3d{1.0, 1.0, 1.0}));

  EXPECT_TRUE(coord_1.almost_equal(coord_2));
}
