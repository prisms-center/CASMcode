#include "gtest/gtest.h"

/// What is being tested:
#include "casm/container/Counter.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/misc/CASM_Eigen_math.hh"

using namespace CASM;
using namespace xtal;

void coordinate_constructor_test() {
  double tol = 1e-5;
  Lattice lat(Lattice::hexagonal());
  Coordinate::vector_type tvec(0.1, 0.2, 0.3);

  {
    // Make sure that new Coordinate gets correct default values
    // could introduce bugs if this changes
    Coordinate tcoord(lat);
    EXPECT_TRUE(almost_zero(tcoord.const_frac(), tol));
    EXPECT_TRUE(almost_zero(tcoord.const_cart(), tol));
    EXPECT_EQ(&(tcoord.home()), &lat);
  }

  {
    Coordinate tcoord(tvec, lat, CART);
    EXPECT_TRUE(almost_equal(tcoord.const_cart(), tvec, tol));
    EXPECT_TRUE(almost_equal(tcoord.const_frac(),
                             lat.inv_lat_column_mat() * tvec, tol));
    EXPECT_EQ(&(tcoord.home()), &lat);
  }

  {
    Coordinate tcoord(tvec, lat, FRAC);
    EXPECT_TRUE(almost_equal(tcoord.const_frac(), tvec, tol));
    EXPECT_TRUE(
        almost_equal(tcoord.const_cart(), lat.lat_column_mat() * tvec, tol));
    EXPECT_EQ(&(tcoord.home()), &lat);
  }

  {
    Coordinate tcoord(0.1, 0.2, 0.3, lat, CART);
    EXPECT_TRUE(almost_equal(tcoord.const_cart(), tvec, tol));
    EXPECT_TRUE(almost_equal(tcoord.const_frac(),
                             lat.inv_lat_column_mat() * tvec, tol));
    EXPECT_EQ(&(tcoord.home()), &lat);
  }

  {
    Coordinate tcoord(0.1, 0.2, 0.3, lat, FRAC);
    EXPECT_TRUE(almost_equal(tcoord.const_frac(), tvec, tol));
    EXPECT_TRUE(
        almost_equal(tcoord.const_cart(), lat.lat_column_mat() * tvec, tol));
    EXPECT_EQ(&(tcoord.home()), &lat);
  }
}

void coordinate_operation_test() {
  double tol = 1e-5;
  Lattice lat(Lattice::hexagonal());

  Coordinate::vector_type vec1(0.1, 0.2, 0.3);
  Coordinate::vector_type vec2(0.6, 0.4, 0.2);
  /// non-const frac access
  {
    Coordinate coord1(lat);
    coord1.frac() = vec2;
    EXPECT_TRUE(almost_equal(coord1.const_frac(), vec2, tol));
    EXPECT_TRUE(
        almost_equal(coord1.const_cart(), lat.lat_column_mat() * vec2, tol));
  }

  /// non-const frac access
  {
    Coordinate coord1(lat);
    coord1.frac(0) = vec2[0];
    coord1.frac(1) = vec2[1];
    coord1.frac(2) = vec2[2];
    EXPECT_TRUE(almost_equal(coord1.const_frac(), vec2, tol));
    EXPECT_TRUE(
        almost_equal(coord1.const_cart(), lat.lat_column_mat() * vec2, tol));
  }

  /// non-const cart access
  {
    Coordinate coord1(lat);
    coord1.cart() = vec2;

    EXPECT_TRUE(almost_equal(coord1.const_cart(), vec2, tol));
    EXPECT_TRUE(almost_equal(coord1.const_frac(),
                             lat.inv_lat_column_mat() * vec2, tol));
  }

  /// non-const cart access
  {
    Coordinate coord1(lat);
    coord1.cart(0) = vec2[0];
    coord1.cart(1) = vec2[1];
    coord1.cart(2) = vec2[2];

    EXPECT_TRUE(almost_equal(coord1.const_cart(), vec2, tol));
    EXPECT_TRUE(almost_equal(coord1.const_frac(),
                             lat.inv_lat_column_mat() * vec2, tol));
  }

  Coordinate coord1(vec1, lat, CART);
  Coordinate coord2(vec2, lat, CART);

  /// frac assignment
  {
    Coordinate tcoord(coord1);
    tcoord.frac() = coord2.frac();
    EXPECT_TRUE(almost_equal(tcoord.const_cart(), coord2.const_cart(), tol));
    EXPECT_TRUE(almost_equal(tcoord.const_frac(), coord2.const_frac(), tol));
  }

  /// frac assignment
  {
    Coordinate tcoord(coord1);
    tcoord.frac(0) = coord2.frac(0);
    tcoord.frac(1) = coord2.frac(1);
    tcoord.frac(2) = coord2.frac(2);
    EXPECT_TRUE(almost_equal(tcoord.const_cart(), coord2.const_cart(), tol));
    EXPECT_TRUE(almost_equal(tcoord.const_frac(), coord2.const_frac(), tol));
  }

  /// cart assignment
  {
    Coordinate tcoord(coord1);
    tcoord.cart(0) = coord2.cart(0);
    tcoord.cart(1) = coord2.cart(1);
    tcoord.cart(2) = coord2.cart(2);

    EXPECT_TRUE(almost_equal(tcoord.const_cart(), coord2.const_cart(), tol));
    EXPECT_TRUE(almost_equal(tcoord.const_frac(), coord2.const_frac(), tol));
  }

  /// unary minus
  {
    Coordinate tcoord(lat);
    tcoord = -coord2;
    EXPECT_TRUE(almost_equal(tcoord.const_cart(), -coord2.const_cart(), tol));
    EXPECT_TRUE(almost_equal(tcoord.const_frac(), -coord2.const_frac(), tol));
  }

  /// plus-assignment
  {
    Coordinate tcoord(coord1);
    tcoord += coord2;
    EXPECT_TRUE(almost_equal(tcoord.const_cart(), vec1 + vec2, tol));
    EXPECT_TRUE(almost_equal(tcoord.const_frac(),
                             lat.inv_lat_column_mat() * (vec1 + vec2), tol));
  }

  /// subtract-assignment
  {
    Coordinate tcoord(coord1);
    tcoord -= coord2;
    EXPECT_TRUE(almost_equal(tcoord.const_cart(), vec1 - vec2, tol));
    EXPECT_TRUE(almost_equal(tcoord.const_frac(),
                             lat.inv_lat_column_mat() * (vec1 - vec2), tol));
  }

  /// equality
  {
    Coordinate tcoord(coord1);
    EXPECT_TRUE(tcoord == coord1);
    EXPECT_TRUE(!(tcoord == coord2));

    EXPECT_TRUE(tcoord != coord2);
    EXPECT_TRUE(!(tcoord != coord1));
  }
}

void coordinate_periodicity_test() {
  double tol = 1e-5;
  Lattice lat(Lattice::hexagonal());

  /// min_dist
  {
    Coordinate::vector_type vec1(0.1, 0.2, 0.3);
    Coordinate::vector_type vec2(0.5, 0.35, 0.2);
    Coordinate coord1(vec1, lat, FRAC);
    Coordinate coord2(vec2, lat, FRAC);
    // 27 translations
    EigenCounter<Eigen::Vector3i> trans_count(Eigen::Vector3i(-2, -2, -2),
                                              Eigen::Vector3i(2, 2, 2),
                                              Eigen::Vector3i(2, 2, 2));
    Coordinate transcoord1(lat);
    Coordinate nearest_trans(lat);
    for (; trans_count.valid(); ++trans_count) {
      transcoord1.frac() = trans_count().cast<double>() + coord1.const_frac();

      // commutativity of min_dist from transcoord1 to coord2
      EXPECT_TRUE(almost_equal(transcoord1.min_dist(coord2),
                               coord2.min_dist(transcoord1), tol));

      // minimality of min_dist
      EXPECT_TRUE((transcoord1.min_dist(coord2) -
                   (coord1.const_cart() - coord2.const_cart()).norm()) < tol);

      // min_dist transcoord1 to coord1
      EXPECT_TRUE(almost_zero(transcoord1.min_dist(coord1), tol));
      EXPECT_TRUE(almost_zero(coord1.min_dist(transcoord1), tol));

      // min_dist transcoord1 to coord1 w/ translation check
      nearest_trans = transcoord1.min_translation(coord1);
      EXPECT_TRUE(almost_zero(transcoord1.min_dist(coord1), tol));
      EXPECT_TRUE(almost_zero(nearest_trans.const_frac(), tol));

      // min_dist coord1 to transcoord1 w/ translation check
      nearest_trans = coord1.min_translation(transcoord1);
      EXPECT_TRUE(almost_zero(coord1.min_dist(transcoord1), tol));
      EXPECT_TRUE(almost_zero(nearest_trans.const_frac(), tol));
    }
  }

  /// voronoi min_dist
  {
    Coordinate::vector_type vec1(1.1, 2.1, 3.1);
    Coordinate coord1(vec1, lat, FRAC);

    // 180 points inside the primitive cell
    EigenCounter<Eigen::Vector3d> point_count(Eigen::Vector3d(0., 0., 0.),
                                              Eigen::Vector3d(1., 1., 1.),
                                              Eigen::Vector3d(0.2, 0.2, 0.25));
    for (; point_count.valid(); ++point_count) {
      Coordinate coord2(point_count(), lat, FRAC);
      // optimality of robust_min_dist
      EXPECT_TRUE(coord1.robust_min_dist(coord2) <
                  (coord1.min_dist(coord2) + tol));

      // commutativity of robust_min_dist
      EXPECT_TRUE(almost_equal(coord2.robust_min_dist(coord1),
                               coord1.robust_min_dist(coord2), tol));

      // equivalence of both versions of robust_min_dist
      EXPECT_TRUE(almost_equal(coord2.robust_min_dist(coord1),
                               coord1.robust_min_dist(coord2), tol));
      EXPECT_TRUE(almost_equal(coord2.robust_min_dist(coord1),
                               coord1.robust_min_dist(coord2), tol));
    }

    Eigen::MatrixXd const &vtable = lat.voronoi_table();
    Coordinate zero_coord(Eigen::Vector3d::Zero(), lat, FRAC);

    // check robustness of edge cases
    for (Index i = 0; i < vtable.rows(); i++) {
      {
        Coordinate coord2((1 + tol / 2.) * vtable.row(i).transpose() /
                              vtable.row(i).squaredNorm(),
                          lat, CART);
        EXPECT_TRUE(coord2.robust_min_dist(zero_coord) <=
                    coord2.const_cart().norm());
      }
      {
        Coordinate coord2((1 - tol / 2.) * vtable.row(i).transpose() /
                              vtable.row(i).squaredNorm(),
                          lat, CART);
        EXPECT_TRUE(almost_equal(coord2.robust_min_dist(zero_coord),
                                 coord2.const_cart().norm(), tol / 2));
      }
    }
  }
}

TEST(CoordinateTest, ConstructorTest) { coordinate_constructor_test(); }

TEST(CoordinateTest, OperationTest) { coordinate_operation_test(); }

TEST(CoordinateTest, PeriodicityTest) { coordinate_periodicity_test(); }
