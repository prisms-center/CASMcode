#include "casm/crystallography/IntegralCoordinateWithin.hh"

#include <stdexcept>
#include <vector>

#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/UnitCellCoord.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/eigen.hh"
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

xtal::Lattice hcp_lattice();

xtal::IntegralCoordinateWithin_f::matrix_type transformation_matrix() {
  Eigen::Matrix3l transformation_matrix;
  transformation_matrix << 1, 0, 3, 1, 1, -2, 1, 2, 0;
  return transformation_matrix;
}

}  // namespace

TEST(LatticePointWithinTest, construct_via_int_transformation) {
  auto trans_mat = transformation_matrix();
  xtal::IntegralCoordinateWithin_f bring_within(trans_mat);
}

TEST(LatticePointWithinTest, construct_via_transformation) {
  Eigen::Matrix3l trans_mat;
  trans_mat << 1, 0, 3, 1, 1, -2, 1, 2, 0;

  xtal::IntegralCoordinateWithin_f bring_within(trans_mat);
}

TEST(LatticePointWithinTest, construct_via_bad_transformation) {
  Eigen::Matrix3l trans_mat;
  trans_mat << 0, 0, 0, 1, 1, -2, 1, 2, 0;

  bool good_catch = false;
  try {
    xtal::IntegralCoordinateWithin_f bring_within(trans_mat);
  } catch (const std::runtime_error &e) {
    good_catch = true;
  }

  EXPECT_TRUE(good_catch) << "Determinant is " << trans_mat.determinant();
}

TEST(LatticePointWithinTest, lattice_point_within_doest_change) {
  auto trans_mat = transformation_matrix();
  xtal::IntegralCoordinateWithin_f bring_within(trans_mat);

  std::vector<Eigen::Vector3l> all_sites_within;
  all_sites_within.emplace_back(0, 0, 0);
  all_sites_within.emplace_back(1, 1, 2);
  all_sites_within.emplace_back(1, 0, 1);
  all_sites_within.emplace_back(2, 0, 1);
  all_sites_within.emplace_back(2, 0, 2);
  all_sites_within.emplace_back(3, 0, 2);
  all_sites_within.emplace_back(3, -1, 1);

  for (const auto &site : all_sites_within) {
    auto brought_within = bring_within(site);
    EXPECT_EQ(site, brought_within)
        << site.transpose() << "    vs    " << brought_within.transpose();
  }

  // Make sure it works if you do it with UnitCell instead of raw vectors
  for (const auto &raw_site : all_sites_within) {
    xtal::UnitCell site(raw_site);
    xtal::UnitCell brought_within = bring_within(site);
    EXPECT_EQ(site, brought_within)
        << site.transpose() << "    vs    " << brought_within.transpose();
  }
}

TEST(LatticePointWithinTest, bring_within_consistent_with_coordinate) {
  auto trans_mat = transformation_matrix();
  auto fcc_lat = ::fcc_lattice();
  auto fcc_superlattice = xtal::make_superlattice(fcc_lat, trans_mat);

  xtal::IntegralCoordinateWithin_f bring_within(trans_mat);

  std::vector<Eigen::Vector3l> sites_outside;
  sites_outside.emplace_back(35, 28, -84);
  sites_outside.emplace_back(-4, 86, 0);
  sites_outside.emplace_back(6, -1000, 7);

  for (const auto &site : sites_outside) {
    xtal::Coordinate coord(site(0), site(1), site(2), fcc_lat, FRAC);
    coord.set_lattice(fcc_superlattice, CART);
    coord.within();
    coord.set_lattice(fcc_lat, CART);

    auto slow_value = round(coord.const_frac());
    auto fast_value = bring_within(site);

    EXPECT_EQ(fast_value, slow_value.cast<long>());
  }
}
