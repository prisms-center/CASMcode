#include "gtest/gtest.h"
#include <algorithm>
#include <cmath>
#include <vector>

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/external/Eigen/src/Core/Matrix.h"
#include "casm/global/definitions.hh"

using namespace CASM;

namespace {
  Eigen::Matrix3d flat_rotation_matrix(double angle) {
    const auto PI = std::acos(-1);
    Eigen::Matrix3d rotate_matrix;
    double radians = angle * 2 * PI / 360;
    rotate_matrix << std::cos(radians), -std::sin(angle), 0, std::sin(angle), std::cos(angle), 0, 0, 0, 1;

    return rotate_matrix;
  }

  xtal::Lattice hcp_lattice() {
    Eigen::Matrix3d hcp_matrix;
    hcp_matrix << 3.1837728024, -1.5918864012, 0.0000000000,
               0.0000000000, 2.7572281267, 0.0000000000,
               0.0000000000, 0.0000000000, 5.2494397163;

    return xtal::Lattice(hcp_matrix);
  }

} // namespace

TEST(SymOpTest, construct) {
  Eigen::Matrix3d m;
  m << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  Eigen::Vector3d t(11, 22, 33);

  xtal::SymOp random_op_true(m, t, true);
  xtal::SymOp random_op_false(m, t, false);
  EXPECT_EQ(xtal::get_matrix(random_op_true), m);
  EXPECT_EQ(xtal::get_translation(random_op_true), t);
  EXPECT_TRUE(xtal::get_time_reversal(random_op_true));
  EXPECT_FALSE(xtal::get_time_reversal(random_op_false));

  xtal::SymOp identity_op = xtal::SymOp::identity();
  EXPECT_EQ(xtal::get_matrix(identity_op), Eigen::Matrix3d::Identity());
  EXPECT_EQ(xtal::get_translation(identity_op), Eigen::Vector3d::Zero());
  EXPECT_FALSE(xtal::get_time_reversal(identity_op));

  xtal::SymOp time_reversal_op = xtal::SymOp::time_reversal();
  EXPECT_EQ(xtal::get_matrix(time_reversal_op), Eigen::Matrix3d::Identity());
  EXPECT_EQ(xtal::get_translation(time_reversal_op), Eigen::Vector3d::Zero());
  EXPECT_TRUE(xtal::get_time_reversal(time_reversal_op));

  xtal::SymOp translation_op = xtal::SymOp::translation_operation(t);
  EXPECT_EQ(xtal::get_matrix(translation_op), t);
  EXPECT_EQ(xtal::get_translation(translation_op), Eigen::Vector3d::Zero());
  EXPECT_FALSE(xtal::get_time_reversal(translation_op));

  xtal::SymOp point_op = xtal::SymOp::point_operation(m);
  EXPECT_EQ(xtal::get_matrix(point_op), m);
  EXPECT_EQ(xtal::get_translation(point_op), Eigen::Vector3d::Zero());
  EXPECT_FALSE(xtal::get_time_reversal(point_op));
}

TEST(SymOpTest, multiplication) {
  Eigen::Matrix3d translate_up_matrix(1.1, 2.2, 3.3);
  Eigen::Matrix3d translate_down_matrix = -translate_up_matrix;

  xtal::SymOp translate_up = xtal::SymOp::translation_operation(translate_up_matrix);
  xtal::SymOp translate_down = xtal::SymOp::translation_operation(translate_down_matrix);
  xtal::SymOp net_zero_translation = translate_up * translate_down;

  xtal::SymOpCompare_f is_identity(xtal::SymOp::identity(), CASM::TOL);
  EXPECT_TRUE(is_identity(net_zero_translation));

  xtal::SymOp rotate_60 = xtal::SymOp::point_operation(::flat_rotation_matrix(60));
  xtal::SymOp rotate_360 = rotate_60 * rotate_60 * rotate_60 * rotate_60 * rotate_60 * rotate_60;
  EXPECT_TRUE(is_identity(rotate_360));
}

TEST(SymOpTest, point_group_closure) {
  xtal::SymOp rotate_72 = xtal::SymOp::point_operation(::flat_rotation_matrix(72));
  xtal::SymOp identity = xtal::SymOp::identity();

  std::vector<xtal::SymOp> rotation5_group{identity, rotate_72};
  xtal::close_group<xtal::SymOpCompare_f>(&rotation5_group, TOL);

  EXPECT_EQ(rotation5_group.size(), 5);
}

TEST(SymOpTest, periodic_group_closure) {
  xtal::Lattice hcp_lattice =::hcp_lattice();
  Eigen::Vector3d translate_6th = hcp_lattice[2] * (1.0 / 6);

  xtal::SymOp screw_60 = xtal::SymOp(::flat_rotation_matrix(60), translate_6th, false);
  std::vector<xtal::SymOp> screw_60_group{screw_60};

  xtal::close_group<xtal::SymOpPeriodicCompare_f>(&screw_60_group, hcp_lattice, TOL);
  EXPECT_EQ(screw_60_group.size(), 6);
}

//TODO: Break this into separate tests?
TEST(SymOpTest, comparisons) {
  xtal::Lattice hcp_lattice =::hcp_lattice();
  Eigen::Vector3d translate_6th = hcp_lattice[2] * (1.0 / 6);

  xtal::SymOp identity = xtal::SymOp::identity();
  xtal::SymOp time_reversal = xtal::SymOp::time_reversal();
  xtal::SymOp rotate_60 = xtal::SymOp::point_operation(::flat_rotation_matrix(60));
  xtal::SymOp screw_60 = xtal::SymOp(rotate_60.matrix, translate_6th, false);

  xtal::SymOp rotate_72 = xtal::SymOp::point_operation(::flat_rotation_matrix(72));
  xtal::SymOp screw_60_equivalent = xtal::SymOp(rotate_60.matrix, 3 * hcp_lattice[2] + translate_6th, false);

  std::vector<xtal::SymOp> fake_group{identity, time_reversal, rotate_60, screw_60};

  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpCompare_f(identity, TOL)), 1);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpCompare_f(time_reversal, TOL)), 1);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpCompare_f(rotate_60, TOL)), 1);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpCompare_f(screw_60, TOL)), 1);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpCompare_f(rotate_72, TOL)), 0);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpCompare_f(screw_60_equivalent, TOL)), 0);

  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpPeriodicCompare_f(identity, hcp_lattice, TOL)), 1);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpPeriodicCompare_f(time_reversal, hcp_lattice, TOL)), 1);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpPeriodicCompare_f(rotate_60, hcp_lattice, TOL)), 1);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpPeriodicCompare_f(screw_60, hcp_lattice, TOL)), 2);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpPeriodicCompare_f(rotate_72, hcp_lattice, TOL)), 0);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpPeriodicCompare_f(screw_60_equivalent, hcp_lattice, TOL)), 2);

  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpMatrixCompare_f(identity, TOL)), 2);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpMatrixCompare_f(time_reversal, TOL)), 2);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpMatrixCompare_f(rotate_60, TOL)), 3);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpMatrixCompare_f(screw_60, TOL)), 3);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpMatrixCompare_f(rotate_72, TOL)), 0);
  EXPECT_EQ(std::count_if(fake_group.begin(), fake_group.end(), xtal::SymOpMatrixCompare_f(screw_60_equivalent, TOL)), 3);
}
