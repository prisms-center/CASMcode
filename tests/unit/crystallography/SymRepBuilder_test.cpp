#include "casm/crystallography/SymRepBuilder.hh"

#include <ostream>

#include "casm/global/eigen.hh"
#include "casm/misc/KroneckerTensorProduct.h"
#include "gtest/gtest.h"

using namespace CASM;

TEST(SymRepBuilderTest, TimeReversalSwapSymRepBuilderConstruction) {
  TimeReversalSwapSymRepBuilder builder = TimeReversalSwapSymRepBuilder();
}

TEST(SymRepBuilderTest, TimeReversalSwapSymRepBuilderTimeReversal) {
  TimeReversalSwapSymRepBuilder builder = TimeReversalSwapSymRepBuilder();
  Eigen::MatrixXd result1 =
      builder.symop_to_matrix(Eigen::Matrix3d(), Eigen::Vector3d(), false, 2);
  Eigen::MatrixXd result2 =
      builder.symop_to_matrix(Eigen::Matrix3d(), Eigen::Vector3d(), true, 2);
  Eigen::MatrixXd exchange(2, 2);
  exchange << 0, 1, 1, 0;
  EXPECT_TRUE(result1.isApprox(Eigen::Matrix2d::Identity(), TOL));
  EXPECT_TRUE(result2.isApprox(exchange, TOL));
}

TEST(SymRepBuilderTest, dOrbitalOccupationConstruction) {
  dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
}

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderIdentity) {
  dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
  Eigen::MatrixXd result = builder.symop_to_matrix(
      Eigen::Matrix3d::Identity(), Eigen::Vector3d(), false, 15);
  EXPECT_TRUE(result.isApprox(Eigen::MatrixXd::Identity(15, 15), TOL));
}

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderInversion) {
  dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
  Eigen::MatrixXd inversion = -1 * Eigen::Matrix3d::Identity();
  Eigen::MatrixXd result =
      builder.symop_to_matrix(inversion, Eigen::Vector3d(), false, 15);
  EXPECT_TRUE(result.isApprox(Eigen::MatrixXd::Identity(15, 15), TOL));
}

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderFourfoldRotation) {
  dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
  Eigen::VectorXd d_xy = Eigen::VectorXd::Zero(15);
  d_xy(0) = 1;
  Eigen::VectorXd d_yz = Eigen::VectorXd::Zero(15);
  d_yz(1) = 1;
  Eigen::VectorXd d_z2 = Eigen::VectorXd::Zero(15);
  d_z2(2) = 1;
  Eigen::VectorXd d_xz = Eigen::VectorXd::Zero(15);
  d_xz(3) = 1;
  Eigen::VectorXd d_x2_y2 = Eigen::VectorXd::Zero(15);
  d_x2_y2(4) = 1;

  Eigen::Matrix3d C4_z;
  C4_z << 0, -1, 0, 1, 0, 0, 0, 0, 1;
  Eigen::MatrixXd symrep_C4_z =
      builder.symop_to_matrix(C4_z, Eigen::Vector3d(), false, 15);
  EXPECT_TRUE(d_xy.isApprox(symrep_C4_z * d_xy, TOL));
  EXPECT_TRUE(d_xz.isApprox(symrep_C4_z * d_yz, TOL));
  EXPECT_TRUE(d_z2.isApprox(symrep_C4_z * d_z2, TOL));
  EXPECT_TRUE(d_yz.isApprox(symrep_C4_z * d_xz, TOL));
  EXPECT_TRUE(d_x2_y2.isApprox(symrep_C4_z * d_x2_y2, TOL));

  Eigen::Matrix3d C4_x;
  C4_x << 1, 0, 0, 0, 0, -1, 0, 1, 0;
  Eigen::MatrixXd symrep_C4_x =
      builder.symop_to_matrix(C4_x, Eigen::Vector3d(), false, 15);
  Eigen::VectorXd d_y2 = Eigen::VectorXd::Zero(15);
  d_y2(2) = 0.25;
  d_y2(4) = 0.75;
  d_y2(6) = sqrt(2) * -0.5 * -sqrt(3) / 2;
  Eigen::VectorXd d_x2_z2 = Eigen::VectorXd::Zero(15);
  d_x2_z2(2) = 0.75;
  d_x2_z2(4) = 0.25;
  d_x2_z2(6) = sqrt(2) * -sqrt(3) / 2 * 0.5;
  EXPECT_TRUE(d_xz.isApprox(symrep_C4_x * d_xy, TOL));
  EXPECT_TRUE(d_yz.isApprox(symrep_C4_x * d_yz, TOL));
  EXPECT_TRUE(d_y2.isApprox(symrep_C4_x * d_z2, TOL));
  EXPECT_TRUE(d_xy.isApprox(symrep_C4_x * d_xz, TOL));
  EXPECT_TRUE(d_x2_z2.isApprox(symrep_C4_x * d_x2_y2, TOL));

  Eigen::Matrix3d C4_y;
  C4_y << 0, 0, 1, 0, 1, 0, -1, 0, 0;
  Eigen::MatrixXd symrep_C4_y =
      builder.symop_to_matrix(C4_y, Eigen::Vector3d(), false, 15);
  Eigen::VectorXd d_x2 = Eigen::VectorXd::Zero(15);
  d_x2(2) = 0.25;
  d_x2(4) = 0.75;
  d_x2(6) = sqrt(2) * -0.5 * sqrt(3) / 2;
  Eigen::VectorXd d_z2_y2 = Eigen::VectorXd::Zero(15);
  d_z2_y2(2) = 0.75;
  d_z2_y2(4) = 0.25;
  d_z2_y2(6) = sqrt(2) * 0.5 * sqrt(3) / 2;
  EXPECT_TRUE(d_yz.isApprox(symrep_C4_y * d_xy, TOL));
  EXPECT_TRUE(d_xy.isApprox(symrep_C4_y * d_yz, TOL));
  EXPECT_TRUE(d_x2.isApprox(symrep_C4_y * d_z2, TOL));
  EXPECT_TRUE(d_xz.isApprox(symrep_C4_y * d_xz, TOL));
  EXPECT_TRUE(d_z2_y2.isApprox(symrep_C4_y * d_x2_y2, TOL));
}

TEST(SymRepBuilderTest,
     dOrbitalOccupationSpinPolarizedSymRepBuilderConstruction) {
  dOrbitalOccupationSpinPolarizedSymRepBuilder builder =
      dOrbitalOccupationSpinPolarizedSymRepBuilder();
}

TEST(SymRepBuilderTest, dOrbitalOccupationSpinPolarizedSymRepBuilderIdentity) {
  dOrbitalOccupationSpinPolarizedSymRepBuilder builder =
      dOrbitalOccupationSpinPolarizedSymRepBuilder();
  Eigen::MatrixXd result = builder.symop_to_matrix(
      Eigen::Matrix3d::Identity(), Eigen::Vector3d(), false, 30);
  EXPECT_TRUE(result.isApprox(Eigen::MatrixXd::Identity(30, 30), TOL));
}

TEST(SymRepBuilderTest,
     dOrbitalOccupationSpinPolarizedSymRepBuilderTimeReversal) {
  dOrbitalOccupationSpinPolarizedSymRepBuilder builder =
      dOrbitalOccupationSpinPolarizedSymRepBuilder();
  Eigen::MatrixXd result = builder.symop_to_matrix(Eigen::Matrix3d::Identity(),
                                                   Eigen::Vector3d(), true, 30);
  Eigen::MatrixXd flip(2, 2);
  flip << 0, 1, 1, 0;
  Eigen::MatrixXd expected(30, 30);
  Eigen::kroneckerProduct(flip, Eigen::MatrixXd::Identity(15, 15), expected);
  EXPECT_TRUE(result.isApprox(expected, TOL));
}
