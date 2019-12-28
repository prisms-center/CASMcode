#include "casm/crystallography/SymRepBuilder.hh"
#include "casm/external/Eigen/src/unsupported/KroneckerTensorProduct.h"
#include "casm/global/eigen.hh"
#include "gtest/gtest.h"
#include <ostream>

using namespace CASM;

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderConstruction)
{
    dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
}

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderDimension)
{
    dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
    Eigen::MatrixXd result = builder.symop_to_matrix(Eigen::Matrix3d::Identity(), Eigen::Vector3d(), false, 30);
    EXPECT_EQ(result.rows(), 30);
    EXPECT_EQ(result.cols(), 30);
}

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderIdentity)
{
    dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
    Eigen::MatrixXd result = builder.symop_to_matrix(Eigen::Matrix3d::Identity(), Eigen::Vector3d(), false, 30);
    EXPECT_TRUE(result.isApprox(Eigen::MatrixXd::Identity(30, 30), TOL));
}

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderTimeReversal)
{
    dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
    Eigen::MatrixXd result = builder.symop_to_matrix(Eigen::Matrix3d::Identity(), Eigen::Vector3d(), true, 30);
    Eigen::MatrixXd flip(2, 2);
    flip << 0, 1, 1, 0;
    Eigen::MatrixXd expected(30, 30);
    Eigen::kroneckerProduct(flip, Eigen::MatrixXd::Identity(15, 15), expected);
    EXPECT_TRUE(result.isApprox(expected, TOL));
}

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderInversion)
{
    dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
    Eigen::MatrixXd inversion = -1 * Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd result = builder.symop_to_matrix(inversion, Eigen::Vector3d(), false, 30);
    EXPECT_TRUE(result.isApprox(Eigen::MatrixXd::Identity(30, 30), TOL));
}

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderOrthogonal)
{
    dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();
    Eigen::MatrixXd id30 = Eigen::MatrixXd::Identity(30, 30);
    Eigen::Matrix3d orth_1, orth_2, orth_3, orth_4, orth_5, orth_6;
    orth_1 << 0.00167153, 0.738058, 0.674735, 0.998724, 0.032826, -0.0383808, -0.050476, 0.673938, -0.737061;
    orth_2 << 0.416724, 0.586135, 0.694828, 0.627959, 0.36705, -0.686252, 0.657273, -0.722301, 0.215111;
    orth_3 << 0.962225, 0.234022, 0.139131, -0.271482, 0.786203, 0.555141, 0.0205302, -0.571942, 0.820037;
    orth_4 << 0.574223, 0.521022, 0.63151, 0.662831, -0.74862, 0.0149398, 0.480545, 0.410005, -0.775224;
    orth_5 << 0.355972, 0.900777, 0.248764, -0.327, -0.129304, 0.936136, 0.875417, -0.414585, 0.248525;
    orth_6 << 0.130944, 0.05601, 0.989806, 0.534307, 0.837007, -0.118049, -0.835087, 0.544318, 0.0796748;
    Eigen::MatrixXd symrep_1 = builder.symop_to_matrix(orth_1, Eigen::Vector3d(), false, 30); 
    Eigen::MatrixXd symrep_2 = builder.symop_to_matrix(orth_2, Eigen::Vector3d(), false, 30); 
    Eigen::MatrixXd symrep_3 = builder.symop_to_matrix(orth_3, Eigen::Vector3d(), false, 30); 
    Eigen::MatrixXd symrep_4 = builder.symop_to_matrix(orth_4, Eigen::Vector3d(), false, 30); 
    Eigen::MatrixXd symrep_5 = builder.symop_to_matrix(orth_5, Eigen::Vector3d(), false, 30); 
    Eigen::MatrixXd symrep_6 = builder.symop_to_matrix(orth_6, Eigen::Vector3d(), false, 30); 
    EXPECT_TRUE(id30.isApprox(symrep_1 * symrep_1.transpose(), TOL));
    EXPECT_TRUE(id30.isApprox(symrep_2 * symrep_2.transpose(), TOL));
    EXPECT_TRUE(id30.isApprox(symrep_3 * symrep_3.transpose(), TOL));
    EXPECT_TRUE(id30.isApprox(symrep_4 * symrep_4.transpose(), TOL));
    EXPECT_TRUE(id30.isApprox(symrep_5 * symrep_5.transpose(), TOL));
    EXPECT_TRUE(id30.isApprox(symrep_6 * symrep_6.transpose(), TOL));
}

TEST(SymRepBuilderTest, dOrbitalOccupationSymRepBuilderFourfoldRotation)
{
    dOrbitalOccupationSymRepBuilder builder = dOrbitalOccupationSymRepBuilder();    
    Eigen::VectorXd d_xy = Eigen::VectorXd::Zero(30);
    d_xy(0) = 1;
    Eigen::VectorXd d_yz = Eigen::VectorXd::Zero(30);
    d_yz(1) = 1;
    Eigen::VectorXd d_z2 = Eigen::VectorXd::Zero(30);
    d_z2(2) = 1;
    Eigen::VectorXd d_xz = Eigen::VectorXd::Zero(30);
    d_xz(3) = 1;
    Eigen::VectorXd d_x2_y2 = Eigen::VectorXd::Zero(30);
    d_x2_y2(4) = 1;

    Eigen::Matrix3d C4_z;
    C4_z << 0, -1, 0, 1, 0, 0, 0, 0, 1;
    Eigen::MatrixXd symrep_C4_z = builder.symop_to_matrix(C4_z, Eigen::Vector3d(), false, 30);
    EXPECT_TRUE(d_xy.isApprox(symrep_C4_z * d_xy, TOL));
    EXPECT_TRUE(d_xz.isApprox(symrep_C4_z * d_yz, TOL));
    EXPECT_TRUE(d_z2.isApprox(symrep_C4_z * d_z2, TOL));
    EXPECT_TRUE(d_yz.isApprox(symrep_C4_z * d_xz, TOL));
    EXPECT_TRUE(d_x2_y2.isApprox(symrep_C4_z * d_x2_y2, TOL));

    Eigen::Matrix3d C4_x;
    C4_x << 1, 0, 0, 0, 0, -1, 0, 1, 0;
    Eigen::MatrixXd symrep_C4_x = builder.symop_to_matrix(C4_x, Eigen::Vector3d(), false, 30);
    Eigen::VectorXd d_y2 = Eigen::VectorXd::Zero(30);
    d_y2(2) = 0.25;
    d_y2(4) = 0.75; 
    d_y2(6) = sqrt(2) * -0.5 * -sqrt(3)/2;
    Eigen::VectorXd d_x2_z2 = Eigen::VectorXd::Zero(30);
    d_x2_z2(2) = 0.75;
    d_x2_z2(4) = 0.25;
    d_x2_z2(6) = sqrt(2) * -sqrt(3)/2 * 0.5;
    EXPECT_TRUE(d_xz.isApprox(symrep_C4_x * d_xy, TOL));
    EXPECT_TRUE(d_yz.isApprox(symrep_C4_x * d_yz, TOL));
    EXPECT_TRUE(d_y2.isApprox(symrep_C4_x * d_z2, TOL));
    EXPECT_TRUE(d_xy.isApprox(symrep_C4_x * d_xz, TOL));
    EXPECT_TRUE(d_x2_z2.isApprox(symrep_C4_x * d_x2_y2, TOL));

    Eigen::Matrix3d C4_y;
    C4_y << 0, 0, 1, 0, 1, 0, -1, 0, 0;
    Eigen::MatrixXd symrep_C4_y = builder.symop_to_matrix(C4_y, Eigen::Vector3d(), false, 30);
    Eigen::VectorXd d_x2 = Eigen::VectorXd::Zero(30);
    d_x2(2) = 0.25;
    d_x2(4) = 0.75;
    d_x2(6) = sqrt(2) * -0.5 * sqrt(3)/2;
    Eigen::VectorXd d_z2_y2 = Eigen::VectorXd::Zero(30);
    d_z2_y2(2) = 0.75;
    d_z2_y2(4) = 0.25;
    d_z2_y2(6) = sqrt(2) * 0.5 * sqrt(3)/2;
    EXPECT_TRUE(d_yz.isApprox(symrep_C4_y * d_xy, TOL));
    EXPECT_TRUE(d_xy.isApprox(symrep_C4_y * d_yz, TOL)); 
    EXPECT_TRUE(d_x2.isApprox(symrep_C4_y * d_z2, TOL));
    EXPECT_TRUE(d_xz.isApprox(symrep_C4_y * d_xz, TOL));
    EXPECT_TRUE(d_z2_y2.isApprox(symrep_C4_y * d_x2_y2, TOL));
}

