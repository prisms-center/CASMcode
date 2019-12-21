#include "casm/crystallography/SymRepBuilder.hh"
#include "casm/external/Eigen/src/unsupported/KroneckerTensorProduct.h"
#include "casm/global/eigen.hh"
#include "gtest/gtest.h"

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
    Eigen::MatrixXd expected(30,30);
    Eigen::kroneckerProduct(flip, Eigen::MatrixXd::Identity(15, 15), expected);
    EXPECT_TRUE(result.isApprox(expected, TOL));
}
