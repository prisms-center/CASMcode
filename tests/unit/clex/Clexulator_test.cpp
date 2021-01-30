#include "casm/clex/Clexulator.hh"

#include "ProjectBaseTest.hh"
#include "casm/clex/Supercell.hh"
#include "gtest/gtest.h"

using namespace CASM;

//
namespace {
Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}
}  // namespace

/// Note: this builds on FCCTernaryProjectClexBasisTest from PrimClex_test.cpp
class OccupationClexulatorTest : protected ProjectBaseTest,
                                 protected testing::Test {
 protected:
  OccupationClexulatorTest()
      : ProjectBaseTest(test::FCC_ternary_prim(), "OccupationClexulatorTest",
                        test::FCC_ternary_clex_basis_specs_str_ex0()),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {}

  // conventional FCC unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
};

TEST_F(OccupationClexulatorTest, UseClexulator) {
  // Zeros configuration, conventional FCC unit cell
  CASM::Configuration configuration{shared_supercell};

  // Check clexulator
  EXPECT_EQ(configuration.size(), 4);
}
