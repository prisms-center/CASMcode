#include "casm/clex/Clexulator.hh"

#include "ProjectBaseTest.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/Supercell.hh"
#include "clex/TestClexBasisSpecs.hh"
#include "crystallography/TestStructures.hh"
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
class OccupationClexulatorTest : public test::ProjectBaseTest {
 protected:
  OccupationClexulatorTest()
      : test::ProjectBaseTest(
            test::FCC_ternary_prim(), "OccupationClexulatorTest",
            jsonParser::parse(test::FCC_ternary_clex_basis_specs_str_ex0())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {
    shared_supercell->set_primclex(primclex_ptr.get());
  }

  // conventional FCC unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
};

TEST_F(OccupationClexulatorTest, UseClexulator) {
  // Zeros configuration, conventional FCC unit cell
  CASM::Configuration configuration{shared_supercell};

  // Check configuration
  EXPECT_EQ(configuration.size(), 4);

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  // Check clexulator
  EXPECT_EQ(clexulator.name(), "OccupationClexulatorTest_Clexulator");
  EXPECT_EQ(clexulator.nlist_size(), 176);
  EXPECT_EQ(clexulator.corr_size(), 75);
  EXPECT_EQ(clexulator.neighborhood().size(), 75);

  Eigen::VectorXd corr = correlations(configuration, clexulator);
  EXPECT_EQ(corr.size(), 75);
}
