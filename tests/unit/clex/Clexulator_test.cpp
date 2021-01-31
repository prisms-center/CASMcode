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

/// Note: this builds on FCCTernaryProjectClexBasisTest from PrimClex_test.cpp
class StrainClexulatorTest : public test::ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  StrainClexulatorTest()
      : test::ProjectBaseTest(test::SimpleCubic_GLstrain_prim(),
                              "StrainClexulatorTest",
                              jsonParser::parse(clex_basis_specs_str())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, Eigen::Matrix3l::Identity())) {
    shared_supercell->set_primclex(primclex_ptr.get());
  }

  // conventional FCC unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
};

std::string StrainClexulatorTest::clex_basis_specs_str() {
  return R"({
"basis_function_specs" : {
"global_max_poly_order": 3
},
"cluster_specs": {
"method": "periodic_max_length",
"params": {
  "orbit_branch_specs": {
  }
}
}
})";
}

TEST_F(StrainClexulatorTest, UseClexulator) {
  // Zeros configuration, simple cubic unit cell
  CASM::Configuration configuration{shared_supercell};

  // Check configuration
  EXPECT_EQ(configuration.size(), 1);

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  // Check clexulator
  EXPECT_EQ(clexulator.name(), "StrainClexulatorTest_Clexulator");
  EXPECT_EQ(clexulator.nlist_size(), 176);
  EXPECT_EQ(clexulator.corr_size(), 75);
  EXPECT_EQ(clexulator.neighborhood().size(), 75);

  Eigen::VectorXd corr = correlations(configuration, clexulator);
  EXPECT_EQ(corr.size(), 75);
}
