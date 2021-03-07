#include "ProjectBaseTest.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/stream/ClexBasis_stream_io.hh"
#include "clex/TestClexBasisSpecs.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

class StrainClexulatorTest : public test::ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  StrainClexulatorTest()
      : test::ProjectBaseTest(test::SimpleCubic_GLstrain_prim(),
                              "StrainClexulatorTest",
                              jsonParser::parse(clex_basis_specs_str())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, Eigen::Matrix3l::Identity())) {
    this->write_basis_set_data();
    this->make_clexulator();
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
  "orbit_branch_specs": {}
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

  bool align = false;
  print_basis_functions(log(), *primclex_ptr, basis_set_name, align);

  // Check clexulator
  EXPECT_EQ(clexulator.name(), "StrainClexulatorTest_Clexulator");
  EXPECT_EQ(clexulator.nlist_size(), 0);
  EXPECT_EQ(clexulator.corr_size(), 11);
  EXPECT_EQ(clexulator.neighborhood().size(), 0);

  Eigen::VectorXd corr = correlations(configuration, clexulator);
  EXPECT_EQ(corr.size(), 11);
}
