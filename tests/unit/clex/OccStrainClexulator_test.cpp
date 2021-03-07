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

namespace {
Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}
}  // namespace

class OccStrainClexulatorTest : public test::ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  OccStrainClexulatorTest()
      : test::ProjectBaseTest(test::FCC_ternary_GLstrain_prim(),
                              "OccStrainClexulatorTest",
                              jsonParser::parse(clex_basis_specs_str())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {
    this->write_basis_set_data();
    this->make_clexulator();
    shared_supercell->set_primclex(primclex_ptr.get());
  }

  // conventional FCC unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
};

std::string OccStrainClexulatorTest::clex_basis_specs_str() {
  return R"({
"basis_function_specs" : {
"global_max_poly_order": 3,
"dof_specs": {
  "occ": {
    "site_basis_functions" : "occupation"
  }
}
},
"cluster_specs": {
"method": "periodic_max_length",
"params": {
  "orbit_branch_specs": {
    "2" : {"max_length" : 2.83}
  }
}
}
})";
}

TEST_F(OccStrainClexulatorTest, UseClexulator) {
  // Zeros configuration, conventional FCC unit cell
  CASM::Configuration configuration{shared_supercell};

  // Check configuration
  EXPECT_EQ(configuration.size(), 4);

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  bool align = false;
  print_basis_functions(log(), *primclex_ptr, basis_set_name, align);

  // Check clexulator
  EXPECT_EQ(clexulator.name(), "OccStrainClexulatorTest_Clexulator");
  EXPECT_EQ(clexulator.nlist_size(), 13);
  EXPECT_EQ(clexulator.corr_size(), 33);
  EXPECT_EQ(clexulator.neighborhood().size(), 13);

  Eigen::VectorXd corr = correlations(configuration, clexulator);
  EXPECT_EQ(corr.size(), 33);
}
