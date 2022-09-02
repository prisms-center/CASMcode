#include "ProjectBaseTest.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/random_alloy_correlations.hh"
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

class FCCChebychevRandomAlloyTest : public test::ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  FCCChebychevRandomAlloyTest()
      : test::ProjectBaseTest(test::FCC_ternary_prim(), "OccClexulatorTest",
                              jsonParser::parse(clex_basis_specs_str())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {
    // this->write_basis_set_data();
    // this->make_clexulator();
    shared_supercell->set_primclex(primclex_ptr.get());
  }

  // conventional FCC unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
};

TEST_F(FCCChebychevRandomAlloyTest, MakeRandomAlloyCorrTest) {
  std::shared_ptr<Structure const> const &shared_prim =
      primclex_ptr->shared_prim();
  ClexBasisSpecs const &basis_set_specs =
      primclex_ptr->basis_set_specs(basis_set_name);
  RandomAlloyCorrCalculator f{shared_prim, basis_set_specs, CASM::log()};

  std::vector<Eigen::VectorXd> sublattice_prob(1);
  sublattice_prob[0].resize(3);
  sublattice_prob[0] << 1. / 3., 1. / 3., 1. / 3.;

  Eigen::VectorXd random_alloy_corr = f(sublattice_prob);
  // std::cout << "random_alloy_corr: " << random_alloy_corr.transpose() <<
  // std::endl;
  EXPECT_TRUE(almost_equal(random_alloy_corr(0), 1.0));
  for (int i = 1; i < random_alloy_corr.size(); ++i) {
    EXPECT_TRUE(almost_equal(random_alloy_corr(i), 0.0));
  }
}

class FCCOccupationRandomAlloyTest : public test::ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  FCCOccupationRandomAlloyTest()
      : test::ProjectBaseTest(test::FCC_ternary_prim(), "OccClexulatorTest",
                              jsonParser::parse(clex_basis_specs_str())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {
    // this->write_basis_set_data();
    // this->make_clexulator();
    shared_supercell->set_primclex(primclex_ptr.get());
  }

  // conventional FCC unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
};

TEST_F(FCCOccupationRandomAlloyTest, MakeRandomAlloyCorrTest) {
  std::shared_ptr<Structure const> const &shared_prim =
      primclex_ptr->shared_prim();
  ClexBasisSpecs const &basis_set_specs =
      primclex_ptr->basis_set_specs(basis_set_name);
  RandomAlloyCorrCalculator f{shared_prim, basis_set_specs, CASM::log()};

  std::vector<Eigen::VectorXd> sublattice_prob(1);
  sublattice_prob[0].resize(3);
  sublattice_prob[0] << 1., 0., 0.;

  Eigen::VectorXd random_alloy_corr = f(sublattice_prob);
  // std::cout << "random_alloy_corr: " << random_alloy_corr.transpose()
  //           << std::endl;
  EXPECT_TRUE(almost_equal(random_alloy_corr(0), 1.0));
  for (int i = 1; i < random_alloy_corr.size(); ++i) {
    EXPECT_TRUE(almost_equal(random_alloy_corr(i), 0.0));
  }
}

std::string FCCChebychevRandomAlloyTest::clex_basis_specs_str() {
  return R"({
"basis_function_specs" : {
"dof_specs": {
  "occ": {
    "site_basis_functions" : "chebychev"
  }
}
},
"cluster_specs": {
"method": "periodic_max_length",
"params": {
  "orbit_branch_specs" : {
    "2" : {"max_length" : 4.01},
    "3" : {"max_length" : 3.01}
  },
  "orbit_specs" : [
    {
      "coordinate_mode" : "Direct",
      "prototype" : [
        [ 0.000000000000, 0.000000000000, 0.000000000000 ],
        [ 1.000000000000, 0.000000000000, 0.000000000000 ],
        [ 2.000000000000, 0.000000000000, 0.000000000000 ],
        [ 3.000000000000, 0.000000000000, 0.000000000000 ]
      ],
      "include_subclusters" : true
    },
    {
      "coordinate_mode" : "Direct",
      "prototype" : [
        [ 0.000000000000, 0.000000000000, 0.000000000000 ],
        [ 0.000000000000, 1.000000000000, 0.000000000000 ],
        [ 0.000000000000, 0.000000000000, 1.000000000000 ],
        [ 1.000000000000, 1.000000000000, 1.000000000000 ]
      ],
      "include_subclusters" : true
    }
  ]
}
}
})";
}

std::string FCCOccupationRandomAlloyTest::clex_basis_specs_str() {
  return R"({
"basis_function_specs" : {
"dof_specs": {
  "occ": {
    "site_basis_functions" : "occupation"
  }
}
},
"cluster_specs": {
"method": "periodic_max_length",
"params": {
  "orbit_branch_specs" : {
    "2" : {"max_length" : 4.01},
    "3" : {"max_length" : 3.01}
  },
  "orbit_specs" : [
    {
      "coordinate_mode" : "Direct",
      "prototype" : [
        [ 0.000000000000, 0.000000000000, 0.000000000000 ],
        [ 1.000000000000, 0.000000000000, 0.000000000000 ],
        [ 2.000000000000, 0.000000000000, 0.000000000000 ],
        [ 3.000000000000, 0.000000000000, 0.000000000000 ]
      ],
      "include_subclusters" : true
    },
    {
      "coordinate_mode" : "Direct",
      "prototype" : [
        [ 0.000000000000, 0.000000000000, 0.000000000000 ],
        [ 0.000000000000, 1.000000000000, 0.000000000000 ],
        [ 0.000000000000, 0.000000000000, 1.000000000000 ],
        [ 1.000000000000, 1.000000000000, 1.000000000000 ]
      ],
      "include_subclusters" : true
    }
  ]
}
}
})";
}
