#include "ProjectBaseTest.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/ClexBasisInfo.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ConfigCorrelations.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/ProtoFuncsPrinter_impl.hh"
#include "casm/clex/io/stream/ClexBasis_stream_io.hh"
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

/// Assert relations between different correlation calculating methods
void assert_correlations_equivalence(PrimClex const &primclex,
                                     std::string const &basis_set_name,
                                     Configuration const &configuration,
                                     Clexulator const &clexulator);

class OccClexulatorFCCTest : public test::ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  OccClexulatorFCCTest()
      : test::ProjectBaseTest(test::FCC_ternary_prim(), "OccClexulatorTest",
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

TEST_F(OccClexulatorFCCTest, UseClexulator) {
  // Zeros configuration, conventional FCC unit cell
  CASM::Configuration configuration{shared_supercell};

  // Check configuration
  EXPECT_EQ(configuration.size(), 4);

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  bool align = false;
  print_basis_functions(log(), *primclex_ptr, basis_set_name, align);

  // Check clexulator
  EXPECT_EQ(clexulator.name(), "OccClexulatorTest_Clexulator");
  EXPECT_EQ(clexulator.nlist_size(), 176);
  EXPECT_EQ(clexulator.corr_size(), 75);
  EXPECT_EQ(clexulator.neighborhood().size(), 75);

  Eigen::VectorXd corr = correlations(configuration, clexulator);
  EXPECT_EQ(corr.size(), 75);
}

TEST_F(OccClexulatorFCCTest, CorrelationsTest) {
  // Ones configuration, conventional FCC unit cell
  CASM::Configuration configuration{shared_supercell};
  configuration.set_occ(0, 1);
  configuration.set_occ(1, 1);
  configuration.set_occ(2, 1);
  configuration.set_occ(3, 1);

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  assert_correlations_equivalence(*primclex_ptr, basis_set_name, configuration,
                                  clexulator);
}

class OccClexulatorZrOTest : public test::ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  OccClexulatorZrOTest()
      : test::ProjectBaseTest(test::ZrO_prim(), "OccClexulatorZrOTest",
                              jsonParser::parse(clex_basis_specs_str())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, Eigen::Matrix3l::Identity())) {
    this->write_basis_set_data();
    this->make_clexulator();
    shared_supercell->set_primclex(primclex_ptr.get());
  }

  // conventional ZrO unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
};

TEST_F(OccClexulatorZrOTest, UseClexulator) {
  // Zeros configuration, single unit cell
  CASM::Configuration configuration{shared_supercell};

  // Check configuration
  EXPECT_EQ(configuration.size(), 4);

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  bool align = false;
  print_basis_functions(log(), *primclex_ptr, basis_set_name, align);

  // Check clexulator
  EXPECT_EQ(clexulator.name(), "OccClexulatorZrOTest_Clexulator");
  EXPECT_EQ(clexulator.nlist_size(), 53);
  EXPECT_EQ(clexulator.corr_size(), 16);
  EXPECT_EQ(clexulator.neighborhood().size(), 27);

  Eigen::VectorXd corr = correlations(configuration, clexulator);
  EXPECT_EQ(corr.size(), 16);
}

TEST_F(OccClexulatorZrOTest, CorrelationsTest) {
  // Configuration w/ 1 Va, 1 O, single unit cell
  CASM::Configuration configuration{shared_supercell};
  configuration.set_occ(2, 1);  // O
  configuration.set_occ(3, 0);  // Va

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  assert_correlations_equivalence(*primclex_ptr, basis_set_name, configuration,
                                  clexulator);
}

class LocalOccClexulatorZrOTest : public test::ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  LocalOccClexulatorZrOTest()
      : test::ProjectBaseTest(test::ZrO_prim(), "LocalOccClexulatorZrOTest",
                              jsonParser::parse(clex_basis_specs_str())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, Eigen::Matrix3l::Identity() * 6)) {
    this->write_basis_set_data();
    this->make_clexulator();
    shared_supercell->set_primclex(primclex_ptr.get());
  }

  // conventional ZrO unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
};

TEST_F(LocalOccClexulatorZrOTest, UseClexulator) {
  // Zeros configuration, 6x6x6 supercell
  CASM::Configuration configuration{shared_supercell};

  // Check configuration
  EXPECT_EQ(configuration.size(), 4 * 216);

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  bool align = false;
  print_basis_functions(log(), *primclex_ptr, basis_set_name, align);

  // Check clexulator
  EXPECT_EQ(clexulator.name(), "LocalOccClexulatorZrOTest_Clexulator");
  EXPECT_EQ(clexulator.nlist_size(), 53);
  EXPECT_EQ(clexulator.corr_size(), 33);
  EXPECT_EQ(clexulator.neighborhood().size(), 26);

  Eigen::VectorXd corr = correlations(configuration, clexulator);
  EXPECT_EQ(corr.size(), 33);
}

TEST_F(LocalOccClexulatorZrOTest, CorrelationsTest) {
  // Configuration w/ 1 O, 6x6x6 supercell
  CASM::Configuration configuration{shared_supercell};
  configuration.set_occ(2 * 216, 1);  // O

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  assert_correlations_equivalence(*primclex_ptr, basis_set_name, configuration,
                                  clexulator);
}

std::string OccClexulatorFCCTest::clex_basis_specs_str() {
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

std::string OccClexulatorZrOTest::clex_basis_specs_str() {
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
    "2" : {"max_length" : 6.01},
    "3" : {"max_length" : 4.14},
    "4" : {"max_length" : 4.14}
  }
}
}
})";
}

std::string LocalOccClexulatorZrOTest::clex_basis_specs_str() {
  return R"({
"basis_function_specs" : {
"dof_specs": {
  "occ": {
    "site_basis_functions" : "occupation"
  }
}
},
"cluster_specs": {
"method": "local_max_length",
"params": {
  "generating_group" : [ 0, 3, 4, 9, 10, 11, 12, 13, 14, 15, 21, 22 ],
  "phenomenal" : {
    "coordinate_mode" : "Integral",
    "sites" : [
      [ 2, 0, 0, 0 ],
      [ 3, 0, 0, 0 ]
    ]
  },
  "orbit_branch_specs" : {
    "1" : {"cutoff_radius": 6.01},
    "2" : {"cutoff_radius": 6.01, "max_length" : 6.01}
  }
}
}
})";
}

/// Assert relations between different correlation calculating methods
///
/// Assert relations between:
/// - correlations(configuration, clexulator)
/// - all_corr_contribution(configuration, clexulator)
/// - all_point_corr(configuration, clexulator)
void assert_correlations_equivalence(PrimClex const &primclex,
                                     std::string const &basis_set_name,
                                     Configuration const &configuration,
                                     Clexulator const &clexulator) {
  std::shared_ptr<Structure const> const &shared_prim = primclex.shared_prim();
  ClexBasisSpecs const &basis_set_specs =
      primclex.basis_set_specs(basis_set_name);

  ClexBasisInfo clex_basis_info =
      make_clex_basis_info(shared_prim, basis_set_specs, primclex.nlist());
  ASSERT_EQ(clex_basis_info.basis_function_info.size(), clexulator.corr_size());

  // jsonParser basis_json;
  // write_clex_basis(shared_prim, basis_set_specs, basis_json);
  // std::cout << "basis: \n" << basis_json << std::endl;

  Eigen::VectorXd C_corr = correlations(configuration, clexulator);

  Index volume = configuration.supercell().volume();

  Eigen::VectorXd sum;
  Eigen::VectorXd mean;

  // confirm mean over unit cells of all_corr_contributions = correlations
  Eigen::MatrixXd C_all_corr_contribution =
      all_corr_contribution(configuration, clexulator);
  ASSERT_EQ(C_all_corr_contribution.rows(), volume);
  ASSERT_EQ(C_all_corr_contribution.cols(), clexulator.corr_size());
  sum = Eigen::VectorXd::Zero(clexulator.corr_size());
  for (Index i = 0; i < C_all_corr_contribution.rows(); ++i) {
    sum += C_all_corr_contribution.row(i);
  }
  mean = sum / (1.0 * C_all_corr_contribution.rows());
  ASSERT_TRUE(almost_equal(mean, C_corr));

  // confirm:
  // - all_point_corr for null cluster is 0.0
  // - mean over unit cells of all_point_corr / cluster size = correlations
  Eigen::MatrixXd C_all_point_corr = all_point_corr(configuration, clexulator);
  ASSERT_EQ(C_all_point_corr.rows(), clex_basis_info.n_point_corr * volume);
  ASSERT_EQ(C_all_point_corr.cols(), clexulator.corr_size());
  sum = Eigen::VectorXd::Zero(clexulator.corr_size());
  for (Index i = 0; i < C_all_point_corr.rows(); ++i) {
    sum += C_all_point_corr.row(i);
  }
  mean = sum / (1.0 * C_all_point_corr.rows() / clex_basis_info.n_point_corr);
  for (Index j = 0; j < clexulator.corr_size(); ++j) {
    double cluster_size =
        clex_basis_info.basis_function_info[j].prototype.size();
    if (j == 0) {
      EXPECT_TRUE(almost_equal(mean(j), 0.));
    } else {
      EXPECT_TRUE(almost_equal(mean(j) / cluster_size, C_corr(j)));
    }
  }
}
