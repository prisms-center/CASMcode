#include "ProjectBaseTest.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/ClexBasisFunctionInfo_impl.hh"
#include "casm/clex/ClexBasis_impl.hh"
#include "casm/clex/Clexulator.hh"
#include "casm/clex/ConfigCorrelations.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/NeighborhoodInfo_impl.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/ProtoFuncsPrinter_impl.hh"
#include "casm/clex/io/stream/ClexBasis_stream_io.hh"
#include "casm/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

namespace {
Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}

struct MakeClexBasisInfo {
  MakeClexBasisInfo(PrimNeighborList &_prim_neighbor_list,
                    std::vector<ClexBasisFunctionInfo> &_basis_function_info,
                    std::unique_ptr<NeighborhoodInfo> &_neighborhood_info)
      : prim_neighbor_list(_prim_neighbor_list),
        basis_function_info(_basis_function_info),
        neighborhood_info(_neighborhood_info) {}

  PrimNeighborList &prim_neighbor_list;
  std::vector<ClexBasisFunctionInfo> &basis_function_info;
  std::unique_ptr<NeighborhoodInfo> &neighborhood_info;

  template <typename OrbitVecType>
  void operator()(ClexBasis const &clex_basis,
                  OrbitVecType const &orbits) const {
    neighborhood_info = notstd::make_unique<NeighborhoodInfo>(
        prim_neighbor_list, orbits.begin(), orbits.end());
    basis_function_info =
        make_clex_basis_function_info(clex_basis, orbits.begin(), orbits.end());
  }
};

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
  EXPECT_EQ(clexulator.name(), "OccClexulatorTest_Clexulator_default");
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
  EXPECT_EQ(clexulator.name(), "OccClexulatorZrOTest_Clexulator_default");
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

    // Uncomment to preserve project
    // tmp_dir.do_not_remove_on_destruction();
  }

  void use_clexulator_test();
  void eval_correlations_test();
  void correlations_correctness_test();

  Eigen::VectorXd local_corr(xtal::UnitCell unitcell,
                             Configuration const &config,
                             Clexulator &clexulator) {
    auto const &scel_sym_info = config.supercell().sym_info();
    return corr_contribution(scel_sym_info.unitcell_index_converter()(unitcell),
                             config, clexulator);
  }

  // conventional ZrO unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
};

void LocalOccClexulatorZrOTest::use_clexulator_test() {
  // Zeros configuration, 6x6x6 supercell
  CASM::Configuration configuration{shared_supercell};

  // Check configuration
  EXPECT_EQ(configuration.size(), 4 * 216);

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  bool align = false;
  print_basis_functions(log(), *primclex_ptr, basis_set_name, align);

  // Check clexulator
  EXPECT_EQ(clexulator.name(), "LocalOccClexulatorZrOTest_Clexulator_default");
  EXPECT_EQ(clexulator.nlist_size(), 53);
  EXPECT_EQ(clexulator.corr_size(), 33);
  EXPECT_EQ(clexulator.neighborhood().size(), 26);

  Eigen::VectorXd corr = correlations(configuration, clexulator);
  EXPECT_EQ(corr.size(), 33);
}

void LocalOccClexulatorZrOTest::eval_correlations_test() {
  // Configuration w/ 1 O, 6x6x6 supercell
  CASM::Configuration configuration{shared_supercell};
  configuration.set_occ(2 * 216, 1);  // O

  Clexulator clexulator = primclex_ptr->clexulator(basis_set_name);

  assert_correlations_equivalence(*primclex_ptr, basis_set_name, configuration,
                                  clexulator);
}

void LocalOccClexulatorZrOTest::correlations_correctness_test() {
  std::vector<Clexulator> clexulator =
      primclex_ptr->local_clexulator(basis_set_name);
  EXPECT_EQ(clexulator.size(), 2);

  std::vector<IntegralCluster> phenom_orbit = this->make_phenom_orbit();
  EXPECT_EQ(phenom_orbit.size(), 2);

  // {
  //     // To see the phenomenal clusters:
  //     // - phenom_orbit[0]: [{2, 0, 0, 0}, {3, 0, 0, 0}]
  //     // - phenom_orbit[1]: [{3, 0, 0, 0}, {2, 0, 0, 1}]
  //     jsonParser json;
  //     to_json(phenom_orbit, json);
  //     std::cout << json << std::endl;
  // }

  // {
  //   // To see the sites that make up the local orbits:
  //   jsonParser json = jsonParser::array();
  //   for (Index e = 0; e < clexulator.size(); ++e) {
  //     for (Index i = 0; i < clexulator[e].corr_size(); ++i) {
  //       jsonParser tjson;
  //       tjson["equivalent_clex"] = e;
  //       tjson["linear_orbit_index"] = i;
  //       tjson["sites"] = clexulator[e].site_neighborhood(i);
  //       json.push_back(tjson);
  //     }
  //   }
  //   std::cout << json << std::endl;
  // }

  // default corr = constant term & zeros
  Eigen::VectorXd default_corr = Eigen::VectorXd::Zero(33);
  default_corr(0) = 1.0;

  Eigen::VectorXd corr;
  Eigen::VectorXd expected_corr;
  {
    // 6x6x6 supercell, w/ 1 O on site (orbit 1)

    // local correlations associated with unitcell {0, 0, 0}
    // local_corr[0] == 1 0 ... (O is on a phenom site)
    // local_corr[1] == 1 1/2 0 ...  (O is on orbit 1 point site w/ mult 2)
    CASM::Configuration configuration{shared_supercell};
    occ(configuration, {2, 0, 0, 0}) = 1;  // O

    corr = local_corr({0, 0, 0}, configuration, clexulator[0]);
    expected_corr = default_corr;
    EXPECT_TRUE(almost_equal(corr, expected_corr));

    corr = local_corr({0, 0, 0}, configuration, clexulator[1]);
    expected_corr = default_corr;
    expected_corr(1) = 1. / 2.;
    EXPECT_TRUE(almost_equal(corr, expected_corr));
  }

  {
    // 6x6x6 supercell, w/ 1 O on site (orbit 3; orbit 4)

    // local correlations associated with unitcell {0, 0, 0}
    // local_corr[0] == 1 0 0 1/12 0 ...
    // local_corr[1] == 1 0 0 0 1/12 0 ...

    CASM::Configuration configuration{shared_supercell};
    occ(configuration, {2, -1, -1, 0}) = 1;  // O

    corr = local_corr({0, 0, 0}, configuration, clexulator[0]);
    expected_corr = default_corr;
    expected_corr(3) = 1. / 12.;
    EXPECT_TRUE(almost_equal(corr, expected_corr));

    corr = local_corr({0, 0, 0}, configuration, clexulator[1]);
    expected_corr = default_corr;
    expected_corr(4) = 1. / 12.;
    EXPECT_TRUE(almost_equal(corr, expected_corr));
  }

  {
    // 6x6x6 supercell, w/ 2 O

    CASM::Configuration configuration{shared_supercell};
    occ(configuration, {2, -1, -1, 0}) = 1;  // O
    occ(configuration, {3, -1, -1, 0}) = 1;  // O

    corr = local_corr({0, 0, 0}, configuration, clexulator[0]);
    expected_corr = default_corr;
    expected_corr(3) = 2. / 12.;  // both points in same point orbit
    expected_corr(8) = 1. / 6.;   // pair orbit
    EXPECT_TRUE(almost_equal(corr, expected_corr));

    corr = local_corr({0, 0, 0}, configuration, clexulator[1]);
    expected_corr = default_corr;
    expected_corr(3) = 1. / 12.;  // points in different point orbits
    expected_corr(4) = 1. / 12.;
    expected_corr(9) = 1. / 12.;  // pair orbit
    EXPECT_TRUE(almost_equal(corr, expected_corr));
  }
}

TEST_F(LocalOccClexulatorZrOTest, RunAllTests) {
  // to avoid re-compiling repeatedly
  use_clexulator_test();
  eval_correlations_test();
  correlations_correctness_test();
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

void assert_equal(
    NeighborhoodInfo const &neighborhood_info,
    std::vector<ClexBasisFunctionInfo> const &clex_basis_function_info,
    Eigen::VectorXd const &corr_A, std::string corr_A_name,
    Eigen::VectorXd const &corr_B, std::string corr_B_name) {
  // std::cout << corr_A_name << " size: " << corr_A.size() << std::endl;
  // std::cout << corr_B_name << " size: " << corr_B.size() << std::endl;
  ASSERT_EQ(corr_A.size(), corr_B.size());
  for (Index i = 0; i < corr_A.size(); ++i) {
    // std::cout << i << ": " << corr_A(i) << " " << corr_B(i) << std::endl;
    ASSERT_TRUE(almost_equal(corr_A(i), corr_B(i)));
  }
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
  ConfigDoF const &configdof = configuration.configdof();
  SuperNeighborList const &supercell_neighbor_list =
      configuration.supercell().nlist();
  ClexBasisSpecs const &basis_set_specs =
      primclex.basis_set_specs(basis_set_name);
  Index n_unitcells = supercell_neighbor_list.n_unitcells();

  std::vector<ClexBasisFunctionInfo> clex_basis_function_info;
  std::unique_ptr<NeighborhoodInfo> neighborhood_info;
  MakeClexBasisInfo f{primclex.nlist(), clex_basis_function_info,
                      neighborhood_info};
  for_clex_basis_and_orbits(shared_prim, basis_set_specs, CASM::log(), f);
  ASSERT_EQ(clex_basis_function_info.size(), clexulator.corr_size());
  ASSERT_EQ(neighborhood_info->n_point_corr, clexulator.n_point_corr());

  // jsonParser basis_json;
  // write_clex_basis(shared_prim, basis_set_specs, basis_json);
  // std::cout << "basis: \n" << basis_json << std::endl;
  // std::cout << "n_unitcells: " << n_unitcells << std::endl;
  // std::cout << "n_corr: " << clexulator.corr_size() << std::endl;
  // std::cout << "functions: " << std::endl;
  // for (Index i=0; i<clex_basis_function_info.size(); ++i) {
  //   auto const &info = clex_basis_function_info[i];
  //   std::cout << "c: " << info.linear_orbit_index << " f: " <<
  //   info.linear_function_index << " m: " << info.multiplicity << " g: " <<
  //   info.cluster_invariant_group_indices.size() << std::endl;
  // }

  Eigen::VectorXd C_corr;
  correlations(C_corr, configdof, supercell_neighbor_list, clexulator);
  ASSERT_EQ(C_corr.size(), clexulator.corr_size());

  // -- check extensive correlations --
  Eigen::VectorXd C_ext_corr;
  extensive_correlations(C_ext_corr, configdof, supercell_neighbor_list,
                         clexulator);
  Eigen::VectorXd C_ext_corr_mean = (C_ext_corr / n_unitcells);
  ASSERT_EQ(C_ext_corr_mean.size(), clexulator.corr_size());

  assert_equal(*neighborhood_info, clex_basis_function_info, C_corr, "C_corr",
               C_ext_corr_mean, "C_ext_corr_mean");

  // -- check restricted correlations --
  std::vector<unsigned int> correlation_indices;
  std::vector<unsigned int> other_indices;
  for (Index i = 0; i < clexulator.corr_size(); i++) {
    if (i % 2 == 0) {
      correlation_indices.push_back(i);
    } else {
      other_indices.push_back(i);
    }
  }
  Eigen::VectorXd C_corr_restricted =
      Eigen::VectorXd::Zero(clexulator.corr_size());
  restricted_correlations(C_corr_restricted, configdof, supercell_neighbor_list,
                          clexulator, correlation_indices.data(),
                          end_ptr(correlation_indices));
  for (Index i : correlation_indices) {
    ASSERT_TRUE(almost_equal(C_corr_restricted(i), C_corr(i)));
  }
  for (Index i : other_indices) {
    ASSERT_TRUE(almost_equal(C_corr_restricted(i), 0.));
  }

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
  assert_equal(*neighborhood_info, clex_basis_function_info, C_corr, "C_corr",
               mean, "C_all_corr_contribution_mean");

  // confirm:
  // - all_point_corr for null cluster is 0.0
  // - mean over unit cells of all_point_corr / cluster size = correlations
  Eigen::MatrixXd C_all_point_corr = all_point_corr(configuration, clexulator);
  ASSERT_EQ(C_all_point_corr.rows(), clexulator.n_point_corr() * volume);
  ASSERT_EQ(C_all_point_corr.cols(), clexulator.corr_size());
  sum = Eigen::VectorXd::Zero(clexulator.corr_size());
  for (Index i = 0; i < C_all_point_corr.rows(); ++i) {
    sum += C_all_point_corr.row(i);
  }
  mean = sum / n_unitcells;

  Eigen::VectorXd C_all_point_corr_expected(C_corr.size());
  for (Index j = 0; j < clexulator.corr_size(); ++j) {
    if (j == 0) {
      // null cluster point corr evaluates as 0.:
      C_all_point_corr_expected(0) = 0.0;
    } else {
      // multiple counting functions at each cluster site:
      double cluster_size = clex_basis_function_info[j].prototype.size();
      C_all_point_corr_expected(j) = C_corr(j) * cluster_size;
    }
  }

  assert_equal(*neighborhood_info, clex_basis_function_info,
               C_all_point_corr_expected, "C_all_point_corr_expected", mean,
               "C_all_point_corr_mean");
}
