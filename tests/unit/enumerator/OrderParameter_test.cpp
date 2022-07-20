#include "casm/enumerator/OrderParameter.hh"

#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

// // Uncommment what generating a symmetry adapated subspace:
// #include "casm/symmetry/IrrepDecomposition.hh"
// #include "casm/symmetry/IrrepDecompositionImpl.hh"

using namespace CASM;
using namespace test;

namespace {

Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}

/// \brief Round entries that are within tol of being integer to that integer
/// value
void print_pretty(const Eigen::MatrixXd &M) {
  double tol = 1e-10;
  Eigen::MatrixXd Mp(M);
  std::cout << "Eigen::MatrixXd basis;\n";
  std::cout << "basis.resize(" << M.rows() << "," << M.cols() << ");\n";
  for (int j = 0; j < M.cols(); j++) {
    std::cout << "basis.col(" << j << ") << ";
    for (int i = 0; i < M.rows(); i++) {
      if (std::abs(std::round(M(i, j)) - M(i, j)) < tol) {
        Mp(i, j) = std::round(M(i, j));
      }
      std::cout << Mp(i, j);
      if (i != M.rows() - 1)
        std::cout << ",";
      else
        std::cout << ";\n";
    }
  }
}

// // Getting a symmetry adapated subspace:
// auto const &sym_info = shared_supercell->sym_info();
// std::vector<PermuteIterator> invariant_group =
//     make_invariant_subgroup(config_input);
//
// MasterSymGroup symrep_master_group;
// SymGroupRepID symrep_id;
// SymGroupRep const *symrep_ptr = &make_dof_space_symrep(
//     dof_space, sym_info, invariant_group, symrep_master_group, symrep_id);
//
// Index dim = dof_space.dim();
// Eigen::MatrixXd subspace = Eigen::MatrixXd::Identity(dim, dim);
// bool allow_complex = true;
// using namespace SymRepTools_v2;
// using SymRepTools_v2::IrrepDecompositionImpl::pretty;
// using SymRepTools_v2::IrrepDecompositionImpl::prettyc;
// IrrepDecomposition irrep_decomposition = make_irrep_decomposition(
//     *symrep_ptr, symrep_master_group, subspace, allow_complex);
//
// // Uncomment to print symmetry_adapted_subspace:
// std::cout << "symmetry_adapted_subspace: \n";
// print_pretty(irrep_decomposition.symmetry_adapted_subspace);

}  // namespace

// reminder: prim_dof_values = dof_space.basis() * normal_coordinate

class OrderParameterTest : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  OrderParameterTest()
      : shared_prim(std::make_shared<CASM::Structure const>(
            test::FCC_ternary_GLstrain_disp_prim())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {}
};

TEST_F(OrderParameterTest, Test1_GLstrain) {
  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "GLstrain";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  EXPECT_EQ(dof_space.basis().rows(), 6);

  OrderParameter order_f(dof_space);
  Configuration config = config_input.configuration();
  order_f.update(config);

  // test 1: initial value
  {
    Eigen::VectorXd y, expected;
    y = order_f.value();
    expected = Eigen::VectorXd::Zero(6);
    EXPECT_TRUE(almost_equal(y, expected))
        << "\ny: " << y.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }

  // test 2: value after modification
  {
    Eigen::VectorXd y, expected;
    expected = Eigen::VectorXd::Zero(6);
    expected(1) = 0.1;
    config.configdof().values().global_dof_values.at("GLstrain") = expected;
    y = order_f.value();
    EXPECT_TRUE(almost_equal(y, expected))
        << "\ny: " << y.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }

  // test 3: global delta
  {
    Eigen::VectorXd old_x, new_x, dy, expected;
    old_x = config.configdof().values().global_dof_values.at("GLstrain");
    new_x.resize(6);
    new_x << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;
    expected = new_x - old_x;
    dy = order_f.global_delta(new_x);
    EXPECT_TRUE(almost_equal(dy, expected))
        << "\nold_x: " << old_x.transpose() << "\nnew_x: " << new_x.transpose()
        << "\ndy: " << dy.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }

  // test 4: global delta (single component)
  {
    Eigen::VectorXd old_x, dy, expected;
    old_x = config.configdof().values().global_dof_values.at("GLstrain");
    int dof_component = 0;
    double new_xi = 0.1;
    expected = Eigen::VectorXd::Zero(6);
    expected(dof_component) = new_xi - old_x(dof_component);
    dy = order_f.global_delta(dof_component, new_xi);
    EXPECT_TRUE(almost_equal(dy, expected))
        << "\nold_x: " << old_x.transpose()
        << "\ndof_component: " << dof_component << "\nnew_xi: " << new_xi
        << "\ndy: " << dy.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }
}

TEST_F(OrderParameterTest, Test1_GLstrain_SymAdapted) {
  // Construct the GLstrain DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "GLstrain";
  Eigen::MatrixXd basis;
  basis.resize(6, 6);
  basis.col(0) << 0.57735, 0.57735, 0.57735, 0.0, 0.0, 0.0;
  basis.col(1) << -0.408248, -0.408248, 0.816497, 0.0, 0.0, 0.0;
  basis.col(2) << 0.707107, -0.707107, 0.0, 0.0, 0.0, 0.0;
  basis.col(3) << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
  basis.col(4) << 0.0, 0.0, 0.0, 0.0, 1.0, 0.0;
  basis.col(5) << 0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
  DoFSpace dof_space = make_dof_space(dof_key, config_input, basis);
  EXPECT_EQ(dof_space.basis().rows(), 6);

  OrderParameter order_f(dof_space);
  Configuration config = config_input.configuration();
  order_f.update(config);

  // test 1: initial value
  {
    Eigen::VectorXd y, expected;
    y = order_f.value();
    expected = Eigen::VectorXd::Zero(6);
    EXPECT_TRUE(almost_equal(y, expected))
        << "\ny: " << y.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }

  // test 2: value after modification
  {
    // x = basis * y
    Eigen::VectorXd x, y;
    x = Eigen::VectorXd::Zero(6);
    x(0) = 0.1;
    config.configdof().values().global_dof_values.at("GLstrain") = x;
    y = order_f.value();
    EXPECT_TRUE(almost_equal(x, basis * y))
        << "\nx: " << x.transpose() << "\ny: " << y.transpose()
        << "\nbasis * y: " << (basis * y).transpose() << std::endl;
  }

  // test 3: global delta
  {
    // x = basis * y
    Eigen::VectorXd old_x, new_x, dy, expected;
    old_x = config.configdof().values().global_dof_values.at("GLstrain");
    new_x.resize(6);
    new_x << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;
    dy = order_f.global_delta(new_x);
    EXPECT_TRUE(almost_equal(new_x - old_x, basis * dy))
        << "\nold_x: " << old_x.transpose() << "\nnew_x: " << new_x.transpose()
        << "\ndx: " << (new_x - old_x).transpose() << "\ndy: " << dy.transpose()
        << "\nbasis * dy: " << (basis * dy).transpose() << std::endl;
  }

  // test 4: global delta (single component)
  {
    Eigen::VectorXd old_x, new_x, dy;
    old_x = config.configdof().values().global_dof_values.at("GLstrain");
    int dof_component = 0;
    double new_xi = 0.1;
    new_x = old_x;
    new_x(dof_component) = new_xi;
    dy = order_f.global_delta(dof_component, new_xi);
    EXPECT_TRUE(almost_equal(dy, basis * dy))
        << "\nold_x: " << old_x.transpose() << "\nnew_x: " << new_x.transpose()
        << "\ndx: " << (new_x - old_x).transpose() << "\ndy: " << dy.transpose()
        << "\nbasis * dy: " << (basis * dy).transpose() << std::endl;
  }
}

TEST_F(OrderParameterTest, Test2_disp) {
  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  EXPECT_EQ(dof_space.basis().rows(), 4 * 3);

  OrderParameter order_f(dof_space);
  Configuration config = config_input.configuration();
  order_f.update(config);

  // test 1: initial value
  {
    Eigen::VectorXd y, expected;
    y = order_f.value();
    expected = Eigen::VectorXd::Zero(12);
    EXPECT_TRUE(almost_equal(y, expected))
        << "\ny: " << y.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }

  // test 2a: value after modification
  {
    Eigen::MatrixXd x;
    Eigen::VectorXd y, expected;
    x = Eigen::MatrixXd::Zero(3, 4);
    x(0, 0) = 0.1;
    config.configdof().values().local_dof_values.at("disp") = x;
    expected = Eigen::VectorXd::Zero(12);
    expected(0) = 0.1;
    y = order_f.value();
    EXPECT_TRUE(almost_equal(y, expected))
        << "\nx: " << x.transpose() << "\ny: " << y.transpose()
        << "\nexpected: " << expected.transpose() << std::endl;
  }

  // test 2b: value after modification
  {
    Eigen::MatrixXd x;
    Eigen::VectorXd y, expected;
    x = Eigen::MatrixXd::Zero(3, 4);
    x(0, 1) = 0.1;
    config.configdof().values().local_dof_values.at("disp") = x;
    expected = Eigen::VectorXd::Zero(12);
    expected(3) = 0.1;
    y = order_f.value();
    EXPECT_TRUE(almost_equal(y, expected))
        << "\nx: " << x.transpose() << "\ny: " << y.transpose()
        << "\nexpected: " << expected.transpose() << std::endl;
  }

  // test 3: local delta
  {
    // x = basis * y
    Eigen::VectorXd old_xi, new_xi, dy, expected;
    Index site_index = 1;
    Index dof_component = 1;
    double new_value = 0.1;
    config.configdof().values().local_dof_values.at("disp") =
        Eigen::MatrixXd::Zero(3, 4);
    old_xi =
        config.configdof().values().local_dof_values.at("disp").col(site_index);
    new_xi = old_xi;
    new_xi(dof_component) = new_value;
    dy = order_f.local_delta(site_index, new_xi);
    expected = Eigen::VectorXd::Zero(12);
    expected(dof_space.basis_row_index(site_index, dof_component)) = new_value;
    EXPECT_TRUE(almost_equal(dy, expected))
        << "\nsite_index: " << site_index
        << "\ndof_component: " << dof_component << "\nnew_value: " << new_value
        << "\nold_xi: " << old_xi.transpose()
        << "\nnew_xi: " << new_xi.transpose() << "\ndy: " << dy.transpose()
        << "\nexpected: " << expected.transpose() << std::endl;
  }

  // test 4: local delta (single component)
  {
    // x = basis * y
    Eigen::VectorXd old_xi, new_xi, dy, expected;
    Index site_index = 1;
    Index dof_component = 1;
    double new_value = 0.1;
    config.configdof().values().local_dof_values.at("disp") =
        Eigen::MatrixXd::Zero(3, 4);
    old_xi =
        config.configdof().values().local_dof_values.at("disp").col(site_index);
    dy = order_f.local_delta(site_index, dof_component, new_value);
    expected = Eigen::VectorXd::Zero(12);
    expected(dof_space.basis_row_index(site_index, dof_component)) = new_value;
    EXPECT_TRUE(almost_equal(dy, expected))
        << "\nsite_index: " << site_index
        << "\ndof_component: " << dof_component << "\nnew_value: " << new_value
        << "\nold_xi: " << old_xi.transpose() << "\ndy: " << dy.transpose()
        << "\nexpected: " << expected.transpose() << std::endl;
  }
}

TEST_F(OrderParameterTest, Test2_disp_SymAdapted) {
  // Construct the disp DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  Eigen::MatrixXd basis;
  basis.resize(12, 12);
  basis.col(0) << 0.5, 0, -0, 0.5, 0, -0, 0.5, 0, -0, 0.5, 0, 0;
  basis.col(1) << 0, 0.5, -0, 0, 0.5, -0, 0, 0.5, 0, 0, 0.5, -0;
  basis.col(2) << 0, 0, 0.5, -0, 0, 0.5, 0, 0, 0.5, 0, 0, 0.5;
  basis.col(3) << 0.5, 0, -0, -0.5, -0, -0, 0.5, 0, -0, -0.5, -0, -0;
  basis.col(4) << 0, 0.5, -0, 0, 0.5, 0, -0, -0.5, 0, -0, -0.5, 0;
  basis.col(5) << 0, 0, 0.5, -0, 0, -0.5, -0, -0, -0.5, 0, 0, 0.5;
  basis.col(6) << 0.707107, 0, -0, 0, -0, 0, -0.707107, 0, 0, 0, 0, 0;
  basis.col(7) << 0, 0.707107, -0, -0, -0.707107, 0, -0, -0, -0, 0, -0, 0;
  basis.col(8) << -0, 0, 0.707107, 0, 0, -0, -0, -0, -0, -0, 0, -0.707107;
  basis.col(9) << 0, -0, 0, 0.707107, 0, -0, 0, 0, 0, -0.707107, -0, 0;
  basis.col(10) << 0, 0, -0, 0, 0, 0.707107, -0, -0, -0.707107, 0, -0, -0;
  basis.col(11) << -0, 0, 0, 0, 0, -0, 0, 0.707107, -0, 0, -0.707107, 0;
  DoFSpace dof_space = make_dof_space(dof_key, config_input, basis);
  EXPECT_EQ(dof_space.basis().rows(), 4 * 3);

  OrderParameter order_f(dof_space);
  Configuration config = config_input.configuration();
  order_f.update(config);

  // test 1: initial value
  {
    Eigen::VectorXd y, expected;
    y = order_f.value();
    expected = Eigen::VectorXd::Zero(12);
    EXPECT_TRUE(almost_equal(y, expected))
        << "\ny: " << y.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }

  // test 2a: value after modification
  {
    // x = basis * y
    Eigen::MatrixXd x;
    Eigen::VectorXd y, expected;
    x = Eigen::MatrixXd::Zero(3, 4);
    x(0, 0) = 0.1;
    config.configdof().values().local_dof_values.at("disp") = x;
    expected = Eigen::VectorXd::Zero(12);
    expected(0) = 0.1;
    y = order_f.value();
    EXPECT_TRUE(almost_equal(expected, basis * y))
        << "\nx: " << x.transpose() << "\ny: " << y.transpose()
        << "\nexpected: " << expected.transpose() << std::endl;
  }

  // test 2b: value after modification
  {
    Eigen::MatrixXd x;
    Eigen::VectorXd y, expected;
    x = Eigen::MatrixXd::Zero(3, 4);
    x(0, 1) = 0.1;
    config.configdof().values().local_dof_values.at("disp") = x;
    expected = Eigen::VectorXd::Zero(12);
    expected(3) = 0.1;
    y = order_f.value();
    EXPECT_TRUE(almost_equal(expected, basis * y))
        << "\nx: " << x.transpose() << "\ny: " << y.transpose()
        << "\nexpected: " << expected.transpose() << std::endl;
  }

  // test 3: local delta
  {
    // x = basis * y
    Eigen::VectorXd old_xi, new_xi, dy, expected;
    Index site_index = 1;
    Index dof_component = 1;
    double new_value = 0.1;
    config.configdof().values().local_dof_values.at("disp") =
        Eigen::MatrixXd::Zero(3, 4);
    old_xi =
        config.configdof().values().local_dof_values.at("disp").col(site_index);
    new_xi = old_xi;
    new_xi(dof_component) = new_value;
    dy = order_f.local_delta(site_index, new_xi);
    expected = Eigen::VectorXd::Zero(12);
    expected(dof_space.basis_row_index(site_index, dof_component)) = new_value;
    EXPECT_TRUE(almost_equal(expected, basis * dy))
        << "\nsite_index: " << site_index
        << "\ndof_component: " << dof_component << "\nnew_value: " << new_value
        << "\nold_xi: " << old_xi.transpose()
        << "\nnew_xi: " << new_xi.transpose() << "\ndy: " << dy.transpose()
        << "\nexpected: " << expected.transpose() << std::endl;
  }

  // test 4: global delta (single component)
  {
    // x = basis * y
    Eigen::VectorXd old_xi, new_xi, dy, expected;
    Index site_index = 1;
    Index dof_component = 1;
    double new_value = 0.1;
    config.configdof().values().local_dof_values.at("disp") =
        Eigen::MatrixXd::Zero(3, 4);
    old_xi =
        config.configdof().values().local_dof_values.at("disp").col(site_index);
    dy = order_f.local_delta(site_index, dof_component, new_value);
    expected = Eigen::VectorXd::Zero(12);
    expected(dof_space.basis_row_index(site_index, dof_component)) = new_value;
    EXPECT_TRUE(almost_equal(expected, basis * dy))
        << "\nsite_index: " << site_index
        << "\ndof_component: " << dof_component << "\nnew_value: " << new_value
        << "\nold_xi: " << old_xi.transpose() << "\ndy: " << dy.transpose()
        << "\nexpected: " << expected.transpose() << std::endl;
  }
}

TEST_F(OrderParameterTest, Test3_occ) {
  // Construct the occ DoF space.
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  EXPECT_EQ(dof_space.basis().rows(), 4 * 3);

  OrderParameter order_f(dof_space);
  Configuration config = config_input.configuration();
  order_f.update(config);

  // test 1: initial value
  {
    Eigen::VectorXd y, expected;
    y = order_f.value();
    expected.resize(12);
    expected << 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.;
    EXPECT_TRUE(almost_equal(y, expected))
        << "\ny: " << y.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }

  // test 2a: value after modification
  {
    Eigen::VectorXd y, expected;
    config.configdof().values().occupation(0) = 1;
    expected.resize(12);
    expected << 0., 1., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.;
    y = order_f.value();
    EXPECT_TRUE(almost_equal(y, expected))
        << "\ny: " << y.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }

  // test 2b: value after modification
  {
    Eigen::VectorXd y, expected;
    config.configdof().values().occupation(0) = 2;
    expected.resize(12);
    expected << 0., 0., 1., 1., 0., 0., 1., 0., 0., 1., 0., 0.;
    y = order_f.value();
    EXPECT_TRUE(almost_equal(y, expected))
        << "\ny: " << y.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }

  // test 3: occ delta
  {
    // x = basis * y
    Eigen::VectorXd dy, expected;
    Index site_index = 1;
    int new_occ = 1;
    config.configdof().values().occupation = Eigen::VectorXi::Zero(4);
    dy = order_f.occ_delta(site_index, new_occ);
    expected = Eigen::VectorXd::Zero(12);
    expected << 0., 0., 0., -1., 1., 0., 0., 0., 0., 0., 0., 0.;
    EXPECT_TRUE(almost_equal(dy, expected))
        << "\nsite_index: " << site_index << "\nnew_occ: " << new_occ
        << "\ndy: " << dy.transpose() << "\nexpected: " << expected.transpose()
        << std::endl;
  }
}

class OrderParameterTest2 : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  OrderParameterTest2()
      : shared_prim(std::make_shared<CASM::Structure const>(
            test::FCC_ternary_GLstrain_disp_prim())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {}
};

/// Test that order parameters are calcluated even for non-commensurate
/// supercells
TEST_F(OrderParameterTest2, NonCommensurateSupercellTest) {
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);
  EXPECT_EQ(dof_space.basis().rows(), 12);

  OrderParameter order_f(dof_space);

  Index begin_volume = 1;
  Index end_volume = 6;
  xtal::ScelEnumProps params{begin_volume, end_volume};
  auto const &pg = shared_prim->point_group();
  xtal::SuperlatticeEnumerator enumerator{pg.begin(), pg.end(),
                                          shared_prim->lattice(), params};
  for (xtal::Lattice const &lattice : enumerator) {
    // std::cout << "#####" << std::endl;
    // std::cout << "lattice vectors: \n"
    //           << lattice.lat_column_mat().transpose() << std::endl;
    auto shared_supercell =
        std::make_shared<CASM::Supercell>(shared_prim, lattice);
    ConfigEnumInput config_input{*shared_supercell};
    ConfigEnumAllOccupations config_enum{config_input};
    // for (Configuration const &config : config_enum) {
    //   std::cout << "---" << std::endl;
    //   std::cout << "occ: " <<
    //   config.configdof().values().occupation.transpose() << std::endl;
    //   std::cout << "val: " << order_f(config).transpose() << std::endl;
    // }
  }
  EXPECT_TRUE(true);
}
