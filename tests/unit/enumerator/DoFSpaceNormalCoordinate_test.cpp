
#include "casm/casm_io/Log.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/enumerator/ConfigEnumInput_impl.hh"
#include "casm/enumerator/DoFSpace.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "casm/symmetry/VectorSpaceSymReport.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;
using namespace test;

namespace {

Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}

}  // namespace

class NormalCoordinateTest0 : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  std::optional<Configuration> config;
  std::optional<DoFSpace> dof_space;

  NormalCoordinateTest0()
      : shared_prim(
            std::make_shared<CASM::Structure const>(test::FCC_binary_prim())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {
    config = Configuration(shared_supercell);

    // occ DoF space, with all sites
    ConfigEnumInput config_input{*config};
    DoFKey dof_key = "occ";
    DoFSpace tmp_dof_space = make_dof_space(dof_key, config_input);

    // symmetry adapted DoF space
    auto const &sym_info = config_input.configuration().supercell().sym_info();
    std::vector<PermuteIterator> group = make_invariant_subgroup(config_input);
    bool calc_wedges = false;
    std::optional<SymRepTools_v2::VectorSpaceSymReport> sym_report;
    dof_space = make_symmetry_adapted_dof_space_v2(
        tmp_dof_space, sym_info, group, calc_wedges, sym_report);

    // std::cout << "symmetry adapted basis: \n" <<
    // symmetry_adapted_dof_space.basis() << std::endl; jsonParser json;
    // to_json(symmetry_adapted_dof_space, json, "test", config_input,
    // sym_report); std::cout << json << std::endl;
  }
};

TEST_F(NormalCoordinateTest0, SymmetryAdaptedTest0) {
  Index eta1_i = 0;
  Index eta2_i = 5;
  Index eta3_i = 6;
  Index eta4_i = 7;

  {
    Configuration tconfig = *config;
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta0: " << eta.transpose() << std::endl;
    // std::cout << "x0: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], sqrt(2.)));
    EXPECT_TRUE(almost_equal(eta[eta2_i], 0.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], 0.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], 0.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(0, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta1_1: " << eta.transpose() << std::endl;
    // std::cout << "x1_1: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta2_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], -sqrt(2.) / 2.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(1, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta1_2: " << eta.transpose() << std::endl;
    // std::cout << "x1_2: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta2_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], sqrt(2.) / 2.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(2, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta1_3: " << eta.transpose() << std::endl;
    // std::cout << "x1_3: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta2_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], sqrt(2.) / 2.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(3, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta1_4: " << eta.transpose() << std::endl;
    // std::cout << "x1_4: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta2_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], -sqrt(2.) / 2.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(0, 1);
    tconfig.set_occ(1, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta2: " << eta.transpose() << std::endl;
    // std::cout << "x2: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], 0.));
    EXPECT_TRUE(almost_equal(eta[eta2_i], -sqrt(2.)));
    EXPECT_TRUE(almost_equal(eta[eta3_i], 0.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], 0.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(0, 1);
    tconfig.set_occ(1, 1);
    tconfig.set_occ(2, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta3_1: " << eta.transpose() << std::endl;
    // std::cout << "x3_1: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta2_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], sqrt(2.) / 2.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(0, 1);
    tconfig.set_occ(1, 1);
    tconfig.set_occ(3, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta3_2: " << eta.transpose() << std::endl;
    // std::cout << "x3_2: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta2_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], -sqrt(2.) / 2.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(0, 1);
    tconfig.set_occ(2, 1);
    tconfig.set_occ(3, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta3_3: " << eta.transpose() << std::endl;
    // std::cout << "x3_3: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta2_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], -sqrt(2.) / 2.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(1, 1);
    tconfig.set_occ(2, 1);
    tconfig.set_occ(3, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta3_4: " << eta.transpose() << std::endl;
    // std::cout << "x3_4: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], -sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta2_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], sqrt(2.) / 2.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], sqrt(2.) / 2.));
  }

  {
    Configuration tconfig = *config;
    tconfig.set_occ(0, 1);
    tconfig.set_occ(1, 1);
    tconfig.set_occ(2, 1);
    tconfig.set_occ(3, 1);
    Eigen::VectorXd eta = get_normal_coordinate(tconfig, *dof_space);
    Eigen::VectorXd x = dof_space->basis() * eta;
    // std::cout << "eta4: " << eta.transpose() << std::endl;
    // std::cout << "x4: " << x.transpose() << std::endl << std::endl;
    EXPECT_TRUE(almost_equal(eta[eta1_i], -sqrt(2.)));
    EXPECT_TRUE(almost_equal(eta[eta2_i], 0.));
    EXPECT_TRUE(almost_equal(eta[eta3_i], 0.));
    EXPECT_TRUE(almost_equal(eta[eta4_i], 0.));
  }
}

class NormalCoordinateTest1 : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  NormalCoordinateTest1()
      : shared_prim(std::make_shared<CASM::Structure const>(
            test::FCC_ternary_GLstrain_disp_prim())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())) {}
};

TEST_F(NormalCoordinateTest1, Test1) {
  // default config
  Configuration config{shared_supercell};

  // occ DoF space, with all sites
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  Eigen::VectorXd x = get_normal_coordinate(config, dof_space);

  Eigen::VectorXd x_expected = Eigen::VectorXd::Zero(12);
  x_expected << 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.;
  EXPECT_EQ(x.size(), x_expected.size());
  EXPECT_TRUE(almost_equal(x, x_expected, TOL))
      << "x: " << x.transpose() << "\n"
      << "x_expected: " << x_expected.transpose() << std::endl;
}

TEST_F(NormalCoordinateTest1, Test2) {
  // non-default config
  Configuration config{shared_supercell};
  config.set_occ(0, 2);

  // occ DoFSpace, with all sites
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  Eigen::VectorXd x = get_normal_coordinate(config, dof_space);

  Eigen::VectorXd x_expected = Eigen::VectorXd::Zero(12);
  x_expected << 0., 0., 1., 1., 0., 0., 1., 0., 0., 1., 0., 0.;
  EXPECT_EQ(x.size(), x_expected.size());
  EXPECT_TRUE(almost_equal(x, x_expected, TOL))
      << "x: " << x.transpose() << "\n"
      << "x_expected: " << x_expected.transpose() << std::endl;
}

TEST_F(NormalCoordinateTest1, Test3) {
  // non-default config
  Configuration config{shared_supercell};
  config.set_occ(1, 2);
  config.set_occ(2, 1);

  // occ DoFSpace, with subset of sites
  ConfigEnumInput config_input{*shared_supercell, {1, 2}};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  Eigen::VectorXd x = get_normal_coordinate(config, dof_space);

  Eigen::VectorXd x_expected = Eigen::VectorXd::Zero(6);
  x_expected << 0., 0., 1., 0., 1., 0.;
  EXPECT_EQ(x.size(), x_expected.size());
  EXPECT_TRUE(almost_equal(x, x_expected, TOL))
      << "x: " << x.transpose() << "\n"
      << "x_expected: " << x_expected.transpose() << std::endl;
}

TEST_F(NormalCoordinateTest1, Test4) {
  // default config - different supercell than DoF space
  auto shared_prim_supercell = std::make_shared<CASM::Supercell>(
      shared_prim, Eigen::Matrix3l::Identity());
  Configuration config{shared_prim_supercell};

  // occ DoF space, with all sites
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // `get_normal_coordinate` throws if config and dof_space have different
  // supercells
  EXPECT_THROW(get_normal_coordinate(config, dof_space), std::runtime_error);
}

TEST_F(NormalCoordinateTest1, Test5) {
  // default config - different (smaller) supercell than DoF space
  auto shared_prim_supercell = std::make_shared<CASM::Supercell>(
      shared_prim, Eigen::Matrix3l::Identity());
  Configuration config{shared_prim_supercell};

  // occ DoF space, with all sites
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // `get_mean_normal_coordinate`
  Eigen::VectorXd x_mean = get_mean_normal_coordinate(config, dof_space);

  Eigen::VectorXd x_mean_expected = Eigen::VectorXd::Zero(12);
  x_mean_expected << 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.;
  EXPECT_EQ(x_mean.size(), x_mean_expected.size());
  EXPECT_TRUE(almost_equal(x_mean, x_mean_expected, TOL))
      << "x_mean: " << x_mean.transpose() << "\n"
      << "x_mean_expected: " << x_mean_expected.transpose() << std::endl;
}

TEST_F(NormalCoordinateTest1, Test6) {
  // default config - different (larger) supercell than DoF space
  auto conventional_fcc_vol8 = std::make_shared<CASM::Supercell>(
      shared_prim, _fcc_conventional_transf_mat() * 2);
  Configuration config{conventional_fcc_vol8};

  // occ DoF space, with all sites
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // `get_mean_normal_coordinate`
  Eigen::VectorXd x_mean = get_mean_normal_coordinate(config, dof_space);

  Eigen::VectorXd x_mean_expected = Eigen::VectorXd::Zero(12);
  x_mean_expected << 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.;
  EXPECT_EQ(x_mean.size(), x_mean_expected.size());
  EXPECT_TRUE(almost_equal(x_mean, x_mean_expected, TOL))
      << "x_mean: " << x_mean.transpose() << "\n"
      << "x_mean_expected: " << x_mean_expected.transpose() << std::endl;
}

TEST_F(NormalCoordinateTest1, Test7) {
  // non-default config - different (larger) supercell than DoF space
  auto conventional_fcc_vol8 = std::make_shared<CASM::Supercell>(
      shared_prim, _fcc_conventional_transf_mat() * 2);
  Configuration config{shared_supercell};
  config.set_occ(0, 2);
  Configuration super_config = fill_supercell(config, conventional_fcc_vol8);

  // occ DoF space, with all sites
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // `get_mean_normal_coordinate`
  Eigen::VectorXd x_mean = get_mean_normal_coordinate(config, dof_space);

  Eigen::VectorXd x_mean_expected = Eigen::VectorXd::Zero(12);
  x_mean_expected << 0., 0., 1., 1., 0., 0., 1., 0., 0., 1., 0., 0.;
  EXPECT_EQ(x_mean.size(), x_mean_expected.size());
  EXPECT_TRUE(almost_equal(x_mean, x_mean_expected, TOL))
      << "x_mean: " << x_mean.transpose()
      << "\nx_mean_expected: " << x_mean_expected.transpose() << std::endl;
}

TEST_F(NormalCoordinateTest1, Test8) {
  auto supercell_1x1x1 = std::make_shared<CASM::Supercell>(
      shared_prim, _fcc_conventional_transf_mat());
  Lattice const &fcc_conventional_lattice = supercell_1x1x1->lattice();

  Eigen::Matrix3l T;
  T << 2., 0., 0., 0., 1., 0., 0., 0., 1.;
  auto supercell_2x1x1 = std::make_shared<CASM::Supercell>(
      shared_prim, _fcc_conventional_transf_mat() * T);

  T << 3., 0., 0., 0., 1., 0., 0., 0., 1.;
  auto supercell_3x1x1 = std::make_shared<CASM::Supercell>(
      shared_prim, _fcc_conventional_transf_mat() * T);

  // default config,
  // commensurate supercell (6x1x1) is different than config supercell
  Configuration config{supercell_3x1x1};

  // occ DoF space, with all sites
  ConfigEnumInput config_input{*supercell_2x1x1};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // check 1) default config

  // `get_mean_normal_coordinate`
  Eigen::VectorXd x_mean = get_mean_normal_coordinate(config, dof_space);

  // occupation on conventional cell origin in unit cells along x:
  //                 config: [0,  0,  0] [0,  0,  0]
  // commensurate supercell: [0,  0,  0,  0,  0,  0]
  //              dof_space: [-,  -] [-,  -] [-,  -]

  Eigen::VectorXd x_mean_expected = Eigen::VectorXd::Zero(24);
  x_mean_expected << 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0., 1., 0., 0.,
      1., 0., 0., 1., 0., 0., 1., 0., 0.;
  EXPECT_EQ(x_mean.size(), x_mean_expected.size());
  EXPECT_TRUE(almost_equal(x_mean, x_mean_expected, TOL))
      << "x_mean: " << x_mean.transpose() << "\n"
      << "x_mean_expected: " << x_mean_expected.transpose() << std::endl;

  // check 2) non-default config

  DoFSpaceIndexConverter f{config, dof_space};

  {
    xtal::Coordinate coord{0., 0., 0., fcc_conventional_lattice, FRAC};
    Index l = f.config_site_index(coord, TOL);
    config.set_occ(l, 2);
  }

  x_mean = get_mean_normal_coordinate(config, dof_space);

  // occupation on conventional cell origin in unit cells along x:
  //                 config: [2,  0,  0] [2,  0,  0]
  // commensurate supercell: [0,  0,  0,  0,  0,  0]
  //              dof_space: [-,  -] [-,  -] [-,  -]

  // x_mean_expected ( frac(0., 0., 0.) ) = 2./3., 0., 1./3.
  {
    xtal::Coordinate coord{0., 0., 0., fcc_conventional_lattice, FRAC};
    Index l = f.dof_space_site_index(coord, TOL);
    x_mean_expected(dof_space.basis_row_index(l, 0)) = 2. / 3.;
    x_mean_expected(dof_space.basis_row_index(l, 1)) = 0.;
    x_mean_expected(dof_space.basis_row_index(l, 2)) = 1. / 3.;
  }

  // x_mean_expected ( frac(1., 0., 0.) ) = 2./3., 0., 1./3.
  {
    xtal::Coordinate coord{1., 0., 0., fcc_conventional_lattice, FRAC};
    Index l = f.dof_space_site_index(coord, TOL);
    x_mean_expected(dof_space.basis_row_index(l, 0)) = 2. / 3.;
    x_mean_expected(dof_space.basis_row_index(l, 1)) = 0.;
    x_mean_expected(dof_space.basis_row_index(l, 2)) = 1. / 3.;
  }

  EXPECT_EQ(x_mean.size(), x_mean_expected.size());
  EXPECT_TRUE(almost_equal(x_mean, x_mean_expected, TOL))
      << "x_mean: " << x_mean.transpose() << "\n"
      << "x_mean_expected: " << x_mean_expected.transpose() << std::endl;
}

TEST_F(NormalCoordinateTest1, Test9) {
  // test loop over supercells and configurations... not sure what to check
  // besides running without error

  // occ DoF space, conventional FCC with all sites
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // lattice enumeration
  int begin_volume = 1;
  int end_volume = 9;
  std::string dirs = "abc";
  Eigen::Matrix3i generating_matrix = Eigen::Matrix3i::Identity();
  xtal::ScelEnumProps enumeration_params{begin_volume, end_volume, dirs,
                                         generating_matrix};
  CASM::ScelEnumByProps scel_enumerator{shared_prim, enumeration_params};

  for (auto const supercell : scel_enumerator) {
    // std::cout << "\n\n--- --- --- --- --- --- ---\n" << std::endl;

    shared_supercell = std::make_shared<CASM::Supercell>(supercell);
    ConfigEnumAllOccupations config_enumerator{*shared_supercell};

    // std::cout << "Lattice: \n" <<
    // shared_supercell->lattice().lat_column_mat() << std::endl << std::endl;

    for (auto const config : config_enumerator) {
      Eigen::VectorXd x_mean = get_mean_normal_coordinate(config, dof_space);
      // std::cout << "occ: " << config.occupation().transpose() << std::endl;
      // std::cout << "x_mean: " << x_mean.transpose() << std::endl <<
      // std::endl;
    }
  }
}

TEST_F(NormalCoordinateTest1, SymmetryAdaptedTest1) {
  // default config
  Configuration config{shared_supercell};

  // occ DoF space, with all sites
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "occ";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  // symmetry adapted DoF space
  auto const &sym_info = config_input.configuration().supercell().sym_info();
  std::vector<PermuteIterator> group = make_invariant_subgroup(config_input);
  bool calc_wedges = false;
  std::optional<SymRepTools_v2::VectorSpaceSymReport> sym_report;
  DoFSpace symmetry_adapted_dof_space = make_symmetry_adapted_dof_space_v2(
      dof_space, sym_info, group, calc_wedges, sym_report);

  // std::cout << "symmetry adapted basis: \n" <<
  // symmetry_adapted_dof_space.basis().transpose() << std::endl; jsonParser
  // json; to_json(symmetry_adapted_dof_space, json, "test", config_input,
  // sym_report); std::cout << json << std::endl;

  Index eta1_i = 0;
  Index eta2_i = 1;
  Index eta3_i = 2;

  {
    {
      Configuration tconfig = config;
      Eigen::VectorXd eta =
          get_normal_coordinate(tconfig, symmetry_adapted_dof_space);
      Eigen::VectorXd x = symmetry_adapted_dof_space.basis() * eta;
      // std::cout << "eta0_0: " << eta.transpose() << std::endl;
      // std::cout << "x0_0: " << x.transpose() << std::endl << std::endl;
      EXPECT_TRUE(almost_equal(eta[eta1_i], 0.));
      EXPECT_TRUE(almost_equal(eta[eta2_i], sqrt(2.)));
      EXPECT_TRUE(almost_equal(eta[eta3_i], sqrt(2.)));
    }

    {
      Configuration tconfig = config;
      tconfig.set_occ(0, 1);
      tconfig.set_occ(1, 1);
      tconfig.set_occ(2, 1);
      tconfig.set_occ(3, 1);
      Eigen::VectorXd eta =
          get_normal_coordinate(tconfig, symmetry_adapted_dof_space);
      Eigen::VectorXd x = symmetry_adapted_dof_space.basis() * eta;
      // std::cout << "eta0_1: " << eta.transpose() << std::endl;
      // std::cout << "x0_1: " << x.transpose() << std::endl << std::endl;
      EXPECT_TRUE(almost_equal(eta[eta1_i], 2.));
      EXPECT_TRUE(almost_equal(eta[eta2_i], 0.));
      EXPECT_TRUE(almost_equal(eta[eta3_i], 0.));
    }

    {
      Configuration tconfig = config;
      tconfig.set_occ(0, 2);
      tconfig.set_occ(1, 2);
      tconfig.set_occ(2, 2);
      tconfig.set_occ(3, 2);
      Eigen::VectorXd eta =
          get_normal_coordinate(tconfig, symmetry_adapted_dof_space);
      Eigen::VectorXd x = symmetry_adapted_dof_space.basis() * eta;
      // std::cout << "eta0_2: " << eta.transpose() << std::endl;
      // std::cout << "x0_2: " << x.transpose() << std::endl << std::endl;
      EXPECT_TRUE(almost_equal(eta[eta1_i], 0.));
      EXPECT_TRUE(almost_equal(eta[eta2_i], -sqrt(2.)));
      EXPECT_TRUE(almost_equal(eta[eta3_i], sqrt(2.)));
    }
  }
}

TEST_F(NormalCoordinateTest1, DispTest1) {
  // default config
  Configuration config{shared_supercell};

  // disp DoF space, with all sites
  ConfigEnumInput config_input{*shared_supercell};
  DoFKey dof_key = "disp";
  DoFSpace dof_space = make_dof_space(dof_key, config_input);

  Eigen::VectorXd x = get_normal_coordinate(config, dof_space);

  Eigen::VectorXd x_expected = Eigen::VectorXd::Zero(12);
  x_expected << 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.;
  EXPECT_EQ(x.size(), x_expected.size());
  EXPECT_TRUE(almost_equal(x, x_expected, TOL))
      << "x: " << x.transpose() << "\n"
      << "x_expected: " << x_expected.transpose() << std::endl;
}
