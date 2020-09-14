#include "gtest/gtest.h"
#include "autotools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ClusterSitesSelector_impl.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "crystallography/TestStructures.hh"

using namespace CASM;
using namespace test;

class ClusterSitesSelectorTest : public testing::Test {
protected:

  std::shared_ptr<CASM::Structure const> shared_prim;
  Configuration configuration_L12;
  Configuration configuration_L12_3x3x3;

  ConfigEnumInputTest():
    shared_prim(std::make_shared<CASM::Structure const>(test::FCC_ternary_prim())),
    shared_conventional_fcc_supercell(std::make_shared<CASM::Supercell>(shared_prim, _fcc_conventional_transf_mat())),
    configuration_L12(_make_configuration_L12()),
    configuration_L12_3x3x3(_make_configuration_L12_3x3x3()) {}

  Eigen::Matrix3i _fcc_conventional_transf_mat() {
    Eigen::Matrix3i transf_mat;
    transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    return transf_mat;
  }

  Configuration _make_configuration_L12() {
    auto conventional_fcc = std::make_shared<Supercell>(this->shared_prim, _fcc_conventional_transf_mat());
    Configuration configuration {conventional_fcc};
    configuration.configdof().occ(0) = 1;
    return configuration;
  }

  Configuration _make_configuration_L12_3x3x3() {
    auto conventional_fcc_3x3x3 = std::make_shared<Supercell>(this->shared_prim, 3 * _fcc_conventional_transf_mat());
    FillSupercell filler {conventional_fcc_3x3x3, this->config_L12, TOL};
    return filler(this->config_L12);
  }

};

TEST_F(ClusterSitesSelectorTest, Test1) {
  EXPECT_EQ(1, 1);
}
