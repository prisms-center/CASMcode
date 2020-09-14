#include "gtest/gtest.h"
#include "autotools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "crystallography/TestStructures.hh"

using namespace CASM;
using namespace test;

namespace {
  Eigen::Matrix3i _fcc_conventional_transf_mat() {
    Eigen::Matrix3i transf_mat;
    transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    return transf_mat;
  }
}

class ConfigEnumInputTest : public testing::Test {
protected:

  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  ConfigEnumInputTest():
    shared_prim(std::make_shared<CASM::Structure const>(test::FCC_ternary_prim())),
    shared_supercell(std::make_shared<CASM::Supercell>(shared_prim, _fcc_conventional_transf_mat())) {}

};

TEST_F(ConfigEnumInputTest, ConstructorTest1) {

  // Constructor, by supercell, no explicit site selection, expect:
  // - zeros ConfigDoF
  // - all sites selected
  CASM::ConfigEnumInput config_enum_input {*shared_supercell};
  EXPECT_EQ(config_enum_input.sites().size(), 4);
  ConfigDoF const &configdof = config_enum_input.configuration().configdof();
  for(int i = 0; i < configdof.size(); ++i) {
    EXPECT_EQ(configdof.occ(i), 0);
  }

  std::vector<PermuteIterator> invariant_group = make_invariant_group(config_enum_input);
  EXPECT_EQ(invariant_group.size(), 48 * 4);
}

TEST_F(ConfigEnumInputTest, ConstructorTest2) {

  CASM::Configuration configuration {shared_supercell};
  configuration.configdof().occ(0) = 1; // L12: A3B1

  // Constructor, by configuration, no explicit site selection, expect:
  // - matching ConfigDoF
  // - all sites selected
  CASM::ConfigEnumInput config_enum_input {configuration};
  EXPECT_EQ(config_enum_input.sites().size(), 4);

  ConfigDoF const &configdof = config_enum_input.configuration().configdof();
  for(int i = 0; i < configdof.size(); ++i) {
    EXPECT_EQ(configdof.occ(i), configuration.configdof().occ(i));
  }

  std::vector<PermuteIterator> invariant_group = make_invariant_group(config_enum_input);
  EXPECT_EQ(invariant_group.size(), 48);
}

TEST_F(ConfigEnumInputTest, ConstructorTest3) {

  // Constructor, by supercell, no explicit site selection, expect:
  // - zeros ConfigDoF
  // - 1 site selected
  CASM::ConfigEnumInput config_enum_input {*shared_supercell, {0}};
  EXPECT_EQ(config_enum_input.sites().size(), 1);
  ConfigDoF const &configdof = config_enum_input.configuration().configdof();
  for(int i = 0; i < configdof.size(); ++i) {
    EXPECT_EQ(configdof.occ(i), 0);
  }

  std::vector<PermuteIterator> invariant_group = make_invariant_group(config_enum_input);
  EXPECT_EQ(invariant_group.size(), 48);
}

TEST_F(ConfigEnumInputTest, ConstructorTest4) {

  CASM::Configuration configuration {shared_supercell};
  configuration.configdof().occ(0) = 1; // L12: A3B1

  // Constructor, by configuration, no explicit site selection, expect:
  // - matching ConfigDoF
  // - 1 site selected (not the B site)
  CASM::ConfigEnumInput config_enum_input {configuration, {1}};
  EXPECT_EQ(config_enum_input.sites().size(), 1);

  ConfigDoF const &configdof = config_enum_input.configuration().configdof();
  for(int i = 0; i < configdof.size(); ++i) {
    EXPECT_EQ(configdof.occ(i), configuration.configdof().occ(i));
  }

  std::vector<PermuteIterator> invariant_group = make_invariant_group(config_enum_input);
  EXPECT_EQ(invariant_group.size(), 16);
}
