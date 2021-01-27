#include "casm/clex/ConfigDoF.hh"

#include "casm/clex/ConfigDoFTools.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

TEST(ConfigDoFTest, ConstructorFail) {
  // fail construction: shared prim structure with no basis sites: cannot
  // construct supercell or configdof

  auto shared_prim = std::make_shared<Structure const>(test::no_basis_prim());
  ASSERT_ANY_THROW(std::make_shared<Supercell const>(
      shared_prim, Eigen::Matrix3l::Identity()));
}

TEST(ConfigDoFTest, Constructor0) {
  // fail construction: shared prim structure with basis sites but no occupant
  // DoF

  auto shared_prim = std::make_shared<Structure const>(test::no_dof_prim());
  auto shared_supercell = std::make_shared<Supercell const>(
      shared_prim, Eigen::Matrix3l::Identity());

  ConfigDoF configdof = make_configdof(*shared_supercell);

  EXPECT_EQ(
      configdof.has_occupation(),
      true);  // TODO: note this is inconsisten with Site having no occupant DoF
  EXPECT_EQ(configdof.global_dofs().size(), 0);
  EXPECT_EQ(configdof.local_dofs().size(), 0);
}

TEST(ConfigDoFTest, Constructor1) {
  // construction: shared prim structure with occupation DoF only

  auto shared_prim =
      std::make_shared<Structure const>(test::FCC_ternary_prim());
  auto shared_supercell = std::make_shared<Supercell const>(
      shared_prim, Eigen::Matrix3l::Identity());

  ConfigDoF configdof = make_configdof(*shared_supercell);

  EXPECT_EQ(configdof.has_occupation(), true);
  EXPECT_EQ(configdof.occupation(),
            Eigen::VectorXi::Zero(shared_supercell->num_sites()));
  EXPECT_EQ(configdof.global_dofs().size(), 0);
  EXPECT_EQ(configdof.local_dofs().size(), 0);
}

TEST(ConfigDoFTest, Constructor2) {
  // construction: shared prim structure with GLstrain DoF only

  auto shared_prim =
      std::make_shared<Structure const>(test::SimpleCubic_GLstrain_prim());
  auto shared_supercell = std::make_shared<Supercell const>(
      shared_prim, Eigen::Matrix3l::Identity());

  ConfigDoF configdof = make_configdof(*shared_supercell);

  EXPECT_EQ(configdof.has_occupation(), true);
  EXPECT_EQ(configdof.global_dofs().size(), 1);
  EXPECT_EQ(configdof.global_dofs().count("GLstrain"), 1);
  EXPECT_EQ(configdof.local_dofs().size(), 0);
}

TEST(ConfigDoFTest, Constructor3) {
  // construction: shared prim structure with GLstrain DoF only

  auto shared_prim =
      std::make_shared<Structure const>(test::SimpleCubic_GLstrain_prim());
  auto shared_supercell = std::make_shared<Supercell const>(
      shared_prim, Eigen::Matrix3l::Identity());

  ConfigDoF configdof = make_configdof(*shared_supercell);

  EXPECT_EQ(configdof.has_occupation(), true);
  EXPECT_EQ(configdof.global_dofs().size(), 1);
  EXPECT_EQ(configdof.global_dofs().count("GLstrain"), 1);
  EXPECT_EQ(configdof.local_dofs().size(), 0);
}

TEST(ConfigDoFTest, Constructor4) {
  // construction: shared prim structure with disp DoF only

  auto shared_prim =
      std::make_shared<Structure const>(test::SimpleCubic_disp_prim());
  auto shared_supercell = std::make_shared<Supercell const>(
      shared_prim, Eigen::Matrix3l::Identity());

  ConfigDoF configdof = make_configdof(*shared_supercell);

  EXPECT_EQ(configdof.has_occupation(), true);
  EXPECT_EQ(configdof.global_dofs().size(), 0);
  EXPECT_EQ(configdof.local_dofs().size(), 1);
  EXPECT_EQ(configdof.local_dofs().count("disp"), 1);
}

TEST(ConfigDoFTest, Constructor5) {
  // construction: shared prim structure with ternary occupation, GLstrain, and
  // disp DoF

  auto shared_prim =
      std::make_shared<Structure const>(test::FCC_ternary_strain_disp_prim());
  auto shared_supercell = std::make_shared<Supercell const>(
      shared_prim, Eigen::Matrix3l::Identity());

  ConfigDoF configdof = make_configdof(*shared_supercell);

  EXPECT_EQ(configdof.has_occupation(), true);
  EXPECT_EQ(configdof.global_dofs().size(), 1);
  EXPECT_EQ(configdof.global_dofs().count("GLstrain"), 1);
  EXPECT_EQ(configdof.local_dofs().size(), 1);
  EXPECT_EQ(configdof.local_dofs().count("disp"), 1);
}

TEST(ConfigDoFTest, Constructor6) {
  auto shared_prim =
      std::make_shared<Structure const>(test::FCC_ternary_strain_disp_prim());
  auto shared_supercell = std::make_shared<Supercell const>(
      shared_prim, Eigen::Matrix3l::Identity());

  ConfigDoF configdof = Configuration::zeros(*shared_supercell).configdof();

  EXPECT_EQ(configdof.has_occupation(), true);
  EXPECT_EQ(configdof.global_dofs().size(), 1);
  EXPECT_EQ(configdof.global_dofs().count("GLstrain"), 1);
  EXPECT_EQ(configdof.local_dofs().size(), 1);
  EXPECT_EQ(configdof.local_dofs().count("disp"), 1);
}
