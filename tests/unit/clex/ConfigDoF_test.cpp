#include "gtest/gtest.h"
#include "crystallography/TestStructures.hh"

#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"

using namespace CASM;

TEST(ConfigDoFTest, TestEmpty) {

  auto shared_prim = std::make_shared<Structure const>(test::empty_prim());
  auto shared_supercell = std::make_shared<Supercell const>(shared_prim, Eigen::Matrix3l::Identity());

  ConfigDoF configdof = shared_supercell->zero_configdof(shared_prim->lattice().tol());

  EXPECT_EQ(configdof.has_occupation(), false);
  EXPECT_EQ(configdof.global_dofs().size(), 0);
  EXPECT_EQ(configdof.local_dofs().size(), 0);
}
