#include "gtest/gtest.h"
#include "autotools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/DoFSpace_impl.hh"
#include "casm/enumerator/io/json/DoFSpace.hh"
#include "crystallography/TestStructures.hh"

using namespace CASM;
using namespace test;

namespace {
  Eigen::Matrix3l _fcc_conventional_transf_mat() {
    Eigen::Matrix3l transf_mat;
    transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    return transf_mat;
  }
}

class DoFSpaceTest : public testing::Test {
protected:

  std::shared_ptr<CASM::Structure const> shared_prim;
  std::shared_ptr<CASM::Supercell> shared_supercell;  // conventional unit cell

  DoFSpaceTest():
    shared_prim(std::make_shared<CASM::Structure const>(test::FCC_ternary_strain_disp_prim())),
    shared_supercell(std::make_shared<CASM::Supercell>(shared_prim, _fcc_conventional_transf_mat())) {}

};

TEST_F(DoFSpaceTest, ConstructorTest1) {

  // Construct the GLstrain DoF space.
  std::cout << "begin DoFSpaceTest" << std::endl;
  ConfigEnumInput config_input {*shared_supercell};
  std::cout << "DoFSpaceTest 0" << std::endl;
  DoFKey dof_key = "GLstrain";
  std::cout << "DoFSpaceTest 1" << std::endl;
  get_dof_space_dimension(dof_key, config_input.configuration(), config_input.sites());
  std::cout << "DoFSpaceTest 2" << std::endl;
  DoFSpace dof_space {config_input, dof_key};
  std::cout << "DoFSpaceTest 3" << std::endl;

  // Uncomment to print dof_space:
  jsonParser dof_space_json;
  to_json(dof_space, dof_space_json);
  std::cout << "DoFSpace:\n" << dof_space_json << std::endl;

  std::cout << "end DoFSpaceTest" << std::endl;
}
