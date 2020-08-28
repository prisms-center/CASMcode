#include "gtest/gtest.h"
#include "autotools.hh"
#include "Common.hh"

#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/clex/ConfigEnumStrain.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Superlattice.hh"
#include "casm/symmetry/SymRepTools.hh"

#include "crystallography/TestStructures.hh" // for test::ZrO_prim

// This test fixture class constructs a CASM project for enumeration examples
//
// Note:
// - Currently, in order to enumerate Supercells, Configurations, etc. a CASM project must be given
//   a project directory, even if they are never written to disk. In the future this requirement
//   will be removed.
// - To create a project directory for example and testing purposes use the `test::proj_dir`
//   function as shown below. It will automatically create a new directory by appending ".N" to the
//   path argument, where N is an incremented integer to ensure a unique new directory.
//
class ExampleEnumerationSimpleCubicConfigEnumStrain : public testing::Test {
protected:

  std::string title;
  std::shared_ptr<CASM::Structure const> shared_prim;
  CASM::fs::path root_dir;
  CASM::ProjectSettings project_settings;
  CASM::PrimClex primclex;

  ExampleEnumerationSimpleCubicConfigEnumStrain():
    title("ExampleEnumerationSimpleCubicConfigEnumStrain"),
    shared_prim(std::make_shared<CASM::Structure const>(test::SimpleCubicGLstrain())),
    root_dir(test::proj_dir(autotools::abs_srcdir() + "/tests/unit/test_projects/SimpleCubicConfigEnumStrain")),
    project_settings(make_default_project_settings(*shared_prim, title, root_dir)),
    primclex(project_settings, shared_prim) {

    int begin_volume {1};
    int end_volume {2};
    std::string dirs {"abc"};
    Eigen::Matrix3i generating_matrix {Eigen::Matrix3i::Identity()};
    CASM::xtal::ScelEnumProps enumeration_params {begin_volume, end_volume, dirs, generating_matrix};
    bool existing_only = false;
    CASM::ScelEnumByProps enumerator {primclex, enumeration_params, existing_only};

    // Currently, when CASM::ScelEnumByProps enumerates a Supercell it is automatically added
    //  to the supercell database available at `primclex.db<Supercell>()`.
    int count = std::distance(enumerator.begin(), enumerator.end());
    EXPECT_EQ(count, 1);
  }

};

namespace enumeration_test_impl {

  SymGroupRep const &global_dof_symrep(Structure const &prim, DoFKey global_dof_key) {
    return prim.factor_group().representation(prim.global_dof_symrep_ID(global_dof_key));
  }

  SymGroupRep const &global_dof_symrep(ConfigEnumInput const &config_input, DoFKey global_dof_key) {
    return global_dof_symrep(config_input.supercell().prim(), global_dof_key);
  }

  SymGroup make_point_group(ConfigEnumInput const &config_input) {
    return make_point_group(
             config_input.group(),
             config_input.supercell().sym_info().supercell_lattice());
  }
}

TEST_F(ExampleEnumerationSimpleCubicConfigEnumStrain, Example1) {

  using namespace enumeration_test_impl; // for global_dof_symrep, make_point_group

  std::vector<CASM::Configuration> configurations;
  for(auto const &scel : primclex.db<Supercell>()) {

    // Inputs needed for CASM::ConfigEnumStrain enumerator
    CASM::ConfigEnumInput config_input {scel};
    std::vector<SymRepTools::SubWedge> wedges;
    Eigen::VectorXd min_value;
    Eigen::VectorXd max_value;
    Eigen::VectorXd increment_value;
    bool auto_range = false;
    bool trim_corners = true;
    DoFKey strain_key = "GLstrain";

    // Make wedges
    std::cout << "make wedges" << std::endl;
    SymGroup point_group = make_point_group(config_input);
    SymGroupRep const &strain_symrep = global_dof_symrep(config_input, strain_key);
    wedges = SymRepTools::symrep_subwedges(strain_symrep, point_group);

    // Make counter parameters
    std::cout << "make counter parameters" << std::endl;
    // For this case, use the full space, so subspace rank == symrep dimensionality
    int dim = strain_symrep.dim();
    Index rank = dim;
    min_value = Eigen::VectorXd::Constant(rank, 0.0);
    max_value = Eigen::VectorXd::Constant(rank, 1.01);
    increment_value = Eigen::VectorXd::Constant(rank, 0.5);

    std::cout << "make ConfigEnumStrain enumerator" << std::endl;
    CASM::ConfigEnumStrain enumerator {config_input, wedges, min_value, max_value, increment_value,
                                       strain_key, auto_range, trim_corners};

    std::cout << "begin enumeration" << std::endl;
    std::copy(enumerator.begin(), enumerator.end(), std::back_inserter(configurations));
    std::cout << "enumeration complete" << std::endl;
  }
  EXPECT_EQ(configurations.size(), pow(3, 6));
}
