#include "gtest/gtest.h"
#include "autotools.hh"
#include "Common.hh"

#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Superlattice.hh"

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
class ExampleEnumerationZrOConfigEnumAllOccupations : public testing::Test {
protected:

  std::string title;
  std::shared_ptr<CASM::Structure const> shared_prim;
  CASM::fs::path root_dir;
  CASM::ProjectSettings project_settings;
  CASM::PrimClex primclex;

  ExampleEnumerationZrOConfigEnumAllOccupations():
    title("ExampleEnumerationZrOConfigEnumAllOccupations"),
    shared_prim(std::make_shared<CASM::Structure const>(test::ZrO_prim())),
    root_dir(test::proj_dir(autotools::abs_srcdir() + "/tests/unit/test_projects/ZrOConfigEnumAllOccupations")),
    project_settings(make_default_project_settings(*shared_prim, title, root_dir)),
    primclex(project_settings, shared_prim) {

    int begin_volume {1};
    int end_volume {5};
    std::string dirs {"abc"};
    Eigen::Matrix3i generating_matrix {Eigen::Matrix3i::Identity()};
    CASM::xtal::ScelEnumProps enumeration_params {begin_volume, end_volume, dirs, generating_matrix};

    CASM::ScelEnumByProps enumerator {primclex, enumeration_params, existing_only};

    // Currently, when CASM::ScelEnumByProps enumerates a Supercell it is automatically added
    //  to the supercell database available at `primclex.db<Supercell>()`.
    int count = std::distance(enumerator.begin(), enumerator.end());
    EXPECT_EQ(count, 20);
  }

};


TEST_F(ExampleEnumerationZrOConfigEnumAllOccupations, Example1) {

  std::vector<CASM::Configuration> configurations;
  for(auto const &scel : primclex.db<Supercell>()) {
    CASM::ConfigEnumInput input {scel};
    CASM::ConfigEnumAllOccupations enumerator {input};
    configurations.push_back(enumerator.begin(), enumerator.end());
  }
  EXPECT_EQ(configurations.size(), 336);
}
