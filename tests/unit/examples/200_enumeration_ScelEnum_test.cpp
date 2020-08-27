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
class ExampleEnumerationZrOScelEnum : public testing::Test {
protected:

  std::string title;
  std::shared_ptr<CASM::Structure const> shared_prim;
  CASM::fs::path root_dir;
  CASM::ProjectSettings project_settings;
  CASM::PrimClex primclex;

  int begin_volume;
  int end_volume;
  std::string dirs;
  Eigen::Matrix3i generating_matrix;
  CASM::xtal::ScelEnumProps enumeration_params;

  ExampleEnumerationZrOScelEnum():
    title("ExampleEnumerationZrOScelEnum"),
    shared_prim(std::make_shared<CASM::Structure const>(test::ZrO_prim())),
    root_dir(test::proj_dir(autotools::abs_srcdir() + "/tests/unit/test_projects/ZrOScelEnum")),
    project_settings(make_default_project_settings(*shared_prim, title, root_dir)),
    primclex(project_settings, shared_prim),
    begin_volume(1),
    end_volume(5),
    dirs("abc"),
    generating_matrix(Eigen::Matrix3i::Identity()),
    enumeration_params(begin_volume, end_volume, dirs, generating_matrix) {}

};


TEST_F(ExampleEnumerationZrOScelEnum, Example1) {

  bool existing_only = false;

  CASM::ScelEnumByProps enumerator {primclex, enumeration_params, existing_only};

  std::vector<CASM::Supercell> supercells {enumerator.begin(), enumerator.end()};

  EXPECT_EQ(supercells.size(), 20);
}
