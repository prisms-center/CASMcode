#include "gtest/gtest.h"

#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/Superlattice.hh"

#include "crystallography/TestStructures.hh" // for test::ZrO_prim

class ExampleEnumerationZrOScelEnum : public testing::Test {
protected:

  std::string title;
  std::shared_ptr<CASM::Structure const> shared_prim;
  CASM::ProjectSettings project_settings;
  CASM::PrimClex primclex;

  int begin_volume;
  int end_volume;
  std::string dirs;
  Eigen::Matrix3i generating_matrix;
  CASM::xtal::ScelEnumProps enumeration_params;

  ExampleEnumerationScelEnum():
    title("ExampleZrOScelEnum"),
    shared_prim(std::make_shared<CASM::Structure const>(test::ZrO_prim())),
    project_settings(make_default_project_settings(*shared_prim, title)),
    primclex(project_settings, shared_prim),
    begin_volume(1),
    end_volume(5),
    dirs("abc"),
    generating_matrix(Eigen::Matrix3i::Identity()),
    enumeration_params(begin_volume, end_volume, dirs, generating_matrix) {}

};

TEST_F(ExampleEnumerationScelEnum, ScelEnumByPropsExample) {

  bool existing_only = false;

  CASM::ScelEnumByProps enumerator {primclex, enumeration_params, existing_only};

  std::vector<CASM::Supercell> supercells {enumerator.begin(), enumerator.end()};

  EXPECT_EQ(supercells.size(), 16);
}
