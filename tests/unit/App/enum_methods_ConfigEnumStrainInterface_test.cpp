#include "App/TestEnumeratorInterface.hh"
#include "Common.hh"
#include "autotools.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/enum.hh"
#include "casm/app/enum/methods/ConfigEnumStrainInterface.hh"
#include "casm/app/enum/methods/ScelEnumInterface.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/ConfigEnumStrain.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/completer/Handlers.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/database/ConfigDatabaseTools_impl.hh"
#include "casm/database/Database.hh"
#include "casm/database/ScelDatabaseTools_impl.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

// This test fixture class constructs a CASM project for enumeration examples
class enum_methods_ConfigEnumStrainInterfaceTest : public testing::Test {
 protected:
  std::string title;
  test::TmpDir tmp_dir;
  fs::path proj_dir;
  std::shared_ptr<CASM::Structure const> shared_prim;
  CASM::ProjectSettings project_settings;
  bool build_project_placeholder;
  CASM::PrimClex primclex;

  enum_methods_ConfigEnumStrainInterfaceTest();
};

enum_methods_ConfigEnumStrainInterfaceTest::
    enum_methods_ConfigEnumStrainInterfaceTest()
    : title("ConfigEnumStrainInterfaceTest"),
      tmp_dir(),
      proj_dir(tmp_dir.path()),
      shared_prim(std::make_shared<CASM::Structure const>(
          test::FCC_ternary_GLstrain_disp_prim())),
      project_settings(
          make_default_project_settings(*shared_prim, title, proj_dir)),
      build_project_placeholder((build_project(project_settings, *shared_prim),
                                 commit(project_settings), true)),
      primclex(project_settings, shared_prim) {
  EXPECT_EQ(primclex.prim().basis().size(), 1);

  // Enumerate supercells
  int begin_volume{1};
  int end_volume{5};
  std::string dirs{"abc"};
  Eigen::Matrix3i generating_matrix{Eigen::Matrix3i::Identity()};
  CASM::xtal::ScelEnumProps enumeration_params{begin_volume, end_volume, dirs,
                                               generating_matrix};

  CASM::ScelEnumByProps enumerator{shared_prim, enumeration_params};
  for (Supercell const &supercell : enumerator) {
    /// Use while transitioning Supercell to no longer need a `PrimClex const *`
    supercell.set_primclex(&primclex);

    make_canonical_and_insert(enumerator, supercell, primclex.db<Supercell>());
  }
  primclex.db<Supercell>().commit();
  EXPECT_EQ(primclex.db<Supercell>().size(), 13);
  primclex.db<Supercell>().close();
  primclex.db<Supercell>().open();
  EXPECT_EQ(primclex.db<Supercell>().size(), 13);
}

TEST_F(enum_methods_ConfigEnumStrainInterfaceTest, Test1) {
  // ScopedNullLogging logging;
  CASM::log().set_verbosity(Log::debug);

  {
    EXPECT_EQ(primclex.db<Configuration>().size(), 0);
    fs::path test_dir = test::proj_dir(primclex.dir().root_dir() /
                                       "ConfigEnumStrainInterfaceTest");
    std::string cli_str = "casm enum --method ConfigEnumStrain";
    jsonParser json_options;
    json_options["max"] = 0.11;
    json_options["increment"] = 0.1;
    json_options["trim_corners"] = false;
    json_options["scelnames"] = std::vector<std::string>{"SCEL1_1_1_1_0_0_0"};
    json_options["output_configurations"] = true;
    json_options["output_dir"] = test_dir.string();
    test::run_enum_interface<ConfigEnumStrainInterface>(cli_str, primclex,
                                                        json_options);
    EXPECT_EQ(primclex.db<Configuration>().size(), 20);
  }
}
