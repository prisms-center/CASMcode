#include "gtest/gtest.h"
#include "autotools.hh"

/// What is being tested:
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/ConfigEnumAllOccupations.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "casm/enumerator/io/json/ConfigEnumInput_json_io.hh"
#include "casm/crystallography/Structure.hh"
#include "crystallography/TestStructures.hh"

using namespace CASM;
using namespace test;

namespace {
  Eigen::Matrix3i _fcc_conventional_transf_mat() {
    Eigen::Matrix3i transf_mat;
    transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
    return transf_mat;
  }

  fs::path enumerator_test_directory() {
    return autotools::abs_srcdir() + "/tests/unit/enumerator";
  }
}

class ConfigEnumInputjsonIOTest : public testing::Test {
protected:

  std::string title;
  std::shared_ptr<CASM::Structure const> shared_prim;
  CASM::ProjectSettings project_settings;
  CASM::PrimClex primclex;

  ConfigEnumInputjsonIOTest():
    title("ConfigEnumInputjsonIOTest"),
    shared_prim(std::make_shared<CASM::Structure const>(test::FCC_ternary_prim())),
    project_settings(make_default_project_settings(*shared_prim, title)),
    primclex(project_settings, shared_prim) {

    int begin_volume {1};
    int end_volume {5};
    std::string dirs {"abc"};
    Eigen::Matrix3i generating_matrix {Eigen::Matrix3i::Identity()};
    CASM::xtal::ScelEnumProps enumeration_params {begin_volume, end_volume, dirs, generating_matrix};
    bool existing_only = false;

    // The ScelEnumByProps variant that accepts a PrimClex in the constructor inserts Supercells into
    // the Supercell database available at `primclex.db<Supercell>()` as it constructs them.
    CASM::ScelEnumByProps enumerator {primclex, enumeration_params, existing_only};

    // Increments the enumerator iterators to construct all Supercell
    int count = std::distance(enumerator.begin(), enumerator.end());
    EXPECT_EQ(count, 13);
    EXPECT_EQ(primclex.db<Supercell>().size(), 13);

    // Also add Configurations
    for(auto const &supercell : primclex.db<Supercell>()) {
      CASM::ConfigEnumAllOccupations enumerator {supercell};
      for(auto const &configuration : enumerator) {
        primclex.db<Configuration>().insert(configuration);
      }
    }

    EXPECT_EQ(primclex.db<Configuration>().size(), 126);
  }

};

// TEST_F(ConfigEnumInputjsonIOTest, Test0) {
//   for(Supercell const &supercell : primclex.db<Supercell>()) {
//     std::cout << supercell.name() << std::endl;
//   }
//   for(auto const &configuration: primclex.db<Configuration>()) {
//     std::cout << configuration.name() << std::endl;
//   }
//   EXPECT_EQ(1, 1);
// }

TEST_F(ConfigEnumInputjsonIOTest, Test1) {
  jsonParser json = CASM::jsonParser::parse(std::string(R"({"supercells": { "min":1, "max":4 } })"));

  InputParser<std::vector<ConfigEnumInput>> parser {json,
                                                    shared_prim,
                                                    primclex.db<Supercell>(),
                                                    primclex.db<Configuration>()};

  std::runtime_error error_if_invalid {"Failed to parse ConfigEnumInput JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  EXPECT_EQ(parser.value->size(), 13);
}

TEST_F(ConfigEnumInputjsonIOTest, Test2) {
  jsonParser json = CASM::jsonParser::parse(std::string(R"({
    "scelnames": [
      "SCEL1_1_1_1_0_0_0",
      "SCEL2_2_1_1_0_1_1",
      "SCEL2_2_1_1_0_0_1",
      "SCEL3_3_1_1_0_2_2",
      "SCEL3_3_1_1_0_0_2",
      "SCEL3_3_1_1_0_2_1"
    ]
  })"));

  InputParser<std::vector<ConfigEnumInput>> parser {json,
                                                    shared_prim,
                                                    primclex.db<Supercell>(),
                                                    primclex.db<Configuration>()};

  std::runtime_error error_if_invalid {"Failed to parse ConfigEnumInput JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  EXPECT_EQ(parser.value->size(), 6);
}

TEST_F(ConfigEnumInputjsonIOTest, Test3) {

  DB::Selection<Supercell> supercell_selection {primclex.db<Supercell>(), "NONE"};
  for(Supercell const &supercell : primclex.db<Supercell>()) {
    if(supercell.volume() == 4) {
      supercell_selection.data()[supercell.name()] = true;
    }
  }
  fs::path out_path = enumerator_test_directory() / "test3_selection.json";
  supercell_selection.write(make_dictionary<Supercell>(),
                            true, // force
                            out_path,
                            true, // write_json
                            true); // only_selected
  EXPECT_EQ(fs::exists(out_path), true);

  jsonParser json = CASM::jsonParser::parse(std::string(R"({
    "scelnames": [
      "SCEL1_1_1_1_0_0_0",
      "SCEL2_2_1_1_0_1_1",
      "SCEL2_2_1_1_0_0_1",
      "SCEL3_3_1_1_0_2_2",
      "SCEL3_3_1_1_0_0_2",
      "SCEL3_3_1_1_0_2_1"
    ]
  })"));
  json["supercell_selection"] = out_path.string();

  InputParser<std::vector<ConfigEnumInput>> parser {json,
                                                    shared_prim,
                                                    primclex.db<Supercell>(),
                                                    primclex.db<Configuration>()};

  std::runtime_error error_if_invalid {"Failed to parse ConfigEnumInput JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  EXPECT_EQ(parser.value->size(), 13);
}

TEST_F(ConfigEnumInputjsonIOTest, Test4) {
  jsonParser json = CASM::jsonParser::parse(std::string(R"({
    "confignames": [
      "SCEL1_1_1_1_0_0_0/0",
      "SCEL1_1_1_1_0_0_0/1",
      "SCEL1_1_1_1_0_0_0/2"
    ]
  })"));

  InputParser<std::vector<ConfigEnumInput>> parser {json,
                                                    shared_prim,
                                                    primclex.db<Supercell>(),
                                                    primclex.db<Configuration>()};

  std::runtime_error error_if_invalid {"Failed to parse ConfigEnumInput JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  EXPECT_EQ(parser.value->size(), 3);
}

TEST_F(ConfigEnumInputjsonIOTest, Test5) {
  DB::Selection<Configuration> configuration_selection {primclex.db<Configuration>(), "NONE"};
  for(Configuration const &configuration : primclex.db<Configuration>()) {
    if(configuration.supercell().volume() > 1) {
      configuration_selection.data()[configuration.name()] = true;
    }
  }
  fs::path out_path = enumerator_test_directory() / "test5_selection.json";
  configuration_selection.write(make_dictionary<Configuration>(),
                                true, // force
                                out_path,
                                true, // write_json
                                true); // only_selected
  EXPECT_EQ(fs::exists(out_path), true);

  jsonParser json = CASM::jsonParser::parse(std::string(R"({
    "confignames": [
      "SCEL1_1_1_1_0_0_0/0",
      "SCEL1_1_1_1_0_0_0/1",
      "SCEL1_1_1_1_0_0_0/2"
    ]
  })"));
  json["config_selection"] = out_path.string();

  InputParser<std::vector<ConfigEnumInput>> parser {json,
                                                    shared_prim,
                                                    primclex.db<Supercell>(),
                                                    primclex.db<Configuration>()};

  std::runtime_error error_if_invalid {"Failed to parse ConfigEnumInput JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  EXPECT_EQ(parser.value->size(), 126);
}

TEST_F(ConfigEnumInputjsonIOTest, Test6) {
  jsonParser json = CASM::jsonParser::parse(std::string(R"({
    "supercells": {
      "min":1,
      "max":4
    },
    "sites": [
      [0, 0, 0, 0],
      [0, 1, 0, 0]
    ]
  })"));

  std::vector<UnitCellCoord> sites {
    UnitCellCoord {0, UnitCell {0, 0, 0} },
    UnitCellCoord {0, UnitCell {1, 0, 0} }
  };

  InputParser<std::vector<ConfigEnumInput>> parser {json,
                                                    shared_prim,
                                                    primclex.db<Supercell>(),
                                                    primclex.db<Configuration>()};

  std::runtime_error error_if_invalid {"Failed to parse ConfigEnumInput JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  EXPECT_EQ(parser.value->size(), 13);
  for(ConfigEnumInput const &input : *parser.value) {
    EXPECT_EQ(input.sites().count(input.configuration().linear_index(sites[0])), 1);
    EXPECT_EQ(input.sites().count(input.configuration().linear_index(sites[1])), 1);
  }
}

namespace {

  std::string _pos_string(Configuration const &_config, std::string name = "none") {
    std::stringstream ss;
    VaspIO::PrintPOSCAR p(make_simple_structure(_config.supercell(), _config.configdof()), name);
    p.sort();
    p.print(ss);
    return ss.str();
  }
}

TEST_F(ConfigEnumInputjsonIOTest, Test7) {

  std::cout << "here 0" << std::endl;
  auto conventional_fcc = std::make_shared<Supercell>(shared_prim, _fcc_conventional_transf_mat());
  std::cout << "here 1" << std::endl;
  auto conventional_fcc_3x3x3 = std::make_shared<Supercell>(shared_prim, 3 * _fcc_conventional_transf_mat());
  std::cout << "here 2" << std::endl;
  Configuration config_L12 {conventional_fcc};
  std::cout << "here 3" << std::endl;
  config_L12.configdof().occ(0) = 1;
  std::cout << "here 4" << std::endl;
  FillSupercell filler {conventional_fcc_3x3x3, config_L12, TOL};
  std::cout << "here 5" << std::endl;
  Configuration config_L12_3x3x3 = filler(config_L12);

  std::cout << "here 6" << std::endl;
  std::cout << _pos_string(config_L12) << std::endl;
  std::cout << _pos_string(config_L12_3x3x3) << std::endl;

  EXPECT_EQ(1, 1);
}
