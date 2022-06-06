#include "Common.hh"
#include "casm/app/CLIParse.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/completer/Handlers.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/BasicStructureIO.hh"
#include "casm/database/ConfigImport.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "gtest/gtest.h"

namespace test {
jsonParser::const_iterator find_mapped(jsonParser const &report,
                                       std::string to_configname) {
  for (auto it = report.begin(); it != report.end(); ++it) {
    auto it_b = it->find("to_configname");
    if (it_b != it->end() && it_b->is_string() &&
        it_b->get<std::string>() == to_configname) {
      return it;
    }
  }
  return report.end();
}
}  // namespace test

// PrimClex without directory for saving project data
class ImportTest : public testing::Test {
 protected:
  void SetUp() override { tmp_dir.do_not_remove_on_destruction(); }

  // copy files to tmp_dir & build the CASM project & construct `primclex`
  void build() {
    // files to copy should be in `data_dir`
    // - expects a "prim.json" file at minimum

    for (auto const &file : data_files) {
      fs::copy_file(data_dir / file, tmp_dir.path() / file);
    }

    jsonParser prim_json{tmp_dir.path() / "prim.json"};
    shared_prim = std::make_shared<Structure const>(read_prim(prim_json, TOL));
    project_settings_ptr = notstd::make_unique<ProjectSettings>(
        make_default_project_settings(*shared_prim, title, tmp_dir.path()));
    build_project(*project_settings_ptr, *shared_prim);
    primclex = notstd::make_unique<PrimClex>(tmp_dir.path());
  }

  // run import
  void import(std::string cli_str, jsonParser const &json_options) {
    // - expects a "import.json" to read JSON options from

    Completer::ImportOption import_opt;
    parse_args(import_opt, cli_str);

    fs::path curr = fs::current_path();
    fs::current_path(tmp_dir.path());
    DB::Import<Configuration>::run(*primclex, json_options, import_opt);
    fs::current_path(curr);
  }

  std::string title;
  fs::path data_dir;
  test::TmpDir tmp_dir;
  std::shared_ptr<Structure const> shared_prim;
  std::unique_ptr<ProjectSettings> project_settings_ptr;
  std::unique_ptr<PrimClex> primclex;
  std::vector<fs::path> data_files;
};

/// This test includes:
/// - "import_properties": false
/// - "copy_structure_files": false
/// - "primitive_only": false
TEST_F(ImportTest, Test1) {
  // ## setup
  title = "ImportTest_test1";
  data_dir = test::data_dir("database") / "import_test1";
  data_files = std::vector<fs::path>(
      {"AB_Ordering_large_supercell.json", "import_list.txt", "import.json",
       "prim.json", "pure_A_large_supercell.json"});
  build();

  for (auto const &file : data_files) {
    EXPECT_TRUE(fs::exists(tmp_dir / file));
  }
  EXPECT_TRUE(fs::exists(tmp_dir.path() / ".casm"));

  // ## pre-condition tests
  EXPECT_EQ(primclex->db<Supercell>().size(), 0);
  EXPECT_EQ(primclex->db<Configuration>().size(), 0);
  EXPECT_EQ(primclex->db_props<Configuration>("default").size(), 0);

  // ## run import
  std::string cli_str =
      "casm import --batch " + (tmp_dir.path() / "import_list.txt").string();
  jsonParser json_options{tmp_dir.path() / "import.json"};
  import(cli_str, json_options);

  // ## post-condition tests

  fs::path report_dir = tmp_dir.path() / "reports" / "import_report.0";
  fs::path map_fail_report = report_dir / "map_fail.json";
  fs::path map_success_report = report_dir / "map_success.json";
  fs::path import_data_fail_report = report_dir / "import_data_fail.json";
  fs::path import_conflict_report = report_dir / "import_conflict.json";

  // should have a "map_sucess.json" file
  EXPECT_TRUE(fs::exists(map_success_report));

  // should not have "map_fail.json", "import_data_fail.json", or
  // "import_conflict.json"
  EXPECT_FALSE(fs::exists(map_fail_report));
  EXPECT_FALSE(fs::exists(import_data_fail_report));
  EXPECT_FALSE(fs::exists(import_conflict_report));

  // should have 3 configurations (original 2, and primitive of pure A)
  jsonParser map_success_json{map_success_report};
  EXPECT_EQ(map_success_json.size(), 3);

  auto it = test::find_mapped(map_success_json, "SCEL54_6_3_3_0_3_3/0");
  EXPECT_TRUE(it != map_success_json.end());

  it = test::find_mapped(map_success_json, "SCEL54_6_3_3_0_3_3/1");
  EXPECT_TRUE(it != map_success_json.end());

  it = test::find_mapped(map_success_json, "SCEL1_1_1_1_0_0_0/0");
  EXPECT_TRUE(it != map_success_json.end());

  EXPECT_EQ(primclex->db<Supercell>().size(), 2);
  EXPECT_EQ(primclex->db<Configuration>().size(), 3);

  // should not import any properties
  EXPECT_EQ(primclex->db_props<Configuration>("default").size(), 0);
}

/// This test includes:
/// - "import_properties": true
/// - "copy_structure_files":true
/// - "primitive_only": false
TEST_F(ImportTest, Test2) {
  // ## setup
  title = "ImportTest_test2";
  data_dir = test::data_dir("database") / "import_test2";
  data_files = std::vector<fs::path>(
      {"AB_Ordering_large_supercell.json", "import_list.txt", "import.json",
       "prim.json", "pure_A_large_supercell.json"});
  build();

  for (auto const &file : data_files) {
    EXPECT_TRUE(fs::exists(tmp_dir / file));
  }
  EXPECT_TRUE(fs::exists(tmp_dir.path() / ".casm"));

  // ## pre-condition tests
  EXPECT_EQ(primclex->db<Supercell>().size(), 0);
  EXPECT_EQ(primclex->db<Configuration>().size(), 0);
  EXPECT_EQ(primclex->db_props<Configuration>("default").size(), 0);

  // ## run import
  std::string cli_str =
      "casm import --batch " + (tmp_dir.path() / "import_list.txt").string();
  jsonParser json_options{tmp_dir.path() / "import.json"};
  import(cli_str, json_options);

  // ## post-condition tests

  fs::path report_dir = tmp_dir.path() / "reports" / "import_report.0";
  fs::path map_fail_report = report_dir / "map_fail.json";
  fs::path map_success_report = report_dir / "map_success.json";
  fs::path import_data_fail_report = report_dir / "import_data_fail.json";
  fs::path import_conflict_report = report_dir / "import_conflict.json";

  // should have a "map_sucess.json" file
  EXPECT_TRUE(fs::exists(map_success_report));

  // should not have "map_fail.json", "import_data_fail.json", or
  // "import_conflict.json"
  EXPECT_FALSE(fs::exists(map_fail_report));
  EXPECT_FALSE(fs::exists(import_data_fail_report));
  EXPECT_FALSE(fs::exists(import_conflict_report));

  // should have 3 configurations (original 2, and primitive of pure A)
  jsonParser map_success_json{map_success_report};
  EXPECT_EQ(map_success_json.size(), 3);

  auto it = test::find_mapped(map_success_json, "SCEL54_6_3_3_0_3_3/0");
  EXPECT_TRUE(it != map_success_json.end());

  it = test::find_mapped(map_success_json, "SCEL54_6_3_3_0_3_3/1");
  EXPECT_TRUE(it != map_success_json.end());

  it = test::find_mapped(map_success_json, "SCEL1_1_1_1_0_0_0/0");
  EXPECT_TRUE(it != map_success_json.end());

  EXPECT_EQ(primclex->db<Supercell>().size(), 2);
  EXPECT_EQ(primclex->db<Configuration>().size(), 3);

  // should import properties for original size configurations only
  EXPECT_EQ(primclex->db_props<Configuration>("default").size(), 2);
}

/// This test includes:
/// - "import_properties": true
/// - "copy_structure_files":true
/// - "primitive_only": true
TEST_F(ImportTest, Test3) {
  // ## setup
  title = "ImportTest_test3";
  data_dir = test::data_dir("database") / "import_test3";
  data_files = std::vector<fs::path>(
      {"AB_Ordering_large_supercell.json", "import_list.txt", "import.json",
       "prim.json", "pure_A_large_supercell.json"});
  build();

  for (auto const &file : data_files) {
    EXPECT_TRUE(fs::exists(tmp_dir / file));
  }
  EXPECT_TRUE(fs::exists(tmp_dir.path() / ".casm"));

  // ## pre-condition tests
  EXPECT_EQ(primclex->db<Supercell>().size(), 0);
  EXPECT_EQ(primclex->db<Configuration>().size(), 0);
  EXPECT_EQ(primclex->db_props<Configuration>("default").size(), 0);

  // ## run import
  std::string cli_str =
      "casm import --batch " + (tmp_dir.path() / "import_list.txt").string();
  jsonParser json_options{tmp_dir.path() / "import.json"};
  import(cli_str, json_options);

  // ## post-condition tests

  fs::path report_dir = tmp_dir.path() / "reports" / "import_report.0";
  fs::path map_fail_report = report_dir / "map_fail.json";
  fs::path map_success_report = report_dir / "map_success.json";
  fs::path import_data_fail_report = report_dir / "import_data_fail.json";
  fs::path import_conflict_report = report_dir / "import_conflict.json";

  // should have a "map_sucess.json" file
  EXPECT_TRUE(fs::exists(map_success_report));

  // should not have "map_fail.json", "import_data_fail.json", or
  // "import_conflict.json"
  EXPECT_FALSE(fs::exists(map_fail_report));
  EXPECT_FALSE(fs::exists(import_data_fail_report));
  EXPECT_FALSE(fs::exists(import_conflict_report));

  // should have 2 configurations (primitive ordered and primitive of pure A)
  jsonParser map_success_json{map_success_report};
  EXPECT_EQ(map_success_json.size(), 2);

  auto it = test::find_mapped(map_success_json, "SCEL54_6_3_3_0_3_3/0");
  EXPECT_TRUE(it != map_success_json.end());

  it = test::find_mapped(map_success_json, "SCEL1_1_1_1_0_0_0/0");
  EXPECT_TRUE(it != map_success_json.end());

  EXPECT_EQ(primclex->db<Supercell>().size(), 2);
  EXPECT_EQ(primclex->db<Configuration>().size(), 2);

  // should import properties for original size configurations only
  EXPECT_EQ(primclex->db_props<Configuration>("default").size(), 1);
}
