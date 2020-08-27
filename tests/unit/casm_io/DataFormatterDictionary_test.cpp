#include "gtest/gtest.h"

/// What is being tested:
#include "casm/casm_io/dataformatter/DataFormatter.hh"

#include <sstream>
#include "autotools.hh"
#include "Common.hh"
#include "casm/casm_io/json/jsonFile.hh"
#include "casm/database/Selection.hh"
#include "casm/database/Selected.hh"
#include "casm/database/Database.hh"

// #include "FCCTernaryProj.hh"

// #include "casm/clex/PrimClex.hh"
// #include "casm/clex/ScelEnum.hh"
// #include "casm/clex/ConfigEnumAllOccupations_impl.hh"
// #include "casm/clex/ConfigIO.hh"
// #include "casm/app/ProjectSettings.hh"
// #include "casm/app/QueryHandler.hh"
// #include "casm/app/ProjectSettings.hh"

// #include "casm/enumerator/Enumerator_impl.hh"
// #include "casm/app/enum.hh"


using namespace CASM;

// TEST(DataFormatterDictionaryTest, Test0) {
//   DataFormatterDictionary<Configuration> dict;
//   dict.insert(ConfigIO::configname(), ConfigIO::scel_size(), alias_or_name<Configuration>());
//
//   test::FCCTernaryProj proj;
//   proj.check_init();
//
//   PrimClex primclex(proj.dir, null_log());
//   const Structure &prim(primclex.prim());
//   primclex.settings().set_crystallography_tol(1e-5);
//
//   Eigen::Vector3d a, b, c;
//   std::tie(a, b, c) = prim.lattice().vectors();
//
//   // format one config
//   {
//     Supercell scel(&primclex, Lattice(a, b, c));
//     std::vector<Configuration> config_vec;
//     for(Index i = 0; i < 3; ++i) {
//       config_vec.emplace_back(scel);
//       config_vec.back().init_occupation();
//     }
//
//     DataFormatter<Configuration> formatter {dict.parse("configname scel_size alias_or_name")};
//
//     {
//       std::stringstream ss;
//       ss << formatter(config_vec[0]);
//
//       std::string result =
//         "#               configname    scel_size             alias_or_name\n"
//         "    SCEL1_1_1_1_0_0_0/none            1    SCEL1_1_1_1_0_0_0/none\n";
//       ASSERT_EQ(ss.str(), result);
//     }
//
//     {
//       std::stringstream ss;
//       ss << formatter(config_vec.begin(), config_vec.end());
//
//       std::string result =
//         "#               configname    scel_size             alias_or_name\n"
//         "    SCEL1_1_1_1_0_0_0/none            1    SCEL1_1_1_1_0_0_0/none\n"
//         "    SCEL1_1_1_1_0_0_0/none            1    SCEL1_1_1_1_0_0_0/none\n"
//         "    SCEL1_1_1_1_0_0_0/none            1    SCEL1_1_1_1_0_0_0/none\n";
//       ASSERT_EQ(ss.str(), result);
//     }
//
//     {
//       jsonParser json;
//       json = formatter(config_vec[0]);
//
//       jsonParser result;
//       result["configname"] = "SCEL1_1_1_1_0_0_0/none";
//       result["scel_size"] = 1;
//       result["alias_or_name"] = "SCEL1_1_1_1_0_0_0/none";
//       ASSERT_EQ(json, result);
//     }
//
//     {
//       jsonParser json;
//       json = formatter(config_vec.begin(), config_vec.end());
//
//       jsonParser result = jsonParser::array(3, jsonParser::object());
//       for(int i = 0; i < 3; ++i) {
//         result[i]["configname"] = "SCEL1_1_1_1_0_0_0/none";
//         result[i]["scel_size"] = 1;
//         result[i]["alias_or_name"] = "SCEL1_1_1_1_0_0_0/none";
//       }
//       ASSERT_EQ(json, result);
//     }
//   }
// }
//
// TEST(DataFormatterDictionaryTest, Test1) {
//
//   ASSERT_EQ(true, true);
//   test::FCCTernaryProj proj;
//   proj.check_init();
//   proj.check_composition();
//
//   PrimClex primclex(proj.dir, null_log());
//   const Structure &prim(primclex.prim());
//   primclex.settings().set_crystallography_tol(1e-5);
//
//   Completer::EnumOption enum_opt;
//   enum_opt.desc();
//
//   Eigen::Vector3d a, b, c;
//   std::tie(a, b, c) = prim.lattice().vectors();
//
//   // -- Generate Supercell & Configuration --
//
//   ScelEnumByProps enum_scel(primclex, ScelEnumProps(1, 5));
//   ASSERT_EQ(true, true);
//
//   ConfigEnumAllOccupations::run(primclex, enum_scel.begin(), enum_scel.end());
//   ASSERT_EQ(true, true);
//
//   // -- Check Supercell --
//
//   {
//     ASSERT_EQ(primclex.generic_db<Supercell>().size(), 13);
//     DB::Selection<Supercell> selection(primclex, "ALL");
//
//     auto &qh = primclex.settings().query_handler<Supercell>();
//     qh.set_selected(selection);
//     auto &dict = qh.dict();
//
//     std::vector<std::string> all_col = {"name", "selected", "scel_size"};
//     auto formatter = dict.parse(all_col);
//
//     std::stringstream ss;
//     ss << formatter(selection.all().begin(), selection.all().end()) << std::endl;
//     ASSERT_EQ(true, true);
//   }
//
//   // -- Check Configuration --
//
//   {
//     ASSERT_EQ(primclex.generic_db<Configuration>().size(), 126);
//
//     DB::Selection<Configuration> selection(primclex, "ALL");
//
//     auto &qh = primclex.settings().query_handler<Configuration>();
//     qh.set_selected(selection);
//     auto &dict = qh.dict();
//
//     std::vector<std::string> all_col = {"name", "selected", "comp"};
//     auto formatter = dict.parse(all_col);
//
//     std::stringstream ss;
//     ss << formatter(selection.all().begin(), selection.all().end()) << std::endl;
//     ASSERT_EQ(true, true);
//   }
// }
