#include "casm/clex/ConfigEnumAllOccupations.hh"

#include "casm/clex/ScelEnum.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/enumerator/ConfigEnumInput.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

namespace test {

struct CompareVec {
  bool operator()(Eigen::VectorXi const &A, Eigen::VectorXi const &B) const {
    return std::lexicographical_compare(A.data(), A.data() + A.size(), B.data(),
                                        B.data() + B.size());
  }
};

// This stores a histogram of composition (number of each type of molecule) -> #
// of occurances
typedef std::map<Eigen::VectorXi, int, CompareVec> OccupationHistogram;

// Add a configuration to the histogram
void insert(Configuration const &config, OccupationHistogram &occ_histogram) {
  auto num_each_molecule = config.num_each_molecule();
  auto it = occ_histogram.find(num_each_molecule);
  if (it == occ_histogram.end()) {
    occ_histogram.emplace(num_each_molecule, 1);
  } else {
    it->second++;
  }
}

template <typename ConfigIterator>
OccupationHistogram make_occ_histogram(ConfigIterator begin,
                                       ConfigIterator end) {
  OccupationHistogram occ_histogram;
  for (auto it = begin; it != end; ++it) {
    insert(*it, occ_histogram);
  }
  return occ_histogram;
}

std::ostream &print_histogram(std::ostream &sout, Eigen::VectorXi const &vector,
                              std::vector<std::string> const &molecule_name) {
  for (int i = 0; i < vector.size(); ++i) {
    sout << molecule_name[i] << ": " << vector[i] << " ";
  }
  return sout;
}

Eigen::VectorXi to_VectorXi(std::vector<int> const &in) {
  Eigen::VectorXi out;
  out.resize(in.size());
  std::copy(in.begin(), in.end(), out.data());
  return out;
}
}  // namespace test

/// ConfigEnumAllOccupations:
///
/// Notes:
/// - If enumerating on all sites in a supercell: guaranteed primitive &
/// canonical, no duplicates
/// - If enumerating on selected sites and not all sites: guaranteed primitive,
/// may be
///   non-canonical, no duplicates
class ConfigEnumAllOccupationsTest : public testing::Test {
 protected:
  std::shared_ptr<CASM::Structure const> shared_prim;

  ConfigEnumAllOccupationsTest()
      : shared_prim(std::make_shared<CASM::Structure const>(test::ZrO_prim())) {
  }

  // allows checking number of configurations enumerated and histogram of
  // composition (using number of each type of occupant in the entire
  // configuration)
  void check(ConfigEnumInput const &initial_state,
             Index expected_configurations_size,
             test::OccupationHistogram const &expected_occ_histogram) {
    ConfigEnumAllOccupations enumerator{initial_state};
    std::vector<Configuration> configurations{enumerator.begin(),
                                              enumerator.end()};

    auto molecule_name = xtal::struc_molecule_name(*shared_prim);
    auto occ_histogram =
        test::make_occ_histogram(configurations.begin(), configurations.end());

    // // Uncomment to see what was actually enumerated:
    // for(auto const &count: occ_histogram) {
    //   std::cout << "{ occupation: { ";
    //   test::print_histogram(std::cout, count.first, molecule_name);
    //   std::cout << "},  occurances: " << count.second << " } " << std::endl;
    // }

    EXPECT_EQ(configurations.size(), expected_configurations_size);
    for (auto const &expected : expected_occ_histogram) {
      std::stringstream ss;
      ss << "expected: { occupation: { ";
      test::print_histogram(ss, expected.first, molecule_name);
      ss << "},  occurances: " << expected.second << " } " << std::endl;

      auto found = occ_histogram.find(expected.first);
      if (found == occ_histogram.end()) {
        EXPECT_TRUE(found != occ_histogram.end()) << ss.str();
      } else {
        ss << " found: { occupation: { ";
        test::print_histogram(ss, found->first, molecule_name);
        ss << "},  occurances: " << found->second << " } " << std::endl;

        EXPECT_EQ(occ_histogram[expected.first], expected.second) << ss.str();
      }
    }
    EXPECT_EQ(occ_histogram.size(), expected_occ_histogram.size());
  }
};

TEST_F(ConfigEnumAllOccupationsTest, Test1) {
  using namespace test;

  // Expect 3 unique configurations for T==I: (2Va,0O), (1Va,1O), (0Va,2O)
  // Note: no equivalents, only primitives

  Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
  Supercell supercell{shared_prim, T};

  Index expected_configurations_size = 3;

  // each entry: { {n_Zr, n_Va, n_O}, #occurances }
  test::OccupationHistogram expected_occ_histogram = {
      {to_VectorXi({2, 0, 2}), 1},
      {to_VectorXi({2, 1, 1}), 1},
      {to_VectorXi({2, 2, 0}), 1}};

  ConfigEnumInput initial_state{supercell};
  check(initial_state, expected_configurations_size, expected_occ_histogram);
}

TEST_F(ConfigEnumAllOccupationsTest, Test2) {
  using namespace test;

  Eigen::Matrix3l T;
  T << 1, 0, 0, 0, 1, 0, 0, 0, 2;
  Supercell supercell{shared_prim, T};

  Index expected_configurations_size = 3;

  // each entry: { {n_Zr, n_Va, n_O}, #occurances }
  test::OccupationHistogram expected_occ_histogram = {
      {to_VectorXi({4, 1, 3}), 1},
      {to_VectorXi({4, 2, 2}), 1},
      {to_VectorXi({4, 3, 1}), 1}};

  ConfigEnumInput initial_state{supercell};
  check(initial_state, expected_configurations_size, expected_occ_histogram);
}

TEST_F(ConfigEnumAllOccupationsTest, Test3) {
  using namespace test;

  Eigen::Matrix3l T;
  T << 2, 0, 0, 0, 2, 0, 0, 0, 2;
  Supercell supercell{shared_prim, T};

  Index expected_configurations_size = 585;

  // each entry: { {n_Zr, n_Va, n_O}, #occurances }
  test::OccupationHistogram expected_occ_histogram = {
      {to_VectorXi({16, 1, 15}), 1},  {to_VectorXi({16, 2, 14}), 2},
      {to_VectorXi({16, 3, 13}), 9},  {to_VectorXi({16, 4, 12}), 20},
      {to_VectorXi({16, 5, 11}), 43}, {to_VectorXi({16, 6, 10}), 70},
      {to_VectorXi({16, 7, 9}), 95},  {to_VectorXi({16, 8, 8}), 105},
      {to_VectorXi({16, 9, 7}), 95},  {to_VectorXi({16, 10, 6}), 70},
      {to_VectorXi({16, 11, 5}), 43}, {to_VectorXi({16, 12, 4}), 20},
      {to_VectorXi({16, 13, 3}), 9},  {to_VectorXi({16, 14, 2}), 2},
      {to_VectorXi({16, 15, 1}), 1},
  };

  ConfigEnumInput initial_state{supercell};
  check(initial_state, expected_configurations_size, expected_occ_histogram);
}

TEST_F(ConfigEnumAllOccupationsTest, TestClusterSites1) {
  using namespace test;

  Eigen::Matrix3l T = Eigen::Matrix3l::Identity();
  auto shared_supercell = std::make_shared<Supercell const>(shared_prim, T);

  Configuration background_configuration{shared_supercell};

  // only enumerate occupation on one site, site index=2 is fixed to occ value 0
  // (Va)
  std::set<Index> site_indices{2};

  Index expected_configurations_size = 2;

  // each entry: { {n_Zr, n_Va, n_O}, #occurances }
  test::OccupationHistogram expected_occ_histogram = {
      {to_VectorXi({2, 1, 1}), 1}, {to_VectorXi({2, 2, 0}), 1}};

  ConfigEnumInput initial_state{background_configuration, site_indices};
  check(initial_state, expected_configurations_size, expected_occ_histogram);
}

TEST_F(ConfigEnumAllOccupationsTest, TestClusterSites2) {
  using namespace test;

  Eigen::Matrix3l T;
  T << 2, 0, 0, 0, 2, 0, 0, 0, 2;
  auto shared_supercell = std::make_shared<Supercell const>(shared_prim, T);

  Configuration background_configuration{shared_supercell};
  std::set<Index> site_indices{16, 17};

  Index expected_configurations_size = 2;

  // TODO: currently skips non-primitive, outputs non-canonical; change to
  // output all? set via params?

  // each entry: { {n_Zr, n_Va, n_O}, #occurances }
  test::OccupationHistogram expected_occ_histogram = {
      {to_VectorXi({16, 15, 1}), 2},
  };

  ConfigEnumInput initial_state{background_configuration, site_indices};
  check(initial_state, expected_configurations_size, expected_occ_histogram);
}

// // TODO: move to database/ConfigDatabase integration tests
// TEST(ConfigEnumAllOccupationsTest, TestWithDatabase) {
//
//   // ConfigEnumAllOccupations: ZrO enumeration, w/ database
//
//   auto shared_prim = std::make_shared<Structure const>(test::ZrO_prim());
//   auto title = shared_prim->structure().title();
//   auto project_settings = make_default_project_settings(*shared_prim, title);
//   PrimClex primclex {project_settings, shared_prim};
//
//   xtal::ScelEnumProps scel_enum_props {1, 5};
//   ScelEnumByProps supercell_enumerator {shared_prim, scel_enum_props};
//   for(auto const &supercell : supercell_enumerator) {
//     make_canonical_and_insert(supercell_enumerator, supercell,
//     primclex.db<Supercell>());
//   }
//
//   // order should not change as long as supercell comparison is not changed
//   std::vector<Index> number_expected {3, 3, 4, 3, 16, 10, 10, 10, 10, 18, 24,
//   24, 21, 15, 27, 27, 42, 24, 24, 21}; std::vector<Index> number_enumerated;
//   bool primitive_only = true; //
//   for(auto const &supercell : primclex.db<Supercell>()) {
//     ConfigEnumAllOccupations enumerator {supercell};
//     for(auto const &configuration : primclex.db<Configuration>()) {
//       make_canonical_and_insert(enumerator,
//                                 configuration,
//                                 primclex.db<Supercell>(),
//                                 primclex.db<Supercell>(),
//                                 primitive_only);
//     }
//     auto range = primclex.db<Configuration>().scel_range(supercell.name());
//     number_enumerated.push_back(std::distance(range.begin(), range.end()));
//   }
//   EXPECT_EQ(number_enumerated.size(), number_expected.size());
//   EXPECT_EQ(number_enumerated, number_expected);
// }
//
// // TODO: remove
// TEST(ConfigEnumAllOccupationsTest, TestOld) {
//
//   // tests ConfigEnumAllOccupations
//   ScopedNullLogging logging;
//
//   // read test file
//   fs::path test_cases_path(autotools::abs_srcdir() +
//   "/tests/unit/clex/ConfigEnumAllOccupations_test_cases.json"); jsonParser
//   tests(test_cases_path); double tol = TOL; fs::path
//   test_proj_dir(autotools::abs_srcdir() + "/tests/unit/clex/test_proj");
//
//   for(auto test_it = tests.begin(); test_it != tests.end(); ++test_it) {
//
//     // input and expected output data
//     jsonParser &j = *test_it;
//
//     // if false: print calculated results if no test data; if true: suppress
//     bool quiet = false;
//     j.get_else(quiet, "quiet", false);
//
//     EXPECT_TRUE(j.contains("title")) << "test case 'title' is required";
//     EXPECT_TRUE(j.contains("prim")) << "test case 'prim' is required";
//     EXPECT_TRUE(j.contains("min_vol")) << "test case 'min_vol' is required";
//     EXPECT_TRUE(j.contains("max_vol")) << "test case 'max_vol' is required";
//
//     // generate prim
//     Structure prim(read_prim(j["prim"], tol));
//
//     // clean up test proj
//     if(fs::exists(test_proj_dir / ".casm")) {
//       fs::remove_all(test_proj_dir);
//     }
//
//     fs::create_directory(test_proj_dir);
//
//     // build a project
//     auto project_name = j["title"].get<std::string>();
//     auto project_settings = make_default_project_settings(prim, project_name,
//     test_proj_dir); build_project(project_settings, prim);
//
//     // read primclex
//     PrimClex primclex(test_proj_dir);
//
//     // generate supercells
//     xtal::ScelEnumProps enum_props(j["min_vol"].get<int>(),
//     j["max_vol"].get<int>() + 1); ScelEnumByProps
//     scel_enum(primclex.shared_prim(), enum_props); for(const auto &scel :
//     scel_enum) {
//       (void) scel;
//     }
//
//     // run checks:
//     jsonParser json_scel;
//     json_scel = (Index) primclex.generic_db<Supercell>().size();
//     check("Nscel", j, json_scel, test_cases_path, quiet);
//
//     // generate configurations
//     jsonParser json = jsonParser::array();
//     for(auto &scel : primclex.generic_db<Supercell>()) {
//       ConfigEnumAllOccupations e(scel);
//       json.push_back(std::distance(e.begin(), e.end()));
//     }
//
//     // run checks:
//     check("Nconfigs", j, json, test_cases_path, quiet);
//
//     // ... add more here ...
//
//     // clean up test proj
//     if(fs::exists(test_proj_dir / ".casm")) {
//       fs::remove_all(test_proj_dir);
//     }
//   }
// }
//
// // TODO: move to app_enum_methods_ConfigEnumAllOccupationsInterface_test
// TEST(ConfigEnumAllOccupationsTest, InterfaceTest1) {
//
//   // create a project
//   test::FCCTernaryProj proj;
//   proj.check_init();
//
//   // construct PrimClex
//   ScopedNullLogging logging;
//   PrimClex primclex(proj.dir);
//
//   // --dry-run test
//   {
//     std::string cli_str = "casm enum --method ScelEnum --max 4 --dry-run";
//     run_enum_interface<ScelEnumInterface>(cli_str, primclex);
//   }
//   EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
//   primclex.generic_db<Supercell>().close();
//   primclex.generic_db<Supercell>().open();
//   EXPECT_EQ(primclex.generic_db<Supercell>().size(), 0);
//
//
//   {
//     std::string cli_str = "casm enum --method ScelEnum --max 4";
//     run_enum_interface<ScelEnumInterface>(cli_str, primclex);
//   }
//   EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
//   primclex.generic_db<Supercell>().close();
//   primclex.generic_db<Supercell>().open();
//   EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
//
//   // --dry-run test
//   {
//     std::string cli_str = "casm enum --method ConfigEnumAllOccupations -a
//     --dry-run";
//     run_enum_interface<ConfigEnumAllOccupationsInterface>(cli_str, primclex);
//   }
//   EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
//   EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);
//   primclex.generic_db<Configuration>().close();
//   primclex.generic_db<Configuration>().open();
//   EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
//   EXPECT_EQ(primclex.generic_db<Configuration>().size(), 0);
//
//   {
//     std::string cli_str = "casm enum --method ConfigEnumAllOccupations -a";
//     run_enum_interface<ConfigEnumAllOccupationsInterface>(cli_str, primclex);
//   }
//   EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
//   EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);
//   primclex.generic_db<Configuration>().close();
//   primclex.generic_db<Configuration>().open();
//   EXPECT_EQ(primclex.generic_db<Supercell>().size(), 13);
//   EXPECT_EQ(primclex.generic_db<Configuration>().size(), 126);
//
// }
