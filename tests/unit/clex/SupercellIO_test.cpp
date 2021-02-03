#include "ProjectBaseTest.hh"
#include "casm/casm_io/Log.hh"
#include "casm/casm_io/dataformatter/DataFormatterTools_impl.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/SupercellIO_impl.hh"
#include "casm/database/DatabaseTypes_impl.hh"
#include "casm/database/ScelDatabaseTools_impl.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

namespace {
Eigen::Matrix3l _fcc_conventional_transf_mat() {
  Eigen::Matrix3l transf_mat;
  transf_mat << -1, 1, 1, 1, -1, 1, 1, 1, -1;
  return transf_mat;
}
}  // namespace

class SupercellFormatterTest : public test::ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  SupercellFormatterTest()
      : test::ProjectBaseTest(test::FCC_ternary_prim(),
                              "SupercellFormatterTest",
                              jsonParser::parse(clex_basis_specs_str())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())),
        dict(make_attribute_dictionary<Supercell>()),
        supercell_db(primclex_ptr->template db<Supercell>()) {
    this->enumerate_supercells();

    shared_supercell->set_primclex(primclex_ptr.get());
  }

  void enumerate_supercells();

  // conventional FCC unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;
  DataFormatterDictionary<Supercell> dict;
  DB::Database<Supercell> &supercell_db;
};

std::string SupercellFormatterTest::clex_basis_specs_str() {
  return R"({
"basis_function_specs" : {
  "dof_specs": {
    "occ": {
      "site_basis_functions" : "occupation"
    }
  }
},
"cluster_specs": {
  "method": "periodic_max_length",
  "params": {
    "orbit_branch_specs" : {
      "2" : {"max_length" : 4.01},
      "3" : {"max_length" : 3.01}
    }
  }
}
})";
}

void SupercellFormatterTest::enumerate_supercells() {
  xtal::ScelEnumProps enumeration_params{1, 5, "abc",
                                         Eigen::Matrix3i::Identity()};

  ScelEnumByProps enumerator{shared_prim, enumeration_params};

  for (Supercell const &supercell : enumerator) {
    supercell.set_primclex(primclex_ptr.get());
    auto result =
        make_canonical_and_insert(enumerator, supercell, supercell_db);
  }
}

TEST_F(SupercellFormatterTest, IsSupercellOfTest1) {
  // conventional FCC is a supercell of the prim unit cell

  ScelIO::IsSupercellOf f{};
  EXPECT_TRUE(f.parse_args("SCEL1_1_1_1_0_0_0"));
  EXPECT_TRUE(f.init(*shared_supercell));
  EXPECT_TRUE(f.validate(*shared_supercell));
  EXPECT_TRUE(f.evaluate(*shared_supercell));

  jsonParser json;
  f.to_json(*shared_supercell, json);
  EXPECT_TRUE(json.is_bool());
  EXPECT_EQ(json.get<bool>(), true);
}

TEST_F(SupercellFormatterTest, IsSupercellOfTest2) {
  // conventional FCC is not a supercell of the prim unit cell

  ScelIO::IsSupercellOf f{};
  EXPECT_TRUE(f.parse_args(shared_supercell->name()));
  EXPECT_TRUE(f.init(*shared_supercell));
  EXPECT_TRUE(f.validate(*shared_supercell));
  EXPECT_TRUE(f.evaluate(*shared_supercell));

  jsonParser json;
  auto scel_it = supercell_db.find("SCEL1_1_1_1_0_0_0");
  f.to_json(*scel_it, json);
  EXPECT_TRUE(json.is_bool());
  EXPECT_EQ(json.get<bool>(), false);
}

TEST_F(SupercellFormatterTest, IsSupercellOfFromDictTest) {
  DataFormatter<Supercell> formatter =
      dict.parse("is_supercell_of(SCEL1_1_1_1_0_0_0)");

  jsonParser json;
  formatter.to_json(*shared_supercell, json);
  EXPECT_TRUE(json.contains("is_supercell_of"));
  EXPECT_TRUE(json["is_supercell_of"].is_bool());
  EXPECT_EQ(json["is_supercell_of"].get<bool>(), true);
}

TEST_F(SupercellFormatterTest, IsSupercellOfInvalidArgumentTest) {
  /// The error will occur in "SupercellCheckBase::init, when initialized with
  /// the first supercell to be formatted

  DataFormatter<Supercell> formatter =
      dict.parse("is_supercell_of(not_a_supercell_name)");

  jsonParser json;
  EXPECT_THROW(formatter.to_json(*shared_supercell, json), std::runtime_error);
}

// --- lattice_params ---

TEST_F(SupercellFormatterTest, LatticeParamsTest) {
  ScelIO::GenericVectorXdScelFormatter f = ScelIO::lattice_params();
  EXPECT_TRUE(f.parse_args(""));
  EXPECT_TRUE(f.init(*shared_supercell));
  EXPECT_TRUE(f.validate(*shared_supercell));
  Eigen::VectorXd result = f.evaluate(*shared_supercell);

  jsonParser json;
  f.to_json(*shared_supercell, json);
  EXPECT_TRUE(json.is_array());
}

TEST_F(SupercellFormatterTest, LatticeParamsFromDictTest) {
  DataFormatter<Supercell> formatter = dict.parse("lattice_params");

  jsonParser json;
  formatter.to_json(*shared_supercell, json);
  EXPECT_TRUE(json.contains("lattice_params"));
  EXPECT_TRUE(json["lattice_params"].is_array());
  EXPECT_EQ(json["lattice_params"].size(), 6);
}

// quick and dirty test each standard Supercell formatter, individually
// - any that require args or are otherwise invalid should fail gracefully
// - use to quickly identify which queries are problematic
TEST_F(SupercellFormatterTest, QueryEachTest) {
  std::set<std::string> requires_arg{"is_unitcell_of", "is_supercell_of"};

  for (auto it = dict.begin(); it != dict.end(); ++it) {
    std::string query_name = it.base()->first;
    // std::cout << "query: " << query_name << std::endl;
    try {
      DataFormatter<Supercell> formatter = dict.parse(query_name);
      jsonParser json;
      formatter.to_json(*shared_supercell, json);
      // std::cout << json << std::endl;
    } catch (std::exception const &e) {
      EXPECT_TRUE(requires_arg.find(query_name) != requires_arg.end());
      // std::cout << "Caught error: " << e.what() << std::endl;
    }
  }
}

// quick and dirty test all standard Supercell formatters
// - any that require args or are otherwise invalid should fail gracefully by
//   throwing exceptions
TEST_F(SupercellFormatterTest, QueryAllTest) {
  std::set<std::string> requires_arg{"is_unitcell_of", "is_supercell_of"};

  for (auto it = dict.begin(); it != dict.end(); ++it) {
    std::string query_name = it.base()->first;
    // std::cout << "query: " << query_name << std::endl;
    try {
      DataFormatter<Supercell> formatter = dict.parse(query_name);
      jsonParser json;
      formatter(supercell_db.begin(), supercell_db.end()).to_json(json);
      // std::cout << json << std::endl;
    } catch (std::exception const &e) {
      EXPECT_TRUE(requires_arg.find(query_name) != requires_arg.end());
      // std::cout << "Caught error: " << e.what() << std::endl;
    }
  }
}
