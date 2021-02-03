#include "ProjectBaseTest.hh"
#include "casm/casm_io/Log.hh"
#include "casm/clex/PrimClex_impl.hh"
#include "casm/clex/ScelEnum.hh"
#include "casm/clex/Supercell.hh"
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

class SupercellFormatterTest : public ProjectBaseTest {
 protected:
  static std::string clex_basis_specs_str();

  SupercellFormatterTest()
      : test::ProjectBaseTest(test::FCC_ternary_prim(),
                              "SupercellFormatterTest",
                              jsonParser::parse(clex_basis_specs_str())),
        shared_supercell(std::make_shared<CASM::Supercell>(
            shared_prim, _fcc_conventional_transf_mat())),
        dict(make_attribute_dictionary<Supercell>()) {
    this->enumerate_supercells();

    shared_supercell->set_primclex(primclex_ptr.get());
  }

  void enumerate_supercells();

  // conventional FCC unit cell
  std::shared_ptr<CASM::Supercell> shared_supercell;

  DataFormatterDictionary<Supercell> dict;
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
})";
}

void SupercellFormatterTest::enumerate_supercells() {
  CASM::xtal::ScelEnumProps enumeration_params{1, 5, "abc",
                                               Eigen::Matrix3i::Identity()};

  CASM::ScelEnumByProps enumerator{shared_prim, enumeration_params};

  for (Supercell const &supercell : enumerator) {
    auto result = make_canonical_and_insert(enumerator, supercell,
                                            primclex_ptr->db<Supercell>());
  }
}

TEST_F(SupercellFormatterTest, IsSupercellOfTest) {
  DataFormatter<Supercell> formatter =
      dict.parse("is_supercell_of(SCEL1_1_1_1_0_0_0)");

  jsonParser json;
  formatter.to_json(*shared_supercell, json);
  std::cout << json << std::endl;
}
