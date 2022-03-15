#include "casm/basis_set/OccupationDoFTraits.hh"

#include "Common.hh"
#include "FCCTernaryProj.hh"
#include "ZrOProj.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/crystallography/Structure.hh"
#include "gtest/gtest.h"

using namespace CASM;

DoF_impl::OccupationDoFSpecs parse_valid_occ_specs(std::string str,
                                                   Structure const &prim) {
  jsonParser json = jsonParser::parse(str);
  InputParser<DoF_impl::OccupationDoFSpecs> parser{json, prim};
  std::runtime_error error_if_invalid{"Failed to parse OccupationDoFSpecs"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
  return *parser.value;
}

InputParser<DoF_impl::OccupationDoFSpecs> parse_invalid_occ_specs(
    std::string str, Structure const &prim) {
  jsonParser json = jsonParser::parse(str);
  return InputParser<DoF_impl::OccupationDoFSpecs>{json, prim};
}

TEST(OccupationDoFTraitsTest, ParseJSONoccupation) {
  std::string str = R"({
    "site_basis_functions" : "occupation"
  })";

  Structure prim{test::FCC_ternary_prim()};
  DoF_impl::OccupationDoFSpecs occ_specs = parse_valid_occ_specs(str, prim);

  ASSERT_EQ(occ_specs.site_basis_function_type,
            DoF_impl::SITE_BASIS_FUNCTION_TYPE::OCCUPATION);
}

TEST(OccupationDoFTraitsTest, ParseJSONchebychev) {
  std::string str = R"({
    "site_basis_functions" : "chebychev"
  })";

  Structure prim{test::FCC_ternary_prim()};
  DoF_impl::OccupationDoFSpecs occ_specs = parse_valid_occ_specs(str, prim);

  ASSERT_EQ(occ_specs.site_basis_function_type,
            DoF_impl::SITE_BASIS_FUNCTION_TYPE::CHEBYCHEV);
}

TEST(OccupationDoFTraitsTest, ParseJSONcomposition_test1) {
  std::string str = R"({
    "site_basis_functions" : [
      {
        "sublat_indices": [0],
        "composition": {"A": 0.25, "B": 0.25, "C": 0.5}
      }
    ]
  })";

  Structure prim{test::FCC_ternary_prim()};
  DoF_impl::OccupationDoFSpecs occ_specs = parse_valid_occ_specs(str, prim);

  ASSERT_EQ(occ_specs.site_basis_function_type,
            DoF_impl::SITE_BASIS_FUNCTION_TYPE::COMPOSITION);
}

TEST(OccupationDoFTraitsTest, ParseJSONcomposition_test2) {
  // composition values do not need to sum to 1.0 -- they will be normalized
  std::string str = R"({
    "site_basis_functions" : [
      {
        "sublat_indices": [0],
        "composition": {"A": 25, "B": 25, "C": 50}
      }
    ]
  })";

  Structure prim{test::FCC_ternary_prim()};
  DoF_impl::OccupationDoFSpecs occ_specs = parse_valid_occ_specs(str, prim);

  ASSERT_EQ(occ_specs.site_basis_function_type,
            DoF_impl::SITE_BASIS_FUNCTION_TYPE::COMPOSITION);

  std::vector<double> sublat_prob_vec =
      sublat_values_vec(occ_specs.sublat_values, 0, "composition",
                        prim.basis()[0].allowed_occupants(), true);
  ASSERT_EQ(sublat_prob_vec, (std::vector<double>{0.25, 0.25, 0.5}));
}

TEST(OccupationDoFTraitsTest, ParseJSONcomposition_test3) {
  std::string str = R"({
    "site_basis_functions" : [
      {
        "sublat_indices": [2, 3],
        "composition": {"O": 0.2, "Va": 0.8}
      }
    ]
  })";

  Structure prim{test::ZrO_prim()};
  DoF_impl::OccupationDoFSpecs occ_specs = parse_valid_occ_specs(str, prim);

  ASSERT_EQ(occ_specs.site_basis_function_type,
            DoF_impl::SITE_BASIS_FUNCTION_TYPE::COMPOSITION);
  ASSERT_EQ(occ_specs.sublat_values.size(), 1);

  std::vector<double> sublat_prob_vec_2 =
      sublat_values_vec(occ_specs.sublat_values, 2, "composition",
                        prim.basis()[2].allowed_occupants(), true);
  ASSERT_EQ(sublat_prob_vec_2, (std::vector<double>{0.8, 0.2}));

  std::vector<double> sublat_prob_vec_3 =
      sublat_values_vec(occ_specs.sublat_values, 3, "composition",
                        prim.basis()[3].allowed_occupants(), true);
  ASSERT_EQ(sublat_prob_vec_3, (std::vector<double>{0.8, 0.2}));
}
