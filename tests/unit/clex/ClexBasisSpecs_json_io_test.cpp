#include "casm/clex/io/json/ClexBasisSpecs_json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

TEST(ClexBasisSpecsJSONTest, ParseTest1) {
  jsonParser json = jsonParser::parse(std::string(R"({
    "basis_function_specs": {
      "dof_specs": {
        "occ": {
          "site_basis_functions" : "occupation"
        }
      }
    },
    "cluster_specs": {
      "method": "periodic_max_length",
      "params": {
        "orbit_branch_specs": {
          "2" : {"max_length" : 4.01},
          "3" : {"max_length" : 3.01}
        },
        "orbit_specs" : [
          {
            "coordinate_mode" : "Direct",
            "prototype" : [
              [ 0.000000000000, 0.000000000000, 0.000000000000 ],
              [ 1.000000000000, 0.000000000000, 0.000000000000 ],
              [ 2.000000000000, 0.000000000000, 0.000000000000 ],
              [ 3.000000000000, 0.000000000000, 0.000000000000 ]
            ],
            "include_subclusters" : true
          },
          {
            "coordinate_mode" : "Direct",
            "prototype" : [
              [ 0.000000000000, 0.000000000000, 0.000000000000 ],
              [ 0.000000000000, 1.000000000000, 0.000000000000 ],
              [ 0.000000000000, 0.000000000000, 1.000000000000 ],
              [ 1.000000000000, 1.000000000000, 1.000000000000 ]
            ],
            "include_subclusters" : true
          }
        ]
      }
    }
  })"));

  auto shared_prim =
      std::make_shared<Structure const>(test::FCC_ternary_prim());
  ParsingDictionary<DoFType::Traits> const *dof_dict = &DoFType::traits_dict();
  InputParser<ClexBasisSpecs> parser{json, shared_prim, dof_dict};

  // Error to throw if parsing fails
  std::runtime_error error_if_invalid{"Failed to parse clex basis specs JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);
}
