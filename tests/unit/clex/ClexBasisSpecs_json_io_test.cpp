#include "casm/clex/io/json/ClexBasisSpecs_json_io.hh"

#include "casm/basis_set/DoFTraits.hh"
#include "casm/basis_set/OccupationDoFTraits.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ClexBasisSpecs.hh"
#include "casm/crystallography/Structure.hh"
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

  ClexBasisSpecs const &basis_set_specs = *parser.value;

  // Check BasisFunctionSpecs contents
  BasisFunctionSpecs const &basis_function_specs =
      basis_set_specs.basis_function_specs;
  EXPECT_EQ(basis_function_specs.dof_keys.size(), 1);
  EXPECT_EQ(basis_function_specs.dof_keys[0], "occ");

  DoF_impl::OccupationDoFSpecs const &occ_specs =
      get<DoF_impl::OccupationDoFSpecs>("occ", basis_function_specs);
  EXPECT_EQ(occ_specs.site_basis_function_type,
            DoF_impl::SITE_BASIS_FUNCTION_TYPE::OCCUPATION);
  EXPECT_EQ(occ_specs.sublat_values.size(), 0);

  ClusterSpecs const &cluster_specs = *basis_set_specs.cluster_specs;
  EXPECT_EQ(cluster_specs.name(), "periodic_max_length");
  EXPECT_EQ(cluster_specs.periodicity_type(),
            CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC);
}

TEST(ClexBasisSpecsJSONTest, ParseTest2) {
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
        "orbit_branch_specs": {}
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

  ClexBasisSpecs const &basis_set_specs = *parser.value;

  // Check BasisFunctionSpecs contents
  BasisFunctionSpecs const &basis_function_specs =
      basis_set_specs.basis_function_specs;
  EXPECT_EQ(basis_function_specs.dof_keys.size(), 1);
  EXPECT_EQ(basis_function_specs.dof_keys[0], "occ");

  DoF_impl::OccupationDoFSpecs const &occ_specs =
      get<DoF_impl::OccupationDoFSpecs>("occ", basis_function_specs);
  EXPECT_EQ(occ_specs.site_basis_function_type,
            DoF_impl::SITE_BASIS_FUNCTION_TYPE::OCCUPATION);
  EXPECT_EQ(occ_specs.sublat_values.size(), 0);

  ClusterSpecs const &cluster_specs = *basis_set_specs.cluster_specs;
  EXPECT_EQ(cluster_specs.name(), "periodic_max_length");
  EXPECT_EQ(cluster_specs.periodicity_type(),
            CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC);
}

TEST(ClexBasisSpecsJSONTest, ParseTest3) {
  jsonParser json = jsonParser::parse(std::string(R"({
    "basis_function_specs": {
      "global_max_poly_order": 3
    },
    "cluster_specs": {
      "method": "periodic_max_length",
      "params": {
        "orbit_branch_specs": {}
      }
    }
  })"));

  auto shared_prim =
      std::make_shared<Structure const>(test::SimpleCubic_GLstrain_prim());
  ParsingDictionary<DoFType::Traits> const *dof_dict = &DoFType::traits_dict();
  InputParser<ClexBasisSpecs> parser{json, shared_prim, dof_dict};

  // Error to throw if parsing fails
  std::runtime_error error_if_invalid{"Failed to parse clex basis specs JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  ClexBasisSpecs const &basis_set_specs = *parser.value;

  // Check BasisFunctionSpecs contents
  BasisFunctionSpecs const &basis_function_specs =
      basis_set_specs.basis_function_specs;
  EXPECT_EQ(basis_function_specs.dof_keys.size(), 1);
  EXPECT_EQ(basis_function_specs.dof_keys[0], "GLstrain");

  EXPECT_EQ(basis_function_specs.dof_specs.size(), 0);

  // check cluster_specs (generic)
  {
    ClusterSpecs const &cluster_specs = *basis_set_specs.cluster_specs;
    EXPECT_EQ(cluster_specs.name(), "periodic_max_length");
    EXPECT_EQ(cluster_specs.periodicity_type(),
              CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC);
  }

  // check cluster_specs (static_cast<PeriodicMaxLengthClusterSpecs const &>)
  {
    PeriodicMaxLengthClusterSpecs const &cluster_specs =
        static_cast<PeriodicMaxLengthClusterSpecs const &>(
            *basis_set_specs.cluster_specs);
    EXPECT_EQ(cluster_specs.name(), "periodic_max_length");
    EXPECT_EQ(cluster_specs.periodicity_type(),
              CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC);

    // Null cluster orbit is always included
    EXPECT_EQ(cluster_specs.max_length.size(), 1);
  }
}
