#include "casm/clusterography/io/json/ClusterSpecs_json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "crystallography/TestStructures.hh"
#include "gtest/gtest.h"

using namespace CASM;

TEST(ClusterSpecsJSONTest, ParseTest1) {
  jsonParser json = jsonParser::parse(std::string(R"({
    "method": "periodic_max_length",
    "params": {
      "orbit_branch_specs": {
        "2": {"max_length": 1.01}
      }
    }
  })"));

  auto shared_prim =
      std::make_shared<Structure const>(test::SimpleCubic_GLstrain_prim());

  InputParser<ClusterSpecs> parser{json, shared_prim};

  // Error to throw if parsing fails
  std::runtime_error error_if_invalid{"Failed to parse cluster specs JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  // In this example we know to expect prim periodic orbits.
  EXPECT_EQ(parser.value->periodicity_type(),
            CASM::CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC);
  auto const &cspecs =
      static_cast<PeriodicMaxLengthClusterSpecs const &>(*parser.value);
  EXPECT_EQ(cspecs.max_length.size(), 3);
}

TEST(ClusterSpecsJSONTest, ParseTest2) {
  jsonParser json = jsonParser::parse(std::string(R"({
    "method": "periodic_max_length",
    "params": {
      "orbit_branch_specs": {}
    }
  })"));

  auto shared_prim =
      std::make_shared<Structure const>(test::SimpleCubic_GLstrain_prim());

  InputParser<ClusterSpecs> parser{json, shared_prim};

  // Error to throw if parsing fails
  std::runtime_error error_if_invalid{"Failed to parse cluster specs JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  // In this example we know to expect prim periodic orbits.
  EXPECT_EQ(parser.value->periodicity_type(),
            CASM::CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC);
  auto const &cspecs =
      static_cast<PeriodicMaxLengthClusterSpecs const &>(*parser.value);
  EXPECT_EQ(cspecs.max_length.size(), 0);
}
