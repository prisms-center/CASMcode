#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"
#include "casm/clusterography/json_io.hh"
#include "crystallography/TestStructures.hh"  // for test::ZrO_prim
#include "gtest/gtest.h"

namespace example {

// Valid input example
CASM::jsonParser cluster_specs_json_test_valid();

// Invalid input example
CASM::jsonParser cluster_specs_json_test_invalid();

void check_orbits(
    std::shared_ptr<CASM::Structure const> const &shared_prim,
    std::vector<CASM::PrimPeriodicIntegralClusterOrbit> const &orbits);
}  // namespace example

// The InputParser template class helps parse JSON and collect errors and
// warning messages.
//
// The InputParser constructor calls an implementation function, `parse`,
// specialized by type, T, which indicates any extra input required to parse the
// JSON:
//    void parse(InputParser<T> &parser, ... any other required input ...);
//
// The InputParser itself is then constructed like:
//    jsonParser json_input = ...;
//    InputParser<T> parser {json_input, ... any other required input ...};
//
// When parsing is successful, the constructed value of type T will owned by the
// member `value`:
//    std::unique_ptr<T> InputParser<T>::value;
//
// When parsing complex objects, subparsers may be called and they will store
// their errors and warnings in a map of (path to the JSON subobject being
// parsed) : (string with an error or warning message).
//
// To check if parsing was successful and conditionally print errors & warnings,
// with the path to the location in the json where the issue occurred, do:
//     Log& log = ...;
//     if(!parser.valid()) {
//        parser.print_errors(log);
//        log << std::endl << parser.report() << std::endl << std::endl;
//        ... handle error or throw ...
//     }
//     if(parser.all_warnings().size()) {
//         parser.print_warnings(log);
//         log << std::endl << parser.report() << std::endl << std::endl;
//     }
//
// To do the above and throw an exception if the parser has any errors, use:
//     MyExceptionType error {"... message ..."};
//     report_and_throw_if_invalid(parser, log, error);
//
// Get the set of errors and warnings, including all subparsers, use:
//     std::set<std::string> warnings = parser.all_warnings();
//     std::set<std::string> errors = parser.all_errors();
//
// If parsing could not occur, then the `parse` function may leave
// `parser.value` as an empty unique_ptr<T>, otherwise `parser.value` should
// hold the resulting value. So if parsing was successful you can get the
// constructed value with `*parser.value`.

TEST(ExampleClusterographyClusterSpecsIO, ParseIntegralCluster) {
  // Construct a ZrO prim
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // try parsing IntegralCluster from JSON
  CASM::jsonParser cluster_json = CASM::jsonParser::parse(std::string(R"({
    "coordinate_mode" : "Integral",
    "prototype" : [
      [ 2, 0, 0, 0 ],
      [ 2, 1, 0, 0 ],
      [ 2, 1, 1, 0 ]
    ],
    "include_subclusters" : true
  })"));
  CASM::InputParser<CASM::IntegralCluster> cluster_parser{cluster_json,
                                                          *shared_prim};
  std::runtime_error cluster_error_if_invalid{
      "Failed to parse IntegralCluster JSON"};
  report_and_throw_if_invalid(cluster_parser, CASM::log(),
                              cluster_error_if_invalid);
}

TEST(ExampleClusterographyClusterSpecsIO,
     ParseIntegralClusterOrbitGeneratorVec) {
  // Construct a ZrO prim
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  CASM::jsonParser clustergen_json = CASM::jsonParser::parse(std::string(R"([{
    "coordinate_mode" : "Integral",
    "prototype" : [
      [ 2, 0, 0, 0 ],
      [ 2, 1, 0, 0 ],
      [ 2, 1, 1, 0 ]
    ],
    "include_subclusters" : true
  },
  {
    "coordinate_mode" : "Integral",
    "prototype" : [
      [ 2, 0, 0, 0 ],
      [ 3, 0, 0, 0 ],
      [ 2, 1, 1, 0 ]
    ],
    "include_subclusters" : true
  }
  ])"));
  CASM::InputParser<std::vector<CASM::IntegralClusterOrbitGenerator>>
      clustergen_parser{clustergen_json, *shared_prim};
  std::runtime_error clustergen_error_if_invalid{
      "Failed to parse IntegralClusterOrbitGenerator JSON"};
  report_and_throw_if_invalid(clustergen_parser, CASM::log(),
                              clustergen_error_if_invalid);
}

TEST(ExampleClusterographyClusterSpecsIO, PeriodicMaxLengthClusterSpecs) {
  // An example controlling orbit generation from JSON input
  //
  // This demonstrates equivalent orbit generation as that coded in example
  //     `302_clusterography_cluster_specs_test.cpp`
  //
  // Documentation of the JSON format is given in the file:
  //     `src/casm/clusterography/io/json/ClusterSpecs_json_io.cc`
  //
  // It is associated with the function:
  //     `void parse( InputParser<ClusterSpecs> &parser,
  //                  const std::shared_ptr<const Structure> &shared_prim,
  //                  const SymGroup &super_group);`
  //
  // An overload function chooses the shared_prim->factor_group automatically:
  //     `void parse( InputParser<ClusterSpecs> &parser,
  //                  const std::shared_ptr<const Structure> &shared_prim);`
  //
  // TODO: make the JSON format documentation available online

  // Construct a ZrO prim
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // Get JSON input
  CASM::jsonParser cluster_specs_json =
      example::cluster_specs_json_test_valid();

  // Parsing ClusterSpecs requires a shared prim
  // - This uses the `parse` overload that chooses shared_prim->factor_group()
  // automatically
  CASM::InputParser<CASM::ClusterSpecs> parser{cluster_specs_json, shared_prim};

  // Error to throw if parsing fails
  std::runtime_error error_if_invalid{"Failed to parse cluster specs JSON"};
  report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid);

  // In this example we know to expect prim periodic orbits.
  EXPECT_EQ(parser.value->periodicity_type(),
            CASM::CLUSTER_PERIODICITY_TYPE::PRIM_PERIODIC);
  CASM::ClusterSpecs::PeriodicOrbitVec orbits =
      parser.value->make_periodic_orbits(CASM::null_log());

  // Using check_orbits from example `302_clusterography_cluster_specs_test.cpp`
  example::check_orbits(shared_prim, orbits);
}

TEST(ExampleClusterographyClusterSpecsIO,
     InvalidPeriodicMaxLengthClusterSpecs) {
  // Example when JSON input is invalid

  // Construct a ZrO prim
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // Get invalid JSON input (in this case 'Integral' is misspelled as
  // 'IFntegral')
  CASM::jsonParser cluster_specs_json =
      example::cluster_specs_json_test_invalid();

  // Parsing ClusterSpecs requires a shared prim
  // - This uses the `parse` overload that chooses shared_prim->factor_group()
  // automatically
  CASM::InputParser<CASM::ClusterSpecs> parser{cluster_specs_json, shared_prim};

  EXPECT_EQ(parser.valid(), false);

  // Error to throw if parsing fails
  std::runtime_error error_if_invalid{"Failed to parse cluster specs JSON"};
  ASSERT_THROW(
      (report_and_throw_if_invalid(parser, CASM::log(), error_if_invalid)),
      std::runtime_error);

  // View casm_unit_example.log to see the error report
}

namespace example {

// Valid input, matching example `302_clusterography_cluster_specs_test.cpp`
CASM::jsonParser cluster_specs_json_test_valid() {
  return CASM::jsonParser::parse(std::string(R"({
      "method": "periodic_max_length",
      "params": {
        "orbit_branch_specs" : {
          "2" : {"max_length" : 5.17}
        },
        "orbit_specs" : [
          {
            "coordinate_mode" : "Integral",
            "prototype" : [
              [ 2, 0, 0, 0 ],
              [ 2, 1, 0, 0 ],
              [ 2, 1, 1, 0 ]
            ],
            "include_subclusters" : true
          },
          {
            "coordinate_mode" : "Integral",
            "prototype" : [
              [ 2, 0, 0, 0 ],
              [ 3, 0, 0, 0 ],
              [ 2, 1, 1, 0 ]
            ],
            "include_subclusters" : true
          },
          {
            "coordinate_mode" : "Integral",
            "prototype" : [
              [ 2, 0, 0, 0 ],
              [ 2, 0, 0, 1 ],
              [ 2, 0, 0, 2 ]
            ],
            "include_subclusters" : true
          }
        ]
      }
    })"));
}

// Invalid input example
CASM::jsonParser cluster_specs_json_test_invalid() {
  return CASM::jsonParser::parse(std::string(R"({
      "method": "periodic_max_length",
      "params": {
        "orbit_branch_specs" : {
          "2" : {"max_length" : 5.17}
        },
        "orbit_specs" : [
          {
            "coordinate_mode" : "IFntegral",
            "prototype" : [
              [ 2, 0, 0, 0 ],
              [ 2, 1, 0, 0 ],
              [ 2, 1, 1, 0 ]
            ],
            "include_subclusters" : true
          },
          {
            "coordinate_mode" : "Integral",
            "prototype" : [
              [ 2, 0, 0, 0 ],
              [ 3, 0, 0, 0 ],
              [ 2, 1, 1, 0 ]
            ],
            "include_subclusters" : true
          },
          {
            "coordinate_mode" : "Integral",
            "prototype" : [
              [ 2, 0, 0, 0 ],
              [ 2, 0, 0, 1 ],
              [ 2, 0, 0, 2 ]
            ],
            "include_subclusters" : true
          }
        ]
      }
    })"));
}

}  // namespace example
