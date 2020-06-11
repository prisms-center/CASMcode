#include "casm/casm_io/Log.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/casm_io/json/InputParser.hh"
#include "casm/clusterography/io/json/ClusterOrbits_json_io.hh"
#include "casm/clusterography/io/json/IntegralCluster_json_io.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/crystallography/Structure.hh"
// #include "casm/global/enum/json_io.hh"

#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"

namespace CASM {

  /// Write custom orbit specs to JSON
  jsonParser &to_json(const IntegralClusterOrbitGenerator &orbit_generator, jsonParser &json) {
    to_json(orbit_generator.prototype, json);
    json["include_subclusters"] = orbit_generator.include_subclusters;
    return json;
  }

  /// Parse custom orbit specs from JSON
  void parse(
    InputParser<std::vector<IntegralClusterOrbitGenerator>> &parser,
    const std::shared_ptr<Structure const> &shared_prim) {

    // "orbit_specs": [
    //   {
    //     "coordinate_mode" : "Direct",
    //     "prototype" : [
    //       [ 0.000000000000, 0.000000000000, 0.000000000000 ],
    //       [ 1.000000000000, 0.000000000000, 0.000000000000 ],
    //       [ 2.000000000000, 0.000000000000, 0.000000000000 ],
    //       [ 3.000000000000, 0.000000000000, 0.000000000000 ]],
    //     "include_subclusters" : true
    //   },
    //   ...
    // ]

    const jsonParser &json = parser.self;

    if(!json.is_array()) {
      parser.error.insert("Error reading orbit generating clusters: Expected a JSON array");
      return;
    }

    parser.value = notstd::make_unique<std::vector<IntegralClusterOrbitGenerator>>();
    auto &custom_generators = *parser.value;
    try {
      // for each custom orbit
      Index i = 0;
      for(auto it = json.begin(); it != json.end(); ++it) {

        // read orbit generating cluster from JSON
        auto relpath = parser.path / std::to_string(i);
        auto subparser = std::make_shared<InputParser<IntegralCluster>>(parser.input, relpath, true, shared_prim);

        if(subparser->valid()) {
          // check if subclusters should be included (yes by default)
          bool include_subclusters;
          parser.optional_else(include_subclusters, "include_subclusters", true);

          custom_generators.emplace_back(*(subparser->value), include_subclusters);
        }
        else {
          // if can't read cluster, insert errors and stop parsing
          parser.insert(relpath, subparser);
          return;
        }
        ++i;
      }
    }
    catch(std::exception &e) {
      parser.error.insert(std::string("Error: Could not read orbit generating clusters: ") + e.what());
    }
  }
}
