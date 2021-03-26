#include <vector>

#include "casm/casm_io/Log.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ECIContainer.hh"

namespace CASM {

///  Make ECIContainer from JSON (eci.json file)
///
/// Format:
/// \code
/// {
///   "orbits": [
///     {
///       "cluster_functions": [
///         {
///           "\\Phi_{<linear_function_index>}": <formula>,
///           "eci": <number>,
///           "linear_function_index": <integer>
///         },
///         ...
///       ],
///       ...
///     },
///     ...
///   ],
///   ...
/// }
/// \endcode
///
/// Note:
/// - Prior to version 1.0.X, the "orbits" array was named "cluster_functions"
///   and the "linear_orbit_index" was named "linear_cluster_index". For
///   compatibility, that format is also accepted.
///
void parse(InputParser<ECIContainer> &parser) {
  if (parser.self.find("orbits") == parser.self.end()) {
    std::stringstream msg;
    msg << "Error parsing ECI: 'orbits' array not found";
    parser.error.insert(msg.str());
    return;
  }

  std::vector<double> value;
  std::vector<ECIContainer::size_type> index;

  jsonParser const &orbits_json = parser.self["orbits"];
  fs::path orbits_path{"orbits"};
  Index orbit_index = 0;
  for (auto const &orbit_json : orbits_json) {
    fs::path orbit_path = orbits_path / std::to_string(orbit_index);
    if (!orbit_json.contains("cluster_functions")) {
      parser.insert_error(orbit_path, "Error: missing 'cluster_functions'");
      return;
    }
    jsonParser const &cluster_functions_json = orbit_json["cluster_functions"];
    Index function_index = 0;
    for (auto const &function_json : cluster_functions_json) {
      if (function_json.find("eci") != function_json.end()) {
        fs::path func_path =
            orbit_path / "cluster_functions" / std::to_string(function_index);

        double eci_value;
        parser.require(eci_value, func_path / "eci");
        value.push_back(eci_value);

        ECIContainer::size_type linear_function_index;
        parser.require(linear_function_index,
                       func_path / "linear_function_index");
        index.push_back(linear_function_index);
      }
      ++function_index;
    }
    ++orbit_index;
  }
  parser.value = notstd::make_unique<ECIContainer>(value.begin(), value.end(),
                                                   index.begin());
}
}  // namespace CASM
