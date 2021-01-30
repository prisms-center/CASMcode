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
///   "orbits":[  // "cluster_functions", prior to 1.0.X
///     {
///       "eci": <float>,
///       "linear_orbit_index": <int>, // "linear_cluster_index", prior to 1.0.X
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
  std::string orbits_name = "orbits";
  std::string index_name = "linear_orbit_index";

  // v0.X compatibility mode:
  if (parser.self.contains("cluster_functions") &&
      !parser.self.contains("orbits")) {
    err_log() << "WARNING: Detected old ECI format (version < 1.0.X): Reading "
                 "\"cluster_functions\" instead of \"orbits\", and "
                 "\"linear_cluster_index\" instead of \"linear_orbit_index\"."
              << std::endl;

    orbits_name = "cluster_functions";
    index_name = "linear_cluster_functions";
  }

  if (parser.self.find(orbits_name) == parser.self.end()) {
    std::stringstream msg;
    msg << "Error parsing ECI: '" << orbits_name << "' array not found";
    parser.error.insert(msg.str());
    return;
  }

  std::vector<double> value;
  std::vector<ECIContainer::size_type> index;

  jsonParser const &orbits_json = parser.self[orbits_name];
  fs::path orbits_path{orbits_name};
  Index i = 0;
  for (auto const &orbit : orbits_json) {
    if (orbit.find("eci") != orbit.end()) {
      double eci_value;
      ECIContainer::size_type linear_orbit_index;
      fs::path orbit_path = orbits_path / std::to_string(i);
      parser.require(eci_value, orbit_path / "eci");
      parser.require(linear_orbit_index, orbit_path / index_name);
      value.push_back(eci_value);
      index.push_back(linear_orbit_index);
    }
    ++i;
  }
  parser.value = notstd::make_unique<ECIContainer>(value.begin(), value.end(),
                                                   index.begin());
}
}  // namespace CASM
