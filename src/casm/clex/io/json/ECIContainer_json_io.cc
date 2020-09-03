#include <vector>
#include "casm/clex/ECIContainer.hh"
#include "casm/casm_io/json/InputParser_impl.hh"

namespace CASM {

  ///  Make ECIContainer from JSON (eci.json file)
  ///
  /// Format:
  /// \code
  /// {
  ///   "site_functions":[
  ///     {
  ///       "asym_unit": X,
  ///       "sublat_indices: [2, 3],
  ///       "phi_b_0": {"Va":0.0, "O":1.0},
  ///       "phi_b_1": {"Va":0.0, "O":1.0},
  ///        ...
  ///     },
  ///     ...
  ///   ],
  ///   "cluster_functions":[
  ///     {
  ///       "eci": X.XXXXX,
  ///       "prototype_function": "\phi_b_i(s_j)...",
  ///       "orbit": [branch_index, orbit_index],
  ///       "linear_function_index": I,
  ///       "mult": X,
  ///       "prototype": [
  ///         [b, i, j, k],
  ///         ...
  ///       ]
  ///     },
  ///     ...
  ///   ]
  /// }
  /// \endcode
  ///
  void parse(InputParser<ECIContainer> &parser) {

    if(parser.self.find("cluster_functions") == parser.self.end()) {
      parser.error.insert("Error parsing ECI: 'cluster_functions' array not found");
      return;
    }

    std::vector<double> value;
    std::vector<ECIContainer::size_type> index;

    jsonParser const &cluster_functions_json = parser.self["cluster_functions"];
    fs::path cluster_functions_path {"cluster_functions"};
    Index i = 0;
    for(auto const &cluster_json : cluster_functions_json) {
      if(cluster_json.find("eci") != cluster_json.end()) {
        double eci_value;
        ECIContainer::size_type linear_function_index;
        fs::path cluster_path = cluster_functions_path / std::to_string(i);
        parser.require_at(eci_value, cluster_path / "eci");
        parser.require_at(linear_function_index, cluster_path / "linear_function_index");
        value.push_back(eci_value);
        index.push_back(linear_function_index);
      }
      ++i;
    }
    parser.value = notstd::make_unique<ECIContainer>(value.begin(), value.end(), index.begin());
  }
}
