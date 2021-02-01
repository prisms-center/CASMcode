#include "casm/clex/io/json/CompositionAxes_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/CompositionAxes.hh"
#include "casm/clex/io/json/CompositionConverter_json_io.hh"

namespace CASM {

void from_json(CompositionAxes &composition_axes, const jsonParser &json) {
  composition_axes = CompositionAxes();

  json.get_if(composition_axes.all_axes, "possible_axes");

  json.get_if(composition_axes.enumerated, "enumerated");

  // Backward compatibility
  std::map<std::string, CompositionConverter> custom;
  json.get_if(custom, "custom_axes");
  composition_axes.all_axes.insert(custom.begin(), custom.end());

  if (json.contains("standard_axes")) {
    std::map<std::string, CompositionConverter> standard;
    from_json(standard, json["standard_axes"]);
    for (auto const &el : standard) {
      composition_axes.enumerated.insert(el.first);
    }
    standard.insert(standard.begin(), standard.end());
  }

  std::string key;
  json.get_if(key, "current_axes");
  if (!key.empty()) {
    composition_axes.select(key);
  }
}

/// \brief Write CompositionAxes to JSON
jsonParser &to_json(CompositionAxes const &composition_axes, jsonParser &json) {
  json = jsonParser::object();

  if (composition_axes.has_current_axes()) {
    json["current_axes"] = composition_axes.curr_key;
  }

  json["possible_axes"] = composition_axes.all_axes;

  if (composition_axes.enumerated.size()) {
    json["enumerated"] = composition_axes.enumerated;
  }

  return json;
}

}  // namespace CASM
