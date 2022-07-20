#include "casm/monte_carlo/io/json/CorrMatchingPotential_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/monte_carlo/CorrMatchingPotential.hh"

namespace CASM {
namespace Monte {

jsonParser &to_json(CorrMatchingTarget const &target, jsonParser &json) {
  json.put_obj();
  json["index"] = target.index;
  json["value"] = target.value;
  json["weight"] = target.weight;
  return json;
}

void from_json(CorrMatchingTarget &target, jsonParser const &json) {
  from_json(target.index, json["index"]);
  from_json(target.value, json["value"]);
  json.get_else(target.weight, "weight", 1.0);
}

jsonParser &to_json(CorrMatchingParams const &params, jsonParser &json) {
  json.put_obj();
  json["tol"] = params.tol;
  json["exact_matching_weight"] = params.exact_matching_weight;
  json["targets"] = params.targets;
  return json;
}

void from_json(CorrMatchingParams &params, jsonParser const &json) {
  json.get_else(params.tol, "tol", CASM::TOL);
  json.get_else(params.exact_matching_weight, "exact_matching_weight", 0.0);
  params.targets.clear();
  json.get_if(params.targets, "targets");
}

jsonParser &to_json(RandomAlloyCorrMatchingParams const &params,
                    jsonParser &json) {
  json.put_obj();
  json["tol"] = params.tol;
  json["exact_matching_weight"] = params.exact_matching_weight;
  json["sublattice_prob"].put_array();
  for (auto const &prob : params.sublattice_prob) {
    jsonParser tjson;
    to_json(prob, tjson, jsonParser::as_array());
    json["sublattice_prob"].push_back(tjson);
  }
  json["random_alloy_targets"] = params.targets;
  return json;
}

void from_json(RandomAlloyCorrMatchingParams &params, jsonParser const &json) {
  json.get_else(params.tol, "tol", CASM::TOL);
  json.get_else(params.exact_matching_weight, "exact_matching_weight", 0.0);
  from_json(params.sublattice_prob, json["sublattice_prob"]);
  params.update_targets();
}

}  // namespace Monte
}  // namespace CASM
