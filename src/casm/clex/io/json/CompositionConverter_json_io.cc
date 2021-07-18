#include "casm/clex/io/json/CompositionConverter_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/CompositionConverter.hh"

namespace CASM {

/// \brief Serialize CompositionConverter to JSON
jsonParser &to_json(const CompositionConverter &f, jsonParser &json) {
  json = jsonParser::object();
  json["components"] = f.components();
  json["independent_compositions"] = f.independent_compositions();
  to_json(f.origin(), json["origin"], CASM::jsonParser::as_array());
  for (int i = 0; i < f.independent_compositions(); i++) {
    to_json(f.end_member(i), json[CompositionConverter::comp_var(i)],
            CASM::jsonParser::as_array());
  }
  json["mol_formula"] = f.mol_formula();
  json["param_formula"] = f.param_formula();

  return json;
}

/// \brief Deserialize CompositionConverter from JSON
void from_json(CompositionConverter &f, const jsonParser &json) {
  std::vector<std::string> components;
  Eigen::VectorXd origin;

  int independent_compositions;

  from_json(components, json["components"]);
  from_json(origin, json["origin"]);

  from_json(independent_compositions, json["independent_compositions"]);
  Eigen::MatrixXd end_members(components.size(), independent_compositions);
  Eigen::VectorXd tvec;
  for (int i = 0; i < independent_compositions; i++) {
    from_json(tvec, json[CompositionConverter::comp_var(i)]);
    end_members.col(i) = tvec;
  }

  f = CompositionConverter(components.begin(), components.end(), origin,
                           end_members);
  return;
}

}  // namespace CASM
