#include "casm/crystallography/io/DoFSetIO.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/crystallography/DoFSet.hh"

namespace CASM {
jsonParser &to_json(xtal::DoFSet const &_dof, jsonParser &json) {
  json["basis"] = _dof.basis();
  json["axis_names"] = _dof.component_names();
  // json["traits"] = _dof.traits().name();
  return json;
}

jsonParser &to_json(xtal::SiteDoFSet const &_dof, jsonParser &json) {
  to_json(static_cast<const xtal::DoFSet &>(_dof), json);
  json["excluded_occupants"] = _dof.excluded_occupants();
  return json;
}

xtal::DoFSet jsonConstructor<xtal::DoFSet>::from_json(
    const jsonParser &json, AnisoValTraits const &traits) {
  Eigen::MatrixXd basis;
  json.get_if(basis, "basis");

  std::vector<std::string> component_names;
  json.get_if(component_names, "axis_names");

  if (component_names.size()) {
    return xtal::DoFSet(traits, component_names, basis);
  }

  return xtal::DoFSet(traits);
}

xtal::SiteDoFSet jsonConstructor<xtal::SiteDoFSet>::from_json(
    const jsonParser &json, AnisoValTraits const &traits) {
  std::unordered_set<std::string> excluded_occupants;
  json.get_if(excluded_occupants, "excluded_occupants");

  return xtal::SiteDoFSet(json.get<xtal::DoFSet>(traits), excluded_occupants);
}

}  // namespace CASM
