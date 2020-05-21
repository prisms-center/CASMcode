#include "casm/crystallography/io/DoFSetIO.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/crystallography/DoFSet.hh"

namespace CASM {
  jsonParser &to_json(xtal::DoFSet const &_dof, jsonParser &json) {
    json["basis"] = _dof.basis();
    json["axis_names"] = _dof.component_names();
    json["traits"] = _dof.traits().name();
    return json;
  }

  jsonParser &to_json(xtal::SiteDoFSet const &_dof, jsonParser &json) {
    to_json(static_cast<const xtal::DoFSet &>(_dof), json);
    json["excluded_occupants"] = _dof.excluded_occupants();
    return json;
  }

  template <>
  xtal::DoFSet from_json<xtal::DoFSet>(const jsonParser &json) {
    Eigen::MatrixXd basis;
    json.get_if(basis, "basis");

    std::vector<std::string> component_names;
    json.get_if(component_names, "axis_names");

    std::string traits_tag = json["traits"].get<std::string>();

    if(component_names.size()) {
      return xtal::DoFSet(xtal::DoFSet::BasicTraits(traits_tag), component_names, basis);
    }

    return xtal::DoFSet(xtal::DoFSet::BasicTraits(traits_tag));
  }

  template <>
  xtal::SiteDoFSet from_json<xtal::SiteDoFSet>(const jsonParser &json) {
    std::unordered_set<std::string> excluded_occupants;
    json.get_if(excluded_occupants, "excluded_occupants");

    return xtal::SiteDoFSet(from_json<xtal::DoFSet>(json), excluded_occupants);
  }

}
