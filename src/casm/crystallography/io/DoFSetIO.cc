#include "casm/crystallography/io/DoFSetIO.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/crystallography/DoFSet.hh"

namespace CASM {
jsonParser &to_json(xtal::DoFSet const &_dof, jsonParser &json) {
  json["basis"] = _dof.basis().transpose();
  json["axis_names"] = _dof.component_names();
  // json["traits"] = _dof.traits().name();
  return json;
}

jsonParser &to_json(xtal::SiteDoFSet const &_dof, jsonParser &json) {
  to_json(static_cast<const xtal::DoFSet &>(_dof), json);
  if (_dof.excluded_occupants().size()) {
    json["excluded_occupants"] = _dof.excluded_occupants();
  }
  return json;
}

xtal::DoFSet jsonConstructor<xtal::DoFSet>::from_json(
    const jsonParser &json, AnisoValTraits const &traits) {
  if (json.contains("basis") != json.contains("axis_names")) {
    std::stringstream msg;
    msg << "Error reading DoF input: if either is present, 'axis_names' and "
           "'basis' must both be present";
    throw std::runtime_error(msg.str());
  }

  if (!json.contains("axis_names")) {
    return xtal::DoFSet(traits);
  }

  Eigen::MatrixXd row_vector_basis;
  json["basis"].get(row_vector_basis);
  Eigen::MatrixXd basis = row_vector_basis.transpose();

  if (row_vector_basis.cols() != traits.dim()) {
    throw std::runtime_error("Cannot construct DoFSet of type " +
                             traits.name() + ", basis vectors are of size " +
                             std::to_string(row_vector_basis.cols()) +
                             " but must be size " +
                             std::to_string(traits.dim()));
  }

  if (row_vector_basis.rows() > traits.dim()) {
    throw std::runtime_error("Cannot construct DoFSet of type " +
                             traits.name() + ", found " +
                             std::to_string(row_vector_basis.rows()) +
                             " basis vectors, but must be " +
                             std::to_string(traits.dim()) + " or fewer.");
  }

  std::vector<std::string> component_names;
  json["axis_names"].get(component_names);

  if (row_vector_basis.rows() != component_names.size()) {
    throw std::runtime_error(
        "Cannot construct DoFSet of type " + traits.name() + ", found " +
        std::to_string(row_vector_basis.rows()) +
        " basis vectors, but the number of 'axis_names' is " +
        std::to_string(traits.dim()) + ". These numbers must be equal.");
  }

  return xtal::DoFSet(traits, component_names, basis);
}

xtal::SiteDoFSet jsonConstructor<xtal::SiteDoFSet>::from_json(
    const jsonParser &json, AnisoValTraits const &traits) {
  std::unordered_set<std::string> excluded_occupants;
  json.get_if(excluded_occupants, "excluded_occupants");

  return xtal::SiteDoFSet(json.get<xtal::DoFSet>(traits), excluded_occupants);
}

}  // namespace CASM
