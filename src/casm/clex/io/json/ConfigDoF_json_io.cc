#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/io/json/ConfigDoF_json_io.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  jsonParser &to_json(LocalContinuousConfigDoFValues const &_values, jsonParser &_json) {
    to_json(_values.standard_values().transpose(), _json["values"]);
    return _json;
  }

  void from_json(LocalContinuousConfigDoFValues &_values, jsonParser const &_json) {
    Eigen::MatrixXd tval = _json["values"].get<Eigen::MatrixXd>().transpose();
    _values.from_standard_values(tval);
  }

  jsonParser &to_json(LocalDiscreteConfigDoFValues const &_values, jsonParser &_json) {
    to_json_array(_values.values(), _json);
    return _json;
  }

  void from_json(LocalDiscreteConfigDoFValues &_values, jsonParser const &_json) {
    _values.resize_vol(_json.size() / _values.n_sublat());
    _values.values() = _json.get<LocalDiscreteConfigDoFValues::ValueType>();
  }

  jsonParser &to_json(GlobalContinuousConfigDoFValues const &_values, jsonParser &_json) {
    to_json_array(_values.standard_values(), _json["values"]);
    return _json;
  }

  void from_json(GlobalContinuousConfigDoFValues &_values, jsonParser const &_json) {
    _values.from_standard_values(_json["values"].get<Eigen::VectorXd>());
  }


  std::unique_ptr<ConfigDoF> jsonMake<ConfigDoF>::make_from_json(
    const jsonParser &json, Structure const &structure) {

    auto result_ptr = notstd::make_unique<ConfigDoF>((Index) structure.basis().size(), 0,
                                                     global_dof_info(structure), local_dof_info(structure), structure.occupant_symrep_IDs(),
                                                     structure.lattice().tol());
    result_ptr->from_json(json);
    return result_ptr;
  }

  // ConfigDoF jsonConstructor<ConfigDoF>::from_json(
  //   const jsonParser &json, Index NB, std::map<DoFKey, DoFSetInfo> const &global_info,
  //   std::map<DoFKey, std::vector<DoFSetInfo> > const &local_info,
  //   std::vector<SymGroupRepID> const &_occ_symrep_IDs, double _tol) {
  //
  //   ConfigDoF result(NB, 0, global_info, local_info,  _occ_symrep_IDs, _tol);
  //   result.from_json(json);
  //
  //   return result;
  // }

  ConfigDoF jsonConstructor<ConfigDoF>::from_json(const jsonParser &json, Structure const &structure) {

    ConfigDoF result {(Index) structure.basis().size(), 0, global_dof_info(structure),
                      local_dof_info(structure), structure.occupant_symrep_IDs(), structure.lattice().tol()};
    result.from_json(json);
    return result;
  }

  jsonParser &to_json(const ConfigDoF &value, jsonParser &json) {
    return value.to_json(json);
  }

  void from_json(ConfigDoF &value, const jsonParser &json) { //, Index NB) {
    value.from_json(json);
  }

}
